"""Verify the input ZIP and optionally extract selected databases."""
import argparse
import csv
import hashlib
import io
import json
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=ROOT / "data/stepwise_diversity_inputs.zip")
    parser.add_argument("--extract", action="store_true")
    parser.add_argument("--report", type=Path)
    parser.add_argument("--dest", type=Path, default=ROOT / "data")
    parser.add_argument("--datasets", nargs="+", choices=["silva", "ltp", "gg2"], default=["silva", "ltp", "gg2"])
    args = parser.parse_args()
    destination = args.dest.resolve()
    report = {"archive": args.archive.name, "datasets": args.datasets, "files": [], "passed": False}
    inventory = json.loads((ROOT / "data/input_archive.json").read_text())
    archive_sha = hashlib.sha256()
    with args.archive.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            archive_sha.update(chunk)
    if archive_sha.hexdigest() != inventory["sha256"] or args.archive.stat().st_size != inventory["bytes"]:
        raise RuntimeError("Archive SHA-256/size differs from data/input_archive.json")
    report["archive_sha256"] = archive_sha.hexdigest()
    with zipfile.ZipFile(args.archive) as archive:
        if archive.read("input_manifest/files.tsv") != (ROOT / "data/input_manifest/files.tsv").read_bytes():
            raise RuntimeError("Embedded manifest differs from the published manifest")
        if archive.read("input_manifest/runs.json") != (ROOT / "config/runs.json").read_bytes():
            raise RuntimeError("Embedded run settings differ from config/runs.json")
        rows = list(csv.DictReader(io.StringIO(archive.read("input_manifest/files.tsv").decode()), delimiter="\t"))
        expected_members = {r["archive_path"] for r in rows} | {"input_manifest/files.tsv", "input_manifest/runs.json"}
        if set(archive.namelist()) != expected_members:
            raise RuntimeError("ZIP membership differs from its manifest")
        for row in rows:
            if row["dataset"] not in args.datasets:
                continue
            target = (destination / row["archive_path"]).resolve()
            if not target.is_relative_to(destination):
                raise RuntimeError("Unsafe archive path: " + row["archive_path"])
            sha = hashlib.sha256()
            size = 0
            writer = None
            partial = target.with_suffix(target.suffix + ".partial")
            if args.extract:
                target.parent.mkdir(parents=True, exist_ok=True)
                writer = partial.open("wb")
            try:
                with archive.open(row["archive_path"]) as handle:
                    for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                        sha.update(chunk)
                        size += len(chunk)
                        if writer:
                            writer.write(chunk)
            finally:
                if writer:
                    writer.close()
            if sha.hexdigest() != row["sha256"] or size != int(row["bytes"]):
                raise RuntimeError("Input checksum/size mismatch: " + row["archive_path"])
            if args.extract:
                partial.replace(target)
            report["files"].append({"path": row["archive_path"], "bytes": size, "sha256": sha.hexdigest(), "passed": True})
            print(f"Verified {row['archive_path']} ({size:,} bytes)", flush=True)
    report["passed"] = True
    if args.extract:
        destination.mkdir(parents=True, exist_ok=True)
        (destination / "input_verification.json").write_text(json.dumps(report, indent=2) + "\n")
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
