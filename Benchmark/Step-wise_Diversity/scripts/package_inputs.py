"""Create the portable collection from the original Step-wise_Diversity directory."""
import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def write_tsv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def copy_file(source, target):
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    args = parser.parse_args()
    source = args.source.resolve()
    configs = json.loads((ROOT / "config/runs.json").read_text())
    provenance = []
    for path in sorted((source / "pcr_simulation").rglob("*")):
        if not path.is_file() or "__pycache__" in path.parts:
            continue
        if path.suffix not in {".py", ".tsv", ".md"}:
            continue
        target = ROOT / "code" / path.relative_to(source)
        copy_file(path, target)
        provenance.append({"bundle_path": str(target.relative_to(ROOT)),
                           "original_path": str(path), "sha256": digest(target)})
    binary = source / "pcr_simulation/usearch11.0.667_i86linux64"
    copy_file(binary, ROOT / "tools" / binary.name)
    (ROOT / "tools" / binary.name).chmod(0o755)
    provenance.append({"bundle_path": "tools/" + binary.name,
                       "original_path": str(binary), "sha256": digest(binary)})
    for pattern in ("*.tsv", "*.nf", "*.config", "*.yml", "*.sh"):
        for path in sorted(source.glob(pattern)):
            target = ROOT / "metadata/legacy_workflow" / path.name
            copy_file(path, target)
            provenance.append({"bundle_path": str(target.relative_to(ROOT)),
                               "original_path": str(path), "sha256": digest(target)})
    for key, config in configs.items():
        run_source = source / config["source_run"]
        target_dir = ROOT / "results/archive_run" / key
        for path in sorted(run_source.rglob("*")):
            if not path.is_file() or path.suffix not in {".tsv", ".png", ".svg"}:
                continue
            # Raw per-record discard reports and unrelated preview backups are unnecessary.
            if path.stat().st_size > 10_000_000 or "_bkp" in path.name:
                continue
            target = target_dir / path.relative_to(run_source)
            copy_file(path, target)
            provenance.append({"bundle_path": str(target.relative_to(ROOT)),
                               "original_path": str(path), "sha256": digest(target)})
        for path in sorted(run_source.parent.iterdir()):
            if path.is_file() and path.suffix in {".tsv", ".png", ".svg"}:
                target = target_dir / "additional_figures" / path.name
                copy_file(path, target)
                provenance.append({"bundle_path": str(target.relative_to(ROOT)),
                                   "original_path": str(path), "sha256": digest(target)})
    (ROOT / "metadata").mkdir(exist_ok=True)
    git = {
        "source_root": str(source),
        "commit": subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip(),
        "status": subprocess.check_output(["git", "-C", str(source), "status", "--short"], text=True),
        "note": "The working-tree scripts, including local edits, are copied without modification.",
    }
    (ROOT / "metadata/source_snapshot.json").write_text(json.dumps(git, indent=2) + "\n")
    write_tsv(ROOT / "metadata/source_files.tsv", provenance)
    archive = ROOT / "data/stepwise_diversity_inputs.zip"
    archive.parent.mkdir(parents=True, exist_ok=True)
    partial = archive.with_suffix(".zip.partial")
    rows = []
    with zipfile.ZipFile(partial, "w", compression=zipfile.ZIP_DEFLATED,
                         compresslevel=3, allowZip64=True) as output:
        for key, config in configs.items():
            entries = [(config["source_input"], config["input"], "aligned_reference")]
            if "header_list" in config:
                entries.append((config["header_list"], config["header_list"], "taxonomy_header_list"))
            for original, relative, role in entries:
                path = source / original
                sha = hashlib.sha256()
                print(f"Packaging {key}: {path.name} ({path.stat().st_size:,} bytes)", flush=True)
                with path.open("rb") as handle, output.open(relative, "w", force_zip64=True) as member:
                    for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                        sha.update(chunk)
                        member.write(chunk)
                rows.append({"dataset": key, "archive_path": relative, "role": role,
                             "bytes": path.stat().st_size, "sha256": sha.hexdigest(),
                             "original_path": str(path), "source_url": config["source_url"]})
        write_tsv(ROOT / "data/input_manifest/files.tsv", rows)
        output.write(ROOT / "data/input_manifest/files.tsv", "input_manifest/files.tsv")
        output.write(ROOT / "config/runs.json", "input_manifest/runs.json")
    partial.replace(archive)
    sha = digest(archive)
    (archive.parent / (archive.name + ".sha256")).write_text(f"{sha}  {archive.name}\n")
    (ROOT / "data/input_archive.json").write_text(json.dumps({
        "archive": archive.name, "bytes": archive.stat().st_size, "sha256": sha,
        "payload_files": len(rows), "uncompressed_payload_bytes": sum(r["bytes"] for r in rows),
    }, indent=2) + "\n")
    print(f"Created {archive} ({archive.stat().st_size:,} bytes)", flush=True)


if __name__ == "__main__":
    main()
