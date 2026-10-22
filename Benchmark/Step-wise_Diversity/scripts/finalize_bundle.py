"""Collect small verification artifacts and create the code/results upload ZIP."""
import argparse
import csv
import datetime
import hashlib
import json
import shutil
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def collect(source, target):
    for path in source.rglob("*"):
        if not path.is_file() or path.suffix not in {".json", ".tsv", ".png", ".svg", ".log"}:
            continue
        if any(part in {".matplotlib", "collection"} for part in path.relative_to(source).parts):
            continue
        if path.stat().st_size > 10_000_000:
            continue
        destination = target / path.relative_to(source)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, destination)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--test-root", type=Path, required=True)
    args = parser.parse_args()
    test_root = args.test_root.resolve()
    evidence = ROOT / "results/validation_checks"
    evidence.mkdir(parents=True, exist_ok=True)
    checks = {}
    for name, relative in {
        "archive": "inputs/input_verification.json",
        "isolated_commands": "relocated-confirmed/test_report.json",
        "isolated_smoke": "relocated-confirmed/smoke/validation.json",
        "full_ltp": "full_historical/validation.json",
        "reference_populations": "population-historical/reference_population_validation.json",
        "gg2_species": "gg2-species-historical/species_summary_validation.json",
    }.items():
        report = json.loads((test_root / relative).read_text())
        if report.get("passed") is not True:
            raise RuntimeError("Validation did not pass: " + name)
        checks[name] = {"passed": True, "original_report": relative}
    summary = {"checked_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
               "checks": checks, "passed": True,
               "full_window_clustering_not_recomputed": ["silva", "gg2"]}
    (evidence / "validation_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    sources = {"relocated": "relocated-confirmed", "full_ltp": "full_historical", "populations": "population-historical",
               "gg2_species": "gg2-species-historical"}
    for label, relative in sources.items():
        collect(test_root / relative, evidence / label)
    diagnostics = evidence / "diagnostics"
    diagnostics.mkdir(exist_ok=True)
    for original, name in {
        "full/validation.json": "ltp_four_threads_current_taxonomy.json",
        "gg2-species-full/species_summary_validation.json": "gg2_two_threads.json",
        "population-audit/reference_population_validation.json": "uniform_current_taxonomy.json",
        "relocated-historical/test_report.json": "parallel_index_race.json",
    }.items():
        shutil.copy2(test_root / original, diagnostics / name)
    shutil.copy2(test_root / "inputs/input_verification.json", evidence / "input_verification.json")
    manifest = ROOT / "metadata/bundle_files.tsv"
    files = []
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file() or path == manifest:
            continue
        relative = path.relative_to(ROOT)
        if path.suffix in {".zip", ".partial", ".pyc"} or "__pycache__" in relative.parts:
            continue
        if relative.parts[:2] == ("data", "DBs"):
            continue
        if any(part.startswith("reproduced") or part == "smoke_run" for part in relative.parts):
            continue
        files.append(path)
    with manifest.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["path", "bytes", "sha256"])
        for path in files:
            writer.writerow([str(path.relative_to(ROOT)), path.stat().st_size, sha256(path)])
    archive = ROOT.parent / "Step-wise_Diversity_code_and_results.zip"
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as output:
        for path in [*files, manifest]:
            output.write(path, str(Path(ROOT.name) / path.relative_to(ROOT)))
    (archive.parent / (archive.name + ".sha256")).write_text(f"{sha256(archive)}  {archive.name}\n")
    print(json.dumps({"small_archive": str(archive), "bytes": archive.stat().st_size,
                      "sha256": sha256(archive), "files": len(files) + 1}, indent=2))


if __name__ == "__main__":
    main()
