#!/usr/bin/env python3
"""Write source links, small-file provenance and dependency-audit tables."""

import argparse
import ast
import csv
import json
import re
import shutil
import stat
import sys
from collections import Counter, defaultdict
from pathlib import Path

from survey_io import read_tsv, sha256, write_tsv


BUNDLE = Path(__file__).resolve().parents[1]


def audit(root):
    result = []
    for path in sorted(root.rglob("*.py")):
        if "__pycache__" in path.parts:
            continue
        tree = ast.parse(path.read_text())
        imports = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imports.update(alias.name.split(".")[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and not node.level and node.module:
                imports.add(node.module.split(".")[0])
        for name in sorted(imports):
            result.append({"file": str(path.relative_to(root)), "import_root": name,
                           "standard_library": name in sys.stdlib_module_names,
                           "note": "Static audit includes optional imports; see environment/README.md"})
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    args = parser.parse_args()
    workspace = args.workspace.resolve()
    benchmark = BUNDLE.parent / "tacticbench_TIC_Archaea_mock4_mock5_SRP"
    links = {r["community_id"]: r for r in read_tsv(benchmark / "metadata/FASTQ_and_publication_links.tsv")}
    inputs = read_tsv(BUNDLE / "data/input_manifest/inputs.tsv")
    metadata = BUNDLE / "metadata"
    metadata.mkdir(exist_ok=True)
    sources = []
    copied = []

    def copy(source, relative):
        destination = BUNDLE / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists():
            destination.chmod(destination.stat().st_mode | stat.S_IWUSR)
        shutil.copy2(source, destination)
        destination.chmod(destination.stat().st_mode | stat.S_IWUSR)
        copied.append({"bundle_path": str(relative), "original_path": str(source),
                       "size_bytes": source.stat().st_size, "sha256": sha256(destination)})

    for item in inputs:
        source = Path(item["source_path"])
        accessions = sorted(set(re.findall(r"[SED]RR\d+", source.name)))
        parts = source.parts
        mock = next((part for part in parts if part.startswith("mock-")), "")
        if not mock:
            mock = next(("mock-" + part.removeprefix("V4_") for part in parts
                         if part in {"ZIEL_1", "ZIEL_2", "ZYMO", "V4_ZIEL_1", "V4_ZIEL_2", "V4_ZYMO"}), "")
        source_links = "; ".join("https://www.ncbi.nlm.nih.gov/sra/?term=" + acc for acc in accessions)
        if not source_links and mock in links:
            source_links = links[mock]["fastq_source"]
        sources.append({"input_id": item["input_id"], "kind": item["kind"], "community_from_path": mock,
                        "accession_in_filename": "; ".join(accessions), "read_source_urls": source_links,
                        "publication_url": links.get(mock, {}).get("publication_doi_or_url", ""),
                        "source_path": item["source_path"],
                        "note": "Links identify sequencing origin, not a download of the exact processed payload; inherited from filenames/existing metadata, not re-verified online"})

    write_tsv(metadata / "input_sources.tsv", sources)
    copy(benchmark / "metadata/FASTQ_and_publication_links.tsv", Path("metadata/benchmark_community_source_links.tsv"))
    for root in [workspace / "Mocks/Picked/mockrobiota/data", workspace / "Mocks/Picked/Isabel"]:
        for folder in sorted(root.iterdir()):
            if not folder.is_dir():
                continue
            for name in ("dataset-metadata.tsv", "sample-metadata.tsv", "README.md", "SraAccList.txt"):
                path = folder / name
                if path.is_file() and path.stat().st_size < 2 * 1024 * 1024:
                    copy(path, Path("metadata/source_metadata") / folder.name / name)
    hc = workspace / "Mocks/Picked/mockrobiota/data/mock-HC227_V3-V4"
    for name in ("biological_metadata.tsv", "technical_metadata.tsv"):
        if (hc / name).is_file():
            copy(hc / name, Path("metadata/source_metadata/mock-HC227_V3-V4") / name)

    parents = {Path(row["source_path"]).parent for row in inputs}
    context = []
    for parent in sorted(parents):
        for name in ("MANIFEST", "metadata.yml", "TACTICPipeline_args.yml", "spike_stat_mapping_file.csv"):
            source = parent / name
            if source.is_file() and source.stat().st_size < 2 * 1024 * 1024:
                relative = Path("metadata/input_context") / source.relative_to(workspace)
                copy(source, relative)
                context.append({"bundle_path": str(relative), "original_path": str(source)})
    write_tsv(metadata / "copied_metadata_manifest.tsv", copied)
    write_tsv(metadata / "dependency_audit.tsv", audit(BUNDLE / "code/length_survey"))
    write_tsv(benchmark / "metadata/dependency_audit.tsv", audit(benchmark / "code/tacticbench/src"))
    snapshot = []
    for path in sorted((BUNDLE / "code/length_survey").glob("*.py")):
        source = workspace / "length_survey" / path.name
        snapshot.append({"bundle_path": str(path.relative_to(BUNDLE)), "original_path": str(source),
                         "sha256": sha256(path), "identical_to_source": sha256(path) == sha256(source)})
    write_tsv(metadata / "original_scripts_manifest.tsv", snapshot)
    observations = read_tsv(BUNDLE / "data/input_manifest/observations.tsv")
    findings = read_tsv(BUNDLE / "data/input_manifest/findings.tsv")
    recovered = Counter(r["script"] for r in observations)
    missing = Counter(r["script"] for r in findings)
    coverage = [{"script": name, "recovered_observations": recovered[name],
                 "unavailable_observations": missing[name],
                 "recorded_total": recovered[name] + missing[name]}
                for name in sorted(set(recovered) | set(missing))]
    write_tsv(metadata / "survey_coverage.tsv", coverage)
    grouped = defaultdict(list)
    for row in findings:
        analysis = re.search(r"Analysis_[^/;\t]+", row["input"])
        grouped[analysis[0] if analysis else row["input"]].append(row)
    missing_rows = [{"input_identifier": name, "unavailable_observations": len(entries),
                     "scripts": "; ".join(sorted({r["script"] for r in entries})),
                     "evidence": "; ".join(sorted({r["evidence"] for r in entries}))}
                    for name, entries in sorted(grouped.items())]
    write_tsv(metadata / "missing_input_summary.tsv", missing_rows)
    package = []
    for path in sorted(BUNDLE.rglob("*")):
        if not path.is_file() or "__pycache__" in path.parts or "payloads" in path.parts:
            continue
        relative = path.relative_to(BUNDLE)
        if path.name.endswith((".zip", ".partial")) or str(relative).startswith("results/reproduced_run/"):
            continue
        if path.name == "small_files_manifest.tsv":
            continue
        package.append({"bundle_path": str(relative), "bytes": path.stat().st_size, "sha256": sha256(path)})
    write_tsv(metadata / "small_files_manifest.tsv", package)
    print(json.dumps({"input_sources": len(sources), "original_scripts": len(snapshot),
                      "small_metadata_files": len(copied), "unchanged_original_scripts": all(r["identical_to_source"] for r in snapshot)}))


if __name__ == "__main__":
    main()
