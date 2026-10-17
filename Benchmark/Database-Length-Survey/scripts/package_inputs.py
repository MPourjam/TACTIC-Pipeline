#!/usr/bin/env python3
"""Recover the length-survey inputs using saved reports, then build a data ZIP."""

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import threading
import zipfile
from collections import defaultdict
from concurrent.futures import Future, ThreadPoolExecutor, as_completed
from pathlib import Path

from survey_io import expected_matches, scan, sha256, write_tsv


BUNDLE = Path(__file__).resolve().parents[1]
OLD_ROOT = "/media/adminuser/MSD-Storage/Mohsen/TACTIC_Benchmark"
OLD_PICKED = "Mocks/Picked DO NOT ENTER ME DO NOT ENTER ME DO NOT ENTER ME"
MOCKS_ROOT = "Mocks/Benchmark_git/All_Analyses/mocks"


def rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def expected_fastq(row):
    if int(row["reads"]) == 0:
        return {"count": 0, "unique_lengths": 0}
    return {
        "count": int(row["reads"]), "min": int(row["min_length"]),
        "max": int(row["max_length"]), "mean_2dp": row["mean_length"],
        "median_2dp": row["median_length"], "stdev_2dp": row["stdev"],
        "unique_lengths": int(row["unique_lengths"]),
    }


def expected_fasta(row):
    return {
        "count": int(row.get("Sequence_Count", row.get("ZOTU_Count"))),
        "min": int(row.get("Min_Length", row.get("Min_Len"))),
        "max": int(row.get("Max_Length", row.get("Max_Len"))),
        "mean_2dp": row.get("Mean_Length", row.get("Mean_Len")),
        "median": float(row.get("Median_Length", row.get("Median_Len"))),
        "stdev_2dp": row["Stdev"], "q1": int(row["Q1"]), "q3": int(row["Q3"]),
    }


def collect(workspace):
    print("Indexing surviving FASTQ and zOTU paths", flush=True)
    inputs = {}
    observations = []
    findings = []
    cache_file = Path("/tmp/length_survey_input_cache.json")
    recovery_cache = json.loads(cache_file.read_text()) if cache_file.exists() else {}
    fastq_index = defaultdict(list)
    fasta_index = defaultdict(list)
    for current, dirs, files in os.walk(workspace / "Mocks"):
        dirs[:] = [d for d in dirs if d not in {".git", "__pycache__", "reproducibility"}]
        for name in files:
            path = Path(current) / name
            if name.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz")) and path.is_file():
                fastq_index[name].append(path)
            elif name == "ZOTUs-Seqs.fasta":
                fasta_index[path.parent.name].append(path)
    print(f"Indexed {sum(map(len, fastq_index.values()))} FASTQ paths and {sum(map(len, fasta_index.values()))} zOTU paths", flush=True)

    def add(path, script, logical, mock, expected=None, limit=0, evidence="", note=""):
        if path is None:
            findings.append({"severity": "missing", "script": script, "input": logical,
                             "evidence": evidence, "detail": note})
            return
        path = path.absolute()
        key = str(path)
        if key not in inputs:
            inputs[key] = {
                "input_id": f"input_{len(inputs) + 1:04d}", "source_path": key,
                "kind": "fasta" if path.suffix == ".fasta" else "fastq",
                "size_bytes": path.stat().st_size,
            }
        observations.append({
            "observation_id": f"observation_{len(observations) + 1:04d}",
            "script": script, "input_id": inputs[key]["input_id"],
            "original_path_or_identifier": logical, "mock": mock, "max_records": limit,
            "expected_json": json.dumps(expected or {}, sort_keys=True),
            "evidence": evidence, "recovery_note": note,
        })

    def recover_fastq(logical, expected=None, limit=0, preferred=None):
        relative = logical.removeprefix(OLD_ROOT + "/")
        candidates = [workspace / relative, workspace / relative.replace(OLD_PICKED, "Mocks/Picked")]
        if preferred:
            candidates.extend(preferred)
        if "/mocks/raw_data/" in logical:
            suffix = logical.split("/mocks/raw_data/")[1]
            candidates.append(workspace / "Mocks/MOCKROBIOTA/data" / suffix)
            candidates.append(workspace / "Mocks/Picked/All_Analyses" / suffix)
            mock, *parts = Path(suffix).parts
            if mock in {"mock-ZIEL_1", "mock-ZIEL_2", "mock-ZYMO"}:
                candidates.append(workspace / MOCKS_ROOT / "cleaned_up_202605" / mock / "untrimmed_fastq" / Path(suffix).name)
        candidates.extend(sorted(fastq_index[Path(logical).name], key=lambda p: ("cleaned_up" in str(p), len(str(p)), str(p))))
        alternate = Path(logical).name.replace("_1_L001_R1_", "_0_L001_R1_")
        if alternate != Path(logical).name:
            candidates.extend(sorted(fastq_index[alternate], key=lambda p: ("/Picked/All_Analyses/" not in str(p), str(p))))
        if expected and "/mocks/raw_data/" in logical:
            suffix = logical.split("/mocks/raw_data/")[1]
            mock = Path(suffix).parts[0]
            name = Path(logical).name
            mate_match = re.search(r"_R([12])", name) or re.search(r"_([12])\.fastq", name)
            accession = re.search(r"[SED]RR\d+", name)
            sample = re.search(r"Even[12]|Staggered[12]", name, re.I)
            aliases = []
            for candidate_name, paths in fastq_index.items():
                candidate_mate = re.search(r"_R([12])", candidate_name) or re.search(r"_([12])\.fastq", candidate_name)
                if not mate_match or not candidate_mate or mate_match[1] != candidate_mate[1]:
                    continue
                if accession and accession[0] not in candidate_name:
                    continue
                if sample and sample[0].lower() not in candidate_name.lower():
                    continue
                for path in paths:
                    if accession or mock in path.parts:
                        aliases.append(path)
            candidates.extend(sorted(aliases, key=lambda p: ("/trimmed_fastq/" not in str(p),
                              "/cleaned_up_202605/" not in str(p), str(p))))
        available = list(dict.fromkeys(p for p in candidates if p.is_file()))
        if not available:
            return None, "No surviving file with the recorded input name"
        if expected:
            for path in available:
                stat = path.stat()
                key = f"{path.resolve()}|{stat.st_size}|{stat.st_mtime_ns}"
                cached = recovery_cache.get(key)
                if cached:
                    stats = cached["samples"].get(str(limit)) if limit else cached["stats"]
                elif limit:
                    stats = scan(path, "fastq", limit)["stats"]
                else:
                    continue
                if stats and not expected_matches(stats, expected):
                    return path, "Recorded name/relocation plus matching saved length statistics; original byte identity not asserted"
            if limit:
                return None, "Surviving same-name candidates do not match the saved sampled length statistics"
        note = "Recovered by recorded name and source-tree relocation; historical statistics checked separately"
        if available[0].name != Path(logical).name:
            note = "Saved summary names R1 with lane index 1; surviving R1 uses lane index 0; prefix length statistics validated separately"
        return available[0], note

    for report, script in [
        ("mockrobiota_fastq_reports/fastq_file_report_index.csv", "generate_fastq_file_reports"),
        ("raw_data_trimmed_fastq_reports/trimmed_fastq_file_report_index.csv", "generate_trimmed_fastq_only_reports"),
    ]:
        for index, row in enumerate(rows(workspace / "length_survey" / report), 2):
            if index % 25 == 0:
                print(f"Resolving {script}: row {index}", flush=True)
            expected = expected_fastq(row)
            path, note = recover_fastq(row["file"], expected)
            mock = row.get("mock_group", row.get("mock"))
            add(path, script, row["file"], mock, expected, evidence=f"length_survey/{report}:{index}", note=note)
    sampled_report = workspace / "Mocks/Benchmark_git/All_Analyses/mocks/trimmed_fastq_length_detailed.csv"
    for index, row in enumerate(rows(sampled_report), 2):
        mock = row["Mock_Community"]
        logical = f"{OLD_ROOT}/{MOCKS_ROOT}/raw_data/{mock}/trimmed_fastq/{row['Fastq_File']}"
        expected = {"count": int(row["Total_Reads_Sampled"]), "min": int(row["Min_Length"]),
                    "max": int(row["Max_Length"]), "mean_2dp": row["Mean_Length"],
                    "median": float(row["Median_Length"]), "stdev_2dp": row["Stdev"]}
        path, note = recover_fastq(logical, expected, 10000)
        add(path, "analyze_trimmed_fastq_lengths", logical, mock, expected, 10000,
            f"{sampled_report.relative_to(workspace)}:{index}", note)

    untrimmed_report = workspace / "length_survey/untrimmed_fastq_read_length_summary.csv"
    for index, row in enumerate(rows(untrimmed_report), 2):
        for mate in ("r1", "r2"):
            logical = f"{OLD_ROOT}/{MOCKS_ROOT}/raw_data/{row['mock']}/untrimmed_fastq/{row[mate + '_file']}"
            expected = {"count": int(row[mate + "_reads_sampled"]),
                        "min": int(row[mate + "_length_min"]), "max": int(row[mate + "_length_max"])}
            path, note = recover_fastq(logical, expected, 1000)
            add(path, "analyze_untrimmed_fastq_lengths", logical, row["mock"], expected, 1000,
                f"length_survey/untrimmed_fastq_read_length_summary.csv:{index}", note)

    raw_report = workspace / "raw_read_length_analysis.csv"
    raw_expected = {(r["Mock_Community"], r["Fastq_File"]): r for r in rows(raw_report)}
    for mock in ("ZIEL_1", "ZIEL_2", "ZYMO"):
        folder = workspace / "Mocks/Picked/Isabel" / mock / "fastqs"
        for path in sorted(folder.glob("*.fastq.gz")):
            if "__processed" in str(path):
                continue
            logical = str(path).replace(str(workspace / "Mocks/Picked"), f"{OLD_ROOT}/{OLD_PICKED}")
            saved = raw_expected.get((mock, path.name))
            expected = {"count": int(saved["Total_Reads_Sampled"]), "min": int(saved["Min_Length"]),
                        "max": int(saved["Max_Length"]), "unique_lengths": int(saved["Unique_Lengths"])} if saved else {}
            add(path, "analyze_raw_read_lengths", logical, mock, expected, limit=10000,
                evidence="analyze_raw_read_lengths.py:102-114; relocated Picked directory")
            if mock == "ZIEL_1":
                add(path, "analyze_ziel1_comprehensive", logical, mock,
                    evidence="analyze_ziel1_comprehensive.py:38-41; relocated Picked directory")

    fasta_stats = {}
    print("Resolving saved zOTU measurements", flush=True)
    for paths in fasta_index.values():
        for path in paths:
            fasta_stats[path] = scan(path, "fasta")["stats"]

    for report, script in [
        ("Mocks/Benchmark_git/All_Analyses/zotu_comparison_tactic_vs_raw.csv", "compare_zotu_locations"),
        ("Mocks/zotu_length_stats_detailed.csv", "generate_zotu_stats"),
    ]:
        usage = defaultdict(set)
        for index, row in enumerate(rows(workspace / report), 2):
            expected = expected_fasta(row)
            candidates = [p for p in fasta_index[row["Analysis_Dir"]]
                          if not expected_matches(fasta_stats[p], expected)]
            group = row["Mock_Name"]
            if script == "compare_zotu_locations":
                location = row["Location"]
                relative = f"{MOCKS_ROOT}/" + ("raw_data/" if location == "raw_data" else "")
                relative += group + ("/tactic_results" if location == "tactic_results" else "")
                relative += f"/{row['Analysis_Dir']}/ZOTUs-Seqs.fasta"
                logical = f"{OLD_ROOT}/{relative}"
                prefix = "/old_before_202605/" if location == "tactic_results" else "/Picked/All_Analyses/"
                candidates.sort(key=lambda p: (prefix not in str(p), str(p)))
                if location == "raw_data":
                    # Prefer retained May-6 inputs, not an older analysis with the same name.
                    candidates.sort(key=lambda p: ("/cleaned_up_202605_test/" not in str(p), "/mocks-storage/raw_data_202605/" not in str(p), str(p)))
            else:
                logical = f"saved row {index}: group={group}; analysis={row['Analysis_Dir']}"
                prefix = "/Picked/mockrobiota/data/" if group == "Unknown" else (
                    "/Mocks/Benchmark_git/" if group == "mocks" else "/Picked/")
                candidates.sort(key=lambda p: (p in usage[(group, row["Analysis_Dir"])], prefix not in str(p), str(p)))
            path = candidates[0] if candidates else None
            if path:
                usage[(group, row["Analysis_Dir"])].add(path)
            note = "Matched recorded analysis directory and all nine saved length statistics"
            if script == "generate_zotu_stats":
                note += "; original full file path was not saved, so original tree position is not asserted"
            if path is None:
                note = "No surviving zOTU FASTA with the recorded analysis-directory name and saved length statistics"
            add(path, script, logical, group, expected, evidence=f"{report}:{index}", note=note)

    return list(inputs.values()), observations, findings


def copy_small_files(workspace):
    code = BUNDLE / "code/length_survey"
    code.mkdir(parents=True, exist_ok=True)
    for path in (workspace / "length_survey").glob("*.py"):
        shutil.copy2(path, code / path.name)
    outputs = BUNDLE / "results/archive_run/length_survey"
    shutil.copytree(workspace / "length_survey", outputs, dirs_exist_ok=True,
                    ignore=shutil.ignore_patterns("*.py", "__pycache__"))
    extra = [workspace / "raw_read_length_analysis.csv"]
    extra += list((workspace / "Mocks").glob("zotu_length_stats*"))
    extra += list((workspace / "Mocks/Benchmark_git/All_Analyses").glob("zotu_comparison*.csv"))
    extra += list((workspace / "Mocks/Benchmark_git/All_Analyses/mocks").glob("trimmed_fastq_length*.csv"))
    for path in extra:
        if path.is_file():
            dest = BUNDLE / "results/archive_run/other_locations" / path.name
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, dest)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--inventory-only", action="store_true")
    parser.add_argument("--allow-missing", action="store_true", help="Archive surviving inputs and explicitly retain the missing-input inventory")
    args = parser.parse_args()
    workspace = args.workspace.resolve()
    copy_small_files(workspace)
    inputs, observations, findings = collect(workspace)
    manifest = BUNDLE / "data/input_manifest"
    manifest.mkdir(parents=True, exist_ok=True)
    write_tsv(manifest / "observations.tsv", observations)
    write_tsv(manifest / "inputs.tsv", inputs)
    write_tsv(manifest / "findings.tsv", findings, ["severity", "script", "input", "evidence", "detail"])
    print(json.dumps({"source_paths": len(inputs), "observations": len(observations),
                      "logical_bytes": sum(int(i["size_bytes"]) for i in inputs), "missing": len(findings)}), flush=True)
    if args.inventory_only:
        return
    if findings and not args.allow_missing:
        raise RuntimeError("Missing historical inputs; see findings.tsv before packaging")

    cache_path = Path("/tmp/length_survey_input_cache.json")
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    metrics = {}
    known_payloads = {result["sha256"]: result for result in cache.values()}
    pending = {}
    lock = threading.Lock()

    def inspect(item):
        path = Path(item["source_path"])
        stat = path.stat()
        key = f"{path.resolve()}|{stat.st_size}|{stat.st_mtime_ns}"
        if key in cache:
            return item["input_id"], key, cache[key]
        digest = sha256(path)
        with lock:
            if digest in known_payloads:
                return item["input_id"], key, known_payloads[digest]
            owner = digest not in pending
            future = pending.setdefault(digest, Future())
        if owner:
            try:
                process = subprocess.run([sys.executable, str(BUNDLE / "scripts/survey_io.py"),
                                          str(path), item["kind"], digest], text=True, capture_output=True, check=True)
                result = json.loads(process.stdout)
                with lock:
                    known_payloads[digest] = result
                future.set_result(result)
            except BaseException as exc:
                future.set_exception(exc)
                raise
        return item["input_id"], key, future.result()

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(inspect, item) for item in inputs]
        for done, future in enumerate(as_completed(futures), 1):
            input_id, key, result = future.result()
            cache[key] = result
            metrics[input_id] = result
            cache_path.write_text(json.dumps(cache))
            if done % 20 == 0:
                print(f"Read and validated {done}/{len(inputs)} source paths", flush=True)

    payloads = {}
    for item in inputs:
        result = metrics[item["input_id"]]
        item["sha256"] = result["sha256"]
        suffix = ".fasta" if item["kind"] == "fasta" else (
            ".fastq.gz" if item["source_path"].endswith(".gz") else ".fastq")
        item["data_path"] = ("small_inputs/" if item["kind"] == "fasta" else "payloads/") + item["sha256"] + suffix
        payloads.setdefault(item["data_path"], item)
        if item["kind"] == "fasta":
            destination = BUNDLE / "data" / item["data_path"]
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(item["source_path"], destination)

    validation = []
    rejected = set()
    for observation in observations:
        result = metrics[observation["input_id"]]
        limit = str(observation["max_records"])
        actual = result["samples"][limit] if int(limit) else result["stats"]
        expected = json.loads(observation["expected_json"])
        differences = expected_matches(actual, expected)
        if differences:
            rejected.add(observation["observation_id"])
            findings.append({"severity": "unrecoverable", "input": observation["original_path_or_identifier"],
                             "script": observation["script"], "evidence": observation["evidence"],
                             "detail": "Surviving candidate rejected: saved length statistics differ in " + ",".join(differences)})
        validation.append({"observation_id": observation["observation_id"], "script": observation["script"],
                           "input_id": observation["input_id"], "status": "passed" if not differences else "candidate_rejected",
                           "fields_differing": ",".join(differences), "stats_json": json.dumps(actual, sort_keys=True)})
    candidates = [{**o, "candidate_source_path": next(i["source_path"] for i in inputs if i["input_id"] == o["input_id"])}
                  for o in observations if o["observation_id"] in rejected]
    if candidates:
        write_tsv(manifest / "rejected_recovery_candidates.tsv", candidates)
    observations = [o for o in observations if o["observation_id"] not in rejected]
    used = {o["input_id"] for o in observations}
    inputs = [i for i in inputs if i["input_id"] in used]
    payloads = {i["data_path"]: i for i in inputs}
    metrics = {key: value for key, value in metrics.items() if key in used}
    write_tsv(manifest / "inputs.tsv", inputs)
    write_tsv(manifest / "observations.tsv", observations)
    write_tsv(manifest / "findings.tsv", findings, ["severity", "script", "input", "evidence", "detail"])
    write_tsv(manifest / "historical_validation.tsv", validation)
    (manifest / "length_distributions.json").write_text(json.dumps(metrics, indent=2) + "\n")
    if rejected and not args.allow_missing:
        raise RuntimeError("Historical-statistic mismatches: see historical_validation.tsv; no archive finalized")

    archive_path = BUNDLE / "data/length_survey_inputs.zip"
    partial = archive_path.with_suffix(".zip.partial")
    members = {name: Path(item["source_path"]) for name, item in payloads.items() if item["kind"] == "fastq"}
    members.update({"input_manifest/" + p.name: p for p in manifest.iterdir() if p.is_file() and p.name != "SHA256SUMS"})
    checksums = {name: (payloads[name]["sha256"] if name in payloads else sha256(path)) for name, path in members.items()}
    (manifest / "SHA256SUMS").write_text("".join(f"{digest}  {name}\n" for name, digest in sorted(checksums.items())))
    members["input_manifest/SHA256SUMS"] = manifest / "SHA256SUMS"
    with zipfile.ZipFile(partial, "w", allowZip64=True, compresslevel=1) as archive:
        for index, (name, path) in enumerate(sorted(members.items()), 1):
            compression = zipfile.ZIP_STORED if name.endswith(".gz") else zipfile.ZIP_DEFLATED
            archive.write(path, name, compress_type=compression)
            if index % 20 == 0:
                print(f"Archived {index}/{len(members)} files", flush=True)
    from verify_input_archive import verify
    report = verify(partial)
    partial.replace(archive_path)
    digest = sha256(archive_path)
    archive_path.with_suffix(".zip.sha256").write_text(f"{digest}  {archive_path.name}\n")
    report.update({"archive": archive_path.name, "archive_bytes": archive_path.stat().st_size,
                   "archive_sha256": digest, "source_paths": len(inputs), "observations": len(observations),
                   "unique_fastq_payloads": sum(i["kind"] == "fastq" for i in payloads.values()),
                   "unique_fasta_payloads": sum(i["kind"] == "fasta" for i in payloads.values()),
                   "historical_statistic_checks": "all included observations passed; unavailable inputs documented separately",
                   "historical_measurements_not_recovered": len(findings),
                   "recovery_candidates_rejected": len(rejected),
                   "historical_reproduction_complete": not findings})
    (BUNDLE / "data/input_archive_validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
