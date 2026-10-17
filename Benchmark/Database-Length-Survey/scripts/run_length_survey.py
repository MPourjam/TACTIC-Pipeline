#!/usr/bin/env python3
"""Repeat the frozen FASTQ/FASTA length measurements without original absolute paths."""

import argparse
import gzip
import json
import zipfile
from contextlib import ExitStack, contextmanager
from pathlib import Path

from survey_io import expected_matches, read_tsv, scan, write_tsv


BUNDLE = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--script", default="all", help="Original script stem, or all")
    parser.add_argument("--archive", type=Path, help="Read FASTQs inside ZIP instead of extracting")
    parser.add_argument("--data-dir", type=Path, default=BUNDLE / "data")
    parser.add_argument("--output-dir", type=Path, default=BUNDLE / "results/reproduced_run")
    parser.add_argument("--max-inputs", type=int, help="Limit distinct inputs for a smoke test only")
    args = parser.parse_args()
    if args.max_inputs is not None and args.max_inputs <= 0:
        parser.error("--max-inputs must be positive")
    observations = read_tsv(args.data_dir / "input_manifest/observations.tsv")
    available = sorted({row["script"] for row in observations})
    if args.script != "all" and args.script not in available:
        parser.error("Choose all or " + ", ".join(available))
    observations = [r for r in observations if args.script == "all" or r["script"] == args.script]
    required = {r["input_id"] for r in observations}
    inputs = [r for r in read_tsv(args.data_dir / "input_manifest/inputs.tsv") if r["input_id"] in required]
    if args.max_inputs:
        inputs = inputs[:args.max_inputs]
    selected = {r["input_id"] for r in inputs}
    observations = [r for r in observations if r["input_id"] in selected]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    measurements = {}
    distributions = {}
    measured_rows = []

    with ExitStack() as stack:
        archive = stack.enter_context(zipfile.ZipFile(args.archive)) if args.archive else None

        @contextmanager
        def open_payload(path, mode):
            with archive.open(path.as_posix()) as raw:
                if path.suffix == ".gz":
                    with gzip.GzipFile(fileobj=raw) as handle:
                        yield handle
                else:
                    yield raw

        # Identical payloads are measured once, but original observations stay separate.
        cache = {}
        for index, item in enumerate(inputs, 1):
            relative = Path(item["data_path"])
            if relative.is_absolute() or ".." in relative.parts:
                raise ValueError(f"Unsafe data path: {relative}")
            if item["sha256"] not in cache:
                if archive and item["kind"] == "fastq":
                    result = scan(relative, item["kind"], open_file=open_payload)
                else:
                    result = scan(args.data_dir / relative, item["kind"])
                cache[item["sha256"]] = result
            result = cache[item["sha256"]]
            measurements[item["input_id"]] = result
            distributions[item["input_id"]] = result["distribution"]
            measured_rows.append({"input_id": item["input_id"], "data_path": item["data_path"], **result["stats"]})
            if index % 20 == 0:
                print(f"Measured {index}/{len(inputs)} input paths", flush=True)

    results = []
    for row in observations:
        result = measurements[row["input_id"]]
        stats = result["samples"][row["max_records"]] if int(row["max_records"]) else result["stats"]
        expected = json.loads(row["expected_json"])
        differences = expected_matches(stats, expected)
        results.append({"observation_id": row["observation_id"], "script": row["script"],
                        "input_id": row["input_id"], "mock": row["mock"],
                        "max_records": row["max_records"], "historical_check": "passed" if not differences else "mismatch",
                        "expected_available": bool(expected), "fields_differing": ",".join(differences), **stats})
    write_tsv(args.output_dir / "input_statistics.tsv", measured_rows)
    write_tsv(args.output_dir / "survey_observations.tsv", results)
    for script in sorted({r["script"] for r in results}):
        write_tsv(args.output_dir / f"{script}.tsv", [r for r in results if r["script"] == script])
    (args.output_dir / "length_distributions.json").write_text(json.dumps(distributions, indent=2) + "\n")
    findings = read_tsv(args.data_dir / "input_manifest/findings.tsv")
    selected_missing = [r for r in findings if args.script == "all" or r.get("script") == args.script]
    summary = {
        "script_selection": args.script, "smoke_test": args.max_inputs is not None,
        "source_paths_measured": len(inputs), "observations": len(results),
        "saved_historical_checks": sum(bool(json.loads(r["expected_json"])) for r in observations),
        "mismatches": sum(r["historical_check"] != "passed" for r in results),
        "missing_historical_measurements": len(selected_missing),
        "missing_historical_measurements_all_scripts": len(findings),
        "note": "Per-file and sampled statistics are repeated; narrative interpretations and legacy approximate pooled statistics are not recomputed.",
    }
    (args.output_dir / "validation.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    if summary["mismatches"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
