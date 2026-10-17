"""Portable entry point for the archived alignment-cutting and clustering workflow."""
import argparse
import csv
import hashlib
import json
import logging
import math
import multiprocessing
import os
import shutil
import subprocess
import sys
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load_workflow(workdir, threads=19, workers=2, usearch=None):
    workdir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(workdir / ".matplotlib"))
    os.chdir(workdir)
    sys.path.insert(0, str(ROOT / "code"))
    # The original workflow takes 80% of available cores. Limit its pools explicitly.
    multiprocessing.cpu_count = lambda: math.ceil(workers / 0.8)
    from pcr_simulation import silva_algn_cut_and_count as workflow
    from pcr_simulation import usearch as usearch_module
    from pcr_simulation import seqkit as seqkit_module
    usearch_module.USEARCH_11_BIN = str(usearch or ROOT / "tools/usearch11.0.667_i86linux64")
    usearch_module.POOL_SIZE = threads
    workflow.POOL_SIZE = threads
    for handler in logging.getLogger().handlers:
        handler.setLevel(logging.INFO)
    logging.getLogger("matplotlib").setLevel(logging.WARNING)

    def strict_subprocess(command, *args, **kwargs):
        kwargs["check"] = True
        try:
            return subprocess.run(command, *args, **kwargs)
        except subprocess.CalledProcessError as error:
            raise RuntimeError(f"External command failed: {command}") from error

    # Raise a normal exception on failure, so a failed worker cannot silently disappear.
    usearch_module.my_subprocess_run = strict_subprocess
    seqkit_module.my_subprocess_run = strict_subprocess
    if not hasattr(workflow, "_bundle_original_common"):
        workflow._bundle_original_common = workflow.seqkit_common

    def common_with_ready_index(*args, **kwargs):
        result = workflow._bundle_original_common(*args, **kwargs)
        if result[0] == 0:
            # Build the shared index before workers can race while creating it.
            subprocess.run(["seqkit", "subseq", "-r", "1:1", "-o", os.devnull, str(result[1])],
                           check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        return result

    workflow.seqkit_common = common_with_ready_index
    return workflow


def compare_tables(expected, observed):
    with expected.open() as handle:
        expected_rows = list(csv.DictReader(handle, delimiter="\t"))
    with observed.open() as handle:
        observed_rows = list(csv.DictReader(handle, delimiter="\t"))
    expected_by_name = {r["FASTA File"]: r for r in expected_rows}
    observed_by_name = {r["FASTA File"]: r for r in observed_rows}
    differences = []
    for name in sorted(expected_by_name.keys() | observed_by_name.keys()):
        if name not in expected_by_name or name not in observed_by_name:
            differences.append({"row": name, "column": "row_presence", "expected": name in expected_by_name,
                                "observed": name in observed_by_name})
            continue
        for field, value in expected_by_name[name].items():
            if field == "FASTA File":
                continue
            actual = observed_by_name[name].get(field, "")
            try:
                equal = Decimal(value.replace(",", "")) == Decimal(actual.replace(",", ""))
            except Exception:
                equal = value == actual
            if not equal:
                differences.append({"row": name, "column": field, "expected": value, "observed": actual})
    return {"passed": not differences, "expected_rows": len(expected_rows), "observed_rows": len(observed_rows),
            "differences": differences}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["check", "plots", "reproduce"])
    parser.add_argument("--datasets", nargs="+", choices=["silva", "ltp", "gg2"], default=["silva", "ltp", "gg2"])
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data")
    parser.add_argument("--outdir", type=Path, default=ROOT / "results/reproduced_run")
    parser.add_argument("--threads", type=int,
                        help="USEARCH thread count; historical runs used 19 and changing it can change clusters")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--usearch", type=Path)
    parser.add_argument("--smoke-records", type=int, default=0,
                        help="Process only the first N input records; a smoke test, not a historical rerun")
    args = parser.parse_args()
    if args.workers < 1 or (args.threads is not None and args.threads < 1) or args.smoke_records < 0:
        parser.error("Threads/workers must be positive and smoke-records must be nonnegative")
    args.data_dir = args.data_dir.resolve()
    args.outdir = args.outdir.resolve()
    if args.usearch:
        args.usearch = args.usearch.resolve()
    configs = json.loads((ROOT / "config/runs.json").read_text())
    threads_by_dataset = {key: args.threads if args.threads is not None else configs[key]["usearch_threads"]
                          for key in args.datasets}
    workflow = load_workflow(args.outdir, threads_by_dataset[args.datasets[0]], args.workers, args.usearch)
    import matplotlib
    import numpy
    binary = args.usearch or ROOT / "tools/usearch11.0.667_i86linux64"
    checks = {"python": sys.version, "numpy": numpy.__version__, "matplotlib": matplotlib.__version__,
              "seqkit": shutil.which("seqkit"), "usearch_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
              "usearch_threads_by_dataset": threads_by_dataset, "workers": args.workers}
    if not checks["seqkit"]:
        raise RuntimeError("seqkit is missing; activate the supplied Conda environment")
    checks["seqkit_version"] = subprocess.check_output(["seqkit", "version"], text=True).strip()
    checks["usearch_version"] = subprocess.check_output([str(binary), "-version"], text=True).strip()
    (args.outdir / "software.json").write_text(json.dumps(checks, indent=2) + "\n")
    if args.mode == "check":
        print(json.dumps(checks, indent=2))
        return 0
    report = {"mode": args.mode, "smoke_records": args.smoke_records, "datasets": {}, "passed": True}
    for key in args.datasets:
        config = configs[key]
        from historical_taxonomy import configure
        configure(workflow, config["taxonomy_profile"])
        from pcr_simulation import usearch as usearch_module
        usearch_module.POOL_SIZE = threads_by_dataset[key]
        workflow.POOL_SIZE = threads_by_dataset[key]
        archived = ROOT / "results/archive_run" / key
        target = args.outdir / key
        target.mkdir(parents=True, exist_ok=True)
        if args.mode == "reproduce" and list(target.glob("run_*")):
            raise RuntimeError(f"Existing run under {target}; select a fresh --outdir")
        if args.mode == "plots":
            figure = workflow.plot_from_stats_tsv(archived / "silva_amplicon_overall_stats.tsv", target)
            if not figure or not figure.with_suffix(".svg").exists():
                raise RuntimeError("Failed to regenerate feature-length figure for " + key)
            summary = archived / "gg2_top25_species_label_frequency_with_clusters.tsv"
            if summary.exists():
                from plot_species_summary import plot
                plot(summary, target / "top25_species", top_n=25)
                plot(summary, target / "top5_species", top_n=5, vertical=True)
            report["datasets"][key] = {"passed": True, "source": "archived numerical summaries"}
            continue
        input_path = args.data_dir / config["input"]
        if not input_path.is_file():
            raise FileNotFoundError(f"Missing {input_path}; run scripts/prepare_inputs.py --extract first")
        if args.smoke_records:
            subset = target / "smoke_input.fasta"
            selected = 0
            with subset.open("w") as handle, (target / "smoke_selected_headers.tsv").open("w") as headers:
                headers.write("sequence_id\theader\n")
                for header, sequence in workflow._iter_fasta_records(input_path):
                    # Ensure the shortest window has bases, even in small smoke fixtures.
                    start = config["start_pos"] - 1
                    trimmed = sequence[start:config["end_pos"]].replace("-", "").replace(".", "")
                    shortest = sequence[start:start + config["step_size"]].replace("-", "").replace(".", "")
                    if not workflow._validate_species_taxonomy(header)[0]:
                        continue
                    if not config["min_amplicon_length"] <= len(trimmed) <= config["max_amplicon_length"]:
                        continue
                    if len(shortest) < 32:
                        continue
                    handle.write(f">{header}\n{sequence}\n")
                    headers.write(header.split()[0] + "\t" + header + "\n")
                    selected += 1
                    if selected == args.smoke_records:
                        break
            if selected != args.smoke_records:
                raise RuntimeError(f"Only {selected} eligible records found for the smoke test")
            input_path = subset
        workflow.generate_amplicon_fastas(
            silva_algn_fasta=str(input_path), start_pos=config["start_pos"], end_pos=config["end_pos"],
            step_size=config["step_size"], min_amplicon_length=config["min_amplicon_length"],
            max_amplicon_length=config["max_amplicon_length"], output_dir=str(target / "run_"),
            clust_ids=config["clust_ids"], reuse_intermediates=False,
        )
        candidates = list(target.glob("run_*/start*/silva_amplicon_overall_stats.tsv"))
        if len(candidates) != 1:
            raise RuntimeError(f"Expected exactly one fresh output table for {key}; use an empty output directory")
        stats = candidates[0]
        if args.smoke_records:
            with stats.open() as handle:
                observed = list(csv.DictReader(handle, delimiter="\t"))
            expected_rows = 2 + len(range(config["step_size"], config["end_pos"] - config["start_pos"] + 1, config["step_size"]))
            passed = len(observed) == expected_rows and all(int(r["Number of Sequences"]) > 0 for r in observed)
            entry = {"passed": passed, "observed_rows": len(observed), "expected_rows": expected_rows,
                     "input_records": args.smoke_records,
                     "input_subset_sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
                     "source": "first eligible input records; all window lengths and identity thresholds"}
        else:
            entry = compare_tables(archived / "silva_amplicon_overall_stats.tsv", stats)
        entry["output_table"] = str(stats.relative_to(args.outdir))
        entry["taxonomy_profile"] = config["taxonomy_profile"]
        entry["usearch_threads"] = threads_by_dataset[key]
        report["datasets"][key] = entry
        report["passed"] = report["passed"] and entry["passed"]
        (args.outdir / "validation.json").write_text(json.dumps(report, indent=2) + "\n")
    (args.outdir / "validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
