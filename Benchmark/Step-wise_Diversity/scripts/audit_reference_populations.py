"""Check full-database and length-filtered populations against archived counts."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
from run_stepwise import ROOT, load_workflow
from historical_taxonomy import species_name


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    data = args.data_dir.resolve()
    output = args.outdir.resolve()
    workflow = load_workflow(output)
    configs = json.loads((ROOT / "config/runs.json").read_text())
    report = {"datasets": {}, "passed": True,
              "scope": "All input records; taxonomy/length filters and species frequencies. Clustering not rerun here."}
    for key, config in configs.items():
        full_taxa = set()
        filtered_taxa = set()
        full_count = filtered_count = 0
        frequencies = Counter()
        for header, sequence in workflow._iter_fasta_records(data / config["input"]):
            strict, _ = workflow._strict_taxonomy_from_header(header)
            if not strict:
                continue
            taxonomy = (species_name(header) if config["taxonomy_profile"] == "legacy_species_name"
                        else workflow._strict_taxonomy_key_from_header(header.split()[0] + " " + strict))
            full_count += 1
            full_taxa.add(taxonomy)
            full_key = workflow._strict_taxonomy_key_from_header(header)
            if "bacteria" in full_key:
                frequencies[full_key] += 1
            trimmed = sequence[config["start_pos"] - 1:config["end_pos"]].replace("-", "").replace(".", "")
            if config["min_amplicon_length"] <= len(trimmed) <= config["max_amplicon_length"]:
                filtered_count += 1
                filtered_taxa.add(taxonomy)
        with (ROOT / "results/archive_run" / key / "silva_amplicon_overall_stats.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        full = next(r for r in rows if r["FASTA File"].startswith("silva_full_length"))
        filtered = next(r for r in rows if r["FASTA File"].startswith("silva_v1-v9_alignment"))
        checks = {
            "full_reference_records": (full_count, int(full["Number of Sequences"])),
            "full_reference_species_labels": (len(full_taxa), int(full["Number of SILVA Species"])),
            "filtered_v1v9_records": (filtered_count, int(filtered["Number of Sequences"])),
            "filtered_v1v9_species_labels": (len(filtered_taxa), int(filtered["Number of SILVA Species"])),
        }
        if key == "gg2":
            summary = ROOT / "results/archive_run/gg2/gg2_top25_species_label_frequency_with_clusters.tsv"
            with summary.open() as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    if not row["taxonomy_label"].startswith("__"):
                        checks["frequency:" + row["taxonomy_label"]] = (frequencies[row["taxonomy_label"]], int(row["frequency"]))
        entry = {"passed": all(actual == expected for actual, expected in checks.values()),
                 "checks": {name: {"observed": actual, "expected": expected, "passed": actual == expected}
                            for name, (actual, expected) in checks.items()}}
        report["datasets"][key] = entry
        report["passed"] = report["passed"] and entry["passed"]
        (output / "reference_population_validation.json").write_text(json.dumps(report, indent=2) + "\n")
        print(f"{key}: {entry['passed']}; {full_count:,} reference records, {filtered_count:,} filtered records", flush=True)
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
