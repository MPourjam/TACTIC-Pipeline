"""Recompute the GG2 top-25/top-5 species measurements from the archived reference."""
import argparse
import csv
import json
from pathlib import Path
from run_stepwise import ROOT, load_workflow
from plot_species_summary import plot


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--threads", type=int)
    args = parser.parse_args()
    data = args.data_dir.resolve()
    output = args.outdir.resolve()
    config = json.loads((ROOT / "config/runs.json").read_text())["gg2"]
    threads = args.threads if args.threads is not None else config["usearch_threads"]
    if threads < 1:
        parser.error("--threads must be positive")
    workflow = load_workflow(output, threads=threads)
    counts = workflow._count_full_taxonomy_labels_from_header_list(data / config["header_list"])
    top25 = counts.most_common(25)
    selected = {label: set() for label, _ in top25}
    if len(selected) != 25:
        raise RuntimeError("Expected 25 species labels in the archived GG2 header list")
    for header, sequence in workflow._iter_fasta_records(data / config["input"]):
        taxonomy = workflow._strict_taxonomy_key_from_header(header)
        if taxonomy not in selected:
            continue
        trimmed = sequence[config["start_pos"] - 1:config["end_pos"]].replace("-", "").replace(".", "")
        if config["min_amplicon_length"] <= len(trimmed) <= config["max_amplicon_length"]:
            selected[taxonomy].add(trimmed)
    archive = ROOT / "results/archive_run/gg2"
    with (archive / "gg2_top25_species_label_frequency_with_clusters.tsv").open() as handle:
        fields = next(csv.reader(handle, delimiter="\t"))
    cluster_ids = [float(field.removeprefix("Number of Clusters at ")) for field in fields if field.startswith("Number of Clusters at ")]
    rows = {}
    for index, (label, frequency) in enumerate(top25, 1):
        clusters = workflow._cluster_unique_sequences_for_taxonomy(
            sorted(selected[label]), cluster_ids, output / "individual", index, False)
        rows[label] = {"taxonomy_label": label, "frequency": frequency,
                       "unique_v1v9_sequences": len(selected[label]), **clusters}
    report = {"passed": True, "source": config["input"], "cluster_ids": cluster_ids,
              "usearch_threads": threads, "comparisons": {}}
    for top_n, baseline in [(25, archive / "gg2_top25_species_label_frequency_with_clusters.tsv"),
                            (5, archive / "additional_figures/gg2_top5_species_label_frequency_with_clusters.tsv")]:
        subset = top25[:top_n]
        pooled = set().union(*(selected[label] for label, _ in subset))
        clusters = workflow._cluster_unique_sequences_for_taxonomy(
            sorted(pooled), cluster_ids, output / f"combined_{top_n}", 1, False)
        combined = {"taxonomy_label": f"__top_{top_n}_combined__", "frequency": sum(freq for _, freq in subset),
                    "unique_v1v9_sequences": len(pooled), **clusters}
        observed = [combined, *[rows[label] for label, _ in subset]]
        table = output / f"top{top_n}_species_summary.tsv"
        with table.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(observed)
        with baseline.open() as handle:
            expected = list(csv.DictReader(handle, delimiter="\t"))
        expected_by_label = {r["taxonomy_label"]: r for r in expected}
        actual_by_label = {r["taxonomy_label"]: r for r in observed}
        differences = []
        for label in expected_by_label.keys() | actual_by_label.keys():
            if label not in expected_by_label or label not in actual_by_label:
                differences.append({"label": label, "field": "label_presence"})
                continue
            for field in fields[1:]:
                if int(expected_by_label[label][field]) != actual_by_label[label][field]:
                    differences.append({"label": label, "field": field,
                                        "expected": expected_by_label[label][field], "observed": actual_by_label[label][field]})
        result = {"passed": not differences, "rows": len(observed), "differences": differences}
        report["comparisons"][f"top{top_n}"] = result
        report["passed"] = report["passed"] and result["passed"]
        plot(table, output / f"top{top_n}_species", top_n, vertical=(top_n == 5))
        (output / "species_summary_validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
