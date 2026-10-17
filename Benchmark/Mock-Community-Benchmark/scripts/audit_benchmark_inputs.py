#!/usr/bin/env python3
"""Regenerate only TACTICBench input preparation and compare archived outputs."""

import argparse
import csv
import hashlib
import sys
import tempfile
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, default=Path(__file__).resolve().parent.parent)
    args = parser.parse_args()
    bundle = args.bundle.resolve()
    sys.path.insert(0, str(bundle / "code/tacticbench/src"))
    import pandas as pd
    from tacticbench.config import load_config
    from tacticbench.steps import discover, prepare_expected, prepare_methods

    archive = bundle / "results/archive_run"
    with (bundle / "metadata/original_run_manifest.tsv").open(newline="") as handle:
        original_manifest = list(csv.DictReader(handle, delimiter="\t"))
    original_root = str(Path(original_manifest[0]["community_path"]).parent)
    rows = []
    with tempfile.TemporaryDirectory(prefix="tacticbench_input_preparation_") as temp:
        cfg = load_config(bundle / "config/run_config.yaml", overrides={
            "root": str(bundle / "data/communities"), "outdir": temp,
        })
        for step in (discover, prepare_expected, prepare_methods):
            result = step.run(cfg)
            if result:
                raise RuntimeError(f"{step.__name__} returned {result}")
        candidates = sorted(Path(temp).glob("communities/*/expected/*.tsv"))
        candidates += sorted(Path(temp).glob("communities/*/expected/expected_sequences.fasta"))
        candidates += sorted(Path(temp).glob("communities/*/methods/*/table_std_tss1000.tsv"))
        candidates += sorted(Path(temp).glob("communities/*/methods/*/seqs.fasta"))
        candidates += [Path(temp) / "qc/sample_name_mapping.tsv", Path(temp) / "qc/expected_sample_selection.tsv"]
        for generated in candidates:
            relative = generated.relative_to(temp)
            archived = archive / relative
            passed = False
            detail = ""
            if not archived.is_file():
                detail = "No corresponding archived output"
            elif generated.suffix == ".fasta":
                passed = generated.read_bytes() == archived.read_bytes()
                detail = "FASTA bytes identical" if passed else "FASTA bytes differ"
            else:
                try:
                    a = pd.read_csv(generated, sep="\t", index_col=0)
                    b = pd.read_csv(archived, sep="\t", index_col=0)
                    # Relocation changes only the provenance column's absolute prefix.
                    if "source_sequences_file" in a.columns:
                        a["source_sequences_file"] = a["source_sequences_file"].str.replace(str(bundle / "data/communities"), "{COMMUNITIES}", regex=False)
                        b["source_sequences_file"] = b["source_sequences_file"].str.replace(original_root, "{COMMUNITIES}", regex=False)
                    pd.testing.assert_frame_equal(a, b, check_exact=False, check_dtype=False, rtol=1e-10, atol=1e-10)
                    passed = True
                    detail = "Labels/order and values match; numeric tolerance 1e-10; source path prefixes normalized when present"
                except AssertionError as exc:
                    detail = str(exc).replace("\n", " ")[:1000]
            rows.append({"artifact": relative.as_posix(), "passed": passed, "comparison": detail,
                         "archived_sha256": hashlib.sha256(archived.read_bytes()).hexdigest() if archived.is_file() else "",
                         "regenerated_sha256": hashlib.sha256(generated.read_bytes()).hexdigest(),
                         "benchmark_code": "bundled code/tacticbench", "upstream_pipeline_run": False})
    destination = bundle / "data/input_manifest/benchmark_preparation_validation.tsv"
    destination.parent.mkdir(exist_ok=True)
    with destination.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    failed = [r for r in rows if not r["passed"]]
    print(f"Benchmark preparation: {len(rows) - len(failed)}/{len(rows)} artifacts match archived outputs")
    print(destination)
    if failed:
        for row in failed:
            print(row["artifact"], row["comparison"])
        raise SystemExit(1)


if __name__ == "__main__":
    main()
