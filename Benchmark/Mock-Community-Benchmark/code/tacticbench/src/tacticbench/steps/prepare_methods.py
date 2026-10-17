from __future__ import annotations

from pathlib import Path

import pandas as pd

from ..log import get_logger
from ..utils.io import (
    ensure_fasta,
    read_dada2_asv_counts,
    read_expected_taxonomy,
    read_tactic_table,
)
from ..utils.normalize import tss
from ..utils.samplemap import map_by_substring
from ..utils.paths import method_dir, qc_dir
from ..utils.sample_selection import select_and_label_samples


METHOD_SPECS = [
    ("ASV", "ASV_table", "ASV_seqs", "dada2"),
    ("zOTU", "zOTU_table", "zOTU_seqs", "tactic"),
    ("OTU", "OTU_table", "OTU_seqs", "tactic"),
    ("TIC", "TIC_table", "TIC_seqs", "tactic"),
    ("TAC", "TAC_table", "TAC_seqs", "tactic"),
]


def run(cfg) -> int:
    """
    Standardize each method's feature table to:
      - columns = harmonized expected sample labels (mapped by substring using
        original expected sample names)
      - values = TSS normalized to cfg.normalization.tss_total per sample (default 1000)
      - optionally representative sample reduction per community:
          - if both Even and Staggered exist, keep one of each
          - otherwise keep one sample

    Writes per community + method:
      outdir/communities/<cid>/methods/<method>/
        - table_std_tss1000.tsv
        - seqs.fasta (normalized to FASTA when possible)

    Also writes QC:
      outdir/qc/sample_name_mapping.tsv
      outdir/qc/skipped_communities.tsv (if skip enabled and mapping fails)
    """
    log = get_logger()
    outdir = Path(cfg.paths.outdir)
    qc = qc_dir(outdir)
    qc.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(outdir / "manifest.tsv", sep="\t")

    skipped_path = qc / "skipped_communities.tsv"
    skipped = (
        pd.read_csv(skipped_path, sep="\t")
        if skipped_path.exists()
        else pd.DataFrame(columns=["community_id", "reason"])
    )
    # Keep non-sample-mapping skip reasons from previous runs, but recompute
    # sample_mapping_failed fresh each run to avoid stale carry-over.
    if not skipped.empty and "reason" in skipped.columns:
        skipped_static = skipped[skipped["reason"].astype(str) != "sample_mapping_failed"].copy()
    else:
        skipped_static = pd.DataFrame(columns=["community_id", "reason"])

    samplemap_rows = []
    failed_samplemap_communities = set()

    for _, r in manifest.iterrows():
        cid = r["community_id"]
        exp_tax = Path(str(r["expected_taxonomy"]))
        if not exp_tax.exists():
            continue

        expected_raw = read_expected_taxonomy(exp_tax)
        sel = select_and_label_samples(
            community_id=str(cid),
            sample_columns=[str(c) for c in expected_raw.columns],
            keep_one_sample_per_community=bool(getattr(cfg.mapping, "keep_one_sample_per_community", True)),
        )
        expected_samples = list(sel.selected_original)
        expected_label_by_original = dict(sel.label_by_original)
        if not expected_samples:
            log.info(f"[{cid}] no expected sample columns after selection; skipping.")
            continue
        selected_desc = ", ".join([f"{s}->{expected_label_by_original.get(s, s)}" for s in expected_samples])
        log.info(f"[{cid}] selected expected samples for method standardization: {selected_desc}")

        failed = False

        for method, tcol, scol, kind in METHOD_SPECS:
            table_p = Path(str(r.get(tcol, "")))
            seq_p = Path(str(r.get(scol, "")))

            if not table_p.exists():
                continue

            if kind == "dada2":
                df = read_dada2_asv_counts(table_p)  # features x inferred_samples
            else:
                df = read_tactic_table(table_p)  # features x inferred_samples

            inferred_samples = list(df.columns)
            mres = map_by_substring(expected_samples, inferred_samples)

            # record mapping QC
            for e, inf in mres.mapping.items():
                samplemap_rows.append(
                    {
                        "community_id": cid,
                        "method": method,
                        "expected_sample": e,
                        "expected_sample_label": expected_label_by_original.get(e, e),
                        "inferred_sample": inf,
                    }
                )
            for e, ties in mres.ambiguous:
                samplemap_rows.append(
                    {
                        "community_id": cid,
                        "method": method,
                        "expected_sample": e,
                        "expected_sample_label": expected_label_by_original.get(e, e),
                        "inferred_sample": mres.mapping.get(e, ""),
                        "note": "ambiguous:" + ",".join(ties),
                    }
                )

            if mres.missing:
                failed = True
                log.info(f"[{cid}] sample mapping failed for {method}; missing={mres.missing}")
                if cfg.mapping.skip_community_on_samplemap_failure:
                    break

            # build standardized table with expected sample names as columns
            out = pd.DataFrame(0.0, index=df.index, columns=expected_samples)
            for e, inf in mres.mapping.items():
                if e in out.columns and inf in df.columns:
                    out[e] = df[inf].values

            out = out.rename(columns=expected_label_by_original)

            # TSS normalize per sample to 1000
            out = tss(out, total=cfg.normalization.tss_total, axis=0)

            md = method_dir(outdir, cid, method)
            md.mkdir(parents=True, exist_ok=True)
            out.to_csv(md / "table_std_tss1000.tsv", sep="\t")

            if seq_p.exists():
                ensure_fasta(seq_p, md / "seqs.fasta")

        if failed and cfg.mapping.skip_community_on_samplemap_failure:
            failed_samplemap_communities.add(str(cid))

    pd.DataFrame(samplemap_rows).to_csv(qc / "sample_name_mapping.tsv", sep="\t", index=False)
    skipped_dynamic = pd.DataFrame(
        [{"community_id": cid, "reason": "sample_mapping_failed"} for cid in sorted(failed_samplemap_communities)]
    )
    skipped_out = pd.concat([skipped_static, skipped_dynamic], ignore_index=True)
    skipped_out.drop_duplicates(subset=["community_id", "reason"]).to_csv(skipped_path, sep="\t", index=False)

    log.info(f"Wrote {qc/'sample_name_mapping.tsv'}")
    log.info(f"Wrote {skipped_path}")
    return 0
