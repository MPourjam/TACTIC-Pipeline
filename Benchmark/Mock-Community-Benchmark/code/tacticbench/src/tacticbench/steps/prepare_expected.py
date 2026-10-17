from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import pandas as pd

from ..log import get_logger
from ..utils.io import (
    collapse_to_species,
    dedupe_fasta_by_sequence,
    genus_species,
    load_fasta_records,
    load_fasta,
    read_expected_taxonomy,
    tip_to_species_map,
    uniquify_fasta_ids,
    write_fasta_records,
)
from ..utils.normalize import tss
from ..utils.paths import expected_dir, qc_dir
from ..utils.sample_selection import select_and_label_samples


def run(cfg) -> int:
    """
    Prepare expected artifacts per community:

    Outputs under:
      outdir/communities/<community_id>/expected/

    Files:
      - taxonomy_collapsed_species.tsv
      - expected_species_abundance_tss1000.tsv
      - expected_tip_abundance_tss1000.tsv
      - expected_sequences.fasta (copy)
      - tip_to_species.tsv (tip -> genus_species)
      - multi_v4_species_summary.tsv (one-row community summary)
      - multi_v4_species_details.tsv (multi-copy species rows only)
      - qc/multi_v4_species_summary.tsv (one row per selected sample)
      - qc/multi_v4_species_summary_by_community.tsv (one row per community)

    Logic:
      - Optionally reduce per-community expected samples using representative
        selection (cfg.mapping.keep_one_sample_per_community):
          - if both Even and Staggered samples exist, keep one of each
          - otherwise keep one sample
        Labels are harmonized for downstream plotting/tables.
      - Collapse source/taxonomy.tsv to genus_species rows (sum).
      - TSS normalize expected species table to total=cfg.normalization.tss_total per sample.
      - Build tip->species mapping using:
          1) exact match taxonomy row label
          2) substring match (longest taxonomy label)
          3) fallback genus_species(tip_id)
      - Distribute each species' abundance equally among its expected tips,
        producing expected tip-level table, then TSS normalize to total again (1000).
      - Enforce unique FASTA IDs, then deduplicate exact duplicate sequences
        (keep first occurrence).
      - Summarize species that have multiple unique V4 sequences based on the
        source FASTA records and write per-community summary/detail TSVs.
    """
    log = get_logger()
    outdir = Path(cfg.paths.outdir)
    manifest = pd.read_csv(outdir / "manifest.tsv", sep="\t")
    qc = qc_dir(outdir)
    qc.mkdir(parents=True, exist_ok=True)
    duplicate_rows = []
    duplicate_seq_rows = []
    sample_selection_rows = []
    multi_copy_community_rows = []
    multi_copy_sample_rows = []
    expected_sequences_name = str(
        getattr(cfg.discovery, "expected_sequences_filename", "expected-sequences.fasta") or "expected-sequences.fasta"
    ).strip().lower()
    multi_copy_prefix = "multi_v4" if "v4" in expected_sequences_name else "multi_copy"

    def _multi_copy_summary(
        community_id: str,
        source_fasta: Path,
        source_taxonomy: pd.DataFrame,
        species_abundance: pd.DataFrame,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        raw_records = load_fasta_records(source_fasta)
        tip_to_sp = tip_to_species_map(source_fasta, source_taxonomy)

        species_to_sequences: dict[str, set[str]] = defaultdict(set)
        species_to_records: dict[str, list[str]] = defaultdict(list)
        for rid, seq in raw_records:
            spx = tip_to_sp.get(rid, genus_species(rid))
            species_to_sequences[spx].add(str(seq).upper())
            species_to_records[spx].append(str(rid))

        total_abundance = float(species_abundance.to_numpy().sum()) if not species_abundance.empty else 0.0
        detail_rows: list[dict] = []
        for spx in sorted(species_to_sequences):
            n_unique = len(species_to_sequences[spx])
            if n_unique <= 1:
                continue
            abundance_sum = float(species_abundance.loc[spx].sum()) if spx in species_abundance.index else 0.0
            detail_rows.append(
                {
                    "community_id": community_id,
                    "species": spx,
                    "n_v4_unique_sequences": int(n_unique),
                    "n_v4_records": int(len(species_to_records[spx])),
                    "record_ids": ";".join(species_to_records[spx]),
                    "relative_abundance_sum": abundance_sum,
                    "relative_abundance_fraction": (abundance_sum / total_abundance) if total_abundance > 0 else 0.0,
                    "relative_abundance_percent": (100.0 * abundance_sum / total_abundance) if total_abundance > 0 else 0.0,
                }
            )

        detail_df = pd.DataFrame(
            detail_rows,
            columns=[
                "community_id",
                "species",
                "n_v4_unique_sequences",
                "n_v4_records",
                "record_ids",
                "relative_abundance_sum",
                "relative_abundance_fraction",
                "relative_abundance_percent",
            ],
        )
        multi_species = detail_df["species"].astype(str).tolist() if not detail_df.empty else []
        multi_sum = float(detail_df["relative_abundance_sum"].sum()) if not detail_df.empty else 0.0
        summary_df = pd.DataFrame(
            [
                {
                    "community_id": community_id,
                    "source_sequences_file": str(source_fasta),
                    "n_species_total": int(len(species_to_sequences)),
                    "n_species_with_multiple_v4_unique_sequences": int(len(detail_df)),
                    "multi_v4_species": ";".join(multi_species),
                    "multi_v4_species_abundance_sum": multi_sum,
                    "multi_v4_species_abundance_fraction": (multi_sum / total_abundance) if total_abundance > 0 else 0.0,
                    "multi_v4_species_abundance_percent": (100.0 * multi_sum / total_abundance) if total_abundance > 0 else 0.0,
                    "has_multi_v4_species": bool(len(detail_df) > 0),
                }
            ]
        )
        return summary_df, detail_df

    for _, r in manifest.iterrows():
        cid = r["community_id"]
        tax_p = Path(str(r["expected_taxonomy"]))
        seq_p = Path(str(r["expected_sequences"]))

        if not tax_p.exists() or not seq_p.exists():
            log.info(f"[{cid}] missing expected inputs; skipping.")
            continue

        ed = expected_dir(outdir, cid)
        ed.mkdir(parents=True, exist_ok=True)

        raw_all = read_expected_taxonomy(tax_p)
        sel = select_and_label_samples(
            community_id=str(cid),
            sample_columns=[str(c) for c in raw_all.columns],
            keep_one_sample_per_community=bool(getattr(cfg.mapping, "keep_one_sample_per_community", True)),
        )
        if not sel.selected_original:
            log.info(f"[{cid}] no expected sample columns found; skipping.")
            continue

        raw = raw_all[sel.selected_original].copy()
        raw = raw.rename(columns=sel.label_by_original)
        selected_labels = [sel.label_by_original.get(s, s) for s in sel.selected_original]
        selected_desc = ", ".join([f"{s}->{sel.label_by_original.get(s, s)}" for s in sel.selected_original])
        log.info(f"[{cid}] selected expected samples: {selected_desc}")
        for s in sel.selected_original:
            sample_selection_rows.append(
                {
                    "community_id": cid,
                    "original_sample_id": s,
                    "selected_sample_id": sel.label_by_original.get(s, s),
                    "profile": sel.profile_by_original.get(s, "other"),
                    "selected": True,
                }
            )
        sp = collapse_to_species(raw)
        sp = tss(sp, total=cfg.normalization.tss_total, axis=0)

        sp.to_csv(ed / "taxonomy_collapsed_species.tsv", sep="\t")
        sp.to_csv(ed / "expected_species_abundance_tss1000.tsv", sep="\t")

        # Copy expected sequences while enforcing:
        # 1) unique IDs via __CopyN suffixes for repeated IDs
        # 2) unique sequences (exact string) by keeping first occurrence
        expected_fasta = ed / "expected_sequences.fasta"
        raw_records = load_fasta_records(seq_p)
        unique_records, assigned_ids = uniquify_fasta_ids(raw_records, copy_suffix="__Copy")
        dedup_records, collapsed_seq = dedupe_fasta_by_sequence(unique_records)
        write_fasta_records(dedup_records, expected_fasta)

        dup_original_ids = 0
        for original_id, assigned in assigned_ids.items():
            if len(assigned) <= 1:
                continue
            dup_original_ids += 1
            duplicate_rows.append(
                {
                    "community_id": cid,
                    "original_id": original_id,
                    "total_count": len(assigned),
                    "duplicate_count": len(assigned) - 1,
                    "renamed_ids": ";".join(assigned[1:]),
                }
            )
        if dup_original_ids > 0:
            log.info(f"[{cid}] expected FASTA had duplicate IDs; renamed {dup_original_ids} original IDs.")
        if collapsed_seq:
            for rec in collapsed_seq:
                duplicate_seq_rows.append(
                    {
                        "community_id": cid,
                        "dropped_id": rec["dropped_id"],
                        "kept_id": rec["kept_id"],
                        "sequence_length": rec["sequence_length"],
                    }
                )
            log.info(f"[{cid}] expected FASTA had duplicate sequences; collapsed {len(collapsed_seq)} records.")

        # Summarize species with multiple unique V4 sequences before deduplication,
        # so repeated records remain visible in the counts.
        summary_df, detail_df = _multi_copy_summary(
            community_id=str(cid),
            source_fasta=seq_p,
            source_taxonomy=raw,
            species_abundance=sp,
        )
        summary_path = ed / f"{multi_copy_prefix}_species_summary.tsv"
        details_path = ed / f"{multi_copy_prefix}_species_details.tsv"
        summary_df.to_csv(summary_path, sep="\t", index=False)
        detail_df.to_csv(details_path, sep="\t", index=False)
        multi_copy_community_rows.extend(summary_df.to_dict(orient="records"))
        log.info(
            f"[{cid}] wrote {summary_path.name} "
            f"({int(summary_df.iloc[0]['n_species_with_multiple_v4_unique_sequences'])} multi-copy species)."
        )

        multi_species = set(detail_df["species"].astype(str).tolist()) if not detail_df.empty else set()
        n_species_total = int(summary_df.iloc[0]["n_species_total"]) if not summary_df.empty else 0
        n_multi_species = int(summary_df.iloc[0]["n_species_with_multiple_v4_unique_sequences"]) if not summary_df.empty else 0
        for original_sample_id, selected_sample_id in zip(sel.selected_original, selected_labels):
            if selected_sample_id not in sp.columns:
                continue
            sample_series = sp[selected_sample_id]
            total_abundance = float(sample_series.sum())
            multi_abundance_sum = float(sample_series.loc[sample_series.index.astype(str).isin(multi_species)].sum())
            multi_copy_sample_rows.append(
                {
                    "community_id": cid,
                    "original_sample_id": original_sample_id,
                    "selected_sample_id": selected_sample_id,
                    "profile": sel.profile_by_original.get(original_sample_id, "other"),
                    "source_sequences_file": str(seq_p),
                    "n_species_total": n_species_total,
                    "n_species_with_multiple_v4_unique_sequences": n_multi_species,
                    "multi_v4_species": ";".join(sorted(multi_species)),
                    "multi_v4_species_abundance_sum": multi_abundance_sum,
                    "multi_v4_species_abundance_fraction": (multi_abundance_sum / total_abundance) if total_abundance > 0 else 0.0,
                    "multi_v4_species_abundance_percent": (100.0 * multi_abundance_sum / total_abundance) if total_abundance > 0 else 0.0,
                    "has_multi_v4_species": bool(len(multi_species) > 0),
                }
            )

        # tip -> species
        tip_to_sp = tip_to_species_map(expected_fasta, raw)
        pd.DataFrame([{"tip": t, "species": s} for t, s in tip_to_sp.items()]).to_csv(
            ed / "tip_to_species.tsv", sep="\t", index=False
        )

        # distribute species abundance across tips
        seqs = load_fasta(expected_fasta)
        tips = list(seqs.keys())

        sp_to_tips: dict[str, list[str]] = {}
        for tip in tips:
            spx = tip_to_sp.get(tip)
            if spx is None:
                continue
            sp_to_tips.setdefault(spx, []).append(tip)

        tip_df = pd.DataFrame(0.0, index=tips, columns=sp.columns)
        for spx, abund in sp.iterrows():
            tt = sp_to_tips.get(spx, [])
            if not tt:
                continue
            share = abund / float(len(tt))
            for tip in tt:
                tip_df.loc[tip, :] = share.values

        # normalize tip table to total again (ensures per-sample total = 1000)
        tip_df = tss(tip_df, total=cfg.normalization.tss_total, axis=0)
        tip_df.to_csv(ed / "expected_tip_abundance_tss1000.tsv", sep="\t")

        log.info(f"[{cid}] prepared expected.")

    dup_qc = qc / "expected_fasta_duplicate_ids.tsv"
    pd.DataFrame(
        duplicate_rows,
        columns=["community_id", "original_id", "total_count", "duplicate_count", "renamed_ids"],
    ).to_csv(dup_qc, sep="\t", index=False)
    log.info(f"Wrote {dup_qc}")

    dup_seq_qc = qc / "expected_fasta_duplicate_sequences.tsv"
    pd.DataFrame(
        duplicate_seq_rows,
        columns=["community_id", "dropped_id", "kept_id", "sequence_length"],
    ).to_csv(dup_seq_qc, sep="\t", index=False)
    log.info(f"Wrote {dup_seq_qc}")

    sample_sel_qc = qc / "expected_sample_selection.tsv"
    pd.DataFrame(
        sample_selection_rows,
        columns=["community_id", "original_sample_id", "selected_sample_id", "profile", "selected"],
    ).to_csv(sample_sel_qc, sep="\t", index=False)
    log.info(f"Wrote {sample_sel_qc}")

    multi_copy_community_qc = qc / f"{multi_copy_prefix}_species_summary_by_community.tsv"
    pd.DataFrame(multi_copy_community_rows).to_csv(multi_copy_community_qc, sep="\t", index=False)
    log.info(f"Wrote {multi_copy_community_qc}")

    multi_copy_sample_qc = qc / f"{multi_copy_prefix}_species_summary.tsv"
    pd.DataFrame(
        multi_copy_sample_rows,
        columns=[
            "community_id",
            "original_sample_id",
            "selected_sample_id",
            "profile",
            "source_sequences_file",
            "n_species_total",
            "n_species_with_multiple_v4_unique_sequences",
            "multi_v4_species",
            "multi_v4_species_abundance_sum",
            "multi_v4_species_abundance_fraction",
            "multi_v4_species_abundance_percent",
            "has_multi_v4_species",
        ],
    ).to_csv(multi_copy_sample_qc, sep="\t", index=False)
    log.info(f"Wrote {multi_copy_sample_qc}")

    return 0
