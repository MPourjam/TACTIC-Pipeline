from __future__ import annotations

from fnmatch import fnmatch
from pathlib import Path
import pandas as pd

from ..log import get_logger
from ..utils.labels import (
    build_community_display_mapping,
    set_community_display_mapping,
    write_community_display_mapping_table,
)
from ..utils.paths import qc_dir


def _pick_latest(dirs: list[Path]) -> Path | None:
    """
    Pick most recently modified directory from a list.
    """
    dirs = [d for d in dirs if d.is_dir()]
    if not dirs:
        return None
    return sorted(dirs, key=lambda p: p.stat().st_mtime, reverse=True)[0]


def _is_excluded(cid: str, patterns: list[str]) -> bool:
    """
    Check whether a community id should be excluded by configured patterns.

    Match modes:
      - fnmatch glob (e.g., "mock-ZIEL_*")
      - substring fallback for plain-text patterns (e.g., "ZIEL")
    """
    for pat in patterns:
        p = str(pat).strip()
        if not p:
            continue
        if fnmatch(cid, p):
            return True
        if p in cid:
            return True
    return False


def run(cfg) -> int:
    """
    Build manifest.tsv by discovering expected inputs, DADA2 outputs, and TACTIC outputs.

    Rules:
      - Communities can be excluded by name pattern:
          cfg.discovery.exclude_community_name_patterns
      - For each community (folder containing source/taxonomy.tsv):
          - expected: source/taxonomy.tsv + source/expected-sequences.fasta
          - DADA2: dada2_results/{dada_out|dada2_out}/ASV_counts.tsv and ASVs.fasta
          - TACTIC (new): trimmed_fastq/Analysis_*/{ZOTU,OTU,TIC,TAC}_Pipeline/*
              - choose latest Analysis_* by mtime
          - TACTIC (legacy fallback): tactic_results/Analysis*am-<MODE>* (latest mtime)
              - exclude any directory containing 'SampleMinCount_Normalized'
          - de-novo: ZOTUs-Table.tab, ZOTUs-Seqs.*, OTUs-Table.tab, OTUs-Seqs.*
          - TIC:    SOTUs-Table-TIC.tab, SOTUs-Seqs-TIC.*
          - TAC:    OTUs-Table-TAC.tab, OTUs-Seqs-TAC.*
    """
    log = get_logger()
    root = Path(cfg.paths.root)
    outdir = Path(cfg.paths.outdir)
    exclude_patterns = list(getattr(cfg.discovery, "exclude_community_name_patterns", []) or [])
    outdir.mkdir(parents=True, exist_ok=True)
    qc_dir(outdir).mkdir(parents=True, exist_ok=True)

    expected_seq_name = str(getattr(cfg.discovery, "expected_sequences_filename", "expected-sequences.fasta") or "expected-sequences.fasta").strip()
    fallback_to_full = bool(getattr(cfg.discovery, "expected_sequences_fallback_to_full", True))

    rows = []
    for tax_tsv in root.rglob("source/taxonomy.tsv"):
        comm_dir = tax_tsv.parent.parent
        cid = comm_dir.name
        if _is_excluded(cid, exclude_patterns):
            log.info(f"[{cid}] excluded by discovery.exclude_community_name_patterns")
            continue

        expected_fasta = comm_dir / "source" / expected_seq_name
        if (not expected_fasta.exists()) and fallback_to_full:
            fallback_fasta = comm_dir / "source/expected-sequences.fasta"
            if fallback_fasta.exists():
                expected_fasta = fallback_fasta
        tactic_root = comm_dir / "tactic_results"
        trimmed_fastq_root = comm_dir / "trimmed_fastq"

        # DADA2 outputs might be in dada_out or dada2_out
        asv_table = ""
        asv_fasta = ""
        for sub in ("dada_out", "dada2_out"):
            t = comm_dir / "dada2_results" / sub / "ASV_counts.tsv"
            f = comm_dir / "dada2_results" / sub / "ASVs.fasta"
            if t.exists():
                asv_table = str(t)
            if f.exists():
                asv_fasta = str(f)

        # pick latest Analysis dirs for each mode (excluding SampleMinCount_Normalized)
        def _mode_dir(pattern: str) -> Path | None:
            if not tactic_root.exists():
                return None
            hits = [
                p
                for p in tactic_root.glob(pattern)
                if p.is_dir() and "SampleMinCount_Normalized" not in str(p)
            ]
            return _pick_latest(hits)

        denovo_dir = _mode_dir("Analysis*am-de-novo*")
        tic_dir = _mode_dir("Analysis*am-TIC*")
        tac_dir = _mode_dir("Analysis*am-TAC*")

        # New TACTIC output layout:
        #   trimmed_fastq/Analysis_*/<MODE>_Pipeline/<files>
        # Keep legacy support above for older runs.
        def _analysis_dirs(base: Path) -> list[Path]:
            if not base.exists():
                return []
            return [
                p
                for p in base.glob("Analysis_*")
                if p.is_dir() and "SampleMinCount_Normalized" not in str(p)
            ]

        latest_analysis = _pick_latest(_analysis_dirs(trimmed_fastq_root))

        def _first_match(d: Path | None, pat: str) -> str:
            if d is None:
                return ""
            hits = [h for h in d.glob(pat) if "SampleMinCount_Normalized" not in str(h)]
            return str(hits[0]) if hits else ""

        def _pipeline_file(analysis_dir: Path | None, pipeline: str, filename: str) -> str:
            if analysis_dir is None:
                return ""
            p = analysis_dir / pipeline / filename
            return str(p) if p.exists() else ""

        def _prefer(*vals: str) -> str:
            for v in vals:
                if str(v).strip():
                    return v
            return ""

        rows.append(
            {
                "community_id": cid,
                "community_path": str(comm_dir),
                "expected_taxonomy": str(tax_tsv),
                "expected_sequences": str(expected_fasta) if expected_fasta.exists() else "",
                "ASV_table": asv_table,
                "ASV_seqs": asv_fasta,
                "zOTU_table": _prefer(
                    _pipeline_file(latest_analysis, "ZOTU_Pipeline", "ZOTUs-Table.tab"),
                    _first_match(denovo_dir, "ZOTUs-Table.tab"),
                ),
                "zOTU_seqs": _prefer(
                    _pipeline_file(latest_analysis, "ZOTU_Pipeline", "ZOTUs-Seqs.fasta"),
                    _first_match(denovo_dir, "ZOTUs-Seqs.*"),
                ),
                "OTU_table": _prefer(
                    _pipeline_file(latest_analysis, "OTU_Pipeline", "OTUs-Table.tab"),
                    _first_match(denovo_dir, "OTUs-Table.tab"),
                ),
                "OTU_seqs": _prefer(
                    _pipeline_file(latest_analysis, "OTU_Pipeline", "OTUs-Seqs.fasta"),
                    _first_match(denovo_dir, "OTUs-Seqs.*"),
                ),
                "TIC_table": _prefer(
                    _pipeline_file(latest_analysis, "TIC_Pipeline", "SOTUs-Table-TIC.tab"),
                    _first_match(tic_dir, "SOTUs-Table-TIC.tab"),
                ),
                "TIC_seqs": _prefer(
                    _pipeline_file(latest_analysis, "TIC_Pipeline", "SOTUs-Seqs-TIC.fasta"),
                    _first_match(tic_dir, "SOTUs-Seqs-TIC.*"),
                ),
                "TAC_table": _prefer(
                    _pipeline_file(latest_analysis, "TAC_Pipeline", "OTUs-Table-TAC.tab"),
                    _first_match(tac_dir, "OTUs-Table-TAC.tab"),
                ),
                "TAC_seqs": _prefer(
                    _pipeline_file(latest_analysis, "TAC_Pipeline", "OTUs-Seqs-TAC.fasta"),
                    _first_match(tac_dir, "OTUs-Seqs-TAC.*"),
                ),
            }
        )

    mf = pd.DataFrame(rows)
    if "community_id" in mf.columns:
        mf = mf.sort_values("community_id")
    manifest_path = outdir / "manifest.tsv"
    mf.to_csv(manifest_path, sep="\t", index=False)
    log.info(f"Wrote {manifest_path}")

    # Write publication-friendly mapping with contiguous Mock numbering after
    # applying global figure exclusions.
    fig_exclude = list(getattr(cfg.figures, "exclude_community_name_patterns", []) or [])
    cids = mf["community_id"].astype(str).tolist() if "community_id" in mf.columns else []
    mapping = build_community_display_mapping(cids, figure_exclude_patterns=fig_exclude)
    set_community_display_mapping(mapping)
    mapping_path = outdir / "community_name_mapping.tsv"
    write_community_display_mapping_table(mapping_path, mapping)
    log.info(f"Wrote {mapping_path}")

    skipped_path = qc_dir(outdir) / "skipped_communities.tsv"
    if not skipped_path.exists():
        pd.DataFrame(columns=["community_id", "reason"]).to_csv(skipped_path, sep="\t", index=False)

    return 0
