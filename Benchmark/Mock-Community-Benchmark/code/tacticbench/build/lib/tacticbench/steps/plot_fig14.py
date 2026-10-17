from __future__ import annotations

from fnmatch import fnmatch
from itertools import combinations
from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..log import get_logger
from ..utils.blast import filter_hits, read_hits
from ..utils.io import genus_species
from ..utils.labels import display_community_id, display_methods
from ..utils.paths import expected_dir, figures_root, method_dir, qc_dir
from . import plot_fig09

METHOD_ORDER = ["ASV", "zOTU", "OTU", "TIC", "TAC"]
CLASS_ORDER = [
    "v4_low_resolution_specific",
    "high_similarity_full16s",
    "not_explained_by_v4_threshold",
    "missing_similarity_data",
]
CLASS_LABELS = {
    "v4_low_resolution_specific": "V4-limited resolution (V4>=t, full16S<t)",
    "high_similarity_full16s": "High full-length 16S similarity (V4>=t, full16S>=t)",
    "not_explained_by_v4_threshold": "Not supported by V4 similarity (V4<t)",
    "missing_similarity_data": "Unclassifiable (missing similarity data)",
}
CLASS_COLORS = {
    # Color-blind-safe (Okabe-Ito style) palette
    "v4_low_resolution_specific": "#E69F00",   # orange
    "high_similarity_full16s": "#0072B2",      # blue
    "not_explained_by_v4_threshold": "#CC79A7",# reddish purple
    "missing_similarity_data": "#7F7F7F",      # neutral gray
}
_METHOD_ALIASES = {
    "asv": "asv",
    "dada2": "asv",
    "zotu": "zotu",
    "unoise3": "zotu",
    "otu": "otu",
    "uparse": "otu",
    "tic": "tic",
    "tac": "tac",
    "default": "default",
}


def _as_percent_identity(x: float | int | str) -> float:
    v = float(x)
    if 0.0 <= v <= 1.0:
        return v * 100.0
    return v


def _normalize_method_token(raw: str) -> str:
    key = str(raw).strip().lower()
    return _METHOD_ALIASES.get(key, key)


def _fig14_method_min_effective_pident(cfg, method: str) -> float:
    if method in {"ASV", "zOTU"}:
        base = _as_percent_identity(getattr(cfg.figures, "fig14_graph_min_effective_pident_denoising", 100.0))
    elif method in {"TIC", "TAC"}:
        base = _as_percent_identity(getattr(cfg.figures, "fig14_graph_min_effective_pident_tic_tac", 98.7))
    else:
        base = _as_percent_identity(getattr(cfg.figures, "fig14_graph_min_effective_pident_otu", 97.0))

    by_method = getattr(cfg.figures, "fig14_edge_threshold_by_method", {}) or {}
    if not isinstance(by_method, dict) or (not by_method):
        return float(base)

    mkey = _normalize_method_token(method)
    picked = None
    for rk, rv in by_method.items():
        if _normalize_method_token(str(rk)) == mkey:
            picked = rv
            break
    if picked is None and "default" in by_method:
        picked = by_method.get("default")
    if picked is None:
        return float(base)
    try:
        return _as_percent_identity(float(picked))
    except Exception:
        return float(base)


def _fig14_min_query_coverage(cfg, method: str) -> float:
    qcov = getattr(cfg.figures, "fig14_min_query_coverage", None)
    if qcov is None:
        qcov = getattr(cfg.figures, "fig09_min_query_coverage", 0.0)
    try:
        base = _as_percent_identity(float(qcov or 0.0))
    except Exception:
        base = 0.0

    by_method = getattr(cfg.figures, "fig14_min_query_coverage_by_method", {}) or {}
    if not isinstance(by_method, dict) or (not by_method):
        return float(base)

    mkey = _normalize_method_token(method)
    picked = None
    for rk, rv in by_method.items():
        if _normalize_method_token(str(rk)) == mkey:
            picked = rv
            break
    if picked is None and "default" in by_method:
        picked = by_method.get("default")
    if picked is None:
        return float(base)
    try:
        return _as_percent_identity(float(picked))
    except Exception:
        return float(base)


def _is_excluded(cid: str, patterns: list[str]) -> bool:
    for pat in patterns:
        p = str(pat).strip()
        if not p:
            continue
        if fnmatch(cid, p) or (p in cid):
            return True
    return False


def _resolve_mock_source_dir(cfg, community_id: str) -> Path:
    # cfg.paths.root is expected to point to the directory that directly contains
    # community folders (e.g., .../All_Analyses/mocks).
    root = Path(cfg.paths.root)
    direct = root / str(community_id) / "source"
    if direct.exists():
        return direct
    # Backward-compatible fallback for layouts where root is one level above mocks/.
    nested = root / "mocks" / str(community_id) / "source"
    return nested


def _load_similarity_matrices(
    paths: list[Path],
    label_transform,
) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for p in paths:
        try:
            m = pd.read_csv(p, sep="\t", index_col=0)
        except Exception:
            continue
        if m.empty:
            continue

        m.index = [label_transform(str(x)) for x in m.index.astype(str)]
        m.columns = [label_transform(str(x)) for x in m.columns.astype(str)]
        m = m.apply(pd.to_numeric, errors="coerce")
        m = m.groupby(level=0, sort=False).max()
        m = m.T.groupby(level=0, sort=False).max().T
        m = m.clip(lower=0.0, upper=1.0)
        frames.append(m)

    if not frames:
        return pd.DataFrame()

    if len(frames) == 1:
        return frames[0]

    labels = sorted(
        set().union(
            *[
                set(frame.index.astype(str).tolist()) | set(frame.columns.astype(str).tolist())
                for frame in frames
            ]
        )
    )
    merged_arr = np.full((len(labels), len(labels)), np.nan, dtype=float)
    for frame in frames:
        arr = frame.reindex(index=labels, columns=labels).to_numpy(dtype=float)
        merged_arr = np.fmax(merged_arr, arr)
    return pd.DataFrame(merged_arr, index=labels, columns=labels).clip(lower=0.0, upper=1.0)


def _load_species_similarity_matrix(source_dir: Path, region: str) -> pd.DataFrame:
    paths = sorted(source_dir.glob(f"sample_*__species_similarity_matrix__{region}.tsv"))
    if not paths:
        return pd.DataFrame()
    return _load_similarity_matrices(paths, label_transform=genus_species)


def _canonical_tip_id(tip_id: str) -> str:
    return re.sub(r"__Copy\d+$", "", str(tip_id))


def _load_tip_similarity_matrix(source_dir: Path, region: str) -> pd.DataFrame:
    paths = sorted(source_dir.glob(f"sample_*__tip_similarity_matrix__{region}.tsv"))
    if not paths:
        return pd.DataFrame()
    return _load_similarity_matrices(paths, label_transform=_canonical_tip_id)


def _build_species_pair_similarity_from_tips(
    tip_matrix: pd.DataFrame,
    tip_to_species: dict[str, str],
    fallback_species_matrix: pd.DataFrame,
) -> pd.DataFrame:
    # Use max cross-tip similarity between species pairs to match how mapping
    # edges are established at expected-tip level. Falls back to the species
    # matrix where tip-based values are unavailable.
    if tip_matrix.empty:
        return fallback_species_matrix.copy() if not fallback_species_matrix.empty else pd.DataFrame()

    m = tip_matrix.copy()
    labels = sorted(set(m.index.astype(str).tolist()) | set(m.columns.astype(str).tolist()))
    species_to_tips: dict[str, list[str]] = {}
    for lbl in labels:
        sp = tip_to_species.get(lbl)
        if sp is None:
            sp = tip_to_species.get(_canonical_tip_id(lbl))
        if sp is None:
            sp = genus_species(lbl)
        sp = genus_species(str(sp))
        species_to_tips.setdefault(sp, []).append(lbl)

    species = sorted(species_to_tips.keys())
    out = pd.DataFrame(index=species, columns=species, dtype=float)
    for sp_a in species:
        tips_a = [t for t in species_to_tips.get(sp_a, []) if t in m.index]
        if not tips_a:
            continue
        for sp_b in species:
            tips_b = [t for t in species_to_tips.get(sp_b, []) if t in m.columns]
            if not tips_b:
                continue
            sub = m.loc[tips_a, tips_b]
            arr = sub.to_numpy(dtype=float)
            if arr.size == 0:
                continue
            with np.errstate(invalid="ignore"):
                v = np.nanmax(arr)
            if np.isfinite(v):
                out.at[sp_a, sp_b] = float(v)

    out = out.clip(lower=0.0, upper=1.0)
    if fallback_species_matrix.empty:
        return out

    fb = fallback_species_matrix.copy()
    all_species = sorted(set(out.index.astype(str)) | set(out.columns.astype(str)) | set(fb.index.astype(str)) | set(fb.columns.astype(str)))
    merged = pd.DataFrame(index=all_species, columns=all_species, dtype=float)
    for sp_a in all_species:
        for sp_b in all_species:
            vals: list[float] = []
            if sp_a in out.index and sp_b in out.columns:
                v = out.at[sp_a, sp_b]
                if pd.notna(v):
                    vals.append(float(v))
            if sp_a in fb.index and sp_b in fb.columns:
                v = fb.at[sp_a, sp_b]
                if pd.notna(v):
                    vals.append(float(v))
            if vals:
                merged.at[sp_a, sp_b] = max(vals)
    return merged.clip(lower=0.0, upper=1.0)


def _similarity_lookup(m: pd.DataFrame, species_a: str, species_b: str) -> float:
    if m.empty:
        return float("nan")
    a = genus_species(str(species_a))
    b = genus_species(str(species_b))
    try:
        if (a in m.index) and (b in m.columns):
            v = m.at[a, b]
            return float(v) if pd.notna(v) else float("nan")
        if (b in m.index) and (a in m.columns):
            v = m.at[b, a]
            return float(v) if pd.notna(v) else float("nan")
    except Exception:
        return float("nan")
    return float("nan")


def _classify_pair(v4_sim: float, full_sim: float, threshold_fraction: float) -> str:
    if not np.isfinite(v4_sim):
        return "missing_similarity_data"
    if v4_sim < threshold_fraction:
        return "not_explained_by_v4_threshold"
    if not np.isfinite(full_sim):
        return "missing_similarity_data"
    if full_sim < threshold_fraction:
        return "v4_low_resolution_specific"
    return "high_similarity_full16s"


def _build_merge_events_for_method(
    cfg,
    outdir: Path,
    community_id: str,
    method: str,
    tip_to_species: dict[str, str],
    v4_matrix: pd.DataFrame,
    full_matrix: pd.DataFrame,
) -> list[dict[str, object]]:
    md = method_dir(outdir, community_id, method)
    tab_p = md / "table_std_tss1000.tsv"
    hits_p = md / "blast_hits.tsv"
    if (not tab_p.exists()) or (not hits_p.exists()) or hits_p.stat().st_size == 0:
        return []

    try:
        tab = pd.read_csv(tab_p, sep="\t", index_col=0)
    except Exception:
        return []
    if tab.empty:
        return []
    tab = tab.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    tab.index = tab.index.astype(str)

    try:
        hits_all = read_hits(hits_p)
    except Exception:
        return []
    if hits_all.empty:
        return []

    min_ident = _fig14_method_min_effective_pident(cfg, method)
    min_qcov = _fig14_min_query_coverage(cfg, method)
    hits = filter_hits(
        hits_all,
        min_effective_pident=min_ident,
        min_qcovhsp=min_qcov,
        strict_pident_gt=False,
    )
    if hits.empty:
        return []

    hits = hits.copy()
    hits["qseqid"] = hits["qseqid"].astype(str)
    hits["sseqid"] = hits["sseqid"].astype(str)
    hits["target_species"] = hits["sseqid"].map(lambda x: tip_to_species.get(str(x), genus_species(str(x))))
    hit_pairs = hits[["qseqid", "target_species"]].dropna().drop_duplicates()
    feature_to_species: dict[str, set[str]] = {}
    for q, g in hit_pairs.groupby("qseqid", sort=False):
        feature_to_species[str(q)] = set(g["target_species"].astype(str).tolist())

    threshold_fraction = float(min_ident) / 100.0
    out: list[dict[str, object]] = []
    for sample_id in tab.columns:
        present = set(tab.index[tab[sample_id].astype(float) > 0].astype(str).tolist())
        for feat in present:
            species_set = sorted(feature_to_species.get(feat, set()))
            k = len(species_set)
            if k <= 1:
                continue
            merge_weight = float(k - 1)
            pair_list = list(combinations(species_set, 2))
            if not pair_list:
                continue
            pair_weight = merge_weight / float(len(pair_list))

            for sp_a, sp_b in pair_list:
                v4_sim = _similarity_lookup(v4_matrix, sp_a, sp_b)
                full_sim = _similarity_lookup(full_matrix, sp_a, sp_b)
                merge_class = _classify_pair(v4_sim, full_sim, threshold_fraction)
                out.append(
                    {
                        "community_id": str(community_id),
                        "sample_id": str(sample_id),
                        "sample_key": f"{community_id}|{sample_id}",
                        "method": str(method),
                        "inferred_feature": str(feat),
                        "n_merged_species_for_feature": int(k),
                        "feature_merge_weight": float(merge_weight),
                        "pair_species_a": str(sp_a),
                        "pair_species_b": str(sp_b),
                        "pair_weight": float(pair_weight),
                        "threshold_effective_identity_percent": float(min_ident),
                        "threshold_effective_identity_fraction": float(threshold_fraction),
                        "v4_similarity": float(v4_sim) if np.isfinite(v4_sim) else np.nan,
                        "full16s_similarity": float(full_sim) if np.isfinite(full_sim) else np.nan,
                        "merge_class": str(merge_class),
                    }
                )
    return out


def _plot_left_panel(
    summary: pd.DataFrame,
    out_path: Path,
    title: str,
    font_size: float,
    show_missing_similarity_data_category: bool,
) -> None:
    tick_size = max(7.0, font_size - 2.0)
    methods = [m for m in METHOD_ORDER if m in set(summary["method"].astype(str).tolist())]
    if not methods:
        methods = METHOD_ORDER[:]
    labels = display_methods(methods)
    x = np.arange(len(methods), dtype=float)

    fig, ax = plt.subplots(1, 1, figsize=(10.2, 5.4))
    bottoms = np.zeros(len(methods), dtype=float)
    class_order_plot = list(CLASS_ORDER)
    if not show_missing_similarity_data_category:
        class_order_plot = [c for c in class_order_plot if c != "missing_similarity_data"]

    for cls in class_order_plot:
        vals = []
        for m in methods:
            sub = summary[summary["method"].astype(str) == m]
            vals.append(float(sub[f"fraction_{cls}"].iloc[0]) if not sub.empty else 0.0)
        vals = np.array(vals, dtype=float)
        ax.bar(
            x,
            vals,
            bottom=bottoms,
            color=CLASS_COLORS.get(cls, "#777777"),
            edgecolor="#333333",
            linewidth=0.6,
            width=0.78,
            label=CLASS_LABELS.get(cls, cls),
        )
        bottoms += vals

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=tick_size)
    ax.set_ylim(0.0, 1.02)
    ax.set_ylabel("Weighted Fraction of Species-Merge Pairs", fontsize=font_size)
    ax.set_xlabel("Method", fontsize=font_size)
    ax.set_title("Merge Attribution Classes", fontsize=font_size)
    ax.tick_params(axis="y", labelsize=tick_size)
    ax.grid(True, axis="y", alpha=0.25, linewidth=0.6)

    handles, labels_ = ax.get_legend_handles_labels()
    fig.tight_layout(rect=(0.02, 0.20, 0.98, 0.92))
    if handles:
        fig.legend(
            handles,
            labels_,
            loc="lower center",
            bbox_to_anchor=(0.5, 0.03),
            ncol=1,
            frameon=False,
            fontsize=tick_size,
        )
    fig.suptitle(f"{title} - Left Panel", fontsize=max(font_size + 1.0, 10.0), y=0.98)
    fig.savefig(out_path, format="svg", bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


def _plot_right_panel(summary: pd.DataFrame, out_path: Path, title: str, font_size: float) -> None:
    tick_size = max(7.0, font_size - 2.0)
    methods = [m for m in METHOD_ORDER if m in set(summary["method"].astype(str).tolist())]
    if not methods:
        methods = METHOD_ORDER[:]
    labels = display_methods(methods)
    x = np.arange(len(methods), dtype=float)

    fig, ax2 = plt.subplots(1, 1, figsize=(9.8, 5.4))
    explained_vals = []
    specific_vals = []
    for m in methods:
        sub = summary[summary["method"].astype(str) == m]
        explained_vals.append(float(sub["fraction_v4_explained_usable"].iloc[0]) if not sub.empty else 0.0)
        specific_vals.append(float(sub["fraction_v4_specific_usable"].iloc[0]) if not sub.empty else 0.0)
    explained_vals = np.array(explained_vals, dtype=float)
    specific_vals = np.array(specific_vals, dtype=float)

    w = 0.36
    ax2.bar(
        x - (w / 2.0),
        explained_vals,
        width=w,
        color="#0072B2",
        edgecolor="#333333",
        linewidth=0.6,
        label="Explained by sequence similarity (usable pairs)",
    )
    ax2.bar(
        x + (w / 2.0),
        specific_vals,
        width=w,
        color="#E69F00",
        edgecolor="#333333",
        linewidth=0.6,
        label="V4-limited resolution only (usable pairs)",
    )
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels, fontsize=tick_size)
    ax2.set_ylim(0.0, 1.02)
    ax2.set_ylabel("Fraction", fontsize=font_size)
    ax2.set_xlabel("Method", fontsize=font_size)
    ax2.set_title("Usable Merge-Pair Fractions", fontsize=font_size)
    ax2.tick_params(axis="y", labelsize=tick_size)
    ax2.grid(True, axis="y", alpha=0.25, linewidth=0.6)

    # annotate usable pair counts
    for i, m in enumerate(methods):
        sub = summary[summary["method"].astype(str) == m]
        if sub.empty:
            continue
        n_usable = float(sub["usable_pair_weight"].iloc[0])
        ax2.text(
            i,
            1.01,
            f"usable={n_usable:.1f}",
            ha="center",
            va="bottom",
            fontsize=max(6.0, tick_size - 1.0),
            color="#333333",
        )

    handles, labels_ = ax2.get_legend_handles_labels()
    fig.tight_layout(rect=(0.02, 0.20, 0.98, 0.92))
    fig.legend(
        handles,
        labels_,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.03),
        ncol=1,
        frameon=False,
        fontsize=tick_size,
    )
    fig.suptitle(f"{title} - Right Panel", fontsize=max(font_size + 1.0, 10.0), y=0.98)
    fig.savefig(out_path, format="svg", bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


def run(cfg) -> int:
    """
    Figure 14:
      Quantify whether species-level merge events are explained by V4 resolution.

    For each method/community/sample, species-level merge events are reconstructed
    from accepted BLAST hit edges using Figure-14 thresholds:
      effective_identity >= method threshold and qcov >= fig14.min_query_coverage
      (fallback to fig09.min_query_coverage only when fig14 value is unset).
    Each inferred feature that merges k species contributes (k-1) merge weight.
    That weight is distributed across all merged species pairs and each pair is
    classified using species-pair similarity derived from expected tip matrices
    (max cross-tip similarity per species pair; fallback to species matrices):
      - v4_low_resolution_specific: V4 >= threshold, full16S < threshold
      - high_similarity_full16s: V4 >= threshold, full16S >= threshold
      - not_explained_by_v4_threshold: V4 < threshold
      - missing_similarity_data: required matrix value unavailable
    """
    log = get_logger()
    outdir = Path(cfg.paths.outdir)

    mf_p = outdir / "manifest.tsv"
    if not mf_p.exists():
        log.info("Figure 14: manifest.tsv not found; run discover first.")
        return 0
    mf = pd.read_csv(mf_p, sep="\t")

    skipped = set()
    sp = qc_dir(outdir) / "skipped_communities.tsv"
    if sp.exists():
        try:
            skipped = set(pd.read_csv(sp, sep="\t")["community_id"].astype(str).tolist())
        except Exception:
            skipped = set()

    exclude_patterns = list(getattr(cfg.figures, "exclude_community_name_patterns", []) or [])
    exclude_patterns += list(getattr(cfg.figures, "exclude_community_name_patterns_fig14", []) or [])

    event_rows: list[dict[str, object]] = []
    processed_communities = 0
    included_communities: list[str] = []
    community_ids = list(dict.fromkeys(mf["community_id"].astype(str).tolist()))
    for cid in community_ids:
        if cid in skipped or _is_excluded(cid, exclude_patterns):
            continue
        included_communities.append(str(cid))

        # Load expected mappings + similarity matrices once per community.
        tip_to_species = plot_fig09._load_tip_to_species(expected_dir(outdir, cid) / "tip_to_species.tsv")
        source_dir = _resolve_mock_source_dir(cfg, cid)
        v4_species_matrix = _load_species_similarity_matrix(source_dir, region="v4")
        full_species_matrix = _load_species_similarity_matrix(source_dir, region="full")
        v4_tip_matrix = _load_tip_similarity_matrix(source_dir, region="v4")
        full_tip_matrix = _load_tip_similarity_matrix(source_dir, region="full")
        v4_matrix = _build_species_pair_similarity_from_tips(
            tip_matrix=v4_tip_matrix,
            tip_to_species=tip_to_species,
            fallback_species_matrix=v4_species_matrix,
        )
        full_matrix = _build_species_pair_similarity_from_tips(
            tip_matrix=full_tip_matrix,
            tip_to_species=tip_to_species,
            fallback_species_matrix=full_species_matrix,
        )

        for method in METHOD_ORDER:
            rows = _build_merge_events_for_method(
                cfg=cfg,
                outdir=outdir,
                community_id=cid,
                method=method,
                tip_to_species=tip_to_species,
                v4_matrix=v4_matrix,
                full_matrix=full_matrix,
            )
            event_rows.extend(rows)
        processed_communities += 1

    out_dir = figures_root(outdir) / "fig14_merge_resolution_audit"
    out_dir.mkdir(parents=True, exist_ok=True)

    event_df = pd.DataFrame(event_rows)
    title = (
        getattr(cfg.figures, "fig14_title", None)
        or "Merge Resolution Audit: V4 vs Full-Length 16S Similarity"
    )
    font_size = float(getattr(cfg.figures, "fig14_font_size", 11.0))
    show_missing = bool(getattr(cfg.figures, "fig14_show_missing_similarity_data_category", False))

    # Community-level occurrence counts:
    # one occurrence row = one (community_id, sample_id, method) with >=1 species-level merge event.
    occurrence_cols = [
        "community_id",
        "n_species_merge_occurrence_rows",
        "n_samples_with_species_merge",
        "n_methods_with_species_merge",
        "n_merge_pair_events",
        "sum_pair_weight",
    ]
    if event_df.empty:
        occ_df = pd.DataFrame({"community_id": sorted(set(included_communities))})
        occ_df["n_species_merge_occurrence_rows"] = 0
        occ_df["n_samples_with_species_merge"] = 0
        occ_df["n_methods_with_species_merge"] = 0
        occ_df["n_merge_pair_events"] = 0
        occ_df["sum_pair_weight"] = 0.0
        occ_df = occ_df[occurrence_cols]
    else:
        per_occ = (
            event_df[["community_id", "sample_id", "method"]]
            .drop_duplicates()
            .groupby("community_id", as_index=False)
            .agg(
                n_species_merge_occurrence_rows=("sample_id", "count"),
                n_samples_with_species_merge=("sample_id", "nunique"),
                n_methods_with_species_merge=("method", "nunique"),
            )
        )
        per_evt = (
            event_df.groupby("community_id", as_index=False)
            .agg(
                n_merge_pair_events=("pair_species_a", "count"),
                sum_pair_weight=("pair_weight", "sum"),
            )
        )
        occ_df = per_occ.merge(per_evt, on="community_id", how="outer")
        if included_communities:
            all_comm = pd.DataFrame({"community_id": sorted(set(included_communities))})
            occ_df = all_comm.merge(occ_df, on="community_id", how="left")
        for c in [
            "n_species_merge_occurrence_rows",
            "n_samples_with_species_merge",
            "n_methods_with_species_merge",
            "n_merge_pair_events",
        ]:
            occ_df[c] = pd.to_numeric(occ_df[c], errors="coerce").fillna(0.0).astype(int)
        occ_df["sum_pair_weight"] = pd.to_numeric(occ_df["sum_pair_weight"], errors="coerce").fillna(0.0)
        occ_df = occ_df[occurrence_cols]
    occ_df.to_csv(out_dir / "species_merge_occurrence_counts_by_community.tsv", sep="\t", index=False)

    if event_df.empty:
        pd.DataFrame(
            columns=[
                "community_id",
                "sample_id",
                "sample_key",
                "method",
                "inferred_feature",
                "n_merged_species_for_feature",
                "feature_merge_weight",
                "pair_species_a",
                "pair_species_b",
                "pair_weight",
                "threshold_effective_identity_percent",
                "threshold_effective_identity_fraction",
                "v4_similarity",
                "full16s_similarity",
                "merge_class",
            ]
        ).to_csv(out_dir / "merge_pairs_resolution_events.tsv", sep="\t", index=False)
        log.info("Figure 14: no species-level merge pairs found under current settings.")
        return 0

    event_df.to_csv(out_dir / "merge_pairs_resolution_events.tsv", sep="\t", index=False)

    # Method-level weighted summaries.
    class_weight = (
        event_df.groupby(["method", "merge_class"], as_index=False)
        .agg(pair_weight=("pair_weight", "sum"))
        .sort_values(["method", "merge_class"], kind="mergesort")
    )
    total_weight = event_df.groupby("method", as_index=False).agg(total_pair_weight=("pair_weight", "sum"))
    summary = total_weight.copy()

    for cls in CLASS_ORDER:
        sub = class_weight[class_weight["merge_class"].astype(str) == cls][["method", "pair_weight"]].rename(
            columns={"pair_weight": f"weight_{cls}"}
        )
        summary = summary.merge(sub, on="method", how="left")
        summary[f"weight_{cls}"] = pd.to_numeric(summary[f"weight_{cls}"], errors="coerce").fillna(0.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            summary[f"fraction_{cls}"] = np.where(
                summary["total_pair_weight"] > 0,
                summary[f"weight_{cls}"] / summary["total_pair_weight"],
                0.0,
            )

    summary["v4_explained_weight"] = (
        summary["weight_v4_low_resolution_specific"] + summary["weight_high_similarity_full16s"]
    )
    summary["usable_pair_weight"] = summary["total_pair_weight"] - summary["weight_missing_similarity_data"]
    with np.errstate(divide="ignore", invalid="ignore"):
        summary["fraction_v4_explained_total"] = np.where(
            summary["total_pair_weight"] > 0,
            summary["v4_explained_weight"] / summary["total_pair_weight"],
            0.0,
        )
        summary["fraction_v4_explained_usable"] = np.where(
            summary["usable_pair_weight"] > 0,
            summary["v4_explained_weight"] / summary["usable_pair_weight"],
            0.0,
        )
        summary["fraction_v4_specific_usable"] = np.where(
            summary["usable_pair_weight"] > 0,
            summary["weight_v4_low_resolution_specific"] / summary["usable_pair_weight"],
            0.0,
        )

    summary = summary.sort_values(
        by=["method"],
        key=lambda s: pd.Categorical(s, categories=METHOD_ORDER, ordered=True),
    )
    summary.to_csv(out_dir / "merge_resolution_summary_by_method.tsv", sep="\t", index=False)

    # Community-level weighted summaries.
    comm = (
        event_df.groupby(["community_id", "merge_class"], as_index=False)
        .agg(pair_weight=("pair_weight", "sum"))
        .sort_values(["community_id", "merge_class"], kind="mergesort")
    )
    comm_tot = event_df.groupby("community_id", as_index=False).agg(total_pair_weight=("pair_weight", "sum"))
    comm_sum = comm_tot.copy()
    for cls in CLASS_ORDER:
        sub = comm[comm["merge_class"].astype(str) == cls][["community_id", "pair_weight"]].rename(
            columns={"pair_weight": f"weight_{cls}"}
        )
        comm_sum = comm_sum.merge(sub, on="community_id", how="left")
        comm_sum[f"weight_{cls}"] = pd.to_numeric(comm_sum[f"weight_{cls}"], errors="coerce").fillna(0.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            comm_sum[f"fraction_{cls}"] = np.where(
                comm_sum["total_pair_weight"] > 0,
                comm_sum[f"weight_{cls}"] / comm_sum["total_pair_weight"],
                0.0,
            )
    comm_sum.to_csv(out_dir / "merge_resolution_summary_by_community.tsv", sep="\t", index=False)

    _plot_left_panel(
        summary=summary,
        out_path=out_dir / "figure14_merge_resolution_audit_left.svg",
        title=title,
        font_size=font_size,
        show_missing_similarity_data_category=show_missing,
    )
    _plot_right_panel(
        summary=summary,
        out_path=out_dir / "figure14_merge_resolution_audit_right.svg",
        title=title,
        font_size=font_size,
    )

    log.info(
        f"Wrote Figure 14 to {out_dir} "
        f"(communities={processed_communities}, merge-pair events={len(event_df)})"
    )
    return 0
