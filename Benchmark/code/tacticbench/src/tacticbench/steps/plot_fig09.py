from __future__ import annotations

from fnmatch import fnmatch
import hashlib
from itertools import combinations
import math
from pathlib import Path
import shutil
import subprocess
import tempfile
import textwrap

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter, MaxNLocator
import numpy as np
import pandas as pd

from ..log import get_logger
from ..utils.paths import expected_dir, figures_root, method_dir, qc_dir
from ..utils.labels import display_methods, display_sample_key
from ..utils.blast import filter_hits, read_hits
from ..utils.io import genus_species

METHOD_ORDER = ["ASV", "zOTU", "OTU", "TIC", "TAC"]
FIG09_RICHNESS_METHOD_ORDER = ["ASV", "zOTU", "TAC", "TIC", "OTU"]
CLUSTERING = {"OTU", "TIC", "TAC"}
PURITY_METHODS = {"TIC", "TAC"}
METHOD_COLORS = {
    "ASV": "#0072B2",
    "zOTU": "#56B4E9",
    "OTU": "#D55E00",
    "TIC": "#009E73",
    "TAC": "#E69F00",
}
GRAPH_LEVEL_ORDER = ("tip", "species")
FIG09_TRANSITION_CLASS_ORDER = ["Exact match", "Near-exact match", "No hit"]
FIG09_TRANSITION_OUTCOME_ORDER = ["Exact match", "Near-exact match", "No hit", "Unassigned"]


def _normalize_plot_type(raw: str, default: str = "boxplot") -> str:
    v = str(raw or "").strip().lower()
    if v in {"box", "boxplot"}:
        return "boxplot"
    if v in {"violin", "violinplot"}:
        return "violin"
    return default


def _normalize_split_index_mode(raw: str, default: str = "mapping") -> str:
    v = str(raw or "").strip().lower()
    if v in {"mapping", "threshold_hits", "graph"}:
        return v
    return default


def _normalize_split_target_level(raw: str, default: str = "species") -> str:
    v = str(raw or "").strip().lower()
    if v in {"species", "tip"}:
        return v
    return default


def _normalize_graph_target_mode(raw: str, default: str = "single") -> str:
    v = str(raw or "").strip().lower()
    aliases = {
        "single": "single",
        "default": "single",
        "per_metric": "single",
        "tip_species_split_violin": "tip_species_split_violin",
        "tip_species_dual_violin": "tip_species_split_violin",
        "dual_violin": "tip_species_split_violin",
        "split_violin_tip_species": "tip_species_split_violin",
        "tip_species": "tip_species_split_violin",
        "both_levels": "tip_species_split_violin",
        "dual": "tip_species_split_violin",
    }
    return aliases.get(v, default)


def _normalize_dual_level_layout(raw, default: str = "separate") -> str:
    if isinstance(raw, (bool, np.bool_)):
        return "combined" if bool(raw) else "separate"
    v = str(raw or "").strip().lower()
    aliases = {
        "separate": "separate",
        "split": "separate",
        "multi_file": "separate",
        "combined": "combined",
        "stacked": "combined",
        "single_file": "combined",
    }
    return aliases.get(v, default)


def _normalize_split_index_denominator(raw: str, default: str = "mapped_expected_features") -> str:
    v = str(raw or "").strip().lower()
    if v in {"mapped_expected_features", "mapped_inferred_features"}:
        return v
    return default


def _normalize_one_to_one_denominator(raw: str, default: str = "expected_present_features") -> str:
    v = str(raw or "").strip().lower()
    if v in {"expected_present_features", "mapped_inferred_features", "inferred_present_features"}:
        return v
    return default


def _normalize_merge_denominator(raw: str, default: str = "mapped_inferred_features") -> str:
    v = str(raw or "").strip().lower()
    if v in {"mapped_inferred_features", "expected_present_features"}:
        return v
    return default


def _normalize_pairwise_test(raw: str, default: str = "signed_rank_permutation") -> str:
    v = str(raw or "").strip().lower()
    if v in {"signed_rank_permutation", "sign_test_exact"}:
        return v
    return default


def _normalize_stats_model(raw: str, default: str = "mixed_effect_glmmtmb") -> str:
    v = str(raw or "").strip().lower()
    if v in {"mixed_effect_glmmtmb", "mixed_effect", "paired_permutation"}:
        return v
    return default


def _normalize_glmm_count_family(raw: str, default: str = "nbinom2") -> str:
    v = str(raw or "").strip().lower()
    if v in {"nbinom2", "nbinom1", "poisson"}:
        return v
    return default


def _normalize_emmeans_scale(raw: str, default: str = "link") -> str:
    v = str(raw or "").strip().lower()
    if v in {"link", "response"}:
        return v
    return default


def _as_bool_scalar(v) -> bool:
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    s = str(v).strip().lower()
    if s in {"true", "t", "1", "yes", "y"}:
        return True
    if s in {"false", "f", "0", "no", "n"}:
        return False
    return bool(v)


def _method_min_effective_pident(cfg, method: str) -> float:
    base = _as_percent_identity(
        float(
        getattr(
            cfg.mapping,
            "min_effective_pident",
            getattr(cfg.mapping, "min_identity", 98.7),
        ))
    )
    clust = getattr(
        cfg.mapping,
        "min_effective_pident_clustering",
        getattr(cfg.mapping, "min_identity_clustering", None),
    )
    if method in CLUSTERING and clust is not None:
        return _as_percent_identity(float(clust))
    return base


def _as_percent_identity(raw: float) -> float:
    v = float(raw)
    if 0.0 <= v <= 1.0:
        v *= 100.0
    return float(np.clip(v, 0.0, 100.0))


def _normalize_graph_metric_key(raw: str) -> str:
    k = str(raw or "").strip().lower()
    aliases = {
        "split": "split_index",
        "split_index": "split_index",
        "merge": "merge_deficit",
        "merge_index": "merge_deficit",
        "merge_deficit": "merge_deficit",
        "unmapped": "unmapped_feature_fraction",
        "unmapped_feature_fraction": "unmapped_feature_fraction",
        "one_to_one": "one_to_one_mapping_fraction",
        "one_to_one_mapping": "one_to_one_mapping_fraction",
        "one_to_one_mapping_fraction": "one_to_one_mapping_fraction",
        "one_to_one_fraction": "one_to_one_mapping_fraction",
        "cluster_purity": "cluster_purity_fraction",
        "cluster_purity_fraction": "cluster_purity_fraction",
        "purity": "cluster_purity_fraction",
        "no_hit_unmapped_fraction": "no_hit_over_unmapped_fraction",
        "no_hit_fraction": "no_hit_over_unmapped_fraction",
        "no_hit_over_unmapped": "no_hit_over_unmapped_fraction",
        "no_hit_over_unmapped_fraction": "no_hit_over_unmapped_fraction",
        "relative_richness_error": "relative_richness_error",
        "richness_error": "relative_richness_error",
        "richness": "relative_richness_error",
        "panel_d": "relative_richness_error",
    }
    return aliases.get(k, k)


def _normalize_method_token(raw: str) -> str:
    m = str(raw or "").strip().lower()
    aliases = {
        "asv": "asv",
        "dada2": "asv",
        "zotu": "zotu",
        "unoise3": "zotu",
        "otu": "otu",
        "uparse": "otu",
        "tic": "tic",
        "tac": "tac",
    }
    return aliases.get(m, m)


def _normalize_method_order(raw, default: list[str] | tuple[str, ...] = METHOD_ORDER) -> list[str]:
    if not isinstance(raw, (list, tuple)):
        raw = list(default)
    aliases = {_normalize_method_token(method): method for method in METHOD_ORDER}
    out: list[str] = []
    for item in raw:
        token = _normalize_method_token(str(item))
        method = aliases.get(token)
        if method and method not in out:
            out.append(method)
    return out or list(default)


def _normalize_panel_key(raw: str) -> str:
    k = str(raw or "").strip().lower()
    aliases = {
        "a": "a",
        "panel_a": "a",
        "species_split_merge": "a",
        "b": "b",
        "panel_b": "b",
        "tip_split_merge": "b",
        "c": "c",
        "panel_c": "c",
        "mixed_split_merge": "c",
        "d": "d",
        "panel_d": "d",
        "mixed_fractions": "d",
    }
    return aliases.get(k, k)


def _graph_target_level_for_method(cfg, method: str, metric_col: str = "split_index") -> str:
    # Figure 09 single-mode defaults:
    # denoising methods at tip level; clustering methods at species level.
    # Optional override:
    #   figures.fig09.graph_target_level_by_metric.<metric>.<method> = tip|species
    defaults = {
        "ASV": "tip",
        "zOTU": "tip",
        "TIC": "species",
        "TAC": "species",
        "OTU": "species",
    }
    d = defaults.get(method, "species" if method in CLUSTERING else "tip")

    by_metric = getattr(cfg.figures, "fig09_graph_target_level_by_metric", {}) or {}
    if isinstance(by_metric, dict) and by_metric:
        metric_key = _normalize_graph_metric_key(metric_col)
        metric_node = None
        for rk, rv in by_metric.items():
            if _normalize_graph_metric_key(str(rk)) == metric_key:
                metric_node = rv
                break
        if isinstance(metric_node, dict):
            method_key = _normalize_method_token(method)
            picked = None
            for rk, rv in metric_node.items():
                if _normalize_method_token(str(rk)) == method_key:
                    picked = rv
                    break
            if picked is None and "default" in metric_node:
                picked = metric_node.get("default")
            if picked is not None:
                return _normalize_split_target_level(picked, default=d)
        elif isinstance(metric_node, str):
            return _normalize_split_target_level(metric_node, default=d)

    return _normalize_split_target_level(d, default=d)


def _method_graph_min_effective_pident(cfg, method: str, metric_col: str = "split_index") -> float:
    base = np.nan
    if method in {"ASV", "zOTU"}:
        base = _as_percent_identity(
            getattr(
                cfg.figures,
                "fig09_graph_min_effective_pident_denoising",
                getattr(cfg.figures, "fig09_graph_min_identity_denoising", 100.0),
            )
        )
    elif method in {"TIC", "TAC"}:
        base = _as_percent_identity(
            getattr(
                cfg.figures,
                "fig09_graph_min_effective_pident_tic_tac",
                getattr(cfg.figures, "fig09_graph_min_identity_tic_tac", 98.7),
            )
        )
    else:
        base = _as_percent_identity(
            getattr(
                cfg.figures,
                "fig09_graph_min_effective_pident_otu",
                getattr(cfg.figures, "fig09_graph_min_identity_otu", 97.0),
            )
        )

    # Optional metric-specific override:
    # figures.fig09.graph_edge_threshold_by_metric.<metric>.<method|default>
    by_metric = getattr(cfg.figures, "fig09_graph_edge_threshold_by_metric", {}) or {}
    if not isinstance(by_metric, dict) or (not by_metric):
        return float(base)

    metric_key = _normalize_graph_metric_key(metric_col)
    metric_node = None
    for rk, rv in by_metric.items():
        if _normalize_graph_metric_key(str(rk)) == metric_key:
            metric_node = rv
            break
    # Fallback 1: explicit default metric node.
    if metric_node is None and "default" in by_metric:
        metric_node = by_metric.get("default")
    # Fallback 2: compact global method map (no metric nesting).
    if metric_node is None:
        methodish_keys = {"asv", "zotu", "otu", "tic", "tac", "default"}
        has_methodish = any(_normalize_method_token(k) in methodish_keys for k in by_metric.keys())
        if has_methodish:
            metric_node = by_metric
    if metric_node is None:
        return float(base)

    picked = None
    if isinstance(metric_node, dict):
        method_key = _normalize_method_token(method)
        for rk, rv in metric_node.items():
            if _normalize_method_token(str(rk)) == method_key:
                picked = rv
                break
        if picked is None and "default" in metric_node:
            picked = metric_node.get("default")
    elif isinstance(metric_node, (int, float, np.floating, np.integer, str)):
        picked = metric_node

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


def _load_tip_to_species(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    try:
        m = pd.read_csv(path, sep="\t")
        if {"tip", "species"}.issubset(set(m.columns)):
            return dict(zip(m["tip"].astype(str), m["species"].astype(str)))
    except Exception:
        return {}
    return {}


def _collapse_tip_table_to_species(tab: pd.DataFrame, tip_to_species: dict[str, str]) -> pd.DataFrame:
    tmp = tab.copy()
    tmp["__species__"] = [str(tip_to_species.get(str(t), str(t))) for t in tmp.index]
    return tmp.groupby("__species__", sort=False).sum(numeric_only=True)


def _feature_to_species_map(mapping: pd.DataFrame, tip_to_species: dict[str, str]) -> dict[str, str]:
    if mapping.empty:
        return {}
    fcol = "inferred_feature"
    if "target_species" in mapping.columns:
        m = mapping.copy()
        m[fcol] = m[fcol].astype(str)
        m["target_species"] = m["target_species"].astype(str)
        # soft mapping mode may contain multiple rows per inferred feature.
        # For split-index bookkeeping, use the dominant mapped species per feature.
        if "weight" in m.columns:
            m["weight"] = pd.to_numeric(m["weight"], errors="coerce").fillna(0.0)
            sort_cols = [fcol, "weight"]
            asc = [True, False]
            for c, a in (("pident", False), ("bitscore", False), ("evalue", True), ("length", False)):
                if c in m.columns:
                    sort_cols.append(c)
                    asc.append(a)
            m = m.sort_values(sort_cols, ascending=asc).drop_duplicates(subset=[fcol], keep="first")
        return dict(zip(m[fcol], m["target_species"]))
    if "target_tip" in mapping.columns:
        return {
            str(q): str(tip_to_species.get(str(t), str(t)))
            for q, t in zip(mapping[fcol].astype(str), mapping["target_tip"].astype(str))
        }
    return {}


def _normalize_effective_identity_threshold(v) -> float | None:
    try:
        f = float(v)
    except Exception:
        return None
    if not np.isfinite(f):
        return None
    if 0.0 <= f <= 1.0:
        return 100.0 * f
    return float(f)


def _parse_no_hit_effective_identity_thresholds(raw) -> dict[str, float]:
    """Parse method-specific no-hit effective-identity thresholds."""
    if not isinstance(raw, dict):
        return {}

    canon = {
        "asv": "ASV",
        "zotu": "zOTU",
        "otu": "OTU",
        "tic": "TIC",
        "tac": "TAC",
        "default": "default",
    }
    parsed: dict[str, float] = {}
    for k, v in raw.items():
        key = str(k).strip().lower()
        key = "".join(ch for ch in key if ch.isalnum() or ch == "_")
        m = canon.get(key)
        if m is None:
            continue
        thr = _normalize_effective_identity_threshold(v)
        if thr is None:
            continue
        parsed[m] = float(thr)

    out: dict[str, float] = {}
    default_thr = parsed.get("default")
    for m in METHOD_ORDER:
        if m in parsed:
            out[m] = float(parsed[m])
        elif default_thr is not None:
            out[m] = float(default_thr)
    return out


def _fig09_no_hit_effective_identity_thresholds(cfg) -> dict[str, float]:
    """
    Parse method-specific no-hit effective-identity thresholds for Figure 09.
    Priority order:
      1) fig09.no_hit_min_effective_pident_by_method
      2) fig12.no_hit_min_effective_pident_by_method (fallback)
    Returned keys are canonical METHOD_ORDER labels.
    """
    fig12_raw = getattr(cfg.figures, "fig12_no_hit_min_effective_pident_by_method", {}) or {}
    fig09_raw = getattr(cfg.figures, "fig09_no_hit_min_effective_pident_by_method", {}) or {}
    out = _parse_no_hit_effective_identity_thresholds(fig12_raw)
    out.update(_parse_no_hit_effective_identity_thresholds(fig09_raw))
    return out


def _fig09_transition_display_label(
    label: str,
    no_match_label: str = "No\nmatch",
    non_bacterial_label: str = "Non\nbacterial",
) -> str:
    mapping = {
        "Exact match": "Exact",
        "Near-exact match": "Near-\nexact",
        "No hit": str(no_match_label),
        "Unassigned": str(non_bacterial_label),
    }
    return mapping.get(str(label), str(label))


def _fig09_transition_classify_feature(
    best_eff: float,
    near_exact_threshold: float,
    exact_identity: float = 100.0,
) -> str:
    exact_identity = _as_percent_identity(float(exact_identity))
    if np.isfinite(best_eff) and best_eff >= (exact_identity - 1e-9):
        return "Exact match"
    if np.isfinite(best_eff) and (best_eff > float(near_exact_threshold)) and (best_eff < (exact_identity - 1e-9)):
        return "Near-exact match"
    return "No hit"


def _fig09_transition_pick_priority_class(classes: list[str]) -> str:
    if "Exact match" in classes:
        return "Exact match"
    if "Near-exact match" in classes:
        return "Near-exact match"
    return "No hit"


def _fig09_format_pct(value: float) -> str:
    pct = float(value) * 100.0
    if not np.isfinite(pct):
        return "NA"
    if pct == 0.0:
        return "0%"
    if abs(pct) < 0.01:
        return "<0.01%" if pct > 0 else ">-0.01%"
    if abs(pct) < 1.0:
        return f"{pct:.2f}%"
    return f"{pct:.1f}%"


def _invert_cluster_membership_map(cluster_map: dict[str, set[str]]) -> dict[str, str]:
    member_to_cluster: dict[str, str] = {}
    for cluster_id, members in cluster_map.items():
        for member_id in members:
            member_to_cluster[str(member_id)] = str(cluster_id)
    return member_to_cluster


def _best_effective_identity_by_feature_from_hits(
    hits_all: pd.DataFrame,
    min_qcov: float,
) -> dict[str, float]:
    if hits_all is None or hits_all.empty:
        return {}
    hits = filter_hits(
        hits_all.copy(),
        min_effective_pident=None,
        min_qcovhsp=min_qcov,
        strict_pident_gt=False,
    )
    if hits.empty:
        return {}
    hits = hits.copy()
    hits["qseqid"] = hits["qseqid"].astype(str)
    hits["effective_identity"] = pd.to_numeric(hits.get("effective_identity", np.nan), errors="coerce")
    hits["pident"] = pd.to_numeric(hits.get("pident", np.nan), errors="coerce")
    hits["bitscore"] = pd.to_numeric(hits.get("bitscore", np.nan), errors="coerce")
    hits["evalue"] = pd.to_numeric(hits.get("evalue", np.nan), errors="coerce")
    hits["length"] = pd.to_numeric(hits.get("length", np.nan), errors="coerce")
    hits = hits.dropna(subset=["qseqid", "effective_identity", "pident"])
    if hits.empty:
        return {}
    best_hits = (
        hits.sort_values(
            ["qseqid", "effective_identity", "pident", "bitscore", "evalue", "length"],
            ascending=[True, False, False, False, True, False],
        )
        .groupby("qseqid", sort=False)
        .first()
    )
    return {
        str(q): float(eff)
        for q, eff in zip(best_hits.index.astype(str), best_hits["effective_identity"].astype(float))
    }


def _best_hit_info_by_feature_from_hits(
    hits_all: pd.DataFrame,
    min_qcov: float,
) -> dict[str, dict[str, object]]:
    if hits_all is None or hits_all.empty:
        return {}
    hits = filter_hits(
        hits_all.copy(),
        min_effective_pident=None,
        min_qcovhsp=min_qcov,
        strict_pident_gt=False,
    )
    if hits.empty:
        return {}
    hits = hits.copy()
    hits["qseqid"] = hits["qseqid"].astype(str)
    hits["sseqid"] = hits["sseqid"].astype(str)
    hits["effective_identity"] = pd.to_numeric(hits.get("effective_identity", np.nan), errors="coerce")
    hits["pident"] = pd.to_numeric(hits.get("pident", np.nan), errors="coerce")
    hits["bitscore"] = pd.to_numeric(hits.get("bitscore", np.nan), errors="coerce")
    hits["evalue"] = pd.to_numeric(hits.get("evalue", np.nan), errors="coerce")
    hits["length"] = pd.to_numeric(hits.get("length", np.nan), errors="coerce")
    hits = hits.dropna(subset=["qseqid", "effective_identity", "pident"])
    if hits.empty:
        return {}
    best_hits = (
        hits.sort_values(
            ["qseqid", "effective_identity", "pident", "bitscore", "evalue", "length"],
            ascending=[True, False, False, False, True, False],
        )
        .groupby("qseqid", sort=False)
        .first()
    )
    return {
        str(q): {
            "target_id": str(r.get("sseqid", "")),
            "effective_identity": float(r.get("effective_identity", np.nan)),
            "pident": float(r.get("pident", np.nan)),
        }
        for q, r in best_hits.iterrows()
    }


def _build_fig09_transition_flow_summary(flow_rows: list[dict]) -> pd.DataFrame:
    if not flow_rows:
        return pd.DataFrame(
            columns=["transition_method", "zotu_class", "final_fate", "count", "row_total", "row_fraction"]
        )
    flow_df = pd.DataFrame(flow_rows)
    summary = (
        flow_df.groupby(["transition_method", "zotu_class", "final_fate"], as_index=False)
        .size()
        .rename(columns={"size": "count"})
    )
    records: list[dict[str, object]] = []
    for method in ["TIC", "TAC"]:
        sub = summary[summary["transition_method"] == method].copy()
        for origin in FIG09_TRANSITION_CLASS_ORDER:
            row_total = int(pd.to_numeric(sub.loc[sub["zotu_class"] == origin, "count"], errors="coerce").sum())
            for fate in FIG09_TRANSITION_OUTCOME_ORDER:
                hit = sub[(sub["zotu_class"] == origin) & (sub["final_fate"] == fate)]
                count = int(pd.to_numeric(hit["count"], errors="coerce").sum()) if not hit.empty else 0
                records.append(
                    {
                        "transition_method": method,
                        "zotu_class": origin,
                        "final_fate": fate,
                        "count": count,
                        "row_total": row_total,
                        "row_fraction": (float(count) / float(row_total)) if row_total > 0 else np.nan,
                    }
                )
    return pd.DataFrame(records)


def _build_fig09_homogeneity_tables(transition_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    mix_specs = [
        ("pure_exact", "Pure exact", ("Exact match",), True),
        ("pure_near", "Pure near-exact", ("Near-exact match",), True),
        ("pure_no_match", "Pure no match", ("No hit",), True),
        ("exact_near", "Exact + near-exact", ("Exact match", "Near-exact match"), False),
    ]
    cluster_columns = [
        "community_id",
        "sample_id",
        "sample_key",
        "method",
        "cluster_id",
        "n_present_members",
        "total_member_abundance_tss1000",
        "original_community_abundance_tss1000",
        "match_statuses",
        "n_match_statuses",
        "species_labels",
        "known_species_labels",
        "n_species_labels",
        "n_known_species_labels",
        "unknown_member_count",
        "is_species_homogeneous_strict",
        "mix_type",
        "mix_label",
        "taxonomy_status",
    ]
    summary_columns = [
        "mix_type",
        "mix_label",
        "method",
        "taxonomy_status",
        "cluster_count",
        "total_clusters_in_group",
        "cluster_fraction",
        "abundance_tss1000",
        "present_sample_count",
        "present_mean_abundance_fraction",
    ]
    method_columns = [
        "mix_type",
        "mix_label",
        "method",
        "clusters_in_group",
        "single_species_clusters",
        "multi_species_or_unresolved_clusters",
        "single_species_fraction",
        "single_species_present_sample_count",
        "single_species_present_mean_abundance_fraction",
    ]
    if transition_df.empty:
        return (
            pd.DataFrame(columns=cluster_columns),
            pd.DataFrame(columns=summary_columns),
            pd.DataFrame(columns=method_columns),
        )

    df = transition_df.copy()
    required_cols = {"transition_method", "cluster_id", "final_fate", "zotu_class"}
    if not required_cols.issubset(df.columns):
        return (
            pd.DataFrame(columns=cluster_columns),
            pd.DataFrame(columns=summary_columns),
            pd.DataFrame(columns=method_columns),
        )
    df = df[
        df["transition_method"].astype(str).isin(["TIC", "TAC"])
        & df["cluster_id"].astype(str).str.strip().ne("")
        & df["final_fate"].astype(str).ne("Unassigned")
    ].copy()
    if df.empty:
        return (
            pd.DataFrame(columns=cluster_columns),
            pd.DataFrame(columns=summary_columns),
            pd.DataFrame(columns=method_columns),
        )
    df["method"] = df["transition_method"].astype(str)
    df["feature_value"] = pd.to_numeric(df.get("zotu_feature_value_tss1000", np.nan), errors="coerce").fillna(0.0)
    if "best_hit_species" not in df.columns:
        df["best_hit_species"] = ""
    df["species_label"] = np.where(
        df["best_hit_species"].astype(str).str.strip().ne(""),
        df["best_hit_species"].astype(str),
        "Unknown / no best hit",
    )
    original_abundance = (
        transition_df[transition_df["transition_method"].astype(str).isin(["TIC", "TAC"])]
        .assign(
            method=lambda x: x["transition_method"].astype(str),
            feature_value=lambda x: pd.to_numeric(x.get("zotu_feature_value_tss1000", np.nan), errors="coerce").fillna(0.0),
        )
        .groupby(["community_id", "sample_id", "sample_key", "method"], sort=False)["feature_value"]
        .sum()
    )

    cluster_rows: list[dict[str, object]] = []
    for keys, grp in df.groupby(["community_id", "sample_id", "sample_key", "method", "cluster_id"], sort=False):
        community_id, sample_id, sample_key, method, cluster_id = keys
        if int(grp.shape[0]) <= 1:
            continue
        statuses = sorted(set(grp["zotu_class"].astype(str).tolist()))
        species_labels = sorted(set(grp["species_label"].astype(str).tolist()))
        known_species_labels = [label for label in species_labels if label != "Unknown / no best hit"]
        unknown_member_count = int((grp["species_label"].astype(str) == "Unknown / no best hit").sum())
        n_species_labels = len(species_labels)
        n_known_species_labels = len(set(known_species_labels))
        is_species_homogeneous_strict = (
            n_species_labels == 1
            and unknown_member_count == 0
            and n_known_species_labels == 1
        )
        total_member_abundance = float(pd.to_numeric(grp["feature_value"], errors="coerce").fillna(0.0).sum())
        original_total = float(original_abundance.get((community_id, sample_id, sample_key, method), np.nan))

        for mix_type, mix_label, required, require_exact_status_set in mix_specs:
            status_set = set(statuses)
            if require_exact_status_set:
                if status_set != set(required):
                    continue
            elif not all(status in status_set for status in required):
                continue
            taxonomy_status = "Single species" if is_species_homogeneous_strict else "Multi-species or unresolved"
            cluster_rows.append(
                {
                    "community_id": community_id,
                    "sample_id": sample_id,
                    "sample_key": sample_key,
                    "method": method,
                    "cluster_id": cluster_id,
                    "n_present_members": int(grp.shape[0]),
                    "total_member_abundance_tss1000": total_member_abundance,
                    "original_community_abundance_tss1000": original_total,
                    "match_statuses": ";".join(statuses),
                    "n_match_statuses": int(len(statuses)),
                    "species_labels": ";".join(species_labels),
                    "known_species_labels": ";".join(known_species_labels),
                    "n_species_labels": int(n_species_labels),
                    "n_known_species_labels": int(n_known_species_labels),
                    "unknown_member_count": int(unknown_member_count),
                    "is_species_homogeneous_strict": bool(is_species_homogeneous_strict),
                    "mix_type": mix_type,
                    "mix_label": mix_label,
                    "taxonomy_status": taxonomy_status,
                }
            )
    cluster_df = pd.DataFrame(cluster_rows, columns=cluster_columns)
    if cluster_df.empty:
        return cluster_df, pd.DataFrame(columns=summary_columns), pd.DataFrame(columns=method_columns)

    status_order = ["Single species", "Multi-species or unresolved"]
    count_raw = (
        cluster_df.groupby(["mix_type", "mix_label", "method", "taxonomy_status"], as_index=False)
        .size()
        .rename(columns={"size": "cluster_count"})
    )
    abundance_raw = (
        cluster_df.groupby(["mix_type", "mix_label", "method", "taxonomy_status"], as_index=False)[
            "total_member_abundance_tss1000"
        ]
        .sum()
        .rename(columns={"total_member_abundance_tss1000": "abundance_tss1000"})
    )
    present_sample = (
        cluster_df.groupby(
            ["mix_type", "mix_label", "method", "taxonomy_status", "community_id", "sample_id", "sample_key"],
            as_index=False,
        )
        .agg(
            abundance_tss1000=("total_member_abundance_tss1000", "sum"),
            original_community_abundance_tss1000=("original_community_abundance_tss1000", "first"),
        )
    )
    present_sample["sample_abundance_fraction"] = np.where(
        pd.to_numeric(present_sample["original_community_abundance_tss1000"], errors="coerce") > 0,
        pd.to_numeric(present_sample["abundance_tss1000"], errors="coerce")
        / pd.to_numeric(present_sample["original_community_abundance_tss1000"], errors="coerce"),
        np.nan,
    )
    present_raw = (
        present_sample.groupby(["mix_type", "mix_label", "method", "taxonomy_status"], as_index=False)
        .agg(
            present_sample_count=("sample_abundance_fraction", "size"),
            present_mean_abundance_fraction=("sample_abundance_fraction", "mean"),
        )
    )
    totals = (
        cluster_df.groupby(["mix_type", "mix_label", "method"], as_index=False)
        .size()
        .rename(columns={"size": "total_clusters_in_group"})
    )

    summary_records: list[dict[str, object]] = []
    method_records: list[dict[str, object]] = []
    for mix_type, mix_label, _required, _require_exact_status_set in mix_specs:
        for method in ["TIC", "TAC"]:
            total_hit = totals[(totals["mix_type"] == mix_type) & (totals["method"] == method)]
            total = int(total_hit["total_clusters_in_group"].iloc[0]) if not total_hit.empty else 0
            single_count = 0
            multi_count = 0
            single_present_n = 0
            single_present_mean = np.nan
            for status in status_order:
                hit = count_raw[
                    (count_raw["mix_type"] == mix_type)
                    & (count_raw["method"] == method)
                    & (count_raw["taxonomy_status"] == status)
                ]
                count = int(hit["cluster_count"].iloc[0]) if not hit.empty else 0
                abundance_hit = abundance_raw[
                    (abundance_raw["mix_type"] == mix_type)
                    & (abundance_raw["method"] == method)
                    & (abundance_raw["taxonomy_status"] == status)
                ]
                abundance = float(abundance_hit["abundance_tss1000"].iloc[0]) if not abundance_hit.empty else 0.0
                present_hit = present_raw[
                    (present_raw["mix_type"] == mix_type)
                    & (present_raw["method"] == method)
                    & (present_raw["taxonomy_status"] == status)
                ]
                present_n = int(present_hit["present_sample_count"].iloc[0]) if not present_hit.empty else 0
                present_mean = (
                    float(present_hit["present_mean_abundance_fraction"].iloc[0])
                    if not present_hit.empty
                    else np.nan
                )
                if status == "Single species":
                    single_count = count
                    single_present_n = present_n
                    single_present_mean = present_mean
                else:
                    multi_count = count
                summary_records.append(
                    {
                        "mix_type": mix_type,
                        "mix_label": mix_label,
                        "method": method,
                        "taxonomy_status": status,
                        "cluster_count": count,
                        "total_clusters_in_group": total,
                        "cluster_fraction": (float(count) / float(total)) if total > 0 else np.nan,
                        "abundance_tss1000": abundance,
                        "present_sample_count": present_n,
                        "present_mean_abundance_fraction": present_mean,
                    }
                )
            method_records.append(
                {
                    "mix_type": mix_type,
                    "mix_label": mix_label,
                    "method": method,
                    "clusters_in_group": total,
                    "single_species_clusters": single_count,
                    "multi_species_or_unresolved_clusters": multi_count,
                    "single_species_fraction": (float(single_count) / float(total)) if total > 0 else np.nan,
                    "single_species_present_sample_count": single_present_n,
                    "single_species_present_mean_abundance_fraction": single_present_mean,
                }
            )
    return (
        cluster_df,
        pd.DataFrame(summary_records, columns=summary_columns),
        pd.DataFrame(method_records, columns=method_columns),
    )


def _plot_fig09_homogeneity_companion(
    summary_df: pd.DataFrame,
    method_summary_df: pd.DataFrame,
    out_svg: Path,
    cfg,
) -> None:
    raw_font_size = getattr(cfg.figures, "fig09_homogeneity_font_size", None)
    if raw_font_size is None:
        raw_font_size = getattr(cfg.figures, "fig09_font_size", 10.0)
    font_size = max(6.0, float(raw_font_size))
    raw_tick_size = getattr(cfg.figures, "fig09_homogeneity_tick_size", None)
    if raw_tick_size is None:
        raw_tick_size = font_size - 2.0
    tick_size = max(5.0, float(raw_tick_size))
    title = str(
        getattr(
            cfg.figures,
            "fig09_homogeneity_title",
            "Taxonomy homogeneity across cluster match-status groups",
        )
    )
    count_ylabel = str(getattr(cfg.figures, "fig09_homogeneity_count_ylabel", "Cluster count"))
    effect_ylabel = str(getattr(cfg.figures, "fig09_homogeneity_effect_ylabel", "Abundance effect"))
    count_label = str(getattr(cfg.figures, "fig09_homogeneity_count_label", "Count"))
    effect_label = str(getattr(cfg.figures, "fig09_homogeneity_effect_label", "Abundance\neffect"))
    single_label = str(getattr(cfg.figures, "fig09_homogeneity_single_species_label", "Single species"))
    multi_label = str(
        getattr(
            cfg.figures,
            "fig09_homogeneity_multi_species_label",
            "Multi-species or unresolved",
        )
    )
    effect_axis_max = float(getattr(cfg.figures, "fig09_homogeneity_effect_axis_max", 1.0))
    if not np.isfinite(effect_axis_max) or effect_axis_max <= 0:
        effect_axis_max = 1.0
    count_axis_padding = max(1.0, float(getattr(cfg.figures, "fig09_homogeneity_count_axis_padding", 1.25)))
    show_value_labels = bool(getattr(cfg.figures, "fig09_homogeneity_show_value_labels", True))
    show_bar_type_labels = bool(getattr(cfg.figures, "fig09_homogeneity_show_bar_type_labels", True))
    show_legends = bool(getattr(cfg.figures, "fig09_homogeneity_show_legends", True))
    figsize_raw = getattr(cfg.figures, "fig09_homogeneity_figsize", None)
    if isinstance(figsize_raw, (list, tuple)) and len(figsize_raw) >= 2:
        figsize = (float(figsize_raw[0]), float(figsize_raw[1]))
    else:
        figsize = (13.2, 10.2)

    methods = ["TIC", "TAC"]
    status_order = ["Single species", "Multi-species or unresolved"]
    status_labels = {
        "Single species": single_label,
        "Multi-species or unresolved": multi_label,
    }
    colors = {
        "Single species": str(getattr(cfg.figures, "fig09_homogeneity_single_species_color", "#009E73")),
        "Multi-species or unresolved": str(
            getattr(cfg.figures, "fig09_homogeneity_multi_species_color", "#D55E00")
        ),
    }
    mix_order = [
        ("pure_exact", str(getattr(cfg.figures, "fig09_homogeneity_pure_exact_title", "Pure exact"))),
        ("pure_near", str(getattr(cfg.figures, "fig09_homogeneity_pure_near_title", "Pure near-exact"))),
        ("pure_no_match", str(getattr(cfg.figures, "fig09_homogeneity_pure_no_match_title", "Pure no match"))),
        ("exact_near", str(getattr(cfg.figures, "fig09_homogeneity_exact_near_title", "Exact + near-exact"))),
    ]

    fig, axes_grid = plt.subplots(2, 2, figsize=figsize, sharey=False)
    axes = axes_grid.ravel()
    effect_axes = []
    legend_handles = None
    for ax, (mix_type, mix_label) in zip(axes, mix_order):
        ax_effect = ax.twinx()
        effect_axes.append(ax_effect)
        xs = np.arange(len(methods), dtype=float)
        count_xs = xs - 0.17
        effect_xs = xs + 0.17
        count_bottom = np.zeros(len(methods), dtype=float)
        effect_bottom = np.zeros(len(methods), dtype=float)
        for status in status_order:
            count_vals: list[float] = []
            effect_vals: list[float] = []
            counts: list[int] = []
            effect_sample_counts: list[int] = []
            for method in methods:
                hit = summary_df[
                    (summary_df["mix_type"].astype(str) == mix_type)
                    & (summary_df["method"].astype(str) == method)
                    & (summary_df["taxonomy_status"].astype(str) == status)
                ]
                count_vals.append(float(hit["cluster_count"].iloc[0]) if not hit.empty else 0.0)
                effect_vals.append(
                    float(hit["present_mean_abundance_fraction"].iloc[0])
                    if not hit.empty and pd.notna(hit["present_mean_abundance_fraction"].iloc[0])
                    else 0.0
                )
                counts.append(int(hit["cluster_count"].iloc[0]) if not hit.empty else 0)
                effect_sample_counts.append(int(hit["present_sample_count"].iloc[0]) if not hit.empty else 0)
            count_vals_arr = np.asarray(count_vals, dtype=float)
            effect_vals_arr = np.asarray(effect_vals, dtype=float)
            bars = ax.bar(
                count_xs,
                count_vals_arr,
                bottom=count_bottom,
                width=0.30,
                color=colors[status],
                edgecolor="#333333",
                linewidth=0.7,
                label=status_labels[status],
            )
            ax_effect.bar(
                effect_xs,
                effect_vals_arr,
                bottom=effect_bottom,
                width=0.30,
                color=colors[status],
                edgecolor="#333333",
                linewidth=0.7,
                hatch=str(getattr(cfg.figures, "fig09_homogeneity_effect_hatch", "...")),
                alpha=float(getattr(cfg.figures, "fig09_homogeneity_effect_alpha", 0.78)),
            )
            if legend_handles is None:
                legend_handles = bars
            if show_value_labels:
                for x, y0, h, count in zip(count_xs, count_bottom, count_vals_arr, counts):
                    if h <= 0:
                        continue
                    ax.text(
                        x,
                        y0 + (h / 2.0),
                        f"{count}",
                        ha="center",
                        va="center",
                        fontsize=max(5.0, tick_size - 1.8),
                        color="white" if h >= 2 else "#1F1F1F",
                    )
                for x, y0, h, sample_count in zip(effect_xs, effect_bottom, effect_vals_arr, effect_sample_counts):
                    if h <= 0:
                        continue
                    ax_effect.text(
                        x,
                        y0 + (h / 2.0),
                        f"{_fig09_format_pct(h)}\nn={sample_count}",
                        ha="center",
                        va="center",
                        fontsize=max(5.0, tick_size - 2.0),
                        color="white" if h >= 0.18 else "#1F1F1F",
                        linespacing=0.9,
                    )
            count_bottom += count_vals_arr
            effect_bottom += effect_vals_arr

        for x, method in zip(xs, methods):
            hit = method_summary_df[
                (method_summary_df["mix_type"].astype(str) == mix_type)
                & (method_summary_df["method"].astype(str) == method)
            ]
            if hit.empty:
                txt = "n=0"
            else:
                row = hit.iloc[0]
                total = int(row["clusters_in_group"])
                single = int(row["single_species_clusters"])
                frac = float(row["single_species_fraction"]) if pd.notna(row["single_species_fraction"]) else np.nan
                txt = f"n={total}\n{single}/{total}"
                if np.isfinite(frac):
                    txt += f" ({_fig09_format_pct(frac)})"
            ax.text(
                x,
                1.02,
                txt,
                transform=ax.get_xaxis_transform(),
                ha="center",
                va="bottom",
                fontsize=max(5.0, tick_size - 2.0),
                color="#222222",
                linespacing=0.92,
            )

        ax.set_xticks(xs)
        ax.set_xticklabels(methods, fontsize=tick_size)
        if show_bar_type_labels:
            for x in xs:
                ax.text(
                    x - 0.17,
                    -0.075,
                    count_label,
                    ha="center",
                    va="top",
                    fontsize=max(5.0, tick_size - 2.0),
                    transform=ax.get_xaxis_transform(),
                    color="#333333",
                )
                ax.text(
                    x + 0.17,
                    -0.075,
                    effect_label,
                    ha="center",
                    va="top",
                    fontsize=max(5.0, tick_size - 2.0),
                    transform=ax.get_xaxis_transform(),
                    color="#333333",
                )
        count_top = max(1.0, float(np.max(count_bottom)) if count_bottom.size else 1.0)
        ax.set_ylim(0.0, count_top * count_axis_padding)
        ax_effect.set_ylim(0.0, effect_axis_max)
        ax_effect.yaxis.set_major_formatter(FuncFormatter(lambda val, _pos: _fig09_format_pct(val)))
        ax_effect.tick_params(axis="y", labelsize=max(5.0, tick_size - 1.0), colors="#555555")
        ax_effect.spines["right"].set_color("#555555")
        ax.set_title(mix_label, fontsize=font_size)
        ax.tick_params(axis="y", labelsize=tick_size)
        ax.grid(axis="y", color="#E6E6E6", linewidth=0.8)
        ax.set_axisbelow(True)
    axes[0].set_ylabel(count_ylabel, fontsize=font_size)
    axes[2].set_ylabel(count_ylabel, fontsize=font_size)
    effect_axes[1].set_ylabel(effect_ylabel, fontsize=font_size, color="#555555")
    effect_axes[3].set_ylabel(effect_ylabel, fontsize=font_size, color="#555555")
    axes[2].set_xlabel("Method", fontsize=font_size)
    axes[3].set_xlabel("Method", fontsize=font_size)
    if show_legends:
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(
            handles,
            labels,
            loc="lower center",
            bbox_to_anchor=(0.36, 0.01),
            ncol=2,
            frameon=False,
            fontsize=max(5.0, tick_size - 0.8),
        )
        type_handles = [
            Rectangle((0, 0), 1, 1, facecolor="#F2F2F2", edgecolor="#333333", linewidth=0.7),
            Rectangle(
                (0, 0),
                1,
                1,
                facecolor="#F2F2F2",
                edgecolor="#333333",
                linewidth=0.7,
                hatch=str(getattr(cfg.figures, "fig09_homogeneity_effect_hatch", "...")),
            ),
        ]
        fig.legend(
            type_handles,
            [str(getattr(cfg.figures, "fig09_homogeneity_count_legend_label", "Absolute count")), effect_ylabel],
            loc="lower center",
            bbox_to_anchor=(0.70, 0.01),
            ncol=2,
            frameon=False,
            fontsize=max(5.0, tick_size - 0.8),
        )
    fig.suptitle(title, fontsize=font_size + 2.0, y=0.985)
    fig.tight_layout(rect=[0.0, 0.08, 1.0, 0.955], h_pad=2.2, w_pad=1.8)
    out_svg.parent.mkdir(parents=True, exist_ok=True)
    save_kwargs: dict[str, object] = {}
    if bool(getattr(cfg.figures, "fig09_homogeneity_save_bbox_tight", False)):
        save_kwargs["bbox_inches"] = "tight"
        save_kwargs["pad_inches"] = max(0.0, float(getattr(cfg.figures, "fig09_homogeneity_save_pad_inches", 0.02)))
    fig.savefig(out_svg, format="svg", **save_kwargs)
    plt.close(fig)


def _build_fig09_relative_richness_tables(
    diagnostics_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    required = [
        "community_id",
        "sample_id",
        "sample_key",
        "method",
        "inferred_present_features",
        "expected_present_features",
    ]
    missing = [col for col in required if col not in diagnostics_df.columns]
    sample_columns = [
        "community_id",
        "sample_id",
        "sample_key",
        "method",
        "method_label",
        "expected_species_count",
        "inferred_feature_count",
        "richness_difference",
        "relative_richness_error",
        "relative_richness_error_percent",
        "richness_error_direction",
    ]
    method_columns = [
        "method",
        "method_label",
        "n_samples",
        "median_relative_richness_error_percent",
        "mean_relative_richness_error_percent",
        "q1_relative_richness_error_percent",
        "q3_relative_richness_error_percent",
        "min_relative_richness_error_percent",
        "max_relative_richness_error_percent",
        "median_richness_difference",
        "mean_richness_difference",
        "fraction_overestimated",
        "fraction_underestimated",
        "fraction_equal",
    ]
    effect_columns = [
        "community_id",
        "sample_id",
        "sample_key",
        "clustering_method",
        "clustering_method_label",
        "expected_species_count",
        "unoise3_feature_count",
        "clustering_feature_count",
        "unoise3_excess_richness",
        "excess_removed",
        "percent_excess_removed",
        "overcorrected_below_expected",
    ]
    effect_summary_columns = [
        "clustering_method",
        "clustering_method_label",
        "n_samples_with_unoise3_excess",
        "median_percent_excess_removed",
        "mean_percent_excess_removed",
        "q1_percent_excess_removed",
        "q3_percent_excess_removed",
        "min_percent_excess_removed",
        "max_percent_excess_removed",
        "fraction_overcorrected_below_expected",
    ]
    if missing:
        return (
            pd.DataFrame(columns=sample_columns),
            pd.DataFrame(columns=method_columns),
            pd.DataFrame(columns=effect_columns),
            pd.DataFrame(columns=effect_summary_columns),
        )

    sample_df = diagnostics_df[required].drop_duplicates(
        subset=["community_id", "sample_id", "sample_key", "method"]
    ).copy()
    sample_df = sample_df[sample_df["method"].astype(str).isin(METHOD_ORDER)].copy()
    sample_df["method"] = sample_df["method"].astype(str)
    sample_df["method_label"] = sample_df["method"].map(lambda method: display_methods([method])[0])
    sample_df["expected_species_count"] = pd.to_numeric(
        sample_df["expected_present_features"], errors="coerce"
    )
    sample_df["inferred_feature_count"] = pd.to_numeric(
        sample_df["inferred_present_features"], errors="coerce"
    )
    sample_df = sample_df.dropna(subset=["expected_species_count", "inferred_feature_count"])
    sample_df = sample_df[sample_df["expected_species_count"] > 0].copy()
    sample_df["relative_richness_error"] = (
        (sample_df["inferred_feature_count"] - sample_df["expected_species_count"])
        / sample_df["expected_species_count"]
    )
    sample_df["relative_richness_error_percent"] = 100.0 * sample_df["relative_richness_error"]
    sample_df["richness_difference"] = (
        sample_df["inferred_feature_count"] - sample_df["expected_species_count"]
    )
    sample_df["richness_error_direction"] = np.select(
        [
            sample_df["richness_difference"] > 0,
            sample_df["richness_difference"] < 0,
        ],
        ["overestimated", "underestimated"],
        default="equal",
    )
    sample_df = sample_df[sample_columns].copy()

    method_rows: list[dict[str, object]] = []
    for method in METHOD_ORDER:
        sub = sample_df[sample_df["method"] == method].copy()
        vals = pd.to_numeric(sub["relative_richness_error_percent"], errors="coerce").dropna()
        diffs = pd.to_numeric(sub["richness_difference"], errors="coerce").dropna()
        method_rows.append(
            {
                "method": method,
                "method_label": display_methods([method])[0],
                "n_samples": int(vals.size),
                "median_relative_richness_error_percent": float(vals.median()) if vals.size else np.nan,
                "mean_relative_richness_error_percent": float(vals.mean()) if vals.size else np.nan,
                "q1_relative_richness_error_percent": float(vals.quantile(0.25)) if vals.size else np.nan,
                "q3_relative_richness_error_percent": float(vals.quantile(0.75)) if vals.size else np.nan,
                "min_relative_richness_error_percent": float(vals.min()) if vals.size else np.nan,
                "max_relative_richness_error_percent": float(vals.max()) if vals.size else np.nan,
                "median_richness_difference": float(diffs.median()) if diffs.size else np.nan,
                "mean_richness_difference": float(diffs.mean()) if diffs.size else np.nan,
                "fraction_overestimated": float((diffs > 0).mean()) if diffs.size else np.nan,
                "fraction_underestimated": float((diffs < 0).mean()) if diffs.size else np.nan,
                "fraction_equal": float((diffs == 0).mean()) if diffs.size else np.nan,
            }
        )
    method_summary = pd.DataFrame(method_rows, columns=method_columns)

    count_wide = sample_df.pivot_table(
        index=["community_id", "sample_id", "sample_key"],
        columns="method",
        values="inferred_feature_count",
        aggfunc="first",
    )
    expected_wide = sample_df.pivot_table(
        index=["community_id", "sample_id", "sample_key"],
        values="expected_species_count",
        aggfunc="first",
    ).rename(columns={"expected_species_count": "expected_species_count"})
    wide = count_wide.join(expected_wide)
    effect_rows: list[dict[str, object]] = []
    for idx, row in wide.iterrows():
        if "zOTU" not in row.index or not np.isfinite(float(row.get("zOTU", np.nan))):
            continue
        expected_count = float(row.get("expected_species_count", np.nan))
        zotu_count = float(row.get("zOTU", np.nan))
        excess_unoise3 = zotu_count - expected_count
        if not np.isfinite(excess_unoise3) or excess_unoise3 <= 0.0:
            continue
        community_id, sample_id, sample_key = idx
        for method in ("TAC", "TIC"):
            method_count = float(row.get(method, np.nan))
            if not np.isfinite(method_count):
                continue
            excess_removed = zotu_count - method_count
            effect_rows.append(
                {
                    "community_id": community_id,
                    "sample_id": sample_id,
                    "sample_key": sample_key,
                    "clustering_method": method,
                    "clustering_method_label": display_methods([method])[0],
                    "expected_species_count": expected_count,
                    "unoise3_feature_count": zotu_count,
                    "clustering_feature_count": method_count,
                    "unoise3_excess_richness": excess_unoise3,
                    "excess_removed": excess_removed,
                    "percent_excess_removed": 100.0 * excess_removed / excess_unoise3,
                    "overcorrected_below_expected": bool(method_count < expected_count),
                }
            )
    effect_df = pd.DataFrame(effect_rows, columns=effect_columns)

    effect_summary_rows: list[dict[str, object]] = []
    for method in ("TAC", "TIC"):
        sub = effect_df[effect_df["clustering_method"] == method].copy()
        vals = pd.to_numeric(sub["percent_excess_removed"], errors="coerce").dropna()
        effect_summary_rows.append(
            {
                "clustering_method": method,
                "clustering_method_label": display_methods([method])[0],
                "n_samples_with_unoise3_excess": int(vals.size),
                "median_percent_excess_removed": float(vals.median()) if vals.size else np.nan,
                "mean_percent_excess_removed": float(vals.mean()) if vals.size else np.nan,
                "q1_percent_excess_removed": float(vals.quantile(0.25)) if vals.size else np.nan,
                "q3_percent_excess_removed": float(vals.quantile(0.75)) if vals.size else np.nan,
                "min_percent_excess_removed": float(vals.min()) if vals.size else np.nan,
                "max_percent_excess_removed": float(vals.max()) if vals.size else np.nan,
                "fraction_overcorrected_below_expected": float(
                    sub["overcorrected_below_expected"].astype(bool).mean()
                )
                if not sub.empty
                else np.nan,
            }
        )
    effect_summary = pd.DataFrame(effect_summary_rows, columns=effect_summary_columns)
    return sample_df, method_summary, effect_df, effect_summary


def _write_mapping_diagnostics_schematic(out_dir: Path, cfg, log) -> Path | None:
    src = Path(__file__).resolve().parents[3] / "fig09_mapping_diagnostics_schematic.svg"
    dest = out_dir / "fig09_mapping_diagnostics_schematic.svg"
    try:
        if not src.exists():
            log.info(f"Figure 09 schematic source not found: {src}")
            return None
        split_label = _normalize_split_index_denominator(
            getattr(cfg.figures, "fig09_split_index_denominator", "mapped_expected_features")
        )
        merge_label = _normalize_merge_denominator(getattr(cfg.figures, "fig09_merge_denominator", "mapped_inferred_features"))
        one_to_one_label = _normalize_one_to_one_denominator(
            getattr(cfg.figures, "fig09_one_to_one_denominator", "expected_present_features")
        )
        content = src.read_text(encoding="utf-8")
        content = content.replace(
            "SplitIndex = sum_T max(0, deg(T) - 1) / D_split",
            f"SplitIndex = sum_T max(0, deg(T) - 1) / D_split ({split_label})",
        )
        content = content.replace(
            "MergeIndex = sum_I max(0, deg(I) - 1) / D_merge",
            f"MergeIndex = sum_I max(0, deg(I) - 1) / D_merge ({merge_label})",
        )
        content = content.replace(
            "OneToOne = count(I with deg(I)=1 and deg(T(I))=1) / D_1to1",
            f"OneToOne = count(I with deg(I)=1 and deg(T(I))=1) / D_1to1 ({one_to_one_label})",
        )
        dest.write_text(content, encoding="utf-8")
        return dest
    except OSError as err:
        log.info(f"Figure 09 schematic copy failed: {err}")
        return None


def _manifest_path(v) -> Path | None:
    if v is None:
        return None
    try:
        if pd.isna(v):
            return None
    except Exception:
        pass
    s = str(v).strip()
    if (not s) or (s.lower() == "nan"):
        return None
    return Path(s)


def _load_cluster_membership_map(table_path: Path | None, method: str) -> dict[str, set[str]]:
    """
    Build cluster -> member-zOTU mapping for TIC/TAC.
      TIC: Map-SOTU-ZOTU.tab  (SOTU, ZOTU)
      TAC: Map-ZOTUs-OTUs-TAC.tab (ZOTU, OTU)
    """
    if table_path is None:
        return {}
    if method == "TIC":
        map_path = table_path.parent / "Map-SOTU-ZOTU.tab"
    elif method == "TAC":
        map_path = table_path.parent / "Map-ZOTUs-OTUs-TAC.tab"
    else:
        return {}
    if not map_path.exists():
        return {}

    out: dict[str, set[str]] = {}
    try:
        with map_path.open("r", encoding="utf-8") as fh:
            for i, ln in enumerate(fh):
                line = ln.strip()
                if not line:
                    continue
                parts = [p.strip() for p in line.split("\t")]
                if len(parts) < 2:
                    continue
                a, b = str(parts[0]), str(parts[1])
                # Header guards.
                if i == 0:
                    a_l = a.strip().lower().replace(" ", "_")
                    b_l = b.strip().lower().replace(" ", "_")
                    tic_header = a_l in {"sotu", "sotu_id", "#sotu_id"} and b_l in {
                        "zotu",
                        "zotu_id",
                        "#zotu_id",
                    }
                    tac_header = a_l in {"zotu", "zotu_id", "#zotu_id"} and b_l in {
                        "otu",
                        "otu_id",
                        "#otu_id",
                    }
                    if method == "TIC" and tic_header:
                        continue
                    if method == "TAC" and tac_header:
                        continue
                if method == "TIC":
                    cluster, member = a, b
                else:
                    member, cluster = a, b
                if (not cluster) or (not member):
                    continue
                out.setdefault(cluster, set()).add(member)
    except Exception:
        return {}
    return out


def _build_member_target_map(
    hits_all: pd.DataFrame,
    min_ident: float,
    min_qcov: float,
    target_level: str,
    tip_to_species: dict[str, str],
) -> dict[str, set[str]]:
    if hits_all is None or hits_all.empty:
        return {}
    h = filter_hits(
        hits_all.copy(),
        min_effective_pident=min_ident,
        min_qcovhsp=min_qcov,
        strict_pident_gt=False,
    )
    if h.empty:
        return {}
    h["qseqid"] = h["qseqid"].astype(str)
    h["sseqid"] = h["sseqid"].astype(str)
    if _normalize_split_target_level(target_level, default="species") == "species":
        h["target"] = h["sseqid"].map(lambda x: tip_to_species.get(str(x), genus_species(str(x))))
    else:
        h["target"] = h["sseqid"]
    hit_pairs = h[["qseqid", "target"]].dropna().drop_duplicates()
    tmap: dict[str, set[str]] = {}
    for q, g in hit_pairs.groupby("qseqid", sort=False):
        tmap[str(q)] = set(g["target"].astype(str).tolist())
    return tmap


def _feature_has_best_hit(
    method: str,
    feature_id: str,
    best_effective_pident_by_feature: dict[str, float],
    best_pident_by_feature: dict[str, float],
    no_hit_threshold_by_method: dict[str, float],
) -> bool:
    eff = float(best_effective_pident_by_feature.get(str(feature_id), np.nan))
    pid = float(best_pident_by_feature.get(str(feature_id), np.nan))
    score = eff if np.isfinite(eff) else pid
    if not np.isfinite(score):
        return False
    thr = no_hit_threshold_by_method.get(str(method))
    if thr is not None:
        return bool(score >= float(thr))
    return True


def _compute_cluster_purity_summary(
    present_clusters: set[str],
    cluster_to_members: dict[str, set[str]],
    present_members: set[str],
    member_targets: dict[str, set[str]],
) -> dict[str, float]:
    """
    Cluster purity (per sample/method):
    - A cluster is "pure" if all present members map and the union of mapped targets has size 1.
    - Report:
        cluster_purity_fraction = pure_clusters / total_clusters
        cluster_purity_score = mean( mapped_member_fraction * (1 / n_targets_union) )
    """
    total_clusters = 0
    pure_clusters = 0
    score_sum = 0.0

    for cluster_id in sorted([str(c) for c in present_clusters]):
        members_all = set(str(m) for m in cluster_to_members.get(cluster_id, set()))
        if not members_all:
            continue
        members_present = [m for m in members_all if m in present_members]
        members = members_present if members_present else list(members_all)
        if not members:
            continue

        total_clusters += 1
        union_targets: set[str] = set()
        mapped_members = 0
        for mem in members:
            targets = set(str(t) for t in member_targets.get(str(mem), set()))
            if targets:
                mapped_members += 1
                union_targets.update(targets)

        is_pure = (mapped_members == len(members)) and (len(union_targets) == 1)
        if is_pure:
            pure_clusters += 1

        mapped_frac = float(mapped_members / len(members))
        inv_target_multiplicity = (1.0 / float(len(union_targets))) if union_targets else 0.0
        score_sum += float(mapped_frac * inv_target_multiplicity)

    if total_clusters <= 0:
        return {
            "cluster_purity_fraction": np.nan,
            "cluster_purity_score": np.nan,
            "cluster_purity_pure_clusters": 0.0,
            "cluster_purity_total_clusters": 0.0,
        }

    return {
        "cluster_purity_fraction": float(pure_clusters / total_clusters),
        "cluster_purity_score": float(score_sum / total_clusters),
        "cluster_purity_pure_clusters": float(pure_clusters),
        "cluster_purity_total_clusters": float(total_clusters),
    }


def _bh_adjust(pvals: list[float]) -> list[float]:
    """
    Benjamini-Hochberg FDR adjustment.
    """
    out = [float("nan")] * len(pvals)
    finite_idx = [i for i, p in enumerate(pvals) if np.isfinite(p)]
    if not finite_idx:
        return out

    order = sorted(finite_idx, key=lambda i: float(pvals[i]))
    m = len(order)

    # Raw BH adjusted values in sorted p-value order.
    raw = [0.0] * m
    for rank1, idx in enumerate(order, start=1):
        raw[rank1 - 1] = float(pvals[idx]) * m / rank1

    # Enforce monotonicity from largest to smallest rank.
    mono = [0.0] * m
    running_min = 1.0
    for j in range(m - 1, -1, -1):
        running_min = min(running_min, raw[j])
        mono[j] = min(1.0, running_min)

    for j, idx in enumerate(order):
        out[idx] = mono[j]
    return out


def _pvalue_stars(p: float, alpha: float = 0.05) -> str:
    if not np.isfinite(p):
        return "ns"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < alpha:
        return "*"
    return "ns"


def _friedman_statistic(values: np.ndarray) -> float:
    """
    Friedman Q statistic with tie correction.
    values shape: n_blocks x n_methods
    """
    if values.ndim != 2:
        return float("nan")
    n, k = values.shape
    if n < 2 or k < 2:
        return float("nan")

    ranks = pd.DataFrame(values).rank(axis=1, method="average").to_numpy(dtype=float)
    rank_sums = ranks.sum(axis=0)
    q = (12.0 / (n * k * (k + 1.0))) * float(np.sum(rank_sums**2)) - 3.0 * n * (k + 1.0)

    tie_num = 0.0
    for i in range(n):
        _, counts = np.unique(values[i, :], return_counts=True)
        tie_num += float(np.sum(counts**3 - counts))
    tie_den = float(n * (k**3 - k))
    if tie_den > 0:
        c = 1.0 - (tie_num / tie_den)
        if c > 0:
            q = q / c
    return float(q)


def _friedman_permutation_pvalue(values: np.ndarray, n_perm: int = 4000, seed: int = 1) -> tuple[float, float]:
    obs = _friedman_statistic(values)
    if not np.isfinite(obs):
        return float("nan"), float("nan")

    rng = np.random.default_rng(seed)
    ge = 0
    n, k = values.shape
    n_perm = max(1, int(n_perm))

    for _ in range(n_perm):
        permuted = np.empty_like(values)
        for i in range(n):
            permuted[i, :] = values[i, rng.permutation(k)]
        h = _friedman_statistic(permuted)
        if np.isfinite(h) and h >= obs:
            ge += 1
    p = (ge + 1.0) / (n_perm + 1.0)
    return float(obs), float(p)


def _paired_signed_rank_permutation(
    x: np.ndarray,
    y: np.ndarray,
    n_perm: int = 4000,
    seed: int = 1,
) -> tuple[float, float, int, float]:
    """
    Paired signed-rank permutation test (two-sided).
    Returns: W_plus, p_value, n_nonzero_pairs, median(x-y)
    """
    d = np.asarray(x, dtype=float) - np.asarray(y, dtype=float)
    d = d[np.isfinite(d)]
    d = d[d != 0]
    n = int(d.size)
    if n == 0:
        return float("nan"), 1.0, 0, 0.0

    ranks = pd.Series(np.abs(d)).rank(method="average").to_numpy(dtype=float)
    w_plus = float(ranks[d > 0].sum())
    total = float(ranks.sum())
    center = 0.5 * total
    obs_dev = abs(w_plus - center)

    rng = np.random.default_rng(seed)
    ge = 0
    n_perm = max(1, int(n_perm))
    for _ in range(n_perm):
        signs = rng.choice(np.array([-1.0, 1.0]), size=n, replace=True)
        w_perm = float(ranks[signs > 0].sum())
        if abs(w_perm - center) >= (obs_dev - 1e-12):
            ge += 1
    p = (ge + 1.0) / (n_perm + 1.0)
    return w_plus, float(p), n, float(np.median(d))


def _paired_sign_test_exact(
    x: np.ndarray,
    y: np.ndarray,
) -> tuple[float, float, int, int, int, int, float]:
    """
    Exact paired sign test (two-sided), robust to ties/zeros by dropping d==0.
    Returns: n_pos, p_value, n_nonzero_pairs, n_pos, n_neg, n_zero, median(x-y)
    """
    d_all = np.asarray(x, dtype=float) - np.asarray(y, dtype=float)
    d_all = d_all[np.isfinite(d_all)]
    n_zero = int(np.sum(d_all == 0))
    d = d_all[d_all != 0]

    n = int(d.size)
    if n == 0:
        return float("nan"), 1.0, 0, 0, 0, n_zero, 0.0

    n_pos = int(np.sum(d > 0))
    n_neg = int(np.sum(d < 0))
    k = min(n_pos, n_neg)

    tail_count = 0
    for i in range(0, k + 1):
        tail_count += math.comb(n, i)
    p = float((2.0 * tail_count) / (2.0**n))
    p = float(min(1.0, max(0.0, p)))
    return float(n_pos), p, n, n_pos, n_neg, n_zero, float(np.median(d))


def _compute_panel_stats_permutation(
    df: pd.DataFrame,
    metric_col: str,
    method_order: list[str],
    n_perm: int,
    seed: int,
    alpha: float,
    pairwise_test: str,
    multiple_test_correction: bool,
) -> dict:
    wide = df.pivot_table(index="sample_key", columns="method", values=metric_col, aggfunc="first")
    methods = [m for m in method_order if m in wide.columns]
    if not methods:
        return {
            "n_blocks": 0,
            "methods": [],
            "friedman_q": float("nan"),
            "friedman_p": float("nan"),
            "pairwise": [],
            "pairwise_sig": [],
        }

    wide = wide[methods].dropna(axis=0, how="any")
    n_blocks = int(wide.shape[0])
    if n_blocks < 2 or len(methods) < 2:
        return {
            "n_blocks": n_blocks,
            "methods": methods,
            "friedman_q": float("nan"),
            "friedman_p": float("nan"),
            "pairwise": [],
            "pairwise_sig": [],
        }

    values = wide.to_numpy(dtype=float)
    friedman_q, friedman_p = _friedman_permutation_pvalue(values, n_perm=n_perm, seed=seed)

    pairs: list[dict] = []
    for pair_idx, (a, b) in enumerate(combinations(methods, 2), start=1):
        xa = wide[a].to_numpy(dtype=float)
        xb = wide[b].to_numpy(dtype=float)
        if pairwise_test == "sign_test_exact":
            stat, p_raw, n_nonzero, n_pos, n_neg, n_zero, median_diff = _paired_sign_test_exact(xa, xb)
            pairs.append(
                {
                    "method_a": a,
                    "method_b": b,
                    "pairwise_test": "sign_test_exact",
                    "statistic": stat,
                    "w_plus": float("nan"),
                    "p_raw": p_raw,
                    "n_nonzero_pairs": n_nonzero,
                    "n_pos_pairs": n_pos,
                    "n_neg_pairs": n_neg,
                    "n_zero_pairs": n_zero,
                    "median_diff_a_minus_b": median_diff,
                }
            )
        else:
            w_plus, p_raw, n_nonzero, median_diff = _paired_signed_rank_permutation(
                xa,
                xb,
                n_perm=n_perm,
                seed=seed + 1000 + pair_idx,
            )
            pairs.append(
                {
                    "method_a": a,
                    "method_b": b,
                    "pairwise_test": "signed_rank_permutation",
                    "statistic": w_plus,
                    "w_plus": w_plus,
                    "p_raw": p_raw,
                    "n_nonzero_pairs": n_nonzero,
                    "n_pos_pairs": float("nan"),
                    "n_neg_pairs": float("nan"),
                    "n_zero_pairs": float("nan"),
                    "median_diff_a_minus_b": median_diff,
                }
            )

    p_raw_vals = [float(r["p_raw"]) for r in pairs]
    p_adj_bh_vals = _bh_adjust(p_raw_vals)
    p_used_vals = p_adj_bh_vals if multiple_test_correction else p_raw_vals
    for i, p_used in enumerate(p_used_vals):
        p_adj = p_adj_bh_vals[i]
        pairs[i]["p_adj_bh"] = p_adj
        pairs[i]["p_value_used"] = p_used
        pairs[i]["multiple_test_correction"] = bool(multiple_test_correction)
        pairs[i]["p_adjustment"] = "bh" if multiple_test_correction else "none"
        pairs[i]["significant_bh"] = bool(np.isfinite(p_adj) and p_adj < alpha)
        pairs[i]["significant_alpha"] = bool(np.isfinite(p_used) and p_used < alpha)

    sig_pairs = [r for r in pairs if bool(r["significant_alpha"])]
    return {
        "n_blocks": n_blocks,
        "methods": methods,
        "friedman_q": friedman_q,
        "friedman_p": friedman_p,
        "multiple_test_correction": bool(multiple_test_correction),
        "stats_model": "paired_permutation",
        "pairwise": pairs,
        "pairwise_sig": sig_pairs,
    }


def _apply_pairwise_adjustments(pairs: list[dict], alpha: float, multiple_test_correction: bool) -> tuple[list[dict], list[dict]]:
    p_raw_vals = [float(r.get("p_raw", np.nan)) for r in pairs]
    p_adj_bh_vals = _bh_adjust(p_raw_vals)
    p_used_vals = p_adj_bh_vals if multiple_test_correction else p_raw_vals
    for i, p_used in enumerate(p_used_vals):
        p_adj = p_adj_bh_vals[i]
        pairs[i]["p_adj_bh"] = p_adj
        pairs[i]["p_value_used"] = p_used
        pairs[i]["multiple_test_correction"] = bool(multiple_test_correction)
        pairs[i]["p_adjustment"] = "bh" if multiple_test_correction else "none"
        pairs[i]["significant_bh"] = bool(np.isfinite(p_adj) and p_adj < alpha)
        pairs[i]["significant_alpha"] = bool(np.isfinite(p_used) and p_used < alpha)
    sig_pairs = [r for r in pairs if bool(r["significant_alpha"])]
    return pairs, sig_pairs


def _prepare_mixed_effect_input(df: pd.DataFrame, metric_col: str) -> tuple[pd.DataFrame, str]:
    base = df.copy()
    base["community_id"] = base["community_id"].astype(str)
    base["sample_key"] = base["sample_key"].astype(str)
    base["method"] = base["method"].astype(str)

    if metric_col == "unmapped_feature_fraction":
        d = base[
            pd.to_numeric(base["inferred_present_features"], errors="coerce").notna()
            & (pd.to_numeric(base["inferred_present_features"], errors="coerce") > 0)
        ].copy()
        d["success"] = pd.to_numeric(d["unmapped_present_features"], errors="coerce").fillna(0.0).clip(lower=0.0)
        d["trials"] = pd.to_numeric(d["inferred_present_features"], errors="coerce").fillna(0.0).clip(lower=0.0)
        d["success"] = np.minimum(d["success"], d["trials"])
        d["failure"] = d["trials"] - d["success"]
        d = d[(d["trials"] > 0) & np.isfinite(d["success"]) & np.isfinite(d["failure"])].copy()
        return d[["community_id", "sample_key", "method", "success", "failure"]], "binomial_counts"

    if metric_col == "one_to_one_mapping_fraction":
        d = base[
            pd.to_numeric(base["one_to_one_denominator_features"], errors="coerce").notna()
            & (pd.to_numeric(base["one_to_one_denominator_features"], errors="coerce") > 0)
            & pd.to_numeric(base["one_to_one_mapping_fraction"], errors="coerce").notna()
        ].copy()
        d["trials"] = pd.to_numeric(d["one_to_one_denominator_features"], errors="coerce").fillna(0.0).clip(lower=0.0)
        frac = pd.to_numeric(d["one_to_one_mapping_fraction"], errors="coerce").fillna(0.0)
        d["success"] = np.rint((frac * d["trials"]).to_numpy(dtype=float))
        d["success"] = np.minimum(np.maximum(d["success"], 0.0), d["trials"])
        d["failure"] = d["trials"] - d["success"]
        d = d[(d["trials"] > 0) & np.isfinite(d["success"]) & np.isfinite(d["failure"])].copy()
        return d[["community_id", "sample_key", "method", "success", "failure"]], "binomial_counts"

    if metric_col == "no_hit_over_unmapped_fraction":
        d = base[
            pd.to_numeric(base.get("no_hit_n_unmapped", np.nan), errors="coerce").notna()
            & (pd.to_numeric(base.get("no_hit_n_unmapped", np.nan), errors="coerce") > 0)
            & pd.to_numeric(base.get("no_hit_n_no_hit", np.nan), errors="coerce").notna()
        ].copy()
        d["trials"] = pd.to_numeric(d.get("no_hit_n_unmapped", np.nan), errors="coerce").fillna(0.0).clip(lower=0.0)
        # GLMM orientation:
        # success = unmapped features with a best hit
        # failure = no-hit features
        d["failure"] = pd.to_numeric(d.get("no_hit_n_no_hit", np.nan), errors="coerce").fillna(0.0).clip(lower=0.0)
        d["failure"] = np.minimum(d["failure"], d["trials"])
        d["success"] = d["trials"] - d["failure"]
        d["success"] = np.minimum(d["success"], d["trials"])
        d = d[(d["trials"] > 0) & np.isfinite(d["success"]) & np.isfinite(d["failure"])].copy()
        return d[["community_id", "sample_key", "method", "success", "failure"]], "binomial_counts"

    if metric_col == "cluster_purity_fraction":
        d = base[
            pd.to_numeric(base.get("cluster_purity_total_clusters", np.nan), errors="coerce").notna()
            & (pd.to_numeric(base.get("cluster_purity_total_clusters", np.nan), errors="coerce") > 0)
            & pd.to_numeric(base.get("cluster_purity_pure_clusters", np.nan), errors="coerce").notna()
        ].copy()
        d["trials"] = pd.to_numeric(
            d.get("cluster_purity_total_clusters", np.nan), errors="coerce"
        ).fillna(0.0).clip(lower=0.0)
        d["success"] = pd.to_numeric(
            d.get("cluster_purity_pure_clusters", np.nan), errors="coerce"
        ).fillna(0.0).clip(lower=0.0)
        d["success"] = np.minimum(d["success"], d["trials"])
        d["failure"] = d["trials"] - d["success"]
        d = d[(d["trials"] > 0) & np.isfinite(d["success"]) & np.isfinite(d["failure"])].copy()
        return d[["community_id", "sample_key", "method", "success", "failure"]], "binomial_counts"

    if metric_col == "split_index":
        d = base[
            pd.to_numeric(base["split_index_denominator_features"], errors="coerce").notna()
            & (pd.to_numeric(base["split_index_denominator_features"], errors="coerce") > 0)
            & pd.to_numeric(base["split_duplicates"], errors="coerce").notna()
        ].copy()
        d["denom"] = pd.to_numeric(d["split_index_denominator_features"], errors="coerce").fillna(0.0).clip(lower=0.0)
        d["count"] = pd.to_numeric(d["split_duplicates"], errors="coerce").fillna(0.0).clip(lower=0.0)
        d = d[(d["denom"] > 0) & np.isfinite(d["count"]) & np.isfinite(d["denom"])].copy()
        return d[["community_id", "sample_key", "method", "count", "denom"]], "negbin_offset"

    if metric_col == "merge_deficit":
        merge_denom_col = (
            "merge_mapped_present_features"
            if "merge_mapped_present_features" in base.columns
            else "split_mapped_present_features"
        )
        d = base[
            pd.to_numeric(base[merge_denom_col], errors="coerce").notna()
            & (pd.to_numeric(base[merge_denom_col], errors="coerce") > 0)
            & pd.to_numeric(base["graph_merge_duplicates"], errors="coerce").notna()
        ].copy()
        d["denom"] = pd.to_numeric(d[merge_denom_col], errors="coerce").fillna(0.0).clip(lower=0.0)
        d["count"] = pd.to_numeric(d["graph_merge_duplicates"], errors="coerce").fillna(0.0).clip(lower=0.0)
        d = d[(d["denom"] > 0) & np.isfinite(d["count"]) & np.isfinite(d["denom"])].copy()
        return d[["community_id", "sample_key", "method", "count", "denom"]], "negbin_offset"

    return pd.DataFrame(), ""


def _run_mixed_effect_r(
    model_input: pd.DataFrame,
    method_order: list[str],
    family_mode: str,
) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    r_code = textwrap.dedent(
        """
        suppressWarnings(suppressMessages({
          library(MASS)
          library(nlme)
        }))

        args <- commandArgs(trailingOnly = TRUE)
        in_tsv <- args[[1]]
        out_omni <- args[[2]]
        out_pair <- args[[3]]
        methods <- strsplit(args[[4]], ",", fixed = TRUE)[[1]]
        family_mode <- args[[5]]

        d <- read.delim(in_tsv, sep = "\\t", stringsAsFactors = FALSE, check.names = FALSE)
        d$community_id <- factor(d$community_id)
        d$sample_key <- factor(d$sample_key)
        d$method <- factor(d$method, levels = methods)
        d <- d[!is.na(d$method), , drop = FALSE]

        if (nrow(d) < 2 || length(unique(as.character(d$method))) < 2) {
          write.table(data.frame(n_blocks=length(unique(as.character(d$sample_key))),
                                 n_methods=length(unique(as.character(d$method))),
                                 omnibus_stat=NA_real_,
                                 omnibus_df=NA_real_,
                                 omnibus_p=NA_real_,
                                 model_family=family_mode,
                                 model_fit=FALSE,
                                 note="insufficient_data"),
                      file=out_omni, sep="\\t", row.names=FALSE, quote=FALSE)
          write.table(data.frame(method_a=character(), method_b=character(), estimate=numeric(),
                                 se=numeric(), statistic=numeric(), p_raw=numeric(), pairwise_test=character()),
                      file=out_pair, sep="\\t", row.names=FALSE, quote=FALSE)
          quit(save="no", status=0)
        }

        fit <- NULL
        fit_ok <- TRUE
        fit_note <- "ok"

        if (family_mode == "binomial_counts") {
          d$success <- as.numeric(d$success)
          d$failure <- as.numeric(d$failure)
          fit <- tryCatch(
            glmmPQL(cbind(success, failure) ~ method,
                    random = list(community_id = ~1, sample_key = ~1),
                    family = binomial,
                    data = d,
                    verbose = FALSE),
            error = function(e) e
          )
        } else {
          d$count <- as.numeric(d$count)
          d$denom <- as.numeric(d$denom)
          d$log_offset <- log(pmax(d$denom, 1e-9))
          theta <- tryCatch(glm.nb(count ~ method + offset(log_offset), data=d)$theta, error=function(e) NA_real_)
          if (!is.finite(theta) || theta <= 0) theta <- 1.0
          fit <- tryCatch(
            glmmPQL(count ~ method + offset(log_offset),
                    random = list(community_id = ~1, sample_key = ~1),
                    family = negative.binomial(theta),
                    data = d,
                    verbose = FALSE),
            error = function(e) e
          )
        }

        if (inherits(fit, "error")) {
          fit_ok <- FALSE
          fit_note <- conditionMessage(fit)
          write.table(data.frame(n_blocks=length(unique(as.character(d$sample_key))),
                                 n_methods=length(unique(as.character(d$method))),
                                 omnibus_stat=NA_real_,
                                 omnibus_df=NA_real_,
                                 omnibus_p=NA_real_,
                                 model_family=family_mode,
                                 model_fit=FALSE,
                                 note=fit_note),
                      file=out_omni, sep="\\t", row.names=FALSE, quote=FALSE)
          write.table(data.frame(method_a=character(), method_b=character(), estimate=numeric(),
                                 se=numeric(), statistic=numeric(), p_raw=numeric(), pairwise_test=character()),
                      file=out_pair, sep="\\t", row.names=FALSE, quote=FALSE)
          quit(save="no", status=0)
        }

        b <- nlme::fixed.effects(fit)
        V <- as.matrix(fit$varFix)
        meth_terms <- grep("^method", names(b), value=TRUE)
        if (length(meth_terms) > 0) {
          bb <- b[meth_terms]
          VV <- V[meth_terms, meth_terms, drop=FALSE]
          W <- as.numeric(t(bb) %*% MASS::ginv(VV) %*% bb)
          p_omni <- stats::pchisq(W, df=length(meth_terms), lower.tail=FALSE)
          df_omni <- as.numeric(length(meth_terms))
        } else {
          W <- 0.0
          p_omni <- 1.0
          df_omni <- 0.0
        }

        write.table(data.frame(n_blocks=length(unique(as.character(d$sample_key))),
                               n_methods=length(unique(as.character(d$method))),
                               omnibus_stat=W,
                               omnibus_df=df_omni,
                               omnibus_p=p_omni,
                               model_family=family_mode,
                               model_fit=fit_ok,
                               note=fit_note),
                    file=out_omni, sep="\\t", row.names=FALSE, quote=FALSE)

        method_levels <- methods[methods %in% unique(as.character(d$method))]
        ref <- method_levels[1]
        term_name <- function(m) {
          if (m == ref) return(NA_character_)
          paste0("method", m)
        }

        pair_rows <- list()
        idx <- 1
        for (i in seq_len(length(method_levels) - 1)) {
          for (j in seq.int(i + 1, length(method_levels))) {
            a <- method_levels[i]
            b2 <- method_levels[j]
            L <- rep(0.0, length(b))
            names(L) <- names(b)
            ta <- term_name(a)
            tb <- term_name(b2)
            if (!is.na(ta) && ta %in% names(L)) L[ta] <- L[ta] + 1.0
            if (!is.na(tb) && tb %in% names(L)) L[tb] <- L[tb] - 1.0
            est <- as.numeric(sum(L * b))
            se2 <- as.numeric(t(L) %*% V %*% L)
            if (!is.finite(se2) || se2 <= 0) {
              se <- NA_real_
              z <- NA_real_
              p <- NA_real_
            } else {
              se <- sqrt(se2)
              z <- est / se
              p <- 2.0 * stats::pnorm(-abs(z))
            }
            pair_rows[[idx]] <- data.frame(
              method_a = a,
              method_b = b2,
              estimate = est,
              se = se,
              statistic = z,
              p_raw = p,
              pairwise_test = "mixed_glmm_pql_wald",
              stringsAsFactors = FALSE
            )
            idx <- idx + 1
          }
        }
        pair_df <- if (length(pair_rows) > 0) do.call(rbind, pair_rows) else data.frame(
          method_a=character(), method_b=character(), estimate=numeric(), se=numeric(),
          statistic=numeric(), p_raw=numeric(), pairwise_test=character()
        )
        write.table(pair_df, file=out_pair, sep="\\t", row.names=FALSE, quote=FALSE)
        """
    )

    with tempfile.TemporaryDirectory(prefix="tacticbench_fig09_mixed_") as td:
        td_path = Path(td)
        in_tsv = td_path / "fig09_input.tsv"
        out_omni = td_path / "fig09_omnibus.tsv"
        out_pair = td_path / "fig09_pairwise.tsv"
        r_script = td_path / "fit_fig09_mixed.R"
        model_input.to_csv(in_tsv, sep="\t", index=False)
        r_script.write_text(r_code, encoding="utf-8")
        methods_csv = ",".join([str(m) for m in method_order])
        cp = subprocess.run(
            ["Rscript", str(r_script), str(in_tsv), str(out_omni), str(out_pair), methods_csv, family_mode],
            capture_output=True,
            text=True,
            check=False,
        )
        if cp.returncode != 0:
            msg = cp.stderr.strip() or cp.stdout.strip() or f"Rscript exited with code {cp.returncode}"
            return pd.DataFrame(), pd.DataFrame(), msg
        if (not out_omni.exists()) or (not out_pair.exists()):
            return pd.DataFrame(), pd.DataFrame(), "Rscript did not write expected outputs."
        omni_df = pd.read_csv(out_omni, sep="\t")
        pair_df = pd.read_csv(out_pair, sep="\t")
        return omni_df, pair_df, ""


def _run_mixed_effect_r_glmmtmb(
    model_input: pd.DataFrame,
    method_order: list[str],
    family_mode: str,
    glmm_spec: dict[str, str],
) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    r_code = textwrap.dedent(
        """
        suppressWarnings(suppressMessages({
          has_glmmTMB <- requireNamespace("glmmTMB", quietly = TRUE)
          has_emmeans <- requireNamespace("emmeans", quietly = TRUE)
        }))

        args <- commandArgs(trailingOnly = TRUE)
        in_tsv <- args[[1]]
        out_omni <- args[[2]]
        out_pair <- args[[3]]
        methods <- strsplit(args[[4]], ",", fixed = TRUE)[[1]]
        family_mode <- args[[5]]
        bin_formula_txt <- args[[6]]
        bin_null_formula_txt <- args[[7]]
        count_formula_txt <- args[[8]]
        count_null_formula_txt <- args[[9]]
        count_family <- args[[10]]
        emmeans_scale <- args[[11]]

        bin_formula <- stats::as.formula(bin_formula_txt)
        bin_null_formula <- stats::as.formula(bin_null_formula_txt)
        count_formula <- stats::as.formula(count_formula_txt)
        count_null_formula <- stats::as.formula(count_null_formula_txt)

        d <- read.delim(in_tsv, sep = "\\t", stringsAsFactors = FALSE, check.names = FALSE)
        d$community_id <- factor(d$community_id)
        d$sample_key <- factor(d$sample_key)
        d$method <- factor(d$method, levels = methods)
        d <- d[!is.na(d$method), , drop = FALSE]
        n_blocks <- length(unique(as.character(d$sample_key)))
        n_methods <- length(unique(as.character(d$method)))

        write_empty <- function(note, model_fit = FALSE, model_family = family_mode) {
          write.table(data.frame(n_blocks=n_blocks,
                                 n_methods=n_methods,
                                 omnibus_stat=NA_real_,
                                 omnibus_df=NA_real_,
                                 omnibus_p=NA_real_,
                                 model_family=model_family,
                                 model_fit=model_fit,
                                 note=note),
                      file=out_omni, sep="\\t", row.names=FALSE, quote=FALSE)
          write.table(data.frame(method_a=character(), method_b=character(), estimate=numeric(),
                                 se=numeric(), statistic=numeric(), p_raw=numeric(), pairwise_test=character()),
                      file=out_pair, sep="\\t", row.names=FALSE, quote=FALSE)
        }

        if (!has_glmmTMB || !has_emmeans) {
          missing <- c()
          if (!has_glmmTMB) missing <- c(missing, "glmmTMB")
          if (!has_emmeans) missing <- c(missing, "emmeans")
          write_empty(paste0("missing_R_packages:", paste(missing, collapse=",")), FALSE)
          quit(save="no", status=0)
        }

        if (nrow(d) < 2 || n_methods < 2) {
          write_empty("insufficient_data", FALSE)
          quit(save="no", status=0)
        }

        fit <- NULL
        fit0 <- NULL
        fit_ok <- TRUE
        fit_note <- "ok"

        if (family_mode == "binomial_counts") {
          d$success <- as.numeric(d$success)
          d$failure <- as.numeric(d$failure)
          fit <- tryCatch(
            glmmTMB::glmmTMB(bin_formula,
                             family = stats::binomial(),
                             data = d),
            error = function(e) e
          )
          fit0 <- tryCatch(
            glmmTMB::glmmTMB(bin_null_formula,
                             family = stats::binomial(),
                             data = d),
            error = function(e) e
          )
        } else {
          d$count <- as.numeric(d$count)
          d$denom <- as.numeric(d$denom)
          d$log_offset <- log(pmax(d$denom, 1e-9))
          fam_obj <- switch(
            count_family,
            "nbinom1" = glmmTMB::nbinom1(),
            "poisson" = stats::poisson(),
            glmmTMB::nbinom2()
          )
          fit <- tryCatch(
            glmmTMB::glmmTMB(count_formula,
                             family = fam_obj,
                             data = d),
            error = function(e) e
          )
          fit0 <- tryCatch(
            glmmTMB::glmmTMB(count_null_formula,
                             family = fam_obj,
                             data = d),
            error = function(e) e
          )
        }

        if (inherits(fit, "error") || inherits(fit0, "error")) {
          msg <- if (inherits(fit, "error")) conditionMessage(fit) else conditionMessage(fit0)
          write_empty(msg, FALSE)
          quit(save="no", status=0)
        }

        omnibus_stat <- NA_real_
        omnibus_df <- NA_real_
        omnibus_p <- NA_real_
        lrt <- tryCatch(stats::anova(fit0, fit, test = "Chisq"), error = function(e) e)
        if (inherits(lrt, "error")) {
          fit_ok <- FALSE
          fit_note <- conditionMessage(lrt)
        } else if (nrow(lrt) < 2) {
          fit_ok <- FALSE
          fit_note <- "anova_output_unexpected"
        } else {
          lrt_row <- lrt[2, , drop = FALSE]
          chisq_col <- intersect(c("Chisq", "LRT", "Deviance"), colnames(lrt_row))
          df_col <- intersect(c("Df", "df"), colnames(lrt_row))
          p_col <- grep("^Pr\\\\(", colnames(lrt_row), value = TRUE)
          if (length(chisq_col) > 0) omnibus_stat <- as.numeric(lrt_row[[chisq_col[1]]])
          if (length(df_col) > 0) omnibus_df <- as.numeric(lrt_row[[df_col[1]]])
          if (length(p_col) > 0) omnibus_p <- as.numeric(lrt_row[[p_col[1]]])
          if (!is.finite(omnibus_p)) {
            # Keep model fit valid for pairwise inference even when omnibus
            # LRT p-value is unavailable in this R/glmmTMB setup.
            fit_note <- "anova_pvalue_missing"
          }
        }

        method_levels <- methods[methods %in% unique(as.character(d$method))]
        pair_df <- data.frame(method_a=character(), method_b=character(), estimate=numeric(),
                              se=numeric(), statistic=numeric(), p_raw=numeric(), pairwise_test=character(),
                              stringsAsFactors=FALSE)
        if (fit_ok && length(method_levels) >= 2) {
          em <- tryCatch(emmeans::emmeans(fit, specs = ~ method, type = emmeans_scale), error = function(e) e)
          if (inherits(em, "error")) {
            fit_ok <- FALSE
            fit_note <- conditionMessage(em)
          } else {
            pw <- tryCatch(
              as.data.frame(emmeans::contrast(em, method = "pairwise", adjust = "none")),
              error = function(e) e
            )
            if (inherits(pw, "error")) {
              fit_ok <- FALSE
              fit_note <- conditionMessage(pw)
            } else if (nrow(pw) > 0) {
              if (!("contrast" %in% colnames(pw))) {
                fit_ok <- FALSE
                fit_note <- "emmeans_pairs_missing_contrast"
              } else {
                pw$contrast <- as.character(pw$contrast)
                pw$method_a <- trimws(sub(" - .*", "", pw$contrast))
                pw$method_b <- trimws(sub(".* - ", "", pw$contrast))
                stat_col <- if ("z.ratio" %in% colnames(pw)) "z.ratio" else if ("t.ratio" %in% colnames(pw)) "t.ratio" else NA_character_
                se_col <- if ("SE" %in% colnames(pw)) "SE" else if ("se" %in% colnames(pw)) "se" else NA_character_
                est_col <- if ("estimate" %in% colnames(pw)) "estimate" else NA_character_
                p_col <- if ("p.value" %in% colnames(pw)) "p.value" else NA_character_
                pair_df <- data.frame(
                  method_a = pw$method_a,
                  method_b = pw$method_b,
                  estimate = if (!is.na(est_col)) as.numeric(pw[[est_col]]) else NA_real_,
                  se = if (!is.na(se_col)) as.numeric(pw[[se_col]]) else NA_real_,
                  statistic = if (!is.na(stat_col)) as.numeric(pw[[stat_col]]) else NA_real_,
                  p_raw = if (!is.na(p_col)) as.numeric(pw[[p_col]]) else NA_real_,
                  pairwise_test = "mixed_glmmtmb_emmeans",
                  stringsAsFactors = FALSE
                )
                pair_df <- pair_df[
                  pair_df$method_a %in% method_levels & pair_df$method_b %in% method_levels,
                  ,
                  drop = FALSE
                ]
              }
            }
          }
        }

        write.table(data.frame(n_blocks=n_blocks,
                               n_methods=n_methods,
                               omnibus_stat=omnibus_stat,
                               omnibus_df=omnibus_df,
                               omnibus_p=omnibus_p,
                               model_family=family_mode,
                               model_fit=fit_ok,
                               note=fit_note),
                    file=out_omni, sep="\\t", row.names=FALSE, quote=FALSE)
        if (!fit_ok) {
          pair_df <- data.frame(method_a=character(), method_b=character(), estimate=numeric(),
                                se=numeric(), statistic=numeric(), p_raw=numeric(), pairwise_test=character(),
                                stringsAsFactors=FALSE)
        }
        write.table(pair_df, file=out_pair, sep="\\t", row.names=FALSE, quote=FALSE)
        """
    )

    with tempfile.TemporaryDirectory(prefix="tacticbench_fig09_glmmtmb_") as td:
        td_path = Path(td)
        in_tsv = td_path / "fig09_input.tsv"
        out_omni = td_path / "fig09_omnibus.tsv"
        out_pair = td_path / "fig09_pairwise.tsv"
        r_script = td_path / "fit_fig09_glmmtmb.R"
        model_input.to_csv(in_tsv, sep="\t", index=False)
        r_script.write_text(r_code, encoding="utf-8")
        methods_csv = ",".join([str(m) for m in method_order])
        bin_formula = str(glmm_spec.get("binomial_formula", "") or "")
        bin_null_formula = str(glmm_spec.get("binomial_null_formula", "") or "")
        count_formula = str(glmm_spec.get("count_formula", "") or "")
        count_null_formula = str(glmm_spec.get("count_null_formula", "") or "")
        count_family = str(glmm_spec.get("count_family", "nbinom2") or "nbinom2")
        emmeans_scale = str(glmm_spec.get("emmeans_scale", "link") or "link")
        cp = subprocess.run(
            [
                "Rscript",
                str(r_script),
                str(in_tsv),
                str(out_omni),
                str(out_pair),
                methods_csv,
                family_mode,
                bin_formula,
                bin_null_formula,
                count_formula,
                count_null_formula,
                count_family,
                emmeans_scale,
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        if cp.returncode != 0:
            msg = cp.stderr.strip() or cp.stdout.strip() or f"Rscript exited with code {cp.returncode}"
            return pd.DataFrame(), pd.DataFrame(), msg
        if (not out_omni.exists()) or (not out_pair.exists()):
            return pd.DataFrame(), pd.DataFrame(), "Rscript did not write expected outputs."
        omni_df = pd.read_csv(out_omni, sep="\t")
        pair_df = pd.read_csv(out_pair, sep="\t")
        return omni_df, pair_df, ""


def _compute_panel_stats_mixed_effect(
    df: pd.DataFrame,
    metric_col: str,
    method_order: list[str],
    alpha: float,
    multiple_test_correction: bool,
    engine: str,
    glmm_spec: dict[str, str],
    log,
) -> dict | None:
    engine_norm = str(engine or "").strip().lower()
    use_glmmtmb = engine_norm == "glmmtmb"
    stats_model_name = "mixed_effect_glmmtmb" if use_glmmtmb else "mixed_effect"
    pairwise_default = "mixed_glmmtmb_emmeans" if use_glmmtmb else "mixed_glmm_pql_wald"

    model_input, family_mode = _prepare_mixed_effect_input(df, metric_col)
    if model_input.empty or (not family_mode):
        return None

    methods = [m for m in method_order if m in set(model_input["method"].astype(str).tolist())]
    n_blocks = int(model_input["sample_key"].astype(str).nunique())
    if len(methods) < 2 or n_blocks < 2:
        return {
            "n_blocks": n_blocks,
            "methods": methods,
            "friedman_q": float("nan"),
            "friedman_p": float("nan"),
            "multiple_test_correction": bool(multiple_test_correction),
            "stats_model": stats_model_name,
            "pairwise": [],
            "pairwise_sig": [],
        }

    if use_glmmtmb:
        omni_df, pair_df, err = _run_mixed_effect_r_glmmtmb(
            model_input=model_input, method_order=methods, family_mode=family_mode, glmm_spec=glmm_spec
        )
    else:
        omni_df, pair_df, err = _run_mixed_effect_r(model_input=model_input, method_order=methods, family_mode=family_mode)

    if err:
        if log is not None:
            log.info(f"Figure 09 mixed-effect ({engine_norm}) stats failed for {metric_col}: {err}")
        return None
    if omni_df.empty:
        return None

    omni = omni_df.iloc[0]
    fit_ok = _as_bool_scalar(omni.get("model_fit", False))
    if not fit_ok:
        note = str(omni.get("note", "model_fit_false"))
        if log is not None:
            log.info(f"Figure 09 mixed-effect ({engine_norm}) model not fitted for {metric_col}: {note}")
        return None

    pairs: list[dict] = []
    if not pair_df.empty:
        for _, row in pair_df.iterrows():
            pairs.append(
                {
                    "method_a": str(row.get("method_a", "")),
                    "method_b": str(row.get("method_b", "")),
                    "pairwise_test": str(row.get("pairwise_test", pairwise_default)),
                    "statistic": float(row.get("statistic", np.nan)),
                    "w_plus": float("nan"),
                    "p_raw": float(row.get("p_raw", np.nan)),
                    "n_nonzero_pairs": float("nan"),
                    "n_pos_pairs": float("nan"),
                    "n_neg_pairs": float("nan"),
                    "n_zero_pairs": float("nan"),
                    "median_diff_a_minus_b": float(row.get("estimate", np.nan)),
                    "estimate": float(row.get("estimate", np.nan)),
                    "std_error": float(row.get("se", np.nan)),
                    "model_family": str(omni.get("model_family", family_mode)),
                }
            )
    pairs, sig_pairs = _apply_pairwise_adjustments(
        pairs=pairs, alpha=alpha, multiple_test_correction=multiple_test_correction
    )

    return {
        "n_blocks": int(omni.get("n_blocks", n_blocks)),
        "methods": methods,
        "friedman_q": float(omni.get("omnibus_stat", np.nan)),
        "friedman_p": float(omni.get("omnibus_p", np.nan)),
        "multiple_test_correction": bool(multiple_test_correction),
        "stats_model": stats_model_name,
        "pairwise": pairs,
        "pairwise_sig": sig_pairs,
    }


def _compute_panel_stats(
    df: pd.DataFrame,
    metric_col: str,
    method_order: list[str],
    n_perm: int,
    seed: int,
    alpha: float,
    pairwise_test: str,
    multiple_test_correction: bool,
    stats_model: str,
    glmm_spec: dict[str, str],
    log,
) -> dict:
    if stats_model == "mixed_effect_glmmtmb":
        mixed = _compute_panel_stats_mixed_effect(
            df=df,
            metric_col=metric_col,
            method_order=method_order,
            alpha=alpha,
            multiple_test_correction=multiple_test_correction,
            engine="glmmtmb",
            glmm_spec=glmm_spec,
            log=log,
        )
        if mixed is not None:
            return mixed
        if log is not None:
            log.info(
                f"Figure 09: falling back to mixed_effect (PQL) stats for panel '{metric_col}' after glmmTMB failure."
            )
        mixed = _compute_panel_stats_mixed_effect(
            df=df,
            metric_col=metric_col,
            method_order=method_order,
            alpha=alpha,
            multiple_test_correction=multiple_test_correction,
            engine="pql",
            glmm_spec=glmm_spec,
            log=log,
        )
        if mixed is not None:
            return mixed
        if log is not None:
            log.info(f"Figure 09: falling back to paired_permutation stats for panel '{metric_col}'.")
    elif stats_model == "mixed_effect":
        mixed = _compute_panel_stats_mixed_effect(
            df=df,
            metric_col=metric_col,
            method_order=method_order,
            alpha=alpha,
            multiple_test_correction=multiple_test_correction,
            engine="pql",
            glmm_spec=glmm_spec,
            log=log,
        )
        if mixed is not None:
            return mixed
        if log is not None:
            log.info(f"Figure 09: falling back to paired_permutation stats for panel '{metric_col}'.")

    return _compute_panel_stats_permutation(
        df=df,
        metric_col=metric_col,
        method_order=method_order,
        n_perm=n_perm,
        seed=seed,
        alpha=alpha,
        pairwise_test=pairwise_test,
        multiple_test_correction=multiple_test_correction,
    )


def _draw_pairwise_brackets(
    ax,
    pairwise_rows: list[dict],
    method_order: list[str],
    alpha: float = 0.05,
    constrain_to_axis: bool = False,
    placement: str = "top",
    spacing_scale: float = 1.25,
) -> None:
    """
    Draw significance brackets between methods using the configured pairwise
    p-value (BH-adjusted or uncorrected).
    """
    if not pairwise_rows:
        return
    place_bottom = str(placement).strip().lower() == "bottom"
    spacing_scale = max(0.8, float(spacing_scale))

    method_to_x = {m: i for i, m in enumerate(method_order, start=1)}
    intervals: list[tuple[int, int, float]] = []
    for row in pairwise_rows:
        a = str(row.get("method_a", ""))
        b = str(row.get("method_b", ""))
        p = float(row.get("p_value_used", row.get("p_adj_bh", row.get("p_raw", np.nan))))
        if (a not in method_to_x) or (b not in method_to_x) or (not np.isfinite(p)):
            continue
        x1 = method_to_x[a]
        x2 = method_to_x[b]
        if x1 == x2:
            continue
        if x1 > x2:
            x1, x2 = x2, x1
        intervals.append((x1, x2, p))
    if not intervals:
        return

    intervals.sort(key=lambda t: ((t[1] - t[0]), t[2]))

    # Greedy level assignment to reduce overlaps.
    levels: list[list[tuple[int, int]]] = []
    placed: list[tuple[int, int, int, float]] = []  # x1, x2, level, p
    for x1, x2, p in intervals:
        placed_level = None
        for lv, spans in enumerate(levels):
            if all((x2 < a) or (x1 > b) for a, b in spans):
                spans.append((x1, x2))
                placed_level = lv
                break
        if placed_level is None:
            levels.append([(x1, x2)])
            placed_level = len(levels) - 1
        placed.append((x1, x2, placed_level, p))
    n_levels = max(1, len(levels))

    y_min, y_max = ax.get_ylim()
    y_range = float(y_max - y_min)
    if y_range <= 0:
        y_range = max(1.0, abs(y_max))

    if constrain_to_axis:
        # Keep brackets inside current y-limits (used for strict fraction panels).
        if place_bottom:
            band_low = y_min + (0.03 * y_range)
            band_high = y_min + (0.50 * y_range)
        else:
            band_low = y_min + (0.56 * y_range)
            band_high = y_min + (0.90 * y_range)
        band_span = max(1e-9, band_high - band_low)
        raw_step = (band_span / float(max(1, n_levels + 1))) * spacing_scale
        step = min(raw_step, band_span / float(max(1, n_levels - 1))) if n_levels > 1 else 0.0
        base = band_low
        cap_h = 0.012 * y_range
        star_off = 0.014 * y_range
    else:
        if place_bottom:
            base = y_min - (0.06 * y_range)
        else:
            base = y_max + (0.08 * y_range)
        step = 0.12 * y_range * spacing_scale
        cap_h = 0.026 * y_range
        star_off = 0.030 * y_range

    bracket_lw = 1.55
    star_fs = 10.2

    top = y_max
    bottom = y_min
    for x1, x2, level, p in placed:
        y = base + (level * step)
        if place_bottom:
            if constrain_to_axis:
                y = min(max(y, y_min + (0.08 * y_range)), y_max - (0.03 * y_range))
            y_cap = y - cap_h
            if constrain_to_axis:
                y_cap = max(y_cap, y_min + (0.05 * y_range))
            ax.plot([x1, x1, x2, x2], [y, y_cap, y_cap, y], color="#111111", linewidth=bracket_lw, zorder=5)
            star_y = y_cap - star_off
            if constrain_to_axis:
                star_y = max(star_y, y_min + (0.015 * y_range))
            ax.text(
                (x1 + x2) / 2.0,
                star_y,
                _pvalue_stars(p, alpha=alpha),
                ha="center",
                va="top",
                fontsize=star_fs,
                fontweight="bold",
                color="#111111",
                zorder=6,
            )
            bottom = min(bottom, star_y)
        else:
            if constrain_to_axis:
                y = min(max(y, y_min + (0.02 * y_range)), y_max - (0.10 * y_range))
            y_cap = y + cap_h
            if constrain_to_axis:
                y_cap = min(y_cap, y_max - (0.085 * y_range))
            ax.plot([x1, x1, x2, x2], [y, y_cap, y_cap, y], color="#111111", linewidth=bracket_lw, zorder=5)
            star_y = y_cap + star_off
            if constrain_to_axis:
                star_y = min(star_y, y_max - (0.045 * y_range))
            ax.text(
                (x1 + x2) / 2.0,
                star_y,
                _pvalue_stars(p, alpha=alpha),
                ha="center",
                va="bottom",
                fontsize=star_fs,
                fontweight="bold",
                color="#111111",
                zorder=6,
            )
            top = max(top, star_y)

    if not constrain_to_axis:
        if place_bottom:
            ax.set_ylim(bottom - (0.06 * y_range), y_max)
        else:
            ax.set_ylim(y_min, top + (0.09 * y_range))


def run(cfg) -> int:
    """
    Figure 09 (SVG):
      Split/Merge diagnostics by method using per-sample proxies:
        - split_index = duplicated assignments / selected denominator
          denominator selectable via figures.fig09.split_index_denominator:
            * mapped_expected_features (default)
            * mapped_inferred_features
          source selectable via figures.fig09.split_index_mode:
            * mapping        -> use finalized mapping.tsv (legacy behavior)
            * threshold_hits -> use all BLAST hits above method threshold
            * graph          -> BLAST-hit bipartite graph with method-specific thresholds
        - unmapped_feature_fraction = unmapped inferred features / inferred present features
        - one_to_one_mapping_fraction = inferred features in exact 1:1 relation / selected denominator
          denominator selectable via figures.fig09.one_to_one_denominator:
            * expected_present_features (default)
            * mapped_inferred_features
            * inferred_present_features
        - merge_deficit:
            * mapping/threshold_hits: 1 - species_coverage
            * graph: graph merge index = inferred fan-out burden / selected denominator
              denominator selectable via figures.fig09.merge_denominator:
                * mapped_inferred_features (default)
                * expected_present_features
        - cluster_purity_fraction (TIC/TAC):
            pure inferred clusters / inferred clusters with known zOTU members
            where a cluster is pure when all present members map and share one target.
        - no_hit_over_unmapped_fraction:
            no-hit inferred features / unmapped inferred features
            with no-hit evaluated using Figure 09 method-specific best-hit thresholds
            (fallback to Figure 12 thresholds when Figure 09-specific values are unset).
    """
    log = get_logger()
    outdir = Path(cfg.paths.outdir)
    exclude_patterns = list(getattr(cfg.figures, "exclude_community_name_patterns", []) or [])
    exclude_patterns += list(getattr(cfg.figures, "exclude_community_name_patterns_fig09", []) or [])
    split_index_mode = _normalize_split_index_mode(getattr(cfg.figures, "fig09_split_index_mode", "mapping"))
    # Figure 09 single-mode behavior: threshold_hits target is fixed to species.
    split_threshold_target_level = "species"
    split_index_denominator = _normalize_split_index_denominator(
        getattr(cfg.figures, "fig09_split_index_denominator", "mapped_expected_features")
    )
    one_to_one_denominator = _normalize_one_to_one_denominator(
        getattr(cfg.figures, "fig09_one_to_one_denominator", "expected_present_features")
    )
    merge_denominator = _normalize_merge_denominator(
        getattr(cfg.figures, "fig09_merge_denominator", "mapped_inferred_features")
    )
    fig09_min_qcov = float(getattr(cfg.figures, "fig09_min_query_coverage", 0.0))
    # Figure 09 can set a dedicated no-hit/unmapped query-coverage gate.
    # Fallback order: fig09.no_hit_min_query_coverage -> fig12.unmapped_best_hit_min_query_coverage -> fig09.min_query_coverage.
    unmapped_best_hit_min_qcov = getattr(cfg.figures, "fig09_no_hit_min_query_coverage", None)
    if unmapped_best_hit_min_qcov is None:
        unmapped_best_hit_min_qcov = getattr(cfg.figures, "fig12_unmapped_best_hit_min_query_coverage", None)
    if unmapped_best_hit_min_qcov is None:
        unmapped_best_hit_min_qcov = fig09_min_qcov
    else:
        unmapped_best_hit_min_qcov = float(unmapped_best_hit_min_qcov)
    no_hit_threshold_by_method = _fig09_no_hit_effective_identity_thresholds(cfg)
    transition_exact_identity = _as_percent_identity(
        float(getattr(cfg.figures, "fig09_transition_heatmap_exact_identity", 100.0))
    )
    raw_transition_near = getattr(cfg.figures, "fig09_transition_heatmap_near_exact_min_effective_pident", None)
    transition_near_min_effective_pident = (
        None if raw_transition_near is None else _as_percent_identity(float(raw_transition_near))
    )
    # Figure 09 is single-mode only for rendering, but in graph mode we compute
    # both tip and species rows to support panel-specific target-level selection.
    graph_dual_tip_species = bool(split_index_mode == "graph")
    graph_dual_combined = False

    mf_p = outdir / "manifest.tsv"
    if not mf_p.exists():
        log.info("Figure 09: manifest.tsv not found; run discover first.")
        return 0
    mf = pd.read_csv(mf_p, sep="\t")

    skipped = set()
    sp = qc_dir(outdir) / "skipped_communities.tsv"
    if sp.exists():
        try:
            skipped = set(pd.read_csv(sp, sep="\t")["community_id"].astype(str).tolist())
        except Exception:
            skipped = set()

    rows: list[dict] = []
    unmapped_best_hit_rows: list[dict] = []
    denoising_best_hit_rows: list[dict] = []
    transition_heatmap_rows: list[dict] = []

    def _record_unmapped_best_hits(
        community_id: str,
        sample_id: str,
        sample_key: str,
        method: str,
        target_level: str,
        unmapped_features: list[str],
        best_pident_by_feature: dict[str, float],
        best_effective_pident_by_feature: dict[str, float],
        feature_fraction_by_id: dict[str, float],
    ) -> None:
        for feat in unmapped_features:
            best_eff = float(best_effective_pident_by_feature.get(str(feat), np.nan))
            best_pid = float(best_pident_by_feature.get(str(feat), np.nan))
            has_best = bool(np.isfinite(best_eff) or np.isfinite(best_pid))
            frac = float(feature_fraction_by_id.get(str(feat), np.nan))
            unmapped_best_hit_rows.append(
                {
                    "community_id": community_id,
                    "sample_id": sample_id,
                    "sample_key": sample_key,
                    "method": method,
                    "graph_target_level_plot": target_level,
                    "unmapped_feature": str(feat),
                    "best_hit_effective_pident": best_eff,
                    "best_hit_pident": best_pid,
                    "best_hit_status": "hit_found" if has_best else "no_hit_detected",
                    "inferred_feature_fraction": frac,
                }
            )
    for _, r in mf.iterrows():
        cid = str(r["community_id"])
        if cid in skipped:
            continue
        if _is_excluded(cid, exclude_patterns):
            continue

        ed = expected_dir(outdir, cid)
        exp_p = ed / "expected_species_abundance_tss1000.tsv"
        if not exp_p.exists():
            continue
        exp = pd.read_csv(exp_p, sep="\t", index_col=0)
        exp_tip = pd.DataFrame()
        exp_tip_p = ed / "expected_tip_abundance_tss1000.tsv"
        if exp_tip_p.exists():
            try:
                exp_tip = pd.read_csv(exp_tip_p, sep="\t", index_col=0)
            except Exception:
                exp_tip = pd.DataFrame()
        samples = [str(s) for s in exp.columns]
        tip_to_species = _load_tip_to_species(ed / "tip_to_species.tsv")

        # TIC/TAC cluster-purity support:
        # - membership is read from raw TACTIC map files (cluster -> zOTU members)
        # - member hits/presence are taken from zOTU outputs in the benchmark outdir
        zotu_std = pd.DataFrame()
        zotu_std_p = method_dir(outdir, cid, "zOTU") / "table_std_tss1000.tsv"
        if zotu_std_p.exists():
            try:
                zotu_std = pd.read_csv(zotu_std_p, sep="\t", index_col=0)
                zotu_std.index = zotu_std.index.astype(str)
            except Exception:
                zotu_std = pd.DataFrame()
        zotu_hits_all = pd.DataFrame()
        zotu_hits_p = method_dir(outdir, cid, "zOTU") / "blast_hits.tsv"
        if zotu_hits_p.exists() and zotu_hits_p.stat().st_size > 0:
            try:
                zotu_hits_all = read_hits(zotu_hits_p)
            except Exception:
                zotu_hits_all = pd.DataFrame()
        tic_cluster_members = _load_cluster_membership_map(_manifest_path(r.get("TIC_table")), method="TIC")
        tac_cluster_members = _load_cluster_membership_map(_manifest_path(r.get("TAC_table")), method="TAC")
        tic_member_to_cluster = _invert_cluster_membership_map(tic_cluster_members)
        tac_member_to_cluster = _invert_cluster_membership_map(tac_cluster_members)
        zotu_transition_best_eff = _best_effective_identity_by_feature_from_hits(
            zotu_hits_all,
            min_qcov=unmapped_best_hit_min_qcov,
        )
        zotu_transition_best_info = _best_hit_info_by_feature_from_hits(
            zotu_hits_all,
            min_qcov=unmapped_best_hit_min_qcov,
        )
        zotu_transition_threshold = float(
            transition_near_min_effective_pident
            if transition_near_min_effective_pident is not None
            else no_hit_threshold_by_method.get("zOTU", 97.0)
        )

        for method in METHOD_ORDER:
            md = method_dir(outdir, cid, method)
            std_p = md / "table_std_tss1000.tsv"
            if not std_p.exists():
                continue

            table_std = pd.read_csv(std_p, sep="\t", index_col=0)
            table_std.index = table_std.index.astype(str)
            inferred_feature_ids = table_std.index.astype(str).tolist()
            strict_mapped_features: set[str] = set()
            mapping_p = md / "mapping.tsv"
            if mapping_p.exists():
                try:
                    mp = pd.read_csv(mapping_p, sep="\t")
                    if "inferred_feature" in mp.columns:
                        strict_mapped_features = set(mp["inferred_feature"].astype(str).tolist())
                except Exception:
                    strict_mapped_features = set()

            # Exploratory companion plot: for unmapped inferred features, keep
            # the best available BLAST identity after query-coverage gating.
            best_pident_by_feature: dict[str, float] = {}
            best_effective_pident_by_feature: dict[str, float] = {}
            denoising_best_hit_by_feature: dict[str, dict] = {}
            hits_best_path = md / "blast_hits.tsv"
            if hits_best_path.exists() and hits_best_path.stat().st_size > 0:
                try:
                    hits_best_raw = read_hits(hits_best_path)
                    hits_best_all = hits_best_raw.copy()
                    hits_best = filter_hits(
                        hits_best_raw.copy(),
                        min_effective_pident=None,
                        min_qcovhsp=unmapped_best_hit_min_qcov,
                        strict_pident_gt=False,
                    )
                    if method in {"ASV", "zOTU"} and not hits_best_all.empty:
                        hb = hits_best_all.copy()
                        hb["qseqid"] = hb["qseqid"].astype(str)
                        hb["sseqid"] = hb["sseqid"].astype(str)
                        hb["effective_identity"] = pd.to_numeric(
                            hb.get("effective_identity", np.nan),
                            errors="coerce",
                        )
                        hb["pident"] = pd.to_numeric(hb["pident"], errors="coerce")
                        hb["bitscore"] = pd.to_numeric(hb["bitscore"], errors="coerce")
                        hb["evalue"] = pd.to_numeric(hb["evalue"], errors="coerce")
                        hb["length"] = pd.to_numeric(hb["length"], errors="coerce")
                        hb["qcov_exact"] = pd.to_numeric(hb.get("qcov_exact", np.nan), errors="coerce")
                        hb = hb.dropna(subset=["qseqid", "effective_identity", "pident"])
                        if not hb.empty:
                            hb_best = hb.sort_values(
                                ["qseqid", "effective_identity", "pident", "bitscore", "evalue", "length"],
                                ascending=[True, False, False, False, True, False],
                            ).groupby("qseqid", sort=False).first()
                            denoising_best_hit_by_feature = {
                                str(q): {
                                    "best_target": str(r.get("sseqid", "")),
                                    "best_pident": float(r.get("pident", np.nan)),
                                    "best_qcov_exact_percent": float(r.get("qcov_exact", np.nan)) * 100.0
                                    if np.isfinite(float(r.get("qcov_exact", np.nan)))
                                    else np.nan,
                                    "best_effective_identity": float(r.get("effective_identity", np.nan)),
                                }
                                for q, r in hb_best.iterrows()
                            }
                    if not hits_best.empty:
                        hits_best["qseqid"] = hits_best["qseqid"].astype(str)
                        hits_best["effective_identity"] = pd.to_numeric(
                            hits_best.get("effective_identity", np.nan),
                            errors="coerce",
                        )
                        hits_best["pident"] = pd.to_numeric(hits_best["pident"], errors="coerce")
                        hits_best["bitscore"] = pd.to_numeric(hits_best["bitscore"], errors="coerce")
                        hits_best["evalue"] = pd.to_numeric(hits_best["evalue"], errors="coerce")
                        hits_best["length"] = pd.to_numeric(hits_best["length"], errors="coerce")
                        hits_best = hits_best.dropna(subset=["qseqid", "effective_identity", "pident"])
                        if not hits_best.empty:
                            best_hits = hits_best.sort_values(
                                ["qseqid", "effective_identity", "pident", "bitscore", "evalue", "length"],
                                ascending=[True, False, False, False, True, False],
                            ).groupby("qseqid", sort=False).first()
                            best_effective_pident_by_feature = {
                                str(q): float(p)
                                for q, p in zip(
                                    best_hits.index.astype(str),
                                    best_hits["effective_identity"].astype(float),
                                )
                            }
                            best_pident_by_feature = {
                                str(q): float(p)
                                for q, p in zip(best_hits.index.astype(str), best_hits["pident"].astype(float))
                            }
                except Exception as e:
                    raise RuntimeError(
                        f"[{cid}] {method}: strict BLAST parsing failed for unmapped best-hit diagnostics: {e}"
                    ) from e
            if method in {"ASV", "zOTU"}:
                method_label = display_methods([method])[0]
                for feat in inferred_feature_ids:
                    rec = denoising_best_hit_by_feature.get(str(feat), {})
                    denoising_best_hit_rows.append(
                        {
                            "community_id": cid,
                            "method": method,
                            "method_label": method_label,
                            "inferred_feature": str(feat),
                            "best_target": str(rec.get("best_target", "")),
                            "best_pident": float(rec.get("best_pident", np.nan)),
                            "best_qcov_exact_percent": float(rec.get("best_qcov_exact_percent", np.nan)),
                            "best_effective_identity": float(rec.get("best_effective_identity", np.nan)),
                            "accepted_in_strict_mapping": bool(str(feat) in strict_mapped_features),
                        }
                    )

            threshold_targets_by_feature: dict[str, set[str]] = {}
            graph_targets_by_metric: dict[str, dict[str, set[str]]] = {}
            graph_targets_by_level: dict[str, dict[str, set[str]]] = {}
            graph_targets_by_metric_level: dict[str, dict[str, dict[str, set[str]]]] = {}
            graph_min_ident_used = np.nan
            graph_min_ident_by_metric: dict[str, float] = {}
            graph_target_levels = {
                "split_index": _graph_target_level_for_method(cfg, method, "split_index"),
                "merge_deficit": _graph_target_level_for_method(cfg, method, "merge_deficit"),
                "unmapped_feature_fraction": _graph_target_level_for_method(cfg, method, "unmapped_feature_fraction"),
                "one_to_one_mapping_fraction": _graph_target_level_for_method(cfg, method, "one_to_one_mapping_fraction"),
                "cluster_purity_fraction": _graph_target_level_for_method(cfg, method, "cluster_purity_fraction"),
                "no_hit_over_unmapped_fraction": _graph_target_level_for_method(
                    cfg, method, "no_hit_over_unmapped_fraction"
                ),
            }
            cluster_member_targets_cache: dict[tuple[str, float], dict[str, set[str]]] = {}

            if split_index_mode in {"threshold_hits", "graph"}:
                hits_p = md / "blast_hits.tsv"
                if hits_p.exists() and hits_p.stat().st_size > 0:
                    try:
                        if split_index_mode == "graph":
                            hits_all = read_hits(hits_p)
                            for metric_key, chosen_level in graph_target_levels.items():
                                min_ident = _method_graph_min_effective_pident(cfg, method, metric_key)
                                graph_min_ident_by_metric[metric_key] = float(min_ident)
                                hits_metric = filter_hits(
                                    hits_all.copy(),
                                    min_effective_pident=min_ident,
                                    min_qcovhsp=fig09_min_qcov,
                                    strict_pident_gt=False,
                                )
                                metric_level_targets: dict[str, dict[str, set[str]]] = {}
                                if not hits_metric.empty:
                                    hits_metric["qseqid"] = hits_metric["qseqid"].astype(str)
                                    hits_metric["sseqid"] = hits_metric["sseqid"].astype(str)
                                    for level in GRAPH_LEVEL_ORDER:
                                        h2 = hits_metric.copy()
                                        if level == "species":
                                            h2["target"] = h2["sseqid"].map(
                                                lambda x: tip_to_species.get(str(x), genus_species(str(x)))
                                            )
                                        else:
                                            h2["target"] = h2["sseqid"]
                                        hit_pairs = h2[["qseqid", "target"]].dropna().drop_duplicates()
                                        tmap: dict[str, set[str]] = {}
                                        for q, g in hit_pairs.groupby("qseqid", sort=False):
                                            tmap[str(q)] = set(g["target"].astype(str).tolist())
                                        metric_level_targets[level] = tmap
                                else:
                                    for level in GRAPH_LEVEL_ORDER:
                                        metric_level_targets[level] = {}
                                graph_targets_by_metric_level[metric_key] = metric_level_targets
                                graph_targets_by_metric[metric_key] = metric_level_targets.get(chosen_level, {})
                            graph_targets_by_level = graph_targets_by_metric_level.get("split_index", {})
                            graph_min_ident_used = float(graph_min_ident_by_metric.get("split_index", np.nan))
                        else:
                            hits = read_hits(hits_p)
                            min_ident = _method_min_effective_pident(cfg, method)
                            hits = filter_hits(
                                hits,
                                min_effective_pident=min_ident,
                                min_qcovhsp=fig09_min_qcov,
                                # Use inclusive identity thresholds (>=) for all methods.
                                strict_pident_gt=False,
                            )
                            if not hits.empty:
                                hits["qseqid"] = hits["qseqid"].astype(str)
                                hits["sseqid"] = hits["sseqid"].astype(str)
                                if split_threshold_target_level == "species":
                                    hits["target"] = hits["sseqid"].map(
                                        lambda x: tip_to_species.get(str(x), genus_species(str(x)))
                                    )
                                else:
                                    hits["target"] = hits["sseqid"]

                                hit_pairs = hits[["qseqid", "target"]].dropna().drop_duplicates()
                                for q, g in hit_pairs.groupby("qseqid", sort=False):
                                    threshold_targets_by_feature[str(q)] = set(g["target"].astype(str).tolist())
                    except Exception as e:
                        raise RuntimeError(
                            f"[{cid}] {method}: strict BLAST parsing failed while building Figure 09 hit maps: {e}"
                        ) from e

            feat2sp: dict[str, str] = {}
            mapped_sp: pd.DataFrame | None = None

            if split_index_mode == "mapping":
                map_p = md / "mapping.tsv"
                if not map_p.exists():
                    continue
                mapping = pd.read_csv(map_p, sep="\t")
                feat2sp = _feature_to_species_map(mapping, tip_to_species)

                if method in CLUSTERING:
                    msp_p = md / "mapped_species_table_tss1000.tsv"
                    if not msp_p.exists():
                        continue
                    mapped_sp = pd.read_csv(msp_p, sep="\t", index_col=0)
                else:
                    mtip_p = md / "mapped_tip_table_tss1000.tsv"
                    if not mtip_p.exists():
                        continue
                    mapped_tip = pd.read_csv(mtip_p, sep="\t", index_col=0)
                    mapped_sp = _collapse_tip_table_to_species(mapped_tip, tip_to_species)

            elif split_index_mode == "threshold_hits":
                # merge_deficit remains the legacy species coverage proxy in this mode.
                if method in CLUSTERING:
                    msp_p = md / "mapped_species_table_tss1000.tsv"
                    if not msp_p.exists():
                        continue
                    mapped_sp = pd.read_csv(msp_p, sep="\t", index_col=0)
                else:
                    mtip_p = md / "mapped_tip_table_tss1000.tsv"
                    if not mtip_p.exists():
                        continue
                    mapped_tip = pd.read_csv(mtip_p, sep="\t", index_col=0)
                    mapped_sp = _collapse_tip_table_to_species(mapped_tip, tip_to_species)

            for s in samples:
                if s not in table_std.columns:
                    continue

                inferred_present = set(table_std.index[table_std[s].astype(float) > 0].astype(str))
                exp_present = set(exp.index[exp[s].astype(float) > 0].astype(str))
                expected_present_features = int(len(exp_present))
                n_inf = len(inferred_present)
                if n_inf == 0:
                    continue
                sample_vals = pd.to_numeric(table_std[s], errors="coerce").fillna(0.0).clip(lower=0.0)
                sample_total = float(sample_vals.sum())
                if sample_total > 0.0:
                    feature_fraction_by_id = {
                        str(fid): float(val / sample_total)
                        for fid, val in sample_vals.items()
                        if float(val) > 0.0
                    }
                else:
                    feature_fraction_by_id = {}

                if (
                    method in {"TIC", "TAC"}
                    and (not zotu_std.empty)
                    and (s in zotu_std.columns)
                ):
                    cluster_to_members = tic_cluster_members if method == "TIC" else tac_cluster_members
                    member_to_cluster = tic_member_to_cluster if method == "TIC" else tac_member_to_cluster
                    if cluster_to_members:
                        zotu_vals = pd.to_numeric(zotu_std[s], errors="coerce").fillna(0.0)
                        zotu_present = set(zotu_vals.index[zotu_vals.astype(float) > 0].astype(str).tolist())
                        zotu_class_by_feature = {
                            str(zotu_id): _fig09_transition_classify_feature(
                                float(zotu_transition_best_eff.get(str(zotu_id), np.nan)),
                                near_exact_threshold=zotu_transition_threshold,
                                exact_identity=transition_exact_identity,
                            )
                            for zotu_id in zotu_present
                        }
                        cluster_class_by_id: dict[str, str] = {}
                        for cluster_id, members in cluster_to_members.items():
                            present_members = [
                                str(member_id)
                                for member_id in members
                                if str(member_id) in zotu_class_by_feature
                            ]
                            if present_members:
                                cluster_class_by_id[str(cluster_id)] = _fig09_transition_pick_priority_class(
                                    [zotu_class_by_feature[str(member_id)] for member_id in present_members]
                                )

                        for zotu_id in sorted(zotu_present):
                            cluster_id = str(member_to_cluster.get(str(zotu_id), ""))
                            if (not cluster_id) or (cluster_id not in inferred_present):
                                final_fate = "Unassigned"
                                classification_source = "zOTU member not retained in present cluster set"
                            else:
                                final_fate = cluster_class_by_id.get(cluster_id, "No hit")
                                classification_source = "present cluster priority class"
                            best_info = zotu_transition_best_info.get(str(zotu_id), {})
                            best_target_id = str(best_info.get("target_id", ""))
                            transition_heatmap_rows.append(
                                {
                                    "community_id": cid,
                                    "sample_id": s,
                                    "sample_key": f"{cid}|{s}",
                                    "transition_method": method,
                                    "zotu_id": str(zotu_id),
                                    "cluster_id": cluster_id,
                                    "zotu_class": zotu_class_by_feature.get(str(zotu_id), "No hit"),
                                    "final_fate": final_fate,
                                    "best_hit_effective_identity": float(
                                        zotu_transition_best_eff.get(str(zotu_id), np.nan)
                                    ),
                                    "best_hit_target_id": best_target_id,
                                    "best_hit_species": str(tip_to_species.get(best_target_id, "")),
                                    "zotu_feature_value_tss1000": float(zotu_vals.get(str(zotu_id), 0.0)),
                                    "near_exact_threshold": zotu_transition_threshold,
                                    "exact_identity_threshold": transition_exact_identity,
                                    "classification_source": classification_source,
                                }
                            )

                split_duplicates = 0
                split_index = np.nan
                species_coverage = np.nan
                one_to_one_mapping_fraction = np.nan
                one_to_one_denominator_features = 0
                merge_deficit = np.nan
                split_mapped_n = 0
                merge_mapped_n = 0
                split_denominator_features = 0
                graph_edge_count = np.nan
                graph_merge_duplicates = np.nan
                graph_target_level_plot = np.nan

                if split_index_mode == "mapping":
                    mapped_present = [f for f in inferred_present if f in feat2sp]
                    unmapped_present_feats = [f for f in inferred_present if f not in feat2sp]
                    n_mapped = len(mapped_present)
                    n_unmapped = n_inf - n_mapped
                    split_mapped_n = n_mapped

                    species_counts: dict[str, int] = {}
                    for f in mapped_present:
                        sp_name = str(feat2sp.get(f, ""))
                        if sp_name:
                            species_counts[sp_name] = species_counts.get(sp_name, 0) + 1
                    split_duplicates = sum(max(0, c - 1) for c in species_counts.values())
                    split_denominator_features = (
                        len(species_counts) if split_index_denominator == "mapped_expected_features" else split_mapped_n
                    )
                    split_index = (
                        float(split_duplicates / split_denominator_features)
                        if split_denominator_features > 0
                        else np.nan
                    )
                    merge_denominator_features = np.nan
                    unmapped_fraction = float(n_unmapped / n_inf)
                    _record_unmapped_best_hits(
                        community_id=cid,
                        sample_id=s,
                        sample_key=f"{cid}|{s}",
                        method=method,
                        target_level="species",
                        unmapped_features=unmapped_present_feats,
                        best_pident_by_feature=best_pident_by_feature,
                        best_effective_pident_by_feature=best_effective_pident_by_feature,
                        feature_fraction_by_id=feature_fraction_by_id,
                    )

                    if isinstance(mapped_sp, pd.DataFrame) and (s in mapped_sp.columns):
                        rec_present = set(mapped_sp.index[mapped_sp[s].astype(float) > 0].astype(str))
                    else:
                        rec_present = set()
                    if exp_present:
                        species_coverage = float(len(exp_present & rec_present) / len(exp_present))
                        merge_deficit = 1.0 - species_coverage

                    one_to_one_n = 0
                    for f in mapped_present:
                        sp_name = str(feat2sp.get(f, ""))
                        if sp_name and species_counts.get(sp_name, 0) == 1:
                            one_to_one_n += 1
                    if one_to_one_denominator == "expected_present_features":
                        one_to_one_denominator_features = len(exp_present)
                    elif one_to_one_denominator == "inferred_present_features":
                        one_to_one_denominator_features = n_inf
                    else:
                        one_to_one_denominator_features = split_mapped_n
                    one_to_one_mapping_fraction = (
                        float(one_to_one_n / one_to_one_denominator_features)
                        if one_to_one_denominator_features > 0
                        else np.nan
                    )

                elif split_index_mode == "threshold_hits":
                    split_mapped_feats = [f for f in inferred_present if f in threshold_targets_by_feature]
                    unmapped_present_feats = [f for f in inferred_present if f not in threshold_targets_by_feature]
                    split_mapped_n = len(split_mapped_feats)
                    n_mapped = split_mapped_n
                    n_unmapped = n_inf - n_mapped
                    unmapped_fraction = float(n_unmapped / n_inf)
                    _record_unmapped_best_hits(
                        community_id=cid,
                        sample_id=s,
                        sample_key=f"{cid}|{s}",
                        method=method,
                        target_level=split_threshold_target_level,
                        unmapped_features=unmapped_present_feats,
                        best_pident_by_feature=best_pident_by_feature,
                        best_effective_pident_by_feature=best_effective_pident_by_feature,
                        feature_fraction_by_id=feature_fraction_by_id,
                    )

                    target_counts: dict[str, int] = {}
                    edge_count = 0
                    for f in split_mapped_feats:
                        targets = threshold_targets_by_feature.get(f, set())
                        edge_count += len(targets)
                        for tgt in targets:
                            target_counts[tgt] = target_counts.get(tgt, 0) + 1
                    split_duplicates = sum(max(0, c - 1) for c in target_counts.values())
                    split_denominator_features = (
                        len(target_counts) if split_index_denominator == "mapped_expected_features" else split_mapped_n
                    )
                    split_index = (
                        float(split_duplicates / split_denominator_features)
                        if split_denominator_features > 0
                        else np.nan
                    )
                    merge_denominator_features = np.nan

                    if isinstance(mapped_sp, pd.DataFrame) and (s in mapped_sp.columns):
                        rec_present = set(mapped_sp.index[mapped_sp[s].astype(float) > 0].astype(str))
                    else:
                        rec_present = set()
                    if exp_present:
                        species_coverage = float(len(exp_present & rec_present) / len(exp_present))
                        merge_deficit = 1.0 - species_coverage

                    one_to_one_n = 0
                    for f in split_mapped_feats:
                        targets = threshold_targets_by_feature.get(f, set())
                        if len(targets) == 1:
                            tgt = next(iter(targets))
                            if target_counts.get(tgt, 0) == 1:
                                one_to_one_n += 1
                    if one_to_one_denominator == "expected_present_features":
                        if split_threshold_target_level == "tip":
                            if (not exp_tip.empty) and (s in exp_tip.columns):
                                one_to_one_denominator_features = int((exp_tip[s].astype(float) > 0).sum())
                            else:
                                one_to_one_denominator_features = 0
                        else:
                            one_to_one_denominator_features = len(exp_present)
                    elif one_to_one_denominator == "inferred_present_features":
                        one_to_one_denominator_features = n_inf
                    else:
                        one_to_one_denominator_features = split_mapped_n
                    one_to_one_mapping_fraction = (
                        float(one_to_one_n / one_to_one_denominator_features)
                        if one_to_one_denominator_features > 0
                        else np.nan
                    )

                    graph_edge_count = float(edge_count)

                else:
                    if graph_dual_tip_species:
                        exp_present = set(exp.index[exp[s].astype(float) > 0].astype(str))
                        dual_rows: list[dict] = []
                        for dual_level in GRAPH_LEVEL_ORDER:
                            split_targets_level = (
                                graph_targets_by_metric_level.get("split_index", {}).get(dual_level, {})
                            )
                            merge_targets_level = (
                                graph_targets_by_metric_level.get("merge_deficit", {}).get(dual_level, {})
                            )
                            unmapped_targets_level = (
                                graph_targets_by_metric_level.get("unmapped_feature_fraction", {}).get(dual_level, {})
                            )
                            one_to_one_targets_level = (
                                graph_targets_by_metric_level.get("one_to_one_mapping_fraction", {}).get(dual_level, {})
                            )
                            no_hit_targets_level = (
                                graph_targets_by_metric_level.get("no_hit_over_unmapped_fraction", {}).get(
                                    dual_level, {}
                                )
                            )

                            split_mapped_feats = [f for f in inferred_present if f in split_targets_level]
                            split_mapped_n = len(split_mapped_feats)
                            merge_mapped_feats = [f for f in inferred_present if f in merge_targets_level]
                            merge_mapped_n = len(merge_mapped_feats)
                            unmapped_mapped_feats = [f for f in inferred_present if f in unmapped_targets_level]
                            mapped_n = len(unmapped_mapped_feats)
                            unmapped_feats_level = [f for f in inferred_present if f not in unmapped_targets_level]
                            n_unmapped = n_inf - mapped_n
                            unmapped_fraction = float(n_unmapped / n_inf)
                            no_hit_unmapped_feats = [f for f in inferred_present if f not in no_hit_targets_level]
                            no_hit_n_unmapped = int(len(no_hit_unmapped_feats))
                            no_hit_n_no_hit = int(
                                sum(
                                    1
                                    for f in no_hit_unmapped_feats
                                    if not _feature_has_best_hit(
                                        method=method,
                                        feature_id=str(f),
                                        best_effective_pident_by_feature=best_effective_pident_by_feature,
                                        best_pident_by_feature=best_pident_by_feature,
                                        no_hit_threshold_by_method=no_hit_threshold_by_method,
                                    )
                                )
                            )
                            no_hit_over_unmapped_fraction = (
                                float(no_hit_n_no_hit / no_hit_n_unmapped) if no_hit_n_unmapped > 0 else 0.0
                            )
                            _record_unmapped_best_hits(
                                community_id=cid,
                                sample_id=s,
                                sample_key=f"{cid}|{s}",
                                method=method,
                                target_level=dual_level,
                                unmapped_features=unmapped_feats_level,
                                best_pident_by_feature=best_pident_by_feature,
                                best_effective_pident_by_feature=best_effective_pident_by_feature,
                                feature_fraction_by_id=feature_fraction_by_id,
                            )

                            target_counts: dict[str, int] = {}
                            edge_count = 0
                            for f in split_mapped_feats:
                                targets = split_targets_level.get(f, set())
                                edge_count += len(targets)
                                for tgt in targets:
                                    target_counts[tgt] = target_counts.get(tgt, 0) + 1

                            split_duplicates_level = sum(max(0, c - 1) for c in target_counts.values())
                            split_denom_level = (
                                len(target_counts)
                                if split_index_denominator == "mapped_expected_features"
                                else split_mapped_n
                            )
                            split_index_level = (
                                float(split_duplicates_level / split_denom_level)
                                if split_denom_level > 0
                                else np.nan
                            )

                            merge_duplicates_level = 0
                            for f in merge_mapped_feats:
                                deg = len(merge_targets_level.get(f, set()))
                                merge_duplicates_level += max(0, deg - 1)
                            if merge_denominator == "expected_present_features":
                                merge_den_level = len(exp_present)
                            else:
                                merge_den_level = merge_mapped_n
                            merge_index_level = (
                                float(merge_duplicates_level / merge_den_level) if merge_den_level > 0 else np.nan
                            )

                            one_to_one_target_counts: dict[str, int] = {}
                            one_to_one_mapped_feats = [f for f in inferred_present if f in one_to_one_targets_level]
                            one_to_one_mapped_n = len(one_to_one_mapped_feats)
                            for f in one_to_one_mapped_feats:
                                for tgt in one_to_one_targets_level.get(f, set()):
                                    one_to_one_target_counts[tgt] = one_to_one_target_counts.get(tgt, 0) + 1

                            one_to_one_n = 0
                            for f in one_to_one_mapped_feats:
                                targets = one_to_one_targets_level.get(f, set())
                                if len(targets) == 1:
                                    tgt = next(iter(targets))
                                    if one_to_one_target_counts.get(tgt, 0) == 1:
                                        one_to_one_n += 1

                            if one_to_one_denominator == "expected_present_features":
                                if dual_level == "tip":
                                    if (not exp_tip.empty) and (s in exp_tip.columns):
                                        one_to_one_den_level = int((exp_tip[s].astype(float) > 0).sum())
                                    else:
                                        one_to_one_den_level = 0
                                else:
                                    one_to_one_den_level = len(exp_present)
                            elif one_to_one_denominator == "inferred_present_features":
                                one_to_one_den_level = n_inf
                            else:
                                one_to_one_den_level = one_to_one_mapped_n

                            one_to_one_level = (
                                float(one_to_one_n / one_to_one_den_level)
                                if one_to_one_den_level > 0
                                else np.nan
                            )

                            purity_fraction = np.nan
                            purity_score = np.nan
                            purity_pure_clusters = 0.0
                            purity_total_clusters = 0.0
                            if (
                                method in PURITY_METHODS
                                and (not zotu_std.empty)
                                and (s in zotu_std.columns)
                                and (not zotu_hits_all.empty)
                            ):
                                cluster_to_members = tic_cluster_members if method == "TIC" else tac_cluster_members
                                if cluster_to_members:
                                    purity_min_ident = float(
                                        graph_min_ident_by_metric.get("cluster_purity_fraction", np.nan)
                                    )
                                    cache_key = (str(dual_level), float(purity_min_ident))
                                    member_targets = cluster_member_targets_cache.get(cache_key)
                                    if member_targets is None:
                                        if np.isfinite(purity_min_ident):
                                            member_targets = _build_member_target_map(
                                                hits_all=zotu_hits_all,
                                                min_ident=purity_min_ident,
                                                min_qcov=fig09_min_qcov,
                                                target_level=dual_level,
                                                tip_to_species=tip_to_species,
                                            )
                                        else:
                                            member_targets = {}
                                        cluster_member_targets_cache[cache_key] = member_targets
                                    zotu_present = set(
                                        zotu_std.index[zotu_std[s].astype(float) > 0].astype(str).tolist()
                                    )
                                    purity_stats = _compute_cluster_purity_summary(
                                        present_clusters=set(str(f) for f in inferred_present),
                                        cluster_to_members=cluster_to_members,
                                        present_members=zotu_present,
                                        member_targets=member_targets,
                                    )
                                    purity_fraction = float(purity_stats["cluster_purity_fraction"])
                                    purity_score = float(purity_stats["cluster_purity_score"])
                                    purity_pure_clusters = float(purity_stats["cluster_purity_pure_clusters"])
                                    purity_total_clusters = float(purity_stats["cluster_purity_total_clusters"])

                            dual_rows.append(
                                {
                                    "level": dual_level,
                                    "mapped_n": mapped_n,
                                    "unmapped_n": n_unmapped,
                                    "split_mapped_n": split_mapped_n,
                                    "merge_mapped_n": merge_mapped_n,
                                    "merge_denominator_features": merge_den_level,
                                    "split_duplicates": split_duplicates_level,
                                    "split_index": split_index_level,
                                    "split_denominator_features": split_denom_level,
                                    "merge_deficit": merge_index_level,
                                    "graph_edge_count": float(edge_count),
                                    "graph_merge_duplicates": float(merge_duplicates_level),
                                    "unmapped_fraction": unmapped_fraction,
                                    "one_to_one_denominator_features": one_to_one_den_level,
                                    "one_to_one_mapping_fraction": one_to_one_level,
                                    "graph_edge_threshold_split_index": float(
                                        graph_min_ident_by_metric.get("split_index", np.nan)
                                    ),
                                    "graph_edge_threshold_merge_deficit": float(
                                        graph_min_ident_by_metric.get("merge_deficit", np.nan)
                                    ),
                                    "graph_edge_threshold_unmapped_feature_fraction": float(
                                        graph_min_ident_by_metric.get("unmapped_feature_fraction", np.nan)
                                    ),
                                    "graph_edge_threshold_one_to_one_mapping_fraction": float(
                                        graph_min_ident_by_metric.get("one_to_one_mapping_fraction", np.nan)
                                    ),
                                    "graph_edge_threshold_cluster_purity_fraction": float(
                                        graph_min_ident_by_metric.get("cluster_purity_fraction", np.nan)
                                    ),
                                    "graph_edge_threshold_no_hit_over_unmapped_fraction": float(
                                        graph_min_ident_by_metric.get("no_hit_over_unmapped_fraction", np.nan)
                                    ),
                                    "no_hit_over_unmapped_fraction": no_hit_over_unmapped_fraction,
                                    "no_hit_n_no_hit": float(no_hit_n_no_hit),
                                    "no_hit_n_unmapped": float(no_hit_n_unmapped),
                                    "cluster_purity_fraction": purity_fraction,
                                    "cluster_purity_score": purity_score,
                                    "cluster_purity_pure_clusters": purity_pure_clusters,
                                    "cluster_purity_total_clusters": purity_total_clusters,
                                }
                            )

                        for dual in dual_rows:
                            rows.append(
                                {
                                    "community_id": cid,
                                    "sample_id": s,
                                    "sample_key": f"{cid}|{s}",
                                    "method": method,
                                    "inferred_present_features": n_inf,
                                    "expected_present_features": expected_present_features,
                                    "mapped_present_features": dual["mapped_n"],
                                    "unmapped_present_features": dual["unmapped_n"],
                                    "split_mapped_present_features": dual["split_mapped_n"],
                                    "merge_mapped_present_features": dual["merge_mapped_n"],
                                    "merge_denominator": merge_denominator,
                                    "merge_denominator_features": dual["merge_denominator_features"],
                                    "split_index_mode": split_index_mode,
                                    "split_index_threshold_target_level": dual["level"],
                                    "graph_target_level_split_index": dual["level"],
                                    "graph_target_level_merge_deficit": dual["level"],
                                    "graph_target_level_unmapped_feature_fraction": dual["level"],
                                    "graph_target_level_one_to_one_mapping_fraction": dual["level"],
                                    "graph_target_level_cluster_purity_fraction": dual["level"],
                                    "graph_target_level_no_hit_over_unmapped_fraction": dual["level"],
                                    "graph_target_level_plot": dual["level"],
                                    "split_index_denominator": split_index_denominator,
                                    "split_index_denominator_features": dual["split_denominator_features"],
                                    "graph_edge_threshold_pident": dual["graph_edge_threshold_split_index"],
                                    "graph_edge_threshold_split_index": dual["graph_edge_threshold_split_index"],
                                    "graph_edge_threshold_merge_deficit": dual["graph_edge_threshold_merge_deficit"],
                                    "graph_edge_threshold_unmapped_feature_fraction": dual[
                                        "graph_edge_threshold_unmapped_feature_fraction"
                                    ],
                                    "graph_edge_threshold_one_to_one_mapping_fraction": dual[
                                        "graph_edge_threshold_one_to_one_mapping_fraction"
                                    ],
                                    "graph_edge_threshold_cluster_purity_fraction": dual[
                                        "graph_edge_threshold_cluster_purity_fraction"
                                    ],
                                    "graph_edge_threshold_no_hit_over_unmapped_fraction": dual[
                                        "graph_edge_threshold_no_hit_over_unmapped_fraction"
                                    ],
                                    "graph_edge_min_query_coverage": fig09_min_qcov,
                                    "graph_edge_count": dual["graph_edge_count"],
                                    "graph_merge_duplicates": dual["graph_merge_duplicates"],
                                    "split_duplicates": dual["split_duplicates"],
                                    "split_index": dual["split_index"],
                                    "unmapped_feature_fraction": dual["unmapped_fraction"],
                                    "one_to_one_mapping_fraction": dual["one_to_one_mapping_fraction"],
                                    "one_to_one_denominator": one_to_one_denominator,
                                    "one_to_one_denominator_features": dual["one_to_one_denominator_features"],
                                    "cluster_purity_fraction": dual["cluster_purity_fraction"],
                                    "cluster_purity_score": dual["cluster_purity_score"],
                                    "cluster_purity_pure_clusters": dual["cluster_purity_pure_clusters"],
                                    "cluster_purity_total_clusters": dual["cluster_purity_total_clusters"],
                                    "no_hit_over_unmapped_fraction": dual["no_hit_over_unmapped_fraction"],
                                    "no_hit_n_no_hit": dual["no_hit_n_no_hit"],
                                    "no_hit_n_unmapped": dual["no_hit_n_unmapped"],
                                    "species_coverage": np.nan,
                                    "merge_deficit": dual["merge_deficit"],
                                }
                            )
                        continue

                    # Graph mode (single): each metric can use its own target level by method.
                    split_targets = graph_targets_by_metric.get("split_index", {})
                    merge_targets = graph_targets_by_metric.get("merge_deficit", {})
                    unmapped_targets = graph_targets_by_metric.get("unmapped_feature_fraction", split_targets)
                    one_to_one_targets = graph_targets_by_metric.get("one_to_one_mapping_fraction", split_targets)
                    no_hit_targets = graph_targets_by_metric.get("no_hit_over_unmapped_fraction", unmapped_targets)

                    split_mapped_feats = [f for f in inferred_present if f in split_targets]
                    split_mapped_n = len(split_mapped_feats)
                    merge_mapped_feats = [f for f in inferred_present if f in merge_targets]
                    merge_mapped_n = len(merge_mapped_feats)
                    unmapped_mapped_feats = [f for f in inferred_present if f in unmapped_targets]
                    unmapped_present_feats = [f for f in inferred_present if f not in unmapped_targets]
                    n_mapped = len(unmapped_mapped_feats)
                    n_unmapped = n_inf - n_mapped
                    unmapped_fraction = float(n_unmapped / n_inf)
                    no_hit_unmapped_feats = [f for f in inferred_present if f not in no_hit_targets]
                    no_hit_n_unmapped = int(len(no_hit_unmapped_feats))
                    no_hit_n_no_hit = int(
                        sum(
                            1
                            for f in no_hit_unmapped_feats
                            if not _feature_has_best_hit(
                                method=method,
                                feature_id=str(f),
                                best_effective_pident_by_feature=best_effective_pident_by_feature,
                                best_pident_by_feature=best_pident_by_feature,
                                no_hit_threshold_by_method=no_hit_threshold_by_method,
                            )
                        )
                    )
                    no_hit_over_unmapped_fraction = (
                        float(no_hit_n_no_hit / no_hit_n_unmapped) if no_hit_n_unmapped > 0 else 0.0
                    )
                    _record_unmapped_best_hits(
                        community_id=cid,
                        sample_id=s,
                        sample_key=f"{cid}|{s}",
                        method=method,
                        target_level=graph_target_levels.get("unmapped_feature_fraction", "species"),
                        unmapped_features=unmapped_present_feats,
                        best_pident_by_feature=best_pident_by_feature,
                        best_effective_pident_by_feature=best_effective_pident_by_feature,
                        feature_fraction_by_id=feature_fraction_by_id,
                    )

                    split_target_counts: dict[str, int] = {}
                    split_edge_count = 0
                    for f in split_mapped_feats:
                        targets = split_targets.get(f, set())
                        split_edge_count += len(targets)
                        for tgt in targets:
                            split_target_counts[tgt] = split_target_counts.get(tgt, 0) + 1

                    split_duplicates = sum(max(0, c - 1) for c in split_target_counts.values())
                    split_denominator_features = (
                        len(split_target_counts)
                        if split_index_denominator == "mapped_expected_features"
                        else split_mapped_n
                    )
                    split_index = (
                        float(split_duplicates / split_denominator_features)
                        if split_denominator_features > 0
                        else np.nan
                    )

                    merge_duplicates = 0
                    for f in merge_mapped_feats:
                        deg = len(merge_targets.get(f, set()))
                        merge_duplicates += max(0, deg - 1)
                    if merge_denominator == "expected_present_features":
                        exp_present_merge = set(exp.index[exp[s].astype(float) > 0].astype(str))
                        merge_denominator_features = len(exp_present_merge)
                    else:
                        merge_denominator_features = merge_mapped_n
                    merge_deficit = (
                        float(merge_duplicates / merge_denominator_features)
                        if merge_denominator_features > 0
                        else np.nan
                    )

                    one_to_one_n = 0
                    one_to_one_target_counts: dict[str, int] = {}
                    one_to_one_mapped_feats = [f for f in inferred_present if f in one_to_one_targets]
                    one_to_one_mapped_n = len(one_to_one_mapped_feats)
                    for f in one_to_one_mapped_feats:
                        for tgt in one_to_one_targets.get(f, set()):
                            one_to_one_target_counts[tgt] = one_to_one_target_counts.get(tgt, 0) + 1
                    for f in one_to_one_mapped_feats:
                        targets = one_to_one_targets.get(f, set())
                        if len(targets) == 1:
                            tgt = next(iter(targets))
                            if one_to_one_target_counts.get(tgt, 0) == 1:
                                one_to_one_n += 1
                    if one_to_one_denominator == "expected_present_features":
                        if graph_target_levels.get("one_to_one_mapping_fraction", "species") == "tip":
                            if (not exp_tip.empty) and (s in exp_tip.columns):
                                one_to_one_denominator_features = int((exp_tip[s].astype(float) > 0).sum())
                            else:
                                one_to_one_denominator_features = 0
                        else:
                            exp_present = set(exp.index[exp[s].astype(float) > 0].astype(str))
                            one_to_one_denominator_features = len(exp_present)
                    elif one_to_one_denominator == "inferred_present_features":
                        one_to_one_denominator_features = n_inf
                    else:
                        one_to_one_denominator_features = one_to_one_mapped_n
                    one_to_one_mapping_fraction = (
                        float(one_to_one_n / one_to_one_denominator_features)
                        if one_to_one_denominator_features > 0
                        else np.nan
                    )

                    cluster_purity_fraction = np.nan
                    cluster_purity_score = np.nan
                    cluster_purity_pure_clusters = 0.0
                    cluster_purity_total_clusters = 0.0
                    if (
                        method in PURITY_METHODS
                        and (not zotu_std.empty)
                        and (s in zotu_std.columns)
                        and (not zotu_hits_all.empty)
                    ):
                        cluster_to_members = tic_cluster_members if method == "TIC" else tac_cluster_members
                        if cluster_to_members:
                            purity_level = graph_target_levels.get("cluster_purity_fraction", "species")
                            purity_min_ident = float(graph_min_ident_by_metric.get("cluster_purity_fraction", np.nan))
                            cache_key = (str(purity_level), float(purity_min_ident))
                            member_targets = cluster_member_targets_cache.get(cache_key)
                            if member_targets is None:
                                if np.isfinite(purity_min_ident):
                                    member_targets = _build_member_target_map(
                                        hits_all=zotu_hits_all,
                                        min_ident=purity_min_ident,
                                        min_qcov=fig09_min_qcov,
                                        target_level=purity_level,
                                        tip_to_species=tip_to_species,
                                    )
                                else:
                                    member_targets = {}
                                cluster_member_targets_cache[cache_key] = member_targets
                            zotu_present = set(zotu_std.index[zotu_std[s].astype(float) > 0].astype(str).tolist())
                            purity_stats = _compute_cluster_purity_summary(
                                present_clusters=set(str(f) for f in inferred_present),
                                cluster_to_members=cluster_to_members,
                                present_members=zotu_present,
                                member_targets=member_targets,
                            )
                            cluster_purity_fraction = float(purity_stats["cluster_purity_fraction"])
                            cluster_purity_score = float(purity_stats["cluster_purity_score"])
                            cluster_purity_pure_clusters = float(purity_stats["cluster_purity_pure_clusters"])
                            cluster_purity_total_clusters = float(purity_stats["cluster_purity_total_clusters"])

                    graph_edge_count = float(split_edge_count)
                    graph_merge_duplicates = float(merge_duplicates)
                    unique_graph_levels = set(graph_target_levels.values())
                    graph_target_level_plot = (
                        next(iter(unique_graph_levels)) if len(unique_graph_levels) == 1 else "mixed"
                    )
                rows.append(
                    {
                        "community_id": cid,
                        "sample_id": s,
                        "sample_key": f"{cid}|{s}",
                        "method": method,
                        "inferred_present_features": n_inf,
                        "expected_present_features": expected_present_features,
                        "mapped_present_features": n_mapped,
                        "unmapped_present_features": n_unmapped,
                        "split_mapped_present_features": split_mapped_n,
                        "merge_mapped_present_features": (
                            merge_mapped_n if split_index_mode == "graph" else split_mapped_n
                        ),
                        "merge_denominator": (merge_denominator if split_index_mode == "graph" else np.nan),
                        "merge_denominator_features": (
                            merge_denominator_features if split_index_mode == "graph" else np.nan
                        ),
                        "split_index_mode": split_index_mode,
                        "split_index_threshold_target_level": (
                            split_threshold_target_level
                            if split_index_mode == "threshold_hits"
                            else (
                                graph_target_levels.get("split_index", "species")
                                if split_index_mode == "graph"
                                else "species"
                            )
                        ),
                        "graph_target_level_split_index": (
                            graph_target_levels.get("split_index", np.nan) if split_index_mode == "graph" else np.nan
                        ),
                        "graph_target_level_merge_deficit": (
                            graph_target_levels.get("merge_deficit", np.nan) if split_index_mode == "graph" else np.nan
                        ),
                        "graph_target_level_unmapped_feature_fraction": (
                            graph_target_levels.get("unmapped_feature_fraction", np.nan)
                            if split_index_mode == "graph"
                            else np.nan
                        ),
                        "graph_target_level_one_to_one_mapping_fraction": (
                            graph_target_levels.get("one_to_one_mapping_fraction", np.nan)
                            if split_index_mode == "graph"
                            else np.nan
                        ),
                        "graph_target_level_cluster_purity_fraction": (
                            graph_target_levels.get("cluster_purity_fraction", np.nan)
                            if split_index_mode == "graph"
                            else np.nan
                        ),
                        "graph_target_level_no_hit_over_unmapped_fraction": (
                            graph_target_levels.get("no_hit_over_unmapped_fraction", np.nan)
                            if split_index_mode == "graph"
                            else np.nan
                        ),
                        "graph_target_level_plot": graph_target_level_plot,
                        "split_index_denominator": split_index_denominator,
                        "split_index_denominator_features": split_denominator_features,
                        "graph_edge_threshold_pident": (
                            graph_min_ident_by_metric.get("split_index", graph_min_ident_used)
                            if split_index_mode == "graph"
                            else np.nan
                        ),
                        "graph_edge_threshold_split_index": (
                            graph_min_ident_by_metric.get("split_index", np.nan)
                            if split_index_mode == "graph"
                            else np.nan
                        ),
                        "graph_edge_threshold_merge_deficit": (
                            graph_min_ident_by_metric.get("merge_deficit", np.nan)
                            if split_index_mode == "graph"
                            else np.nan
                        ),
                        "graph_edge_threshold_unmapped_feature_fraction": (
                            graph_min_ident_by_metric.get("unmapped_feature_fraction", np.nan)
                            if split_index_mode == "graph"
                            else np.nan
                        ),
                        "graph_edge_threshold_one_to_one_mapping_fraction": (
                            graph_min_ident_by_metric.get("one_to_one_mapping_fraction", np.nan)
                            if split_index_mode == "graph"
                            else np.nan
                        ),
                        "graph_edge_threshold_cluster_purity_fraction": (
                            graph_min_ident_by_metric.get("cluster_purity_fraction", np.nan)
                            if split_index_mode == "graph"
                            else np.nan
                        ),
                        "graph_edge_threshold_no_hit_over_unmapped_fraction": (
                            graph_min_ident_by_metric.get("no_hit_over_unmapped_fraction", np.nan)
                            if split_index_mode == "graph"
                            else np.nan
                        ),
                        "graph_edge_min_query_coverage": (
                            fig09_min_qcov if split_index_mode in {"threshold_hits", "graph"} else np.nan
                        ),
                        "graph_edge_count": graph_edge_count,
                        "graph_merge_duplicates": graph_merge_duplicates,
                        "split_duplicates": split_duplicates,
                        "split_index": split_index,
                        "unmapped_feature_fraction": unmapped_fraction,
                        "one_to_one_mapping_fraction": one_to_one_mapping_fraction,
                        "one_to_one_denominator": one_to_one_denominator,
                        "one_to_one_denominator_features": one_to_one_denominator_features,
                        "cluster_purity_fraction": (
                            cluster_purity_fraction if split_index_mode == "graph" else np.nan
                        ),
                        "cluster_purity_score": (
                            cluster_purity_score if split_index_mode == "graph" else np.nan
                        ),
                        "cluster_purity_pure_clusters": (
                            cluster_purity_pure_clusters if split_index_mode == "graph" else np.nan
                        ),
                        "cluster_purity_total_clusters": (
                            cluster_purity_total_clusters if split_index_mode == "graph" else np.nan
                        ),
                        "no_hit_over_unmapped_fraction": (
                            no_hit_over_unmapped_fraction if split_index_mode == "graph" else np.nan
                        ),
                        "no_hit_n_no_hit": (no_hit_n_no_hit if split_index_mode == "graph" else np.nan),
                        "no_hit_n_unmapped": (no_hit_n_unmapped if split_index_mode == "graph" else np.nan),
                        "species_coverage": species_coverage,
                        "merge_deficit": merge_deficit,
                    }
                )
    if not rows:
        log.info("Figure 09: no rows to plot.")
        return 0

    df = pd.DataFrame(rows)
    out_dir = figures_root(outdir) / "fig09_split_merge_diagnostics"
    out_dir.mkdir(parents=True, exist_ok=True)
    schematic_path = _write_mapping_diagnostics_schematic(out_dir, cfg, log)
    if schematic_path is not None:
        log.info(f"Wrote Figure 09 schematic to {schematic_path}")
    transition_heatmap_df = pd.DataFrame(
        transition_heatmap_rows,
        columns=[
            "community_id",
            "sample_id",
            "sample_key",
            "transition_method",
            "zotu_id",
            "cluster_id",
            "zotu_class",
            "final_fate",
            "best_hit_effective_identity",
            "best_hit_target_id",
            "best_hit_species",
            "zotu_feature_value_tss1000",
            "near_exact_threshold",
            "exact_identity_threshold",
            "classification_source",
        ],
    )
    transition_heatmap_summary = _build_fig09_transition_flow_summary(transition_heatmap_rows)
    transition_heatmap_df.to_csv(out_dir / "fig09_transition_heatmap_rows.tsv", sep="\t", index=False)
    transition_heatmap_summary.to_csv(out_dir / "fig09_transition_heatmap_summary.tsv", sep="\t", index=False)
    if bool(getattr(cfg.figures, "fig09_homogeneity_enabled", True)):
        homogeneity_cluster_df, homogeneity_summary_df, homogeneity_method_summary_df = (
            _build_fig09_homogeneity_tables(transition_heatmap_df)
        )
        homogeneity_cluster_df.to_csv(
            out_dir / "fig09_cluster_match_status_homogeneity_clusters.tsv",
            sep="\t",
            index=False,
        )
        homogeneity_summary_df.to_csv(
            out_dir / "fig09_cluster_match_status_homogeneity_summary.tsv",
            sep="\t",
            index=False,
        )
        homogeneity_method_summary_df.to_csv(
            out_dir / "fig09_cluster_match_status_homogeneity_method_summary.tsv",
            sep="\t",
            index=False,
        )
        _plot_fig09_homogeneity_companion(
            homogeneity_summary_df,
            homogeneity_method_summary_df,
            out_dir / "fig09_cluster_match_status_homogeneity.svg",
            cfg,
        )
    unmapped_best_df = pd.DataFrame(
        unmapped_best_hit_rows,
        columns=[
            "community_id",
            "sample_id",
            "sample_key",
            "method",
            "graph_target_level_plot",
            "unmapped_feature",
            "best_hit_effective_pident",
            "best_hit_pident",
            "best_hit_status",
            "inferred_feature_fraction",
        ],
    )
    if not unmapped_best_df.empty:
        eff = pd.to_numeric(unmapped_best_df["best_hit_effective_pident"], errors="coerce")
        pid = pd.to_numeric(unmapped_best_df["best_hit_pident"], errors="coerce")
        has_best = eff.notna() | pid.notna()
        unmapped_best_df["best_hit_status"] = np.where(has_best, "hit_found", "no_hit_detected")
    unmapped_best_df.to_csv(
        out_dir / "unmapped_feature_best_hit_identity.tsv",
        sep="\t",
        index=False,
    )
    if not unmapped_best_df.empty:
        nohit_summary = (
            unmapped_best_df.assign(
                has_best_hit=lambda x: (pd.to_numeric(x["best_hit_effective_pident"], errors="coerce").notna())
                | (pd.to_numeric(x["best_hit_pident"], errors="coerce").notna())
            )
            .groupby(["graph_target_level_plot", "method"], as_index=False)
            .agg(
                n_unmapped_features=("unmapped_feature", "count"),
                n_with_best_hit=("has_best_hit", "sum"),
            )
        )
        nohit_summary["n_no_hit_detected"] = (
            nohit_summary["n_unmapped_features"] - nohit_summary["n_with_best_hit"]
        )
        with np.errstate(divide="ignore", invalid="ignore"):
            nohit_summary["no_hit_fraction"] = (
                nohit_summary["n_no_hit_detected"] / nohit_summary["n_unmapped_features"]
            )
        nohit_summary.to_csv(
            out_dir / "unmapped_feature_best_hit_nohit_summary.tsv",
            sep="\t",
            index=False,
        )

    # Ensure Figure 09 has a no-hit/unmapped metric for all modes by deriving
    # it from per-feature unmapped best-hit diagnostics.
    if not unmapped_best_df.empty:
        dnh = unmapped_best_df.copy()
        dnh["method"] = dnh["method"].astype(str)
        dnh["_lvl_key"] = dnh["graph_target_level_plot"].astype(str).str.strip().str.lower()
        eff = pd.to_numeric(dnh["best_hit_effective_pident"], errors="coerce")
        pid = pd.to_numeric(dnh["best_hit_pident"], errors="coerce")
        score = eff.where(eff.notna(), pid)
        has_best = score.notna()
        if no_hit_threshold_by_method:
            thr = dnh["method"].map(no_hit_threshold_by_method)
            use_thr = thr.notna()
            pass_thr = score >= pd.to_numeric(thr, errors="coerce")
            has_best = has_best & (~use_thr | pass_thr.fillna(False))
        dnh["has_best_hit"] = has_best.astype(float)
        dnh["inferred_feature_fraction"] = pd.to_numeric(
            dnh.get("inferred_feature_fraction", np.nan), errors="coerce"
        ).fillna(0.0).clip(lower=0.0)
        dnh["weighted_unmapped"] = dnh["inferred_feature_fraction"]
        dnh["weighted_with_best_hit"] = dnh["inferred_feature_fraction"] * dnh["has_best_hit"]

        agg_nohit = (
            dnh.groupby(["community_id", "sample_id", "sample_key", "method", "_lvl_key"], as_index=False)
            .agg(
                no_hit_n_unmapped=("unmapped_feature", "count"),
                n_with_best_hit=("has_best_hit", "sum"),
                no_hit_w_unmapped=("weighted_unmapped", "sum"),
                no_hit_w_with_best_hit=("weighted_with_best_hit", "sum"),
            )
            .copy()
        )
        agg_nohit["no_hit_n_no_hit"] = agg_nohit["no_hit_n_unmapped"] - agg_nohit["n_with_best_hit"]
        agg_nohit["no_hit_w_no_hit"] = agg_nohit["no_hit_w_unmapped"] - agg_nohit["no_hit_w_with_best_hit"]
        with np.errstate(divide="ignore", invalid="ignore"):
            agg_nohit["no_hit_over_unmapped_fraction"] = np.where(
                agg_nohit["no_hit_n_unmapped"] > 0,
                agg_nohit["no_hit_n_no_hit"] / agg_nohit["no_hit_n_unmapped"],
                0.0,
            )
            agg_nohit["no_hit_over_unmapped_fraction_weighted"] = np.where(
                agg_nohit["no_hit_w_unmapped"] > 0,
                agg_nohit["no_hit_w_no_hit"] / agg_nohit["no_hit_w_unmapped"],
                0.0,
            )

        dfx = df.copy()
        dfx["_lvl_key"] = dfx.get("graph_target_level_plot", "species")
        dfx["_lvl_key"] = dfx["_lvl_key"].astype(str).str.strip().str.lower()
        dfx["_lvl_key"] = dfx["_lvl_key"].replace({"nan": "species", "none": "species", "": "species"})

        merged = dfx.merge(
            agg_nohit[
                [
                    "community_id",
                    "sample_id",
                    "sample_key",
                    "method",
                    "_lvl_key",
                    "no_hit_n_unmapped",
                    "no_hit_n_no_hit",
                    "no_hit_over_unmapped_fraction",
                    "no_hit_w_unmapped",
                    "no_hit_w_no_hit",
                    "no_hit_over_unmapped_fraction_weighted",
                ]
            ],
            on=["community_id", "sample_id", "sample_key", "method", "_lvl_key"],
            how="left",
            suffixes=("", "_agg"),
        )
        for col in (
            "no_hit_n_unmapped",
            "no_hit_n_no_hit",
            "no_hit_over_unmapped_fraction",
            "no_hit_w_unmapped",
            "no_hit_w_no_hit",
            "no_hit_over_unmapped_fraction_weighted",
        ):
            agg_col = f"{col}_agg"
            if agg_col not in merged.columns:
                continue
            if col in merged.columns:
                cur = pd.to_numeric(merged[col], errors="coerce")
                add = pd.to_numeric(merged[agg_col], errors="coerce")
                merged[col] = cur.fillna(add)
            else:
                merged[col] = pd.to_numeric(merged[agg_col], errors="coerce")
        merged = merged.drop(columns=[c for c in merged.columns if c.endswith("_agg")], errors="ignore")
        merged["no_hit_n_unmapped"] = pd.to_numeric(merged.get("no_hit_n_unmapped", np.nan), errors="coerce").fillna(0.0)
        merged["no_hit_n_no_hit"] = pd.to_numeric(merged.get("no_hit_n_no_hit", np.nan), errors="coerce").fillna(0.0)
        merged["no_hit_w_unmapped"] = pd.to_numeric(merged.get("no_hit_w_unmapped", np.nan), errors="coerce").fillna(0.0)
        merged["no_hit_w_no_hit"] = pd.to_numeric(merged.get("no_hit_w_no_hit", np.nan), errors="coerce").fillna(0.0)
        frac = pd.to_numeric(merged.get("no_hit_over_unmapped_fraction", np.nan), errors="coerce")
        merged["no_hit_over_unmapped_fraction"] = np.where(
            merged["no_hit_n_unmapped"] > 0,
            frac,
            0.0,
        )
        frac_w = pd.to_numeric(merged.get("no_hit_over_unmapped_fraction_weighted", np.nan), errors="coerce")
        merged["no_hit_over_unmapped_fraction_weighted"] = np.where(
            merged["no_hit_w_unmapped"] > 0,
            frac_w,
            0.0,
        )
        df = merged.drop(columns=["_lvl_key"], errors="ignore")
    else:
        if "no_hit_w_unmapped" not in df.columns:
            df["no_hit_w_unmapped"] = 0.0
        if "no_hit_w_no_hit" not in df.columns:
            df["no_hit_w_no_hit"] = 0.0
        if "no_hit_over_unmapped_fraction_weighted" not in df.columns:
            df["no_hit_over_unmapped_fraction_weighted"] = 0.0

    # Persist the finalized diagnostics table after derived Figure 09 metrics.
    df.to_csv(out_dir / "split_merge_diagnostics.tsv", sep="\t", index=False)
    richness_sample_df, richness_method_summary_df, richness_effect_df, richness_effect_summary_df = (
        _build_fig09_relative_richness_tables(df)
    )
    richness_sample_df.to_csv(
        out_dir / "fig09_panel_d_richness_error_sample_summary.tsv",
        sep="\t",
        index=False,
    )
    richness_method_summary_df.to_csv(
        out_dir / "fig09_panel_d_richness_error_method_summary.tsv",
        sep="\t",
        index=False,
    )
    richness_effect_df.to_csv(
        out_dir / "fig09_panel_d_richness_error_excess_removed.tsv",
        sep="\t",
        index=False,
    )
    richness_effect_summary_df.to_csv(
        out_dir / "fig09_panel_d_richness_error_excess_removed_summary.tsv",
        sep="\t",
        index=False,
    )

    # Extra Figure 09 companion QC: denoising best-hit effective-identity histogram.
    # This captures the distribution of best-hit effective identities for ASV/zOTU
    # features and writes the requested 97-100% percent-scaled panel.
    if denoising_best_hit_rows:
        denoise_qc_dir = qc_dir(outdir) / "denoising_best_hit_effective_identity_histogram"
        denoise_qc_dir.mkdir(parents=True, exist_ok=True)
        hist_font_size = max(6.0, float(getattr(cfg.figures, "fig09_font_size", 10.0)))
        hist_tick_size = max(6.0, hist_font_size - 2.0)

        denoise_df = pd.DataFrame(
            denoising_best_hit_rows,
            columns=[
                "community_id",
                "method",
                "method_label",
                "inferred_feature",
                "best_target",
                "best_pident",
                "best_qcov_exact_percent",
                "best_effective_identity",
                "accepted_in_strict_mapping",
            ],
        )
        denoise_df["best_pident"] = pd.to_numeric(denoise_df["best_pident"], errors="coerce")
        denoise_df["best_qcov_exact_percent"] = pd.to_numeric(
            denoise_df["best_qcov_exact_percent"], errors="coerce"
        )
        denoise_df["best_effective_identity"] = pd.to_numeric(
            denoise_df["best_effective_identity"], errors="coerce"
        )
        denoise_df.to_csv(
            denoise_qc_dir / "denoising_best_hit_effective_identity_features.tsv",
            sep="\t",
            index=False,
        )

        pooled_summary = (
            denoise_df.groupby(["method", "method_label"], as_index=False)
            .agg(
                total_inferred_features=("inferred_feature", "count"),
                features_with_recorded_hit=("best_effective_identity", lambda x: int(pd.to_numeric(x, errors="coerce").notna().sum())),
                features_without_recorded_hit=("best_effective_identity", lambda x: int(pd.to_numeric(x, errors="coerce").isna().sum())),
                features_best_eff_identity_lt_90=("best_effective_identity", lambda x: int((pd.to_numeric(x, errors="coerce") < 90.0).sum())),
                features_best_eff_identity_90_to_100=("best_effective_identity", lambda x: int(((pd.to_numeric(x, errors="coerce") >= 90.0) & (pd.to_numeric(x, errors="coerce") < 100.0)).sum())),
                features_best_eff_identity_eq_100=("best_effective_identity", lambda x: int((pd.to_numeric(x, errors="coerce") >= 100.0).sum())),
                features_accepted_in_mapping=("accepted_in_strict_mapping", lambda x: int(pd.Series(x).fillna(False).astype(bool).sum())),
            )
            .sort_values(["method"], key=lambda s: s.map({m: i for i, m in enumerate(METHOD_ORDER)}).fillna(999))
        )
        pooled_summary.to_csv(
            denoise_qc_dir / "denoising_best_hit_effective_identity_summary_pooled.tsv",
            sep="\t",
            index=False,
        )

        by_community_summary = (
            denoise_df.groupby(["community_id", "method", "method_label"], as_index=False)
            .agg(
                total_inferred_features=("inferred_feature", "count"),
                features_with_recorded_hit=("best_effective_identity", lambda x: int(pd.to_numeric(x, errors="coerce").notna().sum())),
                features_without_recorded_hit=("best_effective_identity", lambda x: int(pd.to_numeric(x, errors="coerce").isna().sum())),
                features_best_eff_identity_lt_90=("best_effective_identity", lambda x: int((pd.to_numeric(x, errors="coerce") < 90.0).sum())),
                features_best_eff_identity_90_to_100=("best_effective_identity", lambda x: int(((pd.to_numeric(x, errors="coerce") >= 90.0) & (pd.to_numeric(x, errors="coerce") < 100.0)).sum())),
                features_best_eff_identity_eq_100=("best_effective_identity", lambda x: int((pd.to_numeric(x, errors="coerce") >= 100.0).sum())),
                features_accepted_in_mapping=("accepted_in_strict_mapping", lambda x: int(pd.Series(x).fillna(False).astype(bool).sum())),
            )
            .sort_values(["community_id", "method"], key=lambda s: s.map({m: i for i, m in enumerate(METHOD_ORDER)}).fillna(999))
        )
        by_community_summary.to_csv(
            denoise_qc_dir / "denoising_best_hit_effective_identity_summary_by_community.tsv",
            sep="\t",
            index=False,
        )

        # Requested panel: x=[97,100], bin=0.1, y=percent with 5-unit ticks.
        bins = np.round(np.arange(97.0, 100.0 + 0.100001, 0.1), 6)
        bins_rows: list[dict] = []
        panel_methods = [m for m in ("ASV", "zOTU") if m in denoise_df["method"].astype(str).unique().tolist()]
        for method in panel_methods:
            sub = denoise_df[denoise_df["method"] == method].copy()
            label = display_methods([method])[0]
            total_n = int(len(sub))
            vals = pd.to_numeric(sub["best_effective_identity"], errors="coerce").to_numpy(dtype=float)
            vals = vals[np.isfinite(vals)]
            vals = vals[(vals >= 97.0) & (vals <= 100.0)]
            counts, edges = np.histogram(vals, bins=bins)
            for i in range(len(counts)):
                c = int(counts[i])
                pct = (100.0 * c / total_n) if total_n > 0 else 0.0
                bins_rows.append(
                    {
                        "method": method,
                        "method_label": label,
                        "bin_left": float(edges[i]),
                        "bin_right": float(edges[i + 1]),
                        "count": c,
                        "total_inferred_features": float(total_n),
                        "percent_of_total_inferred_features": float(pct),
                    }
                )
        bins_df = pd.DataFrame(
            bins_rows,
            columns=[
                "method",
                "method_label",
                "bin_left",
                "bin_right",
                "count",
                "total_inferred_features",
                "percent_of_total_inferred_features",
            ],
        )
        bins_df.to_csv(
            denoise_qc_dir / "denoising_best_hit_effective_identity_bins_0p1_x97_100_percent.tsv",
            sep="\t",
            index=False,
        )

        if not bins_df.empty:
            fig_h = plt.figure(figsize=(10.5, 5.8))
            ax_h = fig_h.add_subplot(1, 1, 1)

            for method in panel_methods:
                sub = bins_df[bins_df["method"] == method].sort_values("bin_left")
                if sub.empty:
                    continue
                x_left = sub["bin_left"].to_numpy(dtype=float)
                x_right = sub["bin_right"].to_numpy(dtype=float)
                y = sub["percent_of_total_inferred_features"].to_numpy(dtype=float)
                if x_left.size <= 0:
                    continue
                x_step = np.append(x_left, x_right[-1])
                y_step = np.append(y, y[-1] if y.size > 0 else 0.0)
                color = METHOD_COLORS.get(method, "#666666")
                label = display_methods([method])[0]
                ax_h.fill_between(
                    x_step,
                    y_step,
                    step="post",
                    alpha=0.18,
                    color=color,
                    edgecolor=color,
                    linewidth=0.8,
                )
                ax_h.step(
                    x_step,
                    y_step,
                    where="post",
                    color=color,
                    linewidth=2.0,
                    label=label,
                )

            ax_h.set_xlim(97.0, 100.0)
            xt = np.arange(97.0, 100.0001, 0.5)
            ax_h.set_xticks(xt)
            ax_h.set_xticklabels([f"{float(t):g}" for t in xt], fontsize=hist_tick_size)

            ymax = float(pd.to_numeric(bins_df["percent_of_total_inferred_features"], errors="coerce").max())
            if not np.isfinite(ymax):
                ymax = 0.0
            ytop = max(5.0, 5.0 * math.ceil(max(0.0, ymax) / 5.0))
            ax_h.set_ylim(0.0, ytop)
            yt = np.arange(0.0, ytop + 0.001, 5.0)
            ax_h.set_yticks(yt)
            ax_h.set_yticklabels([f"{int(v)}" for v in yt], fontsize=hist_tick_size)

            ax_h.grid(True, axis="y", alpha=0.8, linewidth=0.7, color="#dddddd")
            ax_h.grid(True, axis="x", alpha=0.8, linewidth=0.6, color="#eeeeee")
            ax_h.set_xlabel("Best-hit effective identity to expected sequence (%)", fontsize=hist_font_size)
            ax_h.set_ylabel("Inferred features per 0.1% bin (%)", fontsize=hist_font_size)
            ax_h.set_title(
                "Denoising Features Span a Range of Best-Hit Effective Identities",
                fontsize=hist_font_size + 1.0,
            )
            ax_h.legend(frameon=False, fontsize=hist_font_size, loc="upper right")
            fig_h.tight_layout()

            fig_h.savefig(
                denoise_qc_dir / "denoising_best_hit_effective_identity_histogram_x97_100_percent_ticks0p5_y5.png",
                format="png",
                dpi=300,
            )
            fig_h.savefig(
                denoise_qc_dir / "denoising_best_hit_effective_identity_histogram_x97_100_percent_ticks0p5_y5.svg",
                format="svg",
            )
            # Also store this companion Figure 09 panel inside the main fig09 output dir
            # so users can find all Figure 09 graphics together.
            fig_h.savefig(
                out_dir / "denoising_best_hit_effective_identity_histogram_x97_100_percent_ticks0p5_y5.png",
                format="png",
                dpi=300,
            )
            fig_h.savefig(
                out_dir / "denoising_best_hit_effective_identity_histogram_x97_100_percent_ticks0p5_y5.svg",
                format="svg",
            )
            plt.close(fig_h)

    panels = [
        ("split_index", "Split Index"),
        ("merge_deficit", "Merge Index" if split_index_mode == "graph" else "Merge Deficit"),
        ("one_to_one_mapping_fraction", "1-to-1 Mapping Fraction"),
        ("relative_richness_error", "Relative Richness Error"),
    ]
    rng = np.random.default_rng(1)
    title = getattr(cfg.figures, "fig09_title", None) or "Split/Merge Diagnostics"
    plot_type = _normalize_plot_type(getattr(cfg.figures, "fig09_plot_type", "boxplot"), default="boxplot")
    dual_level_output_mode = False
    point_size = max(1.0, float(getattr(cfg.figures, "box_violin_dot_size", 22.0)))
    fig09_font_size = max(6.0, float(getattr(cfg.figures, "fig09_font_size", 10.0)))
    fig09_tick_size = max(6.0, fig09_font_size - 2.0)
    fig09_suptitle_size = fig09_font_size + 2.0
    fig09_stats_text_size = max(6.0, fig09_font_size - 3.0)
    fig09_layout_w_pad = max(0.0, float(getattr(cfg.figures, "fig09_layout_w_pad", 3.8)))
    fig09_layout_rect = [
        float(np.clip(float(getattr(cfg.figures, "fig09_layout_left", 0.02)), 0.0, 1.0)),
        float(np.clip(float(getattr(cfg.figures, "fig09_layout_bottom", 0.02)), 0.0, 1.0)),
        float(np.clip(float(getattr(cfg.figures, "fig09_layout_right", 0.96)), 0.0, 1.0)),
        float(np.clip(float(getattr(cfg.figures, "fig09_layout_top", 0.95)), 0.0, 1.0)),
    ]
    if fig09_layout_rect[2] <= fig09_layout_rect[0]:
        fig09_layout_rect[0], fig09_layout_rect[2] = 0.02, 0.96
    if fig09_layout_rect[3] <= fig09_layout_rect[1]:
        fig09_layout_rect[1], fig09_layout_rect[3] = 0.02, 0.95
    fig09_save_bbox_tight = bool(getattr(cfg.figures, "fig09_save_bbox_tight", False))
    fig09_save_pad_inches = max(0.0, float(getattr(cfg.figures, "fig09_save_pad_inches", 0.02)))
    fig09_write_individual_panels = bool(getattr(cfg.figures, "fig09_write_individual_panels", True))

    def _save_fig09_main(fig, path: Path) -> None:
        save_kwargs: dict[str, object] = {"format": "svg"}
        if fig09_save_bbox_tight:
            save_kwargs["bbox_inches"] = "tight"
            save_kwargs["pad_inches"] = fig09_save_pad_inches
        fig.savefig(path, **save_kwargs)
    stats_enabled = bool(getattr(cfg.figures, "fig09_stats_enabled", True))
    stats_permutations = max(100, int(getattr(cfg.figures, "fig09_stats_permutations", 4000)))
    stats_seed = int(getattr(cfg.figures, "fig09_stats_seed", 1))
    stats_alpha = float(getattr(cfg.figures, "fig09_stats_alpha", 0.05))
    stats_show_pairwise = bool(getattr(cfg.figures, "fig09_stats_show_pairwise", True))
    stats_model = _normalize_stats_model(getattr(cfg.figures, "fig09_stats_model", "mixed_effect_glmmtmb"))
    stats_pairwise_test = _normalize_pairwise_test(
        getattr(cfg.figures, "fig09_stats_pairwise_test", "signed_rank_permutation")
    )
    stats_multiple_test_correction = bool(
        getattr(cfg.figures, "fig09_stats_multiple_test_correction", True)
    )
    glmm_bin_formula = str(
        getattr(
            cfg.figures,
            "fig09_stats_glmm_binomial_formula",
            "cbind(success, failure) ~ method + (1|community_id) + (1|sample_key)",
        )
    )
    glmm_bin_null_formula = str(
        getattr(
            cfg.figures,
            "fig09_stats_glmm_binomial_null_formula",
            "cbind(success, failure) ~ 1 + (1|community_id) + (1|sample_key)",
        )
    )
    glmm_count_formula = str(
        getattr(
            cfg.figures,
            "fig09_stats_glmm_count_formula",
            "count ~ method + offset(log_offset) + (1|community_id) + (1|sample_key)",
        )
    )
    glmm_count_null_formula = str(
        getattr(
            cfg.figures,
            "fig09_stats_glmm_count_null_formula",
            "count ~ 1 + offset(log_offset) + (1|community_id) + (1|sample_key)",
        )
    )
    glmm_count_family = _normalize_glmm_count_family(getattr(cfg.figures, "fig09_stats_glmm_count_family", "nbinom2"))
    glmm_emmeans_scale = _normalize_emmeans_scale(getattr(cfg.figures, "fig09_stats_glmm_emmeans_scale", "link"))
    glmm_spec = {
        "binomial_formula": glmm_bin_formula,
        "binomial_null_formula": glmm_bin_null_formula,
        "count_formula": glmm_count_formula,
        "count_null_formula": glmm_count_null_formula,
        "count_family": glmm_count_family,
        "emmeans_scale": glmm_emmeans_scale,
    }

    def _panel_d_methods() -> list[str]:
        return _normalize_method_order(
            getattr(cfg.figures, "fig09_panel_d_methods", FIG09_RICHNESS_METHOD_ORDER),
            default=FIG09_RICHNESS_METHOD_ORDER,
        )

    def _panel_d_stats_df() -> pd.DataFrame:
        d = richness_sample_df.copy()
        if d.empty:
            return d
        d["relative_richness_error_for_stats"] = pd.to_numeric(
            d.get("relative_richness_error_percent", np.nan),
            errors="coerce",
        )
        return d.dropna(subset=["relative_richness_error_for_stats"]).copy()

    def _fmt_name_num(raw: float) -> str:
        v = float(raw)
        if abs(v - round(v)) < 1e-9:
            return str(int(round(v)))
        return f"{v:.3f}".rstrip("0").rstrip(".").replace(".", "p")

    def _draw_unmapped_best_hit_panel(
        data_df: pd.DataFrame,
        out_name: str,
        fig_title: str,
    ) -> None:
        fig = plt.figure(figsize=(8.2, 4.8))
        ax = fig.add_subplot(1, 1, 1)

        method_data: list[tuple[int, str, np.ndarray]] = []
        all_vals: list[np.ndarray] = []
        for j, method in enumerate(METHOD_ORDER, start=1):
            sub = data_df[data_df["method"] == method]
            vals = sub["best_hit_effective_pident"].to_numpy(dtype=float)
            vals = vals[np.isfinite(vals)]
            if vals.size > 0 and float(np.nanmax(vals)) <= 1.0:
                vals = vals * 100.0
            if vals.size > 0:
                method_data.append((j, method, vals))
                all_vals.append(vals)

        if not method_data:
            ax.text(
                0.5,
                0.5,
                "No unmapped features with best-hit identity available.",
                ha="center",
                va="center",
                fontsize=fig09_font_size,
            )
        else:
            violin_vals: list[np.ndarray] = []
            violin_pos: list[float] = []
            violin_methods: list[str] = []
            for j, method, vals in method_data:
                if vals.size >= 2:
                    violin_vals.append(vals)
                    violin_pos.append(float(j))
                    violin_methods.append(method)
            if violin_vals:
                vp = ax.violinplot(
                    violin_vals,
                    positions=violin_pos,
                    widths=0.64,
                    showmeans=False,
                    showmedians=True,
                    showextrema=False,
                )
                for body, method in zip(vp["bodies"], violin_methods):
                    body.set_facecolor(METHOD_COLORS.get(method, "#999999"))
                    body.set_edgecolor("#333333")
                    body.set_alpha(0.33)
                if "cmedians" in vp:
                    vp["cmedians"].set_color("#111111")
                    vp["cmedians"].set_linewidth(1.2)

            for j, method, vals in method_data:
                xj = j + rng.normal(0, 0.04, size=vals.size)
                ax.scatter(
                    xj,
                    vals,
                    s=point_size,
                    alpha=0.45,
                    color=METHOD_COLORS.get(method, "#666666"),
                    edgecolors="none",
                )

        ax.set_title("Best-Hit Effective Identity of Unmapped Inferred Features", fontsize=fig09_font_size)
        ax.set_ylabel("Best hit effective identity (%)", fontsize=fig09_font_size)
        ax.set_xticks(range(1, len(METHOD_ORDER) + 1))
        ax.set_xticklabels(
            display_methods(METHOD_ORDER),
            rotation=45,
            ha="right",
            rotation_mode="anchor",
            fontsize=fig09_tick_size,
        )
        ax.tick_params(axis="y", labelsize=fig09_tick_size)
        ax.grid(True, axis="y", alpha=0.25, linewidth=0.6)
        if all_vals:
            yvals = np.concatenate(all_vals)
            y_min = float(np.nanmin(yvals))
            y_max = float(np.nanmax(yvals))
            if np.isfinite(y_min) and np.isfinite(y_max):
                if y_max <= y_min:
                    pad = 0.5
                else:
                    pad = max(0.25, 0.08 * (y_max - y_min))
                low = max(0.0, y_min - pad)
                high = min(100.0, y_max + pad)
                if high <= low:
                    high = min(100.0, low + 1.0)
                ax.set_ylim(low, high)
            else:
                ax.set_ylim(0.0, 100.0)
        else:
            ax.set_ylim(0.0, 100.0)
        fig.suptitle(fig_title, fontsize=fig09_suptitle_size)
        fig.tight_layout(rect=[0, 0, 0.88, 0.93])
        fig.savefig(out_dir / out_name, format="svg")
        plt.close(fig)

    def _draw_unmapped_no_hit_fraction_panel(
        data_df: pd.DataFrame,
        out_name: str,
        fig_title: str,
    ) -> None:
        fig = plt.figure(figsize=(8.2, 4.8))
        ax = fig.add_subplot(1, 1, 1)

        bars_x: list[int] = []
        bars_y: list[float] = []
        bar_colors: list[str] = []
        labels: list[str] = []
        annots: list[str] = []
        for j, method in enumerate(METHOD_ORDER, start=1):
            sub = data_df[data_df["method"] == method].copy()
            n_total = int(len(sub))
            if n_total <= 0:
                bars_x.append(j)
                bars_y.append(0.0)
                bar_colors.append(METHOD_COLORS.get(method, "#999999"))
                labels.append(display_methods([method])[0])
                annots.append("0/0")
                continue
            eff = pd.to_numeric(sub["best_hit_effective_pident"], errors="coerce")
            pid = pd.to_numeric(sub["best_hit_pident"], errors="coerce")
            n_with = int((eff.notna() | pid.notna()).sum())
            n_no = int(max(0, n_total - n_with))
            frac = float(n_no / n_total)
            bars_x.append(j)
            bars_y.append(frac)
            bar_colors.append(METHOD_COLORS.get(method, "#999999"))
            labels.append(display_methods([method])[0])
            annots.append(f"{n_no}/{n_total}")

        ax.bar(
            bars_x,
            bars_y,
            color=bar_colors,
            alpha=0.75,
            edgecolor="#333333",
            linewidth=0.8,
            width=0.65,
        )
        for x, y, txt in zip(bars_x, bars_y, annots):
            ax.text(
                x,
                min(1.0, y + 0.02),
                txt,
                ha="center",
                va="bottom",
                fontsize=fig09_tick_size,
                color="#222222",
            )

        ax.set_title("Unmapped Features Without Any Detectable Best Match", fontsize=fig09_font_size)
        ax.set_ylabel("No-match fraction among unmapped features", fontsize=fig09_font_size)
        ax.set_ylim(0.0, 1.0)
        ax.set_xticks(range(1, len(METHOD_ORDER) + 1))
        ax.set_xticklabels(
            labels,
            rotation=45,
            ha="right",
            rotation_mode="anchor",
            fontsize=fig09_tick_size,
        )
        ax.tick_params(axis="y", labelsize=fig09_tick_size)
        ax.grid(True, axis="y", alpha=0.25, linewidth=0.6)

        fig.suptitle(fig_title, fontsize=fig09_suptitle_size)
        fig.tight_layout(rect=[0, 0, 0.88, 0.93])
        fig.savefig(out_dir / out_name, format="svg")
        plt.close(fig)

    def _draw_unmapped_weighted_burden_heatmap(
        data_df: pd.DataFrame,
        ref_df: pd.DataFrame,
        out_name: str,
        fig_title: str,
    ) -> None:
        sample_order = list(dict.fromkeys(ref_df["sample_key"].astype(str).tolist()))
        if not sample_order:
            sample_order = list(dict.fromkeys(data_df["sample_key"].astype(str).tolist()))
        sample_set = set(sample_order)
        pair_set = set(
            (str(sk), str(m))
            for sk, m in zip(ref_df["sample_key"].astype(str), ref_df["method"].astype(str))
        )

        if data_df.empty:
            heat_vals = np.full((len(sample_order), len(METHOD_ORDER)), np.nan, dtype=float)
            for i, sk in enumerate(sample_order):
                for j, m in enumerate(METHOD_ORDER):
                    if (str(sk), str(m)) in pair_set:
                        heat_vals[i, j] = 0.0
            summary_out = pd.DataFrame(
                [
                    {
                        "sample_key": sk,
                        "method": m,
                        "n_unmapped_features": 0 if (str(sk), str(m)) in pair_set else np.nan,
                        "sum_unmapped_feature_fraction": 0.0 if (str(sk), str(m)) in pair_set else np.nan,
                        "weighted_unmapped_burden": 0.0 if (str(sk), str(m)) in pair_set else np.nan,
                    }
                    for sk in sample_order
                    for m in METHOD_ORDER
                ]
            )
        else:
            dd = data_df.copy()
            dd["sample_key"] = dd["sample_key"].astype(str)
            dd["method"] = dd["method"].astype(str)
            dd["inferred_feature_fraction"] = pd.to_numeric(
                dd.get("inferred_feature_fraction", np.nan), errors="coerce"
            ).fillna(0.0)
            agg = (
                dd.groupby(["sample_key", "method"], as_index=False)
                .agg(
                    n_unmapped_features=("unmapped_feature", "count"),
                    sum_unmapped_feature_fraction=("inferred_feature_fraction", "sum"),
                )
            )
            # Requested behavior: do not multiply by n_unmapped features.
            # Heatmap value is the summed unmapped fraction only.
            agg["weighted_unmapped_burden"] = agg["sum_unmapped_feature_fraction"]
            agg_map = {
                (str(r["sample_key"]), str(r["method"])): (
                    float(r["n_unmapped_features"]),
                    float(r["sum_unmapped_feature_fraction"]),
                    float(r["weighted_unmapped_burden"]),
                )
                for _, r in agg.iterrows()
            }

            rows_out: list[dict] = []
            heat_vals = np.full((len(sample_order), len(METHOD_ORDER)), np.nan, dtype=float)
            for i, sk in enumerate(sample_order):
                for j, m in enumerate(METHOD_ORDER):
                    key = (str(sk), str(m))
                    if key in agg_map:
                        n_unm, sum_frac, w_burden = agg_map[key]
                    elif key in pair_set:
                        n_unm, sum_frac, w_burden = 0.0, 0.0, 0.0
                    else:
                        n_unm, sum_frac, w_burden = np.nan, np.nan, np.nan
                    heat_vals[i, j] = w_burden
                    rows_out.append(
                        {
                            "sample_key": sk,
                            "method": m,
                            "n_unmapped_features": n_unm,
                            "sum_unmapped_feature_fraction": sum_frac,
                            "weighted_unmapped_burden": w_burden,
                        }
                    )
            summary_out = pd.DataFrame(rows_out)

        summary_out = summary_out.copy()
        summary_out["sum_unmapped_feature_fraction_percent"] = (
            pd.to_numeric(summary_out["sum_unmapped_feature_fraction"], errors="coerce") * 100.0
        )
        summary_out["weighted_unmapped_burden_percent"] = (
            pd.to_numeric(summary_out["weighted_unmapped_burden"], errors="coerce") * 100.0
        )
        summary_out["sample_label"] = summary_out["sample_key"].astype(str).map(display_sample_key)
        summary_out["method_display"] = summary_out["method"].astype(str).map(lambda x: display_methods([x])[0])
        summary_out.to_csv(
            out_dir / out_name.replace(".svg", ".tsv"),
            sep="\t",
            index=False,
        )

        h = max(4.8, 0.28 * max(1, len(sample_order)) + 2.2)
        fig = plt.figure(figsize=(8.8, h))
        ax = fig.add_subplot(1, 1, 1)

        if heat_vals.size == 0:
            ax.text(
                0.5,
                0.5,
                "No sample-method data available.",
                ha="center",
                va="center",
                fontsize=fig09_font_size,
            )
            fig.suptitle(fig_title, fontsize=fig09_suptitle_size)
            fig.tight_layout(rect=[0, 0, 1, 0.94])
            fig.savefig(out_dir / out_name, format="svg")
            plt.close(fig)
            return

        # Color-blind friendly, light-to-dark for 0->100%
        cmap = plt.cm.cividis_r.copy()
        cmap.set_bad(color="#eeeeee")
        heat_vals_pct = heat_vals * 100.0
        im = ax.imshow(
            heat_vals_pct,
            aspect="auto",
            interpolation="nearest",
            cmap=cmap,
            vmin=0.0,
            vmax=100.0,
        )
        cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
        cbar.set_label("Summed unmapped feature fraction (%)", fontsize=fig09_font_size)
        cbar.ax.tick_params(labelsize=fig09_tick_size)
        cbar.set_ticks([0, 20, 40, 60, 80, 100])

        ax.set_title(
            "Unmapped Fraction Burden Heatmap\n"
            "(sum_unmapped_feature_fraction)",
            fontsize=fig09_font_size,
        )
        ax.set_xticks(np.arange(len(METHOD_ORDER)))
        ax.set_xticklabels(
            display_methods(METHOD_ORDER),
            rotation=45,
            ha="right",
            rotation_mode="anchor",
            fontsize=fig09_tick_size,
        )
        ax.set_yticks(np.arange(len(sample_order)))
        ax.set_yticklabels([display_sample_key(s) for s in sample_order], fontsize=max(6.0, fig09_tick_size - 1.0))
        ax.set_xlabel("Method", fontsize=fig09_font_size)
        ax.set_ylabel("Sample", fontsize=fig09_font_size)

        fig.suptitle(fig_title, fontsize=fig09_suptitle_size)
        fig.tight_layout(rect=[0, 0, 1, 0.94])
        fig.savefig(out_dir / out_name, format="svg")
        plt.close(fig)

    def _graph_target_signature_single() -> str:
        metric_order = [
            "split_index",
            "merge_deficit",
            "one_to_one_mapping_fraction",
            "relative_richness_error",
        ]
        metric_short = {
            "split_index": "si",
            "merge_deficit": "mg",
            "one_to_one_mapping_fraction": "o2",
            "relative_richness_error": "re",
        }
        method_short = {"ASV": "asv", "zOTU": "zo", "OTU": "ot", "TIC": "ti", "TAC": "ta"}
        segs: list[str] = []
        for metric in metric_order:
            if metric == "relative_richness_error":
                segs.append("re-richness")
                continue
            bits = []
            for method in METHOD_ORDER:
                lvl = _graph_target_level_for_method(cfg, method, metric)
                bits.append(f"{method_short[method]}{str(lvl)[0].lower()}")
            segs.append(f"{metric_short[metric]}-{'-'.join(bits)}")
        return "__".join(segs)

    def _figure_stem(level_override: str | None = None) -> str:
        parts = [f"split_merge_diagnostics__mode-{split_index_mode}"]
        if split_index_mode == "graph":
            if (
                level_override in {"tip", "species", "species_and_tip", "panel_abcd"}
                or str(level_override or "").startswith("panel_")
            ):
                parts.append(f"target-{level_override}")
            elif dual_level_output_mode:
                parts.append("target-tip_species")
            else:
                parts.append(f"target-{_graph_target_signature_single()}")
            id_d = _as_percent_identity(
                getattr(
                    cfg.figures,
                    "fig09_graph_min_effective_pident_denoising",
                    getattr(cfg.figures, "fig09_graph_min_identity_denoising", 100.0),
                )
            )
            id_tt = _as_percent_identity(
                getattr(
                    cfg.figures,
                    "fig09_graph_min_effective_pident_tic_tac",
                    getattr(cfg.figures, "fig09_graph_min_identity_tic_tac", 98.7),
                )
            )
            id_o = _as_percent_identity(
                getattr(
                    cfg.figures,
                    "fig09_graph_min_effective_pident_otu",
                    getattr(cfg.figures, "fig09_graph_min_identity_otu", 97.0),
                )
            )
            parts.append(f"idD{_fmt_name_num(id_d)}-idTT{_fmt_name_num(id_tt)}-idO{_fmt_name_num(id_o)}")
            custom_thr = getattr(cfg.figures, "fig09_graph_edge_threshold_by_metric", {}) or {}
            if isinstance(custom_thr, dict) and custom_thr:
                parts.append("idByMetric-custom")
            parts.append(f"qcov{_fmt_name_num(fig09_min_qcov)}")
        elif split_index_mode == "threshold_hits":
            parts.append(f"target-{split_threshold_target_level}")
            parts.append(f"qcov{_fmt_name_num(fig09_min_qcov)}")
        return "__".join(parts)

    def _safe_output_stem(stem: str, max_len: int = 170) -> str:
        """
        Keep output filenames under common FS limits (255 bytes per filename element)
        after appending suffixes like "_stats_pairwise.tsv".
        """
        s = str(stem or "")
        if len(s) <= max_len:
            return s
        digest = hashlib.sha1(s.encode("utf-8")).hexdigest()[:12]
        keep = max(16, max_len - len("__h") - len(digest))
        short = s[:keep].rstrip("_-")
        return f"{short}__h{digest}"

    def _compute_panel_stats_bundle(data_df: pd.DataFrame, stem: str, seed_bias: int = 0) -> dict[str, dict]:
        panel_stats_local: dict[str, dict] = {}
        if not stats_enabled:
            return panel_stats_local
        omnibus_rows: list[dict] = []
        pairwise_rows: list[dict] = []
        for idx, (col, label) in enumerate(panels, start=1):
            if col == "no_hit_over_unmapped_fraction":
                continue
            stats_df = data_df
            stats_metric_col = col
            stats_method_order = METHOD_ORDER
            if col == "relative_richness_error":
                stats_df = _panel_d_stats_df()
                stats_metric_col = "relative_richness_error_for_stats"
                stats_method_order = _panel_d_methods()
                if stats_df.empty:
                    continue
            st = _compute_panel_stats(
                df=stats_df,
                metric_col=stats_metric_col,
                method_order=stats_method_order,
                n_perm=stats_permutations,
                seed=stats_seed + seed_bias + (idx * 10000),
                alpha=stats_alpha,
                pairwise_test=stats_pairwise_test,
                multiple_test_correction=stats_multiple_test_correction,
                stats_model=stats_model,
                glmm_spec=glmm_spec,
                log=log,
            )
            panel_stats_local[col] = st
            omnibus_rows.append(
                {
                    "metric": col,
                    "label": label,
                    "stats_model": st.get("stats_model", stats_model),
                    "pairwise_test": stats_pairwise_test,
                    "multiple_test_correction": bool(st.get("multiple_test_correction", stats_multiple_test_correction)),
                    "n_blocks": st["n_blocks"],
                    "friedman_q": st["friedman_q"],
                    "friedman_p_perm": st["friedman_p"],
                    "friedman_significant_alpha": bool(np.isfinite(st["friedman_p"]) and st["friedman_p"] < stats_alpha),
                    "alpha": stats_alpha,
                    "n_permutations": stats_permutations,
                }
            )
            for row in st["pairwise"]:
                pairwise_rows.append(
                    {
                        "metric": col,
                        "label": label,
                        "method_a": row["method_a"],
                        "method_b": row["method_b"],
                        "stats_model": st.get("stats_model", stats_model),
                        "pairwise_test": row.get("pairwise_test", stats_pairwise_test),
                        "statistic": row.get("statistic", row.get("w_plus", float("nan"))),
                        "w_plus": row["w_plus"],
                        "n_nonzero_pairs": row["n_nonzero_pairs"],
                        "n_pos_pairs": row.get("n_pos_pairs", float("nan")),
                        "n_neg_pairs": row.get("n_neg_pairs", float("nan")),
                        "n_zero_pairs": row.get("n_zero_pairs", float("nan")),
                        "median_diff_a_minus_b": row["median_diff_a_minus_b"],
                        "estimate": row.get("estimate", row["median_diff_a_minus_b"]),
                        "std_error": row.get("std_error", float("nan")),
                        "model_family": row.get("model_family", ""),
                        "p_raw": row["p_raw"],
                        "p_adj_bh": row["p_adj_bh"],
                        "p_value_used": row.get("p_value_used", row["p_adj_bh"]),
                        "p_adjustment": row.get("p_adjustment", "bh"),
                        "multiple_test_correction": bool(
                            row.get("multiple_test_correction", stats_multiple_test_correction)
                        ),
                        "significant_alpha": row.get("significant_alpha", row["significant_bh"]),
                        "significant_bh_alpha": row["significant_bh"],
                        "alpha": stats_alpha,
                        "n_permutations": stats_permutations,
                    }
                )
        pd.DataFrame(omnibus_rows).to_csv(out_dir / f"{stem}_stats_omnibus.tsv", sep="\t", index=False)
        pd.DataFrame(pairwise_rows).to_csv(out_dir / f"{stem}_stats_pairwise.tsv", sep="\t", index=False)
        return panel_stats_local

    def _combined_shared_ylim(level_dfs: dict[str, pd.DataFrame]) -> dict[str, tuple[float, float]]:
        out: dict[str, tuple[float, float]] = {}
        for metric in ("split_index", "merge_deficit"):
            vals_all: list[np.ndarray] = []
            for d in level_dfs.values():
                if metric not in d.columns:
                    continue
                vv = d[metric].to_numpy(dtype=float)
                vv = vv[np.isfinite(vv)]
                if vv.size > 0:
                    vals_all.append(vv)
            if not vals_all:
                continue
            vals = np.concatenate(vals_all)
            ymin = float(np.nanmin(vals))
            ymax = float(np.nanmax(vals))
            if not np.isfinite(ymin) or not np.isfinite(ymax):
                continue
            if ymax <= ymin:
                span = max(1.0, abs(ymax))
                pad = 0.08 * span
            else:
                pad = 0.08 * (ymax - ymin)
            low = max(0.0, ymin - pad)
            high = ymax + pad
            out[metric] = (low, high)
        return out

    def _metric_level_column(metric_col: str) -> str:
        by_metric = {
            "split_index": "graph_target_level_split_index",
            "merge_deficit": "graph_target_level_merge_deficit",
            "unmapped_feature_fraction": "graph_target_level_unmapped_feature_fraction",
            "one_to_one_mapping_fraction": "graph_target_level_one_to_one_mapping_fraction",
            "cluster_purity_fraction": "graph_target_level_cluster_purity_fraction",
            "no_hit_over_unmapped_fraction": "graph_target_level_no_hit_over_unmapped_fraction",
        }
        return by_metric.get(metric_col, "graph_target_level_plot")

    def _select_metric_level_df(
        data_df: pd.DataFrame,
        metric_col: str,
        forced_level: str | None = None,
    ) -> pd.DataFrame:
        if split_index_mode != "graph":
            return data_df.copy()
        d = data_df.copy()
        level_col = _metric_level_column(metric_col)
        if level_col not in d.columns:
            return d

        levels = d[level_col].astype(str).str.lower()
        if forced_level in {"tip", "species"}:
            out = d[levels == str(forced_level).lower()].copy()
            return out if not out.empty else d

        expected_by_method = {
            m: _graph_target_level_for_method(cfg, m, metric_col).lower() for m in METHOD_ORDER
        }
        expected = d["method"].astype(str).map(lambda m: expected_by_method.get(str(m), "species"))
        out = d[levels == expected.astype(str).str.lower()].copy()
        return out if not out.empty else d

    def _shared_ylim_for_metric_sets(metric_sets: dict[str, list[pd.DataFrame]]) -> dict[str, tuple[float, float]]:
        out: dict[str, tuple[float, float]] = {}
        for metric, dfs_for_metric in metric_sets.items():
            vals_all: list[np.ndarray] = []
            for d in dfs_for_metric:
                if metric not in d.columns:
                    continue
                vv = pd.to_numeric(d[metric], errors="coerce").to_numpy(dtype=float)
                vv = vv[np.isfinite(vv)]
                if vv.size > 0:
                    vals_all.append(vv)
            if not vals_all:
                continue
            vals = np.concatenate(vals_all)
            ymin = float(np.nanmin(vals))
            ymax = float(np.nanmax(vals))
            if not np.isfinite(ymin) or not np.isfinite(ymax):
                continue
            if ymax <= ymin:
                span = max(1.0, abs(ymax))
                pad = 0.08 * span
            else:
                pad = 0.08 * (ymax - ymin)
            low = max(0.0, ymin - pad)
            high = ymax + pad
            out[metric] = (low, high)
        # Keep Split and Merge on the same y-axis limits for direct visual comparison.
        if "split_index" in out and "merge_deficit" in out:
            low = min(float(out["split_index"][0]), float(out["merge_deficit"][0]))
            high = max(float(out["split_index"][1]), float(out["merge_deficit"][1]))
            out["split_index"] = (low, high)
            out["merge_deficit"] = (low, high)
        return out

    def _default_panel_title(panel_key: str) -> str:
        defaults = {
            "a": "Split/Merge Index",
            "b": "1-to-1 Mapping",
            "c": "Unmapped and No-match Diagnostics",
        }
        return defaults.get(panel_key, panel_key)

    def _panel_title(panel_key: str) -> str:
        raw = getattr(cfg.figures, "fig09_panel_titles", {}) or {}
        if isinstance(raw, dict):
            for k, v in raw.items():
                if _normalize_panel_key(k) == panel_key:
                    return str(v)
        return _default_panel_title(panel_key)

    def _default_panel_method_levels(panel_key: str, metric_col: str) -> dict[str, str]:
        if panel_key == "a" and metric_col in {"split_index", "merge_deficit"}:
            return {m: "species" for m in METHOD_ORDER}
        if panel_key == "b" and metric_col in {"one_to_one_mapping_fraction"}:
            return {m: "species" for m in METHOD_ORDER}
        if panel_key == "c" and metric_col in {"unmapped_feature_fraction", "no_hit_over_unmapped_fraction"}:
            return {m: _graph_target_level_for_method(cfg, m, "unmapped_feature_fraction") for m in METHOD_ORDER}
        return {m: _graph_target_level_for_method(cfg, m, metric_col) for m in METHOD_ORDER}

    def _panel_method_levels(panel_key: str, metric_col: str) -> dict[str, str]:
        levels = _default_panel_method_levels(panel_key, metric_col)
        raw = getattr(cfg.figures, "fig09_panel_target_level_by_metric", {}) or {}
        if not isinstance(raw, dict) or (not raw):
            return levels
        panel_node = None
        for pk, pv in raw.items():
            if _normalize_panel_key(pk) == panel_key:
                panel_node = pv
                break
        if panel_node is None:
            return levels
        # Shorthand: panel node as a single level string applies to all metrics/methods.
        if isinstance(panel_node, str):
            lvl = _normalize_split_target_level(panel_node, default="species")
            return {m: lvl for m in METHOD_ORDER}
        if not isinstance(panel_node, dict):
            return levels

        metric_node = None
        metric_key = _normalize_graph_metric_key(metric_col)
        for mk, mv in panel_node.items():
            if _normalize_graph_metric_key(mk) == metric_key:
                metric_node = mv
                break
        # Shorthand: panel node as method-map applies to all metrics in the panel.
        if metric_node is None:
            methodish = {"asv", "zotu", "otu", "tic", "tac", "default"}
            if any(_normalize_method_token(k) in methodish for k in panel_node.keys()):
                metric_node = panel_node
        if metric_node is None:
            return levels
        if isinstance(metric_node, str):
            lvl = _normalize_split_target_level(metric_node, default="species")
            return {m: lvl for m in METHOD_ORDER}
        if not isinstance(metric_node, dict):
            return levels

        if "default" in metric_node:
            lvl = _normalize_split_target_level(metric_node.get("default"), default="species")
            levels = {m: lvl for m in METHOD_ORDER}
        for k, v in metric_node.items():
            mk = _normalize_method_token(k)
            for method in METHOD_ORDER:
                if _normalize_method_token(method) == mk:
                    levels[method] = _normalize_split_target_level(v, default=levels.get(method, "species"))
                    break
        return levels

    def _select_df_for_metric_levels(
        data_df: pd.DataFrame,
        metric_col: str,
        method_levels: dict[str, str],
    ) -> pd.DataFrame:
        if split_index_mode != "graph":
            return data_df.copy()
        if data_df.empty:
            return data_df.copy()
        level_col = _metric_level_column(metric_col)
        if level_col not in data_df.columns:
            return data_df.copy()

        chunks: list[pd.DataFrame] = []
        for method in METHOD_ORDER:
            sub = data_df[data_df["method"] == method].copy()
            if sub.empty:
                continue
            desired = _normalize_split_target_level(method_levels.get(method, "species"), default="species")
            got = sub[level_col].astype(str).str.lower()
            chosen = sub[got == desired].copy()
            if chosen.empty:
                # Fallback to any available level for this method.
                chosen = sub
            chunks.append(chosen)
        if not chunks:
            return data_df.copy()
        out = pd.concat(chunks, ignore_index=True)
        return out

    def _metric_text_overrides(
        metric_col: str,
        default_title: str,
        default_xlabel: str = "",
        default_ylabel: str = "",
    ) -> tuple[str, str, str]:
        title = str(default_title)
        xlabel = str(default_xlabel)
        ylabel = str(default_ylabel)
        raw = getattr(cfg.figures, "fig09_graph_labels_by_metric", {}) or {}
        if not isinstance(raw, dict) or (not raw):
            return title, xlabel, ylabel

        metric_key = _normalize_graph_metric_key(metric_col)
        node = None
        for rk, rv in raw.items():
            if _normalize_graph_metric_key(str(rk)) == metric_key:
                node = rv
                break
        if node is None and "default" in raw:
            node = raw.get("default")
        if node is None:
            return title, xlabel, ylabel

        if isinstance(node, str):
            return str(node), xlabel, ylabel
        if isinstance(node, dict):
            if "title" in node:
                title = str(node.get("title", title))
            if "xlabel" in node:
                xlabel = str(node.get("xlabel", xlabel))
            elif "x" in node:
                xlabel = str(node.get("x", xlabel))
            if "ylabel" in node:
                ylabel = str(node.get("ylabel", ylabel))
            elif "y" in node:
                ylabel = str(node.get("y", ylabel))
        return title, xlabel, ylabel

    def _draw_transition_heatmap_panel(
        ax,
        metric_col: str,
        label: str,
    ) -> None:
        title_txt, xlabel_txt, ylabel_txt = _metric_text_overrides(
            metric_col,
            label,
            default_xlabel="Post-clustering outcome",
            default_ylabel="Original zOTU class",
        )
        title_txt = str(getattr(cfg.figures, "fig09_transition_heatmap_title", None) or title_txt)
        xlabel_txt = str(getattr(cfg.figures, "fig09_transition_heatmap_xlabel", None) or xlabel_txt)
        ylabel_txt = str(getattr(cfg.figures, "fig09_transition_heatmap_ylabel", None) or ylabel_txt)
        no_match_label = str(getattr(cfg.figures, "fig09_transition_heatmap_no_match_label", "No\nmatch"))
        non_bacterial_label = str(
            getattr(cfg.figures, "fig09_transition_heatmap_non_bacterial_label", "Non\nbacterial")
        )
        row_labels_raw = FIG09_TRANSITION_CLASS_ORDER
        col_labels_raw = FIG09_TRANSITION_OUTCOME_ORDER
        row_labels = [
            _fig09_transition_display_label(label, no_match_label, non_bacterial_label)
            for label in row_labels_raw
        ]
        col_labels = [
            _fig09_transition_display_label(label, no_match_label, non_bacterial_label)
            for label in col_labels_raw
        ]
        methods = ["TIC", "TAC"]
        show_counts = bool(getattr(cfg.figures, "fig09_transition_heatmap_show_counts", True))
        show_percentages = bool(getattr(cfg.figures, "fig09_transition_heatmap_show_percentages", True))
        show_colorbar = bool(getattr(cfg.figures, "fig09_transition_heatmap_show_colorbar", True))
        x_tick_rotation = float(getattr(cfg.figures, "fig09_transition_heatmap_x_tick_rotation", 0.0))
        row_label_rotation = float(getattr(cfg.figures, "fig09_transition_heatmap_row_label_rotation", 0.0))
        method_label_x_offset = float(
            getattr(cfg.figures, "fig09_transition_heatmap_method_label_x_offset", 0.05)
        )
        colorbar_pad = float(getattr(cfg.figures, "fig09_transition_heatmap_colorbar_pad", 0.18))
        tic_label = str(getattr(cfg.figures, "fig09_transition_heatmap_tic_label", "TIC"))
        tac_label = str(getattr(cfg.figures, "fig09_transition_heatmap_tac_label", "TAC"))
        method_labels = {"TIC": tic_label, "TAC": tac_label}

        if transition_heatmap_summary.empty:
            ax.text(
                0.5,
                0.5,
                "No TIC/TAC transition data available",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize=fig09_font_size,
            )
            ax.set_title(title_txt, fontsize=fig09_font_size)
            ax.axis("off")
            return

        method_matrices: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        for method in methods:
            sub = transition_heatmap_summary[
                transition_heatmap_summary["transition_method"].astype(str) == method
            ].copy()
            heat = np.zeros((len(row_labels_raw), len(col_labels_raw)), dtype=float)
            counts = np.zeros_like(heat, dtype=int)
            row_totals = np.zeros(len(row_labels_raw), dtype=int)
            for i, origin in enumerate(row_labels_raw):
                origin_sub = sub[sub["zotu_class"].astype(str) == origin]
                row_total = int(pd.to_numeric(origin_sub["count"], errors="coerce").fillna(0).sum())
                row_totals[i] = row_total
                for j, fate in enumerate(col_labels_raw):
                    hit = origin_sub[origin_sub["final_fate"].astype(str) == fate]
                    count = int(pd.to_numeric(hit["count"], errors="coerce").fillna(0).sum()) if not hit.empty else 0
                    counts[i, j] = count
                    heat[i, j] = (float(count) / float(row_total)) if row_total > 0 else 0.0
            method_matrices[method] = (heat, counts, row_totals)

        cmap_name = str(getattr(cfg.figures, "fig09_transition_heatmap_cmap", "YlGnBu"))
        try:
            cmap = plt.get_cmap(cmap_name)
        except ValueError:
            cmap = plt.get_cmap("YlGnBu")
        norm = plt.Normalize(vmin=0.0, vmax=1.0)
        half_specs = {"TIC": (0.0, 0.5), "TAC": (0.5, 0.5)}

        for i, _origin in enumerate(row_labels_raw):
            for j, _fate in enumerate(col_labels_raw):
                for method in methods:
                    heat, counts, _row_totals = method_matrices[method]
                    y_offset, height = half_specs[method]
                    frac = float(heat[i, j])
                    count = int(counts[i, j])
                    ax.add_patch(
                        Rectangle(
                            (j, i + y_offset),
                            1.0,
                            height,
                            facecolor=cmap(norm(frac)),
                            edgecolor="white",
                            linewidth=1.0,
                        )
                    )
                    text_parts: list[str] = []
                    if show_counts:
                        text_parts.append(str(count))
                    if show_percentages:
                        text_parts.append(f"{frac * 100:.1f}%")
                    if text_parts:
                        ax.text(
                            j + 0.5,
                            i + y_offset + (height / 2.0),
                            "\n".join(text_parts),
                            ha="center",
                            va="center",
                            fontsize=max(5.2, fig09_tick_size - 2.8),
                            color="white" if frac >= 0.55 else "#222222",
                            linespacing=0.92,
                        )

        for x in range(len(col_labels) + 1):
            ax.axvline(float(x), color="#333333", linewidth=0.6)
        for y in range(len(row_labels) + 1):
            ax.axhline(float(y), color="#333333", linewidth=0.6)
        for y in range(len(row_labels)):
            ax.axhline(float(y) + 0.5, color="white", linewidth=1.0)

        ax.set_xlim(0.0, float(len(col_labels)))
        ax.set_ylim(float(len(row_labels)), 0.0)
        ax.set_aspect("auto")
        ax.set_title(title_txt, fontsize=fig09_font_size)
        ax.set_xticks(np.arange(len(col_labels)) + 0.5)
        ax.set_xticklabels(
            col_labels,
            rotation=x_tick_rotation,
            ha="center" if abs(x_tick_rotation) < 1e-9 else "right",
            rotation_mode="anchor",
            fontsize=max(5.8, fig09_tick_size - 1.0),
        )
        ax.set_yticks(np.arange(len(row_labels)) + 0.5)
        ax.set_yticklabels(
            row_labels,
            rotation=row_label_rotation,
            ha="right",
            va="center",
            rotation_mode="anchor",
            fontsize=max(5.8, fig09_tick_size - 1.0),
        )
        for i in range(len(row_labels_raw)):
            ax.text(
                float(len(col_labels)) + method_label_x_offset,
                i + 0.25,
                method_labels["TIC"],
                ha="left",
                va="center",
                fontsize=max(5.5, fig09_tick_size - 1.8),
                fontweight="normal",
                color="#333333",
                clip_on=False,
            )
            ax.text(
                float(len(col_labels)) + method_label_x_offset,
                i + 0.75,
                method_labels["TAC"],
                ha="left",
                va="center",
                fontsize=max(5.5, fig09_tick_size - 1.8),
                fontweight="normal",
                color="#333333",
                clip_on=False,
            )
        if xlabel_txt.strip():
            ax.set_xlabel(xlabel_txt, fontsize=fig09_font_size)
        if ylabel_txt.strip():
            ax.set_ylabel(ylabel_txt, fontsize=fig09_font_size)
        ax.tick_params(axis="both", length=0)
        note = str(getattr(cfg.figures, "fig09_transition_heatmap_note", "") or "")
        if note.strip():
            ax.text(
                0.5,
                1.02,
                note,
                transform=ax.transAxes,
                ha="center",
                va="bottom",
                fontsize=max(5.8, fig09_tick_size - 1.5),
                color="#555555",
            )
        if show_colorbar:
            sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
            sm.set_array([])
            cbar = ax.figure.colorbar(sm, ax=ax, fraction=0.044, pad=colorbar_pad)
            cbar.ax.set_ylabel(
                str(getattr(cfg.figures, "fig09_transition_heatmap_colorbar_label", "Row fraction of zOTUs")),
                fontsize=max(5.8, fig09_tick_size - 1.3),
            )
            cbar.ax.tick_params(labelsize=max(5.8, fig09_tick_size - 1.5))

    def _draw_relative_richness_error_panel(
        ax,
        label: str,
        panel_stats_local: dict[str, dict],
    ) -> None:
        title_txt, xlabel_txt, ylabel_txt = _metric_text_overrides(
            "relative_richness_error",
            label,
            default_xlabel="Method",
            default_ylabel="Relative richness error (%)",
        )
        panel_title = getattr(cfg.figures, "fig09_panel_d_title", None)
        if panel_title is not None:
            title_txt = str(panel_title)
        panel_xlabel = getattr(cfg.figures, "fig09_panel_d_xlabel", None)
        if panel_xlabel is not None:
            xlabel_txt = str(panel_xlabel)
        panel_ylabel = getattr(cfg.figures, "fig09_panel_d_ylabel", None)
        if panel_ylabel is not None:
            ylabel_txt = str(panel_ylabel)

        panel_font_raw = getattr(cfg.figures, "fig09_panel_d_font_size", None)
        panel_font_size = fig09_font_size if panel_font_raw is None else max(6.0, float(panel_font_raw))
        panel_tick_raw = getattr(cfg.figures, "fig09_panel_d_tick_size", None)
        panel_tick_size = max(5.0, panel_font_size - 2.0) if panel_tick_raw is None else max(5.0, float(panel_tick_raw))
        panel_plot_type = _normalize_plot_type(
            getattr(cfg.figures, "fig09_panel_d_plot_type", plot_type),
            default=plot_type,
        )
        methods = _panel_d_methods()
        xtick_rotation = float(getattr(cfg.figures, "fig09_panel_d_xtick_rotation", 45.0))
        ytick_rotation = float(getattr(cfg.figures, "fig09_panel_d_ytick_rotation", 0.0))
        dot_size_raw = getattr(cfg.figures, "fig09_panel_d_point_size", None)
        dot_size = max(1.0, float(point_size if dot_size_raw is None else dot_size_raw))
        dot_opacity = float(np.clip(float(getattr(cfg.figures, "fig09_panel_d_point_opacity", 0.72)), 0.0, 1.0))
        dot_jitter = max(0.0, float(getattr(cfg.figures, "fig09_panel_d_point_jitter", 0.055)))
        show_zero_line = bool(getattr(cfg.figures, "fig09_panel_d_show_zero_line", True))
        show_reference_note = bool(getattr(cfg.figures, "fig09_panel_d_show_reference_note", True))
        reference_note = str(
            getattr(cfg.figures, "fig09_panel_d_reference_note", "0% = expected species richness")
        )
        show_excess_annotation = bool(
            getattr(cfg.figures, "fig09_panel_d_show_excess_removed_annotation", True)
        )

        plot_df = richness_sample_df[richness_sample_df["method"].astype(str).isin(methods)].copy()
        positions = np.arange(1, len(methods) + 1, dtype=float)
        data = [
            pd.to_numeric(
                plot_df.loc[plot_df["method"].astype(str) == method, "relative_richness_error_percent"],
                errors="coerce",
            )
            .dropna()
            .to_numpy(dtype=float)
            for method in methods
        ]

        if plot_df.empty or not any(vals.size > 0 for vals in data):
            ax.text(
                0.5,
                0.5,
                "No richness-error data available",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize=panel_font_size,
            )
        elif panel_plot_type == "violin":
            violin_vals: list[np.ndarray] = []
            violin_pos: list[float] = []
            violin_methods: list[str] = []
            for x, method, vals in zip(positions, methods, data):
                if vals.size >= 2:
                    violin_vals.append(vals)
                    violin_pos.append(float(x))
                    violin_methods.append(method)
            if violin_vals:
                vp = ax.violinplot(
                    violin_vals,
                    positions=violin_pos,
                    widths=0.64,
                    showmeans=False,
                    showmedians=True,
                    showextrema=False,
                )
                for body, method in zip(vp["bodies"], violin_methods):
                    body.set_facecolor(METHOD_COLORS.get(method, "#999999"))
                    body.set_edgecolor("#333333")
                    body.set_alpha(0.33)
                if "cmedians" in vp:
                    vp["cmedians"].set_color("#111111")
                    vp["cmedians"].set_linewidth(1.2)
        else:
            for x, method, vals in zip(positions, methods, data):
                if vals.size == 0:
                    continue
                ax.boxplot(
                    [vals],
                    positions=[float(x)],
                    widths=0.52,
                    patch_artist=True,
                    showfliers=False,
                    boxprops=dict(
                        facecolor=METHOD_COLORS.get(method, "#999999"),
                        alpha=0.38,
                        edgecolor="#333333",
                    ),
                    medianprops=dict(color="#111111", linewidth=1.2),
                    whiskerprops=dict(color="#333333", linewidth=0.9),
                    capprops=dict(color="#333333", linewidth=0.9),
                )

        for x, method, vals in zip(positions, methods, data):
            if vals.size == 0:
                continue
            xj = float(x) + rng.normal(0.0, dot_jitter, size=vals.size)
            ax.scatter(
                xj,
                vals,
                s=dot_size,
                alpha=dot_opacity,
                color=METHOD_COLORS.get(method, "#666666"),
                edgecolors="#333333",
                linewidths=0.35,
                zorder=3,
            )

        if show_zero_line:
            ax.axhline(0.0, color="#333333", linestyle=(0, (4, 3)), linewidth=1.0)
        if show_reference_note and reference_note.strip():
            ax.text(
                0.99,
                0.02,
                reference_note,
                transform=ax.transAxes,
                ha="right",
                va="bottom",
                fontsize=max(5.5, panel_tick_size - 1.0),
                color="#333333",
            )
        if show_excess_annotation and not richness_effect_summary_df.empty:
            effect_lines: list[str] = []
            for method in ("TAC", "TIC"):
                if method not in methods:
                    continue
                hit = richness_effect_summary_df[
                    richness_effect_summary_df["clustering_method"].astype(str) == method
                ]
                if hit.empty:
                    continue
                val = float(hit["median_percent_excess_removed"].iloc[0])
                n = int(hit["n_samples_with_unoise3_excess"].iloc[0])
                if np.isfinite(val):
                    effect_lines.append(f"{display_methods([method])[0]} median excess removed: {val:.1f}% (n={n})")
            if effect_lines:
                ax.text(
                    0.02,
                    0.98,
                    "\n".join(effect_lines),
                    transform=ax.transAxes,
                    ha="left",
                    va="top",
                    fontsize=max(5.5, panel_tick_size - 1.0),
                    bbox=dict(boxstyle="round,pad=0.28", facecolor="white", edgecolor="#BBBBBB", alpha=0.88),
                    zorder=6,
                )

        ax.set_title(str(title_txt), fontsize=panel_font_size)
        ax.set_xticks(positions)
        ax.set_xticklabels(
            display_methods(methods),
            rotation=xtick_rotation,
            ha="center" if abs(xtick_rotation) < 1e-9 else "right",
            rotation_mode="anchor",
            fontsize=panel_tick_size,
        )
        if str(xlabel_txt).strip():
            ax.set_xlabel(str(xlabel_txt), fontsize=panel_font_size)
        if str(ylabel_txt).strip():
            ax.set_ylabel(str(ylabel_txt), fontsize=panel_font_size)
        ax.tick_params(axis="y", labelsize=panel_tick_size)
        for tick in ax.get_yticklabels():
            tick.set_rotation(ytick_rotation)
            tick.set_va("center")
            tick.set_ha("right")
        ax.grid(axis="y", color="#E6E6E6", linewidth=0.8)
        ax.set_axisbelow(True)

        finite_vals = (
            np.concatenate([vals[np.isfinite(vals)] for vals in data if vals.size > 0])
            if any(vals.size > 0 for vals in data)
            else np.array([], dtype=float)
        )
        y_limits_raw = getattr(cfg.figures, "fig09_panel_d_y_limits", None)
        if isinstance(y_limits_raw, (list, tuple)) and len(y_limits_raw) >= 2:
            ymin = float(y_limits_raw[0])
            ymax = float(y_limits_raw[1])
            if np.isfinite(ymin) and np.isfinite(ymax) and ymax > ymin:
                ax.set_ylim(ymin, ymax)
        elif finite_vals.size:
            ymin = float(np.nanmin(finite_vals))
            ymax = float(np.nanmax(finite_vals))
            span = max(1.0, ymax - ymin)
            ax.set_ylim(ymin - 0.12 * span, ymax + 0.18 * span)

        if stats_enabled and ("relative_richness_error" in panel_stats_local):
            st = panel_stats_local["relative_richness_error"]
            if stats_show_pairwise:
                _draw_pairwise_brackets(
                    ax=ax,
                    pairwise_rows=list(st.get("pairwise_sig", [])),
                    method_order=methods,
                    alpha=stats_alpha,
                    constrain_to_axis=False,
                    placement="top",
                    spacing_scale=1.45,
                )

    def _draw_one_panel(
        ax,
        col: str,
        label: str,
        data_df: pd.DataFrame,
        panel_stats_local: dict[str, dict],
        shared_ylim: dict[str, tuple[float, float]] | None = None,
    ) -> None:
        if col == "no_hit_over_unmapped_fraction":
            _draw_transition_heatmap_panel(ax=ax, metric_col=col, label=label)
            return
        if col == "relative_richness_error":
            _draw_relative_richness_error_panel(
                ax=ax,
                label=label,
                panel_stats_local=panel_stats_local,
            )
            return

        is_fraction_panel = col in {
            "unmapped_feature_fraction",
            "one_to_one_mapping_fraction",
            "cluster_purity_fraction",
        }
        method_data: list[tuple[int, str, np.ndarray]] = []
        for j, method in enumerate(METHOD_ORDER, start=1):
            sub = data_df[data_df["method"] == method]
            vals = sub[col].to_numpy(dtype=float)
            vals = vals[np.isfinite(vals)]
            if vals.size == 0:
                continue
            method_data.append((j, method, vals))

        if plot_type == "violin":
            violin_vals: list[np.ndarray] = []
            violin_pos: list[float] = []
            violin_methods: list[str] = []
            for j, method, vals in method_data:
                if vals.size >= 2:
                    violin_vals.append(vals)
                    violin_pos.append(float(j))
                    violin_methods.append(method)
            if violin_vals:
                vp = ax.violinplot(
                    violin_vals,
                    positions=violin_pos,
                    widths=0.64,
                    showmeans=False,
                    showmedians=True,
                    showextrema=False,
                )
                for body, method in zip(vp["bodies"], violin_methods):
                    body.set_facecolor(METHOD_COLORS.get(method, "#999999"))
                    body.set_edgecolor("#333333")
                    body.set_alpha(0.33)
                if "cmedians" in vp:
                    vp["cmedians"].set_color("#111111")
                    vp["cmedians"].set_linewidth(1.2)
        else:
            for j, method, vals in method_data:
                ax.boxplot(
                    [vals],
                    positions=[j],
                    widths=0.52,
                    patch_artist=True,
                    boxprops=dict(
                        facecolor=METHOD_COLORS.get(method, "#999999"),
                        alpha=0.33,
                        edgecolor="#333333",
                    ),
                    medianprops=dict(color="#111111", linewidth=1.2),
                    whiskerprops=dict(color="#333333", linewidth=0.9),
                    capprops=dict(color="#333333", linewidth=0.9),
                )

        for j, method, vals in method_data:
            xj = j + rng.normal(0, 0.04, size=vals.size)
            ax.scatter(
                xj,
                vals,
                s=point_size,
                alpha=0.45,
                color=METHOD_COLORS.get(method, "#666666"),
                edgecolors="none",
            )

        title_txt, xlabel_txt, ylabel_txt = _metric_text_overrides(
            col,
            label,
            default_xlabel="",
            default_ylabel="",
        )
        ax.set_title(title_txt, fontsize=fig09_font_size)
        ax.set_xticks(range(1, len(METHOD_ORDER) + 1))
        ax.set_xticklabels(
            display_methods(METHOD_ORDER),
            rotation=45,
            ha="right",
            rotation_mode="anchor",
            fontsize=fig09_tick_size,
        )
        if str(xlabel_txt).strip():
            ax.set_xlabel(str(xlabel_txt), fontsize=fig09_font_size)
        if str(ylabel_txt).strip():
            ax.set_ylabel(str(ylabel_txt), fontsize=fig09_font_size)
        ax.tick_params(axis="y", labelsize=fig09_tick_size)
        ax.grid(True, axis="y", alpha=0.25, linewidth=0.6)
        if is_fraction_panel:
            ax.set_ylim(0.0, 1.0)
        if shared_ylim and (col in shared_ylim):
            ax.set_ylim(*shared_ylim[col])

        if stats_enabled and (col in panel_stats_local):
            st = panel_stats_local[col]
            if stats_show_pairwise:
                bracket_placement = "bottom" if col == "one_to_one_mapping_fraction" else "top"
                constrain = bool(is_fraction_panel or (shared_ylim and (col in shared_ylim)))
                _draw_pairwise_brackets(
                    ax=ax,
                    pairwise_rows=list(st.get("pairwise_sig", [])),
                    method_order=METHOD_ORDER,
                    alpha=stats_alpha,
                    constrain_to_axis=constrain,
                    placement=bracket_placement,
                    spacing_scale=1.45,
                )
                if col == "unmapped_feature_fraction":
                    pairwise_label = (
                        f"Pairwise BH threshold (alpha={stats_alpha:.3g})\n*** q<0.001, ** q<0.01, * q<{stats_alpha:.3g}"
                        if stats_multiple_test_correction
                        else f"Pairwise threshold (alpha={stats_alpha:.3g})\n*** p<0.001, ** p<0.01, * p<{stats_alpha:.3g}"
                    )
                    ax.text(
                        0.985,
                        0.985,
                        pairwise_label,
                        transform=ax.transAxes,
                        ha="right",
                        va="top",
                        fontsize=max(fig09_stats_text_size, 7.6),
                        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="#bbbbbb", alpha=0.85),
                        zorder=7,
                    )

    def _draw_panel_set(
        panel_set: list[tuple[str, str]],
        out_name: str,
        fig_title: str,
        data_df: pd.DataFrame,
        panel_stats_local: dict[str, dict],
        shared_ylim: dict[str, tuple[float, float]] | None = None,
    ) -> None:
        if len(panel_set) == 4:
            n_rows = 2
            n_cols = 2
        else:
            n_rows = 1
            n_cols = len(panel_set)

        fig = plt.figure(figsize=(4.0 * n_cols, 4.2 * n_rows))
        for i, (col, label) in enumerate(panel_set, start=1):
            ax = fig.add_subplot(n_rows, n_cols, i)
            _draw_one_panel(
                ax=ax,
                col=col,
                label=label,
                data_df=data_df,
                panel_stats_local=panel_stats_local,
                shared_ylim=shared_ylim,
            )

        fig.suptitle(fig_title, fontsize=fig09_suptitle_size)
        fig.tight_layout(rect=fig09_layout_rect, w_pad=fig09_layout_w_pad)
        _save_fig09_main(fig, out_dir / out_name)
        plt.close(fig)
        if fig09_write_individual_panels:
            base_stem = Path(out_name).with_suffix("").name
            panel_letters = ["a", "b", "c", "d"]
            for panel_letter, (col, label) in zip(panel_letters, panel_set):
                fig_single = plt.figure(figsize=(5.2, 4.6))
                ax_single = fig_single.add_subplot(1, 1, 1)
                _draw_one_panel(
                    ax=ax_single,
                    col=col,
                    label=label,
                    data_df=data_df,
                    panel_stats_local={col: panel_stats_local.get(col, {})} if col in panel_stats_local else {},
                    shared_ylim=shared_ylim,
                )
                fig_single.tight_layout()
                _save_fig09_main(fig_single, out_dir / f"{base_stem}_panel_{panel_letter}.svg")
                plt.close(fig_single)

    if split_index_mode == "graph":
        metric_label = {
            "split_index": "Split Index",
            "merge_deficit": "Merge Index",
            "one_to_one_mapping_fraction": "1-to-1 Mapping Fraction",
            "relative_richness_error": "Relative Richness Error",
        }
        metric_short = {
            "split_index": "si",
            "merge_deficit": "mg",
            "one_to_one_mapping_fraction": "o2",
            "relative_richness_error": "re",
        }
        method_short = {"ASV": "asv", "zOTU": "zo", "OTU": "ot", "TIC": "ti", "TAC": "ta"}

        metric_order = [
            "split_index",
            "merge_deficit",
            "one_to_one_mapping_fraction",
            "relative_richness_error",
        ]
        # Keep panel-specific target-level knobs compatible with old config:
        # a: split/merge, b: one-to-one, d: richness error.
        panel_key_by_metric = {
            "split_index": "a",
            "merge_deficit": "a",
            "one_to_one_mapping_fraction": "b",
            "relative_richness_error": "d",
        }

        def _metric_target_signature(metric_cols_local: list[str]) -> str:
            segs: list[str] = []
            for metric_col_local in metric_cols_local:
                if metric_col_local == "relative_richness_error":
                    segs.append("re-richness")
                    continue
                panel_key_local = panel_key_by_metric.get(metric_col_local, "a")
                levels_local = _panel_method_levels(panel_key_local, metric_col_local)
                bits = []
                for method_local in METHOD_ORDER:
                    lvl = _normalize_split_target_level(
                        levels_local.get(method_local, "species"),
                        default="species",
                    )
                    bits.append(f"{method_short[method_local]}{lvl[0]}")
                segs.append(f"{metric_short[metric_col_local]}-{'-'.join(bits)}")
            return "__".join(segs)

        axis_specs: list[tuple[str, str, str, pd.DataFrame]] = []  # (panel_key_legacy, metric, label, df)
        split_axis_dfs: list[pd.DataFrame] = []
        merge_axis_dfs: list[pd.DataFrame] = []
        for col in metric_order:
            panel_key = panel_key_by_metric.get(col, "a")
            mlevels = _panel_method_levels(panel_key, col)
            axis_df = _select_df_for_metric_levels(df, col, mlevels)
            axis_specs.append((panel_key, col, metric_label[col], axis_df))
            if col == "split_index":
                split_axis_dfs.append(axis_df)
            elif col == "merge_deficit":
                merge_axis_dfs.append(axis_df)

        shared_ylim = _shared_ylim_for_metric_sets(
            {
                "split_index": split_axis_dfs,
                "merge_deficit": merge_axis_dfs,
            }
        )

        stem = _safe_output_stem(
            _figure_stem(
            level_override=f"layout-2x2__{_metric_target_signature(metric_order)}"
            )
        )
        axis_stats: dict[tuple[str, str], dict] = {}
        omnibus_rows: list[dict] = []
        pairwise_rows: list[dict] = []
        if stats_enabled:
            for idx, (panel_key, col, label, axis_df) in enumerate(axis_specs, start=1):
                if col == "no_hit_over_unmapped_fraction":
                    continue
                stats_df = axis_df
                stats_metric_col = col
                stats_method_order = METHOD_ORDER
                if col == "relative_richness_error":
                    stats_df = _panel_d_stats_df()
                    stats_metric_col = "relative_richness_error_for_stats"
                    stats_method_order = _panel_d_methods()
                    if stats_df.empty:
                        continue
                st = _compute_panel_stats(
                    df=stats_df,
                    metric_col=stats_metric_col,
                    method_order=stats_method_order,
                    n_perm=stats_permutations,
                    seed=stats_seed + 2000000 + (idx * 10000),
                    alpha=stats_alpha,
                    pairwise_test=stats_pairwise_test,
                    multiple_test_correction=stats_multiple_test_correction,
                    stats_model=stats_model,
                    glmm_spec=glmm_spec,
                    log=log,
                )
                axis_stats[(panel_key, col)] = st
                omnibus_rows.append(
                    {
                        "panel_group": panel_key,
                        "panel_title": _panel_title(panel_key),
                        "metric": col,
                        "label": label,
                        "stats_model": st.get("stats_model", stats_model),
                        "pairwise_test": stats_pairwise_test,
                        "multiple_test_correction": bool(
                            st.get("multiple_test_correction", stats_multiple_test_correction)
                        ),
                        "n_blocks": st["n_blocks"],
                        "friedman_q": st["friedman_q"],
                        "friedman_p_perm": st["friedman_p"],
                        "friedman_significant_alpha": bool(
                            np.isfinite(st["friedman_p"]) and st["friedman_p"] < stats_alpha
                        ),
                        "alpha": stats_alpha,
                        "n_permutations": stats_permutations,
                    }
                )
                for row in st["pairwise"]:
                    pairwise_rows.append(
                        {
                            "panel_group": panel_key,
                            "panel_title": _panel_title(panel_key),
                            "metric": col,
                            "label": label,
                            "method_a": row["method_a"],
                            "method_b": row["method_b"],
                            "stats_model": st.get("stats_model", stats_model),
                            "pairwise_test": row.get("pairwise_test", stats_pairwise_test),
                            "statistic": row.get("statistic", row.get("w_plus", float("nan"))),
                            "w_plus": row["w_plus"],
                            "n_nonzero_pairs": row["n_nonzero_pairs"],
                            "n_pos_pairs": row.get("n_pos_pairs", float("nan")),
                            "n_neg_pairs": row.get("n_neg_pairs", float("nan")),
                            "n_zero_pairs": row.get("n_zero_pairs", float("nan")),
                            "median_diff_a_minus_b": row["median_diff_a_minus_b"],
                            "estimate": row.get("estimate", row["median_diff_a_minus_b"]),
                            "std_error": row.get("std_error", float("nan")),
                            "model_family": row.get("model_family", ""),
                            "p_raw": row["p_raw"],
                            "p_adj_bh": row["p_adj_bh"],
                            "p_value_used": row.get("p_value_used", row["p_adj_bh"]),
                            "p_adjustment": row.get("p_adjustment", "bh"),
                            "multiple_test_correction": bool(
                                row.get("multiple_test_correction", stats_multiple_test_correction)
                            ),
                            "significant_alpha": row.get("significant_alpha", row["significant_bh"]),
                            "significant_bh_alpha": row["significant_bh"],
                            "alpha": stats_alpha,
                            "n_permutations": stats_permutations,
                        }
                    )
            pd.DataFrame(omnibus_rows).to_csv(out_dir / f"{stem}_stats_omnibus.tsv", sep="\t", index=False)
            pd.DataFrame(pairwise_rows).to_csv(out_dir / f"{stem}_stats_pairwise.tsv", sep="\t", index=False)

        fig = plt.figure(figsize=(9.0, 7.2))
        for i, (panel_key, col, label, axis_df) in enumerate(axis_specs, start=1):
            ax = fig.add_subplot(2, 2, i)
            st_local = axis_stats.get((panel_key, col), {})
            _draw_one_panel(
                ax=ax,
                col=col,
                label=label,
                data_df=axis_df,
                panel_stats_local={col: st_local} if st_local else {},
                shared_ylim=shared_ylim,
            )

        fig.suptitle(title, fontsize=fig09_suptitle_size)
        fig.tight_layout(rect=fig09_layout_rect, w_pad=fig09_layout_w_pad)
        out_path = out_dir / f"{stem}.svg"
        _save_fig09_main(fig, out_path)
        plt.close(fig)
        if fig09_write_individual_panels:
            panel_letters = ["a", "b", "c", "d"]
            for panel_letter, (panel_key, col, label, axis_df) in zip(panel_letters, axis_specs):
                fig_single = plt.figure(figsize=(5.2, 4.6))
                ax_single = fig_single.add_subplot(1, 1, 1)
                st_local = axis_stats.get((panel_key, col), {})
                _draw_one_panel(
                    ax=ax_single,
                    col=col,
                    label=label,
                    data_df=axis_df,
                    panel_stats_local={col: st_local} if st_local else {},
                    shared_ylim=shared_ylim,
                )
                fig_single.tight_layout()
                _save_fig09_main(fig_single, out_dir / f"{stem}_panel_{panel_letter}.svg")
                plt.close(fig_single)

    else:
        stem = _safe_output_stem(_figure_stem(level_override=None))
        panel_stats = _compute_panel_stats_bundle(df, stem=stem, seed_bias=0)
        _draw_panel_set(panels, f"{stem}.svg", title, df, panel_stats)

    log.info(f"Wrote Figure 09 to {out_dir}")
    return 0
