from __future__ import annotations

import copy
from fnmatch import fnmatch
from itertools import combinations
from pathlib import Path
import re
import tempfile
import warnings

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
try:
    import seaborn as sns
except Exception:  # pragma: no cover - optional runtime dependency
    sns = None

from ..log import get_logger
from ..utils.labels import display_community_id, display_method, display_sample_key
from ..utils.paths import figures_root
from . import compute_distances, map_features, project_abundances
from .plot_fig11 import (
    METRIC_LABEL,
    METHOD_COLORS,
    METHOD_ORDER,
    _best_hit_aggregate_projection,
    _bh_adjust,
    _compute_panel_stats_permutation as _compute_base_panel_stats_permutation,
    _draw_pairwise_brackets,
    _normalize_pairwise_test as _normalize_base_pairwise_test,
    _normalize_plot_type,
    _paired_sign_test_exact,
    _paired_signed_rank_permutation,
    _apply_pairwise_adjustments,
    _normalize_projection_mapping_mode,
    _pvalue_stars,
    _prepare_fig11_override_workspace,
)
from .plot_fig19 import _build_distance_matrix, _factor_r2, _pcoa

PHYLO_AGNOSTIC_METRICS = ("jaccard", "weighted_jaccard", "braycurtis")
PHYLO_AWARE_METRICS = ("unifrac", "gunifrac_a0.0", "gunifrac_a0.5", "gunifrac_a1.0")
FIG21_DENOISING_METHODS = ("ASV", "zOTU")
FIG21_CLUSTERING_METHODS = ("OTU", "TIC", "TAC")
METHOD_ALIASES = {
    "asv": "ASV",
    "dada2": "ASV",
    "zotu": "zOTU",
    "unoise3": "zOTU",
    "otu": "OTU",
    "uparse": "OTU",
    "tic": "TIC",
    "tac": "TAC",
}

ENDPOINT_GROUP_ORDER = [
    "denoising_phylogeny_aware",
    "denoising_phylogeny_agnostic",
    "clustering_phylogeny_aware",
    "clustering_phylogeny_agnostic",
]
ENDPOINT_GROUP_LABELS = {
    "denoising_phylogeny_aware": "Denoising | phylo-aware",
    "denoising_phylogeny_agnostic": "Denoising | phylo-agnostic",
    "clustering_phylogeny_aware": "Clustering | phylo-aware",
    "clustering_phylogeny_agnostic": "Clustering | phylo-agnostic",
}
ENDPOINT_GROUP_SHORT_LABELS = {
    "denoising_phylogeny_aware": "Denoising | aware",
    "denoising_phylogeny_agnostic": "Denoising | agnostic",
    "clustering_phylogeny_aware": "Clustering | aware",
    "clustering_phylogeny_agnostic": "Clustering | agnostic",
}
PANEL_A_ENDPOINT_GROUP_ORDER = [
    "denoising_phylogeny_aware",
    "denoising_phylogeny_agnostic",
    "uparse_phylogeny_aware",
    "uparse_phylogeny_agnostic",
    "tactic_phylogeny_aware",
    "tactic_phylogeny_agnostic",
]
PANEL_A_ENDPOINT_GROUP_LABELS = {
    "denoising_phylogeny_aware": "Denoising | phylo-aware",
    "denoising_phylogeny_agnostic": "Denoising | phylo-agnostic",
    "uparse_phylogeny_aware": "UPARSE | phylo-aware",
    "uparse_phylogeny_agnostic": "UPARSE | phylo-agnostic",
    "tactic_phylogeny_aware": "TACTIC | phylo-aware",
    "tactic_phylogeny_agnostic": "TACTIC | phylo-agnostic",
}
PANEL_A_ENDPOINT_GROUP_SHORT_LABELS = {
    "denoising_phylogeny_aware": "Denoising | phylo-aware",
    "denoising_phylogeny_agnostic": "Denoising | phylo-agnostic",
    "uparse_phylogeny_aware": "UPARSE | phylo-aware",
    "uparse_phylogeny_agnostic": "UPARSE | phylo-agnostic",
    "tactic_phylogeny_aware": "TACTIC | phylo-aware",
    "tactic_phylogeny_agnostic": "TACTIC | phylo-agnostic",
}
PAIR_GROUP_PAIRS = [
    ("denoising_phylogeny_aware", "denoising_phylogeny_aware"),
    ("denoising_phylogeny_agnostic", "denoising_phylogeny_agnostic"),
    ("clustering_phylogeny_aware", "clustering_phylogeny_aware"),
    ("clustering_phylogeny_agnostic", "clustering_phylogeny_agnostic"),
    ("denoising_phylogeny_aware", "denoising_phylogeny_agnostic"),
    ("clustering_phylogeny_aware", "clustering_phylogeny_agnostic"),
    ("denoising_phylogeny_aware", "clustering_phylogeny_aware"),
    ("denoising_phylogeny_agnostic", "clustering_phylogeny_agnostic"),
    ("denoising_phylogeny_aware", "clustering_phylogeny_agnostic"),
    ("denoising_phylogeny_agnostic", "clustering_phylogeny_aware"),
]
PAIR_GROUP_COLORS = {
    "denoising_phylogeny_aware__vs__denoising_phylogeny_aware": "#6BAED6",
    "denoising_phylogeny_agnostic__vs__denoising_phylogeny_agnostic": "#FD8D3C",
    "clustering_phylogeny_aware__vs__clustering_phylogeny_aware": "#74C476",
    "clustering_phylogeny_agnostic__vs__clustering_phylogeny_agnostic": "#9E9AC8",
    "denoising_phylogeny_aware__vs__denoising_phylogeny_agnostic": "#E6550D",
    "clustering_phylogeny_aware__vs__clustering_phylogeny_agnostic": "#31A354",
    "denoising_phylogeny_aware__vs__clustering_phylogeny_aware": "#3182BD",
    "denoising_phylogeny_agnostic__vs__clustering_phylogeny_agnostic": "#756BB1",
    "denoising_phylogeny_aware__vs__clustering_phylogeny_agnostic": "#636363",
    "denoising_phylogeny_agnostic__vs__clustering_phylogeny_aware": "#B08D00",
}
PLANNED_PAIR_GROUP_COMPARISONS = [
    (
        "denoising_phylogeny_aware__vs__denoising_phylogeny_agnostic",
        "clustering_phylogeny_aware__vs__clustering_phylogeny_agnostic",
        "metric_family_gap_denoising_vs_clustering",
    ),
    (
        "denoising_phylogeny_aware__vs__clustering_phylogeny_aware",
        "denoising_phylogeny_agnostic__vs__clustering_phylogeny_agnostic",
        "method_family_gap_aware_vs_agnostic",
    ),
]
PAIR_GROUP_ORDER = [f"{a}__vs__{b}" for a, b in PAIR_GROUP_PAIRS]
WITHIN_PAIR_GROUP_ORDER = [
    "denoising_phylogeny_aware__vs__denoising_phylogeny_aware",
    "denoising_phylogeny_agnostic__vs__denoising_phylogeny_agnostic",
    "clustering_phylogeny_aware__vs__clustering_phylogeny_aware",
    "clustering_phylogeny_agnostic__vs__clustering_phylogeny_agnostic",
]
ENDPOINT_GROUP_RANK = {group: idx for idx, group in enumerate(ENDPOINT_GROUP_ORDER)}


def _canonical_method(method: object) -> str | None:
    value = str(method or "").strip()
    if not value:
        return None
    return METHOD_ALIASES.get(value.lower(), value if value in METHOD_ALIASES.values() else None)


def _resolve_method_list(raw: object, available: list[str]) -> list[str]:
    """Return configured canonical methods in their requested order."""
    if not isinstance(raw, list) or not raw:
        return list(available)
    selected: list[str] = []
    for method in raw:
        canonical = _canonical_method(method)
        if canonical in available and canonical not in selected:
            selected.append(canonical)
    return selected


def _resolve_panel_b_method_groups(raw: object) -> tuple[list[str], list[str]]:
    """Resolve configurable denoising/clustering membership for pair-group panels."""
    default = (list(FIG21_DENOISING_METHODS), list(FIG21_CLUSTERING_METHODS))
    if not isinstance(raw, dict):
        return default
    denoising = _resolve_method_list(raw.get("denoising"), list(FIG21_DENOISING_METHODS) + list(FIG21_CLUSTERING_METHODS))
    clustering = _resolve_method_list(raw.get("clustering"), list(FIG21_DENOISING_METHODS) + list(FIG21_CLUSTERING_METHODS))
    if not denoising or not clustering:
        return default
    # A method cannot contribute to both endpoint families.
    clustering = [method for method in clustering if method not in denoising]
    return (denoising, clustering) if clustering else default


def _normalize_panel_f_pairwise_test(raw: object) -> str:
    value = str(raw or "").strip().lower()
    if value in {"signed_rank_permutation", "sign_test_exact", "mean_difference_permutation"}:
        return value
    return "signed_rank_permutation"


def _normalize_panel_f_summary_statistic(raw: object) -> str:
    value = str(raw or "").strip().lower()
    return value if value in {"median", "mean", "both"} else "median"


def _normalize_pairwise_test(raw: object, default: str = "sign_test_exact") -> str:
    """Normalize Figure 21's global paired-test selector.

    Figure 21 additionally supports a paired mean-difference sign-flip test.
    Other figures retain the shared Figure 11 selector and its narrower set of
    accepted values.
    """
    value = str(raw or "").strip().lower()
    if value == "mean_difference_permutation":
        return value
    return _normalize_base_pairwise_test(value, default=default)


def _normalize_permutation_alternative(raw: object) -> str:
    value = str(raw or "").strip().lower().replace("_", "-")
    if value in {"greater", "less", "two-sided", "two.sided", "two sided"}:
        return "two-sided" if value in {"two.sided", "two sided"} else value
    return "two-sided"


def _resolve_panel_f_planned_comparisons(raw: object, methods: list[str]) -> list[tuple[str, str]]:
    """Resolve optional ordered method pairs for Panel F directional tests."""
    if not isinstance(raw, list) or not raw:
        return list(combinations(methods, 2))
    pairs: list[tuple[str, str]] = []
    for pair in raw:
        if not isinstance(pair, list) or len(pair) != 2:
            continue
        method_a = _canonical_method(pair[0])
        method_b = _canonical_method(pair[1])
        if method_a in methods and method_b in methods and method_a != method_b:
            candidate = (method_a, method_b)
            if candidate not in pairs:
                pairs.append(candidate)
    return pairs


def _paired_mean_difference_permutation(
    values_a: np.ndarray,
    values_b: np.ndarray,
    n_perm: int,
    seed: int,
    alternative: str,
    bootstrap_resamples: int,
    ci_level: float,
    practical_difference_delta: float = 0.0,
) -> dict[str, float | int]:
    """Sign-flip test and bootstrap CI for paired mean differences.

    For one-sided alternatives, ``practical_difference_delta`` specifies the
    minimum effect magnitude required to reject the null. The returned mean
    difference and confidence interval remain on the original distance scale.
    """
    delta = np.asarray(values_a, dtype=float) - np.asarray(values_b, dtype=float)
    delta = delta[np.isfinite(delta)]
    practical_difference_delta = max(0.0, float(practical_difference_delta))
    if delta.size == 0:
        return {
            "statistic": float("nan"), "p_raw": float("nan"), "n_pairs": 0,
            "n_nonzero_pairs": 0, "n_pos_pairs": 0, "n_neg_pairs": 0, "n_zero_pairs": 0,
            "median_delta": float("nan"), "std_error": float("nan"),
            "ci_lower": float("nan"), "ci_upper": float("nan"),
        }

    observed = float(np.mean(delta))
    if alternative == "greater":
        test_delta = delta - practical_difference_delta
    elif alternative == "less":
        test_delta = delta + practical_difference_delta
    else:
        if practical_difference_delta > 0.0:
            raise ValueError(
                "Panel F practical_difference_delta requires a one-sided "
                "mean_difference_permutation test (alternative: greater or less)."
            )
        test_delta = delta
    observed_above_delta = float(np.mean(test_delta))
    n_perm = max(1, int(n_perm))
    rng = np.random.default_rng(seed)
    signs = rng.choice(np.array([-1.0, 1.0]), size=(n_perm, test_delta.size))
    null_means = np.mean(signs * test_delta[np.newaxis, :], axis=1)
    tolerance = np.finfo(float).eps * max(1.0, abs(observed_above_delta))
    if alternative == "greater":
        extreme = int(np.sum(null_means >= observed_above_delta - tolerance))
    elif alternative == "less":
        extreme = int(np.sum(null_means <= observed_above_delta + tolerance))
    else:
        extreme = int(np.sum(np.abs(null_means) >= abs(observed_above_delta) - tolerance))
    p_raw = float((extreme + 1) / (n_perm + 1))

    ci_level = min(0.999, max(0.001, float(ci_level)))
    n_boot = max(1, int(bootstrap_resamples))
    boot_indices = rng.integers(0, delta.size, size=(n_boot, delta.size))
    boot_means = np.mean(delta[boot_indices], axis=1)
    tail = (1.0 - ci_level) / 2.0
    ci_lower, ci_upper = np.quantile(boot_means, [tail, 1.0 - tail])
    return {
        "statistic": observed,
        "mean_difference_above_delta": observed_above_delta,
        "practical_difference_delta": practical_difference_delta,
        "p_raw": p_raw,
        "n_pairs": int(delta.size),
        "n_nonzero_pairs": int(np.sum(delta != 0.0)),
        "n_pos_pairs": int(np.sum(delta > 0.0)),
        "n_neg_pairs": int(np.sum(delta < 0.0)),
        "n_zero_pairs": int(np.sum(delta == 0.0)),
        "median_delta": float(np.median(delta)),
        "std_error": float(np.std(delta, ddof=1) / np.sqrt(delta.size)) if delta.size > 1 else 0.0,
        "ci_lower": float(ci_lower),
        "ci_upper": float(ci_upper),
    }


def _compute_panel_stats_permutation(
    df: pd.DataFrame,
    method_order: list[str],
    n_perm: int,
    seed: int,
    alpha: float,
    pairwise_test: str,
    multiple_test_correction: bool,
) -> dict:
    """Compute Figure 21 paired tests, including mean-difference permutations.

    The omnibus Friedman permutation statistic remains unchanged. When the
    global Figure 21 selector requests mean differences, pairwise tests use
    the arithmetic mean of matched sample-level differences and a two-sided
    sign-flip null distribution.
    """
    normalized_test = _normalize_pairwise_test(pairwise_test)
    if normalized_test != "mean_difference_permutation":
        return _compute_base_panel_stats_permutation(
            df=df,
            method_order=method_order,
            n_perm=n_perm,
            seed=seed,
            alpha=alpha,
            pairwise_test=normalized_test,
            multiple_test_correction=multiple_test_correction,
        )

    # Reuse the established Friedman omnibus calculation, then replace the
    # pairwise rank/sign tests with the requested paired mean-difference test.
    base = _compute_base_panel_stats_permutation(
        df=df,
        method_order=method_order,
        n_perm=n_perm,
        seed=seed,
        alpha=alpha,
        pairwise_test="signed_rank_permutation",
        multiple_test_correction=multiple_test_correction,
    )
    wide = df.pivot_table(index="sample_key", columns="method", values="distance", aggfunc="first")
    methods = [method for method in method_order if method in wide.columns]
    wide = wide[methods].dropna(axis=0, how="any") if methods else pd.DataFrame()
    if wide.shape[0] < 2 or len(methods) < 2:
        return base

    pairs: list[dict] = []
    for pair_idx, (method_a, method_b) in enumerate(combinations(methods, 2), start=1):
        result = _paired_mean_difference_permutation(
            wide[method_a].to_numpy(dtype=float),
            wide[method_b].to_numpy(dtype=float),
            n_perm=n_perm,
            seed=seed + 1000 + pair_idx,
            alternative="two-sided",
            bootstrap_resamples=n_perm,
            ci_level=0.95,
        )
        pairs.append(
            {
                "method_a": method_a,
                "method_b": method_b,
                "pairwise_test": "mean_difference_permutation",
                "alternative": "two-sided",
                "statistic": result["statistic"],
                "mean_diff_a_minus_b": result["statistic"],
                "median_diff_a_minus_b": result["median_delta"],
                "estimate": result["statistic"],
                "std_error": result["std_error"],
                "ci_level": 0.95,
                "ci_lower": result["ci_lower"],
                "ci_upper": result["ci_upper"],
                "bootstrap_resamples": n_perm,
                "p_raw": result["p_raw"],
                "n_pairs": result["n_pairs"],
                "n_nonzero_pairs": result["n_nonzero_pairs"],
                "n_pos_pairs": result["n_pos_pairs"],
                "n_neg_pairs": result["n_neg_pairs"],
                "n_zero_pairs": result["n_zero_pairs"],
                "w_plus": float("nan"),
                "model_family": "sign_flip_permutation",
            }
        )

    adjusted_pairs, significant_pairs = _apply_pairwise_adjustments(
        pairs=pairs,
        alpha=alpha,
        multiple_test_correction=multiple_test_correction,
    )
    base["pairwise"] = adjusted_pairs
    base["pairwise_sig"] = significant_pairs
    return base


METRIC_COLORS = {
    "jaccard": "#000000",
    "weighted_jaccard": "#0072B2",
    "braycurtis": "#D55E00",
    "unifrac": "#009E73",
    "gunifrac_a0.0": "#56B4E9",
    "gunifrac_a0.5": "#CC79A7",
    "gunifrac_a1.0": "#E69F00",
}
METRIC_MARKERS = {
    "jaccard": "o",
    "weighted_jaccard": "s",
    "braycurtis": "^",
    "unifrac": "v",
    "gunifrac_a0.0": "D",
    "gunifrac_a0.5": "P",
    "gunifrac_a1.0": "X",
}
SAMPLE_STATUS_MARKERS = {
    "multiple_v4_reference": "o",
    "single_v4_reference": "s",
}
SAMPLE_STATUS_LABELS = {
    "multiple_v4_reference": "Multiple-V4 reference",
    "single_v4_reference": "Single-V4 reference",
}
DISTANCE_DEFLATION_FAMILY_ORDER = ["phylogeny_agnostic", "phylogeny_aware"]
DISTANCE_DEFLATION_FAMILY_LABELS = {
    "phylogeny_agnostic": "Phylo-agnostic",
    "phylogeny_aware": "Phylo-aware",
}


def _is_excluded(cid: str, patterns: list[str]) -> bool:
    for pat in patterns:
        p = str(pat).strip()
        if not p:
            continue
        if fnmatch(cid, p) or (p in cid):
            return True
    return False


def _load_fig21_distances(cfg, log) -> pd.DataFrame:
    outdir = Path(cfg.paths.outdir)
    long_p = outdir / "distances_all_long.tsv"

    target_overrides = getattr(cfg.figures, "fig21_target_level_by_method", {})
    min_ident_overrides = getattr(cfg.figures, "fig21_min_effective_pident_by_method", {})
    projection_mode = _normalize_projection_mapping_mode(
        getattr(cfg.figures, "fig21_projection_mapping_mode", "standard")
    )
    use_overrides = (
        isinstance(target_overrides, dict)
        and bool(target_overrides)
    ) or (
        isinstance(min_ident_overrides, dict)
        and bool(min_ident_overrides)
    ) or (projection_mode != "standard")

    if not use_overrides:
        if not long_p.exists():
            log.info("Figure 21: distances_all_long.tsv not found; run compute-distances first.")
            return pd.DataFrame()
        return pd.read_csv(long_p, sep="\t")

    log.info("Figure 21: using figure-specific mapping/projection overrides; recomputing temporary distances.")
    with tempfile.TemporaryDirectory(prefix="tacticbench_fig21_override_") as td:
        tmp_outdir = Path(td) / "out"
        ok = _prepare_fig11_override_workspace(src_outdir=outdir, tmp_outdir=tmp_outdir)
        if not ok:
            log.info("Figure 21: failed to prepare temporary workspace for figure-specific overrides.")
            return pd.DataFrame()

        tmp_cfg = copy.deepcopy(cfg)
        tmp_cfg.paths.outdir = tmp_outdir
        if isinstance(target_overrides, dict) and target_overrides:
            tmp_cfg.mapping.target_level_by_method = dict(target_overrides)
        if isinstance(min_ident_overrides, dict) and min_ident_overrides:
            tmp_cfg.mapping.min_effective_pident_by_method = dict(min_ident_overrides)

        map_features.run(tmp_cfg)
        if projection_mode == "best_hit_aggregate":
            _best_hit_aggregate_projection(tmp_cfg, tmp_outdir, log, figure_label="Figure 21")
        else:
            project_abundances.run(tmp_cfg)
        compute_distances.run(tmp_cfg)

        tmp_long_p = tmp_outdir / "distances_all_long.tsv"
        if not tmp_long_p.exists():
            log.info("Figure 21: temporary distances were not generated.")
            return pd.DataFrame()
        return pd.read_csv(tmp_long_p, sep="\t")


def _metric_label(metric: str) -> str:
    label = METRIC_LABEL.get(str(metric), str(metric))
    return str(label).replace("alpha=", "α=")


def _safe_filename(value: str) -> str:
    out = re.sub(r"[^A-Za-z0-9._-]+", "__", str(value)).strip("_")
    return out or "sample"


def _load_sample_v4_status(benchmark_outdir: Path) -> pd.DataFrame:
    qc_p = Path(benchmark_outdir) / "qc" / "multi_v4_species_summary.tsv"
    if not qc_p.exists():
        return pd.DataFrame()
    try:
        df = pd.read_csv(qc_p, sep="\t")
    except Exception:
        return pd.DataFrame()
    if df.empty or "community_id" not in df.columns:
        return pd.DataFrame()

    sample_id_col = None
    for candidate in ("selected_sample_id", "sample_id"):
        if candidate in df.columns:
            sample_id_col = candidate
            break
    if sample_id_col is None:
        return pd.DataFrame()

    out = df.copy()
    out["community_id"] = out["community_id"].astype(str)
    out["sample_id"] = out[sample_id_col].astype(str)
    out["sample_key"] = out["community_id"] + "|" + out["sample_id"]
    if "has_multi_v4_species" in out.columns:
        has_multi = out["has_multi_v4_species"].astype(bool)
    else:
        has_multi = pd.Series(False, index=out.index)
    out["v4_reference_status"] = np.where(has_multi, "multiple_v4_reference", "single_v4_reference")
    cols = [
        "community_id",
        "sample_id",
        "sample_key",
        "v4_reference_status",
        "has_multi_v4_species",
        "n_species_total",
        "n_species_with_multiple_v4_unique_sequences",
        "multi_v4_species",
        "multi_v4_species_abundance_sum",
        "multi_v4_species_abundance_fraction",
        "multi_v4_species_abundance_percent",
        "source_sequences_file",
    ]
    return out[[c for c in cols if c in out.columns]].drop_duplicates(subset=["sample_key"], keep="first")


def _attach_sample_v4_status(sample_summary: pd.DataFrame, benchmark_outdir: Path) -> pd.DataFrame:
    if sample_summary.empty:
        return sample_summary
    status = _load_sample_v4_status(benchmark_outdir)
    out = sample_summary.copy()
    if status.empty:
        out["v4_reference_status"] = "single_v4_reference"
        out["has_multi_v4_species"] = False
        return out

    merged = out.merge(
        status[["community_id", "sample_id", "sample_key", "v4_reference_status", "has_multi_v4_species"]],
        on=["community_id", "sample_id", "sample_key"],
        how="left",
    )
    merged["v4_reference_status"] = merged["v4_reference_status"].fillna("single_v4_reference")
    merged["has_multi_v4_species"] = merged["has_multi_v4_species"].fillna(False).astype(bool)
    return merged


def _metric_family_lists(cfg, metrics: list[str]) -> tuple[list[str], list[str]]:
    metric_set = {str(m) for m in metrics}
    raw_agnostic = list(getattr(cfg.figures, "fig21_agnostic_metrics", []) or [])
    raw_aware = list(getattr(cfg.figures, "fig21_aware_metrics", []) or [])

    agnostic = [str(m) for m in raw_agnostic if str(m) in metric_set]
    aware = [str(m) for m in raw_aware if str(m) in metric_set]
    if not agnostic:
        agnostic = [m for m in PHYLO_AGNOSTIC_METRICS if m in metric_set]
    if not aware:
        aware = [m for m in PHYLO_AWARE_METRICS if m in metric_set]
    return agnostic, aware


def _metric_family_sensitivity_points(
    long: pd.DataFrame,
    method_order: list[str],
    agnostic_metrics: list[str],
    aware_metrics: list[str],
) -> pd.DataFrame:
    cols = ["community_id", "sample_id", "sample_key", "method", "metric", "distance"]
    if any(c not in long.columns for c in cols):
        return pd.DataFrame()
    if not agnostic_metrics or not aware_metrics:
        return pd.DataFrame()

    metric_order = agnostic_metrics + aware_metrics
    sub = long[long["metric"].astype(str).isin(metric_order)].copy()
    sub = sub[sub["method"].astype(str).isin(method_order)].copy()
    if sub.empty:
        return pd.DataFrame()

    sub["metric_family"] = np.where(
        sub["metric"].astype(str).isin(agnostic_metrics),
        "phylogeny_agnostic",
        "phylogeny_aware",
    )
    means = (
        sub.groupby(["community_id", "sample_id", "sample_key", "method", "metric_family"], sort=False)["distance"]
        .mean()
        .reset_index()
    )
    wide = means.pivot_table(
        index=["community_id", "sample_id", "sample_key", "method"],
        columns="metric_family",
        values="distance",
        aggfunc="first",
    ).reset_index()
    if "phylogeny_agnostic" not in wide.columns or "phylogeny_aware" not in wide.columns:
        return pd.DataFrame()

    wide = wide.dropna(subset=["phylogeny_agnostic", "phylogeny_aware"]).copy()
    if wide.empty:
        return pd.DataFrame()

    wide["sensitivity"] = wide["phylogeny_agnostic"] - wide["phylogeny_aware"]
    wide["community_label"] = wide["community_id"].map(display_community_id)
    wide["sample_label"] = wide["sample_key"].map(display_sample_key)
    wide["method_label"] = wide["method"].map(display_method)
    wide["agnostic_metrics"] = ",".join(agnostic_metrics)
    wide["aware_metrics"] = ",".join(aware_metrics)
    return wide[
        [
            "community_id",
            "sample_id",
            "sample_key",
            "community_label",
            "sample_label",
            "method",
            "method_label",
            "phylogeny_agnostic",
            "phylogeny_aware",
            "sensitivity",
            "agnostic_metrics",
            "aware_metrics",
        ]
    ]


def _distance_deflation_points(
    endpoints_df: pd.DataFrame,
    agnostic_metrics: list[str],
    aware_metrics: list[str],
) -> pd.DataFrame:
    cols = ["community_id", "sample_id", "sample_key", "method", "metric", "distance"]
    if endpoints_df.empty or any(c not in endpoints_df.columns for c in cols):
        return pd.DataFrame()

    metric_order = list(agnostic_metrics) + list(aware_metrics)
    sub = endpoints_df[endpoints_df["metric"].astype(str).isin(metric_order)].copy()
    sub = sub[sub["method"].astype(str).isin(set(FIG21_DENOISING_METHODS) | set(FIG21_CLUSTERING_METHODS))].copy()
    if sub.empty:
        return pd.DataFrame()

    sub["method_family"] = np.where(
        sub["method"].astype(str).isin(FIG21_DENOISING_METHODS),
        "denoising",
        "clustering",
    )
    sub["metric_family"] = np.where(
        sub["metric"].astype(str).isin(agnostic_metrics),
        "phylogeny_agnostic",
        "phylogeny_aware",
    )
    grouped = (
        sub.groupby(
            ["community_id", "sample_id", "sample_key", "metric_family", "method_family"],
            sort=False,
        )["distance"]
        .mean()
        .reset_index()
    )
    wide = grouped.pivot_table(
        index=["community_id", "sample_id", "sample_key", "metric_family"],
        columns="method_family",
        values="distance",
        aggfunc="first",
    ).reset_index()
    if "denoising" not in wide.columns or "clustering" not in wide.columns:
        return pd.DataFrame()

    wide = wide.dropna(subset=["denoising", "clustering"]).copy()
    if wide.empty:
        return pd.DataFrame()
    wide["delta_denoising_minus_clustering"] = wide["denoising"] - wide["clustering"]
    wide["metric_family_label"] = wide["metric_family"].map(DISTANCE_DEFLATION_FAMILY_LABELS).fillna(
        wide["metric_family"].astype(str)
    )
    wide["community_label"] = wide["community_id"].map(display_community_id)
    wide["sample_label"] = wide["sample_key"].map(display_sample_key)
    wide["agnostic_metrics"] = ",".join(agnostic_metrics)
    wide["aware_metrics"] = ",".join(aware_metrics)
    return wide[
        [
            "community_id",
            "sample_id",
            "sample_key",
            "community_label",
            "sample_label",
            "metric_family",
            "metric_family_label",
            "denoising",
            "clustering",
            "delta_denoising_minus_clustering",
            "agnostic_metrics",
            "aware_metrics",
        ]
    ].sort_values(["metric_family", "community_id", "sample_id", "sample_key"])


def _method_summary(sens: pd.DataFrame) -> pd.DataFrame:
    if sens.empty:
        return pd.DataFrame()
    rows = []
    for method, g in sens.groupby("method", sort=False):
        vals = pd.to_numeric(g["sensitivity"], errors="coerce").dropna().to_numpy(dtype=float)
        if vals.size == 0:
            continue
        rows.append(
            {
                "method": method,
                "method_label": display_method(method),
                "n": int(vals.size),
                "mean_sensitivity": float(np.mean(vals)),
                "median_sensitivity": float(np.median(vals)),
                "sd_sensitivity": float(np.std(vals, ddof=1)) if vals.size > 1 else 0.0,
                "min_sensitivity": float(np.min(vals)),
                "max_sensitivity": float(np.max(vals)),
                "fraction_positive": float(np.mean(vals > 0.0)),
            }
        )
    return pd.DataFrame(rows)


def _write_sensitivity_stats(
    sens: pd.DataFrame,
    out_dir: Path,
    method_order: list[str],
    stats_enabled: bool,
    stats_permutations: int,
    stats_seed: int,
    stats_alpha: float,
    stats_pairwise_test: str,
    stats_multiple_test_correction: bool,
) -> dict | None:
    if not stats_enabled or sens.empty:
        return None
    stats_df = sens.rename(columns={"sensitivity": "distance"}).copy()
    st = _compute_panel_stats_permutation(
        df=stats_df,
        method_order=method_order,
        n_perm=stats_permutations,
        seed=stats_seed + 21001,
        alpha=stats_alpha,
        pairwise_test=stats_pairwise_test,
        multiple_test_correction=stats_multiple_test_correction,
    )
    pd.DataFrame(
        [
            {
                "metric": "metric_family_sensitivity",
                "label": "Metric-family sensitivity",
                "stats_model": st.get("stats_model", "paired_permutation"),
                "pairwise_test": stats_pairwise_test,
                "multiple_test_correction": bool(
                    st.get("multiple_test_correction", stats_multiple_test_correction)
                ),
                "n_blocks": st.get("n_blocks", 0),
                "omnibus_stat": st.get("omnibus_stat", np.nan),
                "omnibus_p": st.get("omnibus_p", np.nan),
                "omnibus_significant_alpha": bool(
                    np.isfinite(float(st.get("omnibus_p", np.nan)))
                    and float(st.get("omnibus_p", np.nan)) < stats_alpha
                ),
                "alpha": stats_alpha,
                "n_permutations": stats_permutations,
            }
        ]
    ).to_csv(out_dir / "sensitivity_by_method_stats_omnibus.tsv", sep="\t", index=False)

    pairwise_rows: list[dict] = []
    for row in st.get("pairwise", []):
        pairwise_rows.append(
            {
                "metric": "metric_family_sensitivity",
                "label": "Metric-family sensitivity",
                "method_a": row.get("method_a", ""),
                "method_b": row.get("method_b", ""),
                "stats_model": st.get("stats_model", "paired_permutation"),
                "pairwise_test": row.get("pairwise_test", stats_pairwise_test),
                "statistic": row.get("statistic", row.get("w_plus", np.nan)),
                "w_plus": row.get("w_plus", np.nan),
                "n_nonzero_pairs": row.get("n_nonzero_pairs", np.nan),
                "n_pos_pairs": row.get("n_pos_pairs", np.nan),
                "n_neg_pairs": row.get("n_neg_pairs", np.nan),
                "n_zero_pairs": row.get("n_zero_pairs", np.nan),
                "median_diff_a_minus_b": row.get("median_diff_a_minus_b", np.nan),
                "estimate": row.get("estimate", row.get("median_diff_a_minus_b", np.nan)),
                "std_error": row.get("std_error", np.nan),
                "model_family": row.get("model_family", ""),
                "p_raw": row.get("p_raw", np.nan),
                "p_adj_bh": row.get("p_adj_bh", np.nan),
                "p_value_used": row.get("p_value_used", row.get("p_adj_bh", np.nan)),
                "p_adjustment": row.get("p_adjustment", "bh"),
                "multiple_test_correction": bool(
                    row.get("multiple_test_correction", stats_multiple_test_correction)
                ),
                "significant_alpha": row.get("significant_alpha", False),
                "significant_bh_alpha": row.get("significant_bh", False),
                "alpha": stats_alpha,
                "n_permutations": stats_permutations,
            }
        )
    pd.DataFrame(pairwise_rows).to_csv(out_dir / "sensitivity_by_method_stats_pairwise.tsv", sep="\t", index=False)
    return st


def _plot_sensitivity_by_method(
    sens: pd.DataFrame,
    out_dir: Path,
    method_order: list[str],
    title: str,
    plot_type: str,
    point_size: float,
    font_size: float,
    agnostic_metrics: list[str],
    aware_metrics: list[str],
    violin_bw_method: float | str | None,
    stats: dict | None,
    stats_enabled: bool,
    stats_show_pairwise: bool,
    stats_alpha: float,
) -> None:
    if sens.empty:
        return

    fig = plt.figure(figsize=(9.8, 5.7))
    ax = fig.add_subplot(1, 1, 1)
    positions = np.arange(1, len(method_order) + 1, dtype=float)
    method_data: list[tuple[int, str, np.ndarray]] = []
    for i, method in enumerate(method_order):
        vals = pd.to_numeric(sens.loc[sens["method"] == method, "sensitivity"], errors="coerce").dropna().to_numpy()
        if vals.size:
            method_data.append((i, method, vals))

    if plot_type == "boxplot":
        for i, method, vals in method_data:
            ax.boxplot(
                [vals],
                positions=[positions[i]],
                widths=0.55,
                patch_artist=True,
                boxprops=dict(facecolor=METHOD_COLORS.get(method, "#999999"), alpha=0.40, edgecolor="#333333"),
                medianprops=dict(color="#111111", linewidth=1.3),
                whiskerprops=dict(color="#333333", linewidth=0.9),
                capprops=dict(color="#333333", linewidth=0.9),
            )
    else:
        violin_vals: list[np.ndarray] = []
        violin_pos: list[float] = []
        violin_methods: list[str] = []
        for i, method, vals in method_data:
            if vals.size >= 2:
                violin_vals.append(vals)
                violin_pos.append(float(positions[i]))
                violin_methods.append(method)
        if violin_vals:
            vp = ax.violinplot(
                violin_vals,
                positions=violin_pos,
                widths=0.78,
                showmeans=False,
                showmedians=True,
                showextrema=False,
                bw_method=violin_bw_method,
            )
            for body, method in zip(vp["bodies"], violin_methods):
                body.set_facecolor(METHOD_COLORS.get(method, "#999999"))
                body.set_edgecolor("#333333")
                body.set_alpha(0.5)
            if "cmedians" in vp:
                vp["cmedians"].set_color("#111111")
                vp["cmedians"].set_linewidth(1.2)

    rng = np.random.default_rng(21021)
    for i, method, vals in method_data:
        xj = positions[i] + rng.normal(0.0, 0.04, size=vals.size)
        ax.scatter(
            xj,
            vals,
            s=point_size,
            alpha=0.45,
            color=METHOD_COLORS.get(method, "#666666"),
            edgecolors="none",
            zorder=3,
        )

    finite_vals = pd.to_numeric(sens["sensitivity"], errors="coerce").dropna().to_numpy(dtype=float)
    if finite_vals.size:
        y_min = float(np.min(finite_vals))
        y_max = float(np.max(finite_vals))
        y_range = y_max - y_min
        if y_range <= 0.0:
            y_range = max(0.05, abs(y_max) * 0.25)
        ax.set_ylim(y_min - 0.18 * y_range, y_max + 0.35 * y_range)

    tick_size = max(6.0, font_size - 2.0)
    stats_text_size = max(6.0, font_size - 3.0)
    pairwise_label_size = max(6.0, font_size - 3.5)
    pairwise_star_size = max(6.0, font_size - 2.8)

    ax.axhline(0.0, color="#333333", linewidth=0.8, linestyle="--", alpha=0.7)
    ax.set_xticks(positions)
    ax.set_xticklabels([display_method(m) for m in method_order], rotation=0, fontsize=tick_size)
    ax.set_xlabel("Method", fontsize=font_size)
    ax.set_ylabel("Mean phylogeny-agnostic distance - mean phylogeny-aware distance", fontsize=font_size)
    ax.tick_params(axis="y", labelsize=tick_size)
    ax.grid(True, axis="y", alpha=0.25, linewidth=0.6)
    ax.set_title(title, fontsize=font_size)

    agnostic_label = ", ".join(_metric_label(m) for m in agnostic_metrics)
    aware_label = ", ".join(_metric_label(m) for m in aware_metrics)
    ax.text(
        0.02,
        0.03,
        f"Agnostic: {agnostic_label}\nPhylogeny-aware: {aware_label}",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=stats_text_size,
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="#bbbbbb", alpha=0.85),
    )

    if stats_enabled and isinstance(stats, dict):
        if stats_show_pairwise:
            _draw_pairwise_brackets(
                ax=ax,
                pairwise_rows=list(stats.get("pairwise_sig", [])),
                method_order=method_order,
                alpha=stats_alpha,
                label_fontsize=pairwise_label_size,
                star_fontsize=pairwise_star_size,
            )
        op = float(stats.get("omnibus_p", np.nan))
        os = float(stats.get("omnibus_stat", np.nan))
        ax.text(
            0.98,
            0.03,
            f"paired_permutation\nomnibus={os:.3g}, p={op:.3g}" if np.isfinite(op) else "paired_permutation\nomnibus=NA, p=NA",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=stats_text_size,
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="#bbbbbb", alpha=0.85),
        )

    fig.tight_layout()
    suffix = "violin" if plot_type == "violin" else "boxplot"
    fig.savefig(out_dir / f"metric_family_sensitivity_by_method_{suffix}.svg", format="svg", dpi=300)
    plt.close(fig)


def _draw_directional_sensitivity_panel(
    ax,
    sens: pd.DataFrame,
    method_order: list[str],
    title: str,
    plot_type: str,
    point_size: float,
    font_size: float,
    violin_bw_method: float | str | None,
    stats: dict | None,
    stats_enabled: bool,
    stats_show_pairwise: bool,
    stats_alpha: float,
    sample_v4_multiple_label: str = SAMPLE_STATUS_LABELS["multiple_v4_reference"],
    sample_v4_single_label: str = SAMPLE_STATUS_LABELS["single_v4_reference"],
    x_axis_label: str = "Analysis method",
    y_axis_label: str = "Mean agnostic distance - mean aware distance",
    dot_size: float | None = None,
    dot_opacity: float = 0.48,
    dot_jitter: float = 0.07,
    show_stats_annotations: bool = True,
) -> None:
    fs = _fig21_text_scale(font_size)
    if sens.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, "Directional sensitivity unavailable", ha="center", va="center", transform=ax.transAxes)
        return

    present_methods = set(sens["method"].astype(str))
    order = [m for m in method_order if m in present_methods]
    if not order:
        ax.axis("off")
        return

    positions = np.arange(1, len(order) + 1, dtype=float)
    panel_df = sens.copy()
    if "v4_reference_status" not in panel_df.columns:
        panel_df["v4_reference_status"] = "single_v4_reference"
    panel_df["v4_reference_status"] = panel_df["v4_reference_status"].fillna("single_v4_reference").astype(str)
    panel_df["v4_reference_status"] = panel_df["v4_reference_status"].where(
        panel_df["v4_reference_status"].isin(SAMPLE_STATUS_MARKERS),
        "single_v4_reference",
    )

    method_data: list[tuple[int, str, pd.DataFrame, np.ndarray]] = []
    for i, method in enumerate(order):
        method_df = panel_df.loc[
            panel_df["method"] == method,
            ["sensitivity", "v4_reference_status"],
        ].copy()
        method_df["sensitivity"] = pd.to_numeric(method_df["sensitivity"], errors="coerce")
        method_df = method_df.dropna(subset=["sensitivity"])
        vals = method_df["sensitivity"].to_numpy(dtype=float)
        if vals.size:
            method_data.append((i, method, method_df, vals))

    family_for_method = {
        m: "denoising" if m in FIG21_DENOISING_METHODS else "clustering" if m in FIG21_CLUSTERING_METHODS else ""
        for m in order
    }
    for family, color in [("denoising", "#D9ECF7"), ("clustering", "#E4F4DF")]:
        idxs = [i for i, m in enumerate(order) if family_for_method.get(m) == family]
        if not idxs:
            continue
        ax.axvspan(min(idxs) + 0.5, max(idxs) + 1.5, color=color, alpha=0.35, zorder=0)
        ax.text(
            (min(idxs) + max(idxs)) / 2.0 + 1.0,
            0.98,
            family.title(),
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="top",
            fontsize=fs["stats"],
            color="#333333",
        )

    if plot_type == "boxplot":
        for i, method, _, vals in method_data:
            ax.boxplot(
                [vals],
                positions=[positions[i]],
                widths=0.55,
                patch_artist=True,
                boxprops=dict(facecolor=METHOD_COLORS.get(method, "#999999"), alpha=0.42, edgecolor="#333333"),
                medianprops=dict(color="#111111", linewidth=1.25),
                whiskerprops=dict(color="#333333", linewidth=0.9),
                capprops=dict(color="#333333", linewidth=0.9),
            )
    else:
        violin_vals = [vals for _, _, _, vals in method_data if vals.size >= 2]
        violin_pos = [float(positions[i]) for i, _, _, vals in method_data if vals.size >= 2]
        violin_methods = [method for _, method, _, vals in method_data if vals.size >= 2]
        if violin_vals:
            vp = ax.violinplot(
                violin_vals,
                positions=violin_pos,
                widths=0.78,
                showmeans=False,
                showmedians=True,
                showextrema=False,
                bw_method=violin_bw_method,
            )
            for body, method in zip(vp["bodies"], violin_methods):
                body.set_facecolor(METHOD_COLORS.get(method, "#999999"))
                body.set_edgecolor("#333333")
                body.set_alpha(0.52)
            if "cmedians" in vp:
                vp["cmedians"].set_color("#111111")
                vp["cmedians"].set_linewidth(1.15)

    rng = np.random.default_rng(21119)
    scatter_size = max(1.0, float(dot_size)) if dot_size is not None else point_size * 1.755
    for i, method, method_df, vals in method_data:
        xj = positions[i] + rng.normal(0.0, max(0.0, float(dot_jitter)), size=vals.size)
        status = method_df["v4_reference_status"].astype(str).to_numpy()
        for status_key, marker in SAMPLE_STATUS_MARKERS.items():
            mask = status == status_key
            if not np.any(mask):
                continue
            ax.scatter(
                xj[mask],
                vals[mask],
                s=scatter_size,
                marker=marker,
                alpha=max(0.0, min(1.0, float(dot_opacity))),
                color=METHOD_COLORS.get(method, "#666666"),
                edgecolors="#111111",
                linewidths=0.45,
                zorder=3,
            )

    finite_vals = pd.to_numeric(panel_df["sensitivity"], errors="coerce").dropna().to_numpy(dtype=float)
    if finite_vals.size:
        y_min = float(np.min(finite_vals))
        y_max = float(np.max(finite_vals))
        y_range = y_max - y_min
        if y_range <= 0.0:
            y_range = max(0.05, abs(y_max) * 0.25)
        ax.set_ylim(y_min - 0.18 * y_range, y_max + 0.38 * y_range)

    ax.axhline(0.0, color="#333333", linewidth=0.85, linestyle="--", alpha=0.75)
    ax.set_xticks(positions)
    ax.set_xticklabels([display_method(m) for m in order], rotation=0, fontsize=fs["tick"])
    ax.set_xlabel(x_axis_label, fontsize=fs["base"])
    ax.set_ylabel(y_axis_label, fontsize=fs["base"])
    ax.tick_params(axis="y", labelsize=fs["tick"])
    ax.grid(True, axis="y", alpha=0.25, linewidth=0.6)
    ax.set_title(title, fontsize=fs["base"])

    if stats_enabled and show_stats_annotations and isinstance(stats, dict):
        _draw_pairwise_brackets(
            ax=ax,
            # Panel D has no pre-specified contrast subset: show the adjusted
            # result for every selected method pair, as in Panel F.
            pairwise_rows=list(stats.get("pairwise", [])),
            method_order=order,
            alpha=stats_alpha,
            label_fontsize=fs["pair_label"],
            star_fontsize=fs["pair_star"],
        )

    status_handles = [
        Line2D(
            [0],
            [0],
            marker=SAMPLE_STATUS_MARKERS["multiple_v4_reference"],
            linestyle="None",
            markerfacecolor="#666666",
            markeredgecolor="#111111",
            markersize=9.75,
            label=sample_v4_multiple_label,
        ),
        Line2D(
            [0],
            [0],
            marker=SAMPLE_STATUS_MARKERS["single_v4_reference"],
            linestyle="None",
            markerfacecolor="#666666",
            markeredgecolor="#111111",
            markersize=9.75,
            label=sample_v4_single_label,
        ),
    ]


def _method_metric_family_distance_points(
    endpoints_df: pd.DataFrame,
    agnostic_metrics: list[str],
    aware_metrics: list[str],
    methods: list[str],
) -> pd.DataFrame:
    """Mean distance per sample, method, and metric family for Panel F."""
    cols = ["community_id", "sample_id", "sample_key", "method", "metric", "distance"]
    if endpoints_df.empty or any(col not in endpoints_df.columns for col in cols):
        return pd.DataFrame()
    sub = endpoints_df[
        endpoints_df["method"].astype(str).isin(methods)
        & endpoints_df["metric"].astype(str).isin(list(agnostic_metrics) + list(aware_metrics))
    ].copy()
    if sub.empty:
        return pd.DataFrame()
    sub["metric_family"] = np.where(
        sub["metric"].astype(str).isin(agnostic_metrics), "phylogeny_agnostic", "phylogeny_aware"
    )
    out = (
        sub.groupby(["community_id", "sample_id", "sample_key", "method", "metric_family"], sort=False)["distance"]
        .mean()
        .reset_index()
    )
    out["community_label"] = out["community_id"].map(display_community_id)
    out["sample_label"] = out["sample_key"].map(display_sample_key)
    out["method_label"] = out["method"].map(display_method)
    out["metric_family_label"] = out["metric_family"].map(DISTANCE_DEFLATION_FAMILY_LABELS).fillna(out["metric_family"])
    return out


def _write_method_metric_family_stats(
    points: pd.DataFrame,
    methods: list[str],
    out_dir: Path,
    stats_enabled: bool,
    stats_permutations: int,
    stats_seed: int,
    stats_alpha: float,
    stats_pairwise_test: str,
    stats_multiple_test_correction: bool,
    panel_f_stats_pairwise_test: str,
    panel_f_stats_alternative: str,
    panel_f_stats_practical_difference_delta: float,
    panel_f_stats_bootstrap_resamples: int,
    panel_f_stats_ci_level: float,
    panel_f_stats_planned_comparisons: list[list[str]],
) -> list[dict]:
    if not stats_enabled or points.empty:
        pd.DataFrame().to_csv(out_dir / "distance_by_method_metric_family_stats.tsv", sep="\t", index=False)
        return []
    rows: list[dict] = []
    panel_f_test = _normalize_panel_f_pairwise_test(panel_f_stats_pairwise_test)
    alternative = _normalize_permutation_alternative(panel_f_stats_alternative)
    for idx, family in enumerate(DISTANCE_DEFLATION_FAMILY_ORDER, start=1):
        df = points.loc[points["metric_family"].astype(str) == family].copy()
        if df.empty:
            continue
        if panel_f_test != "mean_difference_permutation":
            st = _compute_panel_stats_permutation(
                df=df.rename(columns={"distance": "distance"}),
                method_order=methods,
                n_perm=stats_permutations,
                seed=stats_seed + 22100 + idx,
                alpha=stats_alpha,
                pairwise_test=panel_f_test,
                multiple_test_correction=stats_multiple_test_correction,
            )
            for rec in st.get("pairwise", []):
                row = dict(rec)
                row["alternative"] = "two-sided"
                row["metric_family"] = family
                row["metric_family_label"] = DISTANCE_DEFLATION_FAMILY_LABELS.get(family, family)
                rows.append(row)
            continue

        wide = df.pivot_table(index="sample_key", columns="method", values="distance", aggfunc="first")
        active_methods = [method for method in methods if method in wide.columns]
        family_rows: list[dict] = []
        for pair_idx, (method_a, method_b) in enumerate(
            _resolve_panel_f_planned_comparisons(panel_f_stats_planned_comparisons, active_methods), start=1
        ):
            paired = wide[[method_a, method_b]].dropna(axis=0, how="any")
            result = _paired_mean_difference_permutation(
                paired[method_a].to_numpy(dtype=float),
                paired[method_b].to_numpy(dtype=float),
                n_perm=stats_permutations,
                seed=stats_seed + 22100 + (idx * 100) + pair_idx,
                alternative=alternative,
                bootstrap_resamples=panel_f_stats_bootstrap_resamples,
                ci_level=panel_f_stats_ci_level,
                practical_difference_delta=panel_f_stats_practical_difference_delta,
            )
            family_rows.append(
                {
                    "method_a": method_a,
                    "method_b": method_b,
                    "pairwise_test": "mean_difference_permutation",
                    "alternative": alternative,
                    "statistic": result["statistic"],
                    "mean_diff_a_minus_b": result["statistic"],
                    "mean_difference_above_delta": result["mean_difference_above_delta"],
                    "practical_difference_delta": result["practical_difference_delta"],
                    "median_diff_a_minus_b": result["median_delta"],
                    "estimate": result["statistic"],
                    "std_error": result["std_error"],
                    "ci_level": panel_f_stats_ci_level,
                    "ci_lower": result["ci_lower"],
                    "ci_upper": result["ci_upper"],
                    "bootstrap_resamples": panel_f_stats_bootstrap_resamples,
                    "p_raw": result["p_raw"],
                    "n_pairs": result["n_pairs"],
                    "n_nonzero_pairs": result["n_nonzero_pairs"],
                    "n_pos_pairs": result["n_pos_pairs"],
                    "n_neg_pairs": result["n_neg_pairs"],
                    "n_zero_pairs": result["n_zero_pairs"],
                    "model_family": "sign_flip_permutation",
                    "metric_family": family,
                    "metric_family_label": DISTANCE_DEFLATION_FAMILY_LABELS.get(family, family),
                }
            )
        raw_p = [float(row["p_raw"]) for row in family_rows]
        adjusted_p = _bh_adjust(raw_p)
        for row, p_adjusted in zip(family_rows, adjusted_p):
            p_used = p_adjusted if stats_multiple_test_correction else float(row["p_raw"])
            row["p_adj_bh"] = p_adjusted
            row["p_value_used"] = p_used
            row["multiple_test_correction"] = bool(stats_multiple_test_correction)
            row["p_adjustment"] = "bh" if stats_multiple_test_correction else "none"
            row["significant_bh"] = bool(np.isfinite(p_adjusted) and p_adjusted < stats_alpha)
            row["significant_alpha"] = bool(np.isfinite(p_used) and p_used < stats_alpha)
        rows.extend(family_rows)
    pd.DataFrame(rows).to_csv(out_dir / "distance_by_method_metric_family_stats.tsv", sep="\t", index=False)
    return rows
    ax.legend(
        handles=status_handles,
        title="Sample V4 status",
        loc="upper left",
        frameon=True,
        framealpha=0.92,
        facecolor="white",
        edgecolor="#BBBBBB",
        fontsize=fs["legend"],
        title_fontsize=fs["legend_title"],
    )

    if stats_enabled and isinstance(stats, dict):
        if stats_show_pairwise:
            _draw_pairwise_brackets(
                ax=ax,
                pairwise_rows=list(stats.get("pairwise_sig", [])),
                method_order=order,
                alpha=stats_alpha,
                label_fontsize=fs["pair_label"],
                star_fontsize=fs["pair_star"],
            )
        op = float(stats.get("omnibus_p", np.nan))
        os = float(stats.get("omnibus_stat", np.nan))
        ax.text(
            0.98,
            0.04,
            f"paired permutation\nomnibus={os:.3g}, p={op:.3g}" if np.isfinite(op) else "paired permutation\nomnibus=NA, p=NA",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=fs["stats"],
            bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor="#bbbbbb", alpha=0.88),
        )


def _draw_distance_deflation_panel(
    ax,
    deflation: pd.DataFrame,
    stats_rows: list[dict],
    title: str,
    plot_type: str,
    point_size: float,
    font_size: float,
    violin_bw_method: float | str | None,
    stats_enabled: bool,
    stats_alpha: float,
    sample_v4_multiple_label: str = SAMPLE_STATUS_LABELS["multiple_v4_reference"],
    sample_v4_single_label: str = SAMPLE_STATUS_LABELS["single_v4_reference"],
    x_axis_label: str = "Metric family | method family",
    y_axis_label: str = "Mean distance to expected composition",
    dot_size: float | None = None,
    dot_opacity: float = 0.75,
    dot_jitter: float = 0.06,
) -> None:
    fs = _fig21_text_scale(font_size)
    if deflation.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, "Distance-deflation contrast unavailable", ha="center", va="center", transform=ax.transAxes)
        return

    panel_df = deflation.copy()
    if "v4_reference_status" not in panel_df.columns:
        panel_df["v4_reference_status"] = "single_v4_reference"
    panel_df["v4_reference_status"] = panel_df["v4_reference_status"].fillna("single_v4_reference").astype(str)
    panel_df["v4_reference_status"] = panel_df["v4_reference_status"].where(
        panel_df["v4_reference_status"].isin(SAMPLE_STATUS_MARKERS),
        "single_v4_reference",
    )
    order = [f for f in DISTANCE_DEFLATION_FAMILY_ORDER if f in set(panel_df["metric_family"].astype(str))]
    if not order:
        ax.axis("off")
        return

    family_colors = {
        "phylogeny_agnostic": "#FFF1E6",
        "phylogeny_aware": "#EAF4FB",
    }
    method_colors = {
        "denoising": "#0072B2",
        "clustering": "#009E73",
    }
    method_labels = {
        "denoising": "Denoising",
        "clustering": "Clustering",
    }
    positions: dict[tuple[str, str], float] = {}
    x_ticks: list[float] = []
    x_labels: list[str] = []
    for i, family in enumerate(order):
        base = 1.0 + (i * 2.75)
        positions[(family, "denoising")] = base
        positions[(family, "clustering")] = base + 1.0
        x_ticks.extend([base, base + 1.0])
        fam_label = DISTANCE_DEFLATION_FAMILY_LABELS.get(family, family)
        x_labels.extend([f"{fam_label}\nDenoising", f"{fam_label}\nClustering"])
        ax.axvspan(base - 0.45, base + 1.45, color=family_colors.get(family, "#f2f2f2"), alpha=0.55, zorder=0)

    value_rows: list[dict] = []
    line_rows: list[dict] = []
    for _, row in panel_df.iterrows():
        family = str(row.get("metric_family", ""))
        if family not in order:
            continue
        den = pd.to_numeric(pd.Series([row.get("denoising", np.nan)]), errors="coerce").iloc[0]
        clu = pd.to_numeric(pd.Series([row.get("clustering", np.nan)]), errors="coerce").iloc[0]
        if not (np.isfinite(den) and np.isfinite(clu)):
            continue
        status = str(row.get("v4_reference_status", "single_v4_reference"))
        if status not in SAMPLE_STATUS_MARKERS:
            status = "single_v4_reference"
        sample_key = str(row.get("sample_key", ""))
        value_rows.extend(
            [
                {
                    "sample_key": sample_key,
                    "metric_family": family,
                    "method_family": "denoising",
                    "value": float(den),
                    "v4_reference_status": status,
                },
                {
                    "sample_key": sample_key,
                    "metric_family": family,
                    "method_family": "clustering",
                    "value": float(clu),
                    "v4_reference_status": status,
                },
            ]
        )
        line_rows.append(
            {
                "metric_family": family,
                "sample_key": sample_key,
                "denoising": float(den),
                "clustering": float(clu),
                "delta": float(den) - float(clu),
            }
        )
    value_df = pd.DataFrame(value_rows)
    if value_df.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, "Distance-deflation contrast unavailable", ha="center", va="center", transform=ax.transAxes)
        return

    grouped_data: list[tuple[str, str, np.ndarray]] = []
    for family in order:
        for method_family in ("denoising", "clustering"):
            vals = pd.to_numeric(
                value_df.loc[
                    (value_df["metric_family"] == family) & (value_df["method_family"] == method_family),
                    "value",
                ],
                errors="coerce",
            ).dropna().to_numpy(dtype=float)
            if vals.size:
                grouped_data.append((family, method_family, vals))

    for line in line_rows:
        family = str(line["metric_family"])
        x1 = positions[(family, "denoising")]
        x2 = positions[(family, "clustering")]
        color = "#2F7D32" if float(line["delta"]) > 0.0 else "#9B3A32"
        ax.plot(
            [x1, x2],
            [float(line["denoising"]), float(line["clustering"])],
            color=color,
            linewidth=0.85,
            alpha=0.26,
            zorder=1,
        )

    if plot_type == "boxplot":
        for family, method_family, vals in grouped_data:
            ax.boxplot(
                [vals],
                positions=[positions[(family, method_family)]],
                widths=0.48,
                patch_artist=True,
                boxprops=dict(facecolor=method_colors.get(method_family, "#999999"), alpha=0.36, edgecolor="#333333"),
                medianprops=dict(color="#111111", linewidth=1.25),
                whiskerprops=dict(color="#333333", linewidth=0.9),
                capprops=dict(color="#333333", linewidth=0.9),
            )
    else:
        violin_vals = [vals for _, _, vals in grouped_data if vals.size >= 2]
        violin_pos = [float(positions[(family, method_family)]) for family, method_family, vals in grouped_data if vals.size >= 2]
        violin_groups = [(family, method_family) for family, method_family, vals in grouped_data if vals.size >= 2]
        if violin_vals:
            vp = ax.violinplot(
                violin_vals,
                positions=violin_pos,
                widths=0.62,
                showmeans=False,
                showmedians=True,
                showextrema=False,
                bw_method=violin_bw_method,
            )
            for body, (_, method_family) in zip(vp["bodies"], violin_groups):
                body.set_facecolor(method_colors.get(method_family, "#999999"))
                body.set_edgecolor("#333333")
                body.set_alpha(0.40)
            if "cmedians" in vp:
                vp["cmedians"].set_color("#111111")
                vp["cmedians"].set_linewidth(1.15)

    rng = np.random.default_rng(21141)
    scatter_size = max(1.0, float(dot_size)) if dot_size is not None else point_size * 1.62
    for family in order:
        for method_family in ("denoising", "clustering"):
            group_df = value_df.loc[
                (value_df["metric_family"] == family) & (value_df["method_family"] == method_family)
            ].copy()
            if group_df.empty:
                continue
            vals = pd.to_numeric(group_df["value"], errors="coerce").to_numpy(dtype=float)
            xj = positions[(family, method_family)] + rng.normal(0.0, max(0.0, float(dot_jitter)), size=vals.size)
            status = group_df["v4_reference_status"].astype(str).to_numpy()
            for status_key, marker in SAMPLE_STATUS_MARKERS.items():
                mask = status == status_key
                if not np.any(mask):
                    continue
                ax.scatter(
                    xj[mask],
                    vals[mask],
                    s=scatter_size,
                    marker=marker,
                    alpha=max(0.0, min(1.0, float(dot_opacity))),
                    color=method_colors.get(method_family, "#666666"),
                    edgecolors="#111111",
                    linewidths=0.45,
                    zorder=3,
                )

    finite_vals = pd.to_numeric(value_df["value"], errors="coerce").dropna().to_numpy(dtype=float)
    if finite_vals.size:
        y_min = float(np.min(finite_vals))
        y_max = float(np.max(finite_vals))
        y_range = y_max - y_min
        if y_range <= 0.0:
            y_range = max(0.05, abs(y_max) * 0.25)
        ax.set_ylim(max(0.0, y_min - 0.10 * y_range), y_max + 0.34 * y_range)

    ax.set_xlim(min(x_ticks) - 0.55, max(x_ticks) + 0.55)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_labels, fontsize=fs["tick"])
    ax.set_xlabel(x_axis_label, fontsize=fs["base"])
    ax.set_ylabel(y_axis_label, fontsize=fs["base"])
    ax.tick_params(axis="y", labelsize=fs["tick"])
    ax.grid(True, axis="y", alpha=0.25, linewidth=0.6)
    ax.set_title(title, fontsize=fs["base"])

    if stats_enabled and stats_rows:
        y_min, y_max = ax.get_ylim()
        y_range = float(y_max - y_min) if y_max > y_min else 1.0
        stat_by_family = {str(r.get("metric_family", "")): r for r in stats_rows}
        for family in order:
            row = stat_by_family.get(family)
            if not row:
                continue
            p = float(row.get("p_value_used", row.get("p_adj_bh", row.get("p_raw", np.nan))))
            if not np.isfinite(p):
                continue
            x1 = positions[(family, "denoising")]
            x2 = positions[(family, "clustering")]
            y = y_max - 0.18 * y_range
            cap = 0.018 * y_range
            ax.plot([x1, x1, x2, x2], [y, y + cap, y + cap, y], color="#111111", linewidth=0.95, zorder=4)
            ax.text(
                (x1 + x2) / 2.0,
                y + cap + 0.012 * y_range,
                f"p={p:.3g}\n{_pvalue_stars(p, alpha=stats_alpha)}",
                ha="center",
                va="bottom",
                fontsize=fs["stats"],
                color="#111111",
                zorder=5,
            )

    handles = [
        Patch(facecolor=method_colors["denoising"], edgecolor="#333333", alpha=0.40, label="Denoising"),
        Patch(facecolor=method_colors["clustering"], edgecolor="#333333", alpha=0.40, label="Clustering"),
        Line2D(
            [0],
            [0],
            marker=SAMPLE_STATUS_MARKERS["multiple_v4_reference"],
            linestyle="None",
            markerfacecolor="#666666",
            markeredgecolor="#111111",
            markersize=9.75,
            label=sample_v4_multiple_label,
        ),
        Line2D(
            [0],
            [0],
            marker=SAMPLE_STATUS_MARKERS["single_v4_reference"],
            linestyle="None",
            markerfacecolor="#666666",
            markeredgecolor="#111111",
            markersize=9.75,
            label=sample_v4_single_label,
        ),
    ]
    ax.legend(
        handles=handles,
        title="Encoding",
        loc="upper left",
        frameon=True,
        framealpha=0.92,
        facecolor="white",
        edgecolor="#BBBBBB",
        fontsize=fs["legend"],
        title_fontsize=fs["legend_title"],
    )

    # Keep the formal delta interpretation visible without making the reader decode the y-axis.
    ax.text(
        0.99,
        0.02,
        "paired test: Denoising - Clustering; downward paired lines mean clustering reduces distance",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=fs["stats"],
        color="#333333",
        bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor="#bbbbbb", alpha=0.86),
    )


def _plot_sample_method_heatmap(
    sens: pd.DataFrame,
    out_dir: Path,
    method_order: list[str],
    font_size: float,
) -> None:
    if sens.empty:
        return
    matrix = sens.pivot_table(index="sample_key", columns="method", values="sensitivity", aggfunc="first")
    matrix = matrix.reindex(columns=[m for m in method_order if m in matrix.columns])
    sample_meta = (
        sens[["sample_key", "community_id", "sample_id", "sample_label"]]
        .drop_duplicates()
        .sort_values(["community_id", "sample_id", "sample_key"])
    )
    sample_order = [s for s in sample_meta["sample_key"].astype(str).tolist() if s in matrix.index]
    matrix = matrix.reindex(index=sample_order)
    matrix.to_csv(out_dir / "sample_method_sensitivity_matrix.tsv", sep="\t")
    sens.to_csv(out_dir / "sample_method_sensitivity_long.tsv", sep="\t", index=False)

    sample_summary = matrix.copy()
    sample_summary["mean_sensitivity"] = matrix.mean(axis=1, skipna=True)
    sample_summary["median_sensitivity"] = matrix.median(axis=1, skipna=True)
    sample_summary["sd_sensitivity"] = matrix.std(axis=1, skipna=True)
    sample_summary = sample_summary.reset_index().rename(columns={"index": "sample_key"})
    sample_summary = sample_summary.merge(sample_meta, on="sample_key", how="left")
    sample_summary.to_csv(out_dir / "sample_method_sensitivity_sample_summary.tsv", sep="\t", index=False)
    _method_summary(sens).to_csv(out_dir / "sample_method_sensitivity_method_summary.tsv", sep="\t", index=False)

    if matrix.empty:
        return
    arr = matrix.to_numpy(dtype=float)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return
    vmax = float(np.nanmax(np.abs(finite)))
    if vmax <= 0.0:
        vmax = 1.0
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    n_rows = int(matrix.shape[0])
    n_cols = int(matrix.shape[1])
    fig_w = max(7.0, 1.15 * n_cols + 3.2)
    fig_h = max(5.0, min(18.0, 0.34 * n_rows + 2.2))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(arr, aspect="auto", interpolation="nearest", cmap="RdBu_r", norm=norm)
    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels([display_method(m) for m in matrix.columns], rotation=0, fontsize=max(6.0, font_size - 2.0))
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels([display_sample_key(s) for s in matrix.index], fontsize=max(5.0, font_size - 5.0))
    ax.set_xlabel("Method", fontsize=font_size)
    ax.set_ylabel("Sample", fontsize=font_size)
    ax.set_title("Sample-wise Metric-Family Sensitivity", fontsize=font_size + 1.0)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Agnostic mean - phylogeny-aware mean", fontsize=max(6.0, font_size - 2.0))
    cbar.ax.tick_params(labelsize=max(6.0, font_size - 3.0))
    fig.tight_layout()
    fig.savefig(out_dir / "sample_method_sensitivity_heatmap.svg", format="svg", dpi=300, bbox_inches="tight")
    plt.close(fig)


def _plot_per_sample_method_metric_pcoa(
    long: pd.DataFrame,
    out_dir: Path,
    method_order: list[str],
    metric_order: list[str],
    font_size: float,
    point_size: float,
) -> None:
    """
    For each sample, draw one PCoA where each point is a method x metric
    distance-to-expected value. The point-to-point distance matrix is the
    absolute difference between those scalar distance values.
    """
    if long.empty:
        return

    svg_dir = out_dir / "svg"
    matrix_dir = out_dir / "distance_matrices"
    svg_dir.mkdir(parents=True, exist_ok=True)
    matrix_dir.mkdir(parents=True, exist_ok=True)

    value_rows: list[dict] = []
    coord_rows: list[dict] = []
    eigen_rows: list[dict] = []
    summary_rows: list[dict] = []

    sub_all = long[long["method"].isin(method_order) & long["metric"].isin(metric_order)].copy()
    if sub_all.empty:
        return

    sample_meta = (
        sub_all[["community_id", "sample_id", "sample_key"]]
        .drop_duplicates()
        .sort_values(["community_id", "sample_id", "sample_key"])
    )

    for _, smeta in sample_meta.iterrows():
        sample_key = str(smeta["sample_key"])
        sample_rows = sub_all[sub_all["sample_key"].astype(str) == sample_key].copy()
        if sample_rows.empty:
            continue

        # One scalar distance value per method x metric point.
        pivot = (
            sample_rows.pivot_table(
                index=["method", "metric"],
                values="distance",
                aggfunc="first",
            )
            .reset_index()
        )
        order_pairs = [(m, metric) for m in method_order for metric in metric_order]
        order_df = pd.DataFrame(order_pairs, columns=["method", "metric"])
        points = order_df.merge(pivot, on=["method", "metric"], how="left")
        points = points.dropna(subset=["distance"]).copy()
        if points.shape[0] < 3:
            continue

        points["community_id"] = str(smeta["community_id"])
        points["sample_id"] = str(smeta["sample_id"])
        points["sample_key"] = sample_key
        points["community_label"] = points["community_id"].map(display_community_id)
        points["sample_label"] = display_sample_key(sample_key)
        points["method_label"] = points["method"].map(display_method)
        points["metric_label"] = points["metric"].map(_metric_label)
        points["point_id"] = points["method"].astype(str) + "__" + points["metric"].astype(str)

        vals = pd.to_numeric(points["distance"], errors="coerce").to_numpy(dtype=float)
        dmat = np.abs(vals[:, None] - vals[None, :])
        coords, eigvals = _pcoa(dmat)
        if coords.shape[1] == 0:
            coords = np.zeros((points.shape[0], 2), dtype=float)
            eigvals = np.array([0.0, 0.0], dtype=float)
        elif coords.shape[1] == 1:
            coords = np.column_stack([coords[:, 0], np.zeros(points.shape[0], dtype=float)])
            eigvals = np.array([float(eigvals[0]), 0.0], dtype=float)
        elif coords.shape[1] >= 2:
            coords = coords[:, :2]

        eig_sum = float(np.sum(eigvals)) if eigvals.size else 0.0
        pc1_pct = 100.0 * float(eigvals[0] / eig_sum) if eig_sum > 0.0 and eigvals.size >= 1 else 0.0
        pc2_pct = 100.0 * float(eigvals[1] / eig_sum) if eig_sum > 0.0 and eigvals.size >= 2 else 0.0
        safe = _safe_filename(sample_key)

        matrix = pd.DataFrame(dmat, index=points["point_id"], columns=points["point_id"])
        matrix.to_csv(matrix_dir / f"{safe}_method_metric_distance_difference_matrix.tsv", sep="\t")

        for i in range(max(2, len(eigvals))):
            eigen_rows.append(
                {
                    "community_id": str(smeta["community_id"]),
                    "sample_id": str(smeta["sample_id"]),
                    "sample_key": sample_key,
                    "axis": f"PC{i + 1}",
                    "eigenvalue": float(eigvals[i]) if i < len(eigvals) else 0.0,
                    "variance_explained": float((eigvals[i] / eig_sum) if eig_sum > 0.0 and i < len(eigvals) else 0.0),
                    "n_points": int(points.shape[0]),
                }
            )

        for idx, row in points.reset_index(drop=True).iterrows():
            rec = row.to_dict()
            rec["PC1"] = float(coords[idx, 0])
            rec["PC2"] = float(coords[idx, 1])
            rec["PC1_variance_explained"] = pc1_pct / 100.0
            rec["PC2_variance_explained"] = pc2_pct / 100.0
            coord_rows.append(rec)
            value_rows.append(row.to_dict())

        summary_rows.append(
            {
                "community_id": str(smeta["community_id"]),
                "sample_id": str(smeta["sample_id"]),
                "sample_key": sample_key,
                "sample_label": display_sample_key(sample_key),
                "n_points": int(points.shape[0]),
                "n_methods": int(points["method"].nunique()),
                "n_metrics": int(points["metric"].nunique()),
                "min_distance": float(np.min(vals)),
                "max_distance": float(np.max(vals)),
                "mean_distance": float(np.mean(vals)),
                "sd_distance": float(np.std(vals, ddof=1)) if vals.size > 1 else 0.0,
                "pc1_variance_explained": pc1_pct / 100.0,
                "pc2_variance_explained": pc2_pct / 100.0,
                "distance_matrix": str(matrix_dir / f"{safe}_method_metric_distance_difference_matrix.tsv"),
                "svg": str(svg_dir / f"{safe}_method_metric_pcoa.svg"),
            }
        )

        fig, ax = plt.subplots(figsize=(8.8, 5.8))
        for metric in metric_order:
            marker = METRIC_MARKERS.get(metric, "o")
            for method in method_order:
                mask = (points["metric"].astype(str) == metric) & (points["method"].astype(str) == method)
                if not mask.any():
                    continue
                idxs = np.where(mask.to_numpy())[0]
                ax.scatter(
                    coords[idxs, 0],
                    coords[idxs, 1],
                    s=point_size,
                    marker=marker,
                    color=METHOD_COLORS.get(method, "#666666"),
                    alpha=0.86,
                    edgecolors="#111111",
                    linewidths=0.45,
                    zorder=3,
                )

        # If the scalar-difference PCoA is mathematically one-dimensional,
        # give the zero PC2 axis enough breathing room to show all markers.
        y_vals = coords[:, 1]
        if np.nanmax(np.abs(y_vals)) <= 1e-12:
            ax.set_ylim(-0.08, 0.08)

        ax.axhline(0.0, color="#888888", linewidth=0.6, alpha=0.45)
        ax.axvline(0.0, color="#888888", linewidth=0.6, alpha=0.45)
        ax.grid(True, alpha=0.25, linewidth=0.6)
        ax.set_xlabel(f"PC1 ({pc1_pct:.1f}% var. explained)", fontsize=font_size)
        ax.set_ylabel(f"PC2 ({pc2_pct:.1f}% var. explained)", fontsize=font_size)
        ax.tick_params(axis="both", labelsize=max(6.0, font_size - 2.0))
        ax.set_title(f"{display_sample_key(sample_key)}: method x metric distance PCoA", fontsize=font_size + 1.0)

        method_handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                markerfacecolor=METHOD_COLORS.get(method, "#666666"),
                markeredgecolor="#111111",
                markersize=7,
                label=display_method(method),
            )
            for method in method_order
            if method in set(points["method"].astype(str))
        ]
        metric_handles = [
            Line2D(
                [0],
                [0],
                marker=METRIC_MARKERS.get(metric, "o"),
                color="#333333",
                markerfacecolor="#DDDDDD",
                markeredgecolor="#111111",
                linestyle="none",
                markersize=7,
                label=_metric_label(metric),
            )
            for metric in metric_order
            if metric in set(points["metric"].astype(str))
        ]
        leg1 = ax.legend(
            handles=method_handles,
            title="Method",
            loc="upper left",
            bbox_to_anchor=(1.01, 1.00),
            frameon=False,
            fontsize=max(7.0, font_size - 2.0),
            title_fontsize=max(7.5, font_size - 1.5),
        )
        ax.add_artist(leg1)
        ax.legend(
            handles=metric_handles,
            title="Metric",
            loc="upper left",
            bbox_to_anchor=(1.01, 0.50),
            frameon=False,
            fontsize=max(7.0, font_size - 2.0),
            title_fontsize=max(7.5, font_size - 1.5),
        )
        fig.tight_layout()
        fig.savefig(svg_dir / f"{safe}_method_metric_pcoa.svg", format="svg", dpi=300, bbox_inches="tight")
        plt.close(fig)

    if value_rows:
        pd.DataFrame(value_rows).to_csv(out_dir / "per_sample_method_metric_distance_values.tsv", sep="\t", index=False)
    if coord_rows:
        pd.DataFrame(coord_rows).to_csv(out_dir / "per_sample_method_metric_pcoa_coords.tsv", sep="\t", index=False)
    if eigen_rows:
        pd.DataFrame(eigen_rows).to_csv(out_dir / "per_sample_method_metric_pcoa_eigenvalues.tsv", sep="\t", index=False)
    if summary_rows:
        pd.DataFrame(summary_rows).to_csv(out_dir / "per_sample_method_metric_pcoa_sample_summary.tsv", sep="\t", index=False)


def _endpoint_group(
    method: str,
    metric: str,
    agnostic_metrics: list[str],
    aware_metrics: list[str],
    denoising_methods: tuple[str, ...] | list[str] = FIG21_DENOISING_METHODS,
    clustering_methods: tuple[str, ...] | list[str] = FIG21_CLUSTERING_METHODS,
) -> str | None:
    method = str(method)
    metric = str(metric)
    if method in denoising_methods:
        method_family = "denoising"
    elif method in clustering_methods:
        method_family = "clustering"
    else:
        return None

    if metric in aware_metrics:
        metric_family = "phylogeny_aware"
    elif metric in agnostic_metrics:
        metric_family = "phylogeny_agnostic"
    else:
        return None
    return f"{method_family}_{metric_family}"


def _pair_group_id(group_a: str, group_b: str) -> str:
    if ENDPOINT_GROUP_RANK[group_a] <= ENDPOINT_GROUP_RANK[group_b]:
        return f"{group_a}__vs__{group_b}"
    return f"{group_b}__vs__{group_a}"


def _draw_panel_f_stats_brackets(
    ax,
    stats_rows: list[dict],
    positions: dict[tuple[str, str], float],
    alpha: float,
    font_size: float,
    show_non_significant: bool,
) -> None:
    """Draw Panel F method-pair statistics above each metric-family block."""
    intervals: list[tuple[str, float, float, float]] = []
    for row in stats_rows:
        family = str(row.get("metric_family", ""))
        method_a = str(row.get("method_a", ""))
        method_b = str(row.get("method_b", ""))
        p_value = float(row.get("p_value_used", row.get("p_adj_bh", row.get("p_raw", np.nan))))
        key_a = (family, method_a)
        key_b = (family, method_b)
        if key_a not in positions or key_b not in positions or not np.isfinite(p_value):
            continue
        if not show_non_significant and p_value >= alpha:
            continue
        x1, x2 = sorted((positions[key_a], positions[key_b]))
        intervals.append((family, x1, x2, p_value))
    if not intervals:
        return

    levels_by_family: dict[str, list[tuple[float, float]]] = {}
    placed: list[tuple[float, float, int, float]] = []
    for family, x1, x2, p_value in sorted(intervals, key=lambda item: (item[0], item[2] - item[1], item[3])):
        levels = levels_by_family.setdefault(family, [])
        level = 0
        while level < len(levels) and not all((x2 < start) or (x1 > end) for start, end in levels[level]):
            level += 1
        if level == len(levels):
            levels.append([])
        levels[level].append((x1, x2))
        placed.append((x1, x2, level, p_value))

    y_min, y_max = ax.get_ylim()
    y_range = max(1e-6, float(y_max - y_min))
    max_levels = max(len(levels) for levels in levels_by_family.values())
    cap_height = 0.025 * y_range
    step = 0.105 * y_range
    text_offset = 0.018 * y_range
    base = y_max + 0.055 * y_range
    ax.set_ylim(y_min, y_max + (0.13 + max_levels * 0.115) * y_range)
    for x1, x2, level, p_value in placed:
        y = base + level * step
        ax.plot([x1, x1, x2, x2], [y, y + cap_height, y + cap_height, y], color="#111111", linewidth=0.9, zorder=5)
        ax.text(
            (x1 + x2) / 2.0,
            y + cap_height + text_offset,
            f"p={p_value:.3g} {_pvalue_stars(p_value, alpha=alpha)}",
            ha="center",
            va="bottom",
            fontsize=max(6.0, font_size - 4.0),
            color="#111111",
            zorder=6,
        )


def _draw_method_metric_family_distance_panel(
    ax,
    points: pd.DataFrame,
    methods: list[str],
    title: str,
    plot_type: str,
    point_size: float,
    font_size: float,
    violin_bw_method: float | str | None,
    x_axis_label: str,
    y_axis_label: str,
    x_tick_rotation: float = 0.0,
    dot_size: float | None = None,
    dot_opacity: float = 0.75,
    dot_jitter: float = 0.06,
    summary_statistic: str = "median",
    stats_rows: list[dict] | None = None,
    stats_alpha: float = 0.05,
    show_stats_annotations: bool = True,
    show_non_significant_annotations: bool = True,
) -> None:
    """Panel F: selected methods separately, grouped only by metric family."""
    fs = _fig21_text_scale(font_size)
    if points.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, "Method-level distance comparison unavailable", ha="center", va="center", transform=ax.transAxes)
        return
    families = [f for f in DISTANCE_DEFLATION_FAMILY_ORDER if f in set(points["metric_family"].astype(str))]
    methods = [m for m in methods if m in set(points["method"].astype(str))]
    summary_statistic = _normalize_panel_f_summary_statistic(summary_statistic)
    show_median = summary_statistic in {"median", "both"}
    show_mean = summary_statistic in {"mean", "both"}
    if not families or not methods:
        ax.axis("off")
        return

    positions: dict[tuple[str, str], float] = {}
    x_ticks: list[float] = []
    x_labels: list[str] = []
    step = len(methods) + 1.25
    for family_idx, family in enumerate(families):
        start = 1.0 + family_idx * step
        ax.axvspan(start - 0.45, start + len(methods) - 0.55, color=("#FFF1E6" if family == "phylogeny_agnostic" else "#EAF4FB"), alpha=0.55, zorder=0)
        for method_idx, method in enumerate(methods):
            x = start + method_idx
            positions[(family, method)] = x
            x_ticks.append(x)
            x_labels.append(f"{DISTANCE_DEFLATION_FAMILY_LABELS.get(family, family)}\n{display_method(method)}")

    rng = np.random.default_rng(21161)
    scatter_size = max(1.0, float(dot_size)) if dot_size is not None else point_size * 1.62
    for family in families:
        family_df = points.loc[points["metric_family"].astype(str) == family].copy()
        for sample_key, sdf in family_df.groupby("sample_key", sort=False):
            sdf = sdf.set_index("method").reindex(methods).dropna(subset=["distance"]).reset_index()
            if sdf.shape[0] > 1:
                ax.plot(
                    [positions[(family, str(m))] for m in sdf["method"]],
                    pd.to_numeric(sdf["distance"], errors="coerce"),
                    color="#666666", alpha=0.18, linewidth=0.75, zorder=1,
                )
        for method in methods:
            vals = pd.to_numeric(
                family_df.loc[family_df["method"].astype(str) == method, "distance"], errors="coerce"
            ).dropna().to_numpy(dtype=float)
            if not vals.size:
                continue
            pos = positions[(family, method)]
            if plot_type == "boxplot":
                ax.boxplot([vals], positions=[pos], widths=0.54, patch_artist=True,
                           boxprops=dict(facecolor=METHOD_COLORS.get(method, "#999999"), alpha=0.38, edgecolor="#333333"),
                           showmeans=show_mean, meanline=True,
                           meanprops=dict(color="#333333", linewidth=1.2, linestyle="--"),
                           medianprops=dict(color="#111111", linewidth=1.2 if show_median else 0.0),
                           whiskerprops=dict(color="#333333", linewidth=0.9), capprops=dict(color="#333333", linewidth=0.9))
            elif vals.size >= 2:
                vp = ax.violinplot([vals], positions=[pos], widths=0.68, showmeans=show_mean, showmedians=show_median, showextrema=False, bw_method=violin_bw_method)
                vp["bodies"][0].set_facecolor(METHOD_COLORS.get(method, "#999999"))
                vp["bodies"][0].set_edgecolor("#333333")
                vp["bodies"][0].set_alpha(0.46)
                if "cmedians" in vp:
                    vp["cmedians"].set_color("#111111")
                    vp["cmedians"].set_linewidth(1.1)
                if "cmeans" in vp:
                    vp["cmeans"].set_color("#333333")
                    vp["cmeans"].set_linewidth(1.1)
                    vp["cmeans"].set_linestyle("--")
            jitter = pos + rng.normal(0.0, max(0.0, float(dot_jitter)), size=vals.size)
            ax.scatter(jitter, vals, s=scatter_size, color=METHOD_COLORS.get(method, "#666666"), edgecolors="#111111", linewidths=0.45, alpha=max(0.0, min(1.0, float(dot_opacity))), zorder=3)

    ax.set_xticks(x_ticks)
    rotation = float(x_tick_rotation)
    ax.set_xticklabels(
        x_labels,
        fontsize=fs["tick"],
        rotation=rotation,
        ha="center" if rotation % 180.0 == 0.0 else "right",
        rotation_mode="anchor",
    )
    ax.set_xlabel(x_axis_label, fontsize=fs["base"])
    ax.set_ylabel(y_axis_label, fontsize=fs["base"])
    ax.tick_params(axis="y", labelsize=fs["tick"])
    ax.grid(True, axis="y", alpha=0.25, linewidth=0.6)
    ax.set_title(title, fontsize=fs["base"])
    if show_stats_annotations:
        _draw_panel_f_stats_brackets(
            ax=ax,
            stats_rows=list(stats_rows or []),
            positions=positions,
            alpha=stats_alpha,
            font_size=fs["base"],
            show_non_significant=show_non_significant_annotations,
        )


def _panel_a_endpoint_group(method: str, metric: str, agnostic_metrics: list[str], aware_metrics: list[str]) -> str | None:
    method = str(method)
    metric = str(metric)
    if method in FIG21_DENOISING_METHODS:
        method_family = "denoising"
    elif method == "OTU":
        method_family = "uparse"
    elif method in {"TIC", "TAC"}:
        method_family = "tactic"
    else:
        return None

    if metric in aware_metrics:
        metric_family = "phylogeny_aware"
    elif metric in agnostic_metrics:
        metric_family = "phylogeny_agnostic"
    else:
        return None
    return f"{method_family}_{metric_family}"


def _panel_a_pair_group_id(group_a: str, group_b: str) -> str:
    ranks = {g: i for i, g in enumerate(PANEL_A_ENDPOINT_GROUP_ORDER)}
    if ranks[group_a] <= ranks[group_b]:
        return f"{group_a}__vs__{group_b}"
    return f"{group_b}__vs__{group_a}"


def _pair_group_label(pair_group_id: str, endpoint_group_labels: dict[str, str] | None = None) -> str:
    parts = str(pair_group_id).split("__vs__")
    if len(parts) != 2:
        return str(pair_group_id)
    a, b = parts
    # Keep Panel B wording aligned with Panel A style: "Method-group | phylo-*".
    labels = endpoint_group_labels or ENDPOINT_GROUP_LABELS
    la = labels.get(a, a.replace("_", " "))
    lb = labels.get(b, b.replace("_", " "))
    if a == b:
        return f"{la}\nwithin"
    return f"{la}\nvs\n{lb}"


def _pair_group_full_label(pair_group_id: str, endpoint_group_labels: dict[str, str] | None = None) -> str:
    parts = str(pair_group_id).split("__vs__")
    if len(parts) != 2:
        return str(pair_group_id)
    a, b = parts
    labels = endpoint_group_labels or ENDPOINT_GROUP_LABELS
    la = labels.get(a, a.replace("_", " "))
    lb = labels.get(b, b.replace("_", " "))
    if a == b:
        return f"{la} within-group"
    return f"{la} vs {lb}"


def _normalize_violin_bw_method(value):
    if value in (None, ""):
        return None
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in {"none", "null", "default", "auto"}:
            return None
        if lowered in {"scott", "silverman"}:
            return lowered
        try:
            return float(value)
        except ValueError:
            return value
    return value


def _resolve_fig21_panel_font_size(value, fallback: float) -> float:
    if value is None:
        return fallback
    return max(6.0, float(value))


def _fig21_text_scale(font_size: float) -> dict[str, float]:
    """Figure 21 typography scale aligned with Figure 09/11 conventions."""
    base = max(6.0, float(font_size))
    return {
        "base": base,                  # titles and axis labels
        "tick": max(6.0, base - 2.0),  # tick labels
        "legend": max(6.0, base - 2.0),
        "legend_title": max(6.0, base - 2.0),
        "cbar_label": base,
        "cbar_tick": max(6.0, base - 2.0),
        "stats": max(6.0, base - 3.0),
        "pair_label": max(6.0, base - 3.5),
        "pair_star": max(6.0, base - 2.8),
    }


def _resolve_fig21_panel_title(value, fallback: str) -> str:
    if value is None:
        return fallback
    text = str(value).strip()
    return text or fallback


def _normalize_summary_stat(value: str | None) -> str:
    if str(value).strip().lower() == "mean":
        return "mean"
    return "median"


def _summary_stat_colorbar_label(summary_stat: str, summarized: bool) -> str:
    prefix = "Mean" if _normalize_summary_stat(summary_stat) == "mean" else "Median"
    if summarized:
        return f"{prefix} summarized absolute difference"
    return f"{prefix} absolute difference"


def _pair_group_to_panel_a_b_cell(pair_group_id: str) -> tuple[float, float] | None:
    parts = str(pair_group_id).split("__vs__")
    if len(parts) != 2:
        return None
    a, b = parts
    if a not in ENDPOINT_GROUP_RANK or b not in ENDPOINT_GROUP_RANK:
        return None
    ia = int(ENDPOINT_GROUP_RANK[a])
    ib = int(ENDPOINT_GROUP_RANK[b])
    if ia > ib:
        ia, ib = ib, ia
    # imshow cell centers are integer coordinates (x=col, y=row).
    return float(ib), float(ia)


def _draw_panel_a_b_planned_comparison_annotations(
    ax,
    planned_sig_rows: list[dict],
    alpha: float,
    font_size: float,
) -> int:
    if not planned_sig_rows:
        return 0

    fs = _fig21_text_scale(font_size)
    drawn = 0
    n = max(1, len(planned_sig_rows))
    for i, row in enumerate(planned_sig_rows):
        pg1 = str(row.get("method_a", ""))
        pg2 = str(row.get("method_b", ""))
        p = float(row.get("p_value_used", row.get("p_adj_bh", row.get("p_raw", np.nan))))
        if not np.isfinite(p):
            continue
        c1 = _pair_group_to_panel_a_b_cell(pg1)
        c2 = _pair_group_to_panel_a_b_cell(pg2)
        if c1 is None or c2 is None:
            continue
        x1, y1 = c1
        x2, y2 = c2
        if np.isclose(x1, x2) and np.isclose(y1, y2):
            continue

        ax.plot([x1, x2], [y1, y2], color="#111111", linewidth=1.35, alpha=0.95, zorder=7)
        ax.scatter([x1, x2], [y1, y2], s=22.0, color="#111111", zorder=8)

        # Offset labels slightly perpendicular to reduce overlap when multiple lines exist.
        dx = x2 - x1
        dy = y2 - y1
        norm = float(np.hypot(dx, dy))
        if norm > 0.0:
            px, py = (-dy / norm), (dx / norm)
        else:
            px, py = (0.0, 0.0)
        off = (i - (n - 1) / 2.0) * 0.16
        mx = 0.5 * (x1 + x2) + off * px
        my = 0.5 * (y1 + y2) + off * py
        label = f"p={p:.3g}\n{_pvalue_stars(p, alpha=alpha)}"
        ax.text(
            mx,
            my,
            label,
            ha="center",
            va="center",
            fontsize=max(5.5, fs["pair_label"] - 1.0),
            color="#111111",
            zorder=9,
            bbox=dict(boxstyle="round,pad=0.14", facecolor="white", edgecolor="#888888", alpha=0.92),
        )
        drawn += 1
    return drawn


def _panel_b_within_contrast_stat_rows(stats: dict | None, visible_pair_groups: set[str]) -> list[dict]:
    if not isinstance(stats, dict):
        return []
    within = {g for g in WITHIN_PAIR_GROUP_ORDER if g in visible_pair_groups}
    if len(within) < 2:
        return []
    rows = []
    for row in stats.get("pairwise_sig", []):
        method_a = str(row.get("method_a", ""))
        method_b = str(row.get("method_b", ""))
        if method_a in within and method_b in within:
            rows.append(row)
    return rows


def _expand_endpoint_matrix_with_gaps(
    matrix: pd.DataFrame,
    endpoint_meta: pd.DataFrame,
    gap_cells: int,
) -> tuple[np.ndarray, np.ndarray, list[str], list[str], list[str], list[int]]:
    endpoint_ids = endpoint_meta["endpoint_id"].astype(str).tolist()
    endpoint_labels = endpoint_meta["endpoint_label"].astype(str).tolist()
    methods = endpoint_meta["method"].astype(str).tolist()
    base = matrix.reindex(index=endpoint_ids, columns=endpoint_ids).to_numpy(dtype=float)
    if gap_cells <= 0 or len(endpoint_ids) <= 1:
        return base, np.arange(len(endpoint_ids), dtype=int), endpoint_ids, endpoint_labels, methods, []

    boundaries = [i for i in range(1, len(methods)) if methods[i] != methods[i - 1]]
    if not boundaries:
        return base, np.arange(len(endpoint_ids), dtype=int), endpoint_ids, endpoint_labels, methods, []

    gap_cells = int(gap_cells)
    offset = 0
    positions: list[int] = []
    for idx in range(len(endpoint_ids)):
        if idx in boundaries:
            offset += gap_cells
        positions.append(idx + offset)

    expanded_size = len(endpoint_ids) + gap_cells * len(boundaries)
    expanded = np.full((expanded_size, expanded_size), np.nan, dtype=float)
    expanded[np.ix_(positions, positions)] = base
    return expanded, np.asarray(positions, dtype=int), endpoint_ids, endpoint_labels, methods, boundaries


def _pair_group_distance_contrast_tables(
    long: pd.DataFrame,
    method_order: list[str],
    metric_order: list[str],
    agnostic_metrics: list[str],
    aware_metrics: list[str],
    summary_stat: str,
    denoising_methods: tuple[str, ...] | list[str] = FIG21_DENOISING_METHODS,
    clustering_methods: tuple[str, ...] | list[str] = FIG21_CLUSTERING_METHODS,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cols = ["community_id", "sample_id", "sample_key", "method", "metric", "distance"]
    if any(c not in long.columns for c in cols):
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    method_set = set(denoising_methods) | set(clustering_methods)
    sub = long[
        long["method"].astype(str).isin([m for m in method_order if m in method_set])
        & long["metric"].astype(str).isin(metric_order)
    ].copy()
    if sub.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    ordered_endpoints = [
        (method, metric)
        for method in method_order
        if method in method_set
        for metric in metric_order
        if _endpoint_group(
            method, metric, agnostic_metrics, aware_metrics, denoising_methods, clustering_methods
        ) is not None
    ]
    order_df = pd.DataFrame(ordered_endpoints, columns=["method", "metric"])
    if order_df.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    endpoint_rows: list[dict] = []
    pair_rows: list[dict] = []
    sample_meta = (
        sub[["community_id", "sample_id", "sample_key"]]
        .drop_duplicates()
        .sort_values(["community_id", "sample_id", "sample_key"])
    )
    for _, smeta in sample_meta.iterrows():
        sample_key = str(smeta["sample_key"])
        sample_rows = sub[sub["sample_key"].astype(str) == sample_key].copy()
        if sample_rows.empty:
            continue

        pivot = (
            sample_rows.pivot_table(index=["method", "metric"], values="distance", aggfunc="first")
            .reset_index()
        )
        endpoints = order_df.merge(pivot, on=["method", "metric"], how="left").dropna(subset=["distance"]).copy()
        if endpoints.shape[0] < 2:
            continue
        endpoints["community_id"] = str(smeta["community_id"])
        endpoints["sample_id"] = str(smeta["sample_id"])
        endpoints["sample_key"] = sample_key
        endpoints["community_label"] = endpoints["community_id"].map(display_community_id)
        endpoints["sample_label"] = display_sample_key(sample_key)
        endpoints["method_label"] = endpoints["method"].map(display_method)
        endpoints["metric_label"] = endpoints["metric"].map(_metric_label)
        endpoints["endpoint_group_id"] = endpoints.apply(
            lambda r: _endpoint_group(
                str(r["method"]),
                str(r["metric"]),
                agnostic_metrics,
                aware_metrics,
                denoising_methods,
                clustering_methods,
            ),
            axis=1,
        )
        endpoints = endpoints.dropna(subset=["endpoint_group_id"]).reset_index(drop=True)
        if endpoints.shape[0] < 2:
            continue
        endpoints["endpoint_group_label"] = endpoints["endpoint_group_id"].map(ENDPOINT_GROUP_LABELS)
        endpoints["endpoint_id"] = endpoints["method"].astype(str) + "__" + endpoints["metric"].astype(str)
        endpoint_rows.extend(endpoints.to_dict("records"))

        for i, j in combinations(range(endpoints.shape[0]), 2):
            a = endpoints.iloc[i]
            b = endpoints.iloc[j]
            group_a = str(a["endpoint_group_id"])
            group_b = str(b["endpoint_group_id"])
            pair_group = _pair_group_id(group_a, group_b)
            signed = float(a["distance"]) - float(b["distance"])
            pair_rows.append(
                {
                    "community_id": str(smeta["community_id"]),
                    "sample_id": str(smeta["sample_id"]),
                    "sample_key": sample_key,
                    "community_label": display_community_id(str(smeta["community_id"])),
                    "sample_label": display_sample_key(sample_key),
                    "endpoint_a": str(a["endpoint_id"]),
                    "endpoint_b": str(b["endpoint_id"]),
                    "method_a": str(a["method"]),
                    "method_b": str(b["method"]),
                    "metric_a": str(a["metric"]),
                    "metric_b": str(b["metric"]),
                    "endpoint_group_a": group_a,
                    "endpoint_group_b": group_b,
                    "endpoint_group_a_label": ENDPOINT_GROUP_LABELS.get(group_a, group_a),
                    "endpoint_group_b_label": ENDPOINT_GROUP_LABELS.get(group_b, group_b),
                    "pair_group_id": pair_group,
                    "pair_group_label": _pair_group_full_label(pair_group),
                    "distance_a": float(a["distance"]),
                    "distance_b": float(b["distance"]),
                    "signed_difference_a_minus_b": signed,
                    "abs_difference": abs(signed),
                }
            )

    endpoints_df = pd.DataFrame(endpoint_rows)
    pairwise_df = pd.DataFrame(pair_rows)
    if pairwise_df.empty:
        return endpoints_df, pairwise_df, pd.DataFrame()

    summary = (
        pairwise_df.groupby(
            ["community_id", "sample_id", "sample_key", "community_label", "sample_label", "pair_group_id", "pair_group_label"],
            sort=False,
        )
        .agg(
            n_pairs=("abs_difference", "size"),
            mean_abs_difference=("abs_difference", "mean"),
            median_abs_difference=("abs_difference", "median"),
            sd_abs_difference=("abs_difference", "std"),
            min_abs_difference=("abs_difference", "min"),
            max_abs_difference=("abs_difference", "max"),
            mean_signed_difference=("signed_difference_a_minus_b", "mean"),
            median_signed_difference=("signed_difference_a_minus_b", "median"),
        )
        .reset_index()
    )
    summary["sd_abs_difference"] = pd.to_numeric(summary["sd_abs_difference"], errors="coerce").fillna(0.0)
    if summary_stat == "mean":
        summary["summary_abs_difference"] = summary["mean_abs_difference"]
    else:
        summary_stat = "median"
        summary["summary_abs_difference"] = summary["median_abs_difference"]
    summary["summary_stat"] = summary_stat
    return endpoints_df, pairwise_df, summary


def _pair_group_summary(sample_summary: pd.DataFrame) -> pd.DataFrame:
    if sample_summary.empty:
        return pd.DataFrame()
    out = (
        sample_summary.groupby(["pair_group_id", "pair_group_label"], sort=False)
        .agg(
            n_samples=("sample_key", "nunique"),
            mean_summary_abs_difference=("summary_abs_difference", "mean"),
            median_summary_abs_difference=("summary_abs_difference", "median"),
            sd_summary_abs_difference=("summary_abs_difference", "std"),
            min_summary_abs_difference=("summary_abs_difference", "min"),
            max_summary_abs_difference=("summary_abs_difference", "max"),
            mean_n_pairs_per_sample=("n_pairs", "mean"),
        )
        .reset_index()
    )
    out["sd_summary_abs_difference"] = pd.to_numeric(out["sd_summary_abs_difference"], errors="coerce").fillna(0.0)
    out["pair_group_order"] = out["pair_group_id"].map({g: i + 1 for i, g in enumerate(PAIR_GROUP_ORDER)})
    return out.sort_values("pair_group_order")


def _pair_group_heatmap_matrix(sample_summary: pd.DataFrame, summary_stat: str = "median") -> pd.DataFrame:
    matrix = pd.DataFrame(np.nan, index=ENDPOINT_GROUP_ORDER, columns=ENDPOINT_GROUP_ORDER, dtype=float)
    if sample_summary.empty:
        return matrix
    stat = _normalize_summary_stat(summary_stat)
    agg = sample_summary.groupby("pair_group_id", sort=False)["summary_abs_difference"]
    values = agg.mean() if stat == "mean" else agg.median()
    for group_a in ENDPOINT_GROUP_ORDER:
        for group_b in ENDPOINT_GROUP_ORDER:
            pair_group = _pair_group_id(group_a, group_b)
            if pair_group in values.index:
                matrix.loc[group_a, group_b] = float(values.loc[pair_group])
    matrix.index = [ENDPOINT_GROUP_LABELS.get(x, x) for x in matrix.index]
    matrix.columns = [ENDPOINT_GROUP_LABELS.get(x, x) for x in matrix.columns]
    return matrix


def _panel_a_heatmap_matrix_from_endpoints(
    endpoints_df: pd.DataFrame,
    agnostic_metrics: list[str],
    aware_metrics: list[str],
    summary_stat: str = "median",
) -> pd.DataFrame:
    matrix = pd.DataFrame(np.nan, index=PANEL_A_ENDPOINT_GROUP_ORDER, columns=PANEL_A_ENDPOINT_GROUP_ORDER, dtype=float)
    if endpoints_df.empty:
        matrix.index = [PANEL_A_ENDPOINT_GROUP_LABELS.get(x, x) for x in matrix.index]
        matrix.columns = [PANEL_A_ENDPOINT_GROUP_LABELS.get(x, x) for x in matrix.columns]
        return matrix

    sub = endpoints_df.copy()
    sub["panel_a_group_id"] = sub.apply(
        lambda r: _panel_a_endpoint_group(
            str(r.get("method", "")),
            str(r.get("metric", "")),
            agnostic_metrics,
            aware_metrics,
        ),
        axis=1,
    )
    sub = sub.dropna(subset=["panel_a_group_id"]).copy()
    if sub.empty:
        matrix.index = [PANEL_A_ENDPOINT_GROUP_LABELS.get(x, x) for x in matrix.index]
        matrix.columns = [PANEL_A_ENDPOINT_GROUP_LABELS.get(x, x) for x in matrix.columns]
        return matrix

    # Panel A must reflect the method membership used to construct the
    # pair-group data. For example, when UPARSE is excluded from Panel B's
    # method groups, do not retain empty UPARSE rows/columns in Panel A_a.
    group_order = [
        group for group in PANEL_A_ENDPOINT_GROUP_ORDER
        if group in set(sub["panel_a_group_id"].astype(str))
    ]
    matrix = pd.DataFrame(np.nan, index=group_order, columns=group_order, dtype=float)

    pair_rows: list[dict] = []
    for _, sample_df in sub.groupby("sample_key", sort=False):
        sample_df = sample_df.reset_index(drop=True)
        if sample_df.shape[0] < 2:
            continue
        for i, j in combinations(range(sample_df.shape[0]), 2):
            a = sample_df.iloc[i]
            b = sample_df.iloc[j]
            group_a = str(a["panel_a_group_id"])
            group_b = str(b["panel_a_group_id"])
            pair_group = _panel_a_pair_group_id(group_a, group_b)
            pair_rows.append(
                {
                    "pair_group_id": pair_group,
                    "abs_difference": abs(float(a["distance"]) - float(b["distance"])),
                }
            )

    pairwise = pd.DataFrame(pair_rows)
    if pairwise.empty:
        matrix.index = [PANEL_A_ENDPOINT_GROUP_LABELS.get(x, x) for x in matrix.index]
        matrix.columns = [PANEL_A_ENDPOINT_GROUP_LABELS.get(x, x) for x in matrix.columns]
        return matrix

    stat = _normalize_summary_stat(summary_stat)
    agg = pairwise.groupby("pair_group_id", sort=False)["abs_difference"]
    values = agg.mean() if stat == "mean" else agg.median()
    for group_a in group_order:
        for group_b in group_order:
            pair_group = _panel_a_pair_group_id(group_a, group_b)
            if pair_group in values.index:
                matrix.loc[group_a, group_b] = float(values.loc[pair_group])

    matrix.index = [PANEL_A_ENDPOINT_GROUP_LABELS.get(x, x) for x in matrix.index]
    matrix.columns = [PANEL_A_ENDPOINT_GROUP_LABELS.get(x, x) for x in matrix.columns]
    return matrix


def _endpoint_distance_heatmap_matrix(
    endpoints_df: pd.DataFrame,
    method_order: list[str],
    metric_order: list[str],
    agnostic_metrics: list[str],
    aware_metrics: list[str],
    summary_stat: str = "median",
    order_mode: str = "metric_family_first",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if endpoints_df.empty:
        return pd.DataFrame(), pd.DataFrame()

    method_set = set(FIG21_DENOISING_METHODS) | set(FIG21_CLUSTERING_METHODS)
    present = set(endpoints_df["endpoint_id"].astype(str).tolist())
    ordered_rows: list[dict] = []
    method_ordering = [m for m in method_order if m in method_set]
    if order_mode == "method_family_first":
        denoising_methods = [m for m in FIG21_DENOISING_METHODS if m in method_ordering]
        clustering_methods = [m for m in FIG21_CLUSTERING_METHODS if m in method_ordering]
        extras = [m for m in method_ordering if m not in set(denoising_methods) | set(clustering_methods)]
        ordered_methods = denoising_methods + clustering_methods + extras
        for method in ordered_methods:
            for family_metrics in (agnostic_metrics, aware_metrics):
                for metric in family_metrics:
                    endpoint_id = f"{method}__{metric}"
                    if endpoint_id not in present:
                        continue
                    ordered_rows.append(
                        {
                            "endpoint_id": endpoint_id,
                            "method": method,
                            "metric": metric,
                            "method_label": display_method(method),
                            "metric_label": _metric_label(metric),
                            "endpoint_group_id": _endpoint_group(method, metric, agnostic_metrics, aware_metrics),
                            "endpoint_group_label": "",
                            "endpoint_label": f"{_metric_label(metric)} | {display_method(method)}",
                            "endpoint_label_multiline": f"{_metric_label(metric)} | {display_method(method)}",
                        }
                    )
    else:
        for family_metrics in (agnostic_metrics, aware_metrics):
            for method in method_ordering:
                for metric in family_metrics:
                    endpoint_id = f"{method}__{metric}"
                    if endpoint_id not in present:
                        continue
                    ordered_rows.append(
                        {
                            "endpoint_id": endpoint_id,
                            "method": method,
                            "metric": metric,
                            "method_label": display_method(method),
                            "metric_label": _metric_label(metric),
                            "endpoint_group_id": _endpoint_group(method, metric, agnostic_metrics, aware_metrics),
                            "endpoint_group_label": "",
                            "endpoint_label": f"{_metric_label(metric)} | {display_method(method)}",
                            "endpoint_label_multiline": f"{_metric_label(metric)} | {display_method(method)}",
                        }
                    )
    endpoint_meta = pd.DataFrame(ordered_rows)
    if endpoint_meta.empty:
        return pd.DataFrame(), pd.DataFrame()
    endpoint_meta["endpoint_group_id"] = endpoint_meta["endpoint_group_id"].fillna("")
    endpoint_meta["endpoint_group_label"] = endpoint_meta["endpoint_group_id"].map(ENDPOINT_GROUP_LABELS)

    ordered_ids = endpoint_meta["endpoint_id"].astype(str).tolist()
    wide = endpoints_df.pivot_table(
        index="sample_key",
        columns="endpoint_id",
        values="distance",
        aggfunc="first",
    )
    wide = wide.reindex(columns=ordered_ids)

    matrix = pd.DataFrame(np.nan, index=ordered_ids, columns=ordered_ids, dtype=float)
    stat = _normalize_summary_stat(summary_stat)
    for i, endpoint_a in enumerate(ordered_ids):
        xa = pd.to_numeric(wide[endpoint_a], errors="coerce")
        for j in range(i, len(ordered_ids)):
            endpoint_b = ordered_ids[j]
            xb = pd.to_numeric(wide[endpoint_b], errors="coerce")
            diff = (xa - xb).abs()
            vals = diff[np.isfinite(diff.to_numpy(dtype=float))]
            if len(vals) == 0:
                continue
            vals_arr = vals.to_numpy(dtype=float)
            value = float(np.mean(vals_arr)) if stat == "mean" else float(np.median(vals_arr))
            matrix.loc[endpoint_a, endpoint_b] = value
            matrix.loc[endpoint_b, endpoint_a] = value
    return matrix, endpoint_meta


def _draw_pair_group_heatmap_panel(
    ax,
    sample_summary: pd.DataFrame,
    font_size: float,
    title: str,
    matrix_override: pd.DataFrame | None = None,
    short_tick_labels: list[str] | None = None,
    xtick_rotation: float = 35.0,
    ytick_rotation: float = 0.0,
    summary_stat: str = "median",
    vmin: float = 0.0,
    vmax_override: float | None = None,
    x_axis_label: str | None = None,
    y_axis_label: str | None = None,
) -> pd.DataFrame:
    fs = _fig21_text_scale(font_size)
    matrix = (
        matrix_override.copy()
        if isinstance(matrix_override, pd.DataFrame)
        else _pair_group_heatmap_matrix(sample_summary, summary_stat=summary_stat)
    )
    arr = matrix.to_numpy(dtype=float)
    finite = arr[np.isfinite(arr)]
    vmax = float(vmax_override) if vmax_override is not None else (float(np.max(finite)) if finite.size else 1.0)
    if vmax <= 0.0:
        vmax = 1.0
    im = ax.imshow(arr, cmap="YlGnBu", vmin=float(vmin), vmax=vmax, aspect="equal")
    ax.set_xticks(np.arange(len(matrix.columns)))
    ax.set_yticks(np.arange(len(matrix.index)))
    if short_tick_labels is None or len(short_tick_labels) != len(matrix.columns):
        short_tick_labels = [str(x) for x in matrix.columns]
    ax.set_xticklabels(
        short_tick_labels,
        rotation=xtick_rotation,
        ha="right",
        va="top",
        rotation_mode="anchor",
        fontsize=fs["tick"],
    )
    ax.set_yticklabels(
        short_tick_labels,
        rotation=ytick_rotation,
        rotation_mode="anchor",
        fontsize=fs["tick"],
        ha="right",
        va="center_baseline",
    )
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            val = arr[i, j]
            if np.isfinite(val):
                color = "white" if val > 0.62 * vmax else "#111111"
                ax.text(j, i, f"{val:.3g}", ha="center", va="center", fontsize=fs["stats"], color=color)
    ax.set_title(title, fontsize=fs["base"])
    if x_axis_label:
        ax.set_xlabel(str(x_axis_label), fontsize=fs["base"])
    if y_axis_label:
        ax.set_ylabel(str(y_axis_label), fontsize=fs["base"])
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    return matrix


def _draw_pair_group_distribution_panel(
    ax,
    sample_summary: pd.DataFrame,
    font_size: float,
    point_size: float,
    plot_type: str,
    stats: dict | None,
    planned_sig_rows: list[dict],
    stats_enabled: bool,
    stats_show_pairwise: bool,
    show_within_contrast_stats: bool,
    stats_alpha: float,
    title: str,
    violin_bw_method: float | str | None = None,
    sample_v4_multiple_label: str = SAMPLE_STATUS_LABELS["multiple_v4_reference"],
    sample_v4_single_label: str = SAMPLE_STATUS_LABELS["single_v4_reference"],
    xtick_rotation: float = 35.0,
    ytick_rotation: float = 0.0,
    dot_size: float | None = None,
    dot_opacity: float = 0.48,
    dot_jitter: float = 0.065,
    x_axis_label: str = "Endpoint-pair group",
    y_axis_label: str = "Sample-level summarized absolute pairwise difference",
    pair_groups_to_show: list[str] | None = None,
    endpoint_group_labels: dict[str, str] | None = None,
) -> list[str]:
    fs = _fig21_text_scale(font_size)
    panel_df = sample_summary.copy()
    if "v4_reference_status" not in panel_df.columns:
        panel_df["v4_reference_status"] = "single_v4_reference"
    panel_df["v4_reference_status"] = panel_df["v4_reference_status"].fillna("single_v4_reference").astype(str)
    panel_df["v4_reference_status"] = panel_df["v4_reference_status"].where(
        panel_df["v4_reference_status"].isin(SAMPLE_STATUS_MARKERS),
        "single_v4_reference",
    )

    present_pair_groups = set(panel_df["pair_group_id"].astype(str))
    if pair_groups_to_show:
        preferred = [str(g) for g in pair_groups_to_show if str(g).strip()]
        order = [g for g in preferred if g in present_pair_groups]
    else:
        order = [g for g in PAIR_GROUP_ORDER if g in present_pair_groups]
    if not order:
        ax.axis("off")
        return []

    positions = np.arange(1, len(order) + 1, dtype=float)
    group_data: list[tuple[int, str, pd.DataFrame, np.ndarray]] = []
    for idx, group in enumerate(order):
        group_df = panel_df.loc[
            panel_df["pair_group_id"] == group,
            ["summary_abs_difference", "v4_reference_status"],
        ].copy()
        group_df["summary_abs_difference"] = pd.to_numeric(group_df["summary_abs_difference"], errors="coerce")
        group_df = group_df.dropna(subset=["summary_abs_difference"])
        vals = group_df["summary_abs_difference"].to_numpy(dtype=float)
        if vals.size:
            group_data.append((idx, group, group_df, vals))

    if plot_type == "boxplot":
        for idx, group, _, vals in group_data:
            ax.boxplot(
                [vals],
                positions=[positions[idx]],
                widths=0.55,
                patch_artist=True,
                boxprops=dict(facecolor=PAIR_GROUP_COLORS.get(group, "#999999"), alpha=0.45, edgecolor="#333333"),
                medianprops=dict(color="#111111", linewidth=1.2),
                whiskerprops=dict(color="#333333", linewidth=0.9),
                capprops=dict(color="#333333", linewidth=0.9),
            )
    else:
        violin_vals = [vals for _, _, _, vals in group_data if vals.size >= 2]
        violin_pos = [float(positions[idx]) for idx, _, _, vals in group_data if vals.size >= 2]
        violin_groups = [group for _, group, _, vals in group_data if vals.size >= 2]
        if violin_vals:
            vp = ax.violinplot(
                violin_vals,
                positions=violin_pos,
                widths=0.72,
                showmeans=False,
                showmedians=True,
                showextrema=False,
                bw_method=violin_bw_method,
            )
            for body, group in zip(vp["bodies"], violin_groups):
                body.set_facecolor(PAIR_GROUP_COLORS.get(group, "#999999"))
                body.set_edgecolor("#333333")
                body.set_alpha(0.52)
            if "cmedians" in vp:
                vp["cmedians"].set_color("#111111")
                vp["cmedians"].set_linewidth(1.15)

    rng = np.random.default_rng(21103)
    scatter_size = max(1.0, float(dot_size)) if dot_size is not None else point_size * 1.62
    for idx, group, group_df, vals in group_data:
        xj = positions[idx] + rng.normal(0.0, max(0.0, float(dot_jitter)), size=vals.size)
        status = group_df["v4_reference_status"].astype(str).to_numpy()
        for status_key, marker in SAMPLE_STATUS_MARKERS.items():
            mask = status == status_key
            if not np.any(mask):
                continue
            ax.scatter(
                xj[mask],
                vals[mask],
                s=scatter_size,
                marker=marker,
                color=PAIR_GROUP_COLORS.get(group, "#666666"),
                alpha=max(0.0, min(1.0, float(dot_opacity))),
                edgecolors="#111111",
                linewidths=0.45,
                zorder=3,
            )

    finite_vals = pd.to_numeric(panel_df["summary_abs_difference"], errors="coerce").dropna().to_numpy(dtype=float)
    if finite_vals.size:
        y_min = float(np.min(finite_vals))
        y_max = float(np.max(finite_vals))
        y_range = y_max - y_min
        if y_range <= 0.0:
            y_range = max(0.05, abs(y_max) * 0.25)
        ax.set_ylim(max(0.0, y_min - 0.08 * y_range), y_max + 0.45 * y_range)
    ax.set_xticks(positions)
    ax.set_xticklabels(
        [_pair_group_label(g, endpoint_group_labels) for g in order],
        rotation=xtick_rotation,
        ha="right",
        va="top",
        rotation_mode="anchor",
        fontsize=fs["tick"],
    )
    for tick in ax.get_xticklabels():
        tick.set_linespacing(0.9)
    ax.set_xlabel(x_axis_label, fontsize=fs["base"])
    ax.set_ylabel(y_axis_label, fontsize=fs["base"])
    ax.tick_params(axis="y", labelsize=fs["tick"])
    for tick in ax.get_yticklabels():
        tick.set_rotation(float(ytick_rotation))
        tick.set_va("center_baseline")
        tick.set_ha("right")
        tick.set_rotation_mode("anchor")
    ax.grid(True, axis="y", alpha=0.25, linewidth=0.6)
    ax.set_title(title, fontsize=fs["base"])

    if stats_enabled and isinstance(stats, dict):
        visible_set = set(order)
        visible_planned_rows = [
            r
            for r in planned_sig_rows
            if str(r.get("method_a", "")) in visible_set and str(r.get("method_b", "")) in visible_set
        ]
        bracket_rows: list[dict] = []
        if stats_show_pairwise:
            bracket_rows.extend(visible_planned_rows)
        if show_within_contrast_stats:
            bracket_rows.extend(_panel_b_within_contrast_stat_rows(stats, visible_set))
        if bracket_rows:
            _draw_pairwise_brackets(
                ax=ax,
                pairwise_rows=bracket_rows,
                method_order=order,
                alpha=stats_alpha,
                label_fontsize=fs["pair_label"],
                star_fontsize=fs["pair_star"],
            )
        op = float(stats.get("omnibus_p", np.nan))
        os = float(stats.get("omnibus_stat", np.nan))
        ax.text(
            0.98,
            0.03,
            f"paired permutation\nomnibus={os:.3g}, p={op:.3g}" if np.isfinite(op) else "paired permutation\nomnibus=NA, p=NA",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=fs["stats"],
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="#bbbbbb", alpha=0.88),
        )

    status_handles = [
        Line2D(
            [0],
            [0],
            marker=SAMPLE_STATUS_MARKERS["multiple_v4_reference"],
            linestyle="None",
            markerfacecolor="#666666",
            markeredgecolor="#111111",
            markersize=9.75,
            label=sample_v4_multiple_label,
        ),
        Line2D(
            [0],
            [0],
            marker=SAMPLE_STATUS_MARKERS["single_v4_reference"],
            linestyle="None",
            markerfacecolor="#666666",
            markeredgecolor="#111111",
            markersize=9.75,
            label=sample_v4_single_label,
        ),
    ]
    ax.legend(
        handles=status_handles,
        title="Sample V4 status",
        loc="upper left",
        frameon=True,
        framealpha=0.92,
        facecolor="white",
        edgecolor="#BBBBBB",
        fontsize=fs["legend"],
        title_fontsize=fs["legend_title"],
    )
    return order


def _draw_pair_group_raincloud_grid_panel(
    fig,
    parent_spec,
    pairwise_df: pd.DataFrame,
    sample_summary: pd.DataFrame | None,
    font_size: float,
    title: str,
    sample_v4_multiple_label: str = SAMPLE_STATUS_LABELS["multiple_v4_reference"],
    sample_v4_single_label: str = SAMPLE_STATUS_LABELS["single_v4_reference"],
    abs_density_label: str = "Pairwise abs_difference density",
    summary_density_label: str = "Sample summary_abs_difference density",
    show_summary_density: bool = True,
    mean_label: str = "Mean",
    grid_label_font_size: float | None = None,
    grid_tick_font_size: float | None = None,
    grid_text_font_size: float | None = None,
    grid_legend_font_size: float | None = None,
    dot_size: float | None = None,
    dot_opacity: float = 0.95,
    dot_jitter: float = 0.055,
    x_axis_label: str | None = None,
    y_axis_label: str | None = None,
) -> pd.DataFrame:
    fs = _fig21_text_scale(font_size)
    groups = list(ENDPOINT_GROUP_ORDER)
    n = len(groups)
    inner = parent_spec.subgridspec(
        n + 1,
        n,
        height_ratios=[0.48] + [1.0] * n,
        hspace=0.30,
        wspace=0.12,
    )

    label_fs = (
        max(4.0, float(grid_label_font_size))
        if grid_label_font_size is not None
        else fs["tick"]
    )
    tick_fs = (
        max(4.0, float(grid_tick_font_size))
        if grid_tick_font_size is not None
        else fs["tick"]
    )
    text_fs = (
        max(4.0, float(grid_text_font_size))
        if grid_text_font_size is not None
        else fs["stats"]
    )
    legend_fs = (
        max(4.0, float(grid_legend_font_size))
        if grid_legend_font_size is not None
        else fs["legend"]
    )
    sample_dot_size = max(1.0, float(dot_size)) if dot_size is not None else 25.5
    sample_dot_opacity = max(0.0, min(1.0, float(dot_opacity)))
    sample_dot_jitter = max(0.0, float(dot_jitter))
    sample_summary_df = sample_summary.copy() if isinstance(sample_summary, pd.DataFrame) else pd.DataFrame()
    if not sample_summary_df.empty:
        required_cols = {"pair_group_id", "summary_abs_difference"}
        if not required_cols.issubset(set(sample_summary_df.columns)):
            sample_summary_df = pd.DataFrame()
        else:
            if "v4_reference_status" not in sample_summary_df.columns:
                sample_summary_df["v4_reference_status"] = "single_v4_reference"
            sample_summary_df["v4_reference_status"] = sample_summary_df["v4_reference_status"].fillna("single_v4_reference").astype(str)
            sample_summary_df["v4_reference_status"] = sample_summary_df["v4_reference_status"].where(
                sample_summary_df["v4_reference_status"].isin(SAMPLE_STATUS_MARKERS),
                "single_v4_reference",
            )
            sample_summary_df["summary_abs_difference"] = pd.to_numeric(
                sample_summary_df["summary_abs_difference"], errors="coerce"
            )
            sample_summary_df = sample_summary_df.dropna(subset=["summary_abs_difference"])
    cell_values: dict[str, np.ndarray] = {}
    all_vals: list[np.ndarray] = []
    for ga in groups:
        for gb in groups:
            pg = _pair_group_id(ga, gb)
            vals = (
                pd.to_numeric(
                    pairwise_df.loc[pairwise_df["pair_group_id"].astype(str) == pg, "abs_difference"],
                    errors="coerce",
                )
                .dropna()
                .to_numpy(dtype=float)
            )
            vals = vals[np.isfinite(vals)]
            cell_values[pg] = vals
            if vals.size:
                all_vals.append(vals)

    if all_vals:
        pooled = np.concatenate(all_vals)
        pooled = pooled[np.isfinite(pooled)]
        if pooled.size:
            vmax = float(np.quantile(pooled, 0.995))
            vmax = max(vmax, float(np.max(pooled)) * 0.75, 1e-6)
        else:
            vmax = 1.0
    else:
        vmax = 1.0

    header_ax = fig.add_subplot(inner[0, :])
    header_ax.axis("off")
    header_ax.text(
        0.5,
        0.93,
        title,
        transform=header_ax.transAxes,
        ha="center",
        va="top",
        fontsize=fs["base"],
    )

    summary_rows: list[dict] = []
    for i, ga in enumerate(groups):
        for j, gb in enumerate(groups):
            ax = fig.add_subplot(inner[i + 1, j])
            pg = _pair_group_id(ga, gb)
            vals = cell_values.get(pg, np.array([], dtype=float))

            ax.set_xlim(0.0, vmax)
            ax.set_yticks([])
            if i < n - 1:
                ax.set_xticklabels([])
            else:
                ax.tick_params(axis="x", labelsize=tick_fs)
            ax.tick_params(axis="x", length=2.0, pad=1.0)

            if i == 0:
                ax.set_title(
                    ENDPOINT_GROUP_LABELS.get(gb, gb),
                    fontsize=label_fs,
                    pad=3.0,
                )
            if j == 0:
                ax.set_ylabel(
                    ENDPOINT_GROUP_LABELS.get(ga, ga),
                    fontsize=label_fs,
                    rotation=90.0,
                    labelpad=18.0,
                    ha="center",
                    va="center",
                )
            if x_axis_label and i == n - 1 and j == n // 2:
                ax.set_xlabel(str(x_axis_label), fontsize=label_fs)

            if vals.size >= 2:
                cell_color = PAIR_GROUP_COLORS.get(pg, "#9A9A9A")
                summary_vals = np.array([], dtype=float)
                if not sample_summary_df.empty:
                    summary_vals = (
                        pd.to_numeric(
                            sample_summary_df.loc[
                                sample_summary_df["pair_group_id"].astype(str) == pg,
                                "summary_abs_difference",
                            ],
                            errors="coerce",
                        )
                        .dropna()
                        .to_numpy(dtype=float)
                    )
                    summary_vals = summary_vals[np.isfinite(summary_vals)]
                if sns is not None:
                    kde_kwargs = dict(
                        x=vals,
                        fill=True,
                        cut=0,
                        clip=(0.0, vmax),
                        bw_adjust=0.75,
                        color=cell_color,
                        alpha=0.48,
                        linewidth=1.0,
                        ax=ax,
                    )
                    with warnings.catch_warnings():
                        warnings.filterwarnings(
                            "ignore",
                            message=r".*use_inf_as_na option is deprecated.*",
                            category=FutureWarning,
                            module=r"seaborn\._oldcore",
                        )
                        try:
                            sns.kdeplot(**kde_kwargs, warn_singular=False)
                        except TypeError:
                            sns.kdeplot(**kde_kwargs)
                        if show_summary_density and summary_vals.size >= 2:
                            try:
                                sns.kdeplot(
                                    x=summary_vals,
                                    fill=False,
                                    cut=0,
                                    clip=(0.0, vmax),
                                    bw_adjust=0.75,
                                    color="#0b4f8a",
                                    alpha=0.98,
                                    linewidth=1.35,
                                    linestyle="-",
                                    ax=ax,
                                    warn_singular=False,
                                )
                            except TypeError:
                                sns.kdeplot(
                                    x=summary_vals,
                                    fill=False,
                                    cut=0,
                                    clip=(0.0, vmax),
                                    bw_adjust=0.75,
                                    color="#0b4f8a",
                                    alpha=0.98,
                                    linewidth=1.35,
                                    linestyle="-",
                                    ax=ax,
                                )
                else:
                    hist, edges = np.histogram(vals, bins=30, range=(0.0, vmax), density=True)
                    centers = 0.5 * (edges[:-1] + edges[1:])
                    ax.fill_between(centers, hist, color=cell_color, alpha=0.42, linewidth=0.0, zorder=2)
                    ax.plot(centers, hist, color="#333333", linewidth=0.8, zorder=3)
                    if show_summary_density and summary_vals.size >= 2:
                        shist, sedges = np.histogram(summary_vals, bins=24, range=(0.0, vmax), density=True)
                        scenters = 0.5 * (sedges[:-1] + sedges[1:])
                        ax.plot(scenters, shist, color="#0b4f8a", linewidth=1.35, linestyle="-", zorder=3.5)

                y1 = float(ax.get_ylim()[1]) if np.isfinite(ax.get_ylim()[1]) else 1.0
                if y1 <= 0.0:
                    y1 = 1.0
                ax.set_ylim(0.0, y1 * 1.06)

                # Keep boxplot horizontally aligned with KDE to avoid visible
                # axis-start offsets between the two layers.
                box_ax = ax.inset_axes([0.01, 0.02, 0.98, 0.42], sharex=ax)
                box_ax.set_facecolor("none")
                if sns is not None:
                    box_df = pd.DataFrame({"abs_difference": vals})
                    with warnings.catch_warnings():
                        warnings.filterwarnings(
                            "ignore",
                            message=r".*use_inf_as_na option is deprecated.*",
                            category=FutureWarning,
                            module=r"seaborn\._oldcore",
                        )
                        sns.boxplot(
                            data=box_df,
                            x="abs_difference",
                            orient="h",
                            color=cell_color,
                            width=0.75,
                            saturation=0.9,
                            linewidth=1.2,
                            fliersize=0.0,
                            whis=1.5,
                            ax=box_ax,
                        )
                else:
                    box_ax.boxplot(
                        [vals],
                        vert=False,
                        widths=0.675,
                        showfliers=False,
                        patch_artist=True,
                        boxprops=dict(facecolor=cell_color, edgecolor="#2A2A2A", linewidth=1.2, alpha=0.7),
                        medianprops=dict(color="#111111", linewidth=1.5),
                        whiskerprops=dict(color="#2A2A2A", linewidth=1.125),
                        capprops=dict(color="#2A2A2A", linewidth=1.125),
                    )
                box_ax.set_yticks([])
                box_ax.set_ylabel("")
                box_ax.set_xlabel("")
                box_ax.tick_params(axis="x", labelbottom=False, length=0)
                for spine in box_ax.spines.values():
                    spine.set_linewidth(0.4)
                    spine.set_color("#BBBBBB")
                box_ax.grid(True, axis="x", alpha=0.12, linewidth=0.4)
                box_ax.set_xlim(0.0, vmax)

                if not sample_summary_df.empty:
                    sample_points = sample_summary_df.loc[
                        sample_summary_df["pair_group_id"].astype(str) == pg,
                        ["summary_abs_difference", "v4_reference_status"],
                    ].copy()
                    if not sample_points.empty:
                        sv = sample_points["summary_abs_difference"].to_numpy(dtype=float)
                        st = sample_points["v4_reference_status"].astype(str).to_numpy()
                        rng = np.random.default_rng(21127 + i * 97 + j * 13)
                        yj = np.clip(
                            rng.normal(loc=0.0, scale=sample_dot_jitter, size=sv.size),
                            -0.20,
                            0.20,
                        )
                        for status_key, marker in SAMPLE_STATUS_MARKERS.items():
                            mask = st == status_key
                            if not np.any(mask):
                                continue
                            box_ax.scatter(
                                sv[mask],
                                yj[mask],
                                marker=marker,
                                s=sample_dot_size,
                                facecolors="white",
                                edgecolors="#111111",
                                linewidths=0.75,
                                alpha=sample_dot_opacity,
                                zorder=6,
                            )

                mean_v = float(np.mean(vals))
                ax.axvline(mean_v, color="#1f77b4", linewidth=1.0, linestyle="--", zorder=4, alpha=0.95)
                y_marker = float(ax.get_ylim()[1]) * 0.96
                ax.scatter([mean_v], [y_marker], marker="D", s=18.0, color="#1f77b4", edgecolors="white", linewidths=0.3, zorder=5)

                summary_rows.append(
                    {
                        "row_group_id": ga,
                        "col_group_id": gb,
                        "row_group_label": ENDPOINT_GROUP_LABELS.get(ga, ga),
                        "col_group_label": ENDPOINT_GROUP_LABELS.get(gb, gb),
                        "pair_group_id": pg,
                        "pair_group_label": _pair_group_full_label(pg),
                        "n_values": int(vals.size),
                        "mean_abs_difference": float(np.mean(vals)),
                        "median_abs_difference": float(np.median(vals)),
                        "sd_abs_difference": float(np.std(vals, ddof=1)) if vals.size > 1 else 0.0,
                        "q25_abs_difference": float(np.quantile(vals, 0.25)),
                        "q75_abs_difference": float(np.quantile(vals, 0.75)),
                        "min_abs_difference": float(np.min(vals)),
                        "max_abs_difference": float(np.max(vals)),
                    }
                )
            else:
                ax.text(
                    0.5,
                    0.52,
                    "NA" if vals.size == 0 else "n<2",
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                    fontsize=text_fs,
                    color="#666666",
                )
                summary_rows.append(
                    {
                        "row_group_id": ga,
                        "col_group_id": gb,
                        "row_group_label": ENDPOINT_GROUP_LABELS.get(ga, ga),
                        "col_group_label": ENDPOINT_GROUP_LABELS.get(gb, gb),
                        "pair_group_id": pg,
                        "pair_group_label": _pair_group_full_label(pg),
                        "n_values": int(vals.size),
                        "mean_abs_difference": np.nan,
                        "median_abs_difference": np.nan,
                        "sd_abs_difference": np.nan,
                        "q25_abs_difference": np.nan,
                        "q75_abs_difference": np.nan,
                        "min_abs_difference": np.nan,
                        "max_abs_difference": np.nan,
                    }
                )

            for spine in ax.spines.values():
                spine.set_visible(True)
                spine.set_linewidth(0.8)
                spine.set_color("#8F8F8F")
            ax.grid(True, axis="x", alpha=0.14, linewidth=0.45)

    legend_handles = [
        Patch(facecolor="#9A9A9A", edgecolor="#9A9A9A", alpha=0.48, label=abs_density_label),
        Line2D([0], [0], color="#1f77b4", linestyle="--", linewidth=1.1, marker="D", markersize=4.5, label=mean_label),
        Line2D(
            [0],
            [0],
            linestyle="None",
            marker=SAMPLE_STATUS_MARKERS["multiple_v4_reference"],
            markerfacecolor="white",
            markeredgecolor="#111111",
            markersize=5.0,
            label=sample_v4_multiple_label,
        ),
        Line2D(
            [0],
            [0],
            linestyle="None",
            marker=SAMPLE_STATUS_MARKERS["single_v4_reference"],
            markerfacecolor="white",
            markeredgecolor="#111111",
            markersize=5.0,
            label=sample_v4_single_label,
        ),
    ]
    if show_summary_density:
        legend_handles.insert(
            1,
            Line2D([0], [0], color="#0b4f8a", linestyle="-", linewidth=1.35, label=summary_density_label),
        )
    header_ax.legend(
        handles=legend_handles,
        loc="center right",
        bbox_to_anchor=(1.0, 0.22),
        borderaxespad=0.0,
        frameon=True,
        framealpha=0.92,
        fontsize=legend_fs,
        ncol=1,
        handlelength=1.5,
        columnspacing=0.9,
        borderpad=0.22,
    )
    if y_axis_label:
        header_ax.text(
            0.01,
            0.14,
            str(y_axis_label),
            transform=header_ax.transAxes,
            ha="left",
            va="center",
            rotation=90,
            fontsize=label_fs,
        )

    return pd.DataFrame(summary_rows)


def _draw_endpoint_matrix_panel(
    ax,
    endpoints_df: pd.DataFrame,
    method_order: list[str],
    metric_order: list[str],
    agnostic_metrics: list[str],
    aware_metrics: list[str],
    gap_cells: int,
    font_size: float,
    title: str,
    summary_stat: str = "median",
    order_mode: str = "metric_family_first",
    show_xticklabels: bool = True,
    xtick_rotation: float = 90.0,
    ytick_rotation: float = 0.0,
    x_axis_label: str | None = None,
    y_axis_label: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    endpoint_matrix, endpoint_meta = _endpoint_distance_heatmap_matrix(
        endpoints_df=endpoints_df,
        method_order=method_order,
        metric_order=metric_order,
        agnostic_metrics=agnostic_metrics,
        aware_metrics=aware_metrics,
        summary_stat=summary_stat,
        order_mode=order_mode,
    )
    if endpoint_matrix.empty or endpoint_meta.empty:
        ax.axis("off")
        ax.text(
            0.5,
            0.5,
            "Endpoint-level matrix unavailable",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=max(7.0, font_size - 2.0),
        )
        return endpoint_matrix, endpoint_meta

    endpoint_ids = endpoint_meta["endpoint_id"].astype(str).tolist()
    endpoint_id_to_label = dict(
        zip(endpoint_meta["endpoint_id"].astype(str), endpoint_meta["endpoint_label"].astype(str))
    )
    arr_ep, ticks, _, _, method_seq, boundaries = _expand_endpoint_matrix_with_gaps(
        matrix=endpoint_matrix,
        endpoint_meta=endpoint_meta,
        gap_cells=int(gap_cells),
    )
    finite_ep = arr_ep[np.isfinite(arr_ep)]
    vmax_ep = float(np.max(finite_ep)) if finite_ep.size else 1.0
    if vmax_ep <= 0.0:
        vmax_ep = 1.0
    cmap = copy.copy(plt.get_cmap("magma_r"))
    cmap.set_bad(color="white")
    im_ep = ax.imshow(np.ma.masked_invalid(arr_ep), cmap=cmap, vmin=0.0, vmax=vmax_ep, aspect="equal")
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    x_labels = [endpoint_id_to_label.get(x, x) for x in endpoint_ids]
    y_labels = [endpoint_id_to_label.get(x, x) for x in endpoint_ids]
    xrot = float(xtick_rotation) % 360.0
    if np.isclose(xrot, 0.0) or np.isclose(xrot, 180.0):
        x_ha = "center"
    elif 0.0 < xrot < 180.0:
        x_ha = "right"
    else:
        x_ha = "left"
    ax.set_xticklabels(
        x_labels,
        rotation=xtick_rotation,
        ha=x_ha,
        va="top",
        fontsize=max(3.2, font_size - 8.8),
        rotation_mode="anchor",
    )
    ax.set_yticklabels(
        y_labels,
        rotation=ytick_rotation,
        va="center_baseline",
        ha="right",
        fontsize=max(3.2, font_size - 8.8),
        rotation_mode="anchor",
    )
    ax.tick_params(axis="x", pad=1.0)
    ax.tick_params(axis="y", pad=1.5)
    if not show_xticklabels:
        ax.tick_params(axis="x", labelbottom=False)
        for tick in ax.get_xticklabels():
            tick.set_visible(False)
    ax.tick_params(length=0)
    ax.set_title(title, fontsize=font_size)
    if x_axis_label:
        ax.set_xlabel(str(x_axis_label), fontsize=font_size)
    if y_axis_label:
        ax.set_ylabel(str(y_axis_label), fontsize=font_size)

    if int(gap_cells) <= 0:
        draw_boundaries: list[float] = []
        for i in range(1, len(method_seq)):
            if method_seq[i] != method_seq[i - 1]:
                draw_boundaries.append(i - 0.5)
        for b in draw_boundaries:
            ax.axvline(b, color="white", linewidth=0.6, alpha=0.85)
            ax.axhline(b, color="white", linewidth=0.6, alpha=0.85)

    for spine in ax.spines.values():
        spine.set_visible(False)
    return endpoint_matrix, endpoint_meta


def _draw_endpoint_matrix_dual_panel(
    fig,
    parent_spec,
    endpoints_df: pd.DataFrame,
    method_order: list[str],
    metric_order: list[str],
    agnostic_metrics: list[str],
    aware_metrics: list[str],
    gap_cells: int,
    font_size: float,
    title: str,
    summary_stat: str = "median",
    xtick_rotation: float = 90.0,
    ytick_rotation: float = 0.0,
    panel_c_a_title: str | None = None,
    panel_c_b_title: str | None = None,
    panel_c_a_xtick_rotation: float | None = None,
    panel_c_a_ytick_rotation: float | None = None,
    panel_c_b_xtick_rotation: float | None = None,
    panel_c_b_ytick_rotation: float | None = None,
    x_axis_label: str | None = None,
    y_axis_label: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    fs = _fig21_text_scale(font_size)
    inner = parent_spec.subgridspec(
        2,
        2,
        width_ratios=[1.0, 0.045],
        height_ratios=[1.0, 1.0],
        hspace=0.34,
        wspace=0.04,
    )
    ax_top = fig.add_subplot(inner[0, 0])
    ax_bottom = fig.add_subplot(inner[1, 0])
    cax = fig.add_subplot(inner[:, 1])

    top_title = str(panel_c_a_title or f"{title}\nAgnostic-first ordering")
    bottom_title = str(panel_c_b_title or "Denoising-first ordering")
    top_xrot = float(xtick_rotation if panel_c_a_xtick_rotation is None else panel_c_a_xtick_rotation)
    top_yrot = float(ytick_rotation if panel_c_a_ytick_rotation is None else panel_c_a_ytick_rotation)
    bottom_xrot = float(xtick_rotation if panel_c_b_xtick_rotation is None else panel_c_b_xtick_rotation)
    bottom_yrot = float(ytick_rotation if panel_c_b_ytick_rotation is None else panel_c_b_ytick_rotation)
    top_matrix, top_meta = _draw_endpoint_matrix_panel(
        ax_top,
        endpoints_df=endpoints_df,
        method_order=method_order,
        metric_order=metric_order,
        agnostic_metrics=agnostic_metrics,
        aware_metrics=aware_metrics,
        gap_cells=gap_cells,
        font_size=font_size,
        title=top_title,
        summary_stat=summary_stat,
        order_mode="metric_family_first",
        show_xticklabels=True,
        xtick_rotation=top_xrot,
        ytick_rotation=top_yrot,
        x_axis_label=x_axis_label,
        y_axis_label=y_axis_label,
    )
    bottom_matrix, bottom_meta = _draw_endpoint_matrix_panel(
        ax_bottom,
        endpoints_df=endpoints_df,
        method_order=method_order,
        metric_order=metric_order,
        agnostic_metrics=agnostic_metrics,
        aware_metrics=aware_metrics,
        gap_cells=gap_cells,
        font_size=font_size,
        title=bottom_title,
        summary_stat=summary_stat,
        order_mode="method_family_first",
        show_xticklabels=True,
        xtick_rotation=bottom_xrot,
        ytick_rotation=bottom_yrot,
        x_axis_label=x_axis_label,
        y_axis_label=y_axis_label,
    )
    if ax_top.images:
        cbar = fig.colorbar(ax_top.images[0], cax=cax)
        cbar.set_label(
            _summary_stat_colorbar_label(summary_stat=summary_stat, summarized=False),
            fontsize=fs["cbar_label"],
        )
        cbar.ax.tick_params(labelsize=fs["cbar_tick"])
    return top_matrix, top_meta, bottom_matrix, bottom_meta



def _write_pair_group_contrast_stats(
    sample_summary: pd.DataFrame,
    out_dir: Path,
    stats_enabled: bool,
    stats_permutations: int,
    stats_seed: int,
    stats_alpha: float,
    stats_pairwise_test: str,
    stats_multiple_test_correction: bool,
) -> tuple[dict | None, list[dict]]:
    if not stats_enabled or sample_summary.empty:
        return None, []

    methods = [g for g in PAIR_GROUP_ORDER if g in set(sample_summary["pair_group_id"].astype(str))]
    stats_df = sample_summary.rename(
        columns={"pair_group_id": "method", "summary_abs_difference": "distance"}
    ).copy()
    st = _compute_panel_stats_permutation(
        df=stats_df,
        method_order=methods,
        n_perm=stats_permutations,
        seed=stats_seed + 21071,
        alpha=stats_alpha,
        pairwise_test=stats_pairwise_test,
        multiple_test_correction=stats_multiple_test_correction,
    )
    pd.DataFrame(
        [
            {
                "metric": "pair_group_distance_contrast",
                "label": "Pair-group distance contrast",
                "stats_model": st.get("stats_model", "paired_permutation"),
                "pairwise_test": stats_pairwise_test,
                "multiple_test_correction": bool(
                    st.get("multiple_test_correction", stats_multiple_test_correction)
                ),
                "n_blocks": st.get("n_blocks", 0),
                "n_pair_groups": len(methods),
                "omnibus_stat": st.get("omnibus_stat", np.nan),
                "omnibus_p": st.get("omnibus_p", np.nan),
                "omnibus_significant_alpha": bool(
                    np.isfinite(float(st.get("omnibus_p", np.nan)))
                    and float(st.get("omnibus_p", np.nan)) < stats_alpha
                ),
                "alpha": stats_alpha,
                "n_permutations": stats_permutations,
            }
        ]
    ).to_csv(out_dir / "pair_group_distance_contrast_stats_omnibus.tsv", sep="\t", index=False)

    pairwise_rows: list[dict] = []
    for row in st.get("pairwise", []):
        method_a = str(row.get("method_a", ""))
        method_b = str(row.get("method_b", ""))
        pairwise_rows.append(
            {
                "metric": "pair_group_distance_contrast",
                "label": "Pair-group distance contrast",
                "method_a": method_a,
                "method_b": method_b,
                "method_a_label": _pair_group_full_label(method_a),
                "method_b_label": _pair_group_full_label(method_b),
                "stats_model": st.get("stats_model", "paired_permutation"),
                "pairwise_test": row.get("pairwise_test", stats_pairwise_test),
                "statistic": row.get("statistic", row.get("w_plus", np.nan)),
                "w_plus": row.get("w_plus", np.nan),
                "n_nonzero_pairs": row.get("n_nonzero_pairs", np.nan),
                "n_pos_pairs": row.get("n_pos_pairs", np.nan),
                "n_neg_pairs": row.get("n_neg_pairs", np.nan),
                "n_zero_pairs": row.get("n_zero_pairs", np.nan),
                "median_diff_a_minus_b": row.get("median_diff_a_minus_b", np.nan),
                "estimate": row.get("estimate", row.get("median_diff_a_minus_b", np.nan)),
                "std_error": row.get("std_error", np.nan),
                "model_family": row.get("model_family", ""),
                "p_raw": row.get("p_raw", np.nan),
                "p_adj_bh": row.get("p_adj_bh", np.nan),
                "p_value_used": row.get("p_value_used", row.get("p_adj_bh", np.nan)),
                "p_adjustment": row.get("p_adjustment", "bh"),
                "multiple_test_correction": bool(
                    row.get("multiple_test_correction", stats_multiple_test_correction)
                ),
                "significant_alpha": row.get("significant_alpha", False),
                "significant_bh_alpha": row.get("significant_bh", False),
                "alpha": stats_alpha,
                "n_permutations": stats_permutations,
            }
        )
    pairwise_df = pd.DataFrame(pairwise_rows)
    pairwise_df.to_csv(out_dir / "pair_group_distance_contrast_stats_pairwise.tsv", sep="\t", index=False)

    within_rows: list[dict] = []
    if not pairwise_df.empty:
        within_groups = set(WITHIN_PAIR_GROUP_ORDER)
        within_mask = pairwise_df["method_a"].isin(within_groups) & pairwise_df["method_b"].isin(within_groups)
        within_rows = pairwise_df.loc[within_mask].copy().to_dict(orient="records")
    pd.DataFrame(within_rows).to_csv(
        out_dir / "pair_group_distance_contrast_stats_within.tsv",
        sep="\t",
        index=False,
    )

    planned_rows: list[dict] = []
    for method_a, method_b, comparison in PLANNED_PAIR_GROUP_COMPARISONS:
        if pairwise_df.empty:
            continue
        mask = (
            ((pairwise_df["method_a"] == method_a) & (pairwise_df["method_b"] == method_b))
            | ((pairwise_df["method_a"] == method_b) & (pairwise_df["method_b"] == method_a))
        )
        if not mask.any():
            continue
        rec = pairwise_df[mask].iloc[0].to_dict()
        rec["planned_comparison"] = comparison
        rec["planned_comparison_label"] = comparison.replace("_", " ")
        planned_rows.append(rec)
    pd.DataFrame(planned_rows).to_csv(out_dir / "pair_group_distance_contrast_stats_planned.tsv", sep="\t", index=False)
    planned_sig = [
        row
        for row in planned_rows
        if bool(row.get("significant_alpha", False))
        and np.isfinite(float(row.get("p_value_used", np.nan)))
        and float(row.get("p_value_used", np.nan)) < stats_alpha
    ]
    return st, planned_sig


def _write_distance_deflation_stats(
    deflation: pd.DataFrame,
    out_dir: Path,
    stats_enabled: bool,
    stats_permutations: int,
    stats_seed: int,
    stats_alpha: float,
    stats_pairwise_test: str,
    stats_multiple_test_correction: bool,
) -> list[dict]:
    if not stats_enabled or deflation.empty:
        pd.DataFrame().to_csv(out_dir / "distance_deflation_by_metric_family_stats.tsv", sep="\t", index=False)
        return []

    rows: list[dict] = []
    for idx, family in enumerate(DISTANCE_DEFLATION_FAMILY_ORDER, start=1):
        vals = pd.to_numeric(
            deflation.loc[
                deflation["metric_family"].astype(str) == family,
                "delta_denoising_minus_clustering",
            ],
            errors="coerce",
        ).dropna().to_numpy(dtype=float)
        if vals.size == 0:
            continue
        zeros = np.zeros_like(vals, dtype=float)
        d_all = vals[np.isfinite(vals)]
        n_zero = int(np.sum(d_all == 0.0))
        n_pos = int(np.sum(d_all > 0.0))
        n_neg = int(np.sum(d_all < 0.0))
        if stats_pairwise_test == "sign_test_exact":
            stat, p_raw, n_nonzero, n_pos_test, n_neg_test, n_zero_test, median_delta = _paired_sign_test_exact(
                vals,
                zeros,
            )
            w_plus = float("nan")
            n_pos = int(n_pos_test)
            n_neg = int(n_neg_test)
            n_zero = int(n_zero_test)
        else:
            stat, p_raw, n_nonzero, median_delta = _paired_signed_rank_permutation(
                vals,
                zeros,
                n_perm=stats_permutations,
                seed=stats_seed + 22001 + idx,
            )
            w_plus = stat
        rows.append(
            {
                "metric": "distance_deflation_by_metric_family",
                "metric_family": family,
                "metric_family_label": DISTANCE_DEFLATION_FAMILY_LABELS.get(family, family),
                "contrast": "denoising_mean_distance_minus_clustering_mean_distance",
                "null_hypothesis": "delta_equals_zero",
                "stats_model": "one_sample_paired_permutation",
                "pairwise_test": stats_pairwise_test,
                "statistic": stat,
                "w_plus": w_plus,
                "n_samples": int(vals.size),
                "n_nonzero_pairs": int(n_nonzero),
                "n_pos_pairs": n_pos,
                "n_neg_pairs": n_neg,
                "n_zero_pairs": n_zero,
                "mean_delta": float(np.mean(vals)),
                "median_delta": float(median_delta),
                "estimate": float(median_delta),
                "std_error": float(np.std(vals, ddof=1) / np.sqrt(vals.size)) if vals.size > 1 else 0.0,
                "model_family": "permutation",
                "p_raw": p_raw,
                "alpha": stats_alpha,
                "n_permutations": stats_permutations,
            }
        )

    p_raw_vals = [float(r.get("p_raw", np.nan)) for r in rows]
    p_adj_vals = _bh_adjust(p_raw_vals)
    p_used_vals = p_adj_vals if stats_multiple_test_correction else p_raw_vals
    for i, row in enumerate(rows):
        p_adj = p_adj_vals[i]
        p_used = p_used_vals[i]
        row["p_adj_bh"] = p_adj
        row["p_value_used"] = p_used
        row["p_adjustment"] = "bh" if stats_multiple_test_correction else "none"
        row["multiple_test_correction"] = bool(stats_multiple_test_correction)
        row["significant_alpha"] = bool(np.isfinite(p_used) and p_used < stats_alpha)
        row["significant_bh_alpha"] = bool(np.isfinite(p_adj) and p_adj < stats_alpha)

    pd.DataFrame(rows).to_csv(out_dir / "distance_deflation_by_metric_family_stats.tsv", sep="\t", index=False)
    return rows


def _plot_pair_group_distance_contrast(
    sample_summary: pd.DataFrame,
    pairwise_df: pd.DataFrame,
    endpoints_df: pd.DataFrame,
    sens: pd.DataFrame,
    benchmark_outdir: Path,
    out_dir: Path,
    title: str,
    panel_a_title: str,
    panel_a_a_title: str,
    panel_a_b_title: str,
    panel_b_title: str,
    panel_c_title: str,
    panel_c_a_title: str,
    panel_c_b_title: str,
    panel_d_title: str,
    panel_e_title: str,
    panel_f_title: str,
    font_size: float,
    panel_a_font_size: float,
    panel_b_font_size: float,
    panel_c_font_size: float,
    panel_d_font_size: float,
    panel_e_font_size: float,
    panel_f_font_size: float,
    panel_e_grid_label_font_size: float | None,
    panel_e_grid_tick_font_size: float | None,
    panel_e_grid_text_font_size: float | None,
    panel_e_grid_legend_font_size: float | None,
    panel_e_abs_density_label: str,
    panel_e_summary_density_label: str,
    panel_e_show_summary_density: bool,
    panel_e_mean_label: str,
    panel_e_sample_v4_multiple_label: str,
    panel_e_sample_v4_single_label: str,
    panel_a_x_axis_label: str | None,
    panel_a_y_axis_label: str | None,
    panel_b_x_axis_label: str,
    panel_b_y_axis_label: str,
    panel_c_x_axis_label: str | None,
    panel_c_y_axis_label: str | None,
    panel_d_x_axis_label: str,
    panel_d_y_axis_label: str,
    panel_e_x_axis_label: str | None,
    panel_e_y_axis_label: str | None,
    panel_f_x_axis_label: str,
    panel_f_y_axis_label: str,
    panel_f_xtick_rotation: float,
    panel_b_pair_groups_to_show: list[str] | None,
    panel_b_group_labels: dict[str, str],
    show_panel_a: bool,
    show_panel_b: bool,
    show_panel_c: bool,
    show_panel_d: bool,
    show_panel_e: bool,
    show_panel_f: bool,
    sample_v4_multiple_label: str,
    sample_v4_single_label: str,
    panel_a_xtick_rotation: float,
    panel_a_ytick_rotation: float,
    panel_b_xtick_rotation: float,
    panel_b_ytick_rotation: float,
    endpoint_matrix_gap_cells: int,
    endpoint_matrix_xtick_rotation: float,
    endpoint_matrix_ytick_rotation: float,
    panel_c_a_xtick_rotation: float | None,
    panel_c_a_ytick_rotation: float | None,
    panel_c_b_xtick_rotation: float | None,
    panel_c_b_ytick_rotation: float | None,
    pair_group_summary_stat: str,
    panel_a_show_planned_comparison_annotations: bool,
    show_panel_e_raincloud_grid: bool,
    panel_b_dot_size: float | None,
    panel_b_dot_opacity: float,
    panel_b_dot_jitter: float,
    panel_d_dot_size: float | None,
    panel_d_dot_opacity: float,
    panel_d_dot_jitter: float,
    panel_e_dot_size: float | None,
    panel_e_dot_opacity: float,
    panel_e_dot_jitter: float,
    panel_f_dot_size: float | None,
    panel_f_dot_opacity: float,
    panel_f_dot_jitter: float,
    point_size: float,
    plot_type: str,
    method_order: list[str],
    panel_d_methods: list[str],
    panel_f_methods: list[str],
    panel_f_summary_statistic: str,
    panel_f_stats_show_annotations: bool,
    panel_f_stats_show_non_significant_annotations: bool,
    panel_d_stats_show_annotations: bool,
    metric_order: list[str],
    agnostic_metrics: list[str],
    aware_metrics: list[str],
    violin_bw_method: float | str | None,
    stats: dict | None,
    directional_stats: dict | None,
    deflation_points: pd.DataFrame,
    deflation_stats: list[dict],
    planned_sig_rows: list[dict],
    stats_enabled: bool,
    stats_show_pairwise: bool,
    panel_b_show_within_contrast_stats: bool,
    stats_alpha: float,
) -> None:
    if sample_summary.empty:
        return
    fs_panel_a = _fig21_text_scale(panel_a_font_size)
    fs_panel_c = _fig21_text_scale(panel_c_font_size)
    fs_fig = _fig21_text_scale(font_size)
    summary_stat = _normalize_summary_stat(pair_group_summary_stat)

    sample_summary_panel_b = _attach_sample_v4_status(sample_summary, benchmark_outdir)
    sample_summary_panel_b.to_csv(out_dir / "pair_group_distance_contrast_sample_summary_with_v4_status.tsv", sep="\t", index=False)
    sens_panel_d = _attach_sample_v4_status(
        sens[sens["method"].astype(str).isin(panel_d_methods)].copy(), benchmark_outdir
    )
    sens_panel_d.to_csv(out_dir / "directional_sensitivity_method_points_with_v4_status.tsv", sep="\t", index=False)
    deflation_panel_f = _attach_sample_v4_status(deflation_points, benchmark_outdir)
    if not deflation_panel_f.empty:
        deflation_panel_f.to_csv(out_dir / "distance_by_method_metric_family_points_with_v4_status.tsv", sep="\t", index=False)

    order = [g for g in PAIR_GROUP_ORDER if g in set(sample_summary["pair_group_id"].astype(str))]
    if not order:
        return
    matrix_a_a = _panel_a_heatmap_matrix_from_endpoints(
        endpoints_df=endpoints_df,
        agnostic_metrics=agnostic_metrics,
        aware_metrics=aware_metrics,
        summary_stat=summary_stat,
    )
    matrix_a_b = _pair_group_heatmap_matrix(sample_summary=sample_summary, summary_stat=summary_stat)
    # Panel A shares Panel B's denoising/clustering categories. Apply any
    # Panel B label overrides here while retaining A_a's UPARSE/TACTIC labels.
    panel_a_a_group_labels = dict(PANEL_A_ENDPOINT_GROUP_LABELS)
    for group in ("denoising_phylogeny_aware", "denoising_phylogeny_agnostic"):
        panel_a_a_group_labels[group] = panel_b_group_labels.get(group, panel_a_a_group_labels[group])
    panel_a_a_rename = {
        PANEL_A_ENDPOINT_GROUP_LABELS[group]: label
        for group, label in panel_a_a_group_labels.items()
    }
    matrix_a_a = matrix_a_a.rename(index=panel_a_a_rename, columns=panel_a_a_rename)

    panel_a_b_rename = {
        ENDPOINT_GROUP_LABELS[group]: panel_b_group_labels.get(group, ENDPOINT_GROUP_LABELS[group])
        for group in ENDPOINT_GROUP_ORDER
    }
    matrix_a_b = matrix_a_b.rename(index=panel_a_b_rename, columns=panel_a_b_rename)
    # Backward-compatible path keeps the detailed panel A_a matrix.
    matrix_a_a.to_csv(out_dir / "pair_group_distance_contrast_heatmap_matrix.tsv", sep="\t")
    matrix_a_a.to_csv(out_dir / "pair_group_distance_contrast_heatmap_matrix_panel_a_a.tsv", sep="\t")
    matrix_a_b.to_csv(out_dir / "pair_group_distance_contrast_heatmap_matrix_panel_a_b.tsv", sep="\t")
    endpoint_matrix, endpoint_meta = _endpoint_distance_heatmap_matrix(
        endpoints_df=endpoints_df,
        method_order=method_order,
        metric_order=metric_order,
        agnostic_metrics=agnostic_metrics,
        aware_metrics=aware_metrics,
        summary_stat=summary_stat,
        order_mode="metric_family_first",
    )
    endpoint_matrix_method_first, endpoint_meta_method_first = _endpoint_distance_heatmap_matrix(
        endpoints_df=endpoints_df,
        method_order=method_order,
        metric_order=metric_order,
        agnostic_metrics=agnostic_metrics,
        aware_metrics=aware_metrics,
        summary_stat=summary_stat,
        order_mode="method_family_first",
    )
    if not endpoint_matrix.empty:
        endpoint_matrix.to_csv(out_dir / "pair_group_distance_contrast_endpoint_matrix.tsv", sep="\t")
    if not endpoint_meta.empty:
        endpoint_meta.to_csv(out_dir / "pair_group_distance_contrast_endpoint_label_map.tsv", sep="\t", index=False)
    if not endpoint_matrix_method_first.empty:
        endpoint_matrix_method_first.to_csv(
            out_dir / "pair_group_distance_contrast_endpoint_matrix_method_family_first.tsv",
            sep="\t",
        )
    if not endpoint_meta_method_first.empty:
        endpoint_meta_method_first.to_csv(
            out_dir / "pair_group_distance_contrast_endpoint_label_map_method_family_first.tsv",
            sep="\t",
            index=False,
        )

    # Combined figure layout. Panels E and F are optional bottom rows.
    height_ratios = [1.0, 1.70]
    fig_height = 17.6
    if show_panel_e_raincloud_grid:
        height_ratios.append(1.55)
        fig_height += 8.0
    if show_panel_f:
        height_ratios.append(0.95)
        fig_height += 5.2
    fig = plt.figure(figsize=(25.2, fig_height))
    gs = fig.add_gridspec(
        len(height_ratios),
        2,
        width_ratios=[1.0, 1.55],
        height_ratios=height_ratios,
    )
    panel_a_spec = gs[0, 0]
    ax_dist = fig.add_subplot(gs[0, 1])
    ax_endpoint = fig.add_subplot(gs[1, 0])
    ax_direction = fig.add_subplot(gs[1, 1])

    if show_panel_a:
        panel_a_inner = panel_a_spec.subgridspec(
            1,
            3,
            width_ratios=[1.0, 1.0, 0.055],
            wspace=0.28,
        )
        ax_heat_a = fig.add_subplot(panel_a_inner[0, 0])
        ax_heat_b = fig.add_subplot(panel_a_inner[0, 1])
        cax_heat = fig.add_subplot(panel_a_inner[0, 2])

        panel_a_a_tick_labels = [str(label) for label in matrix_a_a.columns]
        panel_a_b_tick_labels = [str(label) for label in matrix_a_b.columns]
        finite_a_a = matrix_a_a.to_numpy(dtype=float)
        finite_a_a = finite_a_a[np.isfinite(finite_a_a)]
        finite_a_b = matrix_a_b.to_numpy(dtype=float)
        finite_a_b = finite_a_b[np.isfinite(finite_a_b)]
        finite_all = np.concatenate([finite_a_a, finite_a_b]) if (finite_a_a.size or finite_a_b.size) else np.array([1.0])
        panel_a_vmax = float(np.max(finite_all)) if finite_all.size else 1.0
        if panel_a_vmax <= 0.0:
            panel_a_vmax = 1.0

        _draw_pair_group_heatmap_panel(
            ax_heat_a,
            sample_summary=sample_summary,
            font_size=panel_a_font_size,
            title=panel_a_a_title,
            matrix_override=matrix_a_a,
            short_tick_labels=panel_a_a_tick_labels,
            xtick_rotation=panel_a_xtick_rotation,
            ytick_rotation=panel_a_ytick_rotation,
            summary_stat=summary_stat,
            vmin=0.0,
            vmax_override=panel_a_vmax,
            x_axis_label=panel_a_x_axis_label,
            y_axis_label=panel_a_y_axis_label,
        )
        _draw_pair_group_heatmap_panel(
            ax_heat_b,
            sample_summary=sample_summary,
            font_size=panel_a_font_size,
            title=panel_a_b_title,
            matrix_override=matrix_a_b,
            short_tick_labels=panel_a_b_tick_labels,
            xtick_rotation=panel_a_xtick_rotation,
            ytick_rotation=panel_a_ytick_rotation,
            summary_stat=summary_stat,
            vmin=0.0,
            vmax_override=panel_a_vmax,
            x_axis_label=panel_a_x_axis_label,
            y_axis_label=panel_a_y_axis_label,
        )
        if panel_a_show_planned_comparison_annotations and planned_sig_rows:
            _draw_panel_a_b_planned_comparison_annotations(
                ax=ax_heat_b,
                planned_sig_rows=planned_sig_rows,
                alpha=stats_alpha,
                font_size=panel_a_font_size,
            )
        if ax_heat_a.images:
            cbar = fig.colorbar(ax_heat_a.images[0], cax=cax_heat)
            cbar.set_label(
                _summary_stat_colorbar_label(summary_stat=summary_stat, summarized=True),
                fontsize=fs_panel_a["cbar_label"],
            )
            cbar.ax.tick_params(labelsize=fs_panel_a["cbar_tick"])
    else:
        ax_panel_a_blank = fig.add_subplot(panel_a_spec)
        ax_panel_a_blank.axis("off")

    if show_panel_b:
        _draw_pair_group_distribution_panel(
            ax_dist,
            sample_summary=sample_summary_panel_b,
            font_size=panel_b_font_size,
            point_size=point_size,
            plot_type=plot_type,
            violin_bw_method=violin_bw_method,
            stats=stats,
            planned_sig_rows=planned_sig_rows,
            stats_enabled=stats_enabled,
            stats_show_pairwise=stats_show_pairwise,
            show_within_contrast_stats=panel_b_show_within_contrast_stats,
            stats_alpha=stats_alpha,
            title=panel_b_title,
            x_axis_label=panel_b_x_axis_label,
            y_axis_label=panel_b_y_axis_label,
            pair_groups_to_show=panel_b_pair_groups_to_show,
            endpoint_group_labels=panel_b_group_labels,
            sample_v4_multiple_label=sample_v4_multiple_label,
            sample_v4_single_label=sample_v4_single_label,
            xtick_rotation=panel_b_xtick_rotation,
            ytick_rotation=panel_b_ytick_rotation,
            dot_size=panel_b_dot_size,
            dot_opacity=panel_b_dot_opacity,
            dot_jitter=panel_b_dot_jitter,
        )
    else:
        ax_dist.axis("off")

    if show_panel_c:
        panel_c_spec = ax_endpoint.get_subplotspec()
        ax_endpoint.remove()
        _draw_endpoint_matrix_dual_panel(
            fig,
            panel_c_spec,
            endpoints_df=endpoints_df,
            method_order=panel_d_methods,
            metric_order=metric_order,
            agnostic_metrics=agnostic_metrics,
            aware_metrics=aware_metrics,
            gap_cells=endpoint_matrix_gap_cells,
            font_size=panel_c_font_size,
            title=panel_c_title,
            summary_stat=summary_stat,
            xtick_rotation=endpoint_matrix_xtick_rotation,
            ytick_rotation=endpoint_matrix_ytick_rotation,
            panel_c_a_title=panel_c_a_title,
            panel_c_b_title=panel_c_b_title,
            panel_c_a_xtick_rotation=panel_c_a_xtick_rotation,
            panel_c_a_ytick_rotation=panel_c_a_ytick_rotation,
            panel_c_b_xtick_rotation=panel_c_b_xtick_rotation,
            panel_c_b_ytick_rotation=panel_c_b_ytick_rotation,
            x_axis_label=panel_c_x_axis_label,
            y_axis_label=panel_c_y_axis_label,
        )
    else:
        ax_endpoint.axis("off")

    if show_panel_d:
        _draw_directional_sensitivity_panel(
            ax_direction,
            sens=sens_panel_d,
            method_order=method_order,
            title=panel_d_title,
            plot_type=plot_type,
            point_size=point_size,
            font_size=panel_d_font_size,
            violin_bw_method=violin_bw_method,
            stats=directional_stats,
            stats_enabled=stats_enabled,
            stats_show_pairwise=stats_show_pairwise,
            stats_alpha=stats_alpha,
            x_axis_label=panel_d_x_axis_label,
            y_axis_label=panel_d_y_axis_label,
            sample_v4_multiple_label=sample_v4_multiple_label,
            sample_v4_single_label=sample_v4_single_label,
            dot_size=panel_d_dot_size,
            dot_opacity=panel_d_dot_opacity,
            dot_jitter=panel_d_dot_jitter,
            show_stats_annotations=panel_d_stats_show_annotations,
        )
    else:
        ax_direction.axis("off")

    panel_e_summary = pd.DataFrame()
    next_bottom_row = 2
    if show_panel_e_raincloud_grid:
        ax_panel_e = fig.add_subplot(gs[next_bottom_row, :])
        panel_e_spec = ax_panel_e.get_subplotspec()
        ax_panel_e.remove()
        panel_e_summary = _draw_pair_group_raincloud_grid_panel(
            fig=fig,
            parent_spec=panel_e_spec,
            pairwise_df=pairwise_df,
            sample_summary=sample_summary_panel_b,
            font_size=panel_e_font_size,
            title=panel_e_title,
            sample_v4_multiple_label=panel_e_sample_v4_multiple_label,
            sample_v4_single_label=panel_e_sample_v4_single_label,
            abs_density_label=panel_e_abs_density_label,
            summary_density_label=panel_e_summary_density_label,
            show_summary_density=panel_e_show_summary_density,
            mean_label=panel_e_mean_label,
            grid_label_font_size=panel_e_grid_label_font_size,
            grid_tick_font_size=panel_e_grid_tick_font_size,
            grid_text_font_size=panel_e_grid_text_font_size,
            grid_legend_font_size=panel_e_grid_legend_font_size,
            dot_size=panel_e_dot_size,
            dot_opacity=panel_e_dot_opacity,
            dot_jitter=panel_e_dot_jitter,
            x_axis_label=panel_e_x_axis_label,
            y_axis_label=panel_e_y_axis_label,
        )
        if not panel_e_summary.empty:
            panel_e_summary.to_csv(
                out_dir / "pair_group_distance_contrast_panel_e_cell_summary.tsv",
                sep="\t",
                index=False,
            )
        next_bottom_row += 1

    if show_panel_f:
        ax_deflation = fig.add_subplot(gs[next_bottom_row, :])
        _draw_method_metric_family_distance_panel(
            ax_deflation,
            points=deflation_panel_f,
            methods=panel_f_methods,
            title=panel_f_title,
            plot_type=plot_type,
            point_size=point_size,
            font_size=panel_f_font_size,
            violin_bw_method=violin_bw_method,
            x_axis_label=panel_f_x_axis_label,
            y_axis_label=panel_f_y_axis_label,
            x_tick_rotation=panel_f_xtick_rotation,
            dot_size=panel_f_dot_size,
            dot_opacity=panel_f_dot_opacity,
            dot_jitter=panel_f_dot_jitter,
            summary_statistic=panel_f_summary_statistic,
            stats_rows=deflation_stats,
            stats_alpha=stats_alpha,
            show_stats_annotations=panel_f_stats_show_annotations,
            show_non_significant_annotations=panel_f_stats_show_non_significant_annotations,
        )

    # Combined figure uses panel-specific titles, so keep the overall title concise.
    fig.suptitle(title, fontsize=fs_fig["base"] + 2.0, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_dir / "fig21_pair_group_distance_contrast.svg", format="svg", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Standalone panel exports.
    if show_panel_a:
        fig_a = plt.figure(figsize=(16.4, 7.0))
        gs_a = fig_a.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 0.055], wspace=0.28)
        ax_a_a = fig_a.add_subplot(gs_a[0, 0])
        ax_a_b = fig_a.add_subplot(gs_a[0, 1])
        cax_a = fig_a.add_subplot(gs_a[0, 2])
        _draw_pair_group_heatmap_panel(
            ax_a_a,
            sample_summary=sample_summary,
            font_size=panel_a_font_size,
            title=panel_a_a_title,
            matrix_override=matrix_a_a,
            short_tick_labels=panel_a_a_tick_labels,
            xtick_rotation=panel_a_xtick_rotation,
            ytick_rotation=panel_a_ytick_rotation,
            summary_stat=summary_stat,
            x_axis_label=panel_a_x_axis_label,
            y_axis_label=panel_a_y_axis_label,
            vmin=0.0,
            vmax_override=panel_a_vmax,
        )
        _draw_pair_group_heatmap_panel(
            ax_a_b,
            sample_summary=sample_summary,
            font_size=panel_a_font_size,
            title=panel_a_b_title,
            matrix_override=matrix_a_b,
            short_tick_labels=panel_a_b_tick_labels,
            xtick_rotation=panel_a_xtick_rotation,
            ytick_rotation=panel_a_ytick_rotation,
            summary_stat=summary_stat,
            x_axis_label=panel_a_x_axis_label,
            y_axis_label=panel_a_y_axis_label,
            vmin=0.0,
            vmax_override=panel_a_vmax,
        )
        if panel_a_show_planned_comparison_annotations and planned_sig_rows:
            _draw_panel_a_b_planned_comparison_annotations(
                ax=ax_a_b,
                planned_sig_rows=planned_sig_rows,
                alpha=stats_alpha,
                font_size=panel_a_font_size,
            )
        if ax_a_a.images:
            cbar = fig_a.colorbar(ax_a_a.images[0], cax=cax_a)
            cbar.set_label(
                _summary_stat_colorbar_label(summary_stat=summary_stat, summarized=True),
                fontsize=fs_panel_a["cbar_label"],
            )
            cbar.ax.tick_params(labelsize=fs_panel_a["cbar_tick"])
        # Manual spacing avoids tight_layout warnings with explicit colorbar axes.
        fig_a.subplots_adjust(left=0.045, right=0.985, top=0.93, bottom=0.10, wspace=0.28)
        fig_a.savefig(out_dir / "fig21_pair_group_distance_contrast_panel_a.svg", format="svg", dpi=300, bbox_inches="tight")
        plt.close(fig_a)

        fig_a_a, ax_a_a_only = plt.subplots(figsize=(8.2, 7.0))
        _draw_pair_group_heatmap_panel(
            ax_a_a_only,
            sample_summary=sample_summary,
            font_size=panel_a_font_size,
            title=panel_a_a_title,
            matrix_override=matrix_a_a,
            short_tick_labels=panel_a_a_tick_labels,
            xtick_rotation=panel_a_xtick_rotation,
            ytick_rotation=panel_a_ytick_rotation,
            summary_stat=summary_stat,
        )
        if ax_a_a_only.images:
            cbar = fig_a_a.colorbar(ax_a_a_only.images[0], ax=ax_a_a_only, fraction=0.046, pad=0.04)
            cbar.set_label(
                _summary_stat_colorbar_label(summary_stat=summary_stat, summarized=True),
                fontsize=fs_panel_a["cbar_label"],
            )
            cbar.ax.tick_params(labelsize=fs_panel_a["cbar_tick"])
        fig_a_a.tight_layout()
        fig_a_a.savefig(
            out_dir / "fig21_pair_group_distance_contrast_panel_a_a.svg",
            format="svg",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close(fig_a_a)

        fig_a_b, ax_a_b_only = plt.subplots(figsize=(8.2, 7.0))
        _draw_pair_group_heatmap_panel(
            ax_a_b_only,
            sample_summary=sample_summary,
            font_size=panel_a_font_size,
            title=panel_a_b_title,
            matrix_override=matrix_a_b,
            short_tick_labels=panel_a_b_tick_labels,
            xtick_rotation=panel_a_xtick_rotation,
            ytick_rotation=panel_a_ytick_rotation,
            summary_stat=summary_stat,
        )
        if panel_a_show_planned_comparison_annotations and planned_sig_rows:
            _draw_panel_a_b_planned_comparison_annotations(
                ax=ax_a_b_only,
                planned_sig_rows=planned_sig_rows,
                alpha=stats_alpha,
                font_size=panel_a_font_size,
            )
        if ax_a_b_only.images:
            cbar = fig_a_b.colorbar(ax_a_b_only.images[0], ax=ax_a_b_only, fraction=0.046, pad=0.04)
            cbar.set_label(
                _summary_stat_colorbar_label(summary_stat=summary_stat, summarized=True),
                fontsize=fs_panel_a["cbar_label"],
            )
            cbar.ax.tick_params(labelsize=fs_panel_a["cbar_tick"])
        fig_a_b.tight_layout()
        fig_a_b.savefig(
            out_dir / "fig21_pair_group_distance_contrast_panel_a_b.svg",
            format="svg",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close(fig_a_b)

    if show_panel_b:
        fig_b, ax_b = plt.subplots(figsize=(11.6, 7.0))
        _draw_pair_group_distribution_panel(
            ax_b,
            sample_summary=sample_summary_panel_b,
            font_size=panel_b_font_size,
            point_size=point_size,
            plot_type=plot_type,
            violin_bw_method=violin_bw_method,
            stats=stats,
            planned_sig_rows=planned_sig_rows,
            stats_enabled=stats_enabled,
            stats_show_pairwise=stats_show_pairwise,
            show_within_contrast_stats=panel_b_show_within_contrast_stats,
            stats_alpha=stats_alpha,
            title=panel_b_title,
            x_axis_label=panel_b_x_axis_label,
            y_axis_label=panel_b_y_axis_label,
            pair_groups_to_show=panel_b_pair_groups_to_show,
            endpoint_group_labels=panel_b_group_labels,
            sample_v4_multiple_label=sample_v4_multiple_label,
            sample_v4_single_label=sample_v4_single_label,
            xtick_rotation=panel_b_xtick_rotation,
            ytick_rotation=panel_b_ytick_rotation,
            dot_size=panel_b_dot_size,
            dot_opacity=panel_b_dot_opacity,
            dot_jitter=panel_b_dot_jitter,
        )
        fig_b.tight_layout()
        fig_b.savefig(out_dir / "fig21_pair_group_distance_contrast_panel_b.svg", format="svg", dpi=300, bbox_inches="tight")
        plt.close(fig_b)

    if show_panel_c:
        fig_c = plt.figure(figsize=(17.5, 21.0))
        panel_c_cspec = fig_c.add_gridspec(1, 1)[0]
        _draw_endpoint_matrix_dual_panel(
            fig_c,
            panel_c_cspec,
            endpoints_df=endpoints_df,
            method_order=method_order,
            metric_order=metric_order,
            agnostic_metrics=agnostic_metrics,
            aware_metrics=aware_metrics,
            gap_cells=endpoint_matrix_gap_cells,
            font_size=panel_c_font_size,
            title=panel_c_title,
            summary_stat=summary_stat,
            xtick_rotation=endpoint_matrix_xtick_rotation,
            ytick_rotation=endpoint_matrix_ytick_rotation,
            panel_c_a_title=panel_c_a_title,
            panel_c_b_title=panel_c_b_title,
            panel_c_a_xtick_rotation=panel_c_a_xtick_rotation,
            panel_c_a_ytick_rotation=panel_c_a_ytick_rotation,
            panel_c_b_xtick_rotation=panel_c_b_xtick_rotation,
            panel_c_b_ytick_rotation=panel_c_b_ytick_rotation,
            x_axis_label=panel_c_x_axis_label,
            y_axis_label=panel_c_y_axis_label,
        )
        fig_c.tight_layout()
        fig_c.savefig(out_dir / "fig21_pair_group_distance_contrast_panel_c.svg", format="svg", dpi=300, bbox_inches="tight")
        plt.close(fig_c)

        # Standalone Panel C subpanels with independent colorbar legends.
        fig_c_a, ax_c_a = plt.subplots(figsize=(14.2, 12.0))
        _draw_endpoint_matrix_panel(
            ax_c_a,
            endpoints_df=endpoints_df,
            method_order=method_order,
            metric_order=metric_order,
            agnostic_metrics=agnostic_metrics,
            aware_metrics=aware_metrics,
            gap_cells=endpoint_matrix_gap_cells,
            font_size=panel_c_font_size,
            title=panel_c_a_title,
            summary_stat=summary_stat,
            order_mode="metric_family_first",
            show_xticklabels=True,
            xtick_rotation=float(
                endpoint_matrix_xtick_rotation if panel_c_a_xtick_rotation is None else panel_c_a_xtick_rotation
            ),
            ytick_rotation=float(
                endpoint_matrix_ytick_rotation if panel_c_a_ytick_rotation is None else panel_c_a_ytick_rotation
            ),
            x_axis_label=panel_c_x_axis_label,
            y_axis_label=panel_c_y_axis_label,
        )
        if ax_c_a.images:
            cbar = fig_c_a.colorbar(ax_c_a.images[0], ax=ax_c_a, fraction=0.046, pad=0.04)
            cbar.set_label(
                _summary_stat_colorbar_label(summary_stat=summary_stat, summarized=False),
                fontsize=fs_panel_c["cbar_label"],
            )
            cbar.ax.tick_params(labelsize=fs_panel_c["cbar_tick"])
        fig_c_a.tight_layout()
        fig_c_a.savefig(out_dir / "fig21_pair_group_distance_contrast_panel_c_a.svg", format="svg", dpi=300, bbox_inches="tight")
        plt.close(fig_c_a)

        fig_c_b, ax_c_b = plt.subplots(figsize=(14.2, 12.0))
        _draw_endpoint_matrix_panel(
            ax_c_b,
            endpoints_df=endpoints_df,
            method_order=method_order,
            metric_order=metric_order,
            agnostic_metrics=agnostic_metrics,
            aware_metrics=aware_metrics,
            gap_cells=endpoint_matrix_gap_cells,
            font_size=panel_c_font_size,
            title=panel_c_b_title,
            summary_stat=summary_stat,
            order_mode="method_family_first",
            show_xticklabels=True,
            xtick_rotation=float(
                endpoint_matrix_xtick_rotation if panel_c_b_xtick_rotation is None else panel_c_b_xtick_rotation
            ),
            ytick_rotation=float(
                endpoint_matrix_ytick_rotation if panel_c_b_ytick_rotation is None else panel_c_b_ytick_rotation
            ),
            x_axis_label=panel_c_x_axis_label,
            y_axis_label=panel_c_y_axis_label,
        )
        if ax_c_b.images:
            cbar = fig_c_b.colorbar(ax_c_b.images[0], ax=ax_c_b, fraction=0.046, pad=0.04)
            cbar.set_label(
                _summary_stat_colorbar_label(summary_stat=summary_stat, summarized=False),
                fontsize=fs_panel_c["cbar_label"],
            )
            cbar.ax.tick_params(labelsize=fs_panel_c["cbar_tick"])
        fig_c_b.tight_layout()
        fig_c_b.savefig(out_dir / "fig21_pair_group_distance_contrast_panel_c_b.svg", format="svg", dpi=300, bbox_inches="tight")
        plt.close(fig_c_b)

    if show_panel_d:
        fig_d, ax_d = plt.subplots(figsize=(11.6, 7.0))
        _draw_directional_sensitivity_panel(
            ax_d,
            sens=sens_panel_d,
            method_order=panel_d_methods,
            title=panel_d_title,
            plot_type=plot_type,
            point_size=point_size,
            font_size=panel_d_font_size,
            violin_bw_method=violin_bw_method,
            stats=directional_stats,
            stats_enabled=stats_enabled,
            stats_show_pairwise=stats_show_pairwise,
            stats_alpha=stats_alpha,
            x_axis_label=panel_d_x_axis_label,
            y_axis_label=panel_d_y_axis_label,
            sample_v4_multiple_label=sample_v4_multiple_label,
            sample_v4_single_label=sample_v4_single_label,
            dot_size=panel_d_dot_size,
            dot_opacity=panel_d_dot_opacity,
            dot_jitter=panel_d_dot_jitter,
            show_stats_annotations=panel_d_stats_show_annotations,
        )
        fig_d.tight_layout()
        fig_d.savefig(out_dir / "fig21_pair_group_distance_contrast_panel_d.svg", format="svg", dpi=300, bbox_inches="tight")
        plt.close(fig_d)

    if show_panel_e_raincloud_grid:
        fig_e = plt.figure(figsize=(18.0, 14.5))
        panel_e_spec_single = fig_e.add_gridspec(1, 1)[0]
        panel_e_summary_single = _draw_pair_group_raincloud_grid_panel(
            fig=fig_e,
            parent_spec=panel_e_spec_single,
            pairwise_df=pairwise_df,
            sample_summary=sample_summary_panel_b,
            font_size=panel_e_font_size,
            title=panel_e_title,
            sample_v4_multiple_label=panel_e_sample_v4_multiple_label,
            sample_v4_single_label=panel_e_sample_v4_single_label,
            abs_density_label=panel_e_abs_density_label,
            summary_density_label=panel_e_summary_density_label,
            show_summary_density=panel_e_show_summary_density,
            mean_label=panel_e_mean_label,
            grid_label_font_size=panel_e_grid_label_font_size,
            grid_tick_font_size=panel_e_grid_tick_font_size,
            grid_text_font_size=panel_e_grid_text_font_size,
            grid_legend_font_size=panel_e_grid_legend_font_size,
            dot_size=panel_e_dot_size,
            dot_opacity=panel_e_dot_opacity,
            dot_jitter=panel_e_dot_jitter,
            x_axis_label=panel_e_x_axis_label,
            y_axis_label=panel_e_y_axis_label,
        )
        if not panel_e_summary_single.empty:
            panel_e_summary_single.to_csv(
                out_dir / "pair_group_distance_contrast_panel_e_cell_summary.tsv",
                sep="\t",
                index=False,
            )
        fig_e.tight_layout()
        fig_e.savefig(
            out_dir / "fig21_pair_group_distance_contrast_panel_e.svg",
            format="svg",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close(fig_e)

    if show_panel_f:
        fig_f, ax_f = plt.subplots(figsize=(11.6, 7.0))
        _draw_method_metric_family_distance_panel(
            ax_f,
            points=deflation_panel_f,
            methods=panel_f_methods,
            title=panel_f_title,
            plot_type=plot_type,
            point_size=point_size,
            font_size=panel_f_font_size,
            violin_bw_method=violin_bw_method,
            x_axis_label=panel_f_x_axis_label,
            y_axis_label=panel_f_y_axis_label,
            x_tick_rotation=panel_f_xtick_rotation,
            dot_size=panel_f_dot_size,
            dot_opacity=panel_f_dot_opacity,
            dot_jitter=panel_f_dot_jitter,
            summary_statistic=panel_f_summary_statistic,
            stats_rows=deflation_stats,
            stats_alpha=stats_alpha,
            show_stats_annotations=panel_f_stats_show_annotations,
            show_non_significant_annotations=panel_f_stats_show_non_significant_annotations,
        )
        fig_f.tight_layout()
        fig_f.savefig(
            out_dir / "fig21_pair_group_distance_contrast_panel_f.svg",
            format="svg",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close(fig_f)

def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf-8")


def _write_fig21_readmes(root_dir: Path) -> None:
    _write_text(
        root_dir / "README.md",
        """# Figure 21 Metric-Family Sensitivity

This folder contains the Figure 21 diagnostics for comparing phylogeny-agnostic and phylogeny-aware distance behavior across methods.

The outputs are:

- `pair_group_distance_contrast`: the main Figure 21 test of endpoint-pair groups.
- `fig21_pair_group_distance_contrast_panel_a.svg`: standalone combined Panel A block (A_a + A_b).
- `fig21_pair_group_distance_contrast_panel_a_a.svg`: standalone Panel A_a (detailed 6x6 grouping).
- `fig21_pair_group_distance_contrast_panel_a_b.svg`: standalone Panel A_b (simplified denoising/clustering grouping).
- `fig21_pair_group_distance_contrast_panel_b.svg`: standalone Panel B pair-group distribution plot. Point shape encodes whether the sample's expected references include multiple unique V4 sequences or only a single V4 reference.
- `fig21_pair_group_distance_contrast_panel_c.svg`: standalone Panel C endpoint-level matrix stack with the top heatmap ordered with all phylogeny-agnostic combinations first and the bottom heatmap ordered with all denoising combinations first, followed by clustering combinations.
- `fig21_pair_group_distance_contrast_panel_c_a.svg`: standalone Panel C_a (agnostic-first ordering) with its own colorbar legend.
- `fig21_pair_group_distance_contrast_panel_c_b.svg`: standalone Panel C_b (denoising-first ordering) with its own colorbar legend.
- `fig21_pair_group_distance_contrast_panel_d.svg`: standalone Panel D directional sensitivity plot. Point shape encodes whether a sample has multiple unique V4 references or only single-V4 references.
- `fig21_pair_group_distance_contrast_panel_e.svg`: standalone Panel E KDE + horizontal boxplot grid across pair-group cells.
- `fig21_pair_group_distance_contrast_panel_f.svg`: standalone Panel F method-level distance plot.
- `sensitivity_by_method`: per-sample sensitivity summaries by method.
- `sample_method_heatmap`: sample x method view of the same sensitivity values.
- `per_sample_method_metric_pcoa`: one PCoA per sample over method x metric distances.
- `method_metric_profile_pcoa`: one PCoA across samples using the full method x metric distance profile.
- `pairwise_difference_pcoa`: one PCoA across samples using pairwise differences between method x metric distances.

The main pair-group contrast starts from every sample's method-metric distances to expected composition.
It computes all pairwise differences among those method-metric endpoints, groups those contrasts by endpoint family, then summarizes each pair-group once per sample before statistical testing.

The older sensitivity diagnostics use:

`mean(agnostic metrics) - mean(aware metrics)`

where the metric families are controlled by the Figure 21 config block.
""",
    )
    _write_text(
        root_dir / "pair_group_distance_contrast" / "README.md",
        """# Pair-Group Distance Contrast

This directory contains the main Figure 21 implementation.

What it shows:

- Each sample starts with one distance-to-expected value for every method x metric endpoint.
- Those endpoints are assigned to four groups: denoising phylogeny-aware, denoising phylogeny-agnostic, clustering phylogeny-aware, and clustering phylogeny-agnostic.
- All endpoint pairs are compared within each sample, producing the full pairwise-difference table.
- For testing and plotting, the raw pairwise differences are summarized once per sample per endpoint-pair group.
- Panel A additionally uses a 6-group split for its heatmap: denoising, UPARSE, and TACTIC each crossed with phylogeny-aware vs phylogeny-agnostic.
- Panel A is split into two subpanels:
  - A_a: detailed 6-group split (denoising, UPARSE, and TACTIC crossed with phylogeny-aware vs phylogeny-agnostic)
  - A_b: simplified 4-group split (denoising vs clustering crossed with phylogeny-aware vs phylogeny-agnostic)

How it is calculated:

For each sample, Figure 21 computes:

`abs_difference = |distance(endpoint_i) - distance(endpoint_j)|`

for every pair of endpoints. The endpoint pair is then assigned to one of the configured pair groups, such as denoising-aware vs denoising-agnostic or clustering-aware vs clustering-agnostic.

The statistical test uses the sample-level summaries, not the raw endpoint pairs. This avoids treating the many within-sample pairwise differences as independent observations.
If `fig21.panel_b.show_within_contrast_stats` is enabled, Panel B also annotates significant paired tests among the four within-group contrast distributions.

In Panel B, point shape encodes sample V4 status:

- circle: multiple unique V4 reference sequences for at least one species in the sample's expected composition
- square: only single-V4 reference species

The two legend labels can be overridden with the Figure 21 config keys
`panel_b_sample_v4_multiple_label` and `panel_b_sample_v4_single_label`.

Panel C now shows two stacked endpoint-level matrices:

- top: all phylogeny-agnostic endpoint combinations first, then all phylogeny-aware combinations
- bottom: all denoising combinations first, then all clustering combinations

Panel D uses the directional contrast:

`mean(phylogeny-agnostic distances) - mean(phylogeny-aware distances)`

for each sample and method. Positive values mean the method appears farther from expected composition under phylogeny-agnostic metrics than under phylogeny-aware metrics. Panel D shows all methods side by side, uses pairwise tests between the analysis methods, and uses the same V4-status point-shape encoding as Panel B.

Panel F shows the raw mean distance to expected composition for each selected method within each metric family. Lines connect the same sample across the displayed methods, making method-level changes visible without collapsing methods into denoising or clustering families.

Important files:

- `fig21_pair_group_distance_contrast.svg`: four-panel plot arranged as 2x2 with
  A) split Panel A block:
     - A_a: 6x6 group-level matrix (Denoising / UPARSE / TACTIC x aware / agnostic),
     - A_b: 4x4 simplified matrix (Denoising / Clustering x aware / agnostic),
  B) pair-group violin/points,
  C) stacked endpoint-level method-metric matrices,
  D) directional method-level sensitivity.
  When Panel E is enabled, the combined figure also includes a bottom-row KDE + horizontal boxplot grid for pair-group cells, with per-sample dots overlaid in the boxplot strip and marker shape encoding sample V4 status.
- `fig21_pair_group_distance_contrast_panel_a_a.svg`: standalone Panel A_a detailed heatmap.
- `fig21_pair_group_distance_contrast_panel_a_b.svg`: standalone Panel A_b simplified heatmap.
- `fig21_pair_group_distance_contrast_panel_c_a.svg`: standalone Panel C_a (agnostic-first endpoint matrix) with its own colorbar legend.
- `fig21_pair_group_distance_contrast_panel_c_b.svg`: standalone Panel C_b (denoising-first endpoint matrix) with its own colorbar legend.
- `pair_group_distance_contrast_endpoint_values.tsv`: the method x metric endpoint distances used as inputs.
- `pair_group_distance_contrast_heatmap_matrix.tsv`: Panel A_a matrix (detailed 6-group version, backward-compatible filename).
- `pair_group_distance_contrast_heatmap_matrix_panel_a_a.tsv`: Panel A_a matrix (detailed 6-group version).
- `pair_group_distance_contrast_heatmap_matrix_panel_a_b.tsv`: Panel A_b matrix (simplified 4-group version).
- `pair_group_distance_contrast_endpoint_matrix.tsv`: endpoint-level absolute-difference matrix (mean or median, per config) for the agnostic-first ordering (typically 30x30 before inserting blank cells for Panel C).
- `pair_group_distance_contrast_endpoint_label_map.tsv`: endpoint ids and labels used in the agnostic-first ordering.
- `pair_group_distance_contrast_endpoint_matrix_method_family_first.tsv`: endpoint-level absolute-difference matrix (mean or median, per config) for the denoising-first ordering.
- `pair_group_distance_contrast_endpoint_label_map_method_family_first.tsv`: endpoint ids and labels used in the denoising-first ordering.
- `fig21_endpoint_matrix_gap_cells`: number of blank cells inserted between method blocks in Panel C.
- `fig21_endpoint_matrix_xtick_rotation`: x-axis label rotation in degrees for Panel C endpoint heatmaps.
- `fig21_endpoint_matrix_ytick_rotation`: y-axis label rotation in degrees for Panel C endpoint heatmaps.
- `fig21_panel_a_show_planned_comparison_annotations`: when true, Panel A_b overlays the planned-significant Panel B comparisons as connecting lines with `p` and `*` annotations.
- `pair_group_distance_contrast_pairwise_differences.tsv`: all raw pairwise endpoint differences.
- `pair_group_distance_contrast_sample_summary.tsv`: one summarized value per sample per pair group; this is the statistical input.
- `pair_group_distance_contrast_sample_summary_with_v4_status.tsv`: the same table with sample V4 status attached for Panel B shapes.
- `pair_group_distance_contrast_stats_planned.tsv`: the planned tests for the biological question.
- `pair_group_distance_contrast_stats_within.tsv`: the paired tests among the four within-group contrast distributions used by optional Panel B within-contrast annotations.
- `directional_sensitivity_method_points.tsv`: per-sample method-level directional sensitivities used for Panel D.
- `directional_sensitivity_method_stats_omnibus.tsv`: omnibus test across all analysis methods for Panel D.
- `directional_sensitivity_method_stats_pairwise.tsv`: pairwise method tests for Panel D.
- `distance_by_method_metric_family_points.tsv`: Panel F sample-level mean distances by method and metric family.
- `distance_by_method_metric_family_points_with_v4_status.tsv`: Panel F values with sample V4 status attached.
- `distance_by_method_metric_family_stats.tsv`: pairwise method tests within each metric family for Panel F.
""",
    )
    _write_text(
        root_dir / "sensitivity_by_method" / "README.md",
        """# Sensitivity By Method

This directory contains the main Figure 21 sensitivity summary for each method.

What it shows:

- Each point is one sample-level sensitivity value for one method.
- Sensitivity is plotted by method, usually as a violin or boxplot.
- The accompanying TSV files contain the raw points, summary statistics, and pairwise test results.

How it is calculated:

For each sample and method, Figure 21 computes:

`sensitivity = mean(distance over phylogeny-agnostic metrics) - mean(distance over phylogeny-aware metrics)`

The agnostic and aware metric groups are defined in the Figure 21 config.
""",
    )
    _write_text(
        root_dir / "sample_method_heatmap" / "README.md",
        """# Sample x Method Heatmap

This directory contains a heatmap view of Figure 21 sensitivity.

What it shows:

- Rows are samples.
- Columns are methods.
- Cell values are the per-sample sensitivity values for each method.

How it is calculated:

The heatmap uses the same sensitivity definition as `sensitivity_by_method`:

`mean(agnostic metrics) - mean(aware metrics)`

The values are reshaped into a matrix so you can compare samples and methods at a glance.
""",
    )
    _write_text(
        root_dir / "per_sample_method_metric_pcoa" / "README.md",
        """# Per-Sample Method x Metric PCoA

This directory contains one PCoA plot per sample.

What it shows:

- Each point is one `method x metric` distance-to-expected value for that sample.
- Point color encodes method.
- Point shape encodes metric.
- The plot is the PCoA of the pairwise absolute differences between those scalar distances.

How it is calculated:

For each sample, Figure 21 builds a vector of distances for all method-metric combinations.
If the sample has values like `d1, d2, ..., dn`, then the point-to-point distance matrix is:

`D(i, j) = |d_i - d_j|`

That distance matrix is then passed to PCoA.

Because each point starts as a single scalar distance, this construction is often effectively one-dimensional.
The `svg/` folder contains the rendered figures, and `distance_matrices/` contains the corresponding pairwise difference tables.
""",
    )
    _write_text(
        root_dir / "per_sample_method_metric_pcoa" / "distance_matrices" / "README.md",
        """# Distance Matrices

This folder stores the per-sample pairwise difference matrices used by the `per_sample_method_metric_pcoa` plots.

Each TSV file corresponds to one sample and one matrix over `method x metric` points.

The matrix entries are absolute differences between scalar distance values:

`D(i, j) = |d_i - d_j|`

These matrices are the direct input to the PCoA in the parent directory.
""",
    )
    _write_text(
        root_dir / "per_sample_method_metric_pcoa" / "svg" / "README.md",
        """# SVG Figures

This folder stores the rendered per-sample Figure 21 PCoA SVGs.

Each SVG corresponds to one sample and visualizes the PCoA computed from that sample's `method x metric` distance matrix.

The plot uses:

- color for method
- shape for metric

The coordinates come from the PCoA of the absolute-difference matrix stored in the sibling `distance_matrices/` folder.
""",
    )
    _write_text(
        root_dir / "method_metric_profile_pcoa" / "README.md",
        """# Method x Metric Profile PCoA

This directory contains the across-sample PCoA of full method-metric distance profiles.

What it shows:

- Each point is one sample.
- The sample is represented by its full vector of `distance__METHOD__METRIC` values.
- Points are colored by `mean_distance`.

How it is calculated:

For each sample, Figure 21 collects all method-metric distances into one vector.
Those vectors are standardized feature-wise, then Euclidean distances are computed between samples.
That distance matrix is passed to PCoA.

This plot asks whether samples cluster by their overall distance profile across methods and metrics.
""",
    )
    _write_text(
        root_dir / "pairwise_difference_pcoa" / "README.md",
        """# Pairwise Difference PCoA

This directory contains the most expanded sample-level Figure 21 view.

What it shows:

- Each point is one sample.
- The sample is represented not by the raw method-metric distances, but by all pairwise differences between those distances.
- Points are colored by `mean_distance`.

How it is calculated:

For a sample with method-metric distance vector `d`, Figure 21 computes all pairwise contrasts:

`d_i - d_j` for every pair `(i, j)`

Those pairwise-difference features form a much larger profile vector, which is standardized and compared with Euclidean distance across samples.
The resulting distance matrix is then used for PCoA.
""",
    )


def _standardize_matrix(x: np.ndarray, enabled: bool) -> tuple[np.ndarray, pd.DataFrame]:
    if not enabled:
        return x, pd.DataFrame()
    mu = np.nanmean(x, axis=0)
    sd = np.nanstd(x, axis=0, ddof=0)
    sd_safe = np.where(sd > 0.0, sd, 1.0)
    z = (x - mu[None, :]) / sd_safe[None, :]
    diag = pd.DataFrame(
        {
            "feature": [],
            "mean_before_standardization": [],
            "sd_before_standardization": [],
        }
    )
    return z, diag


def _axis_diagnostics(
    coords: np.ndarray,
    eigvals: np.ndarray,
    coord_df: pd.DataFrame,
    value_cols: list[str],
) -> pd.DataFrame:
    eig_sum = float(np.sum(eigvals)) if eigvals.size else 0.0
    rows: list[dict] = []
    communities = coord_df["community_id"].astype(str).to_numpy() if "community_id" in coord_df.columns else np.array([])
    for i in range(coords.shape[1]):
        axis = f"PC{i + 1}"
        vals = coords[:, i]
        row = {
            "axis": axis,
            "eigenvalue": float(eigvals[i]) if i < eigvals.size else np.nan,
            "variance_explained": float((eigvals[i] / eig_sum) if eig_sum > 0.0 and i < eigvals.size else 0.0),
            "n_points": int(coords.shape[0]),
            "n_features": int(len(value_cols)),
        }
        if communities.size == vals.size:
            row["r2_community"] = _factor_r2(vals, communities)
        for col in ["mean_sensitivity", "mean_distance", "profile_range", "profile_sd"]:
            if col in coord_df.columns:
                x = pd.to_numeric(coord_df[col], errors="coerce").to_numpy(dtype=float)
                ok = np.isfinite(x) & np.isfinite(vals)
                row[f"corr_{col}"] = float(np.corrcoef(vals[ok], x[ok])[0, 1]) if int(np.sum(ok)) >= 3 else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def _write_pcoa_outputs(
    vectors: pd.DataFrame,
    value_cols: list[str],
    out_dir: Path,
    stem: str,
    title: str,
    cfg,
    font_size: float,
    color_col: str,
) -> None:
    if vectors.empty or len(value_cols) < 2:
        return

    meta_cols = [c for c in vectors.columns if c not in value_cols]
    clean = vectors.dropna(subset=value_cols).copy()
    clean.to_csv(out_dir / f"{stem}_vectors.tsv", sep="\t", index=False)
    if clean.shape[0] < 3:
        return

    feature_summary = []
    for col in value_cols:
        vals = pd.to_numeric(clean[col], errors="coerce").dropna().to_numpy(dtype=float)
        feature_summary.append(
            {
                "feature": col,
                "mean": float(np.mean(vals)) if vals.size else np.nan,
                "sd": float(np.std(vals, ddof=1)) if vals.size > 1 else 0.0,
                "min": float(np.min(vals)) if vals.size else np.nan,
                "max": float(np.max(vals)) if vals.size else np.nan,
                "mean_abs": float(np.mean(np.abs(vals))) if vals.size else np.nan,
            }
        )
    pd.DataFrame(feature_summary).to_csv(out_dir / f"{stem}_feature_summary.tsv", sep="\t", index=False)

    x_raw = clean[value_cols].to_numpy(dtype=float)
    standardize = bool(getattr(cfg.figures, "fig21_standardize_profiles", True))
    x = x_raw.copy()
    if standardize:
        mu = np.mean(x, axis=0)
        sd = np.std(x, axis=0, ddof=0)
        sd_safe = np.where(sd > 0.0, sd, 1.0)
        x = (x - mu[None, :]) / sd_safe[None, :]
        pd.DataFrame(
            {
                "feature": value_cols,
                "mean_before_standardization": mu,
                "sd_before_standardization": sd,
                "sd_used": sd_safe,
            }
        ).to_csv(out_dir / f"{stem}_standardization_diagnostics.tsv", sep="\t", index=False)

    profile_distance = str(getattr(cfg.figures, "fig21_profile_distance", "euclidean")).strip().lower()
    if profile_distance not in {"euclidean", "l2"}:
        profile_distance = "euclidean"
    dmat = _build_distance_matrix(x, metric="euclidean")
    coords, eigvals = _pcoa(dmat)
    if coords.shape[1] < 2:
        return

    coord_df = clean[meta_cols].copy()
    for i in range(coords.shape[1]):
        coord_df[f"PC{i + 1}"] = coords[:, i]

    eig_sum = float(np.sum(eigvals)) if eigvals.size else 0.0
    eig_df = pd.DataFrame(
        {
            "axis": [f"PC{i + 1}" for i in range(len(eigvals))],
            "eigenvalue": eigvals,
            "variance_explained": (eigvals / eig_sum) if eig_sum > 0.0 else np.zeros_like(eigvals),
        }
    )
    coord_df.to_csv(out_dir / f"{stem}_pcoa_coords.tsv", sep="\t", index=False)
    eig_df.to_csv(out_dir / f"{stem}_pcoa_eigenvalues.tsv", sep="\t", index=False)
    _axis_diagnostics(coords, eigvals, coord_df, value_cols).to_csv(
        out_dir / f"{stem}_pcoa_axis_diagnostics.tsv",
        sep="\t",
        index=False,
    )

    dist_df = pd.DataFrame(dmat, index=clean["point_id"].astype(str), columns=clean["point_id"].astype(str))
    dist_df.to_csv(out_dir / f"{stem}_distance_matrix.tsv", sep="\t")

    ve = eigvals / eig_sum if eig_sum > 0.0 else np.zeros_like(eigvals)
    pc1_pct = 100.0 * float(ve[0]) if ve.size >= 1 else 0.0
    pc2_pct = 100.0 * float(ve[1]) if ve.size >= 2 else 0.0

    point_size = max(8.0, float(getattr(cfg.figures, "fig21_point_size", 42.0)))
    color_values = pd.to_numeric(coord_df[color_col], errors="coerce") if color_col in coord_df.columns else None
    fig, ax = plt.subplots(figsize=(10.6, 7.6))
    if color_values is not None and np.isfinite(color_values.to_numpy(dtype=float)).any():
        vals = color_values.to_numpy(dtype=float)
        scatter = ax.scatter(
            coord_df["PC1"],
            coord_df["PC2"],
            c=vals,
            cmap="viridis",
            s=point_size,
            alpha=0.82,
            edgecolors="#111111",
            linewidths=0.35,
            zorder=3,
        )
        cbar = fig.colorbar(scatter, ax=ax, fraction=0.035, pad=0.02)
        cbar.set_label(color_col.replace("_", " "), fontsize=max(7.0, font_size - 2.0))
        cbar.ax.tick_params(labelsize=max(6.0, font_size - 3.0))
    else:
        ax.scatter(
            coord_df["PC1"],
            coord_df["PC2"],
            color="#4D4D4D",
            s=point_size,
            alpha=0.82,
            edgecolors="#111111",
            linewidths=0.35,
            zorder=3,
        )

    if coord_df.shape[0] <= 40 and "sample_label" in coord_df.columns:
        for _, row in coord_df.iterrows():
            ax.text(
                float(row["PC1"]),
                float(row["PC2"]),
                str(row["sample_label"]),
                fontsize=max(5.0, font_size - 7.0),
                ha="left",
                va="bottom",
                alpha=0.75,
            )

    ax.axhline(0.0, color="#888888", linewidth=0.6, alpha=0.45)
    ax.axvline(0.0, color="#888888", linewidth=0.6, alpha=0.45)
    ax.grid(True, alpha=0.25, linewidth=0.6)
    ax.set_xlabel(f"PC1 ({pc1_pct:.1f}% var. explained)", fontsize=font_size)
    ax.set_ylabel(f"PC2 ({pc2_pct:.1f}% var. explained)", fontsize=font_size)
    ax.tick_params(axis="both", labelsize=max(6.0, font_size - 2.0))
    ax.set_title(title, fontsize=font_size + 1.0)
    fig.tight_layout()
    fig.savefig(out_dir / f"{stem}.svg", format="svg", dpi=300, bbox_inches="tight")
    plt.close(fig)


def _sensitivity_profile_vectors(sens: pd.DataFrame, method_order: list[str]) -> tuple[pd.DataFrame, list[str]]:
    if sens.empty:
        return pd.DataFrame(), []
    pivot = sens.pivot_table(
        index=["community_id", "sample_id", "sample_key", "community_label", "sample_label"],
        columns="method",
        values="sensitivity",
        aggfunc="first",
    ).reset_index()
    methods = [m for m in method_order if m in pivot.columns]
    value_cols = []
    for method in methods:
        col = f"sensitivity__{method}"
        pivot[col] = pd.to_numeric(pivot[method], errors="coerce")
        value_cols.append(col)
    pivot = pivot.drop(columns=methods)
    pivot["point_id"] = pivot["sample_key"].astype(str)
    arr = pivot[value_cols].to_numpy(dtype=float)
    pivot["mean_sensitivity"] = np.nanmean(arr, axis=1)
    pivot["profile_range"] = np.nanmax(arr, axis=1) - np.nanmin(arr, axis=1)
    pivot["profile_sd"] = np.nanstd(arr, axis=1, ddof=0)
    return pivot, value_cols


def _method_metric_profile_vectors(
    long: pd.DataFrame,
    method_order: list[str],
    metric_order: list[str],
    agnostic_metrics: list[str],
    aware_metrics: list[str],
) -> tuple[pd.DataFrame, list[str]]:
    sub = long[long["method"].isin(method_order) & long["metric"].isin(metric_order)].copy()
    if sub.empty:
        return pd.DataFrame(), []
    sub["method_metric"] = sub["method"].astype(str) + "__" + sub["metric"].astype(str)
    feature_order = [f"{m}__{metric}" for m in method_order for metric in metric_order]
    pivot = sub.pivot_table(
        index=["community_id", "sample_id", "sample_key"],
        columns="method_metric",
        values="distance",
        aggfunc="first",
    ).reindex(columns=feature_order)
    pivot = pivot.dropna(axis=0, how="any").reset_index()
    if pivot.empty:
        return pd.DataFrame(), []
    pivot["community_label"] = pivot["community_id"].map(display_community_id)
    pivot["sample_label"] = pivot["sample_key"].map(display_sample_key)
    pivot["point_id"] = pivot["sample_key"].astype(str)
    value_cols = []
    for feat in feature_order:
        if feat not in pivot.columns:
            continue
        col = f"distance__{feat}"
        pivot[col] = pd.to_numeric(pivot[feat], errors="coerce")
        value_cols.append(col)
    pivot = pivot.drop(columns=[f for f in feature_order if f in pivot.columns])
    arr = pivot[value_cols].to_numpy(dtype=float)
    pivot["mean_distance"] = np.nanmean(arr, axis=1)
    pivot["profile_range"] = np.nanmax(arr, axis=1) - np.nanmin(arr, axis=1)
    pivot["profile_sd"] = np.nanstd(arr, axis=1, ddof=0)

    ag_cols = [f"distance__{m}__{metric}" for m in method_order for metric in agnostic_metrics if f"distance__{m}__{metric}" in value_cols]
    aw_cols = [f"distance__{m}__{metric}" for m in method_order for metric in aware_metrics if f"distance__{m}__{metric}" in value_cols]
    if ag_cols:
        pivot["mean_agnostic_distance"] = np.nanmean(pivot[ag_cols].to_numpy(dtype=float), axis=1)
    if aw_cols:
        pivot["mean_aware_distance"] = np.nanmean(pivot[aw_cols].to_numpy(dtype=float), axis=1)
    if ag_cols and aw_cols:
        pivot["mean_sensitivity"] = pivot["mean_agnostic_distance"] - pivot["mean_aware_distance"]
    return pivot, value_cols


def _pairwise_difference_vectors(vectors: pd.DataFrame, value_cols: list[str]) -> tuple[pd.DataFrame, list[str]]:
    if vectors.empty or len(value_cols) < 2:
        return pd.DataFrame(), []
    meta_cols = [c for c in vectors.columns if c not in value_cols]
    diff_cols: list[str] = []
    diff_data: dict[str, pd.Series] = {}
    for a, b in combinations(value_cols, 2):
        name_a = a.replace("distance__", "").replace("sensitivity__", "")
        name_b = b.replace("distance__", "").replace("sensitivity__", "")
        col = f"diff__{name_a}__minus__{name_b}"
        diff_data[col] = pd.to_numeric(vectors[a], errors="coerce") - pd.to_numeric(vectors[b], errors="coerce")
        diff_cols.append(col)
    diff_df = pd.DataFrame(diff_data, index=vectors.index)
    out = pd.concat([vectors[meta_cols].copy(), diff_df], axis=1)
    arr = out[diff_cols].to_numpy(dtype=float)
    out["mean_pairwise_diff"] = np.nanmean(arr, axis=1)
    out["mean_abs_pairwise_diff"] = np.nanmean(np.abs(arr), axis=1)
    out["profile_range"] = np.nanmax(arr, axis=1) - np.nanmin(arr, axis=1)
    out["profile_sd"] = np.nanstd(arr, axis=1, ddof=0)
    if "mean_sensitivity" not in out.columns and "mean_agnostic_distance" in out.columns and "mean_aware_distance" in out.columns:
        out["mean_sensitivity"] = out["mean_agnostic_distance"] - out["mean_aware_distance"]
    return out, diff_cols


def run(cfg) -> int:
    """
    Figures 21-26: metric-family sensitivity diagnostics suite.

    Core contrast:
      mean(phylogeny-agnostic distances) - mean(phylogeny-aware distances)

    Outputs:
      - Figure 21: pair-group distance contrast
      - Figure 22: per-method sensitivity violin/boxplot + stats
      - Figure 23: sample x method sensitivity heatmap
      - Figure 24: one PCoA per sample where each point is a method x metric distance
      - Figure 25: PCoA of per-sample method x metric distance profiles
      - Figure 26: PCoA of per-sample pairwise differences among method x metric distances
    """
    log = get_logger()
    outdir = Path(cfg.paths.outdir)
    fr = figures_root(outdir)
    fig21_dir = fr / "fig21_pair_group_distance_contrast"
    fig22_dir = fr / "fig22_sensitivity_by_method"
    fig23_dir = fr / "fig23_sample_method_heatmap"
    fig24_dir = fr / "fig24_per_sample_method_metric_pcoa"
    fig25_dir = fr / "fig25_method_metric_profile_pcoa"
    fig26_dir = fr / "fig26_pairwise_difference_pcoa"

    long = _load_fig21_distances(cfg, log)
    if long.empty:
        log.info("Figure 21: no distance rows to plot.")
        return 0
    required = {"community_id", "sample_id", "method", "metric", "distance"}
    if not required.issubset(set(long.columns)):
        log.info("Figure 21: distances table is missing required columns.")
        return 0

    long["community_id"] = long["community_id"].astype(str)
    long["sample_id"] = long["sample_id"].astype(str)
    long["method"] = long["method"].astype(str)
    long["metric"] = long["metric"].astype(str)
    long["distance"] = pd.to_numeric(long["distance"], errors="coerce")
    long = long[np.isfinite(long["distance"])].copy()
    long["sample_key"] = long["community_id"] + "|" + long["sample_id"]

    exclude_patterns = list(getattr(cfg.figures, "exclude_community_name_patterns", []) or [])
    exclude_patterns += list(getattr(cfg.figures, "exclude_community_name_patterns_fig21", []) or [])
    if exclude_patterns:
        long = long[~long["community_id"].map(lambda cid: _is_excluded(cid, exclude_patterns))].copy()
        if long.empty:
            log.info("Figure 21: no rows left after figure-level community exclusions.")
            return 0

    configured_methods = list(getattr(cfg.distances, "methods", METHOD_ORDER) or METHOD_ORDER)
    method_order = [m for m in METHOD_ORDER if m in configured_methods] + [m for m in configured_methods if m not in METHOD_ORDER]
    method_order = [m for m in method_order if m in set(long["method"].astype(str))]
    if not method_order:
        log.info("Figure 21: no configured methods found in distances.")
        return 0
    long = long[long["method"].isin(method_order)].copy()

    metrics_cfg = [str(m) for m in (getattr(cfg.figures, "fig21_metrics", []) or []) if str(m).strip()]
    configured_metrics = metrics_cfg or [str(m) for m in (getattr(cfg.distances, "metrics", []) or [])]
    metric_order = [m for m in configured_metrics if m in set(long["metric"].astype(str)) and m in METRIC_LABEL]
    if not metric_order:
        log.info("Figure 21: no configured non-Aitchison metrics found in distances.")
        return 0
    long = long[long["metric"].isin(metric_order)].copy()

    agnostic_metrics, aware_metrics = _metric_family_lists(cfg, metric_order)
    sens = _metric_family_sensitivity_points(
        long=long,
        method_order=method_order,
        agnostic_metrics=agnostic_metrics,
        aware_metrics=aware_metrics,
    )
    if sens.empty:
        log.info("Figure 21: sensitivity table is empty; need at least one agnostic and one phylogeny-aware metric.")
        return 0

    panel_b_denoising_methods, panel_b_clustering_methods = _resolve_panel_b_method_groups(
        getattr(cfg.figures, "fig21_panel_b_method_groups", {})
    )
    panel_b_group_labels = dict(ENDPOINT_GROUP_LABELS)
    panel_b_group_labels.update(
        {
            str(group): str(label)
            for group, label in getattr(cfg.figures, "fig21_panel_b_group_labels", {}).items()
            if str(label).strip()
        }
    )
    panel_d_methods = _resolve_method_list(getattr(cfg.figures, "fig21_panel_d_methods", []), method_order)
    panel_f_methods = _resolve_method_list(getattr(cfg.figures, "fig21_panel_f_methods", []), method_order)
    panel_f_summary_statistic = _normalize_panel_f_summary_statistic(
        getattr(cfg.figures, "fig21_panel_f_summary_statistic", "median")
    )
    panel_f_stats_pairwise_test_raw = getattr(cfg.figures, "fig21_panel_f_stats_pairwise_test", None)
    panel_f_stats_alternative = _normalize_permutation_alternative(
        getattr(cfg.figures, "fig21_panel_f_stats_alternative", "two-sided")
    )
    panel_f_stats_practical_difference_delta = max(
        0.0, float(getattr(cfg.figures, "fig21_panel_f_stats_practical_difference_delta", 0.0))
    )
    panel_f_stats_bootstrap_resamples = max(
        1, int(getattr(cfg.figures, "fig21_panel_f_stats_bootstrap_resamples", 4000))
    )
    panel_f_stats_ci_level = float(getattr(cfg.figures, "fig21_panel_f_stats_ci_level", 0.95))
    panel_f_stats_show_annotations = bool(
        getattr(cfg.figures, "fig21_panel_f_stats_show_annotations", True)
    )
    panel_f_stats_show_non_significant_annotations = bool(
        getattr(cfg.figures, "fig21_panel_f_stats_show_non_significant_annotations", True)
    )
    panel_d_stats_show_annotations = bool(
        getattr(cfg.figures, "fig21_panel_d_stats_show_annotations", True)
    )
    panel_f_stats_planned_comparisons = list(
        getattr(cfg.figures, "fig21_panel_f_stats_planned_comparisons", []) or []
    )

    metric_family_def = pd.DataFrame(
        {
            "metric_family": ["phylogeny_agnostic", "phylogeny_aware"],
            "metrics": [",".join(agnostic_metrics), ",".join(aware_metrics)],
            "n_metrics": [len(agnostic_metrics), len(aware_metrics)],
        }
    )

    font_size = max(6.0, float(getattr(cfg.figures, "fig21_font_size", 12.0)))
    panel_a_title = _resolve_fig21_panel_title(
        getattr(cfg.figures, "fig21_panel_a_title", None),
        "A. Group-level contrast",
    )
    panel_a_a_title = _resolve_fig21_panel_title(
        getattr(cfg.figures, "fig21_panel_a_a_title", None),
        f"{panel_a_title} (detailed)",
    )
    panel_a_b_title = _resolve_fig21_panel_title(
        getattr(cfg.figures, "fig21_panel_a_b_title", None),
        f"{panel_a_title} (simplified)",
    )
    panel_b_title = _resolve_fig21_panel_title(
        getattr(cfg.figures, "fig21_panel_b_title", None),
        "B. Pair-group distributions across samples",
    )
    panel_c_title = _resolve_fig21_panel_title(
        getattr(cfg.figures, "fig21_panel_c_title", None),
        "C. Endpoint-level median contrast",
    )
    panel_c_a_title = _resolve_fig21_panel_title(
        getattr(cfg.figures, "fig21_panel_c_a_title", None),
        f"{panel_c_title}\nAgnostic-first ordering",
    )
    panel_c_b_title = _resolve_fig21_panel_title(
        getattr(cfg.figures, "fig21_panel_c_b_title", None),
        "Denoising-first ordering",
    )
    panel_d_title = _resolve_fig21_panel_title(
        getattr(cfg.figures, "fig21_panel_d_title", None),
        "D. Directional method-level sensitivity",
    )
    panel_e_title = _resolve_fig21_panel_title(
        getattr(cfg.figures, "fig21_panel_e_title", None),
        "E. Pair-group KDE + boxplot grid of absolute differences",
    )
    panel_f_title = _resolve_fig21_panel_title(
        getattr(cfg.figures, "fig21_panel_f_title", None),
        "F. Denoising vs clustering distance to expected",
    )
    panel_a_font_size = _resolve_fig21_panel_font_size(
        getattr(cfg.figures, "fig21_panel_a_font_size", None),
        font_size,
    )
    panel_b_font_size = _resolve_fig21_panel_font_size(
        getattr(cfg.figures, "fig21_panel_b_font_size", None),
        font_size,
    )
    panel_c_font_size = _resolve_fig21_panel_font_size(
        getattr(cfg.figures, "fig21_panel_c_font_size", None),
        font_size,
    )
    panel_d_font_size = _resolve_fig21_panel_font_size(
        getattr(cfg.figures, "fig21_panel_d_font_size", None),
        font_size,
    )
    panel_e_font_size = _resolve_fig21_panel_font_size(
        getattr(cfg.figures, "fig21_panel_e_font_size", None),
        font_size,
    )
    panel_f_font_size = _resolve_fig21_panel_font_size(
        getattr(cfg.figures, "fig21_panel_f_font_size", None),
        font_size,
    )
    panel_e_grid_label_font_size = getattr(cfg.figures, "fig21_panel_e_grid_label_font_size", None)
    panel_e_grid_tick_font_size = getattr(cfg.figures, "fig21_panel_e_grid_tick_font_size", None)
    panel_e_grid_text_font_size = getattr(cfg.figures, "fig21_panel_e_grid_text_font_size", None)
    panel_e_grid_legend_font_size = getattr(cfg.figures, "fig21_panel_e_grid_legend_font_size", None)
    panel_e_grid_label_font_size = (
        None if panel_e_grid_label_font_size is None else max(4.0, float(panel_e_grid_label_font_size))
    )
    panel_e_grid_tick_font_size = (
        None if panel_e_grid_tick_font_size is None else max(4.0, float(panel_e_grid_tick_font_size))
    )
    panel_e_grid_text_font_size = (
        None if panel_e_grid_text_font_size is None else max(4.0, float(panel_e_grid_text_font_size))
    )
    panel_e_grid_legend_font_size = (
        None if panel_e_grid_legend_font_size is None else max(4.0, float(panel_e_grid_legend_font_size))
    )
    sample_v4_multiple_label = str(
        getattr(cfg.figures, "fig21_panel_b_sample_v4_multiple_label", SAMPLE_STATUS_LABELS["multiple_v4_reference"])
        or SAMPLE_STATUS_LABELS["multiple_v4_reference"]
    )
    sample_v4_single_label = str(
        getattr(cfg.figures, "fig21_panel_b_sample_v4_single_label", SAMPLE_STATUS_LABELS["single_v4_reference"])
        or SAMPLE_STATUS_LABELS["single_v4_reference"]
    )
    panel_e_abs_density_label = str(
        getattr(cfg.figures, "fig21_panel_e_abs_density_label", "Pairwise abs_difference density")
        or "Pairwise abs_difference density"
    )
    panel_e_summary_density_label = str(
        getattr(cfg.figures, "fig21_panel_e_summary_density_label", "Sample summary_abs_difference density")
        or "Sample summary_abs_difference density"
    )
    panel_e_show_summary_density = bool(
        getattr(cfg.figures, "fig21_panel_e_show_summary_density", True)
    )
    panel_e_mean_label = str(
        getattr(cfg.figures, "fig21_panel_e_mean_label", "Mean")
        or "Mean"
    )
    panel_e_sample_v4_multiple_label = str(
        getattr(cfg.figures, "fig21_panel_e_sample_v4_multiple_label", None)
        or sample_v4_multiple_label
    )
    panel_e_sample_v4_single_label = str(
        getattr(cfg.figures, "fig21_panel_e_sample_v4_single_label", None)
        or sample_v4_single_label
    )
    panel_a_x_axis_label = getattr(cfg.figures, "fig21_panel_a_x_axis_label", None)
    panel_a_y_axis_label = getattr(cfg.figures, "fig21_panel_a_y_axis_label", None)
    panel_b_x_axis_label = str(getattr(cfg.figures, "fig21_panel_b_x_axis_label", "Endpoint-pair group") or "Endpoint-pair group")
    panel_b_y_axis_label = str(
        getattr(cfg.figures, "fig21_panel_b_y_axis_label", "Sample-level summarized absolute pairwise difference")
        or "Sample-level summarized absolute pairwise difference"
    )
    panel_c_x_axis_label = getattr(cfg.figures, "fig21_panel_c_x_axis_label", None)
    panel_c_y_axis_label = getattr(cfg.figures, "fig21_panel_c_y_axis_label", None)
    panel_d_x_axis_label = str(getattr(cfg.figures, "fig21_panel_d_x_axis_label", "Analysis method") or "Analysis method")
    panel_d_y_axis_label = str(
        getattr(cfg.figures, "fig21_panel_d_y_axis_label", "Mean agnostic distance - mean aware distance")
        or "Mean agnostic distance - mean aware distance"
    )
    panel_e_x_axis_label = getattr(cfg.figures, "fig21_panel_e_x_axis_label", None)
    panel_e_y_axis_label = getattr(cfg.figures, "fig21_panel_e_y_axis_label", None)
    panel_f_x_axis_label = str(
        getattr(cfg.figures, "fig21_panel_f_x_axis_label", "Metric family | method family")
        or "Metric family | method family"
    )
    panel_f_y_axis_label = str(
        getattr(
            cfg.figures,
            "fig21_panel_f_y_axis_label",
            "Mean distance to expected composition",
        )
        or "Mean distance to expected composition"
    )
    panel_f_xtick_rotation = float(getattr(cfg.figures, "fig21_panel_f_xtick_rotation", 0.0))
    panel_b_pair_groups_to_show = list(getattr(cfg.figures, "fig21_panel_b_pair_groups_to_show", []) or [])
    if not panel_b_pair_groups_to_show:
        panel_b_pair_groups_to_show = [
            "denoising_phylogeny_aware__vs__denoising_phylogeny_aware",
            "denoising_phylogeny_agnostic__vs__denoising_phylogeny_agnostic",
            "clustering_phylogeny_aware__vs__clustering_phylogeny_aware",
            "clustering_phylogeny_agnostic__vs__clustering_phylogeny_agnostic",
            "denoising_phylogeny_aware__vs__denoising_phylogeny_agnostic",
            "clustering_phylogeny_aware__vs__clustering_phylogeny_agnostic",
            "denoising_phylogeny_aware__vs__clustering_phylogeny_aware",
            "denoising_phylogeny_agnostic__vs__clustering_phylogeny_agnostic",
        ]
    panel_b_show_within_contrast_stats = bool(
        getattr(cfg.figures, "fig21_panel_b_show_within_contrast_stats", False)
    )
    endpoint_matrix_gap_cells = max(0, int(getattr(cfg.figures, "fig21_endpoint_matrix_gap_cells", 1)))
    panel_a_xtick_rotation = float(getattr(cfg.figures, "fig21_panel_a_xtick_rotation", 35.0))
    panel_a_ytick_rotation = float(getattr(cfg.figures, "fig21_panel_a_ytick_rotation", 0.0))
    panel_a_show_planned_comparison_annotations = bool(
        getattr(cfg.figures, "fig21_panel_a_show_planned_comparison_annotations", False)
    )
    show_panel_a = bool(getattr(cfg.figures, "fig21_show_panel_a", True))
    show_panel_b = bool(getattr(cfg.figures, "fig21_show_panel_b", True))
    show_panel_c = bool(getattr(cfg.figures, "fig21_show_panel_c", True))
    show_panel_d = bool(getattr(cfg.figures, "fig21_show_panel_d", True))
    show_panel_e = bool(getattr(cfg.figures, "fig21_show_panel_e", True))
    show_panel_f = bool(getattr(cfg.figures, "fig21_show_panel_f", True))
    show_panel_e_raincloud_grid = bool(getattr(cfg.figures, "fig21_show_panel_e_raincloud_grid", True)) and show_panel_e
    panel_b_xtick_rotation = float(getattr(cfg.figures, "fig21_panel_b_xtick_rotation", 35.0))
    panel_b_ytick_rotation = float(getattr(cfg.figures, "fig21_panel_b_ytick_rotation", 0.0))
    endpoint_matrix_xtick_rotation = float(getattr(cfg.figures, "fig21_endpoint_matrix_xtick_rotation", 90.0))
    endpoint_matrix_ytick_rotation = float(getattr(cfg.figures, "fig21_endpoint_matrix_ytick_rotation", 0.0))
    panel_c_a_xtick_rotation = getattr(cfg.figures, "fig21_panel_c_a_xtick_rotation", None)
    panel_c_a_ytick_rotation = getattr(cfg.figures, "fig21_panel_c_a_ytick_rotation", None)
    panel_c_b_xtick_rotation = getattr(cfg.figures, "fig21_panel_c_b_xtick_rotation", None)
    panel_c_b_ytick_rotation = getattr(cfg.figures, "fig21_panel_c_b_ytick_rotation", None)
    panel_c_a_xtick_rotation = None if panel_c_a_xtick_rotation is None else float(panel_c_a_xtick_rotation)
    panel_c_a_ytick_rotation = None if panel_c_a_ytick_rotation is None else float(panel_c_a_ytick_rotation)
    panel_c_b_xtick_rotation = None if panel_c_b_xtick_rotation is None else float(panel_c_b_xtick_rotation)
    panel_c_b_ytick_rotation = None if panel_c_b_ytick_rotation is None else float(panel_c_b_ytick_rotation)
    fig21_title = str(getattr(cfg.figures, "fig21_title", None) or "Figure 21: Metric-Family Sensitivity Diagnostics")
    fig22_title = str(getattr(cfg.figures, "fig22_title", None) or "Figure 22: Sensitivity By Method")
    fig23_title = str(getattr(cfg.figures, "fig23_title", None) or "Figure 23: Sample x Method Sensitivity Heatmap")
    fig24_title = str(
        getattr(cfg.figures, "fig24_title", None) or "Figure 24: Per-Sample Method x Metric PCoA"
    )
    fig25_title = str(
        getattr(cfg.figures, "fig25_title", None) or "Figure 25: Method x Metric Distance Profiles PCoA"
    )
    fig26_title = str(
        getattr(cfg.figures, "fig26_title", None) or "Figure 26: Pairwise Method x Metric Differences PCoA"
    )
    fig22_font_size = _resolve_fig21_panel_font_size(getattr(cfg.figures, "fig22_font_size", None), font_size)
    fig23_font_size = _resolve_fig21_panel_font_size(getattr(cfg.figures, "fig23_font_size", None), font_size)
    fig24_font_size = _resolve_fig21_panel_font_size(getattr(cfg.figures, "fig24_font_size", None), font_size)
    fig25_font_size = _resolve_fig21_panel_font_size(getattr(cfg.figures, "fig25_font_size", None), font_size)
    fig26_font_size = _resolve_fig21_panel_font_size(getattr(cfg.figures, "fig26_font_size", None), font_size)
    point_size = max(1.0, float(getattr(cfg.figures, "fig21_point_size", 42.0)))
    panel_b_dot_size = getattr(cfg.figures, "fig21_panel_b_dot_size", None)
    panel_d_dot_size = getattr(cfg.figures, "fig21_panel_d_dot_size", None)
    panel_e_dot_size = getattr(cfg.figures, "fig21_panel_e_dot_size", None)
    panel_f_dot_size = getattr(cfg.figures, "fig21_panel_f_dot_size", None)
    panel_b_dot_size = None if panel_b_dot_size is None else max(1.0, float(panel_b_dot_size))
    panel_d_dot_size = None if panel_d_dot_size is None else max(1.0, float(panel_d_dot_size))
    panel_e_dot_size = None if panel_e_dot_size is None else max(1.0, float(panel_e_dot_size))
    panel_f_dot_size = None if panel_f_dot_size is None else max(1.0, float(panel_f_dot_size))
    panel_b_dot_opacity = max(0.0, min(1.0, float(getattr(cfg.figures, "fig21_panel_b_dot_opacity", 0.48))))
    panel_d_dot_opacity = max(0.0, min(1.0, float(getattr(cfg.figures, "fig21_panel_d_dot_opacity", 0.48))))
    panel_e_dot_opacity = max(0.0, min(1.0, float(getattr(cfg.figures, "fig21_panel_e_dot_opacity", 0.95))))
    panel_f_dot_opacity = max(0.0, min(1.0, float(getattr(cfg.figures, "fig21_panel_f_dot_opacity", 0.75))))
    panel_b_dot_jitter = max(0.0, float(getattr(cfg.figures, "fig21_panel_b_dot_jitter", 0.065)))
    panel_d_dot_jitter = max(0.0, float(getattr(cfg.figures, "fig21_panel_d_dot_jitter", 0.07)))
    panel_e_dot_jitter = max(0.0, float(getattr(cfg.figures, "fig21_panel_e_dot_jitter", 0.055)))
    panel_f_dot_jitter = max(0.0, float(getattr(cfg.figures, "fig21_panel_f_dot_jitter", 0.06)))
    plot_type = _normalize_plot_type(getattr(cfg.figures, "fig21_plot_type", "violin"), default="violin")
    violin_bw_method = _normalize_violin_bw_method(getattr(cfg.figures, "fig21_violin_bw_method", None))

    stats_enabled = bool(getattr(cfg.figures, "fig21_stats_enabled", True))
    stats_permutations = max(100, int(getattr(cfg.figures, "fig21_stats_permutations", 4000)))
    stats_seed = int(getattr(cfg.figures, "fig21_stats_seed", 1))
    stats_alpha = float(getattr(cfg.figures, "fig21_stats_alpha", 0.05))
    stats_pairwise_test = _normalize_pairwise_test(getattr(cfg.figures, "fig21_stats_pairwise_test", "sign_test_exact"))
    panel_f_stats_pairwise_test = _normalize_panel_f_pairwise_test(
        panel_f_stats_pairwise_test_raw or stats_pairwise_test
    )
    stats_multiple_test_correction = bool(getattr(cfg.figures, "fig21_stats_multiple_test_correction", True))
    stats_show_pairwise = bool(getattr(cfg.figures, "fig21_stats_show_pairwise", True))
    pair_group_summary_stat = str(getattr(cfg.figures, "fig21_pair_group_summary_stat", "median")).strip().lower()
    if pair_group_summary_stat not in {"mean", "median"}:
        pair_group_summary_stat = "median"
    summary_stat_title = "Mean" if pair_group_summary_stat == "mean" else "Median"
    # Allow panel_b.y_axis_label to include a dynamic placeholder driven by
    # fig21.pair_group_summary_stat. Supported placeholders:
    #   - {{}}                    -> Mean / Median
    #   - {summary_stat_title}    -> Mean / Median
    #   - {summary_stat}          -> mean / median
    panel_b_y_axis_label = (
        str(panel_b_y_axis_label)
        .replace("{{}}", summary_stat_title)
        .replace("{summary_stat_title}", summary_stat_title)
        .replace("{summary_stat}", pair_group_summary_stat)
    )

    enable_fig21 = bool(getattr(cfg.figures, "fig21", True)) and bool(
        getattr(cfg.figures, "fig21_show_pair_group_distance_contrast", True)
    )
    enable_fig22 = bool(getattr(cfg.figures, "fig22", True))
    enable_fig23 = bool(getattr(cfg.figures, "fig23", True))
    enable_fig24 = bool(getattr(cfg.figures, "fig24", True))
    enable_fig25 = bool(getattr(cfg.figures, "fig25", True))
    enable_fig26 = bool(getattr(cfg.figures, "fig26", True))
    if not any([enable_fig21, enable_fig22, enable_fig23, enable_fig24, enable_fig25, enable_fig26]):
        log.info("Figures 21-26: all outputs disabled by config.")
        return 0

    def _write_common_inputs(out_dir: Path) -> None:
        out_dir.mkdir(parents=True, exist_ok=True)
        long.to_csv(out_dir / "distance_points_filtered.tsv", sep="\t", index=False)
        metric_family_def.to_csv(out_dir / "metric_family_definition.tsv", sep="\t", index=False)

    if enable_fig21:
        out = fig21_dir
        _write_common_inputs(out)
        endpoints_df, pairwise_df, sample_summary = _pair_group_distance_contrast_tables(
            long=long,
            method_order=method_order,
            metric_order=metric_order,
            agnostic_metrics=agnostic_metrics,
            aware_metrics=aware_metrics,
            summary_stat=pair_group_summary_stat,
            denoising_methods=panel_b_denoising_methods,
            clustering_methods=panel_b_clustering_methods,
        )
        if not endpoints_df.empty:
            endpoints_df.to_csv(out / "pair_group_distance_contrast_endpoint_values.tsv", sep="\t", index=False)
        if not pairwise_df.empty:
            pairwise_df.to_csv(out / "pair_group_distance_contrast_pairwise_differences.tsv", sep="\t", index=False)
        deflation_points = _method_metric_family_distance_points(
            # Panel F is independent of Panel B family assignment, so start
            # from all configured method/metric distances rather than the
            # Panel-B-filtered endpoint table.
            endpoints_df=long,
            agnostic_metrics=agnostic_metrics,
            aware_metrics=aware_metrics,
            methods=panel_f_methods,
        )
        if not deflation_points.empty:
            deflation_points.to_csv(out / "distance_by_method_metric_family_points.tsv", sep="\t", index=False)
        deflation_stats = _write_method_metric_family_stats(
            points=deflation_points,
            methods=panel_f_methods,
            out_dir=out,
            stats_enabled=stats_enabled,
            stats_permutations=stats_permutations,
            stats_seed=stats_seed,
            stats_alpha=stats_alpha,
            stats_pairwise_test=stats_pairwise_test,
            stats_multiple_test_correction=stats_multiple_test_correction,
            panel_f_stats_pairwise_test=panel_f_stats_pairwise_test,
            panel_f_stats_alternative=panel_f_stats_alternative,
            panel_f_stats_practical_difference_delta=panel_f_stats_practical_difference_delta,
            panel_f_stats_bootstrap_resamples=panel_f_stats_bootstrap_resamples,
            panel_f_stats_ci_level=panel_f_stats_ci_level,
            panel_f_stats_planned_comparisons=panel_f_stats_planned_comparisons,
        )
        directional_stats = None
        if not sample_summary.empty:
            sample_summary.to_csv(out / "pair_group_distance_contrast_sample_summary.tsv", sep="\t", index=False)
            _pair_group_summary(sample_summary).to_csv(
                out / "pair_group_distance_contrast_pair_group_summary.tsv",
                sep="\t",
                index=False,
            )

            if stats_enabled and not sens.empty:
                directional_points = sens[sens["method"].astype(str).isin(panel_d_methods)].copy()
                directional_points.to_csv(out / "directional_sensitivity_method_points.tsv", sep="\t", index=False)
                directional_stats = _compute_panel_stats_permutation(
                    df=directional_points.rename(columns={"sensitivity": "distance"}).copy(),
                    method_order=panel_d_methods,
                    n_perm=stats_permutations,
                    seed=stats_seed + 21041,
                    alpha=stats_alpha,
                    pairwise_test=stats_pairwise_test,
                    multiple_test_correction=stats_multiple_test_correction,
                )
                pd.DataFrame(
                    [
                        {
                            "metric": "directional_method_sensitivity",
                            "label": "Directional method-level sensitivity",
                            "stats_model": directional_stats.get("stats_model", "paired_permutation"),
                            "pairwise_test": stats_pairwise_test,
                            "multiple_test_correction": bool(
                                directional_stats.get("multiple_test_correction", stats_multiple_test_correction)
                            ),
                            "n_blocks": directional_stats.get("n_blocks", 0),
                            "omnibus_stat": directional_stats.get("omnibus_stat", np.nan),
                            "omnibus_p": directional_stats.get("omnibus_p", np.nan),
                            "omnibus_significant_alpha": bool(
                                np.isfinite(float(directional_stats.get("omnibus_p", np.nan)))
                                and float(directional_stats.get("omnibus_p", np.nan)) < stats_alpha
                            ),
                            "alpha": stats_alpha,
                            "n_permutations": stats_permutations,
                        }
                    ]
                ).to_csv(out / "directional_sensitivity_method_stats_omnibus.tsv", sep="\t", index=False)
                directional_pairwise_rows: list[dict] = []
                for row in directional_stats.get("pairwise", []):
                    directional_pairwise_rows.append(
                        {
                            "metric": "directional_method_sensitivity",
                            "label": "Directional method-level sensitivity",
                            "method_a": row.get("method_a", ""),
                            "method_b": row.get("method_b", ""),
                            "method_a_label": display_method(str(row.get("method_a", ""))),
                            "method_b_label": display_method(str(row.get("method_b", ""))),
                            "stats_model": directional_stats.get("stats_model", "paired_permutation"),
                            "pairwise_test": row.get("pairwise_test", stats_pairwise_test),
                            "statistic": row.get("statistic", row.get("w_plus", np.nan)),
                            "w_plus": row.get("w_plus", np.nan),
                            "n_nonzero_pairs": row.get("n_nonzero_pairs", np.nan),
                            "n_pos_pairs": row.get("n_pos_pairs", np.nan),
                            "n_neg_pairs": row.get("n_neg_pairs", np.nan),
                            "n_zero_pairs": row.get("n_zero_pairs", np.nan),
                            "median_diff_a_minus_b": row.get("median_diff_a_minus_b", np.nan),
                            "estimate": row.get("estimate", row.get("median_diff_a_minus_b", np.nan)),
                            "std_error": row.get("std_error", np.nan),
                            "model_family": row.get("model_family", ""),
                            "p_raw": row.get("p_raw", np.nan),
                            "p_adj_bh": row.get("p_adj_bh", np.nan),
                            "p_value_used": row.get("p_value_used", row.get("p_adj_bh", np.nan)),
                            "p_adjustment": row.get("p_adjustment", "bh"),
                            "multiple_test_correction": bool(
                                row.get("multiple_test_correction", stats_multiple_test_correction)
                            ),
                            "significant_alpha": row.get("significant_alpha", False),
                            "significant_bh_alpha": row.get("significant_bh_alpha", False),
                            "alpha": stats_alpha,
                            "n_permutations": stats_permutations,
                        }
                    )
                pd.DataFrame(directional_pairwise_rows).to_csv(
                    out / "directional_sensitivity_method_stats_pairwise.tsv",
                    sep="\t",
                    index=False,
                )

            st, planned_sig = _write_pair_group_contrast_stats(
                sample_summary=sample_summary,
                out_dir=out,
                stats_enabled=stats_enabled,
                stats_permutations=stats_permutations,
                stats_seed=stats_seed,
                stats_alpha=stats_alpha,
                stats_pairwise_test=stats_pairwise_test,
                stats_multiple_test_correction=stats_multiple_test_correction,
            )
            _plot_pair_group_distance_contrast(
                sample_summary=sample_summary,
                pairwise_df=pairwise_df,
                endpoints_df=endpoints_df,
                sens=sens,
                benchmark_outdir=outdir,
                out_dir=out,
                title=fig21_title,
                panel_a_title=panel_a_title,
                panel_a_a_title=panel_a_a_title,
                panel_a_b_title=panel_a_b_title,
                panel_b_title=panel_b_title,
                panel_c_title=panel_c_title,
                panel_c_a_title=panel_c_a_title,
                panel_c_b_title=panel_c_b_title,
                panel_d_title=panel_d_title,
                panel_e_title=panel_e_title,
                panel_f_title=panel_f_title,
                font_size=font_size,
                panel_a_font_size=panel_a_font_size,
                panel_b_font_size=panel_b_font_size,
                panel_c_font_size=panel_c_font_size,
                panel_d_font_size=panel_d_font_size,
                panel_e_font_size=panel_e_font_size,
                panel_f_font_size=panel_f_font_size,
                panel_e_grid_label_font_size=panel_e_grid_label_font_size,
                panel_e_grid_tick_font_size=panel_e_grid_tick_font_size,
                panel_e_grid_text_font_size=panel_e_grid_text_font_size,
                panel_e_grid_legend_font_size=panel_e_grid_legend_font_size,
                panel_e_abs_density_label=panel_e_abs_density_label,
                panel_e_summary_density_label=panel_e_summary_density_label,
                panel_e_show_summary_density=panel_e_show_summary_density,
                panel_e_mean_label=panel_e_mean_label,
                panel_e_sample_v4_multiple_label=panel_e_sample_v4_multiple_label,
                panel_e_sample_v4_single_label=panel_e_sample_v4_single_label,
                panel_a_x_axis_label=panel_a_x_axis_label,
                panel_a_y_axis_label=panel_a_y_axis_label,
                panel_b_x_axis_label=panel_b_x_axis_label,
                panel_b_y_axis_label=panel_b_y_axis_label,
                panel_c_x_axis_label=panel_c_x_axis_label,
                panel_c_y_axis_label=panel_c_y_axis_label,
                panel_d_x_axis_label=panel_d_x_axis_label,
                panel_d_y_axis_label=panel_d_y_axis_label,
                panel_e_x_axis_label=panel_e_x_axis_label,
                panel_e_y_axis_label=panel_e_y_axis_label,
                panel_f_x_axis_label=panel_f_x_axis_label,
                panel_f_y_axis_label=panel_f_y_axis_label,
                panel_f_xtick_rotation=panel_f_xtick_rotation,
                panel_b_pair_groups_to_show=panel_b_pair_groups_to_show,
                panel_b_group_labels=panel_b_group_labels,
                show_panel_a=show_panel_a,
                show_panel_b=show_panel_b,
                show_panel_c=show_panel_c,
                show_panel_d=show_panel_d,
                show_panel_e=show_panel_e,
                show_panel_f=show_panel_f,
                sample_v4_multiple_label=sample_v4_multiple_label,
                sample_v4_single_label=sample_v4_single_label,
                panel_a_xtick_rotation=panel_a_xtick_rotation,
                panel_a_ytick_rotation=panel_a_ytick_rotation,
                panel_b_xtick_rotation=panel_b_xtick_rotation,
                panel_b_ytick_rotation=panel_b_ytick_rotation,
                endpoint_matrix_gap_cells=endpoint_matrix_gap_cells,
                endpoint_matrix_xtick_rotation=endpoint_matrix_xtick_rotation,
                endpoint_matrix_ytick_rotation=endpoint_matrix_ytick_rotation,
                panel_c_a_xtick_rotation=panel_c_a_xtick_rotation,
                panel_c_a_ytick_rotation=panel_c_a_ytick_rotation,
                panel_c_b_xtick_rotation=panel_c_b_xtick_rotation,
                panel_c_b_ytick_rotation=panel_c_b_ytick_rotation,
                pair_group_summary_stat=pair_group_summary_stat,
                panel_a_show_planned_comparison_annotations=panel_a_show_planned_comparison_annotations,
                show_panel_e_raincloud_grid=show_panel_e_raincloud_grid,
                panel_b_dot_size=panel_b_dot_size,
                panel_b_dot_opacity=panel_b_dot_opacity,
                panel_b_dot_jitter=panel_b_dot_jitter,
                panel_d_dot_size=panel_d_dot_size,
                panel_d_dot_opacity=panel_d_dot_opacity,
                panel_d_dot_jitter=panel_d_dot_jitter,
                panel_e_dot_size=panel_e_dot_size,
                panel_e_dot_opacity=panel_e_dot_opacity,
                panel_e_dot_jitter=panel_e_dot_jitter,
                panel_f_dot_size=panel_f_dot_size,
                panel_f_dot_opacity=panel_f_dot_opacity,
                panel_f_dot_jitter=panel_f_dot_jitter,
                point_size=point_size,
                plot_type=plot_type,
                method_order=method_order,
                panel_d_methods=panel_d_methods,
                panel_f_methods=panel_f_methods,
                panel_f_summary_statistic=panel_f_summary_statistic,
                panel_f_stats_show_annotations=panel_f_stats_show_annotations,
                panel_f_stats_show_non_significant_annotations=panel_f_stats_show_non_significant_annotations,
                panel_d_stats_show_annotations=panel_d_stats_show_annotations,
                metric_order=metric_order,
                agnostic_metrics=agnostic_metrics,
                aware_metrics=aware_metrics,
                violin_bw_method=violin_bw_method,
                stats=st,
                directional_stats=directional_stats,
                deflation_points=deflation_points,
                deflation_stats=deflation_stats,
                planned_sig_rows=planned_sig,
                stats_enabled=stats_enabled,
                stats_show_pairwise=stats_show_pairwise,
                panel_b_show_within_contrast_stats=panel_b_show_within_contrast_stats,
                stats_alpha=stats_alpha,
            )
        else:
            log.info("Figure 21: pair-group distance contrast table is empty.")

    if enable_fig22:
        out = fig22_dir
        _write_common_inputs(out)
        sens.to_csv(out / "sensitivity_by_method_points.tsv", sep="\t", index=False)
        _method_summary(sens).to_csv(out / "sensitivity_by_method_summary.tsv", sep="\t", index=False)
        st = _write_sensitivity_stats(
            sens=sens,
            out_dir=out,
            method_order=method_order,
            stats_enabled=stats_enabled,
            stats_permutations=stats_permutations,
            stats_seed=stats_seed,
            stats_alpha=stats_alpha,
            stats_pairwise_test=stats_pairwise_test,
            stats_multiple_test_correction=stats_multiple_test_correction,
        )
        _plot_sensitivity_by_method(
            sens=sens,
            out_dir=out,
            method_order=method_order,
            title=fig22_title,
            plot_type=plot_type,
            point_size=point_size,
            font_size=fig22_font_size,
            agnostic_metrics=agnostic_metrics,
            aware_metrics=aware_metrics,
            violin_bw_method=violin_bw_method,
            stats=st,
            stats_enabled=stats_enabled,
            stats_show_pairwise=stats_show_pairwise,
            stats_alpha=stats_alpha,
        )

    if enable_fig23:
        out = fig23_dir
        _write_common_inputs(out)
        _plot_sample_method_heatmap(sens=sens, out_dir=out, method_order=method_order, font_size=fig23_font_size)

    if enable_fig24:
        out = fig24_dir
        _write_common_inputs(out)
        _plot_per_sample_method_metric_pcoa(
            long=long,
            out_dir=out,
            method_order=method_order,
            metric_order=metric_order,
            font_size=fig24_font_size,
            point_size=point_size,
        )

    method_metric_vectors = pd.DataFrame()
    method_metric_cols: list[str] = []
    if enable_fig25 or enable_fig26:
        method_metric_vectors, method_metric_cols = _method_metric_profile_vectors(
            long=long,
            method_order=method_order,
            metric_order=metric_order,
            agnostic_metrics=agnostic_metrics,
            aware_metrics=aware_metrics,
        )

    if enable_fig25:
        out = fig25_dir
        _write_common_inputs(out)
        _write_pcoa_outputs(
            vectors=method_metric_vectors,
            value_cols=method_metric_cols,
            out_dir=out,
            stem="method_metric_profile_pcoa",
            title=fig25_title,
            cfg=cfg,
            font_size=fig25_font_size,
            color_col="mean_sensitivity",
        )

    if enable_fig26:
        out = fig26_dir
        _write_common_inputs(out)
        diff_vectors, diff_cols = _pairwise_difference_vectors(method_metric_vectors, method_metric_cols)
        _write_pcoa_outputs(
            vectors=diff_vectors,
            value_cols=diff_cols,
            out_dir=out,
            stem="pairwise_difference_pcoa",
            title=fig26_title,
            cfg=cfg,
            font_size=fig26_font_size,
            color_col="mean_sensitivity",
        )

    if enable_fig21:
        _write_text(
            fig21_dir / "README.md",
            "# Figure 21\n\nPair-group distance contrast diagnostics.\n\n"
            "Includes optional Panel E: KDE + horizontal-boxplot distributions for each pair-group cell.\n",
        )
    if enable_fig22:
        _write_text(fig22_dir / "README.md", "# Figure 22\n\nSensitivity-by-method diagnostics.\n")
    if enable_fig23:
        _write_text(fig23_dir / "README.md", "# Figure 23\n\nSample x method sensitivity heatmap.\n")
    if enable_fig24:
        _write_text(fig24_dir / "README.md", "# Figure 24\n\nPer-sample method x metric PCoA diagnostics.\n")
    if enable_fig25:
        _write_text(fig25_dir / "README.md", "# Figure 25\n\nMethod x metric profile PCoA diagnostics.\n")
    if enable_fig26:
        _write_text(fig26_dir / "README.md", "# Figure 26\n\nPairwise-difference profile PCoA diagnostics.\n")

    log.info(
        "Wrote metric-family diagnostics: "
        + ", ".join(
            [
                name
                for name, enabled in [
                    ("fig21_pair_group_distance_contrast", enable_fig21),
                    ("fig22_sensitivity_by_method", enable_fig22),
                    ("fig23_sample_method_heatmap", enable_fig23),
                    ("fig24_per_sample_method_metric_pcoa", enable_fig24),
                    ("fig25_method_metric_profile_pcoa", enable_fig25),
                    ("fig26_pairwise_difference_pcoa", enable_fig26),
                ]
                if enabled
            ]
        )
    )
    return 0
