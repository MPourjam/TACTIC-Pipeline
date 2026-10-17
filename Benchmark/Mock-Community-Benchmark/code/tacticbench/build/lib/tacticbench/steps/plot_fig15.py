from __future__ import annotations

from fnmatch import fnmatch
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..log import get_logger
from ..utils.labels import display_methods
from ..utils.paths import figures_root
from . import plot_fig09

METHOD_ORDER = ["ASV", "zOTU", "OTU", "TIC", "TAC"]
METHOD_COLORS = {
    "ASV": "#0072B2",
    "zOTU": "#56B4E9",
    "OTU": "#D55E00",
    "TIC": "#009E73",
    "TAC": "#E69F00",
}

SEGMENT_COLORS = {
    "mapped": "#009E73",
    "unmapped_with_best_hit": "#E69F00",
    "no_hit": "#D55E00",
}


def _is_excluded(cid: str, patterns: list[str]) -> bool:
    for pat in patterns:
        p = str(pat).strip()
        if not p:
            continue
        if fnmatch(cid, p) or (p in cid):
            return True
    return False


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


def _parse_method_threshold_map(raw) -> dict[str, float]:
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


def _fig15_no_hit_thresholds(cfg) -> dict[str, float]:
    """
    Build Figure 15 no-hit effective-identity thresholds with priority:
      1) fig12.no_hit_min_effective_pident_by_method
      2) fig09.no_hit_min_effective_pident_by_method
      3) fig15.no_hit_min_effective_pident_by_method
    """
    out: dict[str, float] = {}
    out.update(_parse_method_threshold_map(getattr(cfg.figures, "fig12_no_hit_min_effective_pident_by_method", {}) or {}))
    out.update(_parse_method_threshold_map(getattr(cfg.figures, "fig09_no_hit_min_effective_pident_by_method", {}) or {}))
    out.update(_parse_method_threshold_map(getattr(cfg.figures, "fig15_no_hit_min_effective_pident_by_method", {}) or {}))
    return out


def _compute_has_best_hit(
    df: pd.DataFrame,
    effective_identity_threshold_by_method: dict[str, float] | None = None,
) -> pd.Series:
    eff = pd.to_numeric(df.get("best_hit_effective_pident", np.nan), errors="coerce")
    pid = pd.to_numeric(df.get("best_hit_pident", np.nan), errors="coerce")

    # Backward compatibility: some older exports store effective identity in [0,1].
    finite_eff = eff[np.isfinite(eff.to_numpy(dtype=float))]
    if finite_eff.size > 0 and float(np.nanmax(finite_eff.to_numpy(dtype=float))) <= 1.0:
        eff = eff * 100.0

    score = eff.where(eff.notna(), pid)
    has = score.notna()

    thr_map = effective_identity_threshold_by_method or {}
    if thr_map and "method" in df.columns:
        methods = df["method"].astype(str)
        thr = methods.map(thr_map)
        use_thr = thr.notna()
        pass_thr = score >= pd.to_numeric(thr, errors="coerce")
        has = has & (~use_thr | pass_thr.fillna(False))

    return has.astype(float)


def _collapse_to_single_feature_level(df: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse duplicated unmapped-feature rows across graph target levels.
    Preference order: tip > species > other.
    """
    if df.empty:
        return df.copy()
    if "graph_target_level_plot" not in df.columns:
        out = df.copy()
        out["graph_target_level_plot"] = "single"
        return out

    out = df.copy()
    lvl = out["graph_target_level_plot"].astype(str).str.strip().str.lower()
    out["_lvl_rank"] = np.where(lvl == "tip", 0, np.where(lvl == "species", 1, 2))
    key_cols = ["community_id", "sample_id", "sample_key", "method", "unmapped_feature"]
    missing = [c for c in key_cols if c not in out.columns]
    if missing:
        out = out.drop(columns=["_lvl_rank"])
        out["graph_target_level_plot"] = "single"
        return out

    out = out.sort_values(key_cols + ["_lvl_rank"], kind="mergesort")
    out = out.drop_duplicates(subset=key_cols, keep="first").copy()
    out = out.drop(columns=["_lvl_rank"])
    out["graph_target_level_plot"] = "single"
    return out


def _build_sample_method_table(
    level_df: pd.DataFrame,
    ref_pairs: pd.DataFrame,
    no_hit_thresholds: dict[str, float],
) -> pd.DataFrame:
    if ref_pairs.empty:
        ref = (
            level_df[[c for c in ["community_id", "sample_id", "sample_key", "method"] if c in level_df.columns]]
            .drop_duplicates()
            .copy()
        )
    else:
        ref = ref_pairs.copy()

    if level_df.empty:
        out = ref.copy()
        out["n_unmapped_features"] = 0.0
        out["n_unmapped_with_best_hit_features"] = 0.0
        out["n_no_hit_features"] = 0.0
        out["mapped_fraction"] = 1.0
        out["unmapped_with_best_hit_fraction"] = 0.0
        out["no_hit_fraction"] = 0.0
        out["unmapped_total_fraction"] = 0.0
        return out

    dd = level_df.copy()
    dd["sample_key"] = dd["sample_key"].astype(str)
    dd["method"] = dd["method"].astype(str)
    dd["has_best_hit"] = _compute_has_best_hit(
        dd,
        effective_identity_threshold_by_method=no_hit_thresholds,
    )
    dd["inferred_feature_fraction"] = (
        pd.to_numeric(dd.get("inferred_feature_fraction", np.nan), errors="coerce")
        .fillna(0.0)
        .clip(lower=0.0)
    )
    dd["weighted_with_best_hit"] = dd["inferred_feature_fraction"] * dd["has_best_hit"]

    agg = (
        dd.groupby(["sample_key", "method"], as_index=False)
        .agg(
            n_unmapped_features=("unmapped_feature", "count"),
            n_unmapped_with_best_hit_features=("has_best_hit", "sum"),
            unmapped_total_fraction=("inferred_feature_fraction", "sum"),
            unmapped_with_best_hit_fraction=("weighted_with_best_hit", "sum"),
        )
    )
    agg["n_no_hit_features"] = agg["n_unmapped_features"] - agg["n_unmapped_with_best_hit_features"]
    agg["no_hit_fraction"] = agg["unmapped_total_fraction"] - agg["unmapped_with_best_hit_fraction"]
    agg["unmapped_total_fraction"] = agg["unmapped_total_fraction"].clip(lower=0.0, upper=1.0)
    agg["unmapped_with_best_hit_fraction"] = (
        agg["unmapped_with_best_hit_fraction"].clip(lower=0.0).clip(upper=agg["unmapped_total_fraction"])
    )
    agg["no_hit_fraction"] = agg["no_hit_fraction"].clip(lower=0.0).clip(upper=agg["unmapped_total_fraction"])
    agg["mapped_fraction"] = (1.0 - agg["unmapped_total_fraction"]).clip(lower=0.0, upper=1.0)

    out = ref.merge(agg, on=["sample_key", "method"], how="left")
    for c in [
        "n_unmapped_features",
        "n_unmapped_with_best_hit_features",
        "n_no_hit_features",
        "mapped_fraction",
        "unmapped_with_best_hit_fraction",
        "no_hit_fraction",
        "unmapped_total_fraction",
    ]:
        out[c] = pd.to_numeric(out[c], errors="coerce")

    out["n_unmapped_features"] = out["n_unmapped_features"].fillna(0.0)
    out["n_unmapped_with_best_hit_features"] = out["n_unmapped_with_best_hit_features"].fillna(0.0)
    out["n_no_hit_features"] = out["n_no_hit_features"].fillna(0.0)
    out["unmapped_total_fraction"] = out["unmapped_total_fraction"].fillna(0.0).clip(0.0, 1.0)
    out["unmapped_with_best_hit_fraction"] = out["unmapped_with_best_hit_fraction"].fillna(0.0).clip(0.0, 1.0)
    out["no_hit_fraction"] = out["no_hit_fraction"].fillna(0.0).clip(0.0, 1.0)
    out["mapped_fraction"] = out["mapped_fraction"].fillna(1.0).clip(0.0, 1.0)
    return out


def _build_method_summary(sample_method_df: pd.DataFrame) -> pd.DataFrame:
    if sample_method_df.empty:
        return pd.DataFrame(
            columns=[
                "method",
                "n_samples",
                "mean_mapped_fraction",
                "mean_unmapped_with_best_hit_fraction",
                "mean_no_hit_fraction",
                "mean_unmapped_total_fraction",
                "sum_n_unmapped_features",
                "sum_n_no_hit_features",
            ]
        )

    rows: list[dict[str, float | str]] = []
    for method in METHOD_ORDER:
        sub = sample_method_df[sample_method_df["method"].astype(str) == method].copy()
        if sub.empty:
            rows.append(
                {
                    "method": method,
                    "n_samples": 0.0,
                    "mean_mapped_fraction": 0.0,
                    "mean_unmapped_with_best_hit_fraction": 0.0,
                    "mean_no_hit_fraction": 0.0,
                    "mean_unmapped_total_fraction": 0.0,
                    "sum_n_unmapped_features": 0.0,
                    "sum_n_no_hit_features": 0.0,
                }
            )
            continue
        rows.append(
            {
                "method": method,
                "n_samples": float(sub["sample_key"].astype(str).nunique()),
                "mean_mapped_fraction": float(pd.to_numeric(sub["mapped_fraction"], errors="coerce").fillna(0.0).mean()),
                "mean_unmapped_with_best_hit_fraction": float(
                    pd.to_numeric(sub["unmapped_with_best_hit_fraction"], errors="coerce").fillna(0.0).mean()
                ),
                "mean_no_hit_fraction": float(pd.to_numeric(sub["no_hit_fraction"], errors="coerce").fillna(0.0).mean()),
                "mean_unmapped_total_fraction": float(
                    pd.to_numeric(sub["unmapped_total_fraction"], errors="coerce").fillna(0.0).mean()
                ),
                "sum_n_unmapped_features": float(
                    pd.to_numeric(sub["n_unmapped_features"], errors="coerce").fillna(0.0).sum()
                ),
                "sum_n_no_hit_features": float(
                    pd.to_numeric(sub["n_no_hit_features"], errors="coerce").fillna(0.0).sum()
                ),
            }
        )
    return pd.DataFrame(rows)


def _plot_method_stacked_composition(
    summary_df: pd.DataFrame,
    out_path: Path,
    title: str,
    font_size: float,
) -> None:
    tick_size = max(7.0, font_size - 2.0)
    fig = plt.figure(figsize=(10.0, 5.4))
    ax = fig.add_subplot(1, 1, 1)

    x = np.arange(len(METHOD_ORDER))
    mapped = []
    with_hit = []
    no_hit = []
    n_samples = []
    for m in METHOD_ORDER:
        row = summary_df[summary_df["method"].astype(str) == m]
        if row.empty:
            mapped.append(0.0)
            with_hit.append(0.0)
            no_hit.append(0.0)
            n_samples.append(0)
            continue
        r = row.iloc[0]
        mapped.append(float(pd.to_numeric(r.get("mean_mapped_fraction", 0.0), errors="coerce")))
        with_hit.append(float(pd.to_numeric(r.get("mean_unmapped_with_best_hit_fraction", 0.0), errors="coerce")))
        no_hit.append(float(pd.to_numeric(r.get("mean_no_hit_fraction", 0.0), errors="coerce")))
        n_samples.append(int(pd.to_numeric(r.get("n_samples", 0.0), errors="coerce")))

    mapped = np.clip(np.array(mapped, dtype=float), 0.0, 1.0)
    with_hit = np.clip(np.array(with_hit, dtype=float), 0.0, 1.0)
    no_hit = np.clip(np.array(no_hit, dtype=float), 0.0, 1.0)
    total = mapped + with_hit + no_hit
    with np.errstate(invalid="ignore", divide="ignore"):
        scale = np.where(total > 1.0, 1.0 / total, 1.0)
    mapped *= scale
    with_hit *= scale
    no_hit *= scale

    ax.bar(
        x,
        mapped * 100.0,
        color=SEGMENT_COLORS["mapped"],
        edgecolor="#333333",
        linewidth=0.7,
        width=0.72,
        label="Mapped reads",
    )
    ax.bar(
        x,
        with_hit * 100.0,
        bottom=mapped * 100.0,
        color=SEGMENT_COLORS["unmapped_with_best_hit"],
        edgecolor="#333333",
        linewidth=0.7,
        width=0.72,
        label="Unmapped reads\n(best hit found)",
    )
    ax.bar(
        x,
        no_hit * 100.0,
        bottom=(mapped + with_hit) * 100.0,
        color=SEGMENT_COLORS["no_hit"],
        edgecolor="#333333",
        linewidth=0.7,
        width=0.72,
        label="No-hit reads",
    )

    for i in range(len(METHOD_ORDER)):
        unmapped_total = float(with_hit[i] + no_hit[i])
        ax.text(
            x[i],
            min(104.0, 100.0 * (mapped[i] + with_hit[i] + no_hit[i]) + 1.4),
            f"{100.0 * unmapped_total:.1f}%",
            ha="center",
            va="bottom",
            fontsize=max(6.0, tick_size - 1.0),
            color="#222222",
        )
        ax.text(
            x[i],
            -6.0,
            f"n={n_samples[i]}",
            ha="center",
            va="top",
            fontsize=max(6.0, tick_size - 1.0),
            color="#333333",
            clip_on=False,
        )

    ax.set_title("Read-Weighted Mapping Outcomes by Method", fontsize=font_size)
    ax.set_xlabel("Method", fontsize=font_size)
    ax.set_ylabel("Mean read fraction across samples (%)", fontsize=font_size)
    ax.set_xticks(x)
    ax.set_xticklabels(display_methods(METHOD_ORDER), fontsize=tick_size)
    ax.set_ylim(0.0, 108.0)
    ax.set_yticks(np.arange(0.0, 101.0, 10.0))
    ax.tick_params(axis="y", labelsize=tick_size)
    ax.grid(True, axis="y", alpha=0.25, linewidth=0.6)
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.935),
        ncol=3,
        frameon=False,
        fontsize=tick_size,
    )

    fig.suptitle(title, fontsize=max(font_size + 1.0, 10.0), y=0.985)
    fig.subplots_adjust(bottom=0.2, right=0.98, top=0.80)
    fig.savefig(out_path, format="svg")
    plt.close(fig)


def _violin_metric_specs() -> list[tuple[str, str]]:
    return [
        ("mapped_fraction", "Mapped Read Fraction"),
        ("unmapped_with_best_hit_fraction", "Unmapped with Best Hit Fraction"),
        ("no_hit_fraction", "No-Hit Read Fraction"),
    ]


def _build_violin_long_table(sample_method_df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for metric, label in _violin_metric_specs():
        if metric not in sample_method_df.columns:
            continue
        tmp = sample_method_df.copy()
        tmp["value"] = pd.to_numeric(tmp[metric], errors="coerce")
        tmp = tmp[np.isfinite(tmp["value"].to_numpy(dtype=float))].copy()
        if tmp.empty:
            continue
        for _, r in tmp.iterrows():
            rows.append(
                {
                    "community_id": str(r.get("community_id", "")),
                    "sample_id": str(r.get("sample_id", "")),
                    "sample_key": str(r.get("sample_key", "")),
                    "method": str(r.get("method", "")),
                    "metric": metric,
                    "metric_label": label,
                    "value": float(r["value"]),
                    "value_percent": 100.0 * float(r["value"]),
                }
            )
    return pd.DataFrame(rows)


def _compute_violin_stats_bundle(
    cfg,
    long_df: pd.DataFrame,
    out_dir: Path,
    stem: str,
    log,
) -> dict[str, dict]:
    """
    Reuse Figure 09 stats engine for Figure 15 violin metrics.
    """
    stats_enabled = bool(getattr(cfg.figures, "fig09_stats_enabled", True))
    if not stats_enabled or long_df.empty:
        return {}

    stats_permutations = max(100, int(getattr(cfg.figures, "fig09_stats_permutations", 4000)))
    stats_seed = int(getattr(cfg.figures, "fig09_stats_seed", 1))
    stats_alpha = float(getattr(cfg.figures, "fig09_stats_alpha", 0.05))
    stats_model = plot_fig09._normalize_stats_model(getattr(cfg.figures, "fig09_stats_model", "mixed_effect_glmmtmb"))
    stats_pairwise_test = plot_fig09._normalize_pairwise_test(
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
    glmm_count_family = plot_fig09._normalize_glmm_count_family(
        getattr(cfg.figures, "fig09_stats_glmm_count_family", "nbinom2")
    )
    glmm_emmeans_scale = plot_fig09._normalize_emmeans_scale(
        getattr(cfg.figures, "fig09_stats_glmm_emmeans_scale", "link")
    )
    glmm_spec = {
        "binomial_formula": glmm_bin_formula,
        "binomial_null_formula": glmm_bin_null_formula,
        "count_formula": glmm_count_formula,
        "count_null_formula": glmm_count_null_formula,
        "count_family": glmm_count_family,
        "emmeans_scale": glmm_emmeans_scale,
    }

    stats_by_metric: dict[str, dict] = {}
    omnibus_rows: list[dict[str, object]] = []
    pairwise_rows: list[dict[str, object]] = []

    for idx, (metric, label) in enumerate(_violin_metric_specs(), start=1):
        d = long_df[long_df["metric"].astype(str) == metric].copy()
        if d.empty:
            continue
        d = d.rename(columns={"value": metric})
        st = plot_fig09._compute_panel_stats(
            df=d,
            metric_col=metric,
            method_order=METHOD_ORDER,
            n_perm=stats_permutations,
            seed=stats_seed + 3000000 + (idx * 10000),
            alpha=stats_alpha,
            pairwise_test=stats_pairwise_test,
            multiple_test_correction=stats_multiple_test_correction,
            stats_model=stats_model,
            glmm_spec=glmm_spec,
            log=log,
        )
        stats_by_metric[metric] = st

        omnibus_rows.append(
            {
                "metric": metric,
                "label": label,
                "stats_model": st.get("stats_model", stats_model),
                "pairwise_test": stats_pairwise_test,
                "multiple_test_correction": bool(
                    st.get("multiple_test_correction", stats_multiple_test_correction)
                ),
                "n_blocks": st.get("n_blocks", np.nan),
                "friedman_q": st.get("friedman_q", np.nan),
                "friedman_p_perm": st.get("friedman_p", np.nan),
                "friedman_significant_alpha": bool(
                    np.isfinite(float(st.get("friedman_p", np.nan)))
                    and float(st.get("friedman_p", np.nan)) < stats_alpha
                ),
                "alpha": stats_alpha,
                "n_permutations": stats_permutations,
            }
        )
        for row in st.get("pairwise", []):
            pairwise_rows.append(
                {
                    "metric": metric,
                    "label": label,
                    "method_a": row.get("method_a", ""),
                    "method_b": row.get("method_b", ""),
                    "stats_model": st.get("stats_model", stats_model),
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
                    "significant_alpha": row.get("significant_alpha", row.get("significant_bh", False)),
                    "significant_bh_alpha": row.get("significant_bh", False),
                    "alpha": stats_alpha,
                    "n_permutations": stats_permutations,
                }
            )

    pd.DataFrame(omnibus_rows).to_csv(out_dir / f"{stem}_stats_omnibus.tsv", sep="\t", index=False)
    pd.DataFrame(pairwise_rows).to_csv(out_dir / f"{stem}_stats_pairwise.tsv", sep="\t", index=False)
    return stats_by_metric


def _plot_violin_metrics(
    cfg,
    long_df: pd.DataFrame,
    out_path: Path,
    title: str,
    font_size: float,
    stats_by_metric: dict[str, dict],
) -> None:
    tick_size = max(7.0, font_size - 2.0)
    point_size = max(1.0, float(getattr(cfg.figures, "box_violin_dot_size", 22.0)))
    stats_alpha = float(getattr(cfg.figures, "fig09_stats_alpha", 0.05))
    stats_show_pairwise = bool(getattr(cfg.figures, "fig09_stats_show_pairwise", True))
    stats_multiple_test_correction = bool(
        getattr(cfg.figures, "fig09_stats_multiple_test_correction", True)
    )

    specs = _violin_metric_specs()
    fig = plt.figure(figsize=(14.2, 4.8))
    rng = np.random.default_rng(1)

    for i, (metric, label) in enumerate(specs, start=1):
        ax = fig.add_subplot(1, len(specs), i)
        d = long_df[long_df["metric"].astype(str) == metric].copy()

        method_data: list[tuple[int, str, np.ndarray]] = []
        for j, method in enumerate(METHOD_ORDER, start=1):
            sub = d[d["method"].astype(str) == method].copy()
            vals = pd.to_numeric(sub.get("value", np.nan), errors="coerce").to_numpy(dtype=float)
            vals = vals[np.isfinite(vals)]
            if vals.size > 0:
                method_data.append((j, method, vals))

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

        ax.set_title(label, fontsize=font_size)
        ax.set_ylabel("Read fraction", fontsize=font_size)
        ax.set_ylim(0.0, 1.0)
        ax.set_yticks(np.arange(0.0, 1.01, 0.1))
        ax.set_xticks(range(1, len(METHOD_ORDER) + 1))
        ax.set_xticklabels(
            display_methods(METHOD_ORDER),
            rotation=45,
            ha="right",
            rotation_mode="anchor",
            fontsize=tick_size,
        )
        ax.tick_params(axis="y", labelsize=tick_size)
        ax.grid(True, axis="y", alpha=0.25, linewidth=0.6)

        st = stats_by_metric.get(metric, {})
        p_omni = float(st.get("friedman_p", np.nan))
        stats_model_name = str(st.get("stats_model", "paired_permutation"))
        if np.isfinite(p_omni):
            ax.text(
                0.985,
                0.985,
                f"Omnibus p={p_omni:.3g}\nmodel={stats_model_name}",
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=max(6.0, tick_size - 1.0),
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="#bbbbbb", alpha=0.85),
                zorder=7,
            )
        if stats_show_pairwise:
            plot_fig09._draw_pairwise_brackets(
                ax=ax,
                pairwise_rows=list(st.get("pairwise_sig", [])),
                method_order=METHOD_ORDER,
                alpha=stats_alpha,
                constrain_to_axis=True,
                placement="top",
                spacing_scale=1.35,
            )
        if i == len(specs):
            pairwise_label = (
                f"Pairwise BH threshold (alpha={stats_alpha:.3g})"
                if stats_multiple_test_correction
                else f"Pairwise threshold (alpha={stats_alpha:.3g})"
            )
            ax.text(
                0.985,
                0.02,
                pairwise_label,
                transform=ax.transAxes,
                ha="right",
                va="bottom",
                fontsize=max(6.0, tick_size - 1.2),
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="#bbbbbb", alpha=0.85),
                zorder=7,
            )

    fig.suptitle(title, fontsize=max(font_size + 1.0, 10.0))
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(out_path, format="svg")
    plt.close(fig)


def run(cfg) -> int:
    """
    Figure 15:
      Read-weighted mapping outcomes by method.

    Each sample/method is decomposed into:
      - mapped read fraction
      - unmapped read fraction with a usable best hit
      - no-hit read fraction
    """
    log = get_logger()
    outdir = Path(cfg.paths.outdir)
    fig09_dir = figures_root(outdir) / "fig09_split_merge_diagnostics"
    src_unmapped = fig09_dir / "unmapped_feature_best_hit_identity.tsv"
    src_diag = fig09_dir / "split_merge_diagnostics.tsv"

    if not src_unmapped.exists():
        log.info("Figure 15: Figure 09 unmapped diagnostics not found; generating Figure 09 prerequisites.")
        plot_fig09.run(cfg)
    if not src_unmapped.exists():
        log.info("Figure 15: unmapped_feature_best_hit_identity.tsv missing; run plot-fig09 first.")
        return 0

    try:
        unmapped_df = pd.read_csv(src_unmapped, sep="\t")
    except Exception:
        log.info("Figure 15: failed to read Figure 09 unmapped diagnostics TSV.")
        return 0
    if unmapped_df.empty:
        log.info("Figure 15: no unmapped diagnostics rows to summarize.")
        return 0

    diag_df = pd.DataFrame()
    if src_diag.exists():
        try:
            diag_df = pd.read_csv(src_diag, sep="\t")
        except Exception:
            diag_df = pd.DataFrame()

    exclude_patterns = list(getattr(cfg.figures, "exclude_community_name_patterns", []) or [])
    exclude_patterns += list(getattr(cfg.figures, "exclude_community_name_patterns_fig15", []) or [])

    if "community_id" in unmapped_df.columns:
        unmapped_df = unmapped_df[
            ~unmapped_df["community_id"].astype(str).map(lambda x: _is_excluded(x, exclude_patterns))
        ].copy()
    if (not diag_df.empty) and ("community_id" in diag_df.columns):
        diag_df = diag_df[
            ~diag_df["community_id"].astype(str).map(lambda x: _is_excluded(x, exclude_patterns))
        ].copy()
    if unmapped_df.empty:
        log.info("Figure 15: all rows excluded or empty after filtering.")
        return 0

    level_df = _collapse_to_single_feature_level(unmapped_df)
    if level_df.empty:
        log.info("Figure 15: no unmapped rows after collapsing target-level duplicates.")
        return 0

    if (not diag_df.empty) and {"sample_key", "method"}.issubset(diag_df.columns):
        ref_cols = [c for c in ["community_id", "sample_id", "sample_key", "method"] if c in diag_df.columns]
        ref_pairs = diag_df[ref_cols].drop_duplicates().copy()
    else:
        ref_cols = [c for c in ["community_id", "sample_id", "sample_key", "method"] if c in level_df.columns]
        ref_pairs = level_df[ref_cols].drop_duplicates().copy()

    no_hit_thresholds = _fig15_no_hit_thresholds(cfg)
    sample_method_df = _build_sample_method_table(
        level_df=level_df,
        ref_pairs=ref_pairs,
        no_hit_thresholds=no_hit_thresholds,
    )
    if sample_method_df.empty:
        log.info("Figure 15: no sample/method rows to plot.")
        return 0

    sample_method_df["mapped_percent"] = 100.0 * sample_method_df["mapped_fraction"]
    sample_method_df["unmapped_with_best_hit_percent"] = 100.0 * sample_method_df["unmapped_with_best_hit_fraction"]
    sample_method_df["no_hit_percent"] = 100.0 * sample_method_df["no_hit_fraction"]
    sample_method_df["unmapped_total_percent"] = 100.0 * sample_method_df["unmapped_total_fraction"]

    method_summary_df = _build_method_summary(sample_method_df)
    for c in [
        "mean_mapped_fraction",
        "mean_unmapped_with_best_hit_fraction",
        "mean_no_hit_fraction",
        "mean_unmapped_total_fraction",
    ]:
        method_summary_df[c.replace("_fraction", "_percent")] = 100.0 * pd.to_numeric(
            method_summary_df[c], errors="coerce"
        ).fillna(0.0)

    out_dir = figures_root(outdir) / "fig15_read_weighted_mapping_outcomes"
    out_dir.mkdir(parents=True, exist_ok=True)

    sample_out = out_dir / "fig15_read_weighted_mapping_outcomes_by_sample.tsv"
    method_out = out_dir / "fig15_read_weighted_mapping_outcomes_by_method.tsv"
    thr_out = out_dir / "fig15_no_hit_thresholds.tsv"
    fig_out = out_dir / "fig15_read_weighted_mapping_outcomes.svg"
    violin_stem = "fig15_read_weighted_mapping_outcomes_violin"
    violin_data_out = out_dir / f"{violin_stem}_data.tsv"
    violin_fig_out = out_dir / f"{violin_stem}.svg"

    sample_method_df.to_csv(sample_out, sep="\t", index=False)
    method_summary_df.to_csv(method_out, sep="\t", index=False)
    thr_df = pd.DataFrame(
        [{"method": m, "no_hit_min_effective_pident": no_hit_thresholds.get(m, np.nan)} for m in METHOD_ORDER]
    )
    thr_df.to_csv(thr_out, sep="\t", index=False)

    title = getattr(cfg.figures, "fig15_title", None) or "Figure 15: Read-Weighted Mapping Outcomes"
    font_size = max(6.0, float(getattr(cfg.figures, "fig15_font_size", getattr(cfg.figures, "fig09_font_size", 10.0))))
    _plot_method_stacked_composition(
        summary_df=method_summary_df,
        out_path=fig_out,
        title=title,
        font_size=font_size,
    )

    long_df = _build_violin_long_table(sample_method_df)
    if not long_df.empty:
        long_df.to_csv(violin_data_out, sep="\t", index=False)
        stats_by_metric = _compute_violin_stats_bundle(
            cfg=cfg,
            long_df=long_df,
            out_dir=out_dir,
            stem=violin_stem,
            log=log,
        )
        _plot_violin_metrics(
            cfg=cfg,
            long_df=long_df,
            out_path=violin_fig_out,
            title=f"{title} (Violin View)",
            font_size=font_size,
            stats_by_metric=stats_by_metric,
        )
    else:
        log.info("Figure 15: no violin data available after filtering.")

    log.info(f"Wrote Figure 15 to {out_dir}")
    return 0
