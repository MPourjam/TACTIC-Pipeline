from __future__ import annotations

from fnmatch import fnmatch
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..log import get_logger
from ..utils.labels import display_methods, display_sample_key
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
GRAPH_LEVEL_ORDER = ("species", "tip")
FEATURE_BASIS_ORDER = ("unmapped", "no_hit")


def _is_excluded(cid: str, patterns: list[str]) -> bool:
    for pat in patterns:
        p = str(pat).strip()
        if not p:
            continue
        if fnmatch(cid, p) or (p in cid):
            return True
    return False


def _ordered_levels(df: pd.DataFrame) -> list[str]:
    if "graph_target_level_plot" not in df.columns:
        return ["single"]
    vals = [str(v) for v in df["graph_target_level_plot"].dropna().astype(str).tolist() if str(v).strip()]
    if not vals:
        return ["single"]
    uniq = list(dict.fromkeys(vals))
    ordered: list[str] = [lvl for lvl in GRAPH_LEVEL_ORDER if lvl in uniq]
    ordered += [lvl for lvl in uniq if lvl not in ordered]
    return ordered


def _normalize_effective_identity_threshold(v) -> float | None:
    try:
        f = float(v)
    except Exception:
        return None
    if not np.isfinite(f):
        return None
    if 0.0 <= f <= 1.0:
        return 100.0 * f
    return f


def _fig12_effective_identity_thresholds(cfg) -> dict[str, float]:
    """
    Parse Figure 12 method-specific no-hit thresholds from config.
    Returned keys are canonical METHOD_ORDER labels.
    """
    raw = getattr(cfg.figures, "fig12_no_hit_min_effective_pident_by_method", {}) or {}
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


def _compute_has_best_hit(
    df: pd.DataFrame,
    effective_identity_threshold_by_method: dict[str, float] | None = None,
) -> pd.Series:
    """
    Returns a float Series (1.0 / 0.0) indicating whether each unmapped feature
    has a usable best hit under optional method-specific effective-identity thresholds.
    """
    eff = pd.to_numeric(df.get("best_hit_effective_pident", np.nan), errors="coerce")
    pid = pd.to_numeric(df.get("best_hit_pident", np.nan), errors="coerce")

    # Backward compatibility: legacy exports may store effective identity as [0,1].
    finite_eff = eff[np.isfinite(eff.to_numpy(dtype=float))]
    if finite_eff.size > 0 and float(np.nanmax(finite_eff.to_numpy(dtype=float))) <= 1.0:
        eff = eff * 100.0

    score = eff.copy()
    score = score.where(score.notna(), pid)
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
    Collapse duplicated unmapped-feature rows across graph_target_level_plot into one row.
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


def _build_sample_method_summary(
    level_df: pd.DataFrame,
    ref_df: pd.DataFrame,
    effective_identity_threshold_by_method: dict[str, float] | None = None,
) -> pd.DataFrame:
    ref_pairs = ref_df[["sample_key", "method"]].drop_duplicates().copy() if not ref_df.empty else pd.DataFrame()
    if not ref_pairs.empty:
        ref_pairs["sample_key"] = ref_pairs["sample_key"].astype(str)
        ref_pairs["method"] = ref_pairs["method"].astype(str)
    pair_set = set((str(a), str(b)) for a, b in zip(ref_pairs.get("sample_key", []), ref_pairs.get("method", [])))

    if level_df.empty:
        if ref_pairs.empty:
            return pd.DataFrame(
                columns=[
                    "sample_key",
                    "method",
                    "n_unmapped_features",
                    "n_with_best_hit",
                    "n_no_hit_detected",
                    "no_hit_fraction",
                    "sum_unmapped_feature_fraction",
                    "sum_unmapped_feature_fraction_percent",
                ]
            )
        out = ref_pairs.copy()
        out["n_unmapped_features"] = 0.0
        out["n_with_best_hit"] = 0.0
        out["n_no_hit_detected"] = 0.0
        out["no_hit_fraction"] = 0.0
        out["sum_unmapped_feature_fraction"] = 0.0
        out["sum_unmapped_feature_fraction_percent"] = 0.0
        return out

    dd = level_df.copy()
    dd["sample_key"] = dd["sample_key"].astype(str)
    dd["method"] = dd["method"].astype(str)
    dd["has_best_hit"] = _compute_has_best_hit(
        dd,
        effective_identity_threshold_by_method=effective_identity_threshold_by_method,
    )
    dd["inferred_feature_fraction"] = pd.to_numeric(dd.get("inferred_feature_fraction", np.nan), errors="coerce").fillna(0.0)

    agg = (
        dd.groupby(["sample_key", "method"], as_index=False)
        .agg(
            n_unmapped_features=("unmapped_feature", "count"),
            n_with_best_hit=("has_best_hit", "sum"),
            sum_unmapped_feature_fraction=("inferred_feature_fraction", "sum"),
        )
    )
    agg["n_no_hit_detected"] = agg["n_unmapped_features"] - agg["n_with_best_hit"]
    with np.errstate(divide="ignore", invalid="ignore"):
        agg["no_hit_fraction"] = np.where(
            agg["n_unmapped_features"] > 0,
            agg["n_no_hit_detected"] / agg["n_unmapped_features"],
            0.0,
        )
    agg["sum_unmapped_feature_fraction_percent"] = 100.0 * agg["sum_unmapped_feature_fraction"]

    if not ref_pairs.empty:
        all_pairs = ref_pairs.copy()
        out = all_pairs.merge(agg, on=["sample_key", "method"], how="left")
        fill0_cols = [
            "n_unmapped_features",
            "n_with_best_hit",
            "n_no_hit_detected",
            "no_hit_fraction",
            "sum_unmapped_feature_fraction",
            "sum_unmapped_feature_fraction_percent",
        ]
        for c in fill0_cols:
            out[c] = pd.to_numeric(out[c], errors="coerce").fillna(0.0)
        return out

    return agg


def _draw_no_hit_bar(ax, summary_df: pd.DataFrame, font_size: float, tick_size: float) -> None:
    bars_x: list[int] = []
    bars_y: list[float] = []
    bar_colors: list[str] = []
    annots: list[str] = []
    for j, method in enumerate(METHOD_ORDER, start=1):
        sub = summary_df[summary_df["method"].astype(str) == method].copy()
        n_total = int(pd.to_numeric(sub.get("n_unmapped_features", 0.0), errors="coerce").sum()) if not sub.empty else 0
        n_no = int(pd.to_numeric(sub.get("n_no_hit_detected", 0.0), errors="coerce").sum()) if not sub.empty else 0
        frac = float(n_no / n_total) if n_total > 0 else 0.0
        bars_x.append(j)
        bars_y.append(frac)
        bar_colors.append(METHOD_COLORS.get(method, "#999999"))
        annots.append(f"{n_no}/{n_total}")

    ax.bar(
        bars_x,
        bars_y,
        color=bar_colors,
        alpha=0.8,
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
            fontsize=max(6.0, tick_size - 1.0),
            color="#222222",
        )

    ax.set_title("No-Hit Fraction Among Unmapped Features", fontsize=font_size)
    ax.set_ylabel("No-hit fraction", fontsize=font_size)
    ax.set_ylim(0.0, 1.0)
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


def _draw_unmapped_fraction_heatmap(
    ax,
    summary_df: pd.DataFrame,
    sample_order: list[str],
    feature_basis: str,
    font_size: float,
    tick_size: float,
):
    grid = np.full((len(sample_order), len(METHOD_ORDER)), np.nan, dtype=float)
    val_map = {
        (str(r["sample_key"]), str(r["method"])): float(r["sum_unmapped_feature_fraction_percent"])
        for _, r in summary_df.iterrows()
    }
    for i, sk in enumerate(sample_order):
        for j, m in enumerate(METHOD_ORDER):
            key = (str(sk), str(m))
            if key in val_map:
                grid[i, j] = float(val_map[key])

    cmap = plt.cm.cividis_r.copy()
    cmap.set_bad(color="#eeeeee")
    im = ax.imshow(
        grid,
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        vmin=0.0,
        vmax=100.0,
    )
    if feature_basis == "no_hit":
        ax.set_title("Summed No-Hit Feature Fraction Heatmap", fontsize=font_size)
    else:
        ax.set_title("Summed Unmapped-Feature Fraction Heatmap", fontsize=font_size)
    ax.set_xticks(np.arange(len(METHOD_ORDER)))
    ax.set_xticklabels(
        display_methods(METHOD_ORDER),
        rotation=45,
        ha="right",
        rotation_mode="anchor",
        fontsize=tick_size,
    )
    ax.set_yticks(np.arange(len(sample_order)))
    ax.set_yticklabels([display_sample_key(s) for s in sample_order], fontsize=max(6.0, tick_size - 1.0))
    ax.set_xlabel("Method", fontsize=font_size)
    ax.set_ylabel("Sample", fontsize=font_size)
    return im


def _write_unmapped_best_hit_identity_plot(
    data_df: pd.DataFrame,
    out_path: Path,
    fig_title: str,
    font_size: float,
    tick_size: float,
    point_size: float,
) -> None:
    fig = plt.figure(figsize=(8.2, 4.8))
    ax = fig.add_subplot(1, 1, 1)
    rng = np.random.default_rng(1)

    method_data: list[tuple[int, str, np.ndarray]] = []
    all_vals: list[np.ndarray] = []
    for j, method in enumerate(METHOD_ORDER, start=1):
        sub = data_df[data_df["method"].astype(str) == method].copy()
        vals = pd.to_numeric(sub.get("best_hit_effective_pident", np.nan), errors="coerce").to_numpy(dtype=float)
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
            fontsize=font_size,
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

    ax.set_title("Best-Hit Effective Identity of Unmapped Inferred Features", fontsize=font_size)
    ax.set_ylabel("Best hit effective identity (%)", fontsize=font_size)
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

    if all_vals:
        yvals = np.concatenate(all_vals)
        y_min = float(np.nanmin(yvals))
        y_max = float(np.nanmax(yvals))
        if np.isfinite(y_min) and np.isfinite(y_max):
            pad = 0.5 if y_max <= y_min else max(0.25, 0.08 * (y_max - y_min))
            low = max(0.0, y_min - pad)
            high = min(100.0, y_max + pad)
            if high <= low:
                high = min(100.0, low + 1.0)
            ax.set_ylim(low, high)
        else:
            ax.set_ylim(0.0, 100.0)
    else:
        ax.set_ylim(0.0, 100.0)

    fig.suptitle(fig_title, fontsize=font_size + 2.0)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(out_path, format="svg")
    plt.close(fig)


def run(cfg) -> int:
    """
    Figure 12:
      - Unmapped best-hit identity violin/scatter plot
      - No-hit fraction bar plot (unmapped features with no detectable best hit)
      - Two heatmap versions (sample x method), written by default:
        * basis-unmapped: summed fractions of all unmapped features
        * basis-no_hit: summed fractions of no-hit features only

    Uses Figure 09 exploratory table:
      figures/fig09_split_merge_diagnostics/unmapped_feature_best_hit_identity.tsv
    """
    log = get_logger()
    outdir = Path(cfg.paths.outdir)
    fig09_dir = figures_root(outdir) / "fig09_split_merge_diagnostics"
    src_unmapped = fig09_dir / "unmapped_feature_best_hit_identity.tsv"
    src_diag = fig09_dir / "split_merge_diagnostics.tsv"

    if not src_unmapped.exists():
        log.info("Figure 12: Figure 09 unmapped diagnostics not found; generating Figure 09 prerequisites.")
        plot_fig09.run(cfg)

    if not src_unmapped.exists():
        log.info("Figure 12: unmapped_feature_best_hit_identity.tsv still missing; run plot-fig09 first.")
        return 0

    try:
        unmapped_df = pd.read_csv(src_unmapped, sep="\t")
    except Exception:
        log.info("Figure 12: failed to read Figure 09 unmapped diagnostics TSV.")
        return 0
    if unmapped_df.empty:
        log.info("Figure 12: no unmapped diagnostics rows to plot.")
        return 0

    diag_df = pd.DataFrame()
    if src_diag.exists():
        try:
            diag_df = pd.read_csv(src_diag, sep="\t")
        except Exception:
            diag_df = pd.DataFrame()

    exclude_patterns = list(getattr(cfg.figures, "exclude_community_name_patterns", []) or [])
    exclude_patterns += list(getattr(cfg.figures, "exclude_community_name_patterns_fig12", []) or [])
    if "community_id" in unmapped_df.columns:
        unmapped_df = unmapped_df[
            ~unmapped_df["community_id"].astype(str).map(lambda x: _is_excluded(x, exclude_patterns))
        ].copy()
    if (not diag_df.empty) and ("community_id" in diag_df.columns):
        diag_df = diag_df[
            ~diag_df["community_id"].astype(str).map(lambda x: _is_excluded(x, exclude_patterns))
        ].copy()
    if unmapped_df.empty:
        log.info("Figure 12: all rows excluded or empty after filtering.")
        return 0

    out_dir = figures_root(outdir) / "fig12_unmapped_diagnostics"
    out_dir.mkdir(parents=True, exist_ok=True)
    title = getattr(cfg.figures, "fig12_title", None) or "Figure 12: Unmapped-Feature Diagnostics"
    font_size = max(
        6.0,
        float(
            getattr(
                cfg.figures,
                "fig12_font_size",
                getattr(cfg.figures, "fig09_font_size", 10.0),
            )
        ),
    )
    tick_size = max(6.0, font_size - 2.0)
    point_size = max(1.0, float(getattr(cfg.figures, "box_violin_dot_size", 22.0)))
    no_hit_eff_thresholds = _fig12_effective_identity_thresholds(cfg)

    level_df = _collapse_to_single_feature_level(unmapped_df)
    if level_df.empty:
        log.info("Figure 12: no unmapped rows after collapsing target-level duplicates.")
        return 0
    level_diag = diag_df.copy()
    fig_title = title
    stem_base = "fig12_unmapped_diagnostics"

    if (not level_diag.empty) and {"sample_key", "method"}.issubset(level_diag.columns):
        ref_df = level_diag[["sample_key", "method"]].drop_duplicates().copy()
    else:
        ref_df = level_df[["sample_key", "method"]].drop_duplicates().copy()

    _write_unmapped_best_hit_identity_plot(
        data_df=level_df,
        out_path=out_dir / f"{stem_base}__unmapped_best_hit_identity.svg",
        fig_title=fig_title,
        font_size=font_size,
        tick_size=tick_size,
        point_size=point_size,
    )

    summary_nohit_bar = _build_sample_method_summary(
        level_df,
        ref_df,
        effective_identity_threshold_by_method=no_hit_eff_thresholds,
    )
    if summary_nohit_bar.empty:
        log.info(f"Wrote Figure 12 to {out_dir}")
        return 0

    by_method = (
        summary_nohit_bar.groupby("method", as_index=False)
        .agg(
            n_unmapped_features=("n_unmapped_features", "sum"),
            n_no_hit_detected=("n_no_hit_detected", "sum"),
        )
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        by_method["no_hit_fraction"] = np.where(
            by_method["n_unmapped_features"] > 0,
            by_method["n_no_hit_detected"] / by_method["n_unmapped_features"],
            0.0,
        )

    sample_order = list(dict.fromkeys(ref_df["sample_key"].astype(str).tolist()))
    if not sample_order:
        sample_order = list(dict.fromkeys(summary_nohit_bar["sample_key"].astype(str).tolist()))
    for feature_basis in FEATURE_BASIS_ORDER:
        level_df_heat = level_df.copy()
        level_df_heat["has_best_hit"] = _compute_has_best_hit(
            level_df_heat,
            effective_identity_threshold_by_method=no_hit_eff_thresholds,
        )
        if feature_basis == "no_hit":
            level_df_heat = level_df_heat[level_df_heat["has_best_hit"] <= 0.0].copy()

        summary = _build_sample_method_summary(
            level_df_heat,
            ref_df,
            effective_identity_threshold_by_method=no_hit_eff_thresholds,
        ).copy()
        summary["feature_basis"] = feature_basis
        summary["target_level_plot"] = "single"
        summary["sample_label"] = summary["sample_key"].astype(str).map(display_sample_key)
        summary["method_display"] = summary["method"].astype(str).map(lambda x: display_methods([x])[0])

        stem = f"{stem_base}__basis-{feature_basis}"
        summary.to_csv(out_dir / f"{stem}.tsv", sep="\t", index=False)
        by_method.to_csv(out_dir / f"{stem}__no_hit_by_method.tsv", sep="\t", index=False)

        fig_h = max(8.0, 0.28 * max(1, len(sample_order)) + 4.8)
        fig = plt.figure(figsize=(10.0, fig_h))
        gs = fig.add_gridspec(2, 1, height_ratios=[1.0, max(1.8, 0.12 * len(sample_order))], hspace=0.35)
        ax_bar = fig.add_subplot(gs[0, 0])
        ax_heat = fig.add_subplot(gs[1, 0])

        _draw_no_hit_bar(ax_bar, summary_nohit_bar, font_size=font_size, tick_size=tick_size)
        im = _draw_unmapped_fraction_heatmap(
            ax_heat,
            summary,
            sample_order=sample_order,
            feature_basis=feature_basis,
            font_size=font_size,
            tick_size=tick_size,
        )
        cbar = fig.colorbar(im, ax=ax_heat, fraction=0.03, pad=0.02)
        if feature_basis == "no_hit":
            cbar.set_label("Summed no-hit feature fraction (%)", fontsize=font_size)
            fig.suptitle(f"{fig_title} [No-Hit Basis]", fontsize=font_size + 2.0)
        else:
            cbar.set_label("Summed unmapped feature fraction (%)", fontsize=font_size)
            fig.suptitle(f"{fig_title} [Unmapped Basis]", fontsize=font_size + 2.0)
        cbar.ax.tick_params(labelsize=tick_size)
        cbar.set_ticks([0, 20, 40, 60, 80, 100])

        # Avoid tight_layout warning with mixed Axes + colorbar layouts.
        fig.subplots_adjust(top=0.92)
        fig.savefig(out_dir / f"{stem}.svg", format="svg")
        plt.close(fig)

    log.info(f"Wrote Figure 12 to {out_dir}")
    return 0
