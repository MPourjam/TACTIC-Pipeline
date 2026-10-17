from __future__ import annotations

from fnmatch import fnmatch
from pathlib import Path
import string

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import FuncFormatter, MultipleLocator
import numpy as np
import pandas as pd

from ..log import get_logger
from ..utils.labels import display_method
from ..utils.paths import figures_root
from .plot_fig09 import _compute_panel_stats_permutation, _pvalue_stars

METHOD_ORDER = ["ASV", "zOTU", "OTU", "TIC", "TAC"]
METHOD_FAMILY = {
    "ASV": "denoising",
    "zOTU": "denoising",
    "OTU": "clustering",
    "TIC": "clustering",
    "TAC": "clustering",
}
FAMILY_ORDER = ["denoising", "clustering"]

METHOD_COLORS = {
    "ASV": "#0072B2",
    "zOTU": "#56B4E9",
    "OTU": "#D55E00",
    "TIC": "#009E73",
    "TAC": "#E69F00",
}
FAMILY_COLORS = {
    "denoising": "#0072B2",
    "clustering": "#D55E00",
}
METRIC_LABELS = {
    "jaccard": "Jaccard",
    "weighted_jaccard": "Weighted Jaccard",
    "braycurtis": "Bray-Curtis",
    "unifrac": "UniFrac",
    "gunifrac_a0.0": "gUniFrac 0.0",
    "gunifrac_a0.5": "gUniFrac 0.5",
    "gunifrac_a1.0": "gUniFrac 1.0",
}
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
    "braycurtis": "D",
    "unifrac": "^",
    "gunifrac_a0.0": "v",
    "gunifrac_a0.5": "P",
    "gunifrac_a1.0": "X",
}


def _is_excluded(cid: str, patterns: list[str]) -> bool:
    for pat in patterns:
        p = str(pat).strip()
        if not p:
            continue
        if fnmatch(cid, p) or (p in cid):
            return True
    return False


def _normalize_method_view(raw: object) -> str:
    v = str(raw or "").strip().lower()
    if v in {"family", "families", "collapsed", "groups"}:
        return "families"
    return "methods"


def _normalize_lower_panel(raw: object) -> str:
    v = str(raw or "").strip().lower()
    if v in {"bar", "bars"}:
        return "bars"
    return "density"


def _normalize_lower_panel_y(raw: object) -> str:
    v = str(raw or "").strip().lower()
    if v in {"method", "methods"}:
        return "methods"
    return "metrics"


def _normalize_family_aggregation(raw: object) -> str:
    v = str(raw or "").strip().lower()
    if v in {"mean", "min"}:
        return v
    return "median"


def _normalize_stats_annotation_mode(raw: object) -> str:
    v = str(raw or "").strip().lower()
    if v in {"brackets", "stars", "star", "lines"}:
        return "brackets"
    return "letters"


def _family_label(family: str) -> str:
    s = str(family or "").strip().lower()
    if s == "denoising":
        return "Denoising"
    if s == "clustering":
        return "Clustering"
    return s.title()


def _metric_label(metric: str) -> str:
    return METRIC_LABELS.get(str(metric), str(metric))


def _format_xtick(x, _pos) -> str:
    if abs(x) < 1e-12:
        return "0.0"
    return f"{x:.1f}"


def _kde(values: np.ndarray, grid: np.ndarray, bandwidth: float | None) -> np.ndarray:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return np.zeros_like(grid, dtype=float)
    if vals.size == 1:
        bw = float(bandwidth) if bandwidth not in (None, 0) else 0.02
    else:
        std = float(np.std(vals, ddof=1))
        if not np.isfinite(std) or std <= 0.0:
            std = 0.02
        auto_bw = 1.06 * std * (vals.size ** (-1.0 / 5.0))
        bw = float(bandwidth) if bandwidth not in (None, 0) else auto_bw
        if not np.isfinite(bw) or bw <= 0.0:
            bw = max(0.01, auto_bw, std * 0.2)
    bw = max(0.005, float(bw))
    z = (grid[:, None] - vals[None, :]) / bw
    dens = np.exp(-0.5 * (z**2)).sum(axis=1)
    dens /= (vals.size * bw * np.sqrt(2.0 * np.pi))
    return dens


def _load_distance_data(cfg) -> pd.DataFrame:
    path = cfg.paths.outdir / "distances_all_long.tsv"
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return pd.DataFrame()
    df["community_id"] = df["community_id"].astype(str)
    df["sample_id"] = df["sample_id"].astype(str)
    df["method"] = df["method"].astype(str)
    df["metric"] = df["metric"].astype(str)
    df["distance"] = pd.to_numeric(df["distance"], errors="coerce")
    df = df[np.isfinite(df["distance"])].copy()
    df["sample_key"] = df["community_id"] + "|" + df["sample_id"]
    df["method_label"] = df["method"].map(display_method)
    df["method_family"] = df["method"].map(METHOD_FAMILY).fillna("other")
    df["metric_label"] = df["metric"].map(_metric_label)
    return df


def _prepare_plot_data(cfg, df: pd.DataFrame) -> tuple[pd.DataFrame, list[str], list[str], str, str]:
    allowed_metrics = set(METRIC_LABELS.keys())
    available_metrics = set(df["metric"].astype(str).tolist())
    metrics_cfg = [str(m) for m in (getattr(cfg.figures, "fig17_metrics", []) or []) if str(m).strip()]
    if metrics_cfg:
        metric_order = [m for m in metrics_cfg if m in available_metrics and m in allowed_metrics]
    else:
        metric_order = [m for m in getattr(cfg.distances, "metrics", []) if m in available_metrics and m in allowed_metrics]
    if not metric_order:
        metric_order = [m for m in df["metric"].drop_duplicates().astype(str).tolist() if m in allowed_metrics]
    method_candidates = [m for m in METHOD_ORDER if m in set(df["method"].astype(str).tolist())]
    if not method_candidates:
        method_candidates = [m for m in df["method"].drop_duplicates().astype(str).tolist()]

    exclude_patterns = list(getattr(cfg.figures, "exclude_community_name_patterns", []) or [])
    exclude_patterns += list(getattr(cfg.figures, "exclude_community_name_patterns_fig17", []) or [])
    if exclude_patterns:
        df = df[~df["community_id"].map(lambda cid: _is_excluded(cid, exclude_patterns))].copy()

    df = df[df["metric"].isin(metric_order) & df["method"].isin(method_candidates)].copy()

    method_view = _normalize_method_view(getattr(cfg.figures, "fig17_method_view", "methods"))
    family_agg = _normalize_family_aggregation(
        getattr(cfg.figures, "fig17_stats_family_aggregation", "median")
    )
    if method_view == "families":
        agg_name = "median" if family_agg == "median" else family_agg
        grouped = (
            df.groupby(
                ["community_id", "sample_id", "sample_key", "metric", "method_family", "metric_label"],
                as_index=False,
            )["distance"]
            .agg(agg_name)
            .rename(columns={"method_family": "method"})
        )
        grouped["method_label"] = grouped["method"].map(_family_label)
        grouped["method_family"] = grouped["method"]
        df = grouped
        row_order = [f for f in FAMILY_ORDER if f in set(df["method"].astype(str).tolist())]
    else:
        row_order = method_candidates

    return df, metric_order, row_order, method_view, family_agg


def _facet_structure(metric_order: list[str], row_order: list[str], lower_panel_y: str) -> tuple[list[str], list[str], str, str]:
    if lower_panel_y == "methods":
        return row_order, metric_order, "method", "metric"
    return metric_order, row_order, "metric", "method"


def _point_color(method_view: str, row_key: str, metric_key: str, lower_panel_y: str) -> str:
    if lower_panel_y == "methods":
        return METRIC_COLORS.get(metric_key, "#4D4D4D")
    if method_view == "families":
        return FAMILY_COLORS.get(row_key, "#4D4D4D")
    return METHOD_COLORS.get(row_key, "#4D4D4D")


def _row_color(method_view: str, row_key: str, lower_panel_y: str) -> str:
    if lower_panel_y == "methods":
        return METRIC_COLORS.get(row_key, "#4D4D4D")
    if method_view == "families":
        return FAMILY_COLORS.get(row_key, "#4D4D4D")
    return METHOD_COLORS.get(row_key, "#4D4D4D")


def _display_row_label(method_view: str, lower_panel_y: str, row_key: str) -> str:
    if lower_panel_y == "methods":
        return _metric_label(row_key)
    if method_view == "families":
        return _family_label(row_key)
    return display_method(row_key)


def _display_facet_label(method_view: str, lower_panel_y: str, facet_key: str) -> str:
    if lower_panel_y == "methods":
        if method_view == "families":
            return _family_label(facet_key)
        return display_method(facet_key)
    return _metric_label(facet_key)


def _write_points_table(df: pd.DataFrame, out_path: Path, method_view: str) -> None:
    pts = df.copy()
    if method_view == "families":
        pts["method_label"] = pts["method"].map(_family_label)
    else:
        pts["method_label"] = pts["method"].map(display_method)
    pts["method_family"] = pts["method_family"].astype(str)
    pts["metric_label"] = pts["metric"].map(_metric_label)
    cols = [
        "community_id",
        "sample_id",
        "method",
        "metric",
        "distance",
        "sample_key",
        "method_label",
        "method_family",
        "metric_label",
    ]
    pts[cols].to_csv(out_path, sep="\t", index=False)


def _compute_stats_tables(
    df: pd.DataFrame,
    facets: list[str],
    row_order: list[str],
    facet_axis: str,
    row_axis: str,
    lower_panel_y: str,
    cfg,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, list[dict]]]:
    alpha = float(getattr(cfg.figures, "fig17_stats_alpha", 0.05))
    pairwise_test = str(getattr(cfg.figures, "fig17_stats_pairwise_test", "signed_rank_permutation"))
    multiple_corr = bool(getattr(cfg.figures, "fig17_stats_multiple_test_correction", True))
    n_perm = int(getattr(cfg.figures, "fig11_stats_permutations", 4000))
    seed = int(getattr(cfg.figures, "fig11_stats_seed", 1))

    pairwise_rows: list[dict] = []
    summary_rows: list[dict] = []
    sig_map: dict[str, list[dict]] = {}

    for facet in facets:
        sub = df[df[facet_axis] == facet].copy()
        work = sub.rename(columns={row_axis: "method", "distance": "distance_value"})
        stats = _compute_panel_stats_permutation(
            df=work[["sample_key", "method", "distance_value"]],
            metric_col="distance_value",
            method_order=row_order,
            n_perm=n_perm,
            seed=seed,
            alpha=alpha,
            pairwise_test=pairwise_test,
            multiple_test_correction=multiple_corr,
        )
        note = "ok"
        if int(stats.get("n_blocks", 0)) < 2 or len(stats.get("methods", [])) < 2:
            note = "insufficient_data"
        summary_rows.append(
            {
                "facet_axis": facet_axis,
                "facet_value": facet,
                "row_axis": row_axis,
                "n_rows_present": len(stats.get("methods", [])),
                "n_blocks": int(stats.get("n_blocks", 0)),
                "friedman_q": float(stats.get("friedman_q", np.nan)),
                "friedman_p": float(stats.get("friedman_p", np.nan)),
                "stats_model": str(stats.get("stats_model", "paired_permutation")),
                "pairwise_test": pairwise_test,
                "multiple_test_correction": bool(multiple_corr),
                "alpha": alpha,
                "note": note,
            }
        )
        sig_rows = []
        for row in stats.get("pairwise", []):
            rec = {
                "facet_axis": facet_axis,
                "facet_value": facet,
                "row_axis": row_axis,
                "row_a": row.get("method_a"),
                "row_b": row.get("method_b"),
                "pairwise_test": row.get("pairwise_test"),
                "statistic": row.get("statistic"),
                "p_raw": row.get("p_raw"),
                "p_adj_bh": row.get("p_adj_bh"),
                "p_value_used": row.get("p_value_used"),
                "multiple_test_correction": row.get("multiple_test_correction"),
                "significant_alpha": row.get("significant_alpha"),
                "n_nonzero_pairs": row.get("n_nonzero_pairs"),
                "median_diff_a_minus_b": row.get("median_diff_a_minus_b"),
                "alpha": alpha,
            }
            pairwise_rows.append(rec)
            if bool(row.get("significant_alpha")):
                sig_rows.append(rec)
        sig_map[str(facet)] = sig_rows

    pairwise_df = pd.DataFrame(pairwise_rows)
    summary_df = pd.DataFrame(summary_rows)
    # Preserve historical labels in TSVs.
    if not pairwise_df.empty:
        pairwise_df["facet_axis"] = np.where(
            str(lower_panel_y) == "methods", "method", "metric"
        )
        pairwise_df["row_axis"] = np.where(
            str(lower_panel_y) == "methods", "metric", "method"
        )
    if not summary_df.empty:
        summary_df["facet_axis"] = np.where(
            str(lower_panel_y) == "methods", "method", "metric"
        )
        summary_df["row_axis"] = np.where(
            str(lower_panel_y) == "methods", "metric", "method"
        )
    return pairwise_df, summary_df, sig_map


def _draw_upper_panel(
    ax,
    df: pd.DataFrame,
    row_order: list[str],
    metric_order: list[str],
    method_view: str,
    font_size: float,
    dot_size: float,
    xtick_interval: float,
) -> None:
    row_to_y = {row: (len(row_order) - 1 - i) for i, row in enumerate(row_order)}
    for metric in metric_order:
        subm = df[df["metric"] == metric]
        for row in row_order:
            sub = subm[subm["method"] == row]
            if sub.empty:
                continue
            color = _row_color(method_view, row, "metrics")
            is_hollow = (row == "ASV") or (row == "zOTU") or (row == "denoising")
            marker = METRIC_MARKERS.get(metric, "o")
            ax.scatter(
                sub["distance"].to_numpy(dtype=float),
                np.full(sub.shape[0], row_to_y[row], dtype=float),
                s=dot_size,
                marker=marker,
                facecolors="none" if is_hollow else color,
                edgecolors=color,
                linewidths=1.2,
                alpha=0.95,
                zorder=3,
            )
    ax.set_yticks([row_to_y[r] for r in row_order], [_display_row_label(method_view, "metrics", r) for r in row_order])
    ax.set_ylabel("Method", fontsize=font_size)
    ax.set_ylim(-0.75, len(row_order) - 0.25)
    ax.xaxis.set_major_locator(MultipleLocator(max(0.01, float(xtick_interval))))
    ax.xaxis.set_major_formatter(FuncFormatter(_format_xtick))
    ax.grid(axis="x", color="#D9D9D9", linewidth=0.8, alpha=0.7)
    ax.tick_params(axis="both", labelsize=max(font_size - 1.0, 8.0))
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def _draw_lower_bars(ax, facet_df: pd.DataFrame, row_order: list[str], method_view: str, lower_panel_y: str, font_size: float) -> None:
    row_to_y = {row: (len(row_order) - 1 - i) for i, row in enumerate(row_order)}
    for row in row_order:
        sub = facet_df[facet_df["row_key"] == row]
        if sub.empty:
            continue
        color = _row_color(method_view, row, lower_panel_y)
        y0 = row_to_y[row]
        for x in sub["distance"].to_numpy(dtype=float):
            ax.plot([x, x], [y0 - 0.33, y0 + 0.33], color=color, linewidth=1.6, alpha=0.95, solid_capstyle="round")
    ax.set_yticks([row_to_y[r] for r in row_order], [_display_row_label(method_view, lower_panel_y, r) for r in row_order])
    ax.set_ylim(-0.75, len(row_order) - 0.25)
    ax.tick_params(axis="both", labelsize=max(font_size - 1.0, 8.0))


def _draw_density_panel(
    ax,
    facet_df: pd.DataFrame,
    row_order: list[str],
    method_view: str,
    lower_panel_y: str,
    font_size: float,
    bandwidth: float | None,
    density_points: int,
    density_scale: float,
    opacity: float,
) -> None:
    row_to_y = {row: (len(row_order) - 1 - i) for i, row in enumerate(row_order)}
    grid = np.linspace(0.0, 1.0, max(64, int(density_points)))
    scale = max(0.15, float(density_scale))
    peak_label_offset = max(0.008, 0.03 * scale)
    for row in row_order:
        sub = facet_df[facet_df["row_key"] == row]
        vals = sub["distance"].to_numpy(dtype=float)
        if vals.size == 0:
            continue
        color = _row_color(method_view, row, lower_panel_y)
        y0 = row_to_y[row]
        dens = _kde(vals, grid, bandwidth)
        if np.nanmax(dens) > 0:
            dens = dens / float(np.nanmax(dens))
        height = dens * scale
        ax.fill_between(grid, y0, y0 + height, color=color, alpha=max(0.0, min(1.0, opacity)), linewidth=0.0, zorder=1)
        ax.plot(grid, y0 + height, color=color, linewidth=1.8, alpha=0.98, zorder=2)
        median_x = float(np.median(vals))
        ax.plot([median_x, median_x], [y0, y0 + max(0.12, float(np.nanmax(height)))], color="#111111", linewidth=1.0, alpha=0.9, zorder=3)

        # Mark local density modes.
        if dens.size >= 3:
            peaks = np.where((dens[1:-1] > dens[:-2]) & (dens[1:-1] >= dens[2:]))[0] + 1
            if peaks.size == 0:
                peaks = np.array([int(np.argmax(dens))], dtype=int)
            peak_heights = dens[peaks]
            order = peaks[np.argsort(peak_heights)[::-1]]
            keep = []
            for idx in order:
                if dens[idx] < 0.15:
                    continue
                if any(abs(grid[idx] - grid[j]) < 0.04 for j in keep):
                    continue
                keep.append(int(idx))
            if not keep:
                keep = [int(order[0])]
            for idx in sorted(keep):
                px = float(grid[idx])
                py = float(y0 + height[idx])
                ax.scatter([px], [py], s=18.0, color="#111111", zorder=4)
                ax.text(
                    px,
                    py + peak_label_offset,
                    f"{px:.2f}",
                    ha="center",
                    va="bottom",
                    fontsize=max(font_size - 4.0, 6.5),
                    color="#111111",
                    zorder=5,
                )
    ax.set_yticks([row_to_y[r] for r in row_order], [_display_row_label(method_view, lower_panel_y, r) for r in row_order])
    ax.set_ylim(-0.75, len(row_order) - 0.05 + scale)
    ax.tick_params(axis="both", labelsize=max(font_size - 1.0, 8.0))


def _draw_significance_brackets(ax, sig_rows: list[dict], row_order: list[str], lower_panel_y: str, alpha: float, font_size: float) -> None:
    if not sig_rows:
        return
    row_to_y = {row: (len(row_order) - 1 - i) for i, row in enumerate(row_order)}
    intervals = []
    for row in sig_rows:
        a = str(row.get("row_a", ""))
        b = str(row.get("row_b", ""))
        if a not in row_to_y or b not in row_to_y:
            continue
        y1 = row_to_y[a]
        y2 = row_to_y[b]
        if y1 == y2:
            continue
        lo, hi = sorted((y1, y2))
        intervals.append((lo, hi, row))
    if not intervals:
        return
    intervals.sort(key=lambda x: (x[1] - x[0], x[0]))
    levels: list[list[tuple[float, float]]] = []
    placed: list[tuple[float, float, int, dict]] = []
    for lo, hi, row in intervals:
        level = None
        for i, spans in enumerate(levels):
            if all((hi < a) or (lo > b) for a, b in spans):
                spans.append((lo, hi))
                level = i
                break
        if level is None:
            levels.append([(lo, hi)])
            level = len(levels) - 1
        placed.append((lo, hi, level, row))

    x_base = 1.03
    x_step = 0.035
    cap = 0.012
    for lo, hi, level, row in placed:
        x = x_base + (level * x_step)
        ax.plot([x, x], [lo, hi], color="#111111", linewidth=1.2, clip_on=False, zorder=6)
        ax.plot([x - cap, x], [lo, lo], color="#111111", linewidth=1.2, clip_on=False, zorder=6)
        ax.plot([x - cap, x], [hi, hi], color="#111111", linewidth=1.2, clip_on=False, zorder=6)
        p = float(row.get("p_value_used", np.nan))
        ax.text(
            x + 0.006,
            (lo + hi) / 2.0,
            _pvalue_stars(p, alpha=alpha),
            ha="left",
            va="center",
            fontsize=max(font_size - 2.0, 8.0),
            fontweight="bold",
            color="#111111",
            clip_on=False,
            zorder=7,
        )


def _letter_code(idx: int) -> str:
    letters = string.ascii_lowercase
    if idx < 0:
        return "a"
    out = ""
    n = int(idx)
    while True:
        out = letters[n % 26] + out
        n = (n // 26) - 1
        if n < 0:
            break
    return out


def _compute_compact_letter_display(
    row_order: list[str],
    row_present: set[str],
    pairwise_rows: list[dict],
    alpha: float,
) -> dict[str, str]:
    rows = [r for r in row_order if (not row_present or r in row_present)]
    if not rows:
        return {}

    sig_pairs: set[frozenset[str]] = set()
    for rec in pairwise_rows:
        a = str(rec.get("row_a", ""))
        b = str(rec.get("row_b", ""))
        if (a == b) or (a not in rows) or (b not in rows):
            continue
        sig = rec.get("significant_alpha")
        if sig is None:
            p = pd.to_numeric(rec.get("p_value_used"), errors="coerce")
            sig = bool(np.isfinite(p) and float(p) < float(alpha))
        if bool(sig):
            sig_pairs.add(frozenset((a, b)))

    letter_groups: list[set[str]] = []
    for row in rows:
        placed = False
        for grp in letter_groups:
            if all(frozenset((row, other)) not in sig_pairs for other in grp):
                grp.add(row)
                placed = True
        if not placed:
            letter_groups.append({row})

    for i, a in enumerate(rows):
        for b in rows[i + 1 :]:
            pair = frozenset((a, b))
            if pair in sig_pairs:
                continue
            if any((a in grp and b in grp) for grp in letter_groups):
                continue
            grp = {a, b}
            for c in rows:
                if c in grp:
                    continue
                if all(frozenset((c, x)) not in sig_pairs for x in grp):
                    grp.add(c)
            letter_groups.append(grp)

    row_letters: dict[str, list[str]] = {r: [] for r in rows}
    for idx, grp in enumerate(letter_groups):
        lbl = _letter_code(idx)
        for r in rows:
            if (r in grp) and (lbl not in row_letters[r]):
                row_letters[r].append(lbl)

    out: dict[str, str] = {}
    for r in rows:
        codes = row_letters.get(r, [])
        out[r] = "".join(codes) if codes else "a"
    return out


def _draw_significance_letters(
    ax,
    pairwise_rows: list[dict],
    row_order: list[str],
    row_present: set[str],
    lower_panel_y: str,
    alpha: float,
    font_size: float,
) -> None:
    labels = _compute_compact_letter_display(
        row_order=row_order,
        row_present=row_present,
        pairwise_rows=pairwise_rows,
        alpha=alpha,
    )
    if not labels:
        return
    row_to_y = {row: (len(row_order) - 1 - i) for i, row in enumerate(row_order)}
    x = 1.045
    for row in row_order:
        txt = labels.get(row)
        if not txt:
            continue
        y = row_to_y.get(row)
        if y is None:
            continue
        ax.text(
            x,
            y,
            txt,
            ha="left",
            va="center",
            fontsize=max(font_size - 1.8, 8.0),
            fontweight="bold",
            color="#111111",
            clip_on=False,
            zorder=7,
        )


def run(cfg) -> int:
    log = get_logger()
    if not bool(getattr(cfg.figures, "fig17", True)):
        log.info("Figure 17 disabled in config; skipping.")
        return 0

    df = _load_distance_data(cfg)
    if df.empty:
        log.info("Figure 17 skipped: distances_all_long.tsv missing or empty.")
        return 0

    df, metric_order, row_order, method_view, family_agg = _prepare_plot_data(cfg, df)
    if df.empty or not metric_order or not row_order:
        log.info("Figure 17 skipped: no rows remain after filtering.")
        return 0

    lower_panel = _normalize_lower_panel(getattr(cfg.figures, "fig17_lower_panel", "density"))
    lower_panel_y = _normalize_lower_panel_y(getattr(cfg.figures, "fig17_lower_panel_y", "metrics"))
    show_upper = bool(getattr(cfg.figures, "fig17_show_upper_panel", True))
    font_size = float(getattr(cfg.figures, "fig17_font_size", 10.0))
    dot_size = float(getattr(cfg.figures, "fig17_dot_size", 32.0))
    opacity = float(getattr(cfg.figures, "fig17_opacity", 1.0))
    lower_height_ratio = float(getattr(cfg.figures, "fig17_lower_panel_height_ratio", 6.0))
    density_bandwidth = getattr(cfg.figures, "fig17_density_bandwidth", None)
    density_points = int(getattr(cfg.figures, "fig17_density_points", 256))
    density_scale = float(getattr(cfg.figures, "fig17_density_scale", 0.72))
    xtick_interval = float(getattr(cfg.figures, "fig17_xtick_interval", 0.2))
    if xtick_interval <= 0:
        xtick_interval = 0.2
    metric_title_pad = float(getattr(cfg.figures, "fig17_metric_title_pad", 14.0))
    grid_mode = str(getattr(cfg.figures, "fig17_grid_mode", "method_pipeline") or "method_pipeline").strip().lower()
    method_pipeline_grid = grid_mode in {"method_pipeline", "method*pipeline", "method_metric", "matrix"}
    stats_annotation_mode = _normalize_stats_annotation_mode(
        getattr(cfg.figures, "fig17_stats_annotation_mode", "letters")
    )

    figdir = figures_root(cfg.paths.outdir) / "fig17_dot_dash_distances"
    figdir.mkdir(parents=True, exist_ok=True)

    _write_points_table(df, figdir / "dot_dash_distance_points.tsv", method_view=method_view)

    facets, lower_row_order, facet_axis, row_axis = _facet_structure(metric_order, row_order, lower_panel_y)
    lower_df = df.copy().rename(
        columns={
            facet_axis: "facet_key",
            row_axis: "row_key",
        }
    )
    pairwise_df, summary_df, sig_map = _compute_stats_tables(
        df=df,
        facets=facets,
        row_order=lower_row_order,
        facet_axis=facet_axis,
        row_axis=row_axis,
        lower_panel_y=lower_panel_y,
        cfg=cfg,
    )
    pairwise_df.to_csv(figdir / "fig17_stats_pairwise.tsv", sep="\t", index=False)
    summary_df.to_csv(figdir / "fig17_stats_facet_summary.tsv", sep="\t", index=False)
    pairwise_map: dict[str, list[dict]] = {}
    if not pairwise_df.empty:
        for fv, sub in pairwise_df.groupby("facet_value", dropna=False):
            pairwise_map[str(fv)] = sub.to_dict("records")
    cld_map: dict[str, dict[str, str]] = {}
    if bool(getattr(cfg.figures, "fig17_stats_enabled", True)) and stats_annotation_mode == "letters":
        alpha = float(getattr(cfg.figures, "fig17_stats_alpha", 0.05))
        if method_pipeline_grid:
            if lower_panel_y == "metrics":
                facet_keys = list(metric_order)
                row_order_for_cld = list(row_order)
            else:
                facet_keys = list(row_order)
                row_order_for_cld = list(metric_order)
            for fk in facet_keys:
                cld_map[str(fk)] = _compute_compact_letter_display(
                    row_order=row_order_for_cld,
                    row_present=set(row_order_for_cld),
                    pairwise_rows=pairwise_map.get(str(fk), []),
                    alpha=alpha,
                )

    if method_pipeline_grid:
        panel_methods = list(row_order)
        panel_metrics = list(metric_order)
        n_rows = max(1, len(panel_methods))
        n_cols = max(1, len(panel_metrics))
        n_facets = n_rows * n_cols
        panel_size = 2.25
        base_width = max(10.5, panel_size * n_cols + 3.6)
        lower_rows = 1
    else:
        n_facets = max(1, len(facets))
        n_cols = int(np.ceil(np.sqrt(n_facets)))
        n_rows = int(np.ceil(n_facets / float(n_cols)))
        base_width = max(12.0, 2.9 * n_cols + 4.0)
        lower_rows = max(2, len(lower_row_order))
    if show_upper:
        if method_pipeline_grid:
            lower_panel_block_h = max(2.8, panel_size)
        else:
            lower_panel_block_h = max(3.0, 0.78 * lower_rows + (0.55 * lower_height_ratio))
        fig_height = 2.4 + (n_rows * lower_panel_block_h)
        fig = plt.figure(figsize=(base_width, fig_height))
        gs = GridSpec(2, 1, height_ratios=[2.0, max(2.0, lower_height_ratio * n_rows)], hspace=0.18, figure=fig)
        ax_upper = fig.add_subplot(gs[0, 0])
        lower_spec = gs[1, 0].subgridspec(
            n_rows,
            n_cols,
            wspace=0.08 if method_pipeline_grid else 0.18,
            hspace=0.12 if method_pipeline_grid else 0.24,
        )
    else:
        if method_pipeline_grid:
            fig_height = max(4.2, panel_size * n_rows + 1.8)
        else:
            fig_height = max(3.8, n_rows * (0.95 * lower_rows + (0.65 * lower_height_ratio)))
        fig = plt.figure(figsize=(base_width, fig_height))
        ax_upper = None
        lower_spec = GridSpec(
            n_rows,
            n_cols,
            left=0.08,
            right=0.86,
            bottom=0.12,
            top=0.88,
            wspace=0.08 if method_pipeline_grid else 0.18,
            hspace=0.12 if method_pipeline_grid else 0.24,
            figure=fig,
        )

    if show_upper and ax_upper is not None:
        _draw_upper_panel(
            ax_upper,
            df,
            row_order=row_order,
            metric_order=metric_order,
            method_view=method_view,
            font_size=font_size,
            dot_size=dot_size,
            xtick_interval=xtick_interval,
        )
        ax_upper.set_xlabel("Distance", fontsize=font_size)
        ax_upper.set_xlim(0.0, 1.0)
        plt.setp(ax_upper.get_xticklabels(), rotation=45, ha="right")

    lower_axes: list[tuple[plt.Axes, int, int]] = []
    if method_pipeline_grid:
        for r, method_key in enumerate(panel_methods):
            for c, metric_key in enumerate(panel_metrics):
                ax = fig.add_subplot(lower_spec[r, c] if show_upper else lower_spec[r, c])
                lower_axes.append((ax, r, c))
                ax.set_box_aspect(1.0)
                facet_df = df[(df["method"] == method_key) & (df["metric"] == metric_key)].copy()
                facet_df["row_key"] = str(method_key)
                row_order_cell = [str(method_key)]
                if lower_panel == "bars":
                    _draw_lower_bars(
                        ax=ax,
                        facet_df=facet_df,
                        row_order=row_order_cell,
                        method_view=method_view,
                        lower_panel_y="metrics",
                        font_size=font_size,
                    )
                else:
                    _draw_density_panel(
                        ax=ax,
                        facet_df=facet_df,
                        row_order=row_order_cell,
                        method_view=method_view,
                        lower_panel_y="metrics",
                        font_size=font_size,
                        bandwidth=density_bandwidth,
                        density_points=density_points,
                        density_scale=density_scale,
                        opacity=opacity,
                    )
                if r == 0:
                    ax.set_title(_metric_label(metric_key), fontsize=max(font_size - 0.5, 9.0), pad=metric_title_pad)
                else:
                    ax.set_title("")
                ax.set_xlim(0.0, 1.14)
                ax.xaxis.set_major_locator(MultipleLocator(max(0.01, xtick_interval)))
                ax.xaxis.set_major_formatter(FuncFormatter(_format_xtick))
                ax.grid(axis="x", color="#D9D9D9", linewidth=0.8, alpha=0.7)
                if c == 0:
                    ax.set_ylabel(_display_row_label(method_view, "metrics", method_key), fontsize=font_size)
                else:
                    ax.set_ylabel("")
                    ax.set_yticklabels([])
                    ax.tick_params(axis="y", length=0)
                ax.set_yticks([])
                if lower_panel == "bars":
                    ax.set_ylim(-0.36, 0.36)
                else:
                    # Tight single-row density panels: keep baseline close to x-axis
                    # and reduce empty headroom above the curve.
                    upper = max(0.22, float(density_scale) * 1.06)
                    ax.set_ylim(0.0, upper)
                if bool(getattr(cfg.figures, "fig17_stats_enabled", True)) and stats_annotation_mode == "letters":
                    facet_key_for_cld = metric_key if lower_panel_y == "metrics" else method_key
                    row_key_for_cld = method_key if lower_panel_y == "metrics" else metric_key
                    lbl = cld_map.get(str(facet_key_for_cld), {}).get(str(row_key_for_cld), "")
                    if lbl:
                        ax.text(
                            0.97,
                            0.93,
                            lbl,
                            transform=ax.transAxes,
                            ha="right",
                            va="top",
                            fontsize=max(font_size - 1.8, 8.0),
                            fontweight="bold",
                            color="#111111",
                            zorder=8,
                        )
                plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
                for spine in ("top", "right"):
                    ax.spines[spine].set_visible(False)
    else:
        for idx, facet in enumerate(facets):
            r = idx // n_cols
            c = idx % n_cols
            ax = fig.add_subplot(lower_spec[r, c] if show_upper else lower_spec[r, c])
            lower_axes.append((ax, r, c))
            ax.set_box_aspect(1.0)
            facet_df = lower_df[lower_df["facet_key"] == facet].copy()
            if lower_panel == "bars":
                _draw_lower_bars(
                    ax=ax,
                    facet_df=facet_df,
                    row_order=lower_row_order,
                    method_view=method_view,
                    lower_panel_y=lower_panel_y,
                    font_size=font_size,
                )
            else:
                _draw_density_panel(
                    ax=ax,
                    facet_df=facet_df,
                    row_order=lower_row_order,
                    method_view=method_view,
                    lower_panel_y=lower_panel_y,
                    font_size=font_size,
                    bandwidth=density_bandwidth,
                    density_points=density_points,
                    density_scale=density_scale,
                    opacity=opacity,
                )
            ax.set_title(
                _display_facet_label(method_view, lower_panel_y, facet),
                fontsize=max(font_size - 0.5, 9.0),
                pad=metric_title_pad,
            )
            ax.set_xlim(0.0, 1.14)
            ax.xaxis.set_major_locator(MultipleLocator(max(0.01, xtick_interval)))
            ax.xaxis.set_major_formatter(FuncFormatter(_format_xtick))
            ax.grid(axis="x", color="#D9D9D9", linewidth=0.8, alpha=0.7)
            if c == 0:
                ax.set_ylabel("Metric" if lower_panel_y == "methods" else "Method", fontsize=font_size)
            else:
                ax.set_ylabel("")
                ax.set_yticklabels([])
            if c != 0:
                ax.tick_params(axis="y", length=0)
            plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
            if bool(getattr(cfg.figures, "fig17_stats_enabled", True)):
                if stats_annotation_mode == "letters":
                    _draw_significance_letters(
                        ax=ax,
                        pairwise_rows=pairwise_map.get(str(facet), []),
                        row_order=lower_row_order,
                        row_present=set(facet_df["row_key"].astype(str).tolist()),
                        lower_panel_y=lower_panel_y,
                        alpha=float(getattr(cfg.figures, "fig17_stats_alpha", 0.05)),
                        font_size=font_size,
                    )
                else:
                    _draw_significance_brackets(
                        ax=ax,
                        sig_rows=sig_map.get(str(facet), []),
                        row_order=lower_row_order,
                        lower_panel_y=lower_panel_y,
                        alpha=float(getattr(cfg.figures, "fig17_stats_alpha", 0.05)),
                        font_size=font_size,
                    )
            for spine in ("top", "right"):
                ax.spines[spine].set_visible(False)

    for ax, r, _c in lower_axes:
        if method_pipeline_grid and r != (n_rows - 1):
            ax.set_xlabel("")
            ax.set_xticklabels([])
        else:
            ax.set_xlabel("Distance", fontsize=font_size)

    title = str(getattr(cfg.figures, "fig17_title", "") or "").strip()
    if title:
        fig.suptitle(title, fontsize=font_size + 1.0, y=0.98)

    if bool(getattr(cfg.figures, "fig17_stats_enabled", True)) and bool(getattr(cfg.figures, "fig17_stats_show_legend", True)) and (
        (not method_pipeline_grid) or (stats_annotation_mode == "letters")
    ):
        alpha = float(getattr(cfg.figures, "fig17_stats_alpha", 0.05))
        use_q = bool(getattr(cfg.figures, "fig17_stats_multiple_test_correction", True))
        pv = "q" if use_q else "p"
        if stats_annotation_mode == "letters":
            legend_text = "\n".join(
                [
                    "CLD letters",
                    "shared letter: not different",
                    f"no shared letter: {pv} < {alpha:.2f}",
                ]
            )
        else:
            legend_text = "\n".join(
                [
                    f"*** {pv} < 0.001",
                    f"** {pv} < 0.01",
                    f"* {pv} < {alpha:.2f}",
                ]
            )
        fig.text(
            0.885,
            0.94,
            legend_text,
            ha="left",
            va="top",
            fontsize=max(font_size - 2.0, 8.0),
            bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#BFBFBF"},
        )

    fig.subplots_adjust(left=0.11 if method_pipeline_grid else 0.09, right=0.86, bottom=0.10, top=0.90 if not title else 0.92)

    stem = "distance_dot_dash"
    svg_path = figdir / f"{stem}.svg"
    png_path = figdir / f"{stem}.png"
    fig.savefig(svg_path, dpi=300, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    if method_view == "families":
        fig.savefig(figdir / "distance_dot_dash_denoising_vs_clustering.svg", dpi=300, bbox_inches="tight")
    plt.close(fig)

    log.info(f"Figure 17 written to {svg_path}")
    return 0
