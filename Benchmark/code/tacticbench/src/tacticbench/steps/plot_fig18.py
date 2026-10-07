from __future__ import annotations

from fnmatch import fnmatch
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..log import get_logger
from ..utils.labels import display_method
from ..utils.paths import figures_root

METHOD_ORDER = ["ASV", "zOTU", "OTU", "TIC", "TAC"]
METHOD_COLORS = {
    "ASV": "#0072B2",
    "zOTU": "#56B4E9",
    "OTU": "#D55E00",
    "TIC": "#009E73",
    "TAC": "#E69F00",
}


def _normalize_plot_type(raw: str, default: str = "violin") -> str:
    v = str(raw or "").strip().lower()
    if v in {"box", "boxplot"}:
        return "boxplot"
    if v in {"violin", "violinplot"}:
        return "violin"
    return default


def _is_excluded(cid: str, patterns: list[str]) -> bool:
    for pat in patterns:
        p = str(pat).strip()
        if not p:
            continue
        if fnmatch(cid, p) or (p in cid):
            return True
    return False


def run(cfg) -> int:
    """
    Figure 18 (SVG):
      Aitchison distance distributions across methods.

      - x-axis: methods
      - y-axis: Aitchison distance (sample -> expected composition)
      - mark type: violin/boxplot + jittered sample dots
    """
    log = get_logger()
    outdir = Path(cfg.paths.outdir)
    exclude_patterns = list(getattr(cfg.figures, "exclude_community_name_patterns", []) or [])
    exclude_patterns += list(getattr(cfg.figures, "exclude_community_name_patterns_fig18", []) or [])

    long_p = outdir / "distances_all_long.tsv"
    if not long_p.exists():
        log.info("Figure 18: distances_all_long.tsv not found; run compute-distances first.")
        return 0

    long = pd.read_csv(long_p, sep="\t")
    if long.empty:
        log.info("Figure 18: no distance rows to plot.")
        return 0

    long = long[long["metric"].astype(str) == "aitchison"].copy()
    if long.empty:
        log.info(
            "Figure 18 skipped: no 'aitchison' rows found in distances_all_long.tsv. "
            "Add 'aitchison' to distances.metrics and rerun compute-distances."
        )
        return 0

    if exclude_patterns:
        keep = ~long["community_id"].astype(str).map(lambda cid: _is_excluded(cid, exclude_patterns))
        long = long[keep].copy()
        if long.empty:
            log.info("Figure 18: no rows left after figure-level community exclusions.")
            return 0

    long["distance"] = pd.to_numeric(long["distance"], errors="coerce")
    long = long[np.isfinite(long["distance"])].copy()
    if long.empty:
        log.info("Figure 18: all Aitchison rows are non-finite.")
        return 0

    long["sample_key"] = long["community_id"].astype(str) + "|" + long["sample_id"].astype(str)

    method_cfg = [str(m) for m in (getattr(cfg.distances, "methods", METHOD_ORDER) or METHOD_ORDER)]
    method_order = [m for m in method_cfg if m in set(long["method"].astype(str).tolist())]
    if not method_order:
        method_order = [m for m in METHOD_ORDER if m in set(long["method"].astype(str).tolist())]
    if not method_order:
        log.info("Figure 18: no methods available after filtering.")
        return 0

    figdir = figures_root(outdir) / "fig18_aitchison_distance_distributions"
    figdir.mkdir(parents=True, exist_ok=True)

    long[["community_id", "sample_id", "sample_key", "method", "metric", "distance"]].to_csv(
        figdir / "fig18_aitchison_points.tsv", sep="\t", index=False
    )

    summary_rows = []
    for method in method_order:
        vals = long.loc[long["method"] == method, "distance"].to_numpy(dtype=float)
        vals = vals[np.isfinite(vals)]
        if vals.size == 0:
            continue
        summary_rows.append(
            {
                "method": method,
                "n": int(vals.size),
                "min": float(np.min(vals)),
                "max": float(np.max(vals)),
                "mean": float(np.mean(vals)),
                "median": float(np.median(vals)),
            }
        )
    pd.DataFrame(summary_rows).to_csv(figdir / "fig18_aitchison_by_method_summary.tsv", sep="\t", index=False)

    plot_type = _normalize_plot_type(getattr(cfg.figures, "fig18_plot_type", "violin"), default="violin")
    font_size = float(getattr(cfg.figures, "fig18_font_size", 14.0))
    point_size = max(1.0, float(getattr(cfg.figures, "box_violin_dot_size", 22.0)))
    title = getattr(cfg.figures, "fig18_title", None) or "Aitchison Distance by Method"

    positions = np.arange(1, len(method_order) + 1, dtype=float)
    method_data: list[tuple[int, str, np.ndarray]] = []
    for i, method in enumerate(method_order):
        vals = long[long["method"] == method]["distance"].to_numpy(dtype=float)
        vals = vals[np.isfinite(vals)]
        if vals.size > 0:
            method_data.append((i, method, vals))

    if not method_data:
        log.info("Figure 18: no finite Aitchison values after method grouping.")
        return 0

    fig = plt.figure(figsize=(10.5, 6.5))
    ax = fig.add_subplot(1, 1, 1)

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
            )
            for body, method in zip(vp["bodies"], violin_methods):
                body.set_facecolor(METHOD_COLORS.get(method, "#999999"))
                body.set_edgecolor("#333333")
                body.set_alpha(0.5)
            if "cmedians" in vp:
                vp["cmedians"].set_color("#111111")
                vp["cmedians"].set_linewidth(1.2)
        else:
            # Fallback when all methods have one sample only.
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

    rng = np.random.default_rng(1)
    for i, method, vals in method_data:
        xj = positions[i] + rng.normal(0.0, 0.04, size=vals.size)
        ax.scatter(
            xj,
            vals,
            s=point_size,
            alpha=0.40,
            color=METHOD_COLORS.get(method, "#666666"),
            edgecolors="none",
            zorder=3,
        )

    ax.set_xticks(positions)
    ax.set_xticklabels([display_method(m) for m in method_order], rotation=0, fontsize=max(8.0, font_size - 1.0))
    ax.set_xlabel("Method", fontsize=font_size)
    ax.set_ylabel("Aitchison Distance", fontsize=font_size)
    ax.tick_params(axis="y", labelsize=max(8.0, font_size - 1.0))
    ax.grid(True, axis="y", alpha=0.25, linewidth=0.6)
    ax.set_title(title, fontsize=font_size + 2.0)

    yvals = long["distance"].to_numpy(dtype=float)
    yvals = yvals[np.isfinite(yvals)]
    if yvals.size > 0:
        ymin = float(np.min(yvals))
        ymax = float(np.max(yvals))
        if ymax > ymin:
            pad = 0.06 * (ymax - ymin)
            ax.set_ylim(ymin - pad, ymax + pad)

    suffix = "violin" if plot_type == "violin" else "boxplot"
    out_svg = figdir / f"fig18_aitchison_distance_{suffix}.svg"
    fig.tight_layout()
    fig.savefig(out_svg, format="svg", dpi=300, bbox_inches="tight")
    plt.close(fig)

    log.info(f"Wrote Figure 18 to {figdir}")
    return 0

