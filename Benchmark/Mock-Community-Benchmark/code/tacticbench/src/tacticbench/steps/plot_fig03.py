from __future__ import annotations

from fnmatch import fnmatch
from pathlib import Path
import math
import re
from typing import Dict, List

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from ..log import get_logger
from ..utils.paths import figures_root
from ..utils.labels import display_methods, display_sample_key


METHOD_AXES = ["ASV", "zOTU", "OTU", "TIC", "TAC"]
METRICS = ["jaccard", "weighted_jaccard", "braycurtis", "unifrac", "gunifrac_a0.0", "gunifrac_a0.5", "gunifrac_a1.0"]
PHYLO_METRICS = {"unifrac", "gunifrac_a0.0", "gunifrac_a0.5", "gunifrac_a1.0"}

METRIC_COLORS = {
    # Color-blind-friendly palette (Okabe-Ito family hues)
    "jaccard": "#E69F00",
    "weighted_jaccard": "#CC79A7",
    "braycurtis": "#D55E00",
    "unifrac": "#56B4E9",
    "gunifrac_a0.0": "#F0E442",
    "gunifrac_a0.5": "#0072B2",
    "gunifrac_a1.0": "#009E73",
}


def _sanitize_filename(s: str) -> str:
    s = s.replace("|", "__")
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", s)
    return s[:200]


def _is_excluded(cid: str, patterns: list[str]) -> bool:
    for pat in patterns:
        p = str(pat).strip()
        if not p:
            continue
        if fnmatch(cid, p) or (p in cid):
            return True
    return False


def _radar(ax, values: List[float], labels: List[str], rmax: float = 1.0, title: str | None = None) -> None:
    n = len(labels)
    angles = [2 * math.pi * i / n for i in range(n)]
    angles += angles[:1]

    vals = list(values) + [values[0]]
    ax.plot(angles, vals)
    ax.fill(angles, vals, alpha=0.08)

    ax.set_ylim(0, rmax)
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(labels, fontsize=9)

    # radial grid rings at 0.1 ratio increments
    yticks = np.arange(0.0, rmax + 1e-9, 0.1)
    ylabels = [f"{t:.1f}" for t in yticks]
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=7)

    if title:
        ax.set_title(title, fontsize=10)


def _get_values(block: pd.DataFrame, metric: str, methods: List[str]) -> List[float]:
    """
    block: subset for one sample_key (rows: method/metric combos)
    """
    out = []
    for m in methods:
        sub = block[(block["method"] == m) & (block["metric"] == metric)]
        if sub.empty:
            out.append(float("nan"))
        else:
            out.append(float(sub["distance"].iloc[0]))
    return out


def _metric_style(metric: str) -> tuple[str, str]:
    color = METRIC_COLORS.get(metric, "#444444")
    linestyle = "-"
    return color, linestyle


def _radar_multi_metric(ax, block: pd.DataFrame, labels: List[str], title: str | None = None, rmax: float = 1.0) -> None:
    n = len(labels)
    angles = [2 * math.pi * i / n for i in range(n)]
    angles += angles[:1]

    ax.set_ylim(0, rmax)
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(display_methods(labels), fontsize=9)
    yticks = np.arange(0.0, rmax + 1e-9, 0.1)
    ylabels = [f"{t:.1f}" for t in yticks]
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=7)

    for metric in METRICS:
        vals = _get_values(block, metric, labels)
        vals = [min(1.0, max(0.0, v)) if np.isfinite(v) else np.nan for v in vals]
        vals = vals + [vals[0]]
        color, ls = _metric_style(metric)
        ax.plot(angles, vals, color=color, linestyle=ls, linewidth=2.3, label=metric)

    if title:
        ax.set_title(title, fontsize=10)


def _legend_handles() -> List[Line2D]:
    out: List[Line2D] = []
    for metric in METRICS:
        color, ls = _metric_style(metric)
        group = "phylo" if metric in PHYLO_METRICS else "non-phylo"
        out.append(Line2D([0], [0], color=color, linestyle=ls, linewidth=2.6, label=f"{metric} ({group})"))
    return out


def run(cfg) -> int:
    """
    Figure 03 (SVG):
      - Option A (fig03_by_sample): for each sample, generate one multi-line radar SVG
        plus one all-samples SVG (grid of sample radars)
      - Option B (fig03_by_metric): for each metric, generate one (or multiple) SVG grid(s) of radar plots across samples

    Radar:
      - axes: methods (ASV, zOTU, OTU, TIC, TAC)
      - scale: 0..1
    """
    log = get_logger()
    outdir = Path(cfg.paths.outdir)
    exclude_patterns = list(getattr(cfg.figures, "exclude_community_name_patterns", []) or [])
    exclude_patterns += list(getattr(cfg.figures, "exclude_community_name_patterns_fig03", []) or [])

    long_p = outdir / "distances_all_long.tsv"
    if not long_p.exists():
        log.info("Figure 03: distances_all_long.tsv not found; run compute-distances first.")
        return 0

    long = pd.read_csv(long_p, sep="\t")
    if long.empty:
        log.info("Figure 03: no distance rows to plot.")
        return 0
    if exclude_patterns:
        keep = ~long["community_id"].astype(str).map(lambda cid: _is_excluded(cid, exclude_patterns))
        long = long[keep].copy()
        if long.empty:
            log.info("Figure 03: no rows left after figure-level community exclusions.")
            return 0

    # restrict to what Figure 03 expects
    long = long[long["method"].isin(METHOD_AXES) & long["metric"].isin(METRICS)].copy()
    if long.empty:
        log.info("Figure 03: no matching method/metric rows after filtering.")
        return 0

    long["sample_key"] = long["community_id"].astype(str) + "|" + long["sample_id"].astype(str)

    fr = figures_root(outdir)
    fr.mkdir(parents=True, exist_ok=True)
    fig03_by_sample_title = getattr(cfg.figures, "fig03_by_sample_title", None)
    fig03_by_metric_title = getattr(cfg.figures, "fig03_by_metric_title", None)

    # ---------- Option A: one file per sample ----------
    if getattr(cfg.figures, "fig03_by_sample", True):
        out_dir = fr / "fig03_radar_by_sample"
        out_dir.mkdir(parents=True, exist_ok=True)

        sample_blocks: Dict[str, pd.DataFrame] = {}
        for sample_key, block in long.groupby("sample_key", sort=True):
            # ensure unique (method,metric) by taking first
            block = (
                block.sort_values(["method", "metric"])
                .drop_duplicates(subset=["method", "metric"], keep="first")
                .reset_index(drop=True)
            )
            sample_blocks[sample_key] = block

            fig = plt.figure(figsize=(7.4, 7.4))
            ax = fig.add_subplot(1, 1, 1, projection="polar")
            _radar_multi_metric(ax, block, METHOD_AXES, title=display_sample_key(sample_key), rmax=1.0)
            legend_y = 0.90 if fig03_by_sample_title else 0.95
            fig.legend(
                handles=_legend_handles(),
                loc="upper center",
                bbox_to_anchor=(0.5, legend_y),
                ncol=2,
                fontsize=9,
                frameon=False,
            )
            if fig03_by_sample_title:
                fig.suptitle(fig03_by_sample_title, fontsize=13, y=0.985)
                fig.tight_layout(rect=[0, 0, 1, 0.80])
            else:
                fig.tight_layout(rect=[0, 0, 1, 0.86])

            out_svg = out_dir / f"{_sanitize_filename(sample_key)}.svg"
            fig.savefig(out_svg, format="svg", bbox_inches="tight", pad_inches=0.15)
            plt.close(fig)

        # one SVG for all samples (single page grid)
        sample_keys = sorted(sample_blocks.keys())
        if sample_keys:
            C = max(1, int(getattr(cfg.figures, "grid_cols", 6)))
            R = int(math.ceil(len(sample_keys) / C))
            fig = plt.figure(figsize=(C * 3.4, max(1, R) * 3.2))

            for i, sk in enumerate(sample_keys):
                ax = fig.add_subplot(R, C, i + 1, projection="polar")
                _radar_multi_metric(ax, sample_blocks[sk], METHOD_AXES, title=display_sample_key(sk), rmax=1.0)
            for j in range(len(sample_keys), R * C):
                ax = fig.add_subplot(R, C, j + 1, projection="polar")
                ax.axis("off")

            legend_y = 0.93 if fig03_by_sample_title else 0.97
            fig.legend(
                handles=_legend_handles(),
                loc="upper center",
                bbox_to_anchor=(0.5, legend_y),
                ncol=3,
                fontsize=9,
                frameon=False,
            )
            if fig03_by_sample_title:
                fig.suptitle(fig03_by_sample_title, fontsize=13, y=0.985)
                fig.tight_layout(rect=[0, 0, 1, 0.85])
            else:
                fig.tight_layout(rect=[0, 0, 1, 0.90])
            fig.savefig(out_dir / "all_samples.svg", format="svg", bbox_inches="tight", pad_inches=0.15)
            plt.close(fig)

        log.info(f"Wrote Figure 03 (by sample) to {out_dir}")

    # ---------- Option B: one (or multiple pages) per metric ----------
    if getattr(cfg.figures, "fig03_by_metric", True):
        out_dir = fr / "fig03_radar_by_metric"
        out_dir.mkdir(parents=True, exist_ok=True)

        # Dedicated by-metric grid size (optional), with backward-compatible fallback.
        R = int(getattr(cfg.figures, "fig03_by_metric_grid_rows", None) or getattr(cfg.figures, "grid_rows", 5))
        C = int(getattr(cfg.figures, "fig03_by_metric_grid_cols", None) or getattr(cfg.figures, "grid_cols", 6))
        per_page = R * C

        sample_keys = sorted(long["sample_key"].unique().tolist())

        # pre-index for fast lookup: (sample_key, method, metric) -> distance
        idx: Dict[tuple, float] = {}
        for _, row in long.iterrows():
            k = (row["sample_key"], row["method"], row["metric"])
            if k not in idx:
                idx[k] = float(row["distance"])

        def values_for(sk: str, metric: str) -> List[float]:
            vals = []
            for m in METHOD_AXES:
                v = idx.get((sk, m, metric), float("nan"))
                vals.append(min(1.0, max(0.0, v)) if np.isfinite(v) else np.nan)
            return vals

        for metric in METRICS:
            # chunk into pages if needed
            for p0 in range(0, len(sample_keys), per_page):
                chunk = sample_keys[p0 : p0 + per_page]
                page = (p0 // per_page) + 1

                fig = plt.figure(figsize=(C * 3.3, R * 3.0))
                for i in range(per_page):
                    ax = fig.add_subplot(R, C, i + 1, projection="polar")
                    if i >= len(chunk):
                        ax.axis("off")
                        continue
                    sk = chunk[i]
                    vals = values_for(sk, metric)
                    _radar(ax, vals, display_methods(METHOD_AXES), rmax=1.0, title=display_sample_key(sk))

                base = fig03_by_metric_title if fig03_by_metric_title else metric
                fig.suptitle(f"{base} (page {page})", fontsize=13)
                fig.tight_layout(rect=[0, 0, 1, 0.95])

                suffix = f"_p{page:02d}" if len(sample_keys) > per_page else ""
                out_svg = out_dir / f"{metric}{suffix}.svg"
                fig.savefig(out_svg, format="svg")
                plt.close(fig)

        log.info(f"Wrote Figure 03 (by metric) to {out_dir}")

    return 0
