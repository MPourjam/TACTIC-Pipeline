from __future__ import annotations

from fnmatch import fnmatch
from pathlib import Path

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd

from ..log import get_logger
from ..utils.paths import figures_root
from ..utils.labels import display_method, display_sample_key

METHOD_ORDER = ["ASV", "zOTU", "OTU", "TIC", "TAC"]
METRIC_ORDER = [
    "jaccard",
    "weighted_jaccard",
    "braycurtis",
    "unifrac",
    "gunifrac_a0.0",
    "gunifrac_a0.5",
    "gunifrac_a1.0",
]
METRIC_LABEL = {
    "jaccard": "Jaccard",
    "weighted_jaccard": "Weighted Jaccard",
    "braycurtis": "Bray-Curtis",
    "unifrac": "UniFrac",
    "gunifrac_a0.0": "gUniFrac (alpha=0.0)",
    "gunifrac_a0.5": "gUniFrac (alpha=0.5)",
    "gunifrac_a1.0": "gUniFrac (alpha=1)",
}


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
    Figure 02 (SVG):
      Heatmap with:
        - rows: method|metric combinations
        - columns: samples (community_id|sample_id)
        - values: distance (0..1), NaN for missing (e.g., phylo tools absent)
    """
    log = get_logger()
    outdir = Path(cfg.paths.outdir)
    exclude_patterns = list(getattr(cfg.figures, "exclude_community_name_patterns", []) or [])
    exclude_patterns += list(getattr(cfg.figures, "exclude_community_name_patterns_fig02", []) or [])

    long_p = outdir / "distances_all_long.tsv"
    if not long_p.exists():
        log.info("Figure 02: distances_all_long.tsv not found; run compute-distances first.")
        return 0

    long = pd.read_csv(long_p, sep="\t")
    if long.empty:
        log.info("Figure 02: no distance rows to plot.")
        return 0
    if exclude_patterns:
        keep = ~long["community_id"].astype(str).map(lambda cid: _is_excluded(cid, exclude_patterns))
        long = long[keep].copy()
        if long.empty:
            log.info("Figure 02: no rows left after figure-level community exclusions.")
            return 0

    long["sample_key"] = long["community_id"].astype(str) + "|" + long["sample_id"].astype(str)

    mat = long.pivot_table(
        index="sample_key",
        columns=["method", "metric"],
        values="distance",
        aggfunc="first",
    )
    mat = mat.sort_index()

    fr = figures_root(outdir)
    fr.mkdir(parents=True, exist_ok=True)
    out_dir = fr / "fig02_distance_heatmaps"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig02_title = getattr(cfg.figures, "fig02_title", None)

    def _draw_heatmap(
        mat_plot: pd.DataFrame,
        row_labels: list[str],
        out_svg: Path,
        title: str,
        row_label_rotation: float = 0.0,
        row_label_fontsize: int = 22,
        title_fontsize: int = 26,
        method_spacers: list[tuple[str, int]] | None = None,
    ) -> None:
        data = mat_plot.to_numpy(dtype=float)
        data_masked = np.ma.array(data, mask=np.isnan(data))

        cell_size = 0.7
        fig_w = max(12.0, mat_plot.shape[1] * cell_size)
        fig_h = max(12.0, mat_plot.shape[0] * cell_size)
        fig = plt.figure(figsize=(fig_w, fig_h))
        ax = fig.add_subplot(1, 1, 1)

        cmap = plt.get_cmap("cividis").copy()
        cmap.set_bad(color="#f0f0f0")
        im = ax.imshow(data_masked, aspect="auto", cmap=cmap, vmin=0.0, vmax=1.0)

        ax.set_yticks(range(mat_plot.shape[0]))
        ax.set_yticklabels(row_labels, rotation=row_label_rotation, fontsize=row_label_fontsize, ha="right")
        ax.set_xticks(range(mat_plot.shape[1]))
        col_labels = [display_sample_key(str(x)) for x in mat_plot.columns.tolist()]
        ax.set_xticklabels(col_labels, rotation=45, ha="right", rotation_mode="anchor", fontsize=22)
        ax.set_title(title, fontsize=title_fontsize, y=1.05)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.14)
        cbar = fig.colorbar(im, cax=cax)
        cbar.ax.tick_params(labelsize=21)
        if method_spacers:
            mid_x = (mat_plot.shape[1] - 1) / 2.0
            for method, spacer_idx in method_spacers:
                if 0 <= spacer_idx < mat_plot.shape[0]:
                    ax.text(
                        mid_x,
                        spacer_idx,
                        display_method(method),
                        fontsize=26,
                        fontweight="bold",
                        ha="center",
                        va="center",
                        color="black",
                    )
        ax.set_aspect("equal", adjustable="box")
        fig.tight_layout(rect=[0, 0, 1, 0.92])
        fig.savefig(out_svg, format="svg")
        plt.close(fig)

    # combined heatmap: rows grouped as method|metric, columns = samples
    mat_t = mat.T
    present_methods = mat_t.index.get_level_values(0).unique().tolist()
    method_order = [m for m in METHOD_ORDER if m in present_methods] + [m for m in present_methods if m not in METHOD_ORDER]
    ordered_rows = []
    for method in method_order:
        for metric in METRIC_ORDER:
            key = (method, metric)
            if key in mat_t.index:
                ordered_rows.append(key)
    combined_df = None
    combined_labels: list[str] = []
    method_spacers: list[tuple[str, int]] = []
    if ordered_rows:
        mat_t = mat_t.loc[ordered_rows]
        mat_cols = mat_t.columns
        combined_rows: list[np.ndarray] = []
        combined_labels = []
        for method in method_order:
            method_rows = [idx for idx in ordered_rows if idx[0] == method]
            if not method_rows:
                continue
            spacer_idx = len(combined_rows)
            combined_rows.append(np.full(len(mat_cols), np.nan))
            combined_labels.append("")
            method_spacers.append((method, spacer_idx))
            for idx in method_rows:
                combined_rows.append(mat_t.loc[idx].to_numpy(dtype=float))
                combined_labels.append(METRIC_LABEL.get(idx[1], str(idx[1])))
        if combined_rows:
            combined_df = pd.DataFrame(combined_rows, columns=mat_cols)
    if combined_df is not None:
        combined_title = fig02_title if fig02_title else "All Methods"
        _draw_heatmap(
            combined_df,
            combined_labels,
            out_dir / "all_methods.svg",
            combined_title,
            row_label_rotation=45,
            row_label_fontsize=22,
            title_fontsize=32,
            method_spacers=method_spacers,
        )

    # separate heatmap per method: rows = metrics, columns = samples
    for method in method_order:
        if method not in mat.columns.get_level_values(0):
            continue
        sub = mat.xs(method, axis=1, level=0, drop_level=True)
        metric_cols = [k for k in METRIC_ORDER if k in sub.columns]
        if not metric_cols:
            continue
        sub = sub[metric_cols]
        sub_t = sub.T
        row_labels = [METRIC_LABEL.get(k, str(k)) for k in sub_t.index.to_list()]
        out_svg = out_dir / f"{method}_distance_heatmap.svg"
        method_title = f"{fig02_title}: {display_method(method)}" if fig02_title else f"{display_method(method)}"
        _draw_heatmap(sub_t, row_labels, out_svg, method_title)

    log.info(f"Wrote Figure 02 heatmaps to {out_dir}")
    return 0
