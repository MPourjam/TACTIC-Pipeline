from __future__ import annotations

from fnmatch import fnmatch
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse
import numpy as np
import pandas as pd

from ..log import get_logger
from ..utils.labels import display_method
from ..utils.paths import figures_root

METHOD_ORDER = ["ASV", "zOTU", "OTU", "TIC", "TAC"]
METHOD_MARKERS = {
    "ASV": "o",
    "zOTU": "^",
    "OTU": "s",
    "TIC": "D",
    "TAC": "P",
}
METHOD_COLORS = {
    "ASV": "#0072B2",
    "zOTU": "#56B4E9",
    "OTU": "#D55E00",
    "TIC": "#009E73",
    "TAC": "#E69F00",
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


def _is_excluded(cid: str, patterns: list[str]) -> bool:
    for pat in patterns:
        p = str(pat).strip()
        if not p:
            continue
        if fnmatch(cid, p) or (p in cid):
            return True
    return False


def _metric_label(metric: str) -> str:
    return METRIC_LABELS.get(str(metric), str(metric))


def _normalize_profile_distance(raw: object) -> str:
    v = str(raw or "").strip().lower()
    if v in {"euclidean", "l2"}:
        return "euclidean"
    return "euclidean"


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
    if df.empty:
        return pd.DataFrame()
    df["sample_key"] = df["community_id"] + "|" + df["sample_id"]
    df["method_label"] = df["method"].map(display_method)
    df["metric_label"] = df["metric"].map(_metric_label)
    return df


def _build_distance_matrix(x: np.ndarray, metric: str = "euclidean") -> np.ndarray:
    if metric == "euclidean":
        diff = x[:, None, :] - x[None, :, :]
        d2 = np.sum(diff * diff, axis=2)
        np.maximum(d2, 0.0, out=d2)
        return np.sqrt(d2)
    # Future-proof fallback.
    diff = x[:, None, :] - x[None, :, :]
    d2 = np.sum(diff * diff, axis=2)
    np.maximum(d2, 0.0, out=d2)
    return np.sqrt(d2)


def _pcoa(distance_matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Classical PCoA (Torgerson / metric MDS):
      B = -1/2 * J * D^2 * J
      eig(B) -> coordinates = V * sqrt(lambda)
    Returns:
      coords: n x k (k = number of positive eigenvalues)
      eigvals_pos: positive eigenvalues in descending order
    """
    n = distance_matrix.shape[0]
    if n == 0:
        return np.empty((0, 0), dtype=float), np.empty((0,), dtype=float)
    d2 = np.square(distance_matrix, dtype=float)
    j = np.eye(n, dtype=float) - (np.ones((n, n), dtype=float) / float(n))
    b = -0.5 * (j @ d2 @ j)
    eigvals, eigvecs = np.linalg.eigh(b)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    pos = eigvals > 1e-12
    if not np.any(pos):
        return np.zeros((n, 0), dtype=float), np.empty((0,), dtype=float)
    eigvals_pos = eigvals[pos]
    eigvecs_pos = eigvecs[:, pos]
    coords = eigvecs_pos * np.sqrt(eigvals_pos)[None, :]
    return coords, eigvals_pos


def _draw_cov_ellipse(
    ax,
    x: np.ndarray,
    y: np.ndarray,
    color: str,
    n_std: float = 2.0,
    lw: float = 1.3,
    alpha: float = 0.85,
) -> bool:
    """
    Draw covariance ellipse for 2D points.
    Returns True when an ellipse was drawn.
    """
    xx = np.asarray(x, dtype=float)
    yy = np.asarray(y, dtype=float)
    ok = np.isfinite(xx) & np.isfinite(yy)
    xx = xx[ok]
    yy = yy[ok]
    if xx.size < 3:
        return False

    cov = np.cov(xx, yy, ddof=1)
    if cov.shape != (2, 2) or (not np.all(np.isfinite(cov))):
        return False
    # Guard against degenerate covariance.
    if np.linalg.det(cov) <= 0:
        return False

    vals, vecs = np.linalg.eigh(cov)
    order = np.argsort(vals)[::-1]
    vals = vals[order]
    vecs = vecs[:, order]
    if vals[0] <= 0 or vals[1] <= 0:
        return False

    width = 2.0 * n_std * float(np.sqrt(vals[0]))
    height = 2.0 * n_std * float(np.sqrt(vals[1]))
    angle = float(np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0])))
    center = (float(np.mean(xx)), float(np.mean(yy)))

    ell = Ellipse(
        xy=center,
        width=width,
        height=height,
        angle=angle,
        facecolor="none",
        edgecolor=color,
        linewidth=lw,
        alpha=alpha,
        zorder=2,
    )
    ax.add_patch(ell)
    return True


def _factor_r2(values: np.ndarray, groups: np.ndarray) -> float:
    """
    One-way ANOVA-style effect size (eta^2 / R^2) for a categorical factor.
    """
    y = np.asarray(values, dtype=float)
    g = np.asarray(groups, dtype=str)
    ok = np.isfinite(y)
    y = y[ok]
    g = g[ok]
    if y.size <= 1:
        return 0.0
    grand = float(np.mean(y))
    sst = float(np.sum((y - grand) ** 2))
    if sst <= 0.0:
        return 0.0
    ssb = 0.0
    for k in np.unique(g):
        mask = g == k
        if not np.any(mask):
            continue
        yk = y[mask]
        mk = float(np.mean(yk))
        ssb += float(yk.size) * ((mk - grand) ** 2)
    r2 = ssb / sst
    return float(np.clip(r2, 0.0, 1.0))


def _compute_umap_coords(distance_matrix: np.ndarray, cfg, log) -> tuple[np.ndarray | None, str]:
    try:
        import umap  # type: ignore
    except Exception:
        return None, "UMAP skipped: umap-learn is not installed."

    n = int(distance_matrix.shape[0])
    if n < 3:
        return None, "UMAP skipped: fewer than 3 points."
    n_neighbors = int(getattr(cfg.figures, "fig19_umap_n_neighbors", 15))
    n_neighbors = max(2, min(n_neighbors, n - 1))
    min_dist = float(getattr(cfg.figures, "fig19_umap_min_dist", 0.10))
    min_dist = float(np.clip(min_dist, 0.0, 1.0))
    random_state = int(getattr(cfg.figures, "fig19_umap_random_state", 42))

    try:
        reducer = umap.UMAP(
            n_components=2,
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            metric="precomputed",
            random_state=random_state,
        )
        coords = reducer.fit_transform(distance_matrix)
        coords = np.asarray(coords, dtype=float)
        if coords.shape[1] < 2:
            return None, "UMAP skipped: embedding has <2 components."
        return coords[:, :2], ""
    except Exception as e:
        return None, f"UMAP failed: {e}"


def _draw_ordination_panel(
    ax,
    profile_df: pd.DataFrame,
    x_col: str,
    y_col: str,
    method_order: list[str],
    metric_order: list[str],
    point_size: float,
    title: str,
    xlabel: str,
    ylabel: str,
    font_size: float,
) -> None:
    for method in method_order:
        marker = METHOD_MARKERS.get(method, "o")
        for metric in metric_order:
            sub = profile_df[(profile_df["method"] == method) & (profile_df["metric"] == metric)]
            if sub.empty:
                continue
            ax.scatter(
                sub[x_col].to_numpy(dtype=float),
                sub[y_col].to_numpy(dtype=float),
                s=point_size,
                marker=marker,
                c=METRIC_COLORS.get(metric, "#4D4D4D"),
                alpha=0.78,
                edgecolors="#111111",
                linewidths=0.35,
            )

    for metric in metric_order:
        sub = profile_df[profile_df["metric"] == metric]
        if sub.empty:
            continue
        _draw_cov_ellipse(
            ax=ax,
            x=sub[x_col].to_numpy(dtype=float),
            y=sub[y_col].to_numpy(dtype=float),
            color=METRIC_COLORS.get(metric, "#4D4D4D"),
            n_std=2.0,
            lw=1.4,
            alpha=0.9,
        )

    ax.grid(True, alpha=0.25, linewidth=0.6)
    ax.set_xlabel(xlabel, fontsize=font_size)
    ax.set_ylabel(ylabel, fontsize=font_size)
    ax.tick_params(axis="both", labelsize=max(8.0, font_size - 1.0))
    ax.set_title(title, fontsize=font_size + 1.0)


def run(cfg) -> int:
    """
    Figure 19:
      PCoA of (sample, method, metric) points using metric-contrast profiles.

    Point definition:
      each row in the ordination corresponds to one (community_id, sample_id, method, metric).

    Feature embedding (contrast profile):
      for fixed (community_id, sample_id, method) and selected metrics K,
      profile for metric m is:
        [d_m - d_k for k in K]

    Distance between points:
      Euclidean distance on (optionally standardized) contrast profiles.
    """
    log = get_logger()
    outdir = Path(cfg.paths.outdir)

    df = _load_distance_data(cfg)
    if df.empty:
        log.info("Figure 19: no finite distance data found in distances_all_long.tsv.")
        return 0

    # Figure-level exclusions
    exclude_patterns = list(getattr(cfg.figures, "exclude_community_name_patterns", []) or [])
    exclude_patterns += list(getattr(cfg.figures, "exclude_community_name_patterns_fig19", []) or [])
    if exclude_patterns:
        df = df[~df["community_id"].map(lambda cid: _is_excluded(cid, exclude_patterns))].copy()
        if df.empty:
            log.info("Figure 19: no rows left after figure-level community exclusions.")
            return 0

    methods_cfg = [str(m) for m in (getattr(cfg.distances, "methods", METHOD_ORDER) or METHOD_ORDER)]
    methods_present = set(df["method"].astype(str).tolist())
    method_order = [m for m in methods_cfg if m in methods_present]
    if not method_order:
        method_order = [m for m in METHOD_ORDER if m in methods_present]
    if not method_order:
        log.info("Figure 19: no methods available after filtering.")
        return 0

    # Metrics: configurable, but always exclude aitchison for this figure.
    allowed_metrics = set(METRIC_LABELS.keys())
    metrics_cfg = [str(m) for m in (getattr(cfg.figures, "fig19_metrics", []) or []) if str(m).strip()]
    if metrics_cfg:
        metric_order = [m for m in metrics_cfg if m in allowed_metrics]
    else:
        metric_order = [m for m in (getattr(cfg.distances, "metrics", []) or []) if m in allowed_metrics]
    if not metric_order:
        metric_order = [m for m in METRIC_LABELS.keys() if m in set(df["metric"].astype(str).tolist())]

    df = df[df["method"].isin(method_order) & df["metric"].isin(metric_order)].copy()
    if df.empty:
        log.info("Figure 19: no rows left after method/metric filtering.")
        return 0

    if len(metric_order) < 2:
        log.info("Figure 19 skipped: need at least 2 non-Aitchison metrics for contrast profiles.")
        return 0

    # Keep only complete (community, sample, method) groups with all selected metrics.
    grp_cols = ["community_id", "sample_id", "sample_key", "method", "method_label"]
    have = (
        df.groupby(grp_cols, as_index=False)["metric"]
        .nunique()
        .rename(columns={"metric": "n_metrics_present"})
    )
    complete = have[have["n_metrics_present"] == len(metric_order)][grp_cols].copy()
    n_dropped = int(have.shape[0] - complete.shape[0])
    if complete.empty:
        log.info("Figure 19: no complete sample-method groups contain all selected metrics.")
        return 0
    if n_dropped > 0:
        log.info(f"Figure 19: dropped {n_dropped} incomplete sample-method groups.")

    df = df.merge(complete, on=grp_cols, how="inner")

    pivot = (
        df.pivot_table(
            index=grp_cols,
            columns="metric",
            values="distance",
            aggfunc="first",
        )
        .reindex(columns=metric_order)
        .reset_index()
    )
    pivot = pivot.dropna(subset=metric_order).copy()
    if pivot.empty:
        log.info("Figure 19: no complete rows left after pivot.")
        return 0

    # Build one point per (sample, method, metric) with metric-contrast profile.
    profile_rows: list[dict] = []
    profile_vectors: list[np.ndarray] = []

    for _, row in pivot.iterrows():
        base = np.array([float(row[m]) for m in metric_order], dtype=float)
        for midx, metric in enumerate(metric_order):
            vec = base[midx] - base
            profile_vectors.append(vec)
            point_id = f"{row['community_id']}|{row['sample_id']}|{row['method']}|{metric}"
            rec = {
                "point_id": point_id,
                "community_id": str(row["community_id"]),
                "sample_id": str(row["sample_id"]),
                "sample_key": str(row["sample_key"]),
                "method": str(row["method"]),
                "method_label": str(row["method_label"]),
                "metric": metric,
                "metric_label": _metric_label(metric),
                "distance_value": float(base[midx]),
            }
            for j, k in enumerate(metric_order):
                rec[f"contrast_vs_{k}"] = float(vec[j])
            profile_rows.append(rec)

    profile_df = pd.DataFrame(profile_rows)
    x = np.vstack(profile_vectors).astype(float)
    if x.shape[0] < 3:
        log.info("Figure 19: too few points to compute a stable ordination.")
        return 0

    standardize = bool(getattr(cfg.figures, "fig19_standardize_profiles", True))
    if standardize:
        mu = np.mean(x, axis=0)
        sd = np.std(x, axis=0, ddof=0)
        sd_safe = np.where(sd > 0.0, sd, 1.0)
        x = (x - mu[None, :]) / sd_safe[None, :]

    profile_distance = _normalize_profile_distance(getattr(cfg.figures, "fig19_profile_distance", "euclidean"))
    dmat = _build_distance_matrix(x, metric=profile_distance)
    coords, eigvals = _pcoa(dmat)
    if coords.shape[1] < 2:
        log.info("Figure 19: fewer than 2 positive PCoA axes; cannot draw 2D plot.")
        return 0

    profile_df["PC1"] = coords[:, 0]
    profile_df["PC2"] = coords[:, 1]
    if coords.shape[1] >= 3:
        profile_df["PC3"] = coords[:, 2]

    # Variance explained from positive eigenvalues.
    eig_sum = float(np.sum(eigvals)) if eigvals.size else 0.0
    if eig_sum > 0.0:
        ve = eigvals / eig_sum
        pc1_pct = 100.0 * float(ve[0]) if ve.size >= 1 else 0.0
        pc2_pct = 100.0 * float(ve[1]) if ve.size >= 2 else 0.0
    else:
        pc1_pct = 0.0
        pc2_pct = 0.0

    show_umap = bool(getattr(cfg.figures, "fig19_show_umap", True))
    umap_coords = None
    umap_msg = ""
    if show_umap:
        umap_coords, umap_msg = _compute_umap_coords(dmat, cfg=cfg, log=log)
        if umap_coords is not None:
            profile_df["UMAP1"] = umap_coords[:, 0]
            profile_df["UMAP2"] = umap_coords[:, 1]
        elif umap_msg:
            log.info(umap_msg)

    figdir = figures_root(outdir) / "fig19_metric_contrast_pcoa"
    figdir.mkdir(parents=True, exist_ok=True)

    profile_df.to_csv(figdir / "fig19_points_profiles_and_coords.tsv", sep="\t", index=False)
    if umap_coords is not None:
        profile_df[
            ["point_id", "community_id", "sample_id", "sample_key", "method", "metric", "UMAP1", "UMAP2"]
        ].to_csv(figdir / "fig19_umap_coords.tsv", sep="\t", index=False)
    else:
        # Prevent stale files from prior runs when UMAP was enabled.
        stale_umap_coords = figdir / "fig19_umap_coords.tsv"
        if stale_umap_coords.exists():
            stale_umap_coords.unlink()

    pd.DataFrame(
        {
            "axis": [f"PC{i+1}" for i in range(len(eigvals))],
            "eigenvalue": eigvals,
            "variance_explained": (eigvals / eig_sum) if eig_sum > 0.0 else np.zeros_like(eigvals),
        }
    ).to_csv(figdir / "fig19_pcoa_eigenvalues.tsv", sep="\t", index=False)

    # Axis-association diagnostics.
    diag_rows: list[dict] = []
    method_groups = profile_df["method"].astype(str).to_numpy()
    metric_groups = profile_df["metric"].astype(str).to_numpy()
    for i in range(coords.shape[1]):
        axis_name = f"PC{i+1}"
        v = coords[:, i]
        r2_method = _factor_r2(v, method_groups)
        r2_metric = _factor_r2(v, metric_groups)
        if r2_method > r2_metric:
            dominant = "method"
        elif r2_metric > r2_method:
            dominant = "metric"
        else:
            dominant = "tie"
        diag_rows.append(
            {
                "ordination": "PCoA",
                "axis": axis_name,
                "eigenvalue": float(eigvals[i]),
                "variance_explained": float((eigvals[i] / eig_sum) if eig_sum > 0.0 else 0.0),
                "r2_method": float(r2_method),
                "r2_metric": float(r2_metric),
                "dominant_factor": dominant,
                "n_points": int(profile_df.shape[0]),
                "n_methods": int(profile_df["method"].nunique()),
                "n_metrics": int(profile_df["metric"].nunique()),
            }
        )
    if umap_coords is not None:
        for i, axis_name in enumerate(["UMAP1", "UMAP2"]):
            v = profile_df[axis_name].to_numpy(dtype=float)
            r2_method = _factor_r2(v, method_groups)
            r2_metric = _factor_r2(v, metric_groups)
            if r2_method > r2_metric:
                dominant = "method"
            elif r2_metric > r2_method:
                dominant = "metric"
            else:
                dominant = "tie"
            diag_rows.append(
                {
                    "ordination": "UMAP",
                    "axis": axis_name,
                    "eigenvalue": np.nan,
                    "variance_explained": np.nan,
                    "r2_method": float(r2_method),
                    "r2_metric": float(r2_metric),
                    "dominant_factor": dominant,
                    "n_points": int(profile_df.shape[0]),
                    "n_methods": int(profile_df["method"].nunique()),
                    "n_metrics": int(profile_df["metric"].nunique()),
                }
            )
    pd.DataFrame(diag_rows).to_csv(
        figdir / "fig19_axis_association_diagnostics.tsv",
        sep="\t",
        index=False,
    )

    # Full distance matrix output.
    dist_df = pd.DataFrame(dmat, index=profile_df["point_id"], columns=profile_df["point_id"])
    dist_df.to_csv(figdir / "fig19_distance_matrix.tsv", sep="\t")

    # Plot
    font_size = float(getattr(cfg.figures, "fig19_font_size", 14.0))
    point_size = max(4.0, float(getattr(cfg.figures, "fig19_point_size", 36.0)))
    title = getattr(cfg.figures, "fig19_title", None) or "PCoA of Metric-Contrast Profiles"

    use_umap_panel = bool(umap_coords is not None)
    if use_umap_panel:
        fig, axes = plt.subplots(1, 2, figsize=(16.8, 8.0))
        ax_pcoa, ax_umap = axes
    else:
        fig = plt.figure(figsize=(10.8, 8.0))
        ax_pcoa = fig.add_subplot(1, 1, 1)
        ax_umap = None

    _draw_ordination_panel(
        ax=ax_pcoa,
        profile_df=profile_df,
        x_col="PC1",
        y_col="PC2",
        method_order=method_order,
        metric_order=metric_order,
        point_size=point_size,
        title=("PCoA" if use_umap_panel else title),
        xlabel=f"PC1 ({pc1_pct:.1f}% var. explained)",
        ylabel=f"PC2 ({pc2_pct:.1f}% var. explained)",
        font_size=font_size,
    )
    if use_umap_panel and ax_umap is not None:
        _draw_ordination_panel(
            ax=ax_umap,
            profile_df=profile_df,
            x_col="UMAP1",
            y_col="UMAP2",
            method_order=method_order,
            metric_order=metric_order,
            point_size=point_size,
            title="UMAP",
            xlabel="UMAP1",
            ylabel="UMAP2",
            font_size=font_size,
        )
        fig.suptitle(title, fontsize=font_size + 1.5, y=0.98)

    metric_handles = [
        Line2D(
            [0], [0],
            marker="o",
            color="none",
            markerfacecolor=METRIC_COLORS.get(metric, "#4D4D4D"),
            markeredgecolor="#111111",
            markersize=7,
            label=_metric_label(metric),
        )
        for metric in metric_order
    ]
    method_handles = [
        Line2D(
            [0], [0],
            marker=METHOD_MARKERS.get(method, "o"),
            color="none",
            markerfacecolor=METHOD_COLORS.get(method, "#777777"),
            markeredgecolor="#111111",
            markersize=7,
            label=display_method(method),
        )
        for method in method_order
    ]

    legend_host = ax_umap if ax_umap is not None else ax_pcoa
    leg1 = legend_host.legend(
        handles=metric_handles,
        title="Distance Metric",
        loc="upper left",
        bbox_to_anchor=(1.01, 1.00),
        frameon=False,
        fontsize=max(8.0, font_size - 1.5),
        title_fontsize=max(8.5, font_size - 1.0),
    )
    legend_host.add_artist(leg1)
    legend_host.legend(
        handles=method_handles,
        title="Analysis Method",
        loc="upper left",
        bbox_to_anchor=(1.01, 0.52),
        frameon=False,
        fontsize=max(8.0, font_size - 1.5),
        title_fontsize=max(8.5, font_size - 1.0),
    )

    out_svg = figdir / "fig19_metric_contrast_pcoa.svg"
    fig.tight_layout()
    fig.savefig(out_svg, format="svg", dpi=300, bbox_inches="tight")
    if use_umap_panel:
        fig.savefig(figdir / "fig19_metric_contrast_pcoa_umap.svg", format="svg", dpi=300, bbox_inches="tight")
    else:
        stale_umap_svg = figdir / "fig19_metric_contrast_pcoa_umap.svg"
        if stale_umap_svg.exists():
            stale_umap_svg.unlink()
    plt.close(fig)

    log.info(
        f"Wrote Figure 19 to {figdir} "
        f"(points={profile_df.shape[0]}, groups={pivot.shape[0]}, metrics={len(metric_order)}, "
        f"umap={'on' if use_umap_panel else 'off'})."
    )
    return 0
