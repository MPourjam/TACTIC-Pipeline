from __future__ import annotations

from fnmatch import fnmatch
from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..log import get_logger
from ..utils.io import genus_species
from ..utils.labels import display_community_id, display_method, display_methods, display_sample_key
from ..utils.paths import figures_root, method_dir, qc_dir

ALL_METHODS = ["ASV", "zOTU", "OTU", "TIC", "TAC"]
SOURCE_METHODS = ["ASV", "zOTU"]
CLUSTER_METHODS = ["OTU", "TIC", "TAC"]
STATUS_ORDER = ["captured_as_one", "still_split", "not_captured"]
STATUS_LABELS = {
    "captured_as_one": "Captured as one",
    "still_split": "Still split",
    "not_captured": "Not captured",
}
STATUS_COLORS = {
    "captured_as_one": "#009E73",
    "still_split": "#E69F00",
    "not_captured": "#CC79A7",
}


def _is_excluded(cid: str, patterns: list[str]) -> bool:
    for pat in patterns:
        p = str(pat).strip()
        if not p:
            continue
        if fnmatch(cid, p) or (p in cid):
            return True
    return False


def _tip_to_species(target_tip: str) -> str:
    s = str(target_tip).strip()
    if not s:
        return ""
    s = re.sub(r"__Copy\d+$", "", s)
    return genus_species(s)


def _target_species_series(mapping: pd.DataFrame) -> pd.Series:
    if "target_species" in mapping.columns:
        return mapping["target_species"].astype(str).map(genus_species)
    if "target_tip" in mapping.columns:
        return mapping["target_tip"].astype(str).map(_tip_to_species)
    return pd.Series([""] * len(mapping), index=mapping.index, dtype=object)


def _load_method_species_feature_counts(outdir: Path, community_id: str, method: str) -> pd.DataFrame:
    """
    Return per-sample, per-species inferred feature counts for one method.
    Columns: sample_id, target_species, n_features
    """
    md = method_dir(outdir, community_id, method)
    map_p = md / "mapping.tsv"
    tab_p = md / "table_std_tss1000.tsv"
    if not map_p.exists() or not tab_p.exists():
        return pd.DataFrame(columns=["sample_id", "target_species", "n_features"])

    try:
        mapping = pd.read_csv(map_p, sep="\t", dtype=str)
    except Exception:
        return pd.DataFrame(columns=["sample_id", "target_species", "n_features"])
    if mapping.empty or "inferred_feature" not in mapping.columns:
        return pd.DataFrame(columns=["sample_id", "target_species", "n_features"])

    mapping = mapping.copy()
    mapping["inferred_feature"] = mapping["inferred_feature"].astype(str)
    mapping["target_species"] = _target_species_series(mapping)
    mapping = mapping[mapping["target_species"].astype(str).str.strip() != ""].copy()
    if mapping.empty:
        return pd.DataFrame(columns=["sample_id", "target_species", "n_features"])
    mapping = mapping[["inferred_feature", "target_species"]].drop_duplicates()

    try:
        table = pd.read_csv(tab_p, sep="\t", index_col=0)
    except Exception:
        return pd.DataFrame(columns=["sample_id", "target_species", "n_features"])
    if table.empty:
        return pd.DataFrame(columns=["sample_id", "target_species", "n_features"])
    table = table.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    present = (table > 0).stack()
    present = present[present].reset_index()
    if present.empty:
        return pd.DataFrame(columns=["sample_id", "target_species", "n_features"])
    present.columns = ["inferred_feature", "sample_id", "present"]
    present = present.drop(columns=["present"])
    present["inferred_feature"] = present["inferred_feature"].astype(str)
    present["sample_id"] = present["sample_id"].astype(str)

    joined = present.merge(mapping, on="inferred_feature", how="inner")
    if joined.empty:
        return pd.DataFrame(columns=["sample_id", "target_species", "n_features"])

    counts = (
        joined.groupby(["sample_id", "target_species"], as_index=False)
        .agg(n_features=("inferred_feature", "nunique"))
        .sort_values(["sample_id", "target_species"], kind="mergesort")
    )
    return counts


def _status_from_count(n_features: int) -> str:
    if n_features <= 0:
        return "not_captured"
    if n_features == 1:
        return "captured_as_one"
    return "still_split"


def _summary_table(capture_df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for source_method in SOURCE_METHODS:
        for clustering_method in CLUSTER_METHODS:
            sub = capture_df[
                (capture_df["source_method"].astype(str) == source_method)
                & (capture_df["clustering_method"].astype(str) == clustering_method)
            ].copy()
            total = int(sub["event_id"].nunique()) if not sub.empty else 0
            for status in STATUS_ORDER:
                n_events = int((sub["capture_status"].astype(str) == status).sum()) if total > 0 else 0
                frac = float(n_events / total) if total > 0 else 0.0
                rows.append(
                    {
                        "source_method": source_method,
                        "clustering_method": clustering_method,
                        "capture_status": status,
                        "capture_status_label": STATUS_LABELS.get(status, status),
                        "n_events": n_events,
                        "n_total_events": total,
                        "fraction": frac,
                    }
                )
    return pd.DataFrame(rows)


def _plot_summary(
    summary_df: pd.DataFrame,
    out_path: Path,
    title: str,
    font_size: float,
) -> None:
    tick_size = max(7.0, font_size - 2.0)
    fig = plt.figure(figsize=(10.8, 4.8))
    axes = fig.subplots(1, len(SOURCE_METHODS), sharey=True)
    if len(SOURCE_METHODS) == 1:
        axes = [axes]

    for ax, source_method in zip(axes, SOURCE_METHODS):
        sub = summary_df[summary_df["source_method"].astype(str) == source_method].copy()
        x = np.arange(len(CLUSTER_METHODS))
        bottoms = np.zeros(len(CLUSTER_METHODS), dtype=float)

        for status in STATUS_ORDER:
            vals = []
            for cm in CLUSTER_METHODS:
                m = sub[
                    (sub["clustering_method"].astype(str) == cm)
                    & (sub["capture_status"].astype(str) == status)
                ]
                vals.append(float(m["fraction"].iloc[0]) if not m.empty else 0.0)
            vals_arr = np.array(vals, dtype=float)
            ax.bar(
                x,
                vals_arr,
                bottom=bottoms,
                width=0.72,
                color=STATUS_COLORS.get(status, "#999999"),
                edgecolor="#333333",
                linewidth=0.6,
                label=STATUS_LABELS.get(status, status),
            )
            bottoms += vals_arr

        totals = []
        for cm in CLUSTER_METHODS:
            m = sub[sub["clustering_method"].astype(str) == cm]
            totals.append(int(m["n_total_events"].iloc[0]) if not m.empty else 0)
        for i, total in enumerate(totals):
            ax.text(
                i,
                1.02,
                f"n={total}",
                ha="center",
                va="bottom",
                fontsize=max(6.0, tick_size - 1.0),
                color="#222222",
            )

        ax.set_title(f"{display_method(source_method)} split events", fontsize=font_size)
        ax.set_xticks(x)
        ax.set_xticklabels(display_methods(CLUSTER_METHODS), fontsize=tick_size)
        ax.tick_params(axis="y", labelsize=tick_size)
        ax.grid(True, axis="y", alpha=0.25, linewidth=0.6)
        ax.set_ylim(0.0, 1.10)
        ax.set_xlabel("Clustering method", fontsize=font_size)

    axes[0].set_ylabel("Fraction of split events", fontsize=font_size)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.03),
        ncol=3,
        frameon=False,
        fontsize=tick_size,
    )
    fig.suptitle(title, fontsize=max(font_size + 1.0, 10.0), y=1.08)
    fig.tight_layout()
    fig.savefig(out_path, format="svg")
    plt.close(fig)


def run(cfg) -> int:
    """
    Figure 13 (SVG):
      Capture of denoising split events by clustering methods.

    For each community/sample, split events are defined for ASV/zOTU as cases where
    one expected species is represented by >1 inferred feature. For each split event,
    OTU/TIC/TAC are evaluated on the same sample/species and classified as:
      - captured_as_one: exactly one clustered feature for that species
      - still_split: >1 clustered features for that species
      - not_captured: no clustered feature for that species
    """
    log = get_logger()
    outdir = Path(cfg.paths.outdir)

    mf_p = outdir / "manifest.tsv"
    if not mf_p.exists():
        log.info("Figure 13: manifest.tsv not found; run discover first.")
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
    exclude_patterns += list(getattr(cfg.figures, "exclude_community_name_patterns_fig13", []) or [])

    event_rows: list[dict[str, object]] = []
    capture_rows: list[dict[str, object]] = []
    event_id = 0

    for _, r in mf.iterrows():
        cid = str(r["community_id"])
        if cid in skipped or _is_excluded(cid, exclude_patterns):
            continue

        method_counts: dict[str, pd.DataFrame] = {}
        for method in ALL_METHODS:
            method_counts[method] = _load_method_species_feature_counts(outdir, cid, method)

        for source_method in SOURCE_METHODS:
            src = method_counts.get(source_method)
            if src is None or src.empty:
                continue
            splits = src[src["n_features"].astype(int) > 1].copy()
            if splits.empty:
                continue

            for _, ev in splits.iterrows():
                sample_id = str(ev["sample_id"])
                species = str(ev["target_species"])
                n_source = int(ev["n_features"])
                sample_key = f"{cid}|{sample_id}"

                event_id += 1
                event_rows.append(
                    {
                        "event_id": event_id,
                        "community_id": cid,
                        "community_display": display_community_id(cid),
                        "sample_id": sample_id,
                        "sample_key": sample_key,
                        "sample_key_display": display_sample_key(sample_key),
                        "source_method": source_method,
                        "source_method_display": display_method(source_method),
                        "target_species": species,
                        "n_source_features": n_source,
                    }
                )

                for clustering_method in CLUSTER_METHODS:
                    cmp_df = method_counts.get(clustering_method)
                    n_cluster = 0
                    if cmp_df is not None and not cmp_df.empty:
                        m = cmp_df[
                            (cmp_df["sample_id"].astype(str) == sample_id)
                            & (cmp_df["target_species"].astype(str) == species)
                        ]
                        if not m.empty:
                            n_cluster = int(m["n_features"].iloc[0])
                    status = _status_from_count(n_cluster)
                    capture_rows.append(
                        {
                            "event_id": event_id,
                            "community_id": cid,
                            "community_display": display_community_id(cid),
                            "sample_id": sample_id,
                            "sample_key": sample_key,
                            "sample_key_display": display_sample_key(sample_key),
                            "source_method": source_method,
                            "source_method_display": display_method(source_method),
                            "target_species": species,
                            "n_source_features": n_source,
                            "clustering_method": clustering_method,
                            "clustering_method_display": display_method(clustering_method),
                            "n_clustering_features": n_cluster,
                            "capture_status": status,
                            "capture_status_label": STATUS_LABELS.get(status, status),
                        }
                    )

    events_df = pd.DataFrame(event_rows)
    capture_df = pd.DataFrame(capture_rows)
    summary_df = _summary_table(capture_df) if not capture_df.empty else _summary_table(
        pd.DataFrame(columns=["event_id", "source_method", "clustering_method", "capture_status"])
    )

    out_dir = figures_root(outdir) / "fig13_split_capture_by_clustering"
    out_dir.mkdir(parents=True, exist_ok=True)
    events_df.to_csv(out_dir / "split_events_denoising.tsv", sep="\t", index=False)
    capture_df.to_csv(out_dir / "split_capture_by_clustering.tsv", sep="\t", index=False)
    summary_df.to_csv(out_dir / "split_capture_summary.tsv", sep="\t", index=False)

    title = getattr(cfg.figures, "fig13_title", None) or "Capture of Denoising Split Events by Clustering Methods"
    font_size = float(getattr(cfg.figures, "fig13_font_size", 11.0))
    _plot_summary(
        summary_df=summary_df,
        out_path=out_dir / "figure13_split_capture_by_clustering.svg",
        title=title,
        font_size=font_size,
    )

    total_events = int(events_df["event_id"].nunique()) if not events_df.empty else 0
    log.info(f"Wrote Figure 13 to {out_dir} (split events={total_events})")
    return 0

