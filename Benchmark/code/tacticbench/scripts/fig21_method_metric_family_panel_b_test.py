from __future__ import annotations

import argparse
from itertools import combinations, product
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

from tacticbench.config import load_config
from tacticbench.steps.plot_fig11 import (
    METHOD_COLORS,
    METHOD_ORDER,
    _compute_panel_stats_permutation,
    _draw_pairwise_brackets,
    _normalize_pairwise_test,
    _pvalue_stars,
)
from tacticbench.steps.plot_fig21 import (
    ENDPOINT_GROUP_LABELS,
    PAIR_GROUP_COLORS,
    PAIR_GROUP_ORDER,
    PLANNED_PAIR_GROUP_COMPARISONS,
    SAMPLE_STATUS_LABELS,
    SAMPLE_STATUS_MARKERS,
    _attach_sample_v4_status,
    _fig21_text_scale,
    _metric_family_lists,
    _panel_b_within_contrast_stat_rows,
    _metric_label,
    _normalize_summary_stat,
    _normalize_violin_bw_method,
    _pair_group_full_label,
    _resolve_fig21_panel_font_size,
)
from tacticbench.utils.labels import (
    configure_community_display_mapping_from_manifest,
    display_community_id,
    display_method,
    display_sample_key,
)
from tacticbench.utils.paths import figures_root


OUTPUT_DIR_NAME = "fig21_pair_group_distance_contrast-TEST"
OUTPUT_STEM = "fig21_method_metric_family_panel_b_test"
OUTPUT_STEM_PANEL_F = "fig21_hybrid_method_panel_f_test"
OUTPUT_STEM_PANEL_B_FILTERED = "fig21_panel_b_without_dada2_uparse_test"
FAMILY_ORDER = ["phylogeny_aware", "phylogeny_agnostic"]
PANEL_F_FAMILY_ORDER = ["phylogeny_agnostic", "phylogeny_aware"]
FAMILY_LABELS = {
    "phylogeny_aware": "phylo-aware",
    "phylogeny_agnostic": "phylo-agnostic",
}
PANEL_F_METHOD_GROUP_ORDER = ["OTU", "Hybrid", "zOTU"]
DENOISING_METHODS = ("ASV", "zOTU")
CLUSTERING_METHODS = ("OTU", "TIC", "TAC")
HYBRID_METHODS = ("TIC", "TAC")


def _resolve_outdir(cfg, out_dir_arg: str | None) -> Path:
    if out_dir_arg:
        return Path(out_dir_arg)
    return figures_root(Path(cfg.paths.outdir)) / OUTPUT_DIR_NAME


def _resolve_source_tsv(cfg, source_tsv_arg: str | None) -> Path:
    if source_tsv_arg:
        return Path(source_tsv_arg)
    return figures_root(Path(cfg.paths.outdir)) / "fig21_pair_group_distance_contrast" / "distance_points_filtered.tsv"


def _configured_method_order(cfg, long: pd.DataFrame) -> list[str]:
    configured = list(getattr(cfg.distances, "methods", METHOD_ORDER) or METHOD_ORDER)
    order = [m for m in METHOD_ORDER if m in configured] + [m for m in configured if m not in METHOD_ORDER]
    present = set(long["method"].astype(str))
    return [str(m) for m in order if str(m) in present]


def _configured_metric_order(cfg, long: pd.DataFrame) -> list[str]:
    metrics_cfg = [str(m) for m in (getattr(cfg.figures, "fig21_metrics", []) or []) if str(m).strip()]
    configured = metrics_cfg or [str(m) for m in (getattr(cfg.distances, "metrics", []) or []) if str(m).strip()]
    present = set(long["metric"].astype(str))
    return [m for m in configured if m in present]


def _family_for_metric(metric: str, agnostic_metrics: list[str], aware_metrics: list[str]) -> str | None:
    metric = str(metric)
    if metric in aware_metrics:
        return "phylogeny_aware"
    if metric in agnostic_metrics:
        return "phylogeny_agnostic"
    return None


def _group_id(method: str, family: str) -> str:
    return f"{method}__{family}"


def _display_group(method_group: str) -> str:
    if str(method_group) == "Hybrid":
        return "Hybrid"
    return display_method(method_group)


def _group_color(method_group: str) -> str:
    if str(method_group) == "Hybrid":
        return "#4E9F3D"
    return METHOD_COLORS.get(str(method_group), "#888888")


def _group_label(method: str, family: str) -> str:
    return f"{_display_group(method)} | {FAMILY_LABELS.get(family, family)}"


def _grid_label(group_id: str) -> str:
    method, family = str(group_id).split("__", 1)
    return f"{_display_group(method)}\n{FAMILY_LABELS.get(family, family)}"


def _focused_method_groups(method_order: list[str]) -> list[dict]:
    present = set(method_order)
    groups: list[dict] = []
    if "zOTU" in present:
        groups.append(
            {
                "group_id": "zOTU",
                "group_label": display_method("zOTU"),
                "source_methods": ["zOTU"],
                "color_key": "zOTU",
            }
        )
    if "OTU" in present:
        groups.append(
            {
                "group_id": "OTU",
                "group_label": display_method("OTU"),
                "source_methods": ["OTU"],
                "color_key": "OTU",
            }
        )
    hybrid_sources = [m for m in HYBRID_METHODS if m in present]
    if hybrid_sources:
        groups.append(
            {
                "group_id": "Hybrid",
                "group_label": "Hybrid",
                "source_methods": hybrid_sources,
                "color_key": "Hybrid",
            }
        )
    return groups


def _build_method_family_tables(
    long: pd.DataFrame,
    method_groups: list[dict],
    metric_order: list[str],
    agnostic_metrics: list[str],
    aware_metrics: list[str],
    summary_stat: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str]]:
    required = {"community_id", "sample_id", "method", "metric", "distance", "sample_key"}
    missing = required.difference(set(long.columns))
    if missing:
        raise ValueError(f"Distance table is missing required columns: {sorted(missing)}")

    metric_family_by_metric = {
        metric: _family_for_metric(metric, agnostic_metrics=agnostic_metrics, aware_metrics=aware_metrics)
        for metric in metric_order
    }
    metric_family_by_metric = {m: f for m, f in metric_family_by_metric.items() if f is not None}
    method_to_group: dict[str, str] = {}
    group_label_by_id: dict[str, str] = {}
    for group in method_groups:
        group_id = str(group["group_id"])
        group_label_by_id[group_id] = str(group.get("group_label", group_id))
        for method in group.get("source_methods", []):
            method_to_group[str(method)] = group_id

    endpoint_group_order = [
        _group_id(str(group["group_id"]), family)
        for group in method_groups
        for family in FAMILY_ORDER
        if any(_family_for_metric(metric, agnostic_metrics, aware_metrics) == family for metric in metric_order)
    ]

    sub = long[
        long["method"].astype(str).isin(method_to_group.keys())
        & long["metric"].astype(str).isin(metric_family_by_metric)
    ].copy()
    sub["distance"] = pd.to_numeric(sub["distance"], errors="coerce")
    sub = sub[np.isfinite(sub["distance"])].copy()
    if sub.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), endpoint_group_order

    sub["metric_family"] = sub["metric"].map(metric_family_by_metric)
    sub["method_group"] = sub["method"].map(method_to_group)
    sub["method_group_label"] = sub["method_group"].map(group_label_by_id)
    sub["endpoint_group_id"] = sub.apply(lambda r: _group_id(str(r["method_group"]), str(r["metric_family"])), axis=1)
    sub["endpoint_group_label"] = sub.apply(lambda r: _group_label(str(r["method_group"]), str(r["metric_family"])), axis=1)
    sub["endpoint_id"] = sub["method"].astype(str) + "__" + sub["metric"].astype(str)
    sub["community_label"] = sub["community_id"].astype(str).map(display_community_id)
    sub["sample_label"] = sub["sample_key"].astype(str).map(display_sample_key)
    sub["method_label"] = sub["method"].astype(str).map(display_method)
    sub["metric_label"] = sub["metric"].astype(str).map(_metric_label)
    endpoints = sub[
        [
            "community_id",
            "sample_id",
            "sample_key",
            "community_label",
            "sample_label",
            "method",
            "method_label",
            "method_group",
            "method_group_label",
            "metric",
            "metric_label",
            "metric_family",
            "endpoint_group_id",
            "endpoint_group_label",
            "endpoint_id",
            "distance",
        ]
    ].copy()

    pair_rows: list[dict] = []
    sample_meta = (
        endpoints[["community_id", "sample_id", "sample_key", "community_label", "sample_label"]]
        .drop_duplicates()
        .sort_values(["community_id", "sample_id", "sample_key"])
    )
    for smeta in sample_meta.itertuples(index=False):
        sample_df = endpoints[endpoints["sample_key"].astype(str) == str(smeta.sample_key)].copy()
        grouped = {
            group: sample_df[sample_df["endpoint_group_id"].astype(str) == group].copy()
            for group in endpoint_group_order
        }
        for row_group in endpoint_group_order:
            row_df = grouped.get(row_group, pd.DataFrame())
            if row_df.empty:
                continue
            for col_group in endpoint_group_order:
                col_df = grouped.get(col_group, pd.DataFrame())
                if col_df.empty:
                    continue
                if row_group == col_group:
                    pair_iter = combinations(row_df.to_dict("records"), 2)
                else:
                    pair_iter = product(row_df.to_dict("records"), col_df.to_dict("records"))
                for a, b in pair_iter:
                    signed = float(a["distance"]) - float(b["distance"])
                    pair_rows.append(
                        {
                            "community_id": str(smeta.community_id),
                            "sample_id": str(smeta.sample_id),
                            "sample_key": str(smeta.sample_key),
                            "community_label": str(smeta.community_label),
                            "sample_label": str(smeta.sample_label),
                            "row_endpoint_group_id": row_group,
                            "col_endpoint_group_id": col_group,
                            "row_endpoint_group_label": str(a["endpoint_group_label"]),
                            "col_endpoint_group_label": str(b["endpoint_group_label"]),
                            "endpoint_a": str(a["endpoint_id"]),
                            "endpoint_b": str(b["endpoint_id"]),
                            "method_a": str(a["method"]),
                            "method_b": str(b["method"]),
                            "method_group_a": str(a["method_group"]),
                            "method_group_b": str(b["method_group"]),
                            "metric_a": str(a["metric"]),
                            "metric_b": str(b["metric"]),
                            "metric_family_a": str(a["metric_family"]),
                            "metric_family_b": str(b["metric_family"]),
                            "distance_a": float(a["distance"]),
                            "distance_b": float(b["distance"]),
                            "signed_difference_a_minus_b": signed,
                            "abs_difference": abs(signed),
                        }
                    )

    pairwise = pd.DataFrame(pair_rows)
    if pairwise.empty:
        return endpoints, pairwise, pd.DataFrame(), endpoint_group_order

    summary = (
        pairwise.groupby(
            [
                "community_id",
                "sample_id",
                "sample_key",
                "community_label",
                "sample_label",
                "row_endpoint_group_id",
                "col_endpoint_group_id",
                "row_endpoint_group_label",
                "col_endpoint_group_label",
            ],
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
    stat = _normalize_summary_stat(summary_stat)
    summary["summary_abs_difference"] = (
        summary["mean_abs_difference"] if stat == "mean" else summary["median_abs_difference"]
    )
    summary["summary_stat"] = stat
    return endpoints, pairwise, summary, endpoint_group_order


def _write_cell_summary(sample_summary: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    if sample_summary.empty:
        out = pd.DataFrame()
    else:
        out = (
            sample_summary.groupby(
                [
                    "row_endpoint_group_id",
                    "col_endpoint_group_id",
                    "row_endpoint_group_label",
                    "col_endpoint_group_label",
                ],
                sort=False,
            )
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
    out.to_csv(out_dir / f"{OUTPUT_STEM}_cell_summary.tsv", sep="\t", index=False)
    return out


FILTERED_PANEL_B_DENOISING_METHODS = ("zOTU",)
FILTERED_PANEL_B_CLUSTERING_METHODS = ("TIC", "TAC")
FILTERED_PANEL_B_ENDPOINT_GROUP_ORDER = [
    "denoising_phylogeny_aware",
    "denoising_phylogeny_agnostic",
    "clustering_phylogeny_aware",
    "clustering_phylogeny_agnostic",
]
FILTERED_PANEL_B_ENDPOINT_GROUP_RANK = {
    group: idx for idx, group in enumerate(FILTERED_PANEL_B_ENDPOINT_GROUP_ORDER)
}
FILTERED_PANEL_B_WITHIN_ORDER = [
    "denoising_phylogeny_aware__vs__denoising_phylogeny_aware",
    "clustering_phylogeny_aware__vs__clustering_phylogeny_aware",
    "denoising_phylogeny_agnostic__vs__denoising_phylogeny_agnostic",
    "clustering_phylogeny_agnostic__vs__clustering_phylogeny_agnostic",
]
FILTERED_PANEL_B_CROSS_ORDER = [
    "denoising_phylogeny_aware__vs__denoising_phylogeny_agnostic",
    "clustering_phylogeny_aware__vs__clustering_phylogeny_agnostic",
    "denoising_phylogeny_aware__vs__clustering_phylogeny_aware",
    "denoising_phylogeny_agnostic__vs__clustering_phylogeny_agnostic",
]


def _filtered_panel_b_endpoint_group(
    method: str,
    metric: str,
    agnostic_metrics: list[str],
    aware_metrics: list[str],
) -> str | None:
    method = str(method)
    metric = str(metric)
    if method in FILTERED_PANEL_B_DENOISING_METHODS:
        method_family = "denoising"
    elif method in FILTERED_PANEL_B_CLUSTERING_METHODS:
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


def _filtered_panel_b_pair_group_id(group_a: str, group_b: str) -> str:
    if FILTERED_PANEL_B_ENDPOINT_GROUP_RANK[group_a] <= FILTERED_PANEL_B_ENDPOINT_GROUP_RANK[group_b]:
        return f"{group_a}__vs__{group_b}"
    return f"{group_b}__vs__{group_a}"


def _filtered_family_panel_b_tables(
    long: pd.DataFrame,
    metric_order: list[str],
    agnostic_metrics: list[str],
    aware_metrics: list[str],
    summary_stat: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cols = ["community_id", "sample_id", "sample_key", "method", "metric", "distance"]
    missing = set(cols).difference(set(long.columns))
    if missing:
        raise ValueError(f"Distance table is missing required columns: {sorted(missing)}")

    source_methods = list(FILTERED_PANEL_B_DENOISING_METHODS) + list(FILTERED_PANEL_B_CLUSTERING_METHODS)
    sub = long[
        long["method"].astype(str).isin(source_methods)
        & long["metric"].astype(str).isin(metric_order)
    ].copy()
    sub["distance"] = pd.to_numeric(sub["distance"], errors="coerce")
    sub = sub[np.isfinite(sub["distance"])].copy()
    if sub.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    ordered_endpoints = [
        (method, metric)
        for method in source_methods
        for metric in metric_order
        if _filtered_panel_b_endpoint_group(method, metric, agnostic_metrics, aware_metrics) is not None
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
    for smeta in sample_meta.itertuples(index=False):
        sample_key = str(smeta.sample_key)
        sample_rows = sub[sub["sample_key"].astype(str) == sample_key].copy()
        if sample_rows.empty:
            continue
        pivot = sample_rows.pivot_table(index=["method", "metric"], values="distance", aggfunc="first").reset_index()
        endpoints = order_df.merge(pivot, on=["method", "metric"], how="left").dropna(subset=["distance"]).copy()
        if endpoints.shape[0] < 2:
            continue
        endpoints["community_id"] = str(smeta.community_id)
        endpoints["sample_id"] = str(smeta.sample_id)
        endpoints["sample_key"] = sample_key
        endpoints["community_label"] = endpoints["community_id"].map(display_community_id)
        endpoints["sample_label"] = display_sample_key(sample_key)
        endpoints["method_label"] = endpoints["method"].map(display_method)
        endpoints["metric_label"] = endpoints["metric"].map(_metric_label)
        endpoints["endpoint_group_id"] = endpoints.apply(
            lambda r: _filtered_panel_b_endpoint_group(str(r["method"]), str(r["metric"]), agnostic_metrics, aware_metrics),
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
            pair_group = _filtered_panel_b_pair_group_id(group_a, group_b)
            signed = float(a["distance"]) - float(b["distance"])
            pair_rows.append(
                {
                    "community_id": str(smeta.community_id),
                    "sample_id": str(smeta.sample_id),
                    "sample_key": sample_key,
                    "community_label": display_community_id(str(smeta.community_id)),
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
            [
                "community_id",
                "sample_id",
                "sample_key",
                "community_label",
                "sample_label",
                "pair_group_id",
                "pair_group_label",
            ],
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
    stat = _normalize_summary_stat(summary_stat)
    summary["summary_abs_difference"] = summary["mean_abs_difference"] if stat == "mean" else summary["median_abs_difference"]
    summary["summary_stat"] = stat
    return endpoints_df, pairwise_df, summary


def _filtered_family_panel_b_stats(
    sample_summary: pd.DataFrame,
    out_dir: Path,
    stats_enabled: bool,
    n_perm: int,
    seed: int,
    alpha: float,
    pairwise_test: str,
    multiple_test_correction: bool,
) -> tuple[dict | None, list[dict]]:
    if not stats_enabled or sample_summary.empty:
        pd.DataFrame().to_csv(out_dir / f"{OUTPUT_STEM_PANEL_B_FILTERED}_stats_omnibus.tsv", sep="\t", index=False)
        pd.DataFrame().to_csv(out_dir / f"{OUTPUT_STEM_PANEL_B_FILTERED}_stats_pairwise.tsv", sep="\t", index=False)
        pd.DataFrame().to_csv(out_dir / f"{OUTPUT_STEM_PANEL_B_FILTERED}_stats_planned.tsv", sep="\t", index=False)
        return None, []

    methods = [g for g in PAIR_GROUP_ORDER if g in set(sample_summary["pair_group_id"].astype(str))]
    stats_df = sample_summary.rename(columns={"pair_group_id": "method", "summary_abs_difference": "distance"}).copy()
    st = _compute_panel_stats_permutation(
        df=stats_df,
        method_order=methods,
        n_perm=n_perm,
        seed=seed + 52101,
        alpha=alpha,
        pairwise_test=pairwise_test,
        multiple_test_correction=multiple_test_correction,
    )
    pd.DataFrame(
        [
            {
                "metric": "filtered_family_panel_b",
                "label": "Figure 21 Panel B mimic excluding DADA2 and UPARSE",
                "stats_model": st.get("stats_model", "paired_permutation"),
                "pairwise_test": pairwise_test,
                "multiple_test_correction": bool(st.get("multiple_test_correction", multiple_test_correction)),
                "n_blocks": st.get("n_blocks", 0),
                "n_pair_groups": len(methods),
                "omnibus_stat": st.get("omnibus_stat", np.nan),
                "omnibus_p": st.get("omnibus_p", np.nan),
                "omnibus_significant_alpha": bool(
                    np.isfinite(float(st.get("omnibus_p", np.nan)))
                    and float(st.get("omnibus_p", np.nan)) < alpha
                ),
                "alpha": alpha,
                "n_permutations": n_perm,
                "denoising_methods": ",".join(FILTERED_PANEL_B_DENOISING_METHODS),
                "clustering_methods": ",".join(FILTERED_PANEL_B_CLUSTERING_METHODS),
            }
        ]
    ).to_csv(out_dir / f"{OUTPUT_STEM_PANEL_B_FILTERED}_stats_omnibus.tsv", sep="\t", index=False)

    pairwise_rows: list[dict] = []
    for row in st.get("pairwise", []):
        method_a = str(row.get("method_a", ""))
        method_b = str(row.get("method_b", ""))
        rec = dict(row)
        rec.update(
            {
                "metric": "filtered_family_panel_b",
                "method_a_label": _pair_group_full_label(method_a),
                "method_b_label": _pair_group_full_label(method_b),
                "alpha": alpha,
                "n_permutations": n_perm,
                "denoising_methods": ",".join(FILTERED_PANEL_B_DENOISING_METHODS),
                "clustering_methods": ",".join(FILTERED_PANEL_B_CLUSTERING_METHODS),
            }
        )
        pairwise_rows.append(rec)
    pairwise_df = pd.DataFrame(pairwise_rows)
    pairwise_df.to_csv(out_dir / f"{OUTPUT_STEM_PANEL_B_FILTERED}_stats_pairwise.tsv", sep="\t", index=False)

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
    pd.DataFrame(planned_rows).to_csv(out_dir / f"{OUTPUT_STEM_PANEL_B_FILTERED}_stats_planned.tsv", sep="\t", index=False)
    return st, planned_rows


def _filtered_panel_b_tick_label(pair_group_id: str) -> str:
    label = _pair_group_full_label(pair_group_id)
    label = label.replace(" within-group", "\nwithin")
    label = label.replace(" vs ", "\nvs\n")
    return label


def _filtered_panel_b_bracket_rows(
    stats: dict | None,
    planned_rows: list[dict],
    order: list[str],
    stats_show_pairwise: bool,
    show_within_contrast_stats: bool,
) -> list[dict]:
    visible_set = set(order)
    rows: list[dict] = []
    if stats_show_pairwise:
        rows.extend(
            r
            for r in planned_rows
            if str(r.get("method_a", "")) in visible_set and str(r.get("method_b", "")) in visible_set
        )
    if show_within_contrast_stats:
        rows.extend(_panel_b_within_contrast_stat_rows(stats, visible_set))
    return rows


def _draw_filtered_family_panel_b_subplot(
    ax,
    sample_summary: pd.DataFrame,
    stats: dict | None,
    planned_rows: list[dict],
    order: list[str],
    title: str,
    font_size: float,
    violin_bw_method,
    plot_type: str,
    point_size: float,
    point_opacity: float,
    point_jitter: float,
    stats_enabled: bool,
    stats_alpha: float,
    stats_show_pairwise: bool,
    show_within_contrast_stats: bool,
    stat_title: str,
    show_legend: bool = False,
) -> list[str]:
    fs = _fig21_text_scale(font_size)
    panel_df = sample_summary.copy()
    if panel_df.empty:
        ax.axis("off")
        return []
    if "v4_reference_status" not in panel_df.columns:
        panel_df["v4_reference_status"] = "single_v4_reference"
    panel_df["v4_reference_status"] = panel_df["v4_reference_status"].fillna("single_v4_reference").astype(str)
    panel_df["v4_reference_status"] = panel_df["v4_reference_status"].where(
        panel_df["v4_reference_status"].isin(SAMPLE_STATUS_MARKERS),
        "single_v4_reference",
    )

    present_pair_groups = set(panel_df["pair_group_id"].astype(str))
    order = [g for g in order if g in present_pair_groups]
    if not order:
        ax.axis("off")
        return []

    positions = np.arange(1, len(order) + 1, dtype=float)
    plot_df = panel_df[panel_df["pair_group_id"].astype(str).isin(order)].copy()
    plot_df["summary_abs_difference"] = pd.to_numeric(plot_df["summary_abs_difference"], errors="coerce")
    plot_df = plot_df.dropna(subset=["summary_abs_difference"])
    if plot_df.empty:
        ax.axis("off")
        return []

    group_data: list[tuple[int, str, pd.DataFrame, np.ndarray]] = []
    for idx, group in enumerate(order):
        group_df = plot_df[plot_df["pair_group_id"].astype(str) == group].copy()
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
                boxprops=dict(facecolor=PAIR_GROUP_COLORS.get(group, "#999999"), alpha=0.42, edgecolor="#333333"),
                medianprops=dict(color="#111111", linewidth=1.15),
                whiskerprops=dict(color="#333333", linewidth=0.85),
                capprops=dict(color="#333333", linewidth=0.85),
            )
    else:
        violin_vals = [vals for _, _, _, vals in group_data if vals.size >= 2]
        violin_pos = [float(positions[idx]) for idx, _, _, vals in group_data if vals.size >= 2]
        violin_groups = [group for _, group, _, vals in group_data if vals.size >= 2]
        if violin_vals:
            vp = ax.violinplot(
                violin_vals,
                positions=violin_pos,
                widths=0.70,
                showmeans=False,
                showmedians=True,
                showextrema=False,
                bw_method=violin_bw_method,
            )
            for body, group in zip(vp["bodies"], violin_groups):
                body.set_facecolor(PAIR_GROUP_COLORS.get(group, "#999999"))
                body.set_edgecolor("#333333")
                body.set_alpha(0.48)
            if "cmedians" in vp:
                vp["cmedians"].set_color("#111111")
                vp["cmedians"].set_linewidth(1.1)

    sample_keys = sorted(plot_df["sample_key"].astype(str).unique().tolist())
    rng = np.random.default_rng(62891)
    jitter_by_sample = {
        sample_key: float(rng.normal(0.0, max(0.0, float(point_jitter))))
        for sample_key in sample_keys
    }
    scatter_size = max(1.0, float(point_size))
    point_rows: list[dict[str, object]] = []
    for idx, group, group_df, _vals in group_data:
        x_center = float(positions[idx])
        for row in group_df.itertuples(index=False):
            sample_key = str(getattr(row, "sample_key"))
            y = float(getattr(row, "summary_abs_difference"))
            x = x_center + float(jitter_by_sample.get(sample_key, 0.0))
            point_rows.append(
                {
                    "sample_key": sample_key,
                    "pair_group_id": group,
                    "x": x,
                    "y": y,
                    "status": str(getattr(row, "v4_reference_status", "single_v4_reference")),
                }
            )
    point_df = pd.DataFrame(point_rows)

    # Draw sample trajectories first so dots remain readable.
    group_rank = {group: i for i, group in enumerate(order)}
    for sample_key, sdf in point_df.groupby("sample_key", sort=False):
        sdf = sdf.sort_values("pair_group_id", key=lambda s: s.map(group_rank)).copy()
        if sdf.shape[0] < 2:
            continue
        ax.plot(
            sdf["x"].to_numpy(dtype=float),
            sdf["y"].to_numpy(dtype=float),
            color="#3A3A3A",
            alpha=0.16,
            linewidth=0.75,
            zorder=2,
        )

    for group in order:
        gdf = point_df[point_df["pair_group_id"].astype(str) == group].copy()
        if gdf.empty:
            continue
        color = PAIR_GROUP_COLORS.get(group, "#666666")
        for status_key, marker in SAMPLE_STATUS_MARKERS.items():
            sub = gdf[gdf["status"].astype(str) == status_key]
            if sub.empty:
                continue
            ax.scatter(
                sub["x"].to_numpy(dtype=float),
                sub["y"].to_numpy(dtype=float),
                s=scatter_size,
                marker=marker,
                color=color,
                alpha=max(0.0, min(1.0, float(point_opacity))),
                edgecolors="#111111",
                linewidths=0.45,
                zorder=4,
            )

    finite_vals = pd.to_numeric(plot_df["summary_abs_difference"], errors="coerce").dropna().to_numpy(dtype=float)
    if finite_vals.size:
        y_min = float(np.min(finite_vals))
        y_max = float(np.max(finite_vals))
        y_range = y_max - y_min
        if y_range <= 0.0:
            y_range = max(0.05, abs(y_max) * 0.25)
        ax.set_ylim(max(0.0, y_min - 0.08 * y_range), y_max + 0.48 * y_range)

    ax.set_xticks(positions)
    ax.set_xticklabels(
        [_filtered_panel_b_tick_label(g) for g in order],
        rotation=30.0,
        ha="right",
        va="top",
        rotation_mode="anchor",
        fontsize=fs["tick"],
    )
    for tick in ax.get_xticklabels():
        tick.set_linespacing(0.9)
    ax.set_xlabel("Method/metric contrast group", fontsize=fs["base"])
    ax.set_ylabel(f"{stat_title} absolute pairwise difference", fontsize=fs["base"])
    ax.tick_params(axis="y", labelsize=fs["tick"])
    ax.grid(True, axis="y", alpha=0.25, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.set_title(title, fontsize=fs["base"])

    if stats_enabled and isinstance(stats, dict):
        bracket_rows = _filtered_panel_b_bracket_rows(
            stats=stats,
            planned_rows=planned_rows,
            order=order,
            stats_show_pairwise=stats_show_pairwise,
            show_within_contrast_stats=show_within_contrast_stats,
        )
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

    if show_legend:
        status_handles = [
            Line2D(
                [0],
                [0],
                marker=SAMPLE_STATUS_MARKERS["multiple_v4_reference"],
                linestyle="None",
                markerfacecolor="#666666",
                markeredgecolor="#111111",
                markersize=8.5,
                label=SAMPLE_STATUS_LABELS["multiple_v4_reference"],
            ),
            Line2D(
                [0],
                [0],
                marker=SAMPLE_STATUS_MARKERS["single_v4_reference"],
                linestyle="None",
                markerfacecolor="#666666",
                markeredgecolor="#111111",
                markersize=8.5,
                label=SAMPLE_STATUS_LABELS["single_v4_reference"],
            ),
            Line2D([0], [0], color="#3A3A3A", alpha=0.24, linewidth=1.0, label="Same sample"),
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


def _draw_filtered_family_panel_b(
    sample_summary: pd.DataFrame,
    stats: dict | None,
    planned_rows: list[dict],
    out_dir: Path,
    font_size: float,
    violin_bw_method,
    plot_type: str,
    point_size: float,
    point_opacity: float,
    point_jitter: float,
    stats_enabled: bool,
    stats_alpha: float,
    stats_show_pairwise: bool,
    show_within_contrast_stats: bool,
    pair_groups_to_show: list[str],
) -> Path:
    if "summary_stat" in sample_summary.columns and not sample_summary.empty:
        stat = _normalize_summary_stat(str(sample_summary["summary_stat"].dropna().iloc[0]))
    else:
        stat = "median"
    stat_title = "Mean" if stat == "mean" else "Median"

    requested = [str(g) for g in pair_groups_to_show if str(g).strip()]
    if not requested:
        requested = list(FILTERED_PANEL_B_WITHIN_ORDER) + list(FILTERED_PANEL_B_CROSS_ORDER)
    requested_set = set(requested)
    within_order = [g for g in FILTERED_PANEL_B_WITHIN_ORDER if g in requested_set]
    cross_order = [g for g in FILTERED_PANEL_B_CROSS_ORDER if g in requested_set]

    standalone_specs = [
        (
            "within",
            within_order,
            "Within-group contrast distributions",
            out_dir / f"{OUTPUT_STEM_PANEL_B_FILTERED}_within.svg",
        ),
        (
            "cross",
            cross_order,
            "Cross-group contrast distributions",
            out_dir / f"{OUTPUT_STEM_PANEL_B_FILTERED}_cross.svg",
        ),
    ]
    for _kind, order, title, path in standalone_specs:
        fig_single, ax_single = plt.subplots(figsize=(11.6, 6.8))
        _draw_filtered_family_panel_b_subplot(
            ax=ax_single,
            sample_summary=sample_summary,
            stats=stats,
            planned_rows=planned_rows,
            order=order,
            title=title,
            font_size=font_size,
            violin_bw_method=violin_bw_method,
            plot_type=plot_type,
            point_size=point_size,
            point_opacity=point_opacity,
            point_jitter=point_jitter,
            stats_enabled=stats_enabled,
            stats_alpha=stats_alpha,
            stats_show_pairwise=stats_show_pairwise,
            show_within_contrast_stats=show_within_contrast_stats,
            stat_title=stat_title,
            show_legend=True,
        )
        fig_single.tight_layout()
        fig_single.savefig(path, format="svg", dpi=300, bbox_inches="tight")
        plt.close(fig_single)

    fig, axes = plt.subplots(2, 1, figsize=(14.0, 12.4), sharey=True)
    _draw_filtered_family_panel_b_subplot(
        ax=axes[0],
        sample_summary=sample_summary,
        stats=stats,
        planned_rows=planned_rows,
        order=within_order,
        title="Within-group contrast distributions",
        font_size=font_size,
        violin_bw_method=violin_bw_method,
        plot_type=plot_type,
        point_size=point_size,
        point_opacity=point_opacity,
        point_jitter=point_jitter,
        stats_enabled=stats_enabled,
        stats_alpha=stats_alpha,
        stats_show_pairwise=stats_show_pairwise,
        show_within_contrast_stats=show_within_contrast_stats,
        stat_title=stat_title,
        show_legend=True,
    )
    _draw_filtered_family_panel_b_subplot(
        ax=axes[1],
        sample_summary=sample_summary,
        stats=stats,
        planned_rows=planned_rows,
        order=cross_order,
        title="Cross-group contrast distributions",
        font_size=font_size,
        violin_bw_method=violin_bw_method,
        plot_type=plot_type,
        point_size=point_size,
        point_opacity=point_opacity,
        point_jitter=point_jitter,
        stats_enabled=stats_enabled,
        stats_alpha=stats_alpha,
        stats_show_pairwise=stats_show_pairwise,
        show_within_contrast_stats=show_within_contrast_stats,
        stat_title=stat_title,
        show_legend=False,
    )
    fig.suptitle("Panel B mimic excluding DADA2 and UPARSE", fontsize=max(10.0, font_size + 1.0), y=0.995)
    fig.tight_layout(rect=(0.01, 0.01, 0.995, 0.975))
    out_path = out_dir / f"{OUTPUT_STEM_PANEL_B_FILTERED}.svg"
    fig.savefig(out_path, format="svg", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _draw_violin_grid(
    sample_summary: pd.DataFrame,
    endpoint_group_order: list[str],
    out_dir: Path,
    font_size: float,
    violin_bw_method,
    summary_stat: str,
) -> Path:
    n = len(endpoint_group_order)
    if n == 0:
        raise ValueError("No endpoint groups are available for plotting.")

    fs = _fig21_text_scale(font_size)
    all_vals = pd.to_numeric(sample_summary["summary_abs_difference"], errors="coerce").dropna().to_numpy(dtype=float)
    if all_vals.size == 0:
        raise ValueError("No finite sample-level summary values are available for plotting.")
    y_min = max(0.0, float(np.nanmin(all_vals)) - 0.04 * max(0.01, float(np.nanmax(all_vals) - np.nanmin(all_vals))))
    y_max = float(np.nanmax(all_vals))
    pad = 0.08 * max(0.01, y_max - y_min)
    y_max += pad

    fig_width = max(13.5, 1.18 * n + 3.6)
    fig_height = max(13.5, 1.05 * n + 3.7)
    fig, axes = plt.subplots(n, n, figsize=(fig_width, fig_height), sharey=True)
    if n == 1:
        axes = np.asarray([[axes]])

    row_method_colors = {}
    for group in endpoint_group_order:
        method = group.split("__", 1)[0]
        row_method_colors[group] = _group_color(method)

    for i, row_group in enumerate(endpoint_group_order):
        for j, col_group in enumerate(endpoint_group_order):
            ax = axes[i, j]
            cell = sample_summary[
                (sample_summary["row_endpoint_group_id"].astype(str) == row_group)
                & (sample_summary["col_endpoint_group_id"].astype(str) == col_group)
            ].copy()
            vals = pd.to_numeric(cell["summary_abs_difference"], errors="coerce").dropna().to_numpy(dtype=float)
            ax.set_xlim(-0.55, 0.55)
            ax.set_ylim(y_min, y_max)
            ax.set_xticks([])
            if j == 0:
                ax.tick_params(axis="y", labelsize=fs["tick"])
                ax.set_ylabel(_grid_label(row_group), rotation=0, ha="right", va="center", labelpad=50, fontsize=fs["tick"])
            else:
                ax.tick_params(axis="y", labelleft=False, length=0)

            if i == 0:
                ax.set_title(_grid_label(col_group), fontsize=fs["tick"], pad=10)

            if vals.size >= 2:
                vp = ax.violinplot(
                    [vals],
                    positions=[0.0],
                    widths=0.78,
                    showmeans=False,
                    showmedians=True,
                    showextrema=False,
                    bw_method=violin_bw_method,
                )
                color = row_method_colors.get(row_group, "#888888")
                for body in vp["bodies"]:
                    body.set_facecolor(color)
                    body.set_edgecolor("#333333")
                    body.set_alpha(0.52)
                if "cmedians" in vp:
                    vp["cmedians"].set_color("#111111")
                    vp["cmedians"].set_linewidth(1.05)
            elif vals.size == 1:
                ax.scatter([0.0], vals, s=16, color=row_method_colors.get(row_group, "#888888"), edgecolors="#111111", linewidths=0.4)
            else:
                ax.text(0.5, 0.5, "NA", transform=ax.transAxes, ha="center", va="center", fontsize=fs["pair_label"], color="#777777")

            if vals.size:
                ax.text(
                    0.96,
                    0.91,
                    f"n={vals.size}",
                    transform=ax.transAxes,
                    ha="right",
                    va="top",
                    fontsize=max(5.0, fs["pair_label"] - 1.0),
                    color="#333333",
                )
            ax.grid(True, axis="y", linewidth=0.35, alpha=0.20)
            for spine in ax.spines.values():
                spine.set_visible(True)
                spine.set_linewidth(0.45)
                spine.set_color("#CFCFCF")

    stat_label = "mean" if _normalize_summary_stat(summary_stat) == "mean" else "median"
    fig.suptitle(
        "Hybrid-focused metric-family contrast distributions",
        fontsize=fs["base"],
        y=0.985,
    )
    fig.text(
        0.5,
        0.025,
        "Column endpoint group: analysis method x distance-metric family",
        ha="center",
        va="center",
        fontsize=fs["base"],
    )
    fig.text(
        0.018,
        0.5,
        f"Sample-level {stat_label} absolute pairwise difference",
        ha="center",
        va="center",
        rotation=90,
        fontsize=fs["base"],
    )
    fig.subplots_adjust(left=0.16, right=0.995, bottom=0.06, top=0.90, wspace=0.06, hspace=0.08)

    out_path = out_dir / f"{OUTPUT_STEM}.svg"
    fig.savefig(out_path, format="svg", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _comparison_id(parts: list[str] | tuple[str, ...]) -> str:
    return "__".join(str(p).replace(" ", "_").replace("|", "_").replace("-", "_") for p in parts)


def _method_family(method: str) -> str:
    if method == "Hybrid":
        return "hybrid"
    if method in DENOISING_METHODS:
        return "denoising"
    if method in CLUSTERING_METHODS:
        return "clustering"
    return "other"


def _available_methods(method_order: list[str], methods: tuple[str, ...]) -> list[str]:
    method_set = set(methods)
    return [m for m in method_order if m in method_set]


def _facet_specs(method_order: list[str]) -> list[dict]:
    focused = [m for m in ["zOTU", "OTU", "Hybrid"] if m in set(method_order)]
    hybrid = "Hybrid" if "Hybrid" in focused else None
    hybrid_targets = [m for m in ["zOTU", "OTU"] if m in focused and hybrid]
    facets: list[dict] = []

    for family, family_short, facet_label in [
        ("phylogeny_aware", "aware", "Within phylo-aware method groups"),
        ("phylogeny_agnostic", "agnostic", "Within phylo-agnostic method groups"),
    ]:
        within_specs: list[dict] = []
        for method in focused:
            within_specs.append(
                {
                    "comparison_id": _comparison_id(["within", method, family_short]),
                    "row_group": _group_id(method, family),
                    "col_group": _group_id(method, family),
                    "comparison_label": f"{_display_group(method)}\nwithin",
                    "primary_method": method,
                }
            )
        planned_pairs = [
            (_comparison_id(["within", target, family_short]), _comparison_id(["within", hybrid, family_short]))
            for target in hybrid_targets
        ]
        facets.append(
            {
                "facet_id": f"within_endpoint_groups_{family}",
                "facet_label": facet_label,
                "specs": within_specs,
                "planned_pairs": planned_pairs,
                "annotation_mode": "planned",
            }
        )

    family_gap_specs: list[dict] = []
    for method in focused:
        family_gap_specs.append(
            {
                "comparison_id": _comparison_id(["family_gap", method]),
                "row_group": _group_id(method, "phylogeny_aware"),
                "col_group": _group_id(method, "phylogeny_agnostic"),
                "comparison_label": f"{_display_group(method)}\naware vs agnostic",
                "primary_method": method,
            }
        )
    family_gap_planned = [(_comparison_id(["family_gap", target]), _comparison_id(["family_gap", hybrid])) for target in hybrid_targets]
    facets.append(
        {
            "facet_id": "same_method_metric_family_gap",
            "facet_label": "Same method: phylo-aware vs phylo-agnostic",
            "specs": family_gap_specs,
            "planned_pairs": family_gap_planned,
            "annotation_mode": "planned",
        }
    )

    for family, facet_label in [
        ("phylogeny_aware", "Hybrid contrast: phylo-aware metrics"),
        ("phylogeny_agnostic", "Hybrid contrast: phylo-agnostic metrics"),
    ]:
        specs: list[dict] = []
        for target in hybrid_targets:
            specs.append(
                {
                    "comparison_id": _comparison_id(["hybrid_contrast", family, target]),
                    "row_group": _group_id(target, family),
                    "col_group": _group_id(hybrid, family),
                    "comparison_label": f"{_display_group(target)}\nvs\nHybrid",
                    "primary_method": target,
                }
            )
        planned_pairs = []
        if len(specs) >= 2:
            planned_pairs = [(str(specs[0]["comparison_id"]), str(specs[1]["comparison_id"]))]
        facets.append(
            {
                "facet_id": f"hybrid_contrast_{family}",
                "facet_label": facet_label,
                "specs": specs,
                "planned_pairs": planned_pairs,
                "annotation_mode": "planned",
            }
        )
    return facets


def _select_method_family_panel_b_data(
    sample_summary: pd.DataFrame,
    method_order: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame, list[dict]]:
    facets = _facet_specs(method_order)
    plot_rows: list[pd.DataFrame] = []
    spec_rows: list[dict] = []

    for facet_index, facet in enumerate(facets, start=1):
        for comparison_index, spec in enumerate(facet["specs"], start=1):
            row_group = str(spec["row_group"])
            col_group = str(spec["col_group"])
            sub = sample_summary[
                (sample_summary["row_endpoint_group_id"].astype(str) == row_group)
                & (sample_summary["col_endpoint_group_id"].astype(str) == col_group)
            ].copy()
            if sub.empty:
                continue
            sub["facet_id"] = str(facet["facet_id"])
            sub["facet_label"] = str(facet["facet_label"])
            sub["facet_order"] = facet_index
            sub["comparison_id"] = str(spec["comparison_id"])
            sub["comparison_label"] = str(spec["comparison_label"])
            sub["comparison_order"] = comparison_index
            sub["primary_method"] = str(spec.get("primary_method", ""))
            sub["primary_method_family"] = sub["primary_method"].map(_method_family)
            plot_rows.append(sub)
            spec_rows.append(
                {
                    "facet_id": str(facet["facet_id"]),
                    "facet_label": str(facet["facet_label"]),
                    "facet_order": facet_index,
                    "comparison_id": str(spec["comparison_id"]),
                    "comparison_label": str(spec["comparison_label"]).replace("\n", " "),
                    "comparison_order": comparison_index,
                    "row_endpoint_group_id": row_group,
                    "col_endpoint_group_id": col_group,
                    "primary_method": str(spec.get("primary_method", "")),
                    "primary_method_label": display_method(str(spec.get("primary_method", ""))),
                    "annotation_mode": str(facet["annotation_mode"]),
                }
            )

    plot_df = pd.concat(plot_rows, ignore_index=True) if plot_rows else pd.DataFrame()
    spec_df = pd.DataFrame(spec_rows)
    return plot_df, spec_df, facets


def _stats_for_method_family_panel_b(
    plot_df: pd.DataFrame,
    facets: list[dict],
    out_dir: Path,
    n_perm: int,
    seed: int,
    alpha: float,
    pairwise_test: str,
    multiple_test_correction: bool,
) -> tuple[dict[str, dict], list[dict]]:
    stats_by_facet: dict[str, dict] = {}
    omnibus_rows: list[dict] = []
    pairwise_rows: list[dict] = []
    annotated_rows: list[dict] = []

    label_map = (
        plot_df[["comparison_id", "comparison_label", "facet_id", "facet_label"]]
        .drop_duplicates()
        .set_index("comparison_id")
        .to_dict(orient="index")
        if not plot_df.empty
        else {}
    )

    for facet_index, facet in enumerate(facets, start=1):
        facet_id = str(facet["facet_id"])
        sub = plot_df[plot_df["facet_id"].astype(str) == facet_id].copy()
        order = [str(spec["comparison_id"]) for spec in facet["specs"] if str(spec["comparison_id"]) in set(sub["comparison_id"].astype(str))]
        if sub.empty or len(order) < 2:
            continue
        stats_df = sub.rename(columns={"comparison_id": "method", "summary_abs_difference": "distance"}).copy()
        st = _compute_panel_stats_permutation(
            df=stats_df,
            method_order=order,
            n_perm=n_perm,
            seed=seed + 31000 + facet_index,
            alpha=alpha,
            pairwise_test=pairwise_test,
            multiple_test_correction=multiple_test_correction,
        )
        stats_by_facet[facet_id] = st
        omnibus_rows.append(
            {
                "facet_id": facet_id,
                "facet_label": str(facet["facet_label"]),
                "stats_model": st.get("stats_model", "paired_permutation"),
                "pairwise_test": pairwise_test,
                "multiple_test_correction": bool(st.get("multiple_test_correction", multiple_test_correction)),
                "n_blocks": st.get("n_blocks", 0),
                "n_comparisons": len(order),
                "omnibus_stat": st.get("omnibus_stat", np.nan),
                "omnibus_p": st.get("omnibus_p", np.nan),
                "omnibus_significant_alpha": bool(
                    np.isfinite(float(st.get("omnibus_p", np.nan)))
                    and float(st.get("omnibus_p", np.nan)) < alpha
                ),
                "alpha": alpha,
                "n_permutations": n_perm,
            }
        )
        for row in st.get("pairwise", []):
            rec = dict(row)
            method_a = str(rec.get("method_a", ""))
            method_b = str(rec.get("method_b", ""))
            rec.update(
                {
                    "facet_id": facet_id,
                    "facet_label": str(facet["facet_label"]),
                    "method_a_label": str(label_map.get(method_a, {}).get("comparison_label", method_a)).replace("\n", " "),
                    "method_b_label": str(label_map.get(method_b, {}).get("comparison_label", method_b)).replace("\n", " "),
                    "alpha": alpha,
                    "n_permutations": n_perm,
                }
            )
            pairwise_rows.append(rec)

        planned_pairs = {(str(a), str(b)) for a, b in facet.get("planned_pairs", [])}
        planned_pairs |= {(b, a) for a, b in list(planned_pairs)}
        for row in st.get("pairwise", []):
            method_a = str(row.get("method_a", ""))
            method_b = str(row.get("method_b", ""))
            keep = False
            mode = str(facet.get("annotation_mode", "significant"))
            if mode == "planned":
                keep = (method_a, method_b) in planned_pairs
            else:
                keep = bool(row.get("significant_alpha", False) or row.get("significant_bh", False))
            if not keep:
                continue
            ann = dict(row)
            ann.update(
                {
                    "facet_id": facet_id,
                    "facet_label": str(facet["facet_label"]),
                    "method_a_label": str(label_map.get(method_a, {}).get("comparison_label", method_a)).replace("\n", " "),
                    "method_b_label": str(label_map.get(method_b, {}).get("comparison_label", method_b)).replace("\n", " "),
                    "annotation_mode": mode,
                    "alpha": alpha,
                    "n_permutations": n_perm,
                }
            )
            annotated_rows.append(ann)

    pd.DataFrame(omnibus_rows).to_csv(out_dir / f"{OUTPUT_STEM}_faceted_stats_omnibus.tsv", sep="\t", index=False)
    pd.DataFrame(pairwise_rows).to_csv(out_dir / f"{OUTPUT_STEM}_faceted_stats_pairwise.tsv", sep="\t", index=False)
    pd.DataFrame(annotated_rows).to_csv(out_dir / f"{OUTPUT_STEM}_faceted_stats_annotated.tsv", sep="\t", index=False)
    return stats_by_facet, annotated_rows


def _draw_panel_b_style_rows(
    plot_df: pd.DataFrame,
    facets: list[dict],
    annotated_rows: list[dict],
    out_dir: Path,
    font_size: float,
    violin_bw_method,
    summary_stat: str,
    point_size: float,
    point_opacity: float,
    point_jitter: float,
    stats_alpha: float,
) -> Path:
    if plot_df.empty:
        raise ValueError("No selected method-specific Figure 21 panel B data are available for plotting.")

    fs = _fig21_text_scale(font_size)
    visible_facets = [facet for facet in facets if str(facet["facet_id"]) in set(plot_df["facet_id"].astype(str))]
    n_facets = len(visible_facets)
    if n_facets <= 5:
        fig = plt.figure(figsize=(18.0, 11.8))
        gs = fig.add_gridspec(2, 6)
        slots = [
            (0, slice(0, 2)),
            (0, slice(2, 4)),
            (0, slice(4, 6)),
            (1, slice(0, 3)),
            (1, slice(3, 6)),
        ]
        axes = np.asarray([fig.add_subplot(gs[row, col]) for row, col in slots[:n_facets]])
    else:
        fig_height = max(10.0, 3.25 * n_facets)
        fig, axes = plt.subplots(n_facets, 1, figsize=(17.0, fig_height), sharey=False)
        if n_facets == 1:
            axes = np.asarray([axes])

    rng = np.random.default_rng(82117)
    axis_index = 0
    stat_label = "mean" if _normalize_summary_stat(summary_stat) == "mean" else "median"
    for facet in visible_facets:
        facet_id = str(facet["facet_id"])
        sub = plot_df[plot_df["facet_id"].astype(str) == facet_id].copy()
        if sub.empty:
            continue
        ax = axes[axis_index]
        axis_index += 1
        order = [
            str(spec["comparison_id"])
            for spec in facet["specs"]
            if str(spec["comparison_id"]) in set(sub["comparison_id"].astype(str))
        ]
        positions = np.arange(1, len(order) + 1, dtype=float)
        label_map = (
            sub[["comparison_id", "comparison_label", "primary_method"]]
            .drop_duplicates()
            .set_index("comparison_id")
            .to_dict(orient="index")
        )
        group_data: list[tuple[int, str, pd.DataFrame, np.ndarray]] = []
        for idx, comp_id in enumerate(order):
            group_df = sub[sub["comparison_id"].astype(str) == comp_id].copy()
            vals = pd.to_numeric(group_df["summary_abs_difference"], errors="coerce").dropna().to_numpy(dtype=float)
            if vals.size:
                group_data.append((idx, comp_id, group_df, vals))

        violin_vals = [vals for _, _, _, vals in group_data if vals.size >= 2]
        violin_pos = [float(positions[idx]) for idx, _, _, vals in group_data if vals.size >= 2]
        violin_ids = [comp_id for _, comp_id, _, vals in group_data if vals.size >= 2]
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
            for body, comp_id in zip(vp["bodies"], violin_ids):
                method = str(label_map.get(comp_id, {}).get("primary_method", ""))
                body.set_facecolor(_group_color(method))
                body.set_edgecolor("#333333")
                body.set_alpha(0.52)
            if "cmedians" in vp:
                vp["cmedians"].set_color("#111111")
                vp["cmedians"].set_linewidth(1.1)

        scatter_size = max(1.0, point_size)
        for idx, comp_id, group_df, vals in group_data:
            xj = positions[idx] + rng.normal(0.0, max(0.0, point_jitter), size=vals.size)
            method = str(label_map.get(comp_id, {}).get("primary_method", ""))
            color = _group_color(method)
            if "v4_reference_status" in group_df.columns:
                status = group_df["v4_reference_status"].fillna("single_v4_reference").astype(str).to_numpy()
            else:
                status = np.repeat("single_v4_reference", vals.size)
            for status_key, marker in SAMPLE_STATUS_MARKERS.items():
                mask = status == status_key
                if not np.any(mask):
                    continue
                ax.scatter(
                    xj[mask],
                    vals[mask],
                    s=scatter_size,
                    marker=marker,
                    color=color,
                    alpha=max(0.0, min(1.0, point_opacity)),
                    edgecolors="#111111",
                    linewidths=0.45,
                    zorder=3,
                )

        finite_vals = pd.to_numeric(sub["summary_abs_difference"], errors="coerce").dropna().to_numpy(dtype=float)
        y_min = float(np.nanmin(finite_vals)) if finite_vals.size else 0.0
        y_max = float(np.nanmax(finite_vals)) if finite_vals.size else 1.0
        y_range = max(0.03, y_max - y_min)
        bracket_count = sum(1 for r in annotated_rows if str(r.get("facet_id", "")) == facet_id)
        top_pad = (0.42 + 0.08 * min(8, bracket_count)) * y_range
        ax.set_ylim(max(0.0, y_min - 0.08 * y_range), y_max + top_pad)
        ax.set_xticks(positions)
        ax.set_xticklabels(
            [str(label_map.get(comp_id, {}).get("comparison_label", comp_id)) for comp_id in order],
            rotation=35,
            ha="right",
            va="top",
            rotation_mode="anchor",
            fontsize=fs["tick"],
        )
        for tick in ax.get_xticklabels():
            tick.set_linespacing(0.9)
        ax.set_title(str(facet["facet_label"]), fontsize=fs["base"], loc="left")
        ax.set_ylabel(f"{stat_label.capitalize()} abs. difference", fontsize=fs["base"])
        ax.tick_params(axis="y", labelsize=fs["tick"])
        ax.grid(True, axis="y", alpha=0.25, linewidth=0.6)
        rows_for_facet = [r for r in annotated_rows if str(r.get("facet_id", "")) == facet_id]
        if rows_for_facet:
            _draw_pairwise_brackets(
                ax=ax,
                pairwise_rows=rows_for_facet,
                method_order=order,
                alpha=stats_alpha,
                label_fontsize=fs["pair_label"],
                star_fontsize=fs["pair_star"],
            )

    fig.suptitle(
        "Hybrid-focused metric-family contrast distributions",
        fontsize=fs["base"] + 1.0,
        y=0.992,
    )
    fig.text(
        0.5,
        0.012,
        "Each violin shows one summarized absolute pairwise-difference value per sample.",
        ha="center",
        va="center",
        fontsize=fs["base"],
    )
    fig.subplots_adjust(left=0.065, right=0.995, bottom=0.075, top=0.935, hspace=0.50, wspace=0.34)

    out_path = out_dir / f"{OUTPUT_STEM}.svg"
    fig.savefig(out_path, format="svg", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _hybrid_panel_f_points(
    long: pd.DataFrame,
    method_groups: list[dict],
    metric_order: list[str],
    agnostic_metrics: list[str],
    aware_metrics: list[str],
) -> pd.DataFrame:
    metric_family_by_metric = {
        metric: _family_for_metric(metric, agnostic_metrics=agnostic_metrics, aware_metrics=aware_metrics)
        for metric in metric_order
    }
    metric_family_by_metric = {m: f for m, f in metric_family_by_metric.items() if f is not None}
    method_to_group: dict[str, str] = {}
    source_methods_by_group: dict[str, str] = {}
    for group in method_groups:
        group_id = str(group["group_id"])
        sources = [str(m) for m in group.get("source_methods", [])]
        source_methods_by_group[group_id] = ",".join(sources)
        for method in sources:
            method_to_group[method] = group_id

    sub = long[
        long["method"].astype(str).isin(method_to_group.keys())
        & long["metric"].astype(str).isin(metric_family_by_metric)
    ].copy()
    if sub.empty:
        return pd.DataFrame()
    sub["distance"] = pd.to_numeric(sub["distance"], errors="coerce")
    sub = sub[np.isfinite(sub["distance"])].copy()
    sub["metric_family"] = sub["metric"].astype(str).map(metric_family_by_metric)
    sub["method_group"] = sub["method"].astype(str).map(method_to_group)
    grouped = (
        sub.groupby(["community_id", "sample_id", "sample_key", "metric_family", "method_group"], sort=False)
        .agg(
            mean_distance=("distance", "mean"),
            median_distance=("distance", "median"),
            n_endpoint_values=("distance", "size"),
        )
        .reset_index()
    )
    grouped["community_label"] = grouped["community_id"].astype(str).map(display_community_id)
    grouped["sample_label"] = grouped["sample_key"].astype(str).map(display_sample_key)
    grouped["metric_family_label"] = grouped["metric_family"].map(FAMILY_LABELS).fillna(grouped["metric_family"].astype(str))
    grouped["method_group_label"] = grouped["method_group"].map(_display_group)
    grouped["source_methods"] = grouped["method_group"].map(source_methods_by_group).fillna("")
    grouped["agnostic_metrics"] = ",".join(agnostic_metrics)
    grouped["aware_metrics"] = ",".join(aware_metrics)
    return grouped[
        [
            "community_id",
            "sample_id",
            "sample_key",
            "community_label",
            "sample_label",
            "metric_family",
            "metric_family_label",
            "method_group",
            "method_group_label",
            "source_methods",
            "mean_distance",
            "median_distance",
            "n_endpoint_values",
            "agnostic_metrics",
            "aware_metrics",
        ]
    ].sort_values(["metric_family", "method_group", "community_id", "sample_id", "sample_key"])


def _hybrid_panel_f_stats(
    points: pd.DataFrame,
    out_dir: Path,
    n_perm: int,
    seed: int,
    alpha: float,
    pairwise_test: str,
    multiple_test_correction: bool,
) -> tuple[dict[str, dict], list[dict]]:
    stats_by_family: dict[str, dict] = {}
    omnibus_rows: list[dict] = []
    pairwise_rows: list[dict] = []
    annotated_rows: list[dict] = []
    if points.empty:
        pd.DataFrame().to_csv(out_dir / f"{OUTPUT_STEM_PANEL_F}_stats_omnibus.tsv", sep="\t", index=False)
        pd.DataFrame().to_csv(out_dir / f"{OUTPUT_STEM_PANEL_F}_stats_pairwise.tsv", sep="\t", index=False)
        pd.DataFrame().to_csv(out_dir / f"{OUTPUT_STEM_PANEL_F}_stats_annotated.tsv", sep="\t", index=False)
        return stats_by_family, annotated_rows

    label_map = {"OTU": "UPARSE", "Hybrid": "Hybrid", "zOTU": "UNOISE3"}
    planned_pairs = {("OTU", "Hybrid"), ("Hybrid", "OTU"), ("Hybrid", "zOTU"), ("zOTU", "Hybrid")}
    for idx, family in enumerate(PANEL_F_FAMILY_ORDER, start=1):
        sub = points[points["metric_family"].astype(str) == family].copy()
        methods = [m for m in PANEL_F_METHOD_GROUP_ORDER if m in set(sub["method_group"].astype(str))]
        if sub.empty or len(methods) < 2:
            continue
        stats_df = sub.rename(columns={"method_group": "method", "mean_distance": "distance"}).copy()
        st = _compute_panel_stats_permutation(
            df=stats_df,
            method_order=methods,
            n_perm=n_perm,
            seed=seed + 43000 + idx,
            alpha=alpha,
            pairwise_test=pairwise_test,
            multiple_test_correction=multiple_test_correction,
        )
        stats_by_family[family] = st
        omnibus_rows.append(
            {
                "metric": "hybrid_method_panel_f",
                "metric_family": family,
                "metric_family_label": FAMILY_LABELS.get(family, family),
                "stats_model": st.get("stats_model", "paired_permutation"),
                "pairwise_test": pairwise_test,
                "multiple_test_correction": bool(st.get("multiple_test_correction", multiple_test_correction)),
                "n_blocks": st.get("n_blocks", 0),
                "n_method_groups": len(methods),
                "omnibus_stat": st.get("omnibus_stat", np.nan),
                "omnibus_p": st.get("omnibus_p", np.nan),
                "omnibus_significant_alpha": bool(
                    np.isfinite(float(st.get("omnibus_p", np.nan)))
                    and float(st.get("omnibus_p", np.nan)) < alpha
                ),
                "alpha": alpha,
                "n_permutations": n_perm,
            }
        )
        for row in st.get("pairwise", []):
            rec = dict(row)
            a = str(rec.get("method_a", ""))
            b = str(rec.get("method_b", ""))
            rec.update(
                {
                    "metric": "hybrid_method_panel_f",
                    "metric_family": family,
                    "metric_family_label": FAMILY_LABELS.get(family, family),
                    "method_a_label": label_map.get(a, a),
                    "method_b_label": label_map.get(b, b),
                    "alpha": alpha,
                    "n_permutations": n_perm,
                }
            )
            pairwise_rows.append(rec)
            if (a, b) in planned_pairs:
                annotated_rows.append(dict(rec, annotation_mode="planned"))

    pd.DataFrame(omnibus_rows).to_csv(out_dir / f"{OUTPUT_STEM_PANEL_F}_stats_omnibus.tsv", sep="\t", index=False)
    pd.DataFrame(pairwise_rows).to_csv(out_dir / f"{OUTPUT_STEM_PANEL_F}_stats_pairwise.tsv", sep="\t", index=False)
    pd.DataFrame(annotated_rows).to_csv(out_dir / f"{OUTPUT_STEM_PANEL_F}_stats_annotated.tsv", sep="\t", index=False)
    return stats_by_family, annotated_rows


def _draw_hybrid_panel_f(
    points: pd.DataFrame,
    annotated_rows: list[dict],
    out_dir: Path,
    benchmark_outdir: Path,
    font_size: float,
    violin_bw_method,
    plot_type: str,
    point_size: float,
    point_opacity: float,
    point_jitter: float,
    stats_alpha: float,
) -> Path:
    fs = _fig21_text_scale(font_size)
    out_path = out_dir / f"{OUTPUT_STEM_PANEL_F}.svg"
    if points.empty:
        fig, ax = plt.subplots(figsize=(11.6, 7.0))
        ax.axis("off")
        ax.text(0.5, 0.5, "Hybrid-focused Panel F unavailable", ha="center", va="center", transform=ax.transAxes)
        fig.savefig(out_path, format="svg", dpi=300, bbox_inches="tight")
        plt.close(fig)
        return out_path

    panel_df = points.copy()
    panel_df = _attach_sample_v4_status(panel_df, benchmark_outdir)
    if "v4_reference_status" not in panel_df.columns:
        panel_df["v4_reference_status"] = "single_v4_reference"
    panel_df["v4_reference_status"] = panel_df["v4_reference_status"].fillna("single_v4_reference").astype(str)
    panel_df["v4_reference_status"] = panel_df["v4_reference_status"].where(
        panel_df["v4_reference_status"].isin(SAMPLE_STATUS_MARKERS),
        "single_v4_reference",
    )

    fig, ax = plt.subplots(figsize=(13.6, 7.6))
    family_colors = {
        "phylogeny_agnostic": "#FFF1E6",
        "phylogeny_aware": "#EAF4FB",
    }
    positions: dict[tuple[str, str], float] = {}
    x_ticks: list[float] = []
    x_labels: list[str] = []
    for i, family in enumerate(PANEL_F_FAMILY_ORDER):
        base = 1.0 + (i * 4.25)
        for j, method_group in enumerate(PANEL_F_METHOD_GROUP_ORDER):
            x = base + j
            positions[(family, method_group)] = x
            x_ticks.append(x)
            x_labels.append(f"{FAMILY_LABELS.get(family, family)}\n{_display_group(method_group)}")
        ax.axvspan(base - 0.48, base + 2.48, color=family_colors.get(family, "#f2f2f2"), alpha=0.55, zorder=0)

    line_source = panel_df.pivot_table(
        index=["metric_family", "sample_key"],
        columns="method_group",
        values="mean_distance",
        aggfunc="first",
    ).reset_index()
    for row in line_source.itertuples(index=False):
        family = str(getattr(row, "metric_family"))
        vals = []
        xs = []
        for method_group in PANEL_F_METHOD_GROUP_ORDER:
            if method_group not in line_source.columns:
                continue
            val = getattr(row, method_group, np.nan)
            if np.isfinite(val) and (family, method_group) in positions:
                vals.append(float(val))
                xs.append(positions[(family, method_group)])
        if len(vals) >= 2:
            ax.plot(xs, vals, color="#555555", linewidth=0.85, alpha=0.24, zorder=1)

    grouped_data: list[tuple[str, str, np.ndarray]] = []
    for family in PANEL_F_FAMILY_ORDER:
        for method_group in PANEL_F_METHOD_GROUP_ORDER:
            vals = pd.to_numeric(
                panel_df.loc[
                    (panel_df["metric_family"].astype(str) == family)
                    & (panel_df["method_group"].astype(str) == method_group),
                    "mean_distance",
                ],
                errors="coerce",
            ).dropna().to_numpy(dtype=float)
            if vals.size:
                grouped_data.append((family, method_group, vals))

    if plot_type == "boxplot":
        for family, method_group, vals in grouped_data:
            ax.boxplot(
                [vals],
                positions=[positions[(family, method_group)]],
                widths=0.48,
                patch_artist=True,
                boxprops=dict(facecolor=_group_color(method_group), alpha=0.36, edgecolor="#333333"),
                medianprops=dict(color="#111111", linewidth=1.25),
                whiskerprops=dict(color="#333333", linewidth=0.9),
                capprops=dict(color="#333333", linewidth=0.9),
            )
    else:
        violin_vals = [vals for _, _, vals in grouped_data if vals.size >= 2]
        violin_pos = [float(positions[(family, method_group)]) for family, method_group, vals in grouped_data if vals.size >= 2]
        violin_groups = [(family, method_group) for family, method_group, vals in grouped_data if vals.size >= 2]
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
            for body, (_, method_group) in zip(vp["bodies"], violin_groups):
                body.set_facecolor(_group_color(method_group))
                body.set_edgecolor("#333333")
                body.set_alpha(0.40)
            if "cmedians" in vp:
                vp["cmedians"].set_color("#111111")
                vp["cmedians"].set_linewidth(1.15)

    rng = np.random.default_rng(42141)
    for family in PANEL_F_FAMILY_ORDER:
        for method_group in PANEL_F_METHOD_GROUP_ORDER:
            group_df = panel_df.loc[
                (panel_df["metric_family"].astype(str) == family)
                & (panel_df["method_group"].astype(str) == method_group)
            ].copy()
            if group_df.empty:
                continue
            vals = pd.to_numeric(group_df["mean_distance"], errors="coerce").to_numpy(dtype=float)
            xj = positions[(family, method_group)] + rng.normal(0.0, max(0.0, float(point_jitter)), size=vals.size)
            status = group_df["v4_reference_status"].astype(str).to_numpy()
            for status_key, marker in SAMPLE_STATUS_MARKERS.items():
                mask = status == status_key
                if not np.any(mask):
                    continue
                ax.scatter(
                    xj[mask],
                    vals[mask],
                    s=point_size,
                    marker=marker,
                    alpha=max(0.0, min(1.0, float(point_opacity))),
                    color=_group_color(method_group),
                    edgecolors="#111111",
                    linewidths=0.45,
                    zorder=3,
                )

    finite_vals = pd.to_numeric(panel_df["mean_distance"], errors="coerce").dropna().to_numpy(dtype=float)
    if finite_vals.size:
        y_min = float(np.min(finite_vals))
        y_max = float(np.max(finite_vals))
        y_range = max(0.05, y_max - y_min)
        ax.set_ylim(max(0.0, y_min - 0.10 * y_range), y_max + 0.40 * y_range)

    y_min, y_max = ax.get_ylim()
    y_range = max(1e-9, y_max - y_min)
    ann_by_family: dict[str, list[dict]] = {}
    for row in annotated_rows:
        ann_by_family.setdefault(str(row.get("metric_family", "")), []).append(row)
    for family, rows in ann_by_family.items():
        rows = [
            r
            for r in rows
            if str(r.get("method_a", "")) in PANEL_F_METHOD_GROUP_ORDER
            and str(r.get("method_b", "")) in PANEL_F_METHOD_GROUP_ORDER
        ]
        for level, row in enumerate(rows):
            a = str(row.get("method_a", ""))
            b = str(row.get("method_b", ""))
            if (family, a) not in positions or (family, b) not in positions:
                continue
            x1 = positions[(family, a)]
            x2 = positions[(family, b)]
            if x1 > x2:
                x1, x2 = x2, x1
            p = float(row.get("p_value_used", row.get("p_adj_bh", row.get("p_raw", np.nan))))
            if not np.isfinite(p):
                continue
            y = y_max - (0.19 - 0.08 * level) * y_range
            cap = 0.017 * y_range
            ax.plot([x1, x1, x2, x2], [y, y + cap, y + cap, y], color="#111111", linewidth=0.95, zorder=4)
            ax.text(
                (x1 + x2) / 2.0,
                y + cap + 0.010 * y_range,
                f"p={p:.3g}\n{_pvalue_stars(p, alpha=stats_alpha)}",
                ha="center",
                va="bottom",
                fontsize=fs["stats"],
                color="#111111",
                zorder=5,
            )

    ax.set_xlim(min(x_ticks) - 0.58, max(x_ticks) + 0.58)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_labels, fontsize=fs["tick"])
    ax.set_xlabel("Metric family | analysis method category", fontsize=fs["base"])
    ax.set_ylabel("Mean distance to expected composition", fontsize=fs["base"])
    ax.tick_params(axis="y", labelsize=fs["tick"])
    ax.grid(True, axis="y", alpha=0.25, linewidth=0.6)
    ax.set_title("Hybrid-focused distance to expected composition", fontsize=fs["base"])
    handles = [
        Patch(facecolor=_group_color(method_group), edgecolor="#333333", alpha=0.40, label=_display_group(method_group))
        for method_group in PANEL_F_METHOD_GROUP_ORDER
    ]
    handles.append(Line2D([0], [0], color="#555555", linewidth=0.85, alpha=0.45, label="Paired sample trajectory"))
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize=fs["legend"])
    fig.tight_layout()
    fig.savefig(out_path, format="svg", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path


def run(config_path: Path, source_tsv: str | None = None, out_dir_arg: str | None = None) -> Path:
    cfg = load_config(config_path)
    exclude_patterns = list(getattr(cfg.figures, "exclude_community_name_patterns", []) or [])
    exclude_patterns += list(getattr(cfg.figures, "exclude_community_name_patterns_fig21", []) or [])
    configure_community_display_mapping_from_manifest(Path(cfg.paths.outdir), figure_exclude_patterns=exclude_patterns)
    out_dir = _resolve_outdir(cfg, out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)

    source_path = _resolve_source_tsv(cfg, source_tsv)
    if not source_path.exists():
        raise FileNotFoundError(
            f"Could not find Figure 21 distance cache: {source_path}. "
            "Run `PYTHONPATH=src python -m tacticbench.cli plot-fig21 --config ...` first, "
            "or pass --source-tsv."
        )

    long = pd.read_csv(source_path, sep="\t")
    required = {"community_id", "sample_id", "method", "metric", "distance"}
    missing = required.difference(set(long.columns))
    if missing:
        raise ValueError(f"Input distance TSV is missing required columns: {sorted(missing)}")
    long["community_id"] = long["community_id"].astype(str)
    long["sample_id"] = long["sample_id"].astype(str)
    long["method"] = long["method"].astype(str)
    long["metric"] = long["metric"].astype(str)
    long["distance"] = pd.to_numeric(long["distance"], errors="coerce")
    long = long[np.isfinite(long["distance"])].copy()
    if "sample_key" not in long.columns:
        long["sample_key"] = long["community_id"] + "|" + long["sample_id"]
    else:
        long["sample_key"] = long["sample_key"].astype(str)

    method_order = _configured_method_order(cfg, long)
    method_groups = _focused_method_groups(method_order)
    method_group_order = [str(g["group_id"]) for g in method_groups]
    metric_order = _configured_metric_order(cfg, long)
    agnostic_metrics, aware_metrics = _metric_family_lists(cfg, metric_order)
    summary_stat = _normalize_summary_stat(getattr(cfg.figures, "fig21_pair_group_summary_stat", "median"))
    font_size = _resolve_fig21_panel_font_size(
        getattr(cfg.figures, "fig21_panel_b_font_size", None),
        max(6.0, float(getattr(cfg.figures, "fig21_font_size", 12.0))),
    )
    violin_bw_method = _normalize_violin_bw_method(getattr(cfg.figures, "fig21_violin_bw_method", None))
    point_size = getattr(cfg.figures, "fig21_panel_b_dot_size", None)
    point_size = max(1.0, float(point_size)) if point_size is not None else max(1.0, float(getattr(cfg.figures, "fig21_point_size", 42.0)))
    point_opacity = max(0.0, min(1.0, float(getattr(cfg.figures, "fig21_panel_b_dot_opacity", 0.48))))
    point_jitter = max(0.0, float(getattr(cfg.figures, "fig21_panel_b_dot_jitter", 0.065)))
    panel_f_font_size = _resolve_fig21_panel_font_size(
        getattr(cfg.figures, "fig21_panel_f_font_size", None),
        max(6.0, float(getattr(cfg.figures, "fig21_font_size", 12.0))),
    )
    panel_f_dot_size = getattr(cfg.figures, "fig21_panel_f_dot_size", None)
    panel_f_dot_size = max(1.0, float(panel_f_dot_size)) if panel_f_dot_size is not None else point_size
    panel_f_dot_opacity = max(0.0, min(1.0, float(getattr(cfg.figures, "fig21_panel_f_dot_opacity", 0.75))))
    panel_f_dot_jitter = max(0.0, float(getattr(cfg.figures, "fig21_panel_f_dot_jitter", 0.06)))
    plot_type = str(getattr(cfg.figures, "fig21_plot_type", "violin") or "violin").strip().lower()
    if plot_type not in {"violin", "boxplot"}:
        plot_type = "violin"
    stats_enabled = bool(getattr(cfg.figures, "fig21_stats_enabled", True))
    stats_permutations = max(100, int(getattr(cfg.figures, "fig21_stats_permutations", 4000)))
    stats_seed = int(getattr(cfg.figures, "fig21_stats_seed", 1))
    stats_alpha = float(getattr(cfg.figures, "fig21_stats_alpha", 0.05))
    stats_pairwise_test = _normalize_pairwise_test(getattr(cfg.figures, "fig21_stats_pairwise_test", "sign_test_exact"))
    stats_multiple_test_correction = bool(getattr(cfg.figures, "fig21_stats_multiple_test_correction", True))
    stats_show_pairwise = bool(getattr(cfg.figures, "fig21_stats_show_pairwise", True))
    show_within_contrast_stats = bool(getattr(cfg.figures, "fig21_panel_b_show_within_contrast_stats", False))
    pair_groups_to_show = list(getattr(cfg.figures, "fig21_panel_b_pair_groups_to_show", []) or [])
    if not pair_groups_to_show:
        pair_groups_to_show = [
            "denoising_phylogeny_aware__vs__denoising_phylogeny_aware",
            "denoising_phylogeny_agnostic__vs__denoising_phylogeny_agnostic",
            "clustering_phylogeny_aware__vs__clustering_phylogeny_aware",
            "clustering_phylogeny_agnostic__vs__clustering_phylogeny_agnostic",
            "denoising_phylogeny_aware__vs__denoising_phylogeny_agnostic",
            "clustering_phylogeny_aware__vs__clustering_phylogeny_agnostic",
            "denoising_phylogeny_aware__vs__clustering_phylogeny_aware",
            "denoising_phylogeny_agnostic__vs__clustering_phylogeny_agnostic",
        ]

    endpoints, pairwise, sample_summary, endpoint_group_order = _build_method_family_tables(
        long=long,
        method_groups=method_groups,
        metric_order=metric_order,
        agnostic_metrics=agnostic_metrics,
        aware_metrics=aware_metrics,
        summary_stat=summary_stat,
    )
    endpoints.to_csv(out_dir / f"{OUTPUT_STEM}_endpoint_values.tsv", sep="\t", index=False)
    pairwise.to_csv(out_dir / f"{OUTPUT_STEM}_pairwise_differences.tsv", sep="\t", index=False)
    sample_summary.to_csv(out_dir / f"{OUTPUT_STEM}_sample_summary.tsv", sep="\t", index=False)
    _write_cell_summary(sample_summary, out_dir)
    sample_summary = _attach_sample_v4_status(sample_summary, Path(cfg.paths.outdir))

    group_map = []
    source_methods_by_group = {
        str(group["group_id"]): ",".join(str(m) for m in group.get("source_methods", []))
        for group in method_groups
    }
    for group in endpoint_group_order:
        method_group, family = group.split("__", 1)
        group_map.append(
            {
                "endpoint_group_id": group,
                "method_group": method_group,
                "method_group_label": _display_group(method_group),
                "source_methods": source_methods_by_group.get(method_group, ""),
                "metric_family": family,
                "metric_family_label": FAMILY_LABELS.get(family, family),
                "endpoint_group_label": _group_label(method_group, family),
            }
        )
    pd.DataFrame(group_map).to_csv(out_dir / f"{OUTPUT_STEM}_endpoint_group_map.tsv", sep="\t", index=False)

    plot_df, spec_df, facets = _select_method_family_panel_b_data(
        sample_summary=sample_summary,
        method_order=method_group_order,
    )
    plot_df.to_csv(out_dir / f"{OUTPUT_STEM}_faceted_plot_data.tsv", sep="\t", index=False)
    spec_df.to_csv(out_dir / f"{OUTPUT_STEM}_faceted_comparison_map.tsv", sep="\t", index=False)
    if stats_enabled:
        _, annotated_rows = _stats_for_method_family_panel_b(
            plot_df=plot_df,
            facets=facets,
            out_dir=out_dir,
            n_perm=stats_permutations,
            seed=stats_seed,
            alpha=stats_alpha,
            pairwise_test=stats_pairwise_test,
            multiple_test_correction=stats_multiple_test_correction,
        )
    else:
        annotated_rows = []
        pd.DataFrame().to_csv(out_dir / f"{OUTPUT_STEM}_faceted_stats_omnibus.tsv", sep="\t", index=False)
        pd.DataFrame().to_csv(out_dir / f"{OUTPUT_STEM}_faceted_stats_pairwise.tsv", sep="\t", index=False)
        pd.DataFrame().to_csv(out_dir / f"{OUTPUT_STEM}_faceted_stats_annotated.tsv", sep="\t", index=False)

    filtered_endpoints, filtered_pairwise, filtered_sample_summary = _filtered_family_panel_b_tables(
        long=long,
        metric_order=metric_order,
        agnostic_metrics=agnostic_metrics,
        aware_metrics=aware_metrics,
        summary_stat=summary_stat,
    )
    filtered_endpoints.to_csv(out_dir / f"{OUTPUT_STEM_PANEL_B_FILTERED}_endpoint_values.tsv", sep="\t", index=False)
    filtered_pairwise.to_csv(out_dir / f"{OUTPUT_STEM_PANEL_B_FILTERED}_pairwise_differences.tsv", sep="\t", index=False)
    filtered_sample_summary.to_csv(out_dir / f"{OUTPUT_STEM_PANEL_B_FILTERED}_sample_summary.tsv", sep="\t", index=False)
    filtered_sample_summary = _attach_sample_v4_status(filtered_sample_summary, Path(cfg.paths.outdir))
    filtered_stats, filtered_planned_rows = _filtered_family_panel_b_stats(
        sample_summary=filtered_sample_summary,
        out_dir=out_dir,
        stats_enabled=stats_enabled,
        n_perm=stats_permutations,
        seed=stats_seed,
        alpha=stats_alpha,
        pairwise_test=stats_pairwise_test,
        multiple_test_correction=stats_multiple_test_correction,
    )
    filtered_panel_b_path = _draw_filtered_family_panel_b(
        sample_summary=filtered_sample_summary,
        stats=filtered_stats,
        planned_rows=filtered_planned_rows,
        out_dir=out_dir,
        font_size=font_size,
        violin_bw_method=violin_bw_method,
        plot_type=plot_type,
        point_size=point_size,
        point_opacity=point_opacity,
        point_jitter=point_jitter,
        stats_enabled=stats_enabled,
        stats_alpha=stats_alpha,
        stats_show_pairwise=stats_show_pairwise,
        show_within_contrast_stats=show_within_contrast_stats,
        pair_groups_to_show=pair_groups_to_show,
    )

    panel_f_points = _hybrid_panel_f_points(
        long=long,
        method_groups=method_groups,
        metric_order=metric_order,
        agnostic_metrics=agnostic_metrics,
        aware_metrics=aware_metrics,
    )
    panel_f_points.to_csv(out_dir / f"{OUTPUT_STEM_PANEL_F}_points.tsv", sep="\t", index=False)
    if stats_enabled:
        _, panel_f_annotated_rows = _hybrid_panel_f_stats(
            points=panel_f_points,
            out_dir=out_dir,
            n_perm=stats_permutations,
            seed=stats_seed,
            alpha=stats_alpha,
            pairwise_test=stats_pairwise_test,
            multiple_test_correction=stats_multiple_test_correction,
        )
    else:
        panel_f_annotated_rows = []
        pd.DataFrame().to_csv(out_dir / f"{OUTPUT_STEM_PANEL_F}_stats_omnibus.tsv", sep="\t", index=False)
        pd.DataFrame().to_csv(out_dir / f"{OUTPUT_STEM_PANEL_F}_stats_pairwise.tsv", sep="\t", index=False)
        pd.DataFrame().to_csv(out_dir / f"{OUTPUT_STEM_PANEL_F}_stats_annotated.tsv", sep="\t", index=False)
    panel_f_path = _draw_hybrid_panel_f(
        points=panel_f_points,
        annotated_rows=panel_f_annotated_rows,
        out_dir=out_dir,
        benchmark_outdir=Path(cfg.paths.outdir),
        font_size=panel_f_font_size,
        violin_bw_method=violin_bw_method,
        plot_type=plot_type,
        point_size=panel_f_dot_size,
        point_opacity=panel_f_dot_opacity,
        point_jitter=panel_f_dot_jitter,
        stats_alpha=stats_alpha,
    )

    panel_b_path = _draw_panel_b_style_rows(
        plot_df=plot_df,
        facets=facets,
        annotated_rows=annotated_rows,
        out_dir=out_dir,
        font_size=font_size,
        violin_bw_method=violin_bw_method,
        summary_stat=summary_stat,
        point_size=point_size,
        point_opacity=point_opacity,
        point_jitter=point_jitter,
        stats_alpha=stats_alpha,
    )
    print(f"Wrote {filtered_panel_b_path}")
    print(f"Wrote {panel_f_path}")
    return panel_b_path


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a Figure 21 panel-B-style auxiliary plot focused on UNOISE3, UPARSE, "
            "and pooled Hybrid endpoints where Hybrid combines TIC and TAC."
        )
    )
    parser.add_argument("--config", required=True, type=Path, help="TACTICBench YAML config.")
    parser.add_argument(
        "--source-tsv",
        default=None,
        help="Optional distance_points_filtered.tsv to read. Defaults to the main Figure 21 output cache.",
    )
    parser.add_argument(
        "--out-dir",
        default=None,
        help=f"Optional output directory. Defaults to figures/{OUTPUT_DIR_NAME}.",
    )
    args = parser.parse_args()

    out_path = run(config_path=args.config, source_tsv=args.source_tsv, out_dir_arg=args.out_dir)
    print(f"Wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
