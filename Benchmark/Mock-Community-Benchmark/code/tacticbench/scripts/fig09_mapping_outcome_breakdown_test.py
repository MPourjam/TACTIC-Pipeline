from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import PathPatch, Polygon, Rectangle
from matplotlib.path import Path as MplPath
from matplotlib.ticker import FuncFormatter

from tacticbench.config import load_config
from tacticbench.steps.plot_fig09 import METHOD_COLORS, _fig09_no_hit_effective_identity_thresholds
from tacticbench.utils.blast import filter_hits, read_hits
from tacticbench.utils.labels import display_method
from tacticbench.utils.paths import expected_dir, figures_root, method_dir


METHOD_ORDER = ["ASV", "zOTU", "OTU", "TIC", "TAC"]
RICHNESS_METHOD_ORDER = ["ASV", "zOTU", "TAC", "TIC", "OTU"]
OUTPUT_STEM = "fig09_mapping_outcome_breakdown_test_exact100"
TRANSITION_OUTPUT_STEM = "fig09_mapping_transition_absorption_test_exact100"
WAVE_OUTPUT_STEM = "fig09_mapping_transition_wave_test_exact100"
HEATMAP_OUTPUT_STEM = "fig09_mapping_transition_heatmap_test_exact100"
PURITY_OUTPUT_STEM = "fig09_cluster_species_purity_test_exact100"
PURITY_VENN_OUTPUT_STEM = "fig09_cluster_species_purity_venn_test_exact100"
MATCH_SPECIES_OUTPUT_STEM = "fig09_match_status_species_purity_test_exact100"
EXACT_NEAR_SPECIES_OUTPUT_STEM = "fig09_exact_near_species_purity_test_exact100"
RICHNESS_OUTPUT_STEM = "fig09_relative_richness_error_test"
CLASS_ORDER = ["Exact match", "Near-exact match", "No hit"]
CLASS_WITH_UNASSIGNED_ORDER = ["Exact match", "Near-exact match", "No hit", "Unassigned"]


def _wave_display_label(label: str) -> str:
    mapping = {
        "Exact match": "Exact match",
        "Near-exact match": "Near-exact match",
        "No hit": "No match",
        "Unassigned": "Non-bacterial",
    }
    return mapping.get(str(label), str(label))


def _format_pct(value: float) -> str:
    pct = float(value) * 100.0
    if not np.isfinite(pct):
        return "NA"
    if pct == 0.0:
        return "0%"
    if abs(pct) < 0.01:
        return "<0.01%" if pct > 0 else ">-0.01%"
    if abs(pct) < 1.0:
        return f"{pct:.2f}%"
    return f"{pct:.1f}%"


def _format_amount(value: float) -> str:
    val = float(value)
    if not np.isfinite(val):
        return "NA"
    if val == 0.0:
        return "0"
    if abs(val) < 0.01:
        return "<0.01" if val > 0 else ">-0.01"
    if abs(val) < 1.0:
        return f"{val:.2f}"
    return f"{val:.1f}"


def _resolve_no_hit_qcov(cfg) -> float:
    fig09_min_qcov = float(getattr(cfg.figures, "fig09_min_query_coverage", 0.0))
    no_hit_qcov = getattr(cfg.figures, "fig09_no_hit_min_query_coverage", None)
    if no_hit_qcov is None:
        no_hit_qcov = getattr(cfg.figures, "fig12_unmapped_best_hit_min_query_coverage", None)
    if no_hit_qcov is None:
        return fig09_min_qcov
    return float(no_hit_qcov)


def _best_hits_by_feature(
    blast_hits_path: Path,
    min_qcov: float,
) -> pd.DataFrame:
    if (not blast_hits_path.exists()) or blast_hits_path.stat().st_size == 0:
        return pd.DataFrame()

    hits = read_hits(blast_hits_path)
    hits = filter_hits(
        hits,
        min_effective_pident=None,
        min_qcovhsp=min_qcov,
        strict_pident_gt=False,
    )
    if hits.empty:
        return pd.DataFrame()

    hits = hits.copy()
    hits["qseqid"] = hits["qseqid"].astype(str)
    hits["sseqid"] = hits["sseqid"].astype(str)
    hits["effective_identity"] = pd.to_numeric(hits["effective_identity"], errors="coerce")
    hits["pident"] = pd.to_numeric(hits["pident"], errors="coerce")
    hits["bitscore"] = pd.to_numeric(hits["bitscore"], errors="coerce")
    hits["evalue"] = pd.to_numeric(hits["evalue"], errors="coerce")
    hits["length"] = pd.to_numeric(hits["length"], errors="coerce")
    hits = hits.dropna(subset=["qseqid", "effective_identity", "pident"])
    if hits.empty:
        return pd.DataFrame()

    return (
        hits.sort_values(
            ["qseqid", "effective_identity", "pident", "bitscore", "evalue", "length"],
            ascending=[True, False, False, False, True, False],
        )
        .groupby("qseqid", sort=False)
        .first()
    )


def _best_effective_identity_by_feature(
    blast_hits_path: Path,
    min_qcov: float,
) -> dict[str, float]:
    best_hits = _best_hits_by_feature(blast_hits_path, min_qcov=min_qcov)
    if best_hits.empty:
        return {}
    return {
        str(feature_id): float(eff)
        for feature_id, eff in zip(
            best_hits.index.astype(str),
            best_hits["effective_identity"].astype(float),
        )
    }


def _load_tip_to_species(outdir: Path, community_id: str) -> dict[str, str]:
    tip_map_path = expected_dir(outdir, community_id) / "tip_to_species.tsv"
    if not tip_map_path.exists():
        return {}
    tip_map = pd.read_csv(tip_map_path, sep="\t")
    if "tip" not in tip_map.columns or "species" not in tip_map.columns:
        return {}
    return {
        str(tip): str(species)
        for tip, species in zip(
            tip_map["tip"].astype(str),
            tip_map["species"].astype(str),
        )
    }


def _classify_feature(best_eff: float, threshold: float) -> str:
    if np.isfinite(best_eff) and best_eff >= (100.0 - 1e-9):
        return "Exact match"
    if np.isfinite(best_eff) and (best_eff > float(threshold)) and (best_eff < (100.0 - 1e-9)):
        return "Near-exact match"
    return "No hit"


def _load_manifest(outdir: Path) -> pd.DataFrame:
    manifest_path = outdir / "manifest.tsv"
    if not manifest_path.exists():
        raise FileNotFoundError(f"TACTICBench manifest not found: {manifest_path}")
    manifest = pd.read_csv(manifest_path, sep="\t")
    if "community_id" not in manifest.columns:
        raise ValueError(f"TACTICBench manifest is missing required column 'community_id': {manifest_path}")
    manifest["community_id"] = manifest["community_id"].astype(str)
    return manifest


def _membership_path_from_manifest_row(manifest_row: pd.Series, method: str) -> Path:
    if method == "TIC":
        table_key = "TIC_table"
        filename = "Map-SOTU-ZOTU.tab"
    elif method == "TAC":
        table_key = "TAC_table"
        filename = "Map-ZOTUs-OTUs-TAC.tab"
    else:
        raise ValueError(f"Cluster membership is only defined here for TIC/TAC, not {method}")

    table_path = Path(str(manifest_row[table_key]))
    return table_path.parent / filename


def _load_cluster_membership_map(manifest_row: pd.Series, method: str) -> dict[str, set[str]]:
    membership_path = _membership_path_from_manifest_row(manifest_row, method)
    if not membership_path.exists():
        raise FileNotFoundError(f"Missing {method} cluster-membership map: {membership_path}")

    raw = pd.read_csv(membership_path, sep="\t", header=None, dtype=str)
    if raw.empty or raw.shape[1] < 2:
        return {}

    col0 = raw.iloc[:, 0].astype(str).str.strip()
    col1 = raw.iloc[:, 1].astype(str).str.strip()
    if method == "TIC":
        if raw.shape[0] > 0 and col0.iloc[0].lower() == "sotu" and col1.iloc[0].lower() == "zotu":
            raw = raw.iloc[1:, :2].copy()
        else:
            raw = raw.iloc[:, :2].copy()
        raw.columns = ["cluster_id", "member_id"]
    else:
        if raw.shape[0] > 0 and col0.iloc[0].lower() == "zotu" and col1.iloc[0].lower() == "otu":
            raw = raw.iloc[1:, :2].copy()
        else:
            raw = raw.iloc[:, :2].copy()
        raw.columns = ["member_id", "cluster_id"]

    raw["cluster_id"] = raw["cluster_id"].astype(str).str.strip()
    raw["member_id"] = raw["member_id"].astype(str).str.strip()
    raw = raw[(raw["cluster_id"] != "") & (raw["member_id"] != "")]

    cluster_map: dict[str, set[str]] = {}
    for rec in raw.itertuples(index=False):
        cluster_map.setdefault(str(rec.cluster_id), set()).add(str(rec.member_id))
    return cluster_map


def _pick_priority_class(classes: list[str]) -> str:
    if "Exact match" in classes:
        return "Exact match"
    if "Near-exact match" in classes:
        return "Near-exact match"
    return "No hit"


def _invert_cluster_membership_map(cluster_map: dict[str, set[str]]) -> dict[str, str]:
    member_to_cluster: dict[str, str] = {}
    for cluster_id, members in cluster_map.items():
        for member_id in members:
            member_to_cluster[str(member_id)] = str(cluster_id)
    return member_to_cluster


def _load_registry(fig09_dir: Path) -> pd.DataFrame:
    diagnostics_path = fig09_dir / "split_merge_diagnostics.tsv"
    if not diagnostics_path.exists():
        raise FileNotFoundError(
            f"Figure 09 diagnostics table not found: {diagnostics_path}. "
            "Run Figure 09 first so the helper can reuse the same included samples."
        )
    diag = pd.read_csv(diagnostics_path, sep="\t")
    keep_cols = ["community_id", "sample_id", "sample_key", "method", "expected_present_features"]
    missing = [c for c in keep_cols if c not in diag.columns]
    if missing:
        raise ValueError(
            f"Figure 09 diagnostics table is missing required columns: {', '.join(missing)}"
        )
    registry = diag[keep_cols].drop_duplicates(subset=["community_id", "sample_id", "sample_key", "method"]).copy()
    registry["method"] = registry["method"].astype(str)
    registry["expected_present_features"] = pd.to_numeric(
        registry["expected_present_features"], errors="coerce"
    ).fillna(0.0)
    registry = registry[registry["method"].isin(METHOD_ORDER)].copy()
    return registry


def _build_feature_table(cfg, registry: pd.DataFrame) -> pd.DataFrame:
    outdir = Path(cfg.paths.outdir)
    manifest = _load_manifest(outdir)
    no_hit_qcov = _resolve_no_hit_qcov(cfg)
    no_hit_thresholds = _fig09_no_hit_effective_identity_thresholds(cfg)
    default_threshold = 97.0

    feature_rows: list[dict[str, object]] = []
    best_hit_cache: dict[tuple[str, str], dict[str, float]] = {}
    best_hit_row_cache: dict[tuple[str, str], pd.DataFrame] = {}
    std_cache: dict[tuple[str, str], pd.DataFrame] = {}
    cluster_membership_cache: dict[tuple[str, str], dict[str, set[str]]] = {}
    manifest_row_cache: dict[str, pd.Series] = {}
    tip_to_species_cache: dict[str, dict[str, str]] = {}

    for rec in registry.itertuples(index=False):
        community_id = str(rec.community_id)
        sample_id = str(rec.sample_id)
        sample_key = str(rec.sample_key)
        method = str(rec.method)
        threshold = float(no_hit_thresholds.get(method, default_threshold))

        cache_key = (community_id, method)
        if cache_key not in best_hit_cache:
            best_hit_cache[cache_key] = _best_effective_identity_by_feature(
                method_dir(outdir, community_id, method) / "blast_hits.tsv",
                min_qcov=no_hit_qcov,
            )
        if cache_key not in best_hit_row_cache:
            best_hit_row_cache[cache_key] = _best_hits_by_feature(
                method_dir(outdir, community_id, method) / "blast_hits.tsv",
                min_qcov=no_hit_qcov,
            )
        if cache_key not in std_cache:
            std_path = method_dir(outdir, community_id, method) / "table_std_tss1000.tsv"
            if not std_path.exists():
                raise FileNotFoundError(f"Missing inferred-feature table: {std_path}")
            std_df = pd.read_csv(std_path, sep="\t", index_col=0)
            std_df.index = std_df.index.astype(str)
            std_cache[cache_key] = std_df

        std_df = std_cache[cache_key]
        if sample_id not in std_df.columns:
            raise KeyError(
                f"Sample {sample_id!r} not found in {method} inferred-feature table "
                f"for community {community_id}"
            )

        present = pd.to_numeric(std_df[sample_id], errors="coerce").fillna(0.0)
        present = present[present > 0.0]
        best_by_feature = best_hit_cache[cache_key]
        best_rows = best_hit_row_cache[cache_key]
        if community_id not in tip_to_species_cache:
            tip_to_species_cache[community_id] = _load_tip_to_species(outdir, community_id)
        tip_to_species = tip_to_species_cache[community_id]

        zotu_cache_key = (community_id, "zOTU")
        if method in {"TIC", "TAC"}:
            if zotu_cache_key not in best_hit_cache:
                best_hit_cache[zotu_cache_key] = _best_effective_identity_by_feature(
                    method_dir(outdir, community_id, "zOTU") / "blast_hits.tsv",
                    min_qcov=no_hit_qcov,
                )
            if zotu_cache_key not in best_hit_row_cache:
                best_hit_row_cache[zotu_cache_key] = _best_hits_by_feature(
                    method_dir(outdir, community_id, "zOTU") / "blast_hits.tsv",
                    min_qcov=no_hit_qcov,
                )
            if zotu_cache_key not in std_cache:
                zotu_std_path = method_dir(outdir, community_id, "zOTU") / "table_std_tss1000.tsv"
                if not zotu_std_path.exists():
                    raise FileNotFoundError(f"Missing zOTU inferred-feature table: {zotu_std_path}")
                zotu_std_df = pd.read_csv(zotu_std_path, sep="\t", index_col=0)
                zotu_std_df.index = zotu_std_df.index.astype(str)
                std_cache[zotu_cache_key] = zotu_std_df
            if community_id not in manifest_row_cache:
                match = manifest.loc[manifest["community_id"] == community_id]
                if match.empty:
                    raise KeyError(f"Community {community_id!r} not found in manifest.tsv")
                manifest_row_cache[community_id] = match.iloc[0]
            if cache_key not in cluster_membership_cache:
                cluster_membership_cache[cache_key] = _load_cluster_membership_map(
                    manifest_row_cache[community_id],
                    method,
                )

            zotu_std_df = std_cache[zotu_cache_key]
            if sample_id not in zotu_std_df.columns:
                raise KeyError(
                    f"Sample {sample_id!r} not found in zOTU inferred-feature table "
                    f"for community {community_id}"
                )
            zotu_present = pd.to_numeric(zotu_std_df[sample_id], errors="coerce").fillna(0.0)
            zotu_present = zotu_present[zotu_present > 0.0]
            zotu_best_by_feature = best_hit_cache[zotu_cache_key]
            zotu_best_rows = best_hit_row_cache[zotu_cache_key]
            zotu_threshold = float(no_hit_thresholds.get("zOTU", default_threshold))
            zotu_classification_by_feature = {
                str(feature_id): _classify_feature(float(zotu_best_by_feature.get(str(feature_id), np.nan)), zotu_threshold)
                for feature_id in zotu_present.index.astype(str)
            }
            cluster_membership = cluster_membership_cache[cache_key]
            member_to_cluster = _invert_cluster_membership_map(cluster_membership)
            present_clusters = set(present.index.astype(str))

        if method in {"TIC", "TAC"}:
            for member_id, member_abundance in zotu_present.items():
                member_id = str(member_id)
                cluster_id = member_to_cluster.get(member_id, "")
                best_eff = float(zotu_best_by_feature.get(member_id, np.nan))
                classification = zotu_classification_by_feature.get(member_id, "No hit")
                best_hit_target_id = ""
                best_hit_species = ""
                if member_id in zotu_best_rows.index:
                    best_hit_target_id = str(zotu_best_rows.loc[member_id, "sseqid"])
                    best_hit_species = str(tip_to_species.get(best_hit_target_id, ""))
                if (not cluster_id) or (cluster_id not in present_clusters):
                    feature_rows.append(
                        {
                            "community_id": community_id,
                            "sample_id": sample_id,
                            "sample_key": sample_key,
                            "method": method,
                            "feature_id": member_id,
                            "cluster_id": "",
                            "feature_value_tss1000": float(member_abundance),
                            "near_exact_threshold": zotu_threshold,
                            "rule_exact": "best_hit_effective_identity == 100",
                            "rule_near_exact": f"{zotu_threshold:g} < best_hit_effective_identity < 100",
                            "rule_no_hit": f"best_hit_effective_identity <= {zotu_threshold:g} or no recorded hit",
                            "best_hit_effective_identity": best_eff,
                            "best_hit_target_id": best_hit_target_id,
                            "best_hit_species": best_hit_species,
                            "classification": "Unassigned",
                            "classification_source": f"zOTU member not retained in present {method} cluster set",
                            "cluster_member_count_total": np.nan,
                            "cluster_member_count_present_in_sample": np.nan,
                            "cluster_members_all": "",
                            "cluster_members_present_in_sample": "",
                            "cluster_members_used_for_classification": "",
                            "cluster_member_classifications_used": "",
                        }
                    )
                    continue
                members_all = sorted(cluster_membership.get(cluster_id, set()))
                members_present = [mid for mid in members_all if mid in zotu_classification_by_feature]
                feature_rows.append(
                    {
                        "community_id": community_id,
                        "sample_id": sample_id,
                        "sample_key": sample_key,
                        "method": method,
                        "feature_id": member_id,
                        "cluster_id": cluster_id,
                        "feature_value_tss1000": float(member_abundance),
                        "near_exact_threshold": zotu_threshold,
                        "rule_exact": "best_hit_effective_identity == 100",
                        "rule_near_exact": f"{zotu_threshold:g} < best_hit_effective_identity < 100",
                        "rule_no_hit": f"best_hit_effective_identity <= {zotu_threshold:g} or no recorded hit",
                        "best_hit_effective_identity": best_eff,
                        "best_hit_target_id": best_hit_target_id,
                        "best_hit_species": best_hit_species,
                        "classification": classification,
                        "classification_source": f"zOTU member carried through present {method} cluster",
                        "cluster_member_count_total": len(members_all),
                        "cluster_member_count_present_in_sample": len(members_present),
                        "cluster_members_all": ";".join(members_all),
                        "cluster_members_present_in_sample": ";".join(members_present),
                        "cluster_members_used_for_classification": member_id,
                        "cluster_member_classifications_used": classification,
                    }
                )

        else:
            for feature_id, abundance in present.items():
                best_eff = float(best_by_feature.get(str(feature_id), np.nan))
                classification = _classify_feature(best_eff, threshold)
                best_hit_target_id = ""
                best_hit_species = ""
                feature_id = str(feature_id)
                if feature_id in best_rows.index:
                    best_hit_target_id = str(best_rows.loc[feature_id, "sseqid"])
                    best_hit_species = str(tip_to_species.get(best_hit_target_id, ""))
                feature_rows.append(
                    {
                        "community_id": community_id,
                        "sample_id": sample_id,
                        "sample_key": sample_key,
                        "method": method,
                        "feature_id": feature_id,
                        "cluster_id": "",
                        "feature_value_tss1000": float(abundance),
                        "near_exact_threshold": threshold,
                        "rule_exact": "best_hit_effective_identity == 100",
                        "rule_near_exact": f"{threshold:g} < best_hit_effective_identity < 100",
                        "rule_no_hit": f"best_hit_effective_identity <= {threshold:g} or no recorded hit",
                        "best_hit_effective_identity": best_eff,
                        "best_hit_target_id": best_hit_target_id,
                        "best_hit_species": best_hit_species,
                        "classification": classification,
                        "classification_source": "feature best hit",
                        "cluster_member_count_total": np.nan,
                        "cluster_member_count_present_in_sample": np.nan,
                        "cluster_members_all": "",
                        "cluster_members_present_in_sample": "",
                        "cluster_members_used_for_classification": "",
                        "cluster_member_classifications_used": "",
                    }
                )

    if not feature_rows:
        raise RuntimeError("No inferred features were found in the Figure 09 registry.")
    return pd.DataFrame(feature_rows)


def _summarize_by_sample(feature_df: pd.DataFrame, registry: pd.DataFrame) -> pd.DataFrame:
    counts = (
        feature_df.groupby(["community_id", "sample_id", "sample_key", "method", "classification"], as_index=False)
        .size()
        .rename(columns={"size": "count"})
    )
    wide = (
        counts.pivot_table(
            index=["community_id", "sample_id", "sample_key", "method"],
            columns="classification",
            values="count",
            fill_value=0,
        )
        .reset_index()
    )
    for col in CLASS_WITH_UNASSIGNED_ORDER:
        if col not in wide.columns:
            wide[col] = 0
    wide = wide.merge(
        registry,
        on=["community_id", "sample_id", "sample_key", "method"],
        how="left",
    )
    wide = wide.rename(
        columns={
            "Exact match": "exact_match_count",
            "Near-exact match": "near_exact_match_count",
            "No hit": "no_hit_count",
            "Unassigned": "unassigned_count",
            "expected_present_features": "expected_present_features_count",
        }
    )
    wide["inferred_present_features_count"] = (
        wide["exact_match_count"]
        + wide["near_exact_match_count"]
        + wide["no_hit_count"]
        + wide["unassigned_count"]
    )
    for label in ["exact_match", "near_exact_match", "no_hit", "unassigned"]:
        num = pd.to_numeric(wide[f"{label}_count"], errors="coerce").fillna(0.0)
        den = pd.to_numeric(wide["inferred_present_features_count"], errors="coerce").replace(0.0, np.nan)
        wide[f"{label}_fraction"] = num / den
    return wide


def _summarize_by_method(sample_df: pd.DataFrame) -> pd.DataFrame:
    summary = (
        sample_df.groupby("method", as_index=False)[
            [
                "expected_present_features_count",
                "inferred_present_features_count",
                "exact_match_count",
                "near_exact_match_count",
                "no_hit_count",
                "unassigned_count",
            ]
        ]
        .sum()
    )
    summary["exact_match_fraction"] = (
        pd.to_numeric(summary["exact_match_count"], errors="coerce")
        / pd.to_numeric(summary["inferred_present_features_count"], errors="coerce").replace(0.0, np.nan)
    )
    summary["near_exact_match_fraction"] = (
        pd.to_numeric(summary["near_exact_match_count"], errors="coerce")
        / pd.to_numeric(summary["inferred_present_features_count"], errors="coerce").replace(0.0, np.nan)
    )
    summary["no_hit_fraction"] = (
        pd.to_numeric(summary["no_hit_count"], errors="coerce")
        / pd.to_numeric(summary["inferred_present_features_count"], errors="coerce").replace(0.0, np.nan)
    )
    summary["unassigned_fraction"] = (
        pd.to_numeric(summary["unassigned_count"], errors="coerce")
        / pd.to_numeric(summary["inferred_present_features_count"], errors="coerce").replace(0.0, np.nan)
    )
    return summary


def _plot_summary(summary: pd.DataFrame, out_svg: Path, font_size: float) -> None:
    summary = summary.set_index("method").reindex(METHOD_ORDER).reset_index()

    xs = np.arange(1, len(METHOD_ORDER) + 1, dtype=float)
    expected = pd.to_numeric(summary["expected_present_features_count"], errors="coerce").fillna(0.0).to_numpy()
    exact = pd.to_numeric(summary["exact_match_count"], errors="coerce").fillna(0.0).to_numpy()
    near = pd.to_numeric(summary["near_exact_match_count"], errors="coerce").fillna(0.0).to_numpy()
    nohit = pd.to_numeric(summary["no_hit_count"], errors="coerce").fillna(0.0).to_numpy()
    unassigned = pd.to_numeric(summary["unassigned_count"], errors="coerce").fillna(0.0).to_numpy()

    fig, ax = plt.subplots(figsize=(10.5, 6.5))
    left_pos = xs - 0.20
    right_pos = xs + 0.20
    width = 0.34

    c_expected = "#009E73"
    c_exact = "#D9D9D9"
    c_near = "#56B4E9"
    c_nohit = "#D55E00"
    c_unassigned = "#8C8C8C"

    ax.bar(
        left_pos,
        expected,
        color=c_expected,
        edgecolor="#333333",
        linewidth=0.7,
        width=width,
        label="Expected species (E)",
    )
    ax.bar(
        right_pos,
        exact,
        color=c_exact,
        edgecolor="#333333",
        linewidth=0.7,
        width=width,
        label="Exact match",
    )
    ax.bar(
        right_pos,
        near,
        bottom=exact,
        color=c_near,
        edgecolor="#333333",
        linewidth=0.7,
        width=width,
        label="Near-exact match",
    )
    ax.bar(
        right_pos,
        nohit,
        bottom=(exact + near),
        color=c_nohit,
        edgecolor="#333333",
        linewidth=0.7,
        width=width,
        label="No hit",
    )
    ax.bar(
        right_pos,
        unassigned,
        bottom=(exact + near + nohit),
        color=c_unassigned,
        edgecolor="#333333",
        linewidth=0.7,
        width=width,
        label="Unassigned",
    )

    ymax = float(
        max(
            np.max(expected) if expected.size else 0.0,
            np.max(exact + near + nohit + unassigned) if exact.size else 0.0,
        )
    )
    pad = max(2.0, 0.03 * max(1.0, ymax))

    for x_l, x_r, exp_ct, ex_ct, nr_ct, nh_ct, ua_ct in zip(
        left_pos, right_pos, expected, exact, near, nohit, unassigned
    ):
        ax.text(
            x_l,
            float(exp_ct) + pad,
            f"E={int(round(exp_ct))}",
            ha="center",
            va="bottom",
            fontsize=max(7.0, font_size - 2.0),
            color="#222222",
        )
        total = float(ex_ct + nr_ct + nh_ct + ua_ct)
        exact_frac = (float(ex_ct) / total) if total > 0 else np.nan
        near_frac = (float(nr_ct) / total) if total > 0 else np.nan
        ax.text(
            x_r,
            total + pad,
            f"X={int(round(ex_ct))}  N={int(round(nr_ct))}  H={int(round(nh_ct))}  U={int(round(ua_ct))}\n"
            f"exact={exact_frac:.2f}  near={near_frac:.2f}" if np.isfinite(exact_frac) and np.isfinite(near_frac) else
            f"X={int(round(ex_ct))}  N={int(round(nr_ct))}  H={int(round(nh_ct))}  U={int(round(ua_ct))}",
            ha="center",
            va="bottom",
            fontsize=max(6.5, font_size - 3.0),
            color="#222222",
        )

    ax.set_xticks(xs)
    ax.set_xticklabels(METHOD_ORDER, fontsize=max(8.0, font_size - 1.0))
    ax.set_ylabel("Feature count (summed across samples)", fontsize=font_size)
    ax.set_xlabel("Method", fontsize=font_size)
    ax.set_title(
        "Figure 09 test: mapping outcome breakdown\n"
        "Exact = 100%; Near-exact = threshold < best hit < 100%; Unassigned = not retained in TIC/TAC",
        fontsize=font_size + 1.0,
    )
    ax.grid(axis="y", color="#E6E6E6", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize=max(8.0, font_size - 1.0))
    fig.tight_layout(rect=[0.0, 0.0, 0.84, 1.0])
    fig.savefig(out_svg, format="svg", bbox_inches="tight")
    plt.close(fig)


def _build_relative_richness_tables(
    fig09_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    diagnostics_path = fig09_dir / "split_merge_diagnostics.tsv"
    if not diagnostics_path.exists():
        raise FileNotFoundError(
            f"Figure 09 diagnostics table not found: {diagnostics_path}. "
            "Run Figure 09 first so this richness test can reuse the same included samples."
        )

    diag = pd.read_csv(diagnostics_path, sep="\t")
    required = [
        "community_id",
        "sample_id",
        "sample_key",
        "method",
        "inferred_present_features",
        "expected_present_features",
    ]
    missing = [c for c in required if c not in diag.columns]
    if missing:
        raise ValueError(
            "Figure 09 diagnostics table is missing required richness columns: "
            + ", ".join(missing)
        )

    sample_df = diag[required].drop_duplicates(
        subset=["community_id", "sample_id", "sample_key", "method"]
    ).copy()
    sample_df = sample_df[sample_df["method"].astype(str).isin(METHOD_ORDER)].copy()
    sample_df["method"] = sample_df["method"].astype(str)
    sample_df["method_label"] = sample_df["method"].map(display_method)
    sample_df["expected_species_count"] = pd.to_numeric(
        sample_df["expected_present_features"], errors="coerce"
    )
    sample_df["inferred_feature_count"] = pd.to_numeric(
        sample_df["inferred_present_features"], errors="coerce"
    )
    sample_df = sample_df.dropna(subset=["expected_species_count", "inferred_feature_count"])
    sample_df = sample_df[sample_df["expected_species_count"] > 0].copy()
    sample_df["relative_richness_error"] = (
        (sample_df["inferred_feature_count"] - sample_df["expected_species_count"])
        / sample_df["expected_species_count"]
    )
    sample_df["relative_richness_error_percent"] = 100.0 * sample_df["relative_richness_error"]
    sample_df["richness_difference"] = (
        sample_df["inferred_feature_count"] - sample_df["expected_species_count"]
    )
    sample_df["richness_error_direction"] = np.select(
        [
            sample_df["richness_difference"] > 0,
            sample_df["richness_difference"] < 0,
        ],
        ["overestimated", "underestimated"],
        default="equal",
    )
    sample_df = sample_df[
        [
            "community_id",
            "sample_id",
            "sample_key",
            "method",
            "method_label",
            "expected_species_count",
            "inferred_feature_count",
            "richness_difference",
            "relative_richness_error",
            "relative_richness_error_percent",
            "richness_error_direction",
        ]
    ].copy()

    method_rows: list[dict[str, object]] = []
    for method in METHOD_ORDER:
        sub = sample_df[sample_df["method"] == method].copy()
        vals = pd.to_numeric(sub["relative_richness_error_percent"], errors="coerce").dropna()
        diffs = pd.to_numeric(sub["richness_difference"], errors="coerce").dropna()
        method_rows.append(
            {
                "method": method,
                "method_label": display_method(method),
                "n_samples": int(vals.size),
                "median_relative_richness_error_percent": float(vals.median()) if vals.size else np.nan,
                "mean_relative_richness_error_percent": float(vals.mean()) if vals.size else np.nan,
                "q1_relative_richness_error_percent": float(vals.quantile(0.25)) if vals.size else np.nan,
                "q3_relative_richness_error_percent": float(vals.quantile(0.75)) if vals.size else np.nan,
                "min_relative_richness_error_percent": float(vals.min()) if vals.size else np.nan,
                "max_relative_richness_error_percent": float(vals.max()) if vals.size else np.nan,
                "median_richness_difference": float(diffs.median()) if diffs.size else np.nan,
                "mean_richness_difference": float(diffs.mean()) if diffs.size else np.nan,
                "fraction_overestimated": float((diffs > 0).mean()) if diffs.size else np.nan,
                "fraction_underestimated": float((diffs < 0).mean()) if diffs.size else np.nan,
                "fraction_equal": float((diffs == 0).mean()) if diffs.size else np.nan,
            }
        )
    method_summary = pd.DataFrame(method_rows)

    count_wide = sample_df.pivot_table(
        index=["community_id", "sample_id", "sample_key"],
        columns="method",
        values="inferred_feature_count",
        aggfunc="first",
    )
    expected_wide = sample_df.pivot_table(
        index=["community_id", "sample_id", "sample_key"],
        values="expected_species_count",
        aggfunc="first",
    )
    wide = count_wide.join(expected_wide.rename(columns={"expected_species_count": "expected_species_count"}))
    effect_rows: list[dict[str, object]] = []
    for idx, row in wide.iterrows():
        if "zOTU" not in row.index or not np.isfinite(float(row.get("zOTU", np.nan))):
            continue
        expected_count = float(row.get("expected_species_count", np.nan))
        zotu_count = float(row.get("zOTU", np.nan))
        excess_unoise3 = zotu_count - expected_count
        if not np.isfinite(excess_unoise3) or excess_unoise3 <= 0.0:
            continue
        community_id, sample_id, sample_key = idx
        for method in ("TAC", "TIC"):
            method_count = float(row.get(method, np.nan))
            if not np.isfinite(method_count):
                continue
            excess_removed = zotu_count - method_count
            effect_rows.append(
                {
                    "community_id": community_id,
                    "sample_id": sample_id,
                    "sample_key": sample_key,
                    "clustering_method": method,
                    "clustering_method_label": display_method(method),
                    "expected_species_count": expected_count,
                    "unoise3_feature_count": zotu_count,
                    "clustering_feature_count": method_count,
                    "unoise3_excess_richness": excess_unoise3,
                    "excess_removed": excess_removed,
                    "percent_excess_removed": 100.0 * excess_removed / excess_unoise3,
                    "overcorrected_below_expected": bool(method_count < expected_count),
                }
            )
    effect_df = pd.DataFrame(effect_rows)

    effect_summary_rows: list[dict[str, object]] = []
    for method in ("TAC", "TIC"):
        sub = effect_df[effect_df["clustering_method"] == method].copy()
        vals = pd.to_numeric(sub["percent_excess_removed"], errors="coerce").dropna()
        effect_summary_rows.append(
            {
                "clustering_method": method,
                "clustering_method_label": display_method(method),
                "n_samples_with_unoise3_excess": int(vals.size),
                "median_percent_excess_removed": float(vals.median()) if vals.size else np.nan,
                "mean_percent_excess_removed": float(vals.mean()) if vals.size else np.nan,
                "q1_percent_excess_removed": float(vals.quantile(0.25)) if vals.size else np.nan,
                "q3_percent_excess_removed": float(vals.quantile(0.75)) if vals.size else np.nan,
                "min_percent_excess_removed": float(vals.min()) if vals.size else np.nan,
                "max_percent_excess_removed": float(vals.max()) if vals.size else np.nan,
                "fraction_overcorrected_below_expected": float(
                    sub["overcorrected_below_expected"].astype(bool).mean()
                )
                if not sub.empty
                else np.nan,
            }
        )
    effect_summary = pd.DataFrame(effect_summary_rows)
    return sample_df, method_summary, effect_df, effect_summary


def _plot_relative_richness_error_figure(
    sample_df: pd.DataFrame,
    effect_summary_df: pd.DataFrame,
    out_svg: Path,
    font_size: float,
) -> None:
    plot_df = sample_df[sample_df["method"].isin(RICHNESS_METHOD_ORDER)].copy()
    fig, ax = plt.subplots(figsize=(9.8, 6.4))
    rng = np.random.default_rng(7)
    positions = np.arange(1, len(RICHNESS_METHOD_ORDER) + 1, dtype=float)
    data = [
        pd.to_numeric(
            plot_df.loc[plot_df["method"] == method, "relative_richness_error_percent"],
            errors="coerce",
        )
        .dropna()
        .to_numpy(dtype=float)
        for method in RICHNESS_METHOD_ORDER
    ]

    box = ax.boxplot(
        data,
        positions=positions,
        widths=0.52,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="#111111", linewidth=1.5),
        whiskerprops=dict(color="#333333", linewidth=1.0),
        capprops=dict(color="#333333", linewidth=1.0),
    )
    for patch, method in zip(box["boxes"], RICHNESS_METHOD_ORDER):
        patch.set_facecolor(METHOD_COLORS.get(method, "#999999"))
        patch.set_alpha(0.38)
        patch.set_edgecolor("#333333")
        patch.set_linewidth(1.0)

    for x, method, vals in zip(positions, RICHNESS_METHOD_ORDER, data):
        if vals.size == 0:
            continue
        jitter = rng.normal(loc=0.0, scale=0.055, size=vals.size)
        ax.scatter(
            np.full(vals.size, x) + jitter,
            vals,
            s=34,
            alpha=0.72,
            color=METHOD_COLORS.get(method, "#666666"),
            edgecolors="#333333",
            linewidths=0.35,
            zorder=3,
        )

    ax.axhline(0.0, color="#333333", linestyle=(0, (4, 3)), linewidth=1.1)
    ax.text(
        0.99,
        0.02,
        "0% = expected species richness",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=max(7.0, font_size - 2.5),
        color="#333333",
    )

    effect_lines: list[str] = []
    for method in ("TAC", "TIC"):
        hit = effect_summary_df[effect_summary_df["clustering_method"] == method]
        if hit.empty:
            continue
        val = float(hit["median_percent_excess_removed"].iloc[0])
        n = int(hit["n_samples_with_unoise3_excess"].iloc[0])
        if np.isfinite(val):
            effect_lines.append(f"{display_method(method)} median excess removed: {val:.1f}% (n={n})")
    if effect_lines:
        ax.text(
            0.02,
            0.98,
            "\n".join(effect_lines),
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=max(7.2, font_size - 2.2),
            bbox=dict(boxstyle="round,pad=0.32", facecolor="white", edgecolor="#BBBBBB", alpha=0.88),
        )

    ax.set_xticks(positions)
    ax.set_xticklabels([display_method(m) for m in RICHNESS_METHOD_ORDER], fontsize=max(8.0, font_size - 1.0))
    ax.set_ylabel("Relative richness error (%)", fontsize=font_size)
    ax.set_xlabel("Method", fontsize=font_size)
    ax.set_title(
        "Figure 09 test: relative richness error against expected species richness",
        fontsize=font_size + 1.0,
    )
    ax.grid(axis="y", color="#E6E6E6", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", labelsize=max(8.0, font_size - 1.0))

    finite_vals = np.concatenate([v[np.isfinite(v)] for v in data if v.size > 0]) if any(v.size for v in data) else np.array([])
    if finite_vals.size:
        ymin = float(np.nanmin(finite_vals))
        ymax = float(np.nanmax(finite_vals))
        span = max(1.0, ymax - ymin)
        ax.set_ylim(ymin - 0.12 * span, ymax + 0.18 * span)

    fig.tight_layout()
    fig.savefig(out_svg, format="svg", bbox_inches="tight")
    plt.close(fig)


def _build_transition_tables(
    feature_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    z_df = feature_df[feature_df["method"] == "zOTU"].copy()
    z_df = z_df.rename(
        columns={
            "feature_id": "zotu_id",
            "classification": "zotu_class",
            "best_hit_effective_identity": "zotu_best_hit_effective_identity",
        }
    )
    z_df = z_df[
        [
            "community_id",
            "sample_id",
            "sample_key",
            "zotu_id",
            "zotu_class",
            "zotu_best_hit_effective_identity",
            "feature_value_tss1000",
        ]
    ].drop_duplicates()

    pri = {"No hit": 0, "Near-exact match": 1, "Exact match": 2}
    inv = {v: k for k, v in pri.items()}

    transition_rows: list[pd.DataFrame] = []
    near_rows: list[pd.DataFrame] = []
    flow_rows: list[pd.DataFrame] = []

    for method in ["TIC", "TAC"]:
        method_df = feature_df[feature_df["method"] == method].copy()
        if method_df.empty:
            continue
        assigned_df = method_df[
            method_df["cluster_id"].astype(str).str.strip().ne("")
            & method_df["classification"].astype(str).ne("Unassigned")
        ].copy()
        assigned_df = assigned_df.rename(columns={"feature_id": "zotu_id"})
        merged = assigned_df.merge(
            z_df,
            on=["community_id", "sample_id", "sample_key", "zotu_id"],
            how="left",
            suffixes=("", "_zotu"),
        )
        cluster_class = pd.DataFrame(
            columns=["community_id", "sample_id", "sample_key", "cluster_id", "cluster_class"]
        )
        if not merged.empty:
            cluster_class = (
                merged.assign(_score=merged["zotu_class"].map(pri).fillna(0).astype(int))
                .groupby(["community_id", "sample_id", "sample_key", "cluster_id"], as_index=False)["_score"]
                .max()
            )
            cluster_class["cluster_class"] = cluster_class["_score"].map(inv)

            merged = merged.merge(
                cluster_class[["community_id", "sample_id", "sample_key", "cluster_id", "cluster_class"]],
                on=["community_id", "sample_id", "sample_key", "cluster_id"],
                how="left",
            )
            merged["transition"] = merged["zotu_class"].astype(str) + " -> " + merged["cluster_class"].astype(str)
            merged["transition_method"] = method
            transition_rows.append(
                merged[
                    [
                        "community_id",
                        "sample_id",
                        "sample_key",
                        "transition_method",
                        "zotu_id",
                        "cluster_id",
                        "zotu_class",
                        "cluster_class",
                        "transition",
                        "zotu_best_hit_effective_identity",
                        "feature_value_tss1000_zotu",
                        "feature_value_tss1000",
                    ]
                ].rename(
                    columns={
                        "feature_value_tss1000_zotu": "zotu_feature_value_tss1000",
                        "feature_value_tss1000": f"{method.lower()}_cluster_member_value_tss1000",
                    }
                )
            )

        assigned_keys = (
            assigned_df.merge(
                cluster_class[["community_id", "sample_id", "sample_key", "cluster_id", "cluster_class"]],
                on=["community_id", "sample_id", "sample_key", "cluster_id"],
                how="left",
            )[
                ["community_id", "sample_id", "sample_key", "zotu_id", "cluster_class"]
            ]
            .drop_duplicates()
        )
        near_all = z_df[z_df["zotu_class"] == "Near-exact match"][
            [
                "community_id",
                "sample_id",
                "sample_key",
                "zotu_id",
                "zotu_class",
                "zotu_best_hit_effective_identity",
                "feature_value_tss1000",
            ]
        ].copy()
        near_join = near_all.merge(
            assigned_keys,
            on=["community_id", "sample_id", "sample_key", "zotu_id"],
            how="left",
        )
        near_join["transition_method"] = method
        near_join["near_exact_fate"] = near_join["cluster_class"].fillna("Unassigned")
        near_rows.append(
            near_join.rename(columns={"feature_value_tss1000": "zotu_feature_value_tss1000"})
        )

        flow_join = z_df.merge(
            assigned_keys.rename(columns={"cluster_class": "final_fate"}),
            on=["community_id", "sample_id", "sample_key", "zotu_id"],
            how="left",
        )
        flow_join["transition_method"] = method
        flow_join["final_fate"] = flow_join["final_fate"].fillna("Unassigned")
        flow_rows.append(
            flow_join.rename(columns={"feature_value_tss1000": "zotu_feature_value_tss1000"})
        )

    if transition_rows:
        transition_df = pd.concat(transition_rows, ignore_index=True)
    else:
        transition_df = pd.DataFrame(
            columns=[
                "community_id",
                "sample_id",
                "sample_key",
                "transition_method",
                "zotu_id",
                "cluster_id",
                "zotu_class",
                "cluster_class",
                "transition",
                "zotu_best_hit_effective_identity",
                "zotu_feature_value_tss1000",
            ]
        )

    if near_rows:
        near_df = pd.concat(near_rows, ignore_index=True)
    else:
        near_df = pd.DataFrame(
            columns=[
                "community_id",
                "sample_id",
                "sample_key",
                "zotu_id",
                "zotu_class",
                "zotu_best_hit_effective_identity",
                "zotu_feature_value_tss1000",
                "transition_method",
                "cluster_class",
                "near_exact_fate",
            ]
        )

    if flow_rows:
        flow_df = pd.concat(flow_rows, ignore_index=True)
    else:
        flow_df = pd.DataFrame(
            columns=[
                "community_id",
                "sample_id",
                "sample_key",
                "zotu_id",
                "zotu_class",
                "zotu_best_hit_effective_identity",
                "zotu_feature_value_tss1000",
                "transition_method",
                "final_fate",
            ]
        )

    matrix_records: list[dict[str, object]] = []
    for method in ["TIC", "TAC"]:
        sub = transition_df[transition_df["transition_method"] == method].copy()
        for row_class in CLASS_ORDER:
            row_total = int((sub["zotu_class"] == row_class).sum())
            for col_class in CLASS_ORDER:
                count = int(((sub["zotu_class"] == row_class) & (sub["cluster_class"] == col_class)).sum())
                matrix_records.append(
                    {
                        "transition_method": method,
                        "zotu_class": row_class,
                        "cluster_class": col_class,
                        "count": count,
                        "row_total_assigned": row_total,
                        "row_fraction_assigned": (float(count) / float(row_total)) if row_total > 0 else np.nan,
                    }
                )
    matrix_df = pd.DataFrame(matrix_records)

    near_summary = (
        near_df.groupby(["transition_method", "near_exact_fate"], as_index=False)
        .size()
        .rename(columns={"size": "count"})
    )
    total_by_method = (
        near_df.groupby("transition_method", as_index=False)
        .size()
        .rename(columns={"size": "near_exact_total"})
    )
    assigned_by_method = (
        near_df[near_df["near_exact_fate"] != "Unassigned"]
        .groupby("transition_method", as_index=False)
        .size()
        .rename(columns={"size": "near_exact_assigned"})
    )
    near_summary = near_summary.merge(total_by_method, on="transition_method", how="left")
    near_summary = near_summary.merge(assigned_by_method, on="transition_method", how="left")
    near_summary["near_exact_assigned"] = pd.to_numeric(
        near_summary["near_exact_assigned"], errors="coerce"
    ).fillna(0).astype(int)
    near_summary["fraction_of_all_near_exact"] = np.where(
        near_summary["near_exact_total"] > 0,
        pd.to_numeric(near_summary["count"], errors="coerce") / pd.to_numeric(near_summary["near_exact_total"], errors="coerce"),
        np.nan,
    )
    near_summary["fraction_of_assigned_near_exact"] = np.where(
        (near_summary["near_exact_assigned"] > 0) & (near_summary["near_exact_fate"] != "Unassigned"),
        pd.to_numeric(near_summary["count"], errors="coerce") / pd.to_numeric(near_summary["near_exact_assigned"], errors="coerce"),
        np.nan,
    )

    return transition_df, matrix_df, near_summary, flow_df


def _draw_transition_heatmap(ax, matrix_df: pd.DataFrame, method: str, font_size: float) -> None:
    sub = matrix_df[matrix_df["transition_method"] == method].copy()
    heat = np.zeros((len(CLASS_ORDER), len(CLASS_ORDER)), dtype=float)
    counts = np.zeros_like(heat, dtype=int)
    for i, row_class in enumerate(CLASS_ORDER):
        for j, col_class in enumerate(CLASS_ORDER):
            hit = sub[(sub["zotu_class"] == row_class) & (sub["cluster_class"] == col_class)]
            if not hit.empty:
                heat[i, j] = float(hit["row_fraction_assigned"].iloc[0]) if np.isfinite(float(hit["row_fraction_assigned"].iloc[0])) else 0.0
                counts[i, j] = int(hit["count"].iloc[0])

    im = ax.imshow(heat, cmap="YlGnBu", vmin=0.0, vmax=1.0, aspect="auto")
    ax.set_xticks(np.arange(len(CLASS_ORDER)))
    ax.set_yticks(np.arange(len(CLASS_ORDER)))
    ax.set_xticklabels(CLASS_ORDER, rotation=25, ha="right", fontsize=max(7.0, font_size - 2.0))
    ax.set_yticklabels(CLASS_ORDER, fontsize=max(7.0, font_size - 2.0))
    ax.set_xlabel("Final cluster class", fontsize=max(8.0, font_size - 1.0))
    ax.set_ylabel("Original zOTU member class", fontsize=max(8.0, font_size - 1.0))
    ax.set_title(f"{method}: zOTU-to-cluster transitions\n(assigned members only)", fontsize=font_size)

    for i in range(len(CLASS_ORDER)):
        for j in range(len(CLASS_ORDER)):
            frac = heat[i, j]
            count = counts[i, j]
            text = f"{count}\n{frac * 100:.1f}%"
            color = "white" if frac >= 0.55 else "#1F1F1F"
            ax.text(j, i, text, ha="center", va="center", fontsize=max(6.5, font_size - 3.0), color=color)

    ax.set_xticks(np.arange(-0.5, len(CLASS_ORDER), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(CLASS_ORDER), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="minor", bottom=False, left=False)
    return im


def _plot_transition_figure(matrix_df: pd.DataFrame, near_summary: pd.DataFrame, out_svg: Path, font_size: float) -> None:
    fig = plt.figure(figsize=(14.0, 5.8))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.2], wspace=0.5)
    ax_tic = fig.add_subplot(gs[0, 0])
    ax_tac = fig.add_subplot(gs[0, 1])
    ax_bar = fig.add_subplot(gs[0, 2])

    im_tic = _draw_transition_heatmap(ax_tic, matrix_df, "TIC", font_size)
    _draw_transition_heatmap(ax_tac, matrix_df, "TAC", font_size)

    bar_methods = ["TIC", "TAC"]
    fate_order = ["Exact match", "Near-exact match", "No hit", "Unassigned"]
    fate_labels = {
        "Exact match": "Upgraded to exact",
        "Near-exact match": "Stayed near-exact",
        "No hit": "Ended as no hit",
        "Unassigned": "Unassigned to present cluster",
    }
    fate_colors = {
        "Exact match": "#009E73",
        "Near-exact match": "#56B4E9",
        "No hit": "#D55E00",
        "Unassigned": "#BDBDBD",
    }

    xs = np.arange(len(bar_methods), dtype=float)
    bottom = np.zeros(len(bar_methods), dtype=float)
    for fate in fate_order:
        vals = []
        for method in bar_methods:
            hit = near_summary[
                (near_summary["transition_method"] == method) & (near_summary["near_exact_fate"] == fate)
            ]
            frac = float(hit["fraction_of_all_near_exact"].iloc[0]) if not hit.empty else 0.0
            vals.append(frac)
        vals_arr = np.asarray(vals, dtype=float)
        ax_bar.bar(
            xs,
            vals_arr,
            bottom=bottom,
            color=fate_colors[fate],
            edgecolor="#333333",
            linewidth=0.7,
            width=0.62,
            label=fate_labels[fate],
        )
        for x, y0, h, method in zip(xs, bottom, vals_arr, bar_methods):
            if h <= 0:
                continue
            hit = near_summary[
                (near_summary["transition_method"] == method) & (near_summary["near_exact_fate"] == fate)
            ]
            n = int(hit["count"].iloc[0]) if not hit.empty else 0
            ax_bar.text(
                x,
                y0 + (h / 2.0),
                f"{n}\n{h * 100:.1f}%",
                ha="center",
                va="center",
                fontsize=max(6.2, font_size - 3.5),
                color="white" if h >= 0.18 else "#1F1F1F",
            )
        bottom += vals_arr

    for x, method in zip(xs, bar_methods):
        hit = near_summary[near_summary["transition_method"] == method]
        total = int(hit["near_exact_total"].iloc[0]) if not hit.empty else 0
        assigned = int(hit["near_exact_assigned"].iloc[0]) if not hit.empty else 0
        ax_bar.text(
            x,
            1.03,
            f"n={total}\nassigned={assigned}",
            ha="center",
            va="bottom",
            fontsize=max(6.5, font_size - 3.0),
            color="#222222",
        )

    ax_bar.set_xticks(xs)
    ax_bar.set_xticklabels(bar_methods, fontsize=max(8.0, font_size - 1.0))
    ax_bar.set_ylim(0.0, 1.12)
    ax_bar.set_ylabel("Fraction of all near-exact zOTUs", fontsize=max(8.0, font_size - 1.0))
    ax_bar.set_xlabel("Method", fontsize=max(8.0, font_size - 1.0))
    ax_bar.set_title("Fate of near-exact zOTUs", fontsize=font_size)
    ax_bar.grid(axis="y", color="#E6E6E6", linewidth=0.8)
    ax_bar.set_axisbelow(True)
    ax_bar.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False, fontsize=max(7.0, font_size - 2.0))

    cbar = fig.colorbar(im_tic, ax=[ax_tic, ax_tac], fraction=0.03, pad=0.03)
    cbar.ax.set_ylabel("Row fraction among assigned members", fontsize=max(7.0, font_size - 2.0))
    cbar.ax.tick_params(labelsize=max(6.5, font_size - 3.0))

    fig.suptitle(
        "Figure 09 test: absorption of near-exact zOTUs by TIC/TAC clustering",
        fontsize=font_size + 2.0,
        y=1.02,
    )
    fig.tight_layout(rect=[0.0, 0.0, 0.9, 0.97])
    fig.savefig(out_svg, format="svg", bbox_inches="tight")
    plt.close(fig)


def _build_flow_summary(flow_df: pd.DataFrame) -> pd.DataFrame:
    summary = (
        flow_df.groupby(["transition_method", "zotu_class", "final_fate"], as_index=False)
        .size()
        .rename(columns={"size": "count"})
    )
    records: list[dict[str, object]] = []
    for method in ["TIC", "TAC"]:
        sub = summary[summary["transition_method"] == method].copy()
        for origin in CLASS_ORDER:
            row_total = int(sub.loc[sub["zotu_class"] == origin, "count"].sum())
            for fate in CLASS_WITH_UNASSIGNED_ORDER:
                hit = sub[(sub["zotu_class"] == origin) & (sub["final_fate"] == fate)]
                count = int(hit["count"].iloc[0]) if not hit.empty else 0
                records.append(
                    {
                        "transition_method": method,
                        "zotu_class": origin,
                        "final_fate": fate,
                        "count": count,
                        "row_total": row_total,
                        "row_fraction": (float(count) / float(row_total)) if row_total > 0 else np.nan,
                    }
                )
    return pd.DataFrame(records)


def _plot_flow_transition_heatmap_figure(
    flow_summary: pd.DataFrame,
    out_svg: Path,
    font_size: float,
) -> None:
    row_labels_raw = CLASS_ORDER
    col_labels_raw = CLASS_WITH_UNASSIGNED_ORDER
    row_labels = [_wave_display_label(label) for label in row_labels_raw]
    col_labels = [_wave_display_label(label) for label in col_labels_raw]

    methods = ["TIC", "TAC"]
    method_matrices: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for method in methods:
        sub = flow_summary[flow_summary["transition_method"] == method].copy()
        heat = np.zeros((len(row_labels_raw), len(col_labels_raw)), dtype=float)
        counts = np.zeros_like(heat, dtype=int)
        row_totals = np.zeros(len(row_labels_raw), dtype=int)
        for i, origin in enumerate(row_labels_raw):
            origin_sub = sub[sub["zotu_class"] == origin]
            row_total = int(pd.to_numeric(origin_sub["count"], errors="coerce").fillna(0).sum())
            row_totals[i] = row_total
            for j, fate in enumerate(col_labels_raw):
                hit = origin_sub[origin_sub["final_fate"] == fate]
                count = int(pd.to_numeric(hit["count"], errors="coerce").fillna(0).sum()) if not hit.empty else 0
                counts[i, j] = count
                heat[i, j] = (float(count) / float(row_total)) if row_total > 0 else 0.0
        method_matrices[method] = (heat, counts, row_totals)

    fig, ax = plt.subplots(figsize=(10.8, 5.4))
    cmap = plt.get_cmap("YlGnBu")
    norm = plt.Normalize(vmin=0.0, vmax=1.0)
    half_specs = {
        "TIC": (0.0, 0.5, "top half"),
        "TAC": (0.5, 0.5, "bottom half"),
    }

    for i, _origin in enumerate(row_labels_raw):
        for j, _fate in enumerate(col_labels_raw):
            for method in methods:
                heat, counts, _row_totals = method_matrices[method]
                y_offset, height, _half_label = half_specs[method]
                frac = float(heat[i, j])
                count = int(counts[i, j])
                ax.add_patch(
                    Rectangle(
                        (j, i + y_offset),
                        1.0,
                        height,
                        facecolor=cmap(norm(frac)),
                        edgecolor="white",
                        linewidth=1.2,
                    )
                )
                color = "white" if frac >= 0.55 else "#222222"
                ax.text(
                    j + 0.5,
                    i + y_offset + (height / 2.0),
                    f"{count}\n{frac * 100:.1f}%",
                    ha="center",
                    va="center",
                    fontsize=max(6.0, font_size - 4.8),
                    color=color,
                    linespacing=0.95,
                )

    for x in range(len(col_labels) + 1):
        ax.axvline(float(x), color="#333333", linewidth=0.7)
    for y in range(len(row_labels) + 1):
        ax.axhline(float(y), color="#333333", linewidth=0.7)
    for y in range(len(row_labels)):
        ax.axhline(float(y) + 0.5, color="white", linewidth=1.1)

    tic_totals = method_matrices["TIC"][2]
    tac_totals = method_matrices["TAC"][2]
    row_tick_labels = []
    for i, label in enumerate(row_labels):
        if int(tic_totals[i]) == int(tac_totals[i]):
            row_tick_labels.append(f"{label}\n(n={int(tic_totals[i])})")
        else:
            row_tick_labels.append(f"{label}\n(TIC={int(tic_totals[i])}, TAC={int(tac_totals[i])})")

    ax.set_xlim(0.0, float(len(col_labels)))
    ax.set_ylim(float(len(row_labels)), 0.0)
    ax.set_aspect("equal")
    ax.set_xticks(np.arange(len(col_labels)) + 0.5)
    ax.set_xticklabels(col_labels, rotation=35, ha="right", rotation_mode="anchor", fontsize=max(7.0, font_size - 3.0))
    ax.set_yticks(np.arange(len(row_labels)) + 0.5)
    ax.set_yticklabels(row_tick_labels, fontsize=max(7.0, font_size - 3.0))
    for i in range(len(row_labels_raw)):
        ax.text(
            -0.02,
            i + 0.25,
            "TIC",
            transform=ax.get_yaxis_transform(),
            ha="right",
            va="center",
            fontsize=max(6.2, font_size - 4.0),
            fontweight="bold",
            color="#333333",
        )
        ax.text(
            -0.02,
            i + 0.75,
            "TAC",
            transform=ax.get_yaxis_transform(),
            ha="right",
            va="center",
            fontsize=max(6.2, font_size - 4.0),
            fontweight="bold",
            color="#333333",
        )
    ax.set_xlabel("Post-clustering outcome", fontsize=max(8.0, font_size - 2.0))
    ax.set_ylabel("Original zOTU class", fontsize=max(8.0, font_size - 2.0))
    ax.tick_params(axis="both", length=0)

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.03)
    cbar.ax.set_ylabel("Row fraction of zOTUs", fontsize=max(8.0, font_size - 2.0))
    cbar.ax.tick_params(labelsize=max(7.0, font_size - 3.0))
    ax.text(
        0.5,
        1.04,
        "Each cell is split into TIC (top) and TAC (bottom).",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=max(8.0, font_size - 2.0),
        color="#333333",
    )

    fig.suptitle(
        "Figure 09 test: zOTU outcome transitions after TIC/TAC clustering",
        fontsize=font_size + 1.5,
        fontweight="bold",
        y=0.995,
    )
    fig.subplots_adjust(left=0.16, right=0.86, bottom=0.24, top=0.80)
    fig.savefig(out_svg, format="svg", bbox_inches="tight")
    plt.close(fig)


def _plot_wave_transition_figure(flow_summary: pd.DataFrame, out_html: Path, font_size: float) -> None:
    from floweaver import Bundle, Dataset, Partition, ProcessGroup, SankeyDefinition, weave
    from floweaver import CategoricalScale
    from ipywidgets import HTML, HBox, Layout, VBox
    from ipywidgets.embed import embed_minimal_html

    color_by_origin = {
        "Exact match": "#BDBDBD",
        "Near-exact match": "#56B4E9",
        "No hit": "#D55E00",
    }

    def _build_method_widget(method: str):
        sub = flow_summary[flow_summary["transition_method"] == method].copy()
        sub = sub[pd.to_numeric(sub["count"], errors="coerce").fillna(0.0) > 0].copy()
        if sub.empty:
            return VBox([HTML(f"<b>{method}</b>"), HTML("<i>No flow data available.</i>")])

        source_labels_raw = [label for label in CLASS_ORDER if float(sub.loc[sub["zotu_class"] == label, "count"].sum()) > 0]
        target_labels_raw = [
            label
            for label in CLASS_WITH_UNASSIGNED_ORDER
            if float(sub.loc[sub["final_fate"] == label, "count"].sum()) > 0
        ]

        flows_rows: list[dict[str, object]] = []
        for origin in source_labels_raw:
            for fate in target_labels_raw:
                hit = sub[(sub["zotu_class"] == origin) & (sub["final_fate"] == fate)]
                count = float(hit["count"].iloc[0]) if not hit.empty else 0.0
                if count <= 0:
                    continue
                flows_rows.append(
                    {
                        "source": _wave_display_label(origin),
                        "target": _wave_display_label(fate),
                        "value": float(count),
                        "origin_kind": _wave_display_label(origin),
                    }
                )

        if not flows_rows:
            return VBox([HTML(f"<b>{method}</b>"), HTML("<i>No non-zero flows available.</i>")])

        flows = pd.DataFrame(flows_rows, columns=["source", "target", "value", "origin_kind"])
        source_labels = [_wave_display_label(label) for label in source_labels_raw]
        target_labels = [_wave_display_label(label) for label in target_labels_raw]

        nodes = {
            "src": ProcessGroup(
                selection=source_labels,
                partition=Partition.Simple("source", source_labels),
                title="Original zOTU class",
            ),
            "dst": ProcessGroup(
                selection=target_labels,
                partition=Partition.Simple("target", target_labels),
                title="Outcome after clustering",
            ),
        }
        bundles = [Bundle("src", "dst")]
        definition = SankeyDefinition(nodes, bundles, [["src"], ["dst"]])
        color_scale = CategoricalScale(
            "source",
            palette=[
                color_by_origin["Exact match"],
                color_by_origin["Near-exact match"],
                color_by_origin["No hit"],
            ],
        )
        sankey = weave(definition, Dataset(flows), measures="value", link_color=color_scale)
        widget = sankey.to_widget(
            width=720,
            height=500,
            margins={"top": 15, "bottom": 10, "left": 160, "right": 160},
            link_label_format=".0f",
            link_label_min_width=8,
        )
        title_html = HTML(
            value=f"<div style='font-size:{max(12, int(round(font_size)))}px; "
            "font-weight:600; text-align:center; margin-bottom:4px;'>"
            f"{method}</div>"
        )
        return VBox([title_html, widget], layout=Layout(width="760px"))

    title = HTML(
        value=(
            f"<div style='font-size:{max(14, int(round(font_size + 2.0)))}px; "
            "font-weight:700; text-align:center; margin:8px 0 10px 0;'>"
            "Figure 09 test: wave transition view of how TIC/TAC absorb zOTU classes"
            "</div>"
        )
    )
    subtitle = HTML(
        value=(
            f"<div style='font-size:{max(10, int(round(font_size - 1.0)))}px; "
            "text-align:center; margin:0 0 8px 0;'>"
            "TIC includes the Non-bacterial outcome for zOTUs removed before final cluster membership."
            "</div>"
        )
    )
    panel = VBox(
        [
            title,
            subtitle,
            HBox(
                [_build_method_widget("TIC"), _build_method_widget("TAC")],
                layout=Layout(justify_content="center", align_items="flex-start", gap="18px"),
            ),
        ],
        layout=Layout(width="100%"),
    )
    out_html.parent.mkdir(parents=True, exist_ok=True)
    embed_minimal_html(str(out_html), views=[panel], title="Figure 09 test wave transition")


def _add_wave_ribbon(
    ax,
    x_left: float,
    x_right: float,
    y0_bottom: float,
    y0_top: float,
    y1_bottom: float,
    y1_top: float,
    color: str,
    alpha: float,
) -> None:
    ctrl_dx = (x_right - x_left) * 0.42
    verts = [
        (x_left, y0_top),
        (x_left + ctrl_dx, y0_top),
        (x_right - ctrl_dx, y1_top),
        (x_right, y1_top),
        (x_right, y1_bottom),
        (x_right - ctrl_dx, y1_bottom),
        (x_left + ctrl_dx, y0_bottom),
        (x_left, y0_bottom),
        (x_left, y0_top),
    ]
    codes = [
        MplPath.MOVETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.LINETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CLOSEPOLY,
    ]
    patch = PathPatch(
        MplPath(verts, codes),
        facecolor=color,
        edgecolor="none",
        alpha=alpha,
        zorder=1,
    )
    ax.add_patch(patch)


def _add_wave_arrowhead(
    ax,
    x_center: float,
    y_center: float,
    direction: str,
    color: str,
    size_x: float,
    size_y: float,
) -> None:
    if direction == "left":
        points = [
            (x_center - (size_x / 2.0), y_center),
            (x_center + (size_x / 2.0), y_center + (size_y / 2.0)),
            (x_center + (size_x / 2.0), y_center - (size_y / 2.0)),
        ]
    else:
        points = [
            (x_center + (size_x / 2.0), y_center),
            (x_center - (size_x / 2.0), y_center + (size_y / 2.0)),
            (x_center - (size_x / 2.0), y_center - (size_y / 2.0)),
        ]
    ax.add_patch(
        Polygon(
            points,
            closed=True,
            facecolor=color,
            edgecolor="#FFFFFF",
            linewidth=0.45,
            alpha=0.95,
            zorder=2.2,
        )
    )


def _draw_wave_transition_panel(ax, flow_summary: pd.DataFrame, method: str, font_size: float) -> None:
    color_by_origin = {
        "Exact match": "#BDBDBD",
        "Near-exact match": "#56B4E9",
        "No match": "#D55E00",
    }
    display_order = [_wave_display_label(label) for label in CLASS_ORDER]
    target_display_order = [_wave_display_label(label) for label in CLASS_WITH_UNASSIGNED_ORDER]

    sub = flow_summary[flow_summary["transition_method"] == method].copy()
    sub = sub[pd.to_numeric(sub["count"], errors="coerce").fillna(0.0) > 0].copy()
    if sub.empty:
        ax.text(0.5, 0.5, "No flow data available", ha="center", va="center", fontsize=font_size)
        ax.axis("off")
        return

    if "zotu_class_display" not in sub.columns:
        sub["zotu_class_display"] = sub["zotu_class"].astype(str).map(_wave_display_label)
    if "final_fate_display" not in sub.columns:
        sub["final_fate_display"] = sub["final_fate"].astype(str).map(_wave_display_label)

    source_labels = [
        label
        for label in display_order
        if float(sub.loc[sub["zotu_class_display"] == label, "count"].sum()) > 0
    ]
    target_labels = [
        label
        for label in target_display_order
        if float(sub.loc[sub["final_fate_display"] == label, "count"].sum()) > 0
    ]
    total = float(pd.to_numeric(sub["count"], errors="coerce").sum())
    if total <= 0:
        ax.text(0.5, 0.5, "No flow data available", ha="center", va="center", fontsize=font_size)
        ax.axis("off")
        return

    top_margin = 0.94
    bottom_margin = 0.08
    usable_height = top_margin - bottom_margin
    source_gap = 0.02
    target_gap = 0.02
    left_x0, left_x1 = 0.16, 0.22
    right_x0, right_x1 = 0.78, 0.84
    ribbon_x0, ribbon_x1 = left_x1, right_x0

    def _bounds(labels: list[str], column: str, gap: float) -> dict[str, tuple[float, float]]:
        totals = {
            label: float(sub.loc[sub[column] == label, "count"].sum())
            for label in labels
        }
        gap_total = gap * max(0, len(labels) - 1)
        scale = usable_height / total if gap_total < usable_height else 0.0
        y_cursor = top_margin
        bounds: dict[str, tuple[float, float]] = {}
        for idx, label in enumerate(labels):
            h = totals[label] * scale
            y_top = y_cursor
            y_bottom = y_top - h
            bounds[label] = (y_bottom, y_top)
            y_cursor = y_bottom - (gap if idx < len(labels) - 1 else 0.0)
        return bounds

    source_bounds = _bounds(source_labels, "zotu_class_display", source_gap)
    target_bounds = _bounds(target_labels, "final_fate_display", target_gap)
    scale = usable_height / total
    flow_counts: dict[tuple[str, str], float] = {}
    for source_label in source_labels:
        for target_label in target_labels:
            hit = sub[
                (sub["zotu_class_display"] == source_label)
                & (sub["final_fate_display"] == target_label)
            ]
            if hit.empty:
                continue
            count = float(hit["count"].iloc[0])
            if count > 0:
                flow_counts[(source_label, target_label)] = count

    def _ordered_present(preferred: list[str], observed: list[str]) -> list[str]:
        used = [label for label in preferred if label in observed]
        used.extend(label for label in observed if label not in used)
        return used

    source_preference = {
        "Exact match": ["Non-bacterial", "Exact match", "Near-exact match", "No match"],
        "Near-exact match": ["Non-bacterial", "Near-exact match", "Exact match", "No match"],
        "No match": ["No match", "Non-bacterial", "Near-exact match", "Exact match"],
    }
    target_preference = {
        "Exact match": ["Near-exact match", "Exact match", "No match"],
        "Near-exact match": ["Near-exact match", "Exact match", "No match"],
        "No match": ["No match", "Near-exact match", "Exact match"],
        "Non-bacterial": ["Near-exact match", "Exact match", "No match"],
    }

    source_offsets = {label: source_bounds[label][0] for label in source_labels}
    source_segments: dict[tuple[str, str], tuple[float, float]] = {}
    for source_label in source_labels:
        observed_targets = [
            target_label
            for target_label in target_labels
            if (source_label, target_label) in flow_counts
        ]
        ordered_targets = _ordered_present(
            source_preference.get(source_label, target_labels),
            observed_targets,
        )
        for target_label in ordered_targets:
            count = flow_counts[(source_label, target_label)]
            height = count * scale
            y_bottom = source_offsets[source_label]
            y_top = y_bottom + height
            source_segments[(source_label, target_label)] = (y_bottom, y_top)
            source_offsets[source_label] += height

    target_offsets = {label: target_bounds[label][0] for label in target_labels}
    target_segments: dict[tuple[str, str], tuple[float, float]] = {}
    for target_label in target_labels:
        observed_sources = [
            source_label
            for source_label in source_labels
            if (source_label, target_label) in flow_counts
        ]
        ordered_sources = _ordered_present(
            target_preference.get(target_label, source_labels),
            observed_sources,
        )
        for source_label in ordered_sources:
            count = flow_counts[(source_label, target_label)]
            height = count * scale
            y_bottom = target_offsets[target_label]
            y_top = y_bottom + height
            target_segments[(source_label, target_label)] = (y_bottom, y_top)
            target_offsets[target_label] += height

    for (source_label, target_label), count in flow_counts.items():
        y0_bottom, y0_top = source_segments[(source_label, target_label)]
        y1_bottom, y1_top = target_segments[(source_label, target_label)]
        _add_wave_ribbon(
            ax,
            ribbon_x0,
            ribbon_x1,
            y0_bottom,
            y0_top,
            y1_bottom,
            y1_top,
            color=color_by_origin.get(source_label, "#999999"),
            alpha=0.78,
        )

    for label in source_labels:
        y_bottom, y_top = source_bounds[label]
        ax.add_patch(
            Rectangle(
                (left_x0, y_bottom),
                left_x1 - left_x0,
                y_top - y_bottom,
                facecolor="#FFFFFF",
                edgecolor=color_by_origin.get(label, "#444444"),
                linewidth=1.2,
                zorder=3,
            )
        )
        ax.text(
            left_x0 - 0.03,
            (y_bottom + y_top) / 2.0,
            label,
            ha="right",
            va="center",
            fontsize=max(8.0, font_size - 2.0),
        )
        ax.text(
            (left_x0 + left_x1) / 2.0,
            (y_bottom + y_top) / 2.0,
            f"{int(round(float(sub.loc[sub['zotu_class_display'] == label, 'count'].sum())))}",
            ha="center",
            va="center",
            fontsize=max(7.0, font_size - 3.0),
        )

    for label in target_labels:
        y_bottom, y_top = target_bounds[label]
        ax.add_patch(
            Rectangle(
                (right_x0, y_bottom),
                right_x1 - right_x0,
                y_top - y_bottom,
                facecolor="#FFFFFF",
                edgecolor="#555555",
                linewidth=1.0,
                zorder=3,
            )
        )
        ax.text(
            right_x1 + 0.03,
            (y_bottom + y_top) / 2.0,
            label,
            ha="left",
            va="center",
            fontsize=max(8.0, font_size - 2.0),
        )
        ax.text(
            (right_x0 + right_x1) / 2.0,
            (y_bottom + y_top) / 2.0,
            f"{int(round(float(sub.loc[sub['final_fate_display'] == label, 'count'].sum())))}",
            ha="center",
            va="center",
            fontsize=max(7.0, font_size - 3.0),
        )

    ax.text(
        (left_x0 + left_x1) / 2.0,
        0.985,
        "Original zOTU class",
        ha="center",
        va="top",
        fontsize=max(8.0, font_size - 2.0),
        fontweight="bold",
    )
    ax.text(
        (right_x0 + right_x1) / 2.0,
        0.985,
        "Outcome after clustering",
        ha="center",
        va="top",
        fontsize=max(8.0, font_size - 2.0),
        fontweight="bold",
    )
    ax.set_title(method, fontsize=max(11.0, font_size), pad=12.0)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.axis("off")


def _plot_wave_transition_svg(flow_summary: pd.DataFrame, out_svg: Path, font_size: float) -> None:
    color_by_origin = {
        "Exact match": "#BDBDBD",
        "Near-exact match": "#56B4E9",
        "No match": "#D55E00",
    }
    display_order = [_wave_display_label(label) for label in CLASS_ORDER]
    target_display_order = [_wave_display_label(label) for label in CLASS_WITH_UNASSIGNED_ORDER]

    def _ordered_present(preferred: list[str], observed: list[str]) -> list[str]:
        used = [label for label in preferred if label in observed]
        used.extend(label for label in observed if label not in used)
        return used

    def _prep_method(method: str) -> dict[str, object]:
        sub = flow_summary[flow_summary["transition_method"] == method].copy()
        sub = sub[pd.to_numeric(sub["count"], errors="coerce").fillna(0.0) > 0].copy()
        if "zotu_class_display" not in sub.columns:
            sub["zotu_class_display"] = sub["zotu_class"].astype(str).map(_wave_display_label)
        if "final_fate_display" not in sub.columns:
            sub["final_fate_display"] = sub["final_fate"].astype(str).map(_wave_display_label)
        source_labels = [
            label
            for label in display_order
            if float(sub.loc[sub["zotu_class_display"] == label, "count"].sum()) > 0
        ]
        target_labels = [
            label
            for label in target_display_order
            if float(sub.loc[sub["final_fate_display"] == label, "count"].sum()) > 0
        ]
        source_totals = {
            label: float(sub.loc[sub["zotu_class_display"] == label, "count"].sum())
            for label in source_labels
        }
        target_totals = {
            label: float(sub.loc[sub["final_fate_display"] == label, "count"].sum())
            for label in target_labels
        }
        flow_counts: dict[tuple[str, str], float] = {}
        for source_label in source_labels:
            for target_label in target_labels:
                hit = sub[
                    (sub["zotu_class_display"] == source_label)
                    & (sub["final_fate_display"] == target_label)
                ]
                if hit.empty:
                    continue
                count = float(hit["count"].iloc[0])
                if count > 0:
                    flow_counts[(source_label, target_label)] = count
        return {
            "source_labels": source_labels,
            "target_labels": target_labels,
            "source_totals": source_totals,
            "target_totals": target_totals,
            "flow_counts": flow_counts,
        }

    def _bounds_from_totals(
        labels: list[str],
        totals: dict[str, float],
        total_sum: float,
        top_margin: float,
        bottom_margin: float,
        gap: float,
    ) -> tuple[dict[str, tuple[float, float]], float]:
        usable_height = top_margin - bottom_margin
        gap_total = gap * max(0, len(labels) - 1)
        scale = usable_height / total_sum if total_sum > 0 and gap_total < usable_height else 0.0
        y_cursor = top_margin
        bounds: dict[str, tuple[float, float]] = {}
        for idx, label in enumerate(labels):
            height = totals.get(label, 0.0) * scale
            y_top = y_cursor
            y_bottom = y_top - height
            bounds[label] = (y_bottom, y_top)
            y_cursor = y_bottom - (gap if idx < len(labels) - 1 else 0.0)
        return bounds, scale

    def _build_edge_segments(
        source_labels: list[str],
        target_labels: list[str],
        source_bounds: dict[str, tuple[float, float]],
        target_bounds: dict[str, tuple[float, float]],
        flow_counts: dict[tuple[str, str], float],
        scale: float,
    ) -> tuple[dict[tuple[str, str], tuple[float, float]], dict[tuple[str, str], tuple[float, float]]]:
        source_preference = {
            "Exact match": ["Non-bacterial", "Exact match", "Near-exact match", "No match"],
            "Near-exact match": ["Non-bacterial", "Near-exact match", "Exact match", "No match"],
            "No match": ["No match", "Non-bacterial", "Near-exact match", "Exact match"],
        }
        target_preference = {
            "Exact match": ["Near-exact match", "Exact match", "No match"],
            "Near-exact match": ["Near-exact match", "Exact match", "No match"],
            "No match": ["No match", "Near-exact match", "Exact match"],
            "Non-bacterial": ["Near-exact match", "Exact match", "No match"],
        }

        source_offsets = {label: source_bounds[label][0] for label in source_labels}
        source_segments: dict[tuple[str, str], tuple[float, float]] = {}
        for source_label in source_labels:
            observed_targets = [
                target_label
                for target_label in target_labels
                if (source_label, target_label) in flow_counts
            ]
            ordered_targets = _ordered_present(
                source_preference.get(source_label, target_labels),
                observed_targets,
            )
            for target_label in ordered_targets:
                height = flow_counts[(source_label, target_label)] * scale
                y_bottom = source_offsets[source_label]
                y_top = y_bottom + height
                source_segments[(source_label, target_label)] = (y_bottom, y_top)
                source_offsets[source_label] += height

        target_offsets = {label: target_bounds[label][0] for label in target_labels}
        target_segments: dict[tuple[str, str], tuple[float, float]] = {}
        for target_label in target_labels:
            observed_sources = [
                source_label
                for source_label in source_labels
                if (source_label, target_label) in flow_counts
            ]
            ordered_sources = _ordered_present(
                target_preference.get(target_label, source_labels),
                observed_sources,
            )
            for source_label in ordered_sources:
                height = flow_counts[(source_label, target_label)] * scale
                y_bottom = target_offsets[target_label]
                y_top = y_bottom + height
                target_segments[(source_label, target_label)] = (y_bottom, y_top)
                target_offsets[target_label] += height
        return source_segments, target_segments

    tac = _prep_method("TAC")
    tic = _prep_method("TIC")
    center_labels = [
        label
        for label in display_order
        if float(tac["source_totals"].get(label, 0.0)) > 0 or float(tic["source_totals"].get(label, 0.0)) > 0
    ]
    center_totals = {
        label: max(float(tac["source_totals"].get(label, 0.0)), float(tic["source_totals"].get(label, 0.0)))
        for label in center_labels
    }
    total = float(sum(center_totals.values()))

    fig, ax = plt.subplots(figsize=(12.4, 6.2))
    top_margin = 0.79
    bottom_margin = 0.10
    gap = 0.02
    left_x0, left_x1 = 0.18, 0.25
    center_x0, center_x1 = 0.395, 0.535
    right_x0, right_x1 = 0.68, 0.75

    center_bounds, scale = _bounds_from_totals(
        center_labels,
        center_totals,
        total,
        top_margin=top_margin,
        bottom_margin=bottom_margin,
        gap=gap,
    )
    tac_target_bounds, _ = _bounds_from_totals(
        list(tac["target_labels"]),
        dict(tac["target_totals"]),
        total,
        top_margin=top_margin,
        bottom_margin=bottom_margin,
        gap=gap,
    )
    tic_target_bounds, _ = _bounds_from_totals(
        list(tic["target_labels"]),
        dict(tic["target_totals"]),
        total,
        top_margin=top_margin,
        bottom_margin=bottom_margin,
        gap=gap,
    )

    tac_center_segments, tac_target_segments = _build_edge_segments(
        list(tac["source_labels"]),
        list(tac["target_labels"]),
        center_bounds,
        tac_target_bounds,
        dict(tac["flow_counts"]),
        scale,
    )
    tic_center_segments, tic_target_segments = _build_edge_segments(
        list(tic["source_labels"]),
        list(tic["target_labels"]),
        center_bounds,
        tic_target_bounds,
        dict(tic["flow_counts"]),
        scale,
    )

    for (source_label, target_label), _count in dict(tac["flow_counts"]).items():
        y_left_bottom, y_left_top = tac_target_segments[(source_label, target_label)]
        y_center_bottom, y_center_top = tac_center_segments[(source_label, target_label)]
        _add_wave_ribbon(
            ax,
            left_x1,
            center_x0,
            y_left_bottom,
            y_left_top,
            y_center_bottom,
            y_center_top,
            color=color_by_origin.get(source_label, "#999999"),
            alpha=0.78,
        )
        ribbon_height = min(y_left_top - y_left_bottom, y_center_top - y_center_bottom)
        if ribbon_height >= 0.012:
            _add_wave_arrowhead(
                ax,
                x_center=((left_x1 + center_x0) / 2.0) + 0.005,
                y_center=(
                    ((y_left_bottom + y_left_top) / 2.0)
                    + ((y_center_bottom + y_center_top) / 2.0)
                )
                / 2.0,
                direction="left",
                color=color_by_origin.get(source_label, "#999999"),
                size_x=0.02,
                size_y=min(0.018, ribbon_height * 0.75),
            )

    for (source_label, target_label), _count in dict(tic["flow_counts"]).items():
        y_center_bottom, y_center_top = tic_center_segments[(source_label, target_label)]
        y_right_bottom, y_right_top = tic_target_segments[(source_label, target_label)]
        _add_wave_ribbon(
            ax,
            center_x1,
            right_x0,
            y_center_bottom,
            y_center_top,
            y_right_bottom,
            y_right_top,
            color=color_by_origin.get(source_label, "#999999"),
            alpha=0.78,
        )
        ribbon_height = min(y_center_top - y_center_bottom, y_right_top - y_right_bottom)
        if ribbon_height >= 0.012:
            _add_wave_arrowhead(
                ax,
                x_center=((center_x1 + right_x0) / 2.0) - 0.005,
                y_center=(
                    ((y_center_bottom + y_center_top) / 2.0)
                    + ((y_right_bottom + y_right_top) / 2.0)
                )
                / 2.0,
                direction="right",
                color=color_by_origin.get(source_label, "#999999"),
                size_x=0.02,
                size_y=min(0.018, ribbon_height * 0.75),
            )

    for label in list(tac["target_labels"]):
        y_bottom, y_top = tac_target_bounds[label]
        ax.add_patch(
            Rectangle(
                (left_x0, y_bottom),
                left_x1 - left_x0,
                y_top - y_bottom,
                facecolor="#FFFFFF",
                edgecolor="#555555",
                linewidth=1.0,
                zorder=3,
            )
        )
        ax.text(
            left_x0 - 0.03,
            (y_bottom + y_top) / 2.0,
            label,
            ha="right",
            va="center",
            fontsize=max(8.0, font_size - 2.0),
        )
        ax.text(
            (left_x0 + left_x1) / 2.0,
            (y_bottom + y_top) / 2.0,
            f"{int(round(float(tac['target_totals'].get(label, 0.0))))}",
            ha="center",
            va="center",
            fontsize=max(7.0, font_size - 3.0),
        )

    for label in center_labels:
        y_bottom, y_top = center_bounds[label]
        ax.add_patch(
            Rectangle(
                (center_x0, y_bottom),
                center_x1 - center_x0,
                y_top - y_bottom,
                facecolor="#FFFFFF",
                edgecolor=color_by_origin.get(label, "#444444"),
                linewidth=1.2,
                zorder=3,
            )
        )
        ax.text(
            (center_x0 + center_x1) / 2.0,
            (y_bottom + y_top) / 2.0,
            f"{label}\n{int(round(float(center_totals.get(label, 0.0))))}",
            ha="center",
            va="center",
            fontsize=max(7.0, font_size - 3.3),
            linespacing=1.05,
        )

    for label in list(tic["target_labels"]):
        y_bottom, y_top = tic_target_bounds[label]
        ax.add_patch(
            Rectangle(
                (right_x0, y_bottom),
                right_x1 - right_x0,
                y_top - y_bottom,
                facecolor="#FFFFFF",
                edgecolor="#555555",
                linewidth=1.0,
                zorder=3,
            )
        )
        ax.text(
            right_x1 + 0.03,
            (y_bottom + y_top) / 2.0,
            label,
            ha="left",
            va="center",
            fontsize=max(8.0, font_size - 2.0),
        )
        ax.text(
            (right_x0 + right_x1) / 2.0,
            (y_bottom + y_top) / 2.0,
            f"{int(round(float(tic['target_totals'].get(label, 0.0))))}",
            ha="center",
            va="center",
            fontsize=max(7.0, font_size - 3.0),
        )

    header_y = 0.845
    ax.text(
        (left_x0 + left_x1) / 2.0,
        header_y,
        "TAC outcome",
        ha="center",
        va="bottom",
        fontsize=max(8.5, font_size - 2.0),
        fontweight="bold",
    )
    ax.text(
        (center_x0 + center_x1) / 2.0,
        header_y,
        "Original zOTU class",
        ha="center",
        va="bottom",
        fontsize=max(8.5, font_size - 2.0),
        fontweight="bold",
    )
    ax.text(
        (right_x0 + right_x1) / 2.0,
        header_y,
        "TIC outcome",
        ha="center",
        va="bottom",
        fontsize=max(8.5, font_size - 2.0),
        fontweight="bold",
    )
    fig.suptitle(
        "Figure 09 test: compact wave transition view around original zOTU classes",
        fontsize=max(14.0, font_size + 2.0),
        y=0.985,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.942,
        "TIC includes the Non-bacterial outcome for zOTUs removed before final cluster membership.",
        ha="center",
        va="center",
        fontsize=max(10.0, font_size - 1.0),
    )
    fig.text(
        0.5,
        0.915,
        "Read from the middle outward; arrowheads indicate direction toward TAC and TIC outcomes.",
        ha="center",
        va="center",
        fontsize=max(9.0, font_size - 2.0),
        color="#555555",
    )

    legend_handles = [
        Rectangle((0.0, 0.0), 1.0, 1.0, facecolor="#BDBDBD", edgecolor="none", alpha=0.78, label="Exact match"),
        Rectangle((0.0, 0.0), 1.0, 1.0, facecolor="#56B4E9", edgecolor="none", alpha=0.78, label="Near-exact match"),
        Rectangle((0.0, 0.0), 1.0, 1.0, facecolor="#D55E00", edgecolor="none", alpha=0.78, label="No match"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.03),
        ncol=3,
        frameon=False,
        fontsize=max(8.0, font_size - 2.0),
        title="Ribbon color = original zOTU class",
        title_fontsize=max(8.0, font_size - 2.0),
    )
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.axis("off")
    out_svg.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=[0.0, 0.08, 1.0, 0.86])
    fig.savefig(out_svg, format="svg", bbox_inches="tight")
    plt.close(fig)


def _build_species_purity_tables(feature_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    assigned = feature_df[
        feature_df["method"].isin(["TIC", "TAC"])
        & feature_df["cluster_id"].astype(str).str.strip().ne("")
        & feature_df["classification"].astype(str).ne("Unassigned")
    ].copy()
    if assigned.empty:
        empty_cluster = pd.DataFrame(
            columns=[
                "community_id",
                "sample_id",
                "sample_key",
                "method",
                "cluster_id",
                "n_present_members",
                "total_member_abundance_tss1000",
                "dominant_species",
                "dominant_species_fraction_abundance",
                "n_species_labels",
                "unknown_member_count",
                "gini_impurity_abundance",
                "shannon_entropy_abundance",
                "is_pure_species_cluster",
            ]
        )
        empty_sample = pd.DataFrame(
            columns=[
                "community_id",
                "sample_id",
                "sample_key",
                "method",
                "total_clusters",
                "pure_clusters",
                "pure_cluster_fraction",
                "mean_dominant_species_fraction_abundance",
                "mean_gini_impurity_abundance",
                "impure_clusters",
                "mean_impure_gini_impurity_abundance",
                "mean_impure_shannon_entropy_abundance",
            ]
        )
        return empty_cluster, empty_sample

    assigned["species_label"] = np.where(
        assigned["classification"].isin(["Exact match", "Near-exact match"])
        & assigned["best_hit_species"].astype(str).str.strip().ne(""),
        assigned["best_hit_species"].astype(str),
        "Unknown / no-hit",
    )

    cluster_rows: list[dict[str, object]] = []
    for keys, grp in assigned.groupby(["community_id", "sample_id", "sample_key", "method", "cluster_id"], sort=False):
        community_id, sample_id, sample_key, method, cluster_id = keys
        abund = pd.to_numeric(grp["feature_value_tss1000"], errors="coerce").fillna(0.0)
        total_abund = float(abund.sum())
        species_abund = (
            grp.assign(_abund=abund)
            .groupby("species_label", as_index=False)["_abund"]
            .sum()
            .rename(columns={"_abund": "abundance"})
        )
        shares = (
            species_abund["abundance"].astype(float).to_numpy() / total_abund
            if total_abund > 0
            else np.zeros(len(species_abund), dtype=float)
        )
        dominant_idx = int(np.argmax(shares)) if len(shares) > 0 else 0
        dominant_species = str(species_abund["species_label"].iloc[dominant_idx]) if len(species_abund) > 0 else ""
        dominant_frac = float(shares[dominant_idx]) if len(shares) > 0 else np.nan
        gini = float(1.0 - np.sum(shares**2)) if len(shares) > 0 else np.nan
        shannon = float(-np.sum(shares * np.log2(shares))) if len(shares) > 0 else np.nan
        unknown_member_count = int((grp["species_label"] == "Unknown / no-hit").sum())
        known_labels = [label for label in species_abund["species_label"].astype(str).tolist() if label != "Unknown / no-hit"]
        is_pure = (len(set(known_labels)) == 1) and (unknown_member_count == 0) and (len(known_labels) > 0)
        cluster_rows.append(
            {
                "community_id": community_id,
                "sample_id": sample_id,
                "sample_key": sample_key,
                "method": method,
                "cluster_id": cluster_id,
                "n_present_members": int(grp.shape[0]),
                "total_member_abundance_tss1000": total_abund,
                "dominant_species": dominant_species,
                "dominant_species_fraction_abundance": dominant_frac,
                "n_species_labels": int(species_abund.shape[0]),
                "unknown_member_count": unknown_member_count,
                "gini_impurity_abundance": gini,
                "shannon_entropy_abundance": shannon,
                "is_pure_species_cluster": bool(is_pure),
            }
        )

    cluster_df = pd.DataFrame(cluster_rows)
    sample_rows: list[dict[str, object]] = []
    for keys, grp in cluster_df.groupby(["community_id", "sample_id", "sample_key", "method"], sort=False):
        impure = grp[~grp["is_pure_species_cluster"].astype(bool)].copy()
        sample_rows.append(
            {
                "community_id": keys[0],
                "sample_id": keys[1],
                "sample_key": keys[2],
                "method": keys[3],
                "total_clusters": int(grp.shape[0]),
                "pure_clusters": int(grp["is_pure_species_cluster"].astype(bool).sum()),
                "pure_cluster_fraction": float(grp["is_pure_species_cluster"].astype(bool).mean()),
                "mean_dominant_species_fraction_abundance": float(
                    pd.to_numeric(grp["dominant_species_fraction_abundance"], errors="coerce").mean()
                ),
                "mean_gini_impurity_abundance": float(
                    pd.to_numeric(grp["gini_impurity_abundance"], errors="coerce").mean()
                ),
                "impure_clusters": int(impure.shape[0]),
                "mean_impure_gini_impurity_abundance": float(
                    pd.to_numeric(impure["gini_impurity_abundance"], errors="coerce").mean()
                )
                if not impure.empty
                else np.nan,
                "mean_impure_shannon_entropy_abundance": float(
                    pd.to_numeric(impure["shannon_entropy_abundance"], errors="coerce").mean()
                )
                if not impure.empty
                else np.nan,
            }
        )
    sample_df = pd.DataFrame(sample_rows)
    return cluster_df, sample_df


def _build_match_status_species_purity_tables(feature_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    assigned = feature_df[
        feature_df["method"].isin(["TIC", "TAC"])
        & feature_df["cluster_id"].astype(str).str.strip().ne("")
        & feature_df["classification"].astype(str).ne("Unassigned")
    ].copy()
    cluster_columns = [
        "community_id",
        "sample_id",
        "sample_key",
        "method",
        "cluster_id",
        "n_present_members",
        "total_member_abundance_tss1000",
        "original_community_abundance_tss1000",
        "match_statuses",
        "n_match_statuses",
        "is_match_homogeneous",
        "is_match_heterogeneous",
        "species_labels",
        "known_species_labels",
        "n_species_labels",
        "n_known_species_labels",
        "unknown_member_count",
        "is_species_homogeneous_strict",
        "is_species_homogeneous_known_only",
        "category",
    ]
    summary_columns = ["method", "category", "cluster_count", "total_multi_zotu_clusters", "fraction"]
    het_columns = [
        "method",
        "match_heterogeneous_clusters",
        "match_heterogeneous_species_homogeneous_strict",
        "fraction_strict",
        "match_heterogeneous_species_homogeneous_known_only",
        "fraction_known_only",
    ]
    if assigned.empty:
        return (
            pd.DataFrame(columns=cluster_columns),
            pd.DataFrame(columns=summary_columns),
            pd.DataFrame(columns=het_columns),
        )

    assigned["species_label"] = np.where(
        assigned["best_hit_species"].astype(str).str.strip().ne(""),
        assigned["best_hit_species"].astype(str),
        "Unknown / no best hit",
    )
    assigned["match_status"] = assigned["classification"].astype(str)
    original_abundance_by_sample = (
        feature_df[feature_df["method"].isin(["TIC", "TAC"])]
        .groupby(["community_id", "sample_id", "sample_key", "method"], sort=False)[
            "feature_value_tss1000"
        ]
        .sum()
    )

    cluster_rows: list[dict[str, object]] = []
    for keys, grp in assigned.groupby(["community_id", "sample_id", "sample_key", "method", "cluster_id"], sort=False):
        community_id, sample_id, sample_key, method, cluster_id = keys
        abund = pd.to_numeric(grp["feature_value_tss1000"], errors="coerce").fillna(0.0)
        match_statuses = sorted(set(grp["match_status"].astype(str).tolist()))
        species_labels = sorted(set(grp["species_label"].astype(str).tolist()))
        known_species_labels = [label for label in species_labels if label != "Unknown / no best hit"]
        unknown_member_count = int((grp["species_label"].astype(str) == "Unknown / no best hit").sum())
        n_match_statuses = len(match_statuses)
        n_species_labels = len(species_labels)
        n_known_species_labels = len(set(known_species_labels))
        is_match_homogeneous = n_match_statuses == 1
        is_match_heterogeneous = n_match_statuses > 1
        is_species_homogeneous_strict = (
            n_species_labels == 1
            and unknown_member_count == 0
            and n_known_species_labels == 1
        )
        is_species_homogeneous_known_only = n_known_species_labels == 1
        match_part = "Match-homogeneous" if is_match_homogeneous else "Match-heterogeneous"
        species_part = (
            "Species-homogeneous"
            if is_species_homogeneous_strict
            else "Species-heterogeneous/unresolved"
        )
        original_abundance = float(
            original_abundance_by_sample.get((community_id, sample_id, sample_key, method), np.nan)
        )
        cluster_rows.append(
            {
                "community_id": community_id,
                "sample_id": sample_id,
                "sample_key": sample_key,
                "method": method,
                "cluster_id": cluster_id,
                "n_present_members": int(grp.shape[0]),
                "total_member_abundance_tss1000": float(abund.sum()),
                "original_community_abundance_tss1000": original_abundance,
                "match_statuses": ";".join(match_statuses),
                "n_match_statuses": int(n_match_statuses),
                "is_match_homogeneous": bool(is_match_homogeneous),
                "is_match_heterogeneous": bool(is_match_heterogeneous),
                "species_labels": ";".join(species_labels),
                "known_species_labels": ";".join(known_species_labels),
                "n_species_labels": int(n_species_labels),
                "n_known_species_labels": int(n_known_species_labels),
                "unknown_member_count": int(unknown_member_count),
                "is_species_homogeneous_strict": bool(is_species_homogeneous_strict),
                "is_species_homogeneous_known_only": bool(is_species_homogeneous_known_only),
                "category": f"{match_part} | {species_part}",
            }
        )

    cluster_df = pd.DataFrame(cluster_rows, columns=cluster_columns)
    multi_df = cluster_df[pd.to_numeric(cluster_df["n_present_members"], errors="coerce").fillna(0) > 1].copy()
    if multi_df.empty:
        return cluster_df, pd.DataFrame(columns=summary_columns), pd.DataFrame(columns=het_columns)

    category_order = [
        "Match-heterogeneous | Species-homogeneous",
        "Match-heterogeneous | Species-heterogeneous/unresolved",
        "Match-homogeneous | Species-homogeneous",
        "Match-homogeneous | Species-heterogeneous/unresolved",
    ]
    summary = (
        multi_df.groupby(["method", "category"], as_index=False)
        .size()
        .rename(columns={"size": "cluster_count"})
    )
    total = (
        multi_df.groupby("method", as_index=False)
        .size()
        .rename(columns={"size": "total_multi_zotu_clusters"})
    )
    records: list[dict[str, object]] = []
    for method in ["TIC", "TAC"]:
        method_total_hit = total[total["method"] == method]
        method_total = int(method_total_hit["total_multi_zotu_clusters"].iloc[0]) if not method_total_hit.empty else 0
        for category in category_order:
            hit = summary[(summary["method"] == method) & (summary["category"] == category)]
            count = int(hit["cluster_count"].iloc[0]) if not hit.empty else 0
            records.append(
                {
                    "method": method,
                    "category": category,
                    "cluster_count": count,
                    "total_multi_zotu_clusters": method_total,
                    "fraction": (float(count) / float(method_total)) if method_total > 0 else np.nan,
                }
            )
    summary_df = pd.DataFrame(records, columns=summary_columns)

    het_records: list[dict[str, object]] = []
    for method in ["TIC", "TAC"]:
        sub = multi_df[multi_df["method"] == method].copy()
        hetero = sub[sub["is_match_heterogeneous"].astype(bool)].copy()
        n_hetero = int(hetero.shape[0])
        n_strict = int(hetero["is_species_homogeneous_strict"].astype(bool).sum()) if n_hetero > 0 else 0
        n_known = int(hetero["is_species_homogeneous_known_only"].astype(bool).sum()) if n_hetero > 0 else 0
        het_records.append(
            {
                "method": method,
                "match_heterogeneous_clusters": n_hetero,
                "match_heterogeneous_species_homogeneous_strict": n_strict,
                "fraction_strict": (float(n_strict) / float(n_hetero)) if n_hetero > 0 else np.nan,
                "match_heterogeneous_species_homogeneous_known_only": n_known,
                "fraction_known_only": (float(n_known) / float(n_hetero)) if n_hetero > 0 else np.nan,
            }
        )
    hetero_summary_df = pd.DataFrame(het_records, columns=het_columns)
    return cluster_df, summary_df, hetero_summary_df


def _plot_match_status_species_purity_figure(
    summary_df: pd.DataFrame,
    hetero_summary_df: pd.DataFrame,
    out_svg: Path,
    font_size: float,
) -> None:
    methods = ["TIC", "TAC"]
    category_order = [
        "Match-heterogeneous | Species-homogeneous",
        "Match-heterogeneous | Species-heterogeneous/unresolved",
        "Match-homogeneous | Species-homogeneous",
        "Match-homogeneous | Species-heterogeneous/unresolved",
    ]
    category_labels = {
        "Match-heterogeneous | Species-homogeneous": "Match-heterogeneous\nSpecies-homogeneous",
        "Match-heterogeneous | Species-heterogeneous/unresolved": "Match-heterogeneous\nSpecies-heterogeneous/unresolved",
        "Match-homogeneous | Species-homogeneous": "Match-homogeneous\nSpecies-homogeneous",
        "Match-homogeneous | Species-heterogeneous/unresolved": "Match-homogeneous\nSpecies-heterogeneous/unresolved",
    }
    colors = {
        "Match-heterogeneous | Species-homogeneous": "#009E73",
        "Match-heterogeneous | Species-heterogeneous/unresolved": "#D55E00",
        "Match-homogeneous | Species-homogeneous": "#56B4E9",
        "Match-homogeneous | Species-heterogeneous/unresolved": "#E69F00",
    }

    fig, ax = plt.subplots(figsize=(9.4, 5.8))
    xs = np.arange(len(methods), dtype=float)
    bottom = np.zeros(len(methods), dtype=float)
    for category in category_order:
        vals = []
        counts = []
        totals = []
        for method in methods:
            hit = summary_df[(summary_df["method"] == method) & (summary_df["category"] == category)]
            vals.append(float(hit["fraction"].iloc[0]) if not hit.empty and pd.notna(hit["fraction"].iloc[0]) else 0.0)
            counts.append(int(hit["cluster_count"].iloc[0]) if not hit.empty else 0)
            totals.append(int(hit["total_multi_zotu_clusters"].iloc[0]) if not hit.empty else 0)
        vals_arr = np.asarray(vals, dtype=float)
        ax.bar(
            xs,
            vals_arr,
            bottom=bottom,
            width=0.62,
            color=colors[category],
            edgecolor="#333333",
            linewidth=0.7,
            label=category_labels[category],
        )
        for x, y0, h, count in zip(xs, bottom, vals_arr, counts):
            if h <= 0:
                continue
            ax.text(
                x,
                y0 + (h / 2.0),
                f"{count}\n{h * 100:.1f}%",
                ha="center",
                va="center",
                fontsize=max(6.4, font_size - 4.0),
                color="white" if h >= 0.16 else "#1F1F1F",
                linespacing=0.9,
            )
        bottom += vals_arr

    for x, method in zip(xs, methods):
        total_hit = summary_df[summary_df["method"] == method]
        total = int(total_hit["total_multi_zotu_clusters"].max()) if not total_hit.empty else 0
        het_hit = hetero_summary_df[hetero_summary_df["method"] == method]
        if not het_hit.empty:
            n_hetero = int(het_hit["match_heterogeneous_clusters"].iloc[0])
            n_strict = int(het_hit["match_heterogeneous_species_homogeneous_strict"].iloc[0])
            frac = float(het_hit["fraction_strict"].iloc[0]) if pd.notna(het_hit["fraction_strict"].iloc[0]) else np.nan
            txt = f"multi-zOTU n={total}\nheterogeneous: {n_strict}/{n_hetero}"
            if np.isfinite(frac):
                txt += f" ({frac * 100:.1f}%)"
        else:
            txt = f"multi-zOTU n={total}\nheterogeneous: NA"
        ax.text(
            x,
            1.035,
            txt,
            ha="center",
            va="bottom",
            fontsize=max(7.2, font_size - 3.2),
            color="#222222",
        )

    ax.set_xticks(xs)
    ax.set_xticklabels(methods, fontsize=max(8.0, font_size - 1.0))
    ax.set_ylim(0.0, 1.18)
    ax.set_ylabel("Fraction of multi-zOTU clusters", fontsize=max(8.0, font_size - 1.0))
    ax.set_xlabel("Method", fontsize=max(8.0, font_size - 1.0))
    ax.set_title(
        "Match-status heterogeneity versus species homogeneity",
        fontsize=font_size + 1.0,
    )
    ax.grid(axis="y", color="#E6E6E6", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        frameon=False,
        fontsize=max(7.0, font_size - 3.0),
        title="Cluster category",
        title_fontsize=max(7.4, font_size - 2.6),
    )
    fig.suptitle(
        "Figure 09 test: species identity of match-status heterogeneous TIC/TAC clusters",
        fontsize=font_size + 2.0,
        y=1.02,
    )
    fig.tight_layout(rect=[0.0, 0.0, 0.78, 0.96])
    out_svg.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_svg, format="svg", bbox_inches="tight")
    plt.close(fig)


def _build_exact_near_species_purity_tables(cluster_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    mix_specs = [
        ("pure_exact", "Pure exact", ("Exact match",), True),
        ("pure_near", "Pure near-exact", ("Near-exact match",), True),
        ("pure_no_match", "Pure no match", ("No hit",), True),
        ("exact_near", "Exact + near-exact", ("Exact match", "Near-exact match"), False),
    ]
    cluster_columns = list(cluster_df.columns) + [
        "mix_type",
        "mix_label",
        "has_exact",
        "has_near_exact",
        "has_no_match",
        "taxonomy_status",
    ] if not cluster_df.empty else [
        "mix_type",
        "mix_label",
        "community_id",
        "sample_id",
        "sample_key",
        "method",
        "cluster_id",
        "n_present_members",
        "match_statuses",
        "species_labels",
        "known_species_labels",
        "taxonomy_status",
    ]
    summary_columns = [
        "mix_type",
        "mix_label",
        "method",
        "taxonomy_status",
        "cluster_count",
        "total_clusters_in_group",
        "fraction",
        "abundance_tss1000",
        "total_group_abundance_tss1000",
        "group_abundance_fraction",
        "total_original_community_abundance_tss1000",
        "original_community_abundance_fraction",
        "present_sample_count",
        "present_mean_abundance_fraction",
    ]
    method_columns = [
        "mix_type",
        "mix_label",
        "method",
        "clusters_in_group",
        "single_species_clusters",
        "multi_species_or_unresolved_clusters",
        "single_species_fraction",
        "total_group_abundance_tss1000",
        "single_species_abundance_tss1000",
        "multi_species_or_unresolved_abundance_tss1000",
        "single_species_group_abundance_fraction",
        "total_original_community_abundance_tss1000",
        "single_species_original_community_abundance_fraction",
        "group_original_community_abundance_fraction",
        "single_species_present_sample_count",
        "single_species_present_mean_abundance_fraction",
    ]
    if cluster_df.empty:
        return (
            pd.DataFrame(columns=cluster_columns),
            pd.DataFrame(columns=summary_columns),
            pd.DataFrame(columns=method_columns),
        )
    cluster_df = cluster_df[
        pd.to_numeric(cluster_df["n_present_members"], errors="coerce").fillna(0) > 1
    ].copy()
    if cluster_df.empty:
        return (
            pd.DataFrame(columns=cluster_columns),
            pd.DataFrame(columns=summary_columns),
            pd.DataFrame(columns=method_columns),
        )

    rows: list[dict[str, object]] = []
    for rec in cluster_df.to_dict(orient="records"):
        statuses = [s for s in str(rec.get("match_statuses", "")).split(";") if s]
        has_exact = "Exact match" in statuses
        has_near = "Near-exact match" in statuses
        has_no = "No hit" in statuses
        status_flags = {"Exact match": has_exact, "Near-exact match": has_near, "No hit": has_no}
        for mix_type, mix_label, required, require_exact_status_set in mix_specs:
            if require_exact_status_set:
                if set(statuses) != set(required):
                    continue
            elif not all(status_flags.get(status, False) for status in required):
                continue
            taxonomy_status = (
                "Single species"
                if bool(rec.get("is_species_homogeneous_strict", False))
                else "Multi-species or unresolved"
            )
            out = dict(rec)
            out.update(
                {
                    "mix_type": mix_type,
                    "mix_label": mix_label,
                    "has_exact": bool(has_exact),
                    "has_near_exact": bool(has_near),
                    "has_no_match": bool(has_no),
                    "taxonomy_status": taxonomy_status,
                }
            )
            rows.append(out)

    mixed_df = pd.DataFrame(rows)
    if mixed_df.empty:
        return (
            pd.DataFrame(columns=cluster_columns),
            pd.DataFrame(columns=summary_columns),
            pd.DataFrame(columns=method_columns),
        )

    status_order = ["Single species", "Multi-species or unresolved"]
    summary_raw = (
        mixed_df.groupby(["mix_type", "mix_label", "method", "taxonomy_status"], as_index=False)
        .size()
        .rename(columns={"size": "cluster_count"})
    )
    abundance_raw = (
        mixed_df.groupby(["mix_type", "mix_label", "method", "taxonomy_status"], as_index=False)[
            "total_member_abundance_tss1000"
        ]
        .sum()
        .rename(columns={"total_member_abundance_tss1000": "abundance_tss1000"})
    )
    present_sample = (
        mixed_df.groupby(
            [
                "mix_type",
                "mix_label",
                "method",
                "taxonomy_status",
                "community_id",
                "sample_id",
                "sample_key",
            ],
            as_index=False,
        )
        .agg(
            abundance_tss1000=("total_member_abundance_tss1000", "sum"),
            original_community_abundance_tss1000=("original_community_abundance_tss1000", "first"),
        )
    )
    present_sample["sample_abundance_fraction"] = np.where(
        pd.to_numeric(present_sample["original_community_abundance_tss1000"], errors="coerce") > 0,
        pd.to_numeric(present_sample["abundance_tss1000"], errors="coerce")
        / pd.to_numeric(present_sample["original_community_abundance_tss1000"], errors="coerce"),
        np.nan,
    )
    present_raw = (
        present_sample.groupby(["mix_type", "mix_label", "method", "taxonomy_status"], as_index=False)
        .agg(
            present_sample_count=("sample_abundance_fraction", "size"),
            present_mean_abundance_fraction=("sample_abundance_fraction", "mean"),
        )
    )
    totals = (
        mixed_df.groupby(["mix_type", "mix_label", "method"], as_index=False)
        .size()
        .rename(columns={"size": "total_clusters_in_group"})
    )
    abundance_totals = (
        mixed_df.groupby(["mix_type", "mix_label", "method"], as_index=False)[
            "total_member_abundance_tss1000"
        ]
        .sum()
        .rename(columns={"total_member_abundance_tss1000": "total_group_abundance_tss1000"})
    )
    original_totals_base = cluster_df[
        [
            "method",
            "community_id",
            "sample_id",
            "sample_key",
            "original_community_abundance_tss1000",
        ]
    ].drop_duplicates(subset=["method", "community_id", "sample_id", "sample_key"])
    original_totals_by_method = (
        original_totals_base.groupby("method", as_index=False)[
            "original_community_abundance_tss1000"
        ]
        .sum()
        .rename(
            columns={
                "original_community_abundance_tss1000": "total_original_community_abundance_tss1000"
            }
        )
    )
    summary_records: list[dict[str, object]] = []
    method_records: list[dict[str, object]] = []
    for mix_type, mix_label, _required, _require_exact_status_set in mix_specs:
        for method in ["TIC", "TAC"]:
            total_hit = totals[(totals["mix_type"] == mix_type) & (totals["method"] == method)]
            total = int(total_hit["total_clusters_in_group"].iloc[0]) if not total_hit.empty else 0
            abundance_total_hit = abundance_totals[
                (abundance_totals["mix_type"] == mix_type) & (abundance_totals["method"] == method)
            ]
            abundance_total = (
                float(abundance_total_hit["total_group_abundance_tss1000"].iloc[0])
                if not abundance_total_hit.empty
                else 0.0
            )
            original_total_hit = original_totals_by_method[original_totals_by_method["method"] == method]
            original_total = (
                float(original_total_hit["total_original_community_abundance_tss1000"].iloc[0])
                if not original_total_hit.empty
                else 0.0
            )
            single_count = 0
            multi_count = 0
            single_abundance = 0.0
            multi_abundance = 0.0
            single_present_sample_count = 0
            single_present_mean = np.nan
            for status in status_order:
                hit = summary_raw[
                    (summary_raw["mix_type"] == mix_type)
                    & (summary_raw["method"] == method)
                    & (summary_raw["taxonomy_status"] == status)
                ]
                count = int(hit["cluster_count"].iloc[0]) if not hit.empty else 0
                abundance_hit = abundance_raw[
                    (abundance_raw["mix_type"] == mix_type)
                    & (abundance_raw["method"] == method)
                    & (abundance_raw["taxonomy_status"] == status)
                ]
                abundance = (
                    float(abundance_hit["abundance_tss1000"].iloc[0])
                    if not abundance_hit.empty
                    else 0.0
                )
                present_hit = present_raw[
                    (present_raw["mix_type"] == mix_type)
                    & (present_raw["method"] == method)
                    & (present_raw["taxonomy_status"] == status)
                ]
                present_sample_count = (
                    int(present_hit["present_sample_count"].iloc[0])
                    if not present_hit.empty
                    else 0
                )
                present_mean = (
                    float(present_hit["present_mean_abundance_fraction"].iloc[0])
                    if not present_hit.empty
                    else np.nan
                )
                if status == "Single species":
                    single_count = count
                    single_abundance = abundance
                    single_present_sample_count = present_sample_count
                    single_present_mean = present_mean
                else:
                    multi_count = count
                    multi_abundance = abundance
                summary_records.append(
                    {
                        "mix_type": mix_type,
                        "mix_label": mix_label,
                        "method": method,
                        "taxonomy_status": status,
                        "cluster_count": count,
                        "total_clusters_in_group": total,
                        "fraction": (float(count) / float(total)) if total > 0 else np.nan,
                        "abundance_tss1000": abundance,
                        "total_group_abundance_tss1000": abundance_total,
                        "group_abundance_fraction": (
                            float(abundance) / float(abundance_total)
                            if abundance_total > 0
                            else np.nan
                        ),
                        "total_original_community_abundance_tss1000": original_total,
                        "original_community_abundance_fraction": (
                            float(abundance) / float(original_total)
                            if original_total > 0
                            else np.nan
                        ),
                        "present_sample_count": present_sample_count,
                        "present_mean_abundance_fraction": present_mean,
                    }
                )
            method_records.append(
                {
                    "mix_type": mix_type,
                    "mix_label": mix_label,
                    "method": method,
                    "clusters_in_group": total,
                    "single_species_clusters": single_count,
                    "multi_species_or_unresolved_clusters": multi_count,
                    "single_species_fraction": (float(single_count) / float(total)) if total > 0 else np.nan,
                    "total_group_abundance_tss1000": abundance_total,
                    "single_species_abundance_tss1000": single_abundance,
                    "multi_species_or_unresolved_abundance_tss1000": multi_abundance,
                    "single_species_group_abundance_fraction": (
                        float(single_abundance) / float(abundance_total)
                        if abundance_total > 0
                        else np.nan
                    ),
                    "total_original_community_abundance_tss1000": original_total,
                    "single_species_original_community_abundance_fraction": (
                        float(single_abundance) / float(original_total)
                        if original_total > 0
                        else np.nan
                    ),
                    "group_original_community_abundance_fraction": (
                        float(abundance_total) / float(original_total)
                        if original_total > 0
                        else np.nan
                    ),
                    "single_species_present_sample_count": single_present_sample_count,
                    "single_species_present_mean_abundance_fraction": single_present_mean,
                }
            )

    return (
        mixed_df,
        pd.DataFrame(summary_records, columns=summary_columns),
        pd.DataFrame(method_records, columns=method_columns),
    )


def _plot_exact_near_species_purity_figure(
    summary_df: pd.DataFrame,
    method_summary_df: pd.DataFrame,
    out_svg: Path,
    font_size: float,
) -> None:
    methods = ["TIC", "TAC"]
    status_order = ["Single species", "Multi-species or unresolved"]
    colors = {
        "Single species": "#009E73",
        "Multi-species or unresolved": "#D55E00",
    }
    mix_order = [
        ("pure_exact", "Pure exact"),
        ("pure_near", "Pure near-exact"),
        ("pure_no_match", "Pure no match"),
        ("exact_near", "Exact + near-exact"),
    ]
    fig, axes_grid = plt.subplots(2, 2, figsize=(13.2, 10.2), sharey=False)
    axes = axes_grid.ravel()
    legend_handles = None
    effect_axes = []
    for ax, (mix_type, mix_label) in zip(axes, mix_order):
        ax_effect = ax.twinx()
        effect_axes.append(ax_effect)
        xs = np.arange(len(methods), dtype=float)
        count_xs = xs - 0.17
        effect_xs = xs + 0.17
        count_bottom = np.zeros(len(methods), dtype=float)
        effect_bottom = np.zeros(len(methods), dtype=float)
        for status in status_order:
            count_vals = []
            effect_vals = []
            counts = []
            effect_sample_counts = []
            for method in methods:
                hit = summary_df[
                    (summary_df["mix_type"] == mix_type)
                    & (summary_df["method"] == method)
                    & (summary_df["taxonomy_status"] == status)
                ]
                count_vals.append(float(hit["cluster_count"].iloc[0]) if not hit.empty else 0.0)
                effect_vals.append(
                    float(hit["present_mean_abundance_fraction"].iloc[0])
                    if not hit.empty and pd.notna(hit["present_mean_abundance_fraction"].iloc[0])
                    else 0.0
                )
                counts.append(int(hit["cluster_count"].iloc[0]) if not hit.empty else 0)
                effect_sample_counts.append(int(hit["present_sample_count"].iloc[0]) if not hit.empty else 0)
            count_vals_arr = np.asarray(count_vals, dtype=float)
            effect_vals_arr = np.asarray(effect_vals, dtype=float)
            bars = ax.bar(
                count_xs,
                count_vals_arr,
                bottom=count_bottom,
                width=0.30,
                color=colors[status],
                edgecolor="#333333",
                linewidth=0.7,
                label=status,
            )
            ax_effect.bar(
                effect_xs,
                effect_vals_arr,
                bottom=effect_bottom,
                width=0.30,
                color=colors[status],
                edgecolor="#333333",
                linewidth=0.7,
                hatch="...",
                alpha=0.78,
            )
            if legend_handles is None:
                legend_handles = bars
            for x, y0, h, count in zip(count_xs, count_bottom, count_vals_arr, counts):
                if h <= 0:
                    continue
                ax.text(
                    x,
                    y0 + (h / 2.0),
                    f"{count}",
                    ha="center",
                    va="center",
                    fontsize=max(6.6, font_size - 4.2),
                    color="white" if h >= 0.18 else "#1F1F1F",
                    linespacing=0.9,
                )
            for x, y0, h, sample_count in zip(effect_xs, effect_bottom, effect_vals_arr, effect_sample_counts):
                if h <= 0:
                    continue
                ax_effect.text(
                    x,
                    y0 + (h / 2.0),
                    f"{_format_pct(h)}\nn={sample_count}",
                    ha="center",
                    va="center",
                    fontsize=max(5.6, font_size - 5.2),
                    color="white" if h >= 0.18 else "#1F1F1F",
                    linespacing=0.9,
                )
            count_bottom += count_vals_arr
            effect_bottom += effect_vals_arr

        for x, method in zip(xs, methods):
            hit = method_summary_df[
                (method_summary_df["mix_type"] == mix_type)
                & (method_summary_df["method"] == method)
            ]
            if hit.empty:
                txt = "n=0"
            else:
                row = hit.iloc[0]
                total = int(row["clusters_in_group"])
                single = int(row["single_species_clusters"])
                frac = float(row["single_species_fraction"]) if pd.notna(row["single_species_fraction"]) else np.nan
                txt = f"n={total}\n{single}/{total}"
                if np.isfinite(frac):
                    txt += f" ({_format_pct(frac)})"
            ax.text(
                x,
                1.02,
                txt,
                ha="center",
                va="bottom",
                fontsize=max(5.8, font_size - 4.8),
                color="#222222",
                linespacing=0.92,
                transform=ax.get_xaxis_transform(),
            )

        ax.set_xticks(xs)
        ax.set_xticklabels(methods, fontsize=max(8.0, font_size - 1.0))
        for x in xs:
            ax.text(
                x - 0.17,
                -0.075,
                "Count",
                ha="center",
                va="top",
                fontsize=max(6.4, font_size - 4.0),
                transform=ax.get_xaxis_transform(),
                color="#333333",
            )
            ax.text(
                x + 0.17,
                -0.075,
                "Abundance\neffect",
                ha="center",
                va="top",
                fontsize=max(6.4, font_size - 4.0),
                transform=ax.get_xaxis_transform(),
                color="#333333",
            )
        count_top = max(1.0, float(np.max(count_bottom)) if count_bottom.size else 1.0)
        ax.set_ylim(0.0, count_top * 1.25)
        ax_effect.set_ylim(0.0, 1.0)
        ax_effect.yaxis.set_major_formatter(FuncFormatter(lambda val, _pos: _format_pct(val)))
        ax_effect.tick_params(axis="y", labelsize=max(7.0, font_size - 3.2), colors="#555555")
        ax_effect.spines["right"].set_color("#555555")
        ax.set_xlabel("Method", fontsize=max(8.0, font_size - 1.0))
        ax.set_title(mix_label, fontsize=font_size)
        ax.grid(axis="y", color="#E6E6E6", linewidth=0.8)
        ax.set_axisbelow(True)
    axes[0].set_ylabel("Cluster count", fontsize=max(8.0, font_size - 1.0))
    axes[2].set_ylabel("Cluster count", fontsize=max(8.0, font_size - 1.0))
    effect_axes[1].set_ylabel("Abundance effect", fontsize=max(8.0, font_size - 1.0), color="#555555")
    effect_axes[3].set_ylabel("Abundance effect", fontsize=max(8.0, font_size - 1.0), color="#555555")
    axes[2].set_xlabel("Method", fontsize=max(8.0, font_size - 1.0))
    axes[3].set_xlabel("Method", fontsize=max(8.0, font_size - 1.0))
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.36, 0.01),
        ncol=2,
        frameon=False,
        fontsize=max(7.2, font_size - 2.8),
    )
    type_handles = [
        Rectangle((0, 0), 1, 1, facecolor="#F2F2F2", edgecolor="#333333", linewidth=0.7),
        Rectangle((0, 0), 1, 1, facecolor="#F2F2F2", edgecolor="#333333", linewidth=0.7, hatch="..."),
    ]
    fig.legend(
        type_handles,
        ["Absolute count", "Abundance effect"],
        loc="lower center",
        bbox_to_anchor=(0.70, 0.01),
        ncol=2,
        frameon=False,
        fontsize=max(7.2, font_size - 2.8),
    )
    fig.suptitle(
        "Figure 09 test: taxonomy homogeneity across cluster match-status groups",
        fontsize=font_size + 2.0,
        y=0.985,
    )
    fig.tight_layout(rect=[0.0, 0.08, 1.0, 0.955], h_pad=2.2, w_pad=1.8)
    out_svg.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_svg, format="svg", bbox_inches="tight")
    plt.close(fig)


def _build_species_purity_stack_summary(cluster_df: pd.DataFrame) -> pd.DataFrame:
    if cluster_df.empty:
        return pd.DataFrame(
            columns=[
                "community_id",
                "sample_id",
                "sample_key",
                "method",
                "purity_category",
                "cluster_count",
                "total_clusters",
                "fraction",
            ]
        )

    cat = np.where(
        cluster_df["is_pure_species_cluster"].astype(bool),
        "Pure species",
        np.where(
            pd.to_numeric(cluster_df["unknown_member_count"], errors="coerce").fillna(0.0) > 0,
            "Impure + unknown/no-hit",
            "Impure multi-species",
        ),
    )
    base = cluster_df.copy()
    base["purity_category"] = cat
    summary = (
        base.groupby(["community_id", "sample_id", "sample_key", "method", "purity_category"], as_index=False)
        .size()
        .rename(columns={"size": "cluster_count"})
    )
    totals = (
        base.groupby(["community_id", "sample_id", "sample_key", "method"], as_index=False)
        .size()
        .rename(columns={"size": "total_clusters"})
    )
    summary = summary.merge(totals, on=["community_id", "sample_id", "sample_key", "method"], how="left")
    summary["fraction"] = np.where(
        pd.to_numeric(summary["total_clusters"], errors="coerce") > 0,
        pd.to_numeric(summary["cluster_count"], errors="coerce")
        / pd.to_numeric(summary["total_clusters"], errors="coerce"),
        np.nan,
    )
    return summary


def _plot_species_purity_figure(
    cluster_df: pd.DataFrame,
    sample_df: pd.DataFrame,
    out_svg: Path,
    font_size: float,
) -> None:
    stack_df = _build_species_purity_stack_summary(cluster_df)
    order = ["TIC", "TAC"]
    category_order = ["Pure species", "Impure multi-species", "Impure + unknown/no-hit"]
    colors = {
        "Pure species": "#009E73",
        "Impure multi-species": "#E69F00",
        "Impure + unknown/no-hit": "#CC79A7",
    }
    sample_order = sorted(stack_df["sample_key"].astype(str).unique().tolist()) if not stack_df.empty else []
    fig, axes = plt.subplots(1, 2, figsize=(18.0, 6.8), sharey=True, constrained_layout=True)

    for ax, method in zip(axes, order):
        sub = stack_df[stack_df["method"] == method].copy()
        positions = np.arange(len(sample_order), dtype=float)
        bottom = np.zeros(len(sample_order), dtype=float)
        for category in category_order:
            vals = []
            counts = []
            for sample_key in sample_order:
                hit = sub[(sub["sample_key"] == sample_key) & (sub["purity_category"] == category)]
                vals.append(float(hit["fraction"].iloc[0]) if not hit.empty else 0.0)
                counts.append(int(hit["cluster_count"].iloc[0]) if not hit.empty else 0)
            vals_arr = np.asarray(vals, dtype=float)
            ax.bar(
                positions,
                vals_arr,
                bottom=bottom,
                width=0.82,
                color=colors[category],
                edgecolor="#333333",
                linewidth=0.6,
                label=category,
            )
            for x, y0, frac, count in zip(positions, bottom, vals_arr, counts):
                if frac < 0.09:
                    continue
                ax.text(
                    x,
                    y0 + (frac / 2.0),
                    f"{count}",
                    ha="center",
                    va="center",
                    fontsize=max(6.0, font_size - 4.0),
                    color="white" if category != "Impure multi-species" else "#1F1F1F",
                )
            bottom += vals_arr

        pure_mean = pd.to_numeric(sample_df.loc[sample_df["method"] == method, "pure_cluster_fraction"], errors="coerce").mean()
        ax.set_title(
            f"{method}\nmean pure fraction={pure_mean:.2f}" if np.isfinite(pure_mean) else method,
            fontsize=font_size,
        )
        ax.set_xticks(positions)
        ax.set_xticklabels(sample_order, rotation=90, fontsize=max(6.8, font_size - 4.0))
        ax.set_ylim(0.0, 1.02)
        ax.grid(axis="y", color="#E6E6E6", linewidth=0.8)
        ax.set_axisbelow(True)
        ax.set_xlabel("Sample", fontsize=max(8.0, font_size - 1.0))

    axes[0].set_ylabel("Fraction of present clusters", fontsize=max(8.0, font_size - 1.0))
    handles = [Rectangle((0, 0), 1, 1, facecolor=colors[c], edgecolor="#333333", linewidth=0.6, label=c) for c in category_order]
    fig.legend(
        handles=handles,
        loc="upper center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 1.02),
        fontsize=max(7.0, font_size - 2.0),
    )
    fig.suptitle(
        "Figure 09 test: species purity of TIC/TAC clusters",
        fontsize=font_size + 2.0,
        y=1.08,
    )
    fig.savefig(out_svg, format="svg", bbox_inches="tight")
    plt.close(fig)


def _build_species_purity_venn_summary(cluster_df: pd.DataFrame) -> pd.DataFrame:
    if cluster_df.empty:
        return pd.DataFrame(
            columns=[
                "method",
                "total_clusters",
                "multi_zotu_clusters",
                "single_species_clusters",
                "pure_multi_zotu_clusters",
                "impure_multi_zotu_clusters",
                "pure_single_zotu_clusters",
                "other_clusters",
                "pure_fraction_within_multi_zotu",
                "impure_fraction_within_multi_zotu",
            ]
        )

    base = cluster_df.copy()
    base["is_multi_zotu_cluster"] = (
        pd.to_numeric(base["n_present_members"], errors="coerce").fillna(0.0) > 1
    )
    base["is_single_species_cluster"] = base["is_pure_species_cluster"].astype(bool)

    records: list[dict[str, object]] = []
    for method in ["TIC", "TAC"]:
        sub = base[base["method"] == method].copy()
        total = int(sub.shape[0])
        multi_zotu = sub["is_multi_zotu_cluster"].astype(bool)
        single_species = sub["is_single_species_cluster"].astype(bool)
        overlap = int((multi_zotu & single_species).sum())
        left_only = int((multi_zotu & ~single_species).sum())
        right_only = int((~multi_zotu & single_species).sum())
        outside = int((~multi_zotu & ~single_species).sum())
        multi_total = int(multi_zotu.sum())
        records.append(
            {
                "method": method,
                "total_clusters": total,
                "multi_zotu_clusters": multi_total,
                "single_species_clusters": int(single_species.sum()),
                "pure_multi_zotu_clusters": overlap,
                "impure_multi_zotu_clusters": left_only,
                "pure_single_zotu_clusters": right_only,
                "other_clusters": outside,
                "pure_fraction_within_multi_zotu": (float(overlap) / float(multi_total)) if multi_total > 0 else np.nan,
                "impure_fraction_within_multi_zotu": (float(left_only) / float(multi_total)) if multi_total > 0 else np.nan,
            }
        )
    return pd.DataFrame.from_records(records)


def _plot_species_purity_venn_figure(
    venn_df: pd.DataFrame,
    out_svg: Path,
    font_size: float,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 6.2), constrained_layout=True)
    methods = ["TIC", "TAC"]
    left_color = "#56B4E9"
    right_color = "#009E73"

    for ax, method in zip(axes, methods):
        hit = venn_df[venn_df["method"] == method]
        if hit.empty:
            ax.text(0.5, 0.5, "No cluster data available", ha="center", va="center", fontsize=font_size)
            ax.axis("off")
            continue

        row = hit.iloc[0]
        left_only = int(row["impure_multi_zotu_clusters"])
        overlap = int(row["pure_multi_zotu_clusters"])
        right_only = int(row["pure_single_zotu_clusters"])
        outside = int(row["other_clusters"])
        multi_total = int(row["multi_zotu_clusters"])
        pure_frac = float(row["pure_fraction_within_multi_zotu"]) if pd.notna(row["pure_fraction_within_multi_zotu"]) else np.nan
        impure_frac = float(row["impure_fraction_within_multi_zotu"]) if pd.notna(row["impure_fraction_within_multi_zotu"]) else np.nan

        left_center = (0.43, 0.54)
        right_center = (0.61, 0.54)
        radius = 0.22
        ax.add_patch(plt.Circle(left_center, radius, color=left_color, alpha=0.42, ec="#2F2F2F", lw=1.2))
        ax.add_patch(plt.Circle(right_center, radius, color=right_color, alpha=0.42, ec="#2F2F2F", lw=1.2))

        ax.text(left_center[0] - 0.11, 0.54, f"{left_only}", ha="center", va="center", fontsize=font_size + 2.0, fontweight="bold")
        ax.text(0.52, 0.54, f"{overlap}", ha="center", va="center", fontsize=font_size + 2.0, fontweight="bold")
        ax.text(right_center[0] + 0.11, 0.54, f"{right_only}", ha="center", va="center", fontsize=font_size + 2.0, fontweight="bold")

        ax.text(left_center[0] - 0.02, 0.82, "Multi-zOTU\nclusters", ha="center", va="center", fontsize=max(8.0, font_size - 1.5))
        ax.text(right_center[0] + 0.02, 0.82, "Single-species\nclusters", ha="center", va="center", fontsize=max(8.0, font_size - 1.5))

        ax.text(
            0.5,
            0.20,
            (
                f"Pure among multi-zOTU = {overlap}/{multi_total}"
                + (f" ({pure_frac * 100:.1f}%)" if np.isfinite(pure_frac) else "")
            ),
            ha="center",
            va="center",
            fontsize=max(8.0, font_size - 1.5),
        )
        ax.text(
            0.5,
            0.14,
            (
                f"Impure among multi-zOTU = {left_only}/{multi_total}"
                + (f" ({impure_frac * 100:.1f}%)" if np.isfinite(impure_frac) else "")
            ),
            ha="center",
            va="center",
            fontsize=max(8.0, font_size - 1.5),
        )
        ax.text(
            0.5,
            0.08,
            f"Outside circles = {outside} single-zOTU non-pure clusters",
            ha="center",
            va="center",
            fontsize=max(7.5, font_size - 2.0),
            color="#555555",
        )

        ax.set_title(method, fontsize=font_size + 1.0)
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.0)
        ax.set_aspect("equal")
        ax.axis("off")

    fig.suptitle(
        "Figure 09 test: TIC/TAC multi-zOTU species purity as simple Venn summaries",
        fontsize=font_size + 2.0,
        y=1.02,
    )
    fig.savefig(out_svg, format="svg", bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Create a standalone Figure 09 test plot using the trial logic: "
            "Exact match = 100% effective identity; Near-exact match = threshold < best hit < 100%; "
            "No hit = <= threshold or no recorded hit."
        )
    )
    parser.add_argument("--config", required=True, help="Path to tacticbench YAML config")
    parser.add_argument(
        "--fig09-dir",
        default=None,
        help=(
            "Optional Figure 09 diagnostics directory. If omitted, the script uses "
            "figures/fig09_split_merge_diagnostics, falling back to "
            "figures/fig09_split_merge_diagnostics-TEST when available."
        ),
    )
    args = parser.parse_args()

    cfg = load_config(Path(args.config))
    outdir = Path(cfg.paths.outdir)
    if args.fig09_dir:
        fig09_dir = Path(args.fig09_dir)
    else:
        fig09_dir = figures_root(outdir) / "fig09_split_merge_diagnostics"
        test_fig09_dir = figures_root(outdir) / "fig09_split_merge_diagnostics-TEST"
        if (not (fig09_dir / "split_merge_diagnostics.tsv").exists()) and (
            test_fig09_dir / "split_merge_diagnostics.tsv"
        ).exists():
            fig09_dir = test_fig09_dir
    fig09_dir.mkdir(parents=True, exist_ok=True)

    registry = _load_registry(fig09_dir)
    feature_df = _build_feature_table(cfg, registry)
    sample_df = _summarize_by_sample(feature_df, registry)
    method_df = _summarize_by_method(sample_df)

    feature_out = fig09_dir / f"{OUTPUT_STEM}__features.tsv"
    sample_out = fig09_dir / f"{OUTPUT_STEM}__sample_summary.tsv"
    method_out = fig09_dir / f"{OUTPUT_STEM}__method_summary.tsv"
    svg_out = fig09_dir / f"{OUTPUT_STEM}.svg"
    richness_sample_out = fig09_dir / f"{RICHNESS_OUTPUT_STEM}__sample_summary.tsv"
    richness_method_out = fig09_dir / f"{RICHNESS_OUTPUT_STEM}__method_summary.tsv"
    richness_effect_out = fig09_dir / f"{RICHNESS_OUTPUT_STEM}__excess_removed.tsv"
    richness_effect_summary_out = fig09_dir / f"{RICHNESS_OUTPUT_STEM}__excess_removed_summary.tsv"
    richness_svg_out = fig09_dir / f"{RICHNESS_OUTPUT_STEM}.svg"

    feature_df.to_csv(feature_out, sep="\t", index=False)
    sample_df.to_csv(sample_out, sep="\t", index=False)
    method_df.to_csv(method_out, sep="\t", index=False)

    font_size = float(getattr(cfg.figures, "fig09_font_size", 14.0))
    _plot_summary(method_df, svg_out, font_size=font_size)
    richness_sample_df, richness_method_df, richness_effect_df, richness_effect_summary_df = (
        _build_relative_richness_tables(fig09_dir)
    )
    richness_sample_df.to_csv(richness_sample_out, sep="\t", index=False)
    richness_method_df.to_csv(richness_method_out, sep="\t", index=False)
    richness_effect_df.to_csv(richness_effect_out, sep="\t", index=False)
    richness_effect_summary_df.to_csv(richness_effect_summary_out, sep="\t", index=False)
    _plot_relative_richness_error_figure(
        richness_sample_df,
        richness_effect_summary_df,
        richness_svg_out,
        font_size=font_size,
    )

    transition_df, matrix_df, near_summary, flow_df = _build_transition_tables(feature_df)
    transition_rows_out = fig09_dir / f"{TRANSITION_OUTPUT_STEM}__rows.tsv"
    transition_matrix_out = fig09_dir / f"{TRANSITION_OUTPUT_STEM}__matrix.tsv"
    transition_near_out = fig09_dir / f"{TRANSITION_OUTPUT_STEM}__near_exact_summary.tsv"
    transition_svg_out = fig09_dir / f"{TRANSITION_OUTPUT_STEM}.svg"
    flow_rows_out = fig09_dir / f"{WAVE_OUTPUT_STEM}__rows.tsv"
    flow_summary_out = fig09_dir / f"{WAVE_OUTPUT_STEM}__summary.tsv"
    flow_html_out = fig09_dir / f"{WAVE_OUTPUT_STEM}.html"
    flow_svg_out = fig09_dir / f"{WAVE_OUTPUT_STEM}.svg"
    heatmap_svg_out = fig09_dir / f"{HEATMAP_OUTPUT_STEM}.svg"
    purity_cluster_out = fig09_dir / f"{PURITY_OUTPUT_STEM}__clusters.tsv"
    purity_sample_out = fig09_dir / f"{PURITY_OUTPUT_STEM}__samples.tsv"
    purity_stack_out = fig09_dir / f"{PURITY_OUTPUT_STEM}__stacked_summary.tsv"
    purity_svg_out = fig09_dir / f"{PURITY_OUTPUT_STEM}.svg"
    purity_venn_out = fig09_dir / f"{PURITY_VENN_OUTPUT_STEM}.svg"
    purity_venn_summary_out = fig09_dir / f"{PURITY_VENN_OUTPUT_STEM}__summary.tsv"
    match_species_cluster_out = fig09_dir / f"{MATCH_SPECIES_OUTPUT_STEM}__clusters.tsv"
    match_species_summary_out = fig09_dir / f"{MATCH_SPECIES_OUTPUT_STEM}__summary.tsv"
    match_species_heterogeneous_out = fig09_dir / f"{MATCH_SPECIES_OUTPUT_STEM}__heterogeneous_summary.tsv"
    match_species_svg_out = fig09_dir / f"{MATCH_SPECIES_OUTPUT_STEM}.svg"
    exact_near_cluster_out = fig09_dir / f"{EXACT_NEAR_SPECIES_OUTPUT_STEM}__clusters.tsv"
    exact_near_summary_out = fig09_dir / f"{EXACT_NEAR_SPECIES_OUTPUT_STEM}__summary.tsv"
    exact_near_method_summary_out = fig09_dir / f"{EXACT_NEAR_SPECIES_OUTPUT_STEM}__method_summary.tsv"
    exact_near_svg_out = fig09_dir / f"{EXACT_NEAR_SPECIES_OUTPUT_STEM}.svg"

    transition_df.to_csv(transition_rows_out, sep="\t", index=False)
    matrix_df.to_csv(transition_matrix_out, sep="\t", index=False)
    near_summary.to_csv(transition_near_out, sep="\t", index=False)
    _plot_transition_figure(matrix_df, near_summary, transition_svg_out, font_size=font_size)
    flow_df = flow_df.copy()
    if "zotu_class" in flow_df.columns:
        flow_df["zotu_class_display"] = flow_df["zotu_class"].astype(str).map(_wave_display_label)
    if "final_fate" in flow_df.columns:
        flow_df["final_fate_display"] = flow_df["final_fate"].astype(str).map(_wave_display_label)
    flow_df.to_csv(flow_rows_out, sep="\t", index=False)
    flow_summary = _build_flow_summary(flow_df)
    flow_summary = flow_summary.copy()
    if "zotu_class" in flow_summary.columns:
        flow_summary["zotu_class_display"] = flow_summary["zotu_class"].astype(str).map(_wave_display_label)
    if "final_fate" in flow_summary.columns:
        flow_summary["final_fate_display"] = flow_summary["final_fate"].astype(str).map(_wave_display_label)
    flow_summary.to_csv(flow_summary_out, sep="\t", index=False)
    _plot_wave_transition_figure(flow_summary, flow_html_out, font_size=font_size)
    _plot_wave_transition_svg(flow_summary, flow_svg_out, font_size=font_size)
    _plot_flow_transition_heatmap_figure(flow_summary, heatmap_svg_out, font_size=font_size)
    purity_cluster_df, purity_sample_df = _build_species_purity_tables(feature_df)
    purity_stack_df = _build_species_purity_stack_summary(purity_cluster_df)
    purity_venn_df = _build_species_purity_venn_summary(purity_cluster_df)
    purity_cluster_df.to_csv(purity_cluster_out, sep="\t", index=False)
    purity_sample_df.to_csv(purity_sample_out, sep="\t", index=False)
    purity_stack_df.to_csv(purity_stack_out, sep="\t", index=False)
    purity_venn_df.to_csv(purity_venn_summary_out, sep="\t", index=False)
    _plot_species_purity_figure(purity_cluster_df, purity_sample_df, purity_svg_out, font_size=font_size)
    _plot_species_purity_venn_figure(purity_venn_df, purity_venn_out, font_size=font_size)
    match_species_cluster_df, match_species_summary_df, match_species_heterogeneous_df = (
        _build_match_status_species_purity_tables(feature_df)
    )
    match_species_cluster_df.to_csv(match_species_cluster_out, sep="\t", index=False)
    match_species_summary_df.to_csv(match_species_summary_out, sep="\t", index=False)
    match_species_heterogeneous_df.to_csv(match_species_heterogeneous_out, sep="\t", index=False)
    _plot_match_status_species_purity_figure(
        match_species_summary_df,
        match_species_heterogeneous_df,
        match_species_svg_out,
        font_size=font_size,
    )
    exact_near_cluster_df, exact_near_summary_df, exact_near_method_summary_df = (
        _build_exact_near_species_purity_tables(match_species_cluster_df)
    )
    exact_near_cluster_df.to_csv(exact_near_cluster_out, sep="\t", index=False)
    exact_near_summary_df.to_csv(exact_near_summary_out, sep="\t", index=False)
    exact_near_method_summary_df.to_csv(exact_near_method_summary_out, sep="\t", index=False)
    _plot_exact_near_species_purity_figure(
        exact_near_summary_df,
        exact_near_method_summary_df,
        exact_near_svg_out,
        font_size=font_size,
    )

    print(f"Wrote: {svg_out}")
    print(f"Wrote: {method_out}")
    print(f"Wrote: {sample_out}")
    print(f"Wrote: {feature_out}")
    print(f"Wrote: {richness_svg_out}")
    print(f"Wrote: {richness_sample_out}")
    print(f"Wrote: {richness_method_out}")
    print(f"Wrote: {richness_effect_out}")
    print(f"Wrote: {richness_effect_summary_out}")
    print(f"Wrote: {transition_svg_out}")
    print(f"Wrote: {transition_matrix_out}")
    print(f"Wrote: {transition_near_out}")
    print(f"Wrote: {transition_rows_out}")
    print(f"Wrote: {flow_html_out}")
    print(f"Wrote: {flow_svg_out}")
    print(f"Wrote: {heatmap_svg_out}")
    print(f"Wrote: {flow_summary_out}")
    print(f"Wrote: {flow_rows_out}")
    print(f"Wrote: {purity_svg_out}")
    print(f"Wrote: {purity_cluster_out}")
    print(f"Wrote: {purity_sample_out}")
    print(f"Wrote: {purity_stack_out}")
    print(f"Wrote: {purity_venn_out}")
    print(f"Wrote: {purity_venn_summary_out}")
    print(f"Wrote: {match_species_svg_out}")
    print(f"Wrote: {match_species_cluster_out}")
    print(f"Wrote: {match_species_summary_out}")
    print(f"Wrote: {match_species_heterogeneous_out}")
    print(f"Wrote: {exact_near_svg_out}")
    print(f"Wrote: {exact_near_cluster_out}")
    print(f"Wrote: {exact_near_summary_out}")
    print(f"Wrote: {exact_near_method_summary_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
