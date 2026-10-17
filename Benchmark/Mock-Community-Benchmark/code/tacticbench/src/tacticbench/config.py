from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import re
from typing import Any, Dict, Optional

import yaml


@dataclass
class PathsConfig:
    root: Path
    outdir: Path


@dataclass
class NormalizationConfig:
    tss_total: float = 1000.0
    renormalize_after_mapping: bool = False


@dataclass
class MappingConfig:
    # Method thresholds are interpreted as effective percent identity:
    # effective_identity = pident * qcov_exact
    # where qcov_exact = (abs(qend - qstart) + 1) / qlen
    # Accepts percent (e.g., 98.7) or fraction (e.g., 0.987).
    min_effective_pident: float = 98.7
    # Optional dedicated threshold for TIC/TAC only.
    # If None, TIC/TAC use min_effective_pident_clustering (or min_effective_pident).
    min_effective_pident_tic_tac: Optional[float] = None
    min_effective_pident_clustering: Optional[float] = None
    # Optional BLAST output threshold for writing blast_hits.tsv via BLAST -perc_identity.
    # This is a RAW pident cutoff (not effective identity).
    # If None, defaults to the assignment threshold for each method
    # (min_effective_pident / min_effective_pident_clustering).
    # Accepts percent (e.g., 95.0) or fraction (e.g., 0.95).
    blast_record_min_raw_pident: Optional[float] = None
    # Additional BLAST quality gates to suppress short local matches.
    # min_query_coverage accepts percent (e.g., 90.0) or fraction (e.g., 0.9).
    min_query_coverage: float = 0.0
    seed: int = 1
    assignment_mode: str = "one_to_one"  # one_to_one | tie_aware_split | soft_species
    # soft_species mode (ASV/zOTU): weighting field for species candidates
    soft_species_score_field: str = "bitscore"  # bitscore | pident | length | uniform
    # keep species candidates with score >= (max_score * threshold), in [0,1]
    soft_species_min_relative_score: float = 0.0
    skip_community_on_samplemap_failure: bool = True
    sample_match_mode: str = "substring"  # currently only supported mode
    # If true, reduce to representative expected samples per community:
    # - if both Even and Staggered exist, keep one of each (first by order)
    # - otherwise keep one sample (first by order)
    # If false, keep all available expected samples/replicates.
    keep_one_sample_per_community: bool = False
    store_unmapped_mass: bool = True
    # Optional global mapping target level override by method.
    # Values: tip | species
    # Keys: asv, zotu, otu, tic, tac, default
    # If empty, defaults are used (ASV/zOTU=tip; OTU/TIC/TAC=species).
    target_level_by_method: dict[str, Any] = field(default_factory=dict)


@dataclass
class DiscoveryConfig:
    # Exclude communities whose folder name matches any pattern.
    # Pattern matching supports shell-style wildcards, e.g. "mock-ZIEL_*".
    # For plain text patterns (no wildcard), substring matching is also applied.
    exclude_community_name_patterns: list[str] = field(default_factory=list)
    # Expected reference filename to discover under each community source/ directory.
    # Set to "expected-sequences.v4.fasta" to run the full benchmark in V4 space.
    expected_sequences_filename: str = "expected-sequences.fasta"
    # If true and expected_sequences_filename is missing, fall back to
    # source/expected-sequences.fasta.
    expected_sequences_fallback_to_full: bool = True


@dataclass
class TreeConfig:
    build_full_tree: bool = True
    build_centroid_tree: bool = True
    mafft_exe: str = "mafft"
    fasttree_exe: str = "FastTree"
    blastn_exe: str = "blastn"
    makeblastdb_exe: str = "makeblastdb"


@dataclass
class DistancesConfig:
    methods: list[str] = field(default_factory=lambda: ["zOTU", "ASV", "TIC", "TAC", "OTU"])
    metrics: list[str] = field(
        default_factory=lambda: [
            "jaccard",
            "weighted_jaccard",
            "braycurtis",
            "unifrac",
            "gunifrac_a0.0",
            "gunifrac_a0.5",
            "gunifrac_a1.0",
        ]
    )
    penalize_unmapped: bool = False
    unmapped_penalty_lambda: float = 1.0
    penalize_mapping_mismatch: bool = False
    mapping_mismatch_lambda: float = 1.0
    penalize_coverage_deficit: bool = False
    coverage_deficit_lambda: float = 1.0
    penalty_transform: str = "linear"  # linear | exp
    penalty_exp_k: float = 3.0
    phylo_mapping_mode: str = "mapped_projection"  # mapped_projection | augmented_mismatch
    augmented_mismatch_beta: float = 1.0
    augmented_mismatch_eps: float = 1e-6
    augmented_mismatch_include_unmapped: bool = False
    # Zero-replacement pseudocount for Aitchison distance.
    aitchison_pseudocount: float = 1e-6


@dataclass
class FiguresConfig:
    fig01: bool = True
    fig02: bool = True
    fig03_by_sample: bool = True
    fig03_by_metric: bool = True
    fig04: bool = True
    fig05: bool = True
    fig06: bool = True
    fig07: bool = True
    fig08: bool = True
    fig09: bool = True
    fig10: bool = True
    fig11: bool = True
    fig12: bool = True
    fig13: bool = True
    fig14: bool = True
    fig15: bool = True
    fig16: bool = True
    fig17: bool = True
    fig18: bool = True
    fig19: bool = True
    fig20: bool = True
    fig21: bool = True
    fig22: bool = True
    fig23: bool = True
    fig24: bool = True
    fig25: bool = True
    fig26: bool = True
    grid_rows: int = 5
    grid_cols: int = 6
    # Marker area (matplotlib `s`) for jittered points on all box/violin figures.
    box_violin_dot_size: float = 22.0
    # Optional dedicated grid size for Figure 03 by-metric pages.
    # If unset, falls back to grid_rows/grid_cols.
    fig03_by_metric_grid_rows: Optional[int] = None
    fig03_by_metric_grid_cols: Optional[int] = None
    fig01_title: Optional[str] = None
    fig01_title_with_expected_features: Optional[str] = None
    fig02_title: Optional[str] = None
    fig03_by_sample_title: Optional[str] = None
    fig03_by_metric_title: Optional[str] = None
    fig04_title: Optional[str] = None
    fig05_title: Optional[str] = None
    fig06_title: Optional[str] = None
    fig07_title: Optional[str] = None
    fig08_title: Optional[str] = None
    fig09_title: Optional[str] = None
    fig10_title: Optional[str] = None
    fig11_title: Optional[str] = None
    # Optional Figure 11 metric selection override.
    # If empty, Figure 11 falls back to distances.metrics.
    fig11_metrics: list[str] = field(default_factory=list)
    # Optional Figure 11 metric-panel title overrides.
    # Keys use metric ids (e.g., jaccard, weighted_jaccard, braycurtis, unifrac, gunifrac_a0.0, gunifrac_a1.0).
    fig11_panel_title_by_metric: dict[str, str] = field(default_factory=dict)
    # Optional Figure 11 combined-figure suptitle overrides.
    fig11_combined_2x2_title: Optional[str] = None
    fig11_combined_3x2_title: Optional[str] = None
    fig12_title: Optional[str] = None
    fig13_title: Optional[str] = None
    fig14_title: Optional[str] = None
    fig15_title: Optional[str] = None
    fig16_title: Optional[str] = None
    fig17_title: Optional[str] = None
    fig18_title: Optional[str] = None
    fig19_title: Optional[str] = None
    fig20_title: Optional[str] = None
    fig21_title: Optional[str] = None
    fig22_title: Optional[str] = None
    fig23_title: Optional[str] = None
    fig24_title: Optional[str] = None
    fig25_title: Optional[str] = None
    fig26_title: Optional[str] = None
    # Base font size for Figure 11 labels/ticks/titles.
    fig11_font_size: float = 10.0
    # Base font size for Figure 12 labels/ticks/titles.
    fig12_font_size: float = 10.0
    # Base font size for Figure 13 labels/ticks/titles.
    fig13_font_size: float = 10.0
    # Base font size for Figure 14 labels/ticks/titles.
    fig14_font_size: float = 10.0
    # Base font size for Figure 15 labels/ticks/titles.
    fig15_font_size: float = 10.0
    # Base font size for Figure 16 labels/ticks/titles.
    fig16_font_size: float = 10.0
    # Base font size for Figure 17 labels/ticks/titles.
    fig17_font_size: float = 10.0
    # Base font size for Figure 18 labels/ticks/titles.
    fig18_font_size: float = 10.0
    # Base font size for Figure 19 labels/ticks/titles.
    fig19_font_size: float = 10.0
    # Base font size for Figure 20 labels/ticks/titles.
    fig20_font_size: float = 10.0
    # Base font size for Figure 21 labels/ticks/titles.
    fig21_font_size: float = 10.0
    # Optional font-size overrides for Figures 22-26 (fallback to fig21_font_size).
    fig22_font_size: Optional[float] = None
    fig23_font_size: Optional[float] = None
    fig24_font_size: Optional[float] = None
    fig25_font_size: Optional[float] = None
    fig26_font_size: Optional[float] = None
    # Optional Figure 21 per-panel title overrides.
    # Leave as None to use the built-in panel titles.
    fig21_show_panel_a: bool = True
    fig21_show_panel_b: bool = True
    fig21_show_panel_c: bool = True
    fig21_show_panel_d: bool = True
    fig21_show_panel_e: bool = True
    fig21_show_panel_f: bool = True
    fig21_panel_a_title: Optional[str] = None
    fig21_panel_a_a_title: Optional[str] = None
    fig21_panel_a_b_title: Optional[str] = None
    fig21_panel_b_title: Optional[str] = None
    fig21_panel_c_title: Optional[str] = None
    fig21_panel_c_a_title: Optional[str] = None
    fig21_panel_c_b_title: Optional[str] = None
    fig21_panel_d_title: Optional[str] = None
    fig21_panel_e_title: Optional[str] = None
    fig21_panel_f_title: Optional[str] = None
    # Optional Figure 21 per-panel font-size overrides.
    # Leave as None to fall back to fig21_font_size.
    fig21_panel_a_font_size: Optional[float] = None
    fig21_panel_b_font_size: Optional[float] = None
    fig21_panel_c_font_size: Optional[float] = None
    fig21_panel_d_font_size: Optional[float] = None
    fig21_panel_e_font_size: Optional[float] = None
    fig21_panel_f_font_size: Optional[float] = None
    # Optional Figure 21 Panel E (KDE+boxplot grid) font-size overrides.
    # Leave as None to derive from fig21_panel_e_font_size.
    fig21_panel_e_grid_label_font_size: Optional[float] = None
    fig21_panel_e_grid_tick_font_size: Optional[float] = None
    fig21_panel_e_grid_text_font_size: Optional[float] = None
    fig21_panel_e_grid_legend_font_size: Optional[float] = None
    # Optional Figure 21 Panel B legend label overrides for sample V4 status.
    fig21_panel_b_sample_v4_multiple_label: str = "Multiple-V4 reference"
    fig21_panel_b_sample_v4_single_label: str = "Single-V4 reference"
    # Optional Figure 21 sample-dot styling overrides for Panels B, D, and E.
    # Dot sizes are Matplotlib scatter areas. Leave as None to derive from fig21_point_size.
    fig21_panel_b_dot_size: Optional[float] = None
    fig21_panel_b_dot_opacity: float = 0.48
    fig21_panel_b_dot_jitter: float = 0.065
    fig21_panel_d_dot_size: Optional[float] = None
    fig21_panel_d_dot_opacity: float = 0.48
    fig21_panel_d_dot_jitter: float = 0.07
    fig21_panel_d_stats_show_annotations: bool = True
    fig21_panel_e_dot_size: Optional[float] = None
    fig21_panel_e_dot_opacity: float = 0.95
    fig21_panel_e_dot_jitter: float = 0.055
    fig21_panel_f_dot_size: Optional[float] = None
    fig21_panel_f_dot_opacity: float = 0.75
    fig21_panel_f_dot_jitter: float = 0.06
    # Optional Panel E legend label overrides.
    fig21_panel_e_abs_density_label: str = "Pairwise abs_difference density"
    fig21_panel_e_summary_density_label: str = "Sample summary_abs_difference density"
    fig21_panel_e_show_summary_density: bool = True
    fig21_panel_e_mean_label: str = "Mean"
    fig21_panel_e_sample_v4_multiple_label: Optional[str] = None
    fig21_panel_e_sample_v4_single_label: Optional[str] = None
    # Optional Figure 21 axis-title overrides for Panels A-E.
    fig21_panel_a_x_axis_label: Optional[str] = None
    fig21_panel_a_y_axis_label: Optional[str] = None
    fig21_panel_b_x_axis_label: Optional[str] = None
    fig21_panel_b_y_axis_label: Optional[str] = None
    fig21_panel_c_x_axis_label: Optional[str] = None
    fig21_panel_c_y_axis_label: Optional[str] = None
    fig21_panel_d_x_axis_label: Optional[str] = None
    fig21_panel_d_y_axis_label: Optional[str] = None
    fig21_panel_e_x_axis_label: Optional[str] = None
    fig21_panel_e_y_axis_label: Optional[str] = None
    fig21_panel_f_x_axis_label: Optional[str] = None
    fig21_panel_f_y_axis_label: Optional[str] = None
    fig21_panel_f_xtick_rotation: float = 0.0
    # Optional panel-B subset of pair-groups to display (pair_group_id strings).
    fig21_panel_b_pair_groups_to_show: list[str] = field(default_factory=list)
    # If true, annotate Panel B with paired tests among the within-group contrast distributions.
    fig21_panel_b_show_within_contrast_stats: bool = False
    # Rotation in degrees for Figure 21 Panel A x/y tick labels.
    fig21_panel_a_xtick_rotation: float = 35.0
    fig21_panel_a_ytick_rotation: float = 0.0
    # If true, overlay Panel B planned-significance comparisons onto Figure 21 Panel A_b
    # as connecting lines between relevant heatmap cells with p-values and star levels.
    fig21_panel_a_show_planned_comparison_annotations: bool = False
    # If true, include Figure 21 Panel E (KDE + horizontal boxplot grid of pair-group differences).
    fig21_show_panel_e_raincloud_grid: bool = True
    # Rotation in degrees for Figure 21 Panel B x/y tick labels.
    fig21_panel_b_xtick_rotation: float = 35.0
    fig21_panel_b_ytick_rotation: float = 0.0
    # Number of blank cells to insert between method blocks in Figure 21 Panel C.
    fig21_endpoint_matrix_gap_cells: int = 1
    # Rotation in degrees for Figure 21 Panel C endpoint heatmap x tick labels.
    fig21_endpoint_matrix_xtick_rotation: float = 90.0
    # Rotation in degrees for Figure 21 Panel C endpoint heatmap y tick labels.
    fig21_endpoint_matrix_ytick_rotation: float = 0.0
    # Optional Panel C subpanel rotations. Leave as None to use panel C defaults.
    fig21_panel_c_a_xtick_rotation: Optional[float] = None
    fig21_panel_c_a_ytick_rotation: Optional[float] = None
    fig21_panel_c_b_xtick_rotation: Optional[float] = None
    fig21_panel_c_b_ytick_rotation: Optional[float] = None
    # Marker area for Figure 17 sample dots.
    fig17_dot_size: float = 32.0
    # Fill opacity for Figure 17 density ridgeline fills.
    fig17_opacity: float = 1.0
    # If false, hide the upper sample-point panel in Figure 17.
    fig17_show_upper_panel: bool = True
    # Figure 17 method display mode.
    fig17_method_view: str = "methods"  # methods | families
    # Figure 17 lower panel orientation.
    fig17_lower_panel_y: str = "metrics"  # metrics | methods
    # Figure 17 lower panel rendering mode.
    fig17_lower_panel: str = "density"  # density | bars
    # Figure 17 grid arrangement:
    # - method_pipeline: one panel per method x metric combination
    # - facet: one panel per facet value, with rows inside each panel
    fig17_grid_mode: str = "method_pipeline"  # method_pipeline | facet
    # Relative height of the lower panel in Figure 17.
    fig17_lower_panel_height_ratio: float = 6.0
    # Optional Figure 17 KDE controls.
    fig17_density_bandwidth: Optional[float] = None
    fig17_density_points: int = 256
    fig17_density_scale: float = 0.72
    # Figure 17 x-axis major tick interval.
    fig17_xtick_interval: float = 0.2
    # Figure 17 title padding for per-panel metric labels.
    fig17_metric_title_pad: float = 14.0
    # Optional Figure 17 metric selection.
    fig17_metrics: list[str] = field(default_factory=list)
    # Figure 17 pairwise significance settings.
    fig17_stats_enabled: bool = True
    fig17_stats_alpha: float = 0.05
    fig17_stats_pairwise_test: str = "signed_rank_permutation"  # signed_rank_permutation | sign_test_exact
    fig17_stats_multiple_test_correction: bool = True
    fig17_stats_family_aggregation: str = "median"  # median | mean | min
    # Figure 17 significance annotation style.
    fig17_stats_annotation_mode: str = "letters"  # letters | brackets
    fig17_stats_show_legend: bool = True
    # If true, keep "missing_similarity_data" as a separate category in the
    # left panel of Figure 14. If false, hide it from the plot/legend.
    fig14_show_missing_similarity_data_category: bool = False
    # Figure 14-specific edge thresholds (effective identity, percent unless in [0,1]).
    # These are independent from Figure 09 thresholds.
    fig14_graph_min_effective_pident_denoising: float = 100.0
    fig14_graph_min_effective_pident_tic_tac: float = 98.7
    fig14_graph_min_effective_pident_otu: float = 97.0
    # Optional Figure 14 per-method overrides:
    # keys: asv, zotu, otu, tic, tac, default
    # values: percent (98.7) or fraction (0.987)
    fig14_edge_threshold_by_method: dict[str, Any] = field(default_factory=dict)
    # Figure 14-specific min query coverage gate (percent or fraction).
    # If None, falls back to fig09_min_query_coverage for backward compatibility.
    fig14_min_query_coverage: Optional[float] = None
    # Optional Figure 14 per-method overrides for min query coverage.
    # keys: asv, zotu, otu, tic, tac, default
    # values: percent (100.0) or fraction (1.0)
    fig14_min_query_coverage_by_method: dict[str, Any] = field(default_factory=dict)
    # Figure 16 per-community species driver plot:
    # number of species to show, ranked by max absolute deviation across methods.
    fig16_species_driver_top_n: int = 25
    # Figure 16 methods to include. Empty means all methods.
    # Accepted aliases include dada2/asv, unoise3/zotu, uparse/otu, tic, tac.
    fig16_methods: list[str] = field(default_factory=list)
    # Figure 16 per-community species driver plot:
    # aggregate sample-wise deviations per community using mean or median.
    fig16_species_driver_sample_aggregation: str = "mean"  # mean | median
    # Figure 16 per-community species driver plot labels.
    # title_prefix is rendered as: "<prefix> - <community>"
    fig16_species_driver_title_prefix: str = "Species Driver Plot"
    fig16_species_driver_xlabel: str = "Relative Abundance Deviation (Inferred - Expected)"
    # Supports format placeholders {n} or {n_species}.
    fig16_species_driver_ylabel_template: str = "Top {n} Species by |Deviation|"
    # Figure 16 shared tick-label rotation controls for both species and V4 panels.
    fig16_xtick_rotation: float = 45.0
    fig16_ytick_rotation: float = 45.0
    # Figure 16 x-axis limit mode:
    # - "auto": global shared limit from longest plotted range across communities
    # - "unbound": per-community local limit
    # - numeric: fixed percent (e.g., 20 => +/-0.20; 0.2 => +/-0.20)
    fig16_species_driver_xlim: Any = "unbound"
    # Figure 16 abundance scale used to compute deviations:
    # - relative: normalize projected and expected tables per sample before subtraction
    # - raw: subtract raw projected abundances from raw expected abundances
    # - original_relative: use mapped abundance from the original inferred table
    #   divided by the original inferred sample total, and compare that to the
    #   expected relative abundance
    fig16_projection_abundance_mode: str = "relative"  # relative | raw | original_relative
    fig07_plot_type: str = "boxplot"  # boxplot | violin
    # Figure 07 expected-level target by method:
    # keys: asv, zotu, otu, tic, tac, default
    # values: tip | species
    # If empty, defaults are ASV/zOTU=tip and OTU/TIC/TAC=species.
    fig07_target_level_by_method: dict[str, Any] = field(default_factory=dict)
    # Figure 07 stats engine:
    # - mixed_effect_glmmtmb: Gaussian mixed model via glmmTMB + emmeans
    # - paired_permutation: permutation-based Kruskal-style omnibus test
    fig07_stats_model: str = "mixed_effect_glmmtmb"  # mixed_effect_glmmtmb | paired_permutation
    # Figure 07 pairwise test used when fig07_stats_model is paired_permutation.
    # - signed_rank_permutation: permutation-based Wilcoxon signed-rank test
    # - sign_test_exact: exact paired sign test (robust to many zeros/ties)
    fig07_stats_pairwise_test: str = "signed_rank_permutation"  # signed_rank_permutation | sign_test_exact
    # emmeans scale used by Figure 07 GLMM pairwise contrasts.
    fig07_stats_glmm_emmeans_scale: str = "response"  # link | response
    # alpha threshold used when flagging pairwise significance in Figure 07 stats tables.
    fig07_stats_alpha: float = 0.05
    # If true, apply BH correction across Figure 07 pairwise tests per metric.
    fig07_stats_multiple_test_correction: bool = True
    # If true, draw pairwise significance brackets on Figure 07 panels.
    fig07_stats_show_pairwise: bool = True
    # If true, also write the 3-metric Figure 07 variant:
    # observed_species + effective_shannon + effective_simpson.
    fig07_include_richness_effective_variant: bool = True
    # If true, write alpha_diversity_errors_long.tsv beside wide TSV output.
    fig07_write_long_tsv: bool = True
    fig08_plot_type: str = "boxplot"  # boxplot | violin
    fig09_plot_type: str = "boxplot"  # boxplot | violin
    # Base font size for Figure 09 labels/ticks/titles.
    fig09_font_size: float = 10.0
    # Horizontal padding between Figure 09 subplot columns.
    fig09_layout_w_pad: float = 3.8
    # Figure 09 outer layout rectangle passed to tight_layout:
    # [left, bottom, right, top] in figure coordinates.
    fig09_layout_left: float = 0.02
    fig09_layout_bottom: float = 0.02
    fig09_layout_right: float = 0.96
    fig09_layout_top: float = 0.95
    # Optional final SVG cropping.
    fig09_save_bbox_tight: bool = False
    fig09_save_pad_inches: float = 0.02
    # If true, write one SVG per Figure 09 panel beside the combined figure.
    fig09_write_individual_panels: bool = True
    # Split-index source for Figure 09:
    # - mapping: use finalized mapping.tsv assignments (legacy behavior)
    # - threshold_hits: use all BLAST hits above method threshold
    # - graph: graph-based split/merge from BLAST-hit bipartite edges
    fig09_split_index_mode: str = "mapping"  # mapping | threshold_hits | graph
    # target aggregation for threshold_hits mode:
    # - species: collapse expected tips to species before counting splits
    # - tip: operate on expected tip ids directly
    fig09_split_index_threshold_target_level: str = "species"  # species | tip
    # Split-index denominator for Figure 09.
    # - mapped_expected_features: number of expected targets with >=1 mapped edge/assignment
    # - mapped_inferred_features: number of mapped inferred features
    fig09_split_index_denominator: str = "mapped_expected_features"  # mapped_expected_features | mapped_inferred_features
    # 1-to-1 mapping fraction denominator for Figure 09.
    # - expected_present_features: expected targets present in sample
    # - mapped_inferred_features: mapped inferred features in sample
    # - inferred_present_features: inferred present features in sample
    fig09_one_to_one_denominator: str = "expected_present_features"  # expected_present_features | mapped_inferred_features | inferred_present_features
    # Merge index denominator for Figure 09 graph mode.
    # - mapped_inferred_features: mapped inferred features (legacy/current behavior)
    # - expected_present_features: expected species present in sample
    fig09_merge_denominator: str = "mapped_inferred_features"  # mapped_inferred_features | expected_present_features
    # Graph mode edge thresholds on effective identity:
    #   effective_identity = pident * qcov_exact
    #   where qcov_exact = (abs(qend - qstart) + 1) / qlen
    # Values are interpreted as percent identity by default.
    # If a value is in [0,1], it is treated as fraction and multiplied by 100.
    fig09_graph_min_effective_pident_denoising: float = 100.0
    fig09_graph_min_effective_pident_tic_tac: float = 98.7
    fig09_graph_min_effective_pident_otu: float = 97.0
    # Additional hit quality gates for Figure 09 hit-based modes.
    # min_query_coverage accepts percent (e.g., 90.0) or fraction (e.g., 0.9).
    fig09_min_query_coverage: float = 0.0
    # Optional Figure 09-specific query-coverage gate for no-hit/unmapped logic.
    # If None, fallback order is:
    #   fig12_unmapped_best_hit_min_query_coverage -> fig09_min_query_coverage.
    fig09_no_hit_min_query_coverage: Optional[float] = None
    # Optional Figure 09 method-specific no-hit thresholds on effective identity.
    # If a method is provided, no-hit/unmapped classification in Figure 09 uses it.
    # Missing methods fall back to Figure 12 method thresholds.
    # Methods: asv, zotu, otu, tic, tac (+ optional "default" fallback).
    # Values accept percent (98.7) or fraction (0.987).
    fig09_no_hit_min_effective_pident_by_method: dict[str, float] = field(default_factory=dict)
    # Optional dedicated hit-quality gate for Figure 12 diagnostics
    # (used while building exploratory unmapped best-hit identities).
    # If None, falls back to fig09_min_query_coverage.
    fig12_unmapped_best_hit_min_query_coverage: Optional[float] = None
    # Optional method-specific no-hit thresholds on effective identity for Figure 12.
    # If provided for a method, an unmapped feature is counted as "has best hit"
    # only when best_hit_effective_pident >= threshold for that method.
    # Methods: asv, zotu, otu, tic, tac (+ optional "default" fallback).
    # Values accept percent (98.7) or fraction (0.987).
    fig12_no_hit_min_effective_pident_by_method: dict[str, float] = field(default_factory=dict)
    # Optional method-specific no-hit thresholds on effective identity for Figure 15.
    # If provided for a method, weighted no-hit vs unmapped-with-hit split uses this cutoff.
    # Methods: asv, zotu, otu, tic, tac (+ optional "default" fallback).
    # Values accept percent (98.7) or fraction (0.987).
    fig15_no_hit_min_effective_pident_by_method: dict[str, float] = field(default_factory=dict)
    # Graph mode target level per method (tip | species).
    fig09_graph_target_level_asv: str = "tip"
    fig09_graph_target_level_zotu: str = "tip"
    fig09_graph_target_level_tic: str = "species"
    fig09_graph_target_level_tac: str = "species"
    fig09_graph_target_level_otu: str = "species"
    # Optional graph-mode override by metric and method for Figure 09.
    # Shape:
    #   {
    #     "mode": "single" | "tip_species_split_violin",
    #     "dual_level_layout": "separate" | "combined",
    #     "split_index": {"asv":"tip", "zotu":"tip", "tic":"species", ...},
    #     "merge_deficit": {"asv":"species", ...},
    #     "unmapped_feature_fraction": {"asv":"tip", ...},
    #     "one_to_one_mapping_fraction": {"asv":"tip", ...}
    #   }
    # mode behavior:
    # - single: current behavior (one row per sample/method)
    # - tip_species_split_violin: graph mode only; compute rows at both levels
    #   and then:
    #   * dual_level_layout=separate: write two figures (tip and species)
    #   * dual_level_layout=combined: write one combined 8-panel figure
    # Metric keys also accept aliases:
    #   merge_index, merge, one_to_one, one_to_one_mapping, unmapped
    fig09_graph_target_level_by_metric: dict[str, Any] = field(default_factory=dict)
    # Panel-specific target-level overrides for the Figure 09 a/b/c layout.
    # Shape:
    #   {
    #     "a": {"split_index": {"default":"species"}, "merge_deficit": {"default":"species"}},
    #     "b": {"one_to_one_mapping_fraction": {"default":"species"},
    #           "cluster_purity_fraction": {"default":"species"}},
    #     "c": {"unmapped_feature_fraction": {"default":"species"},
    #           "no_hit_over_unmapped_fraction": {"default":"species"}}
    #   }
    # Panel keys also accept aliases: panel_a/panel_b/panel_c.
    # Shorthand supported:
    #   a: species
    #   b: species
    #   c: {asv: tip, zotu: tip, default: species}
    # where a panel-level method map applies to all metrics in that panel.
    fig09_panel_target_level_by_metric: dict[str, Any] = field(default_factory=dict)
    # Optional Figure 09 graph-mode edge-threshold overrides by metric/method.
    # Shape:
    #   {
    #     "split_index": {"default":98.7, "asv":100.0, "zotu":100.0, "otu":97.0, "tic":98.7, "tac":98.7},
    #     "merge_deficit": {...},
    #     "cluster_purity_fraction": {...},
    #     "unmapped_feature_fraction": {...},
    #     "one_to_one_mapping_fraction": {...},
    #     "no_hit_over_unmapped_fraction": {...}
    #   }
    # Metric keys accept aliases:
    #   split, merge, merge_index, one_to_one, one_to_one_mapping, unmapped
    # Method keys accept aliases:
    #   asv/dada2, zotu/unoise3, otu/uparse, tic, tac
    # Shorthand supported:
    #   {asv:100, zotu:100, otu:97, tic:98.7, tac:98.7}
    #   {default:98.7}
    # i.e., a global method-map (or single default) for all metrics.
    # Values accept percent (98.7) or fraction (0.987).
    fig09_graph_edge_threshold_by_metric: dict[str, Any] = field(default_factory=dict)
    # Legacy knob from old 4-panel layout. Ignored by current 3-panel Figure 09.
    fig09_panel_c_enabled: bool = True
    # Optional custom titles for Figure 09 panel rows:
    # keys: a,b,c (or panel_a..panel_c)
    fig09_panel_titles: dict[str, Any] = field(default_factory=dict)
    # Optional per-metric text overrides for Figure 09 2x2 graph layout.
    # Shape:
    #   {
    #     "split_index": {"title":"...", "xlabel":"...", "ylabel":"..."},
    #     "merge_deficit": {"title":"...", "xlabel":"...", "ylabel":"..."},
    #     "one_to_one_mapping_fraction": {"title":"...", "xlabel":"...", "ylabel":"..."},
    #     "no_hit_over_unmapped_fraction": {"title":"...", "xlabel":"...", "ylabel":"..."}
    #   }
    # Shorthand:
    #   split_index: "Custom Title"
    # (sets title only)
    fig09_graph_labels_by_metric: dict[str, Any] = field(default_factory=dict)
    # Figure 09 panel D: relative richness error against expected species richness.
    fig09_panel_d_title: Optional[str] = None
    fig09_panel_d_xlabel: Optional[str] = None
    fig09_panel_d_ylabel: Optional[str] = None
    fig09_panel_d_methods: list[str] = field(default_factory=lambda: ["ASV", "zOTU", "TAC", "TIC", "OTU"])
    fig09_panel_d_plot_type: str = "boxplot"
    fig09_panel_d_font_size: Optional[float] = None
    fig09_panel_d_tick_size: Optional[float] = None
    fig09_panel_d_xtick_rotation: float = 45.0
    fig09_panel_d_ytick_rotation: float = 0.0
    fig09_panel_d_point_size: Optional[float] = None
    fig09_panel_d_point_opacity: float = 0.72
    fig09_panel_d_point_jitter: float = 0.055
    fig09_panel_d_show_zero_line: bool = True
    fig09_panel_d_show_reference_note: bool = True
    fig09_panel_d_reference_note: str = "0% = expected species richness"
    fig09_panel_d_show_excess_removed_annotation: bool = True
    fig09_panel_d_y_limits: list[float] = field(default_factory=list)
    # Replacement panel for no_hit_over_unmapped_fraction in Figure 09:
    # a compact TIC/TAC zOTU transition heatmap.
    fig09_transition_heatmap_title: Optional[str] = None
    fig09_transition_heatmap_xlabel: Optional[str] = None
    fig09_transition_heatmap_ylabel: Optional[str] = None
    fig09_transition_heatmap_note: str = ""
    fig09_transition_heatmap_cmap: str = "YlGnBu"
    fig09_transition_heatmap_show_counts: bool = True
    fig09_transition_heatmap_show_percentages: bool = True
    fig09_transition_heatmap_show_colorbar: bool = True
    fig09_transition_heatmap_x_tick_rotation: float = 0.0
    fig09_transition_heatmap_row_label_rotation: float = 0.0
    fig09_transition_heatmap_method_label_x_offset: float = 0.05
    fig09_transition_heatmap_colorbar_pad: float = 0.18
    fig09_transition_heatmap_tic_label: str = "TIC"
    fig09_transition_heatmap_tac_label: str = "TAC"
    fig09_transition_heatmap_no_match_label: str = "No\nmatch"
    fig09_transition_heatmap_non_bacterial_label: str = "Non\nbacterial"
    fig09_transition_heatmap_colorbar_label: str = "Row fraction of zOTUs"
    fig09_transition_heatmap_exact_identity: float = 100.0
    # If None, the zOTU no-hit threshold is used.
    fig09_transition_heatmap_near_exact_min_effective_pident: Optional[float] = None
    # Companion Figure 09 homogeneity plot written beside the main Figure 09.
    fig09_homogeneity_enabled: bool = True
    fig09_homogeneity_title: str = "Taxonomy homogeneity across cluster match-status groups"
    fig09_homogeneity_font_size: Optional[float] = None
    fig09_homogeneity_tick_size: Optional[float] = None
    fig09_homogeneity_figsize: list[float] = field(default_factory=list)
    fig09_homogeneity_count_ylabel: str = "Cluster count"
    fig09_homogeneity_effect_ylabel: str = "Abundance effect"
    fig09_homogeneity_count_label: str = "Count"
    fig09_homogeneity_effect_label: str = "Abundance\neffect"
    fig09_homogeneity_count_legend_label: str = "Absolute count"
    fig09_homogeneity_single_species_label: str = "Single species"
    fig09_homogeneity_multi_species_label: str = "Multi-species or unresolved"
    fig09_homogeneity_single_species_color: str = "#009E73"
    fig09_homogeneity_multi_species_color: str = "#D55E00"
    fig09_homogeneity_effect_hatch: str = "..."
    fig09_homogeneity_effect_alpha: float = 0.78
    fig09_homogeneity_effect_axis_max: float = 1.0
    fig09_homogeneity_count_axis_padding: float = 1.25
    fig09_homogeneity_show_value_labels: bool = True
    fig09_homogeneity_show_bar_type_labels: bool = True
    fig09_homogeneity_show_legends: bool = True
    fig09_homogeneity_pure_exact_title: str = "Pure exact"
    fig09_homogeneity_pure_near_title: str = "Pure near-exact"
    fig09_homogeneity_pure_no_match_title: str = "Pure no match"
    fig09_homogeneity_exact_near_title: str = "Exact + near-exact"
    fig09_homogeneity_save_bbox_tight: bool = False
    fig09_homogeneity_save_pad_inches: float = 0.02
    fig11_plot_type: str = "violin"  # boxplot | violin
    fig18_plot_type: str = "violin"  # boxplot | violin
    # Figure 19 point size.
    fig19_point_size: float = 36.0
    # Optional Figure 19 metric selection (Aitchison is always excluded for Fig19).
    fig19_metrics: list[str] = field(default_factory=list)
    # Distance metric between contrast profiles for Figure 19.
    fig19_profile_distance: str = "euclidean"  # euclidean
    # If true, z-score contrast profile columns before computing pairwise distances.
    fig19_standardize_profiles: bool = True
    # If true, render a UMAP panel next to PCoA in Figure 19.
    fig19_show_umap: bool = True
    # UMAP hyperparameters for Figure 19 (metric=precomputed on contrast-distance matrix).
    fig19_umap_n_neighbors: int = 15
    fig19_umap_min_dist: float = 0.10
    fig19_umap_random_state: int = 42
    # Figure 20 sample-wise aggregation mode for component relative abundances.
    fig20_sample_aggregation: str = "mean"  # mean | median
    # Figure 20 methods to include. Empty means all methods.
    # Accepted aliases include dada2/asv, unoise3/zotu, uparse/otu, tic, tac.
    fig20_methods: list[str] = field(default_factory=list)
    # Figure 21 metric-family sensitivity diagnostics.
    fig21_plot_type: str = "violin"  # boxplot | violin
    # Violin KDE bandwidth control for Figure 21 pair-group distributions.
    # Use None for the Matplotlib default (Scott's rule).
    # Larger numeric values smooth more; smaller values preserve more local detail.
    fig21_violin_bw_method: float | str | None = None
    fig21_point_size: float = 42.0
    fig21_metrics: list[str] = field(default_factory=list)
    fig21_agnostic_metrics: list[str] = field(default_factory=list)
    fig21_aware_metrics: list[str] = field(default_factory=list)
    fig21_profile_distance: str = "euclidean"  # euclidean
    fig21_standardize_profiles: bool = True
    fig21_show_pair_group_distance_contrast: bool = True
    fig21_pair_group_summary_stat: str = "median"  # median | mean
    # Figure 21 panel-B method-family membership. Keys are denoising and
    # clustering; empty values retain the historical ASV/zOTU and OTU/TIC/TAC groups.
    fig21_panel_b_method_groups: dict[str, list[str]] = field(default_factory=dict)
    # Optional display labels for the four Panel B endpoint groups. Keys are:
    # denoising_phylogeny_aware, denoising_phylogeny_agnostic,
    # clustering_phylogeny_aware, clustering_phylogeny_agnostic.
    fig21_panel_b_group_labels: dict[str, str] = field(default_factory=dict)
    # Figure 21 panels D/F method filters. Empty means all methods in the
    # Figure 21 distance table.
    fig21_panel_d_methods: list[str] = field(default_factory=list)
    fig21_panel_f_methods: list[str] = field(default_factory=list)
    # Panel F summary line(s) over each violin/boxplot.
    fig21_panel_f_summary_statistic: str = "median"  # median | mean | both
    # Panel F pairwise test. mean_difference_permutation is a sign-flip test
    # of the paired arithmetic-mean difference (method_a - method_b).
    # None inherits fig21.stats.pairwise_test; a panel_f.stats value overrides it.
    fig21_panel_f_stats_pairwise_test: str | None = None
    fig21_panel_f_stats_alternative: str = "two-sided"  # two-sided | greater | less
    # Positive minimum paired mean difference required by a one-sided Panel F
    # test. A value of 0.0 retains the conventional null of no difference.
    fig21_panel_f_stats_practical_difference_delta: float = 0.0
    fig21_panel_f_stats_bootstrap_resamples: int = 4000
    fig21_panel_f_stats_ci_level: float = 0.95
    fig21_panel_f_stats_show_annotations: bool = True
    fig21_panel_f_stats_show_non_significant_annotations: bool = True
    # Optional planned pairs of methods for Panel F. Empty tests all pairs.
    fig21_panel_f_stats_planned_comparisons: list[list[str]] = field(default_factory=list)
    fig21_show_sensitivity_by_method: bool = True
    fig21_show_sample_method_heatmap: bool = True
    fig21_show_per_sample_method_metric_pcoa: bool = True
    fig21_show_method_metric_profile_pcoa: bool = True
    fig21_show_pairwise_difference_pcoa: bool = True
    # Figure 11-only projection mapping mode:
    # - standard: reuse module default projection behavior
    # - best_hit_aggregate: use best-hit tie-splitting at tip level with UNKNOWN bin,
    #   and species-hit aggregation at species level.
    fig11_projection_mapping_mode: str = "standard"  # standard | best_hit_aggregate
    # Figure 11-only mapping overrides (do not affect module-wide mapping):
    # target level per method: tip | species
    # keys: asv, zotu, otu, tic, tac, default
    fig11_target_level_by_method: dict[str, Any] = field(default_factory=dict)
    # effective-identity threshold per method (percent or fraction)
    # keys: asv, zotu, otu, tic, tac, default
    fig11_min_effective_pident_by_method: dict[str, Any] = field(default_factory=dict)
    # Figure 21-only projection/mapping overrides, same semantics as Figure 11.
    fig21_projection_mapping_mode: str = "standard"  # standard | best_hit_aggregate
    fig21_target_level_by_method: dict[str, Any] = field(default_factory=dict)
    fig21_min_effective_pident_by_method: dict[str, Any] = field(default_factory=dict)
    fig09_stats_enabled: bool = True
    fig09_stats_permutations: int = 4000
    fig09_stats_seed: int = 1
    fig09_stats_alpha: float = 0.05
    fig09_stats_show_pairwise: bool = True
    # Figure 09 stats engine:
    # - mixed_effect_glmmtmb: GLMMs via R (glmmTMB + emmeans), stronger formal tests
    # - mixed_effect: fallback GLMMs via R (MASS::glmmPQL) per metric panel
    # - paired_permutation: legacy Friedman + paired permutation/sign tests
    fig09_stats_model: str = "mixed_effect_glmmtmb"  # mixed_effect_glmmtmb | mixed_effect | paired_permutation
    # Pairwise test in Figure 09:
    # - signed_rank_permutation: paired signed-rank permutation test
    # - sign_test_exact: exact paired sign test (robust to many zeros/ties)
    fig09_stats_pairwise_test: str = "signed_rank_permutation"  # signed_rank_permutation | sign_test_exact
    # Customizable GLMM (glmmTMB backend) model specs for Figure 09.
    # Binomial panels:
    #   - unmapped_feature_fraction
    #   - one_to_one_mapping_fraction
    fig09_stats_glmm_binomial_formula: str = (
        "cbind(success, failure) ~ method + (1|community_id) + (1|sample_key)"
    )
    fig09_stats_glmm_binomial_null_formula: str = "cbind(success, failure) ~ 1 + (1|community_id) + (1|sample_key)"
    # Count panels:
    #   - split_index
    #   - merge_deficit
    fig09_stats_glmm_count_formula: str = "count ~ method + offset(log_offset) + (1|community_id) + (1|sample_key)"
    fig09_stats_glmm_count_null_formula: str = "count ~ 1 + offset(log_offset) + (1|community_id) + (1|sample_key)"
    fig09_stats_glmm_count_family: str = "nbinom2"  # nbinom2 | nbinom1 | poisson
    fig09_stats_glmm_emmeans_scale: str = "link"  # link | response
    # If true, apply BH multiple-test correction across pairwise comparisons
    # within each Figure 09 panel.
    fig09_stats_multiple_test_correction: bool = True
    fig10_top_n_species: int = 30
    fig07_stats_permutations: int = 2000
    fig07_stats_seed: int = 1
    fig07_show_pvalue_label: bool = True
    fig06_grid_rows: Optional[int] = None
    fig06_grid_cols: Optional[int] = None
    fig11_stats_enabled: bool = True
    fig11_stats_permutations: int = 4000
    fig11_stats_seed: int = 1
    fig11_stats_alpha: float = 0.05
    fig11_stats_show_pairwise: bool = True
    # Figure 11 stats engine:
    # - mixed_effect_glmmtmb_auto: beta-GLMM where suitable, fallback to paired permutation
    # - mixed_effect_glmmtmb_beta: beta-GLMM for all Figure 11 metrics
    # - paired_permutation: Friedman + paired permutation/sign tests
    fig11_stats_model: str = "mixed_effect_glmmtmb_auto"  # mixed_effect_glmmtmb_auto | mixed_effect_glmmtmb_beta | paired_permutation
    # Pairwise test in Figure 11:
    # - signed_rank_permutation
    # - sign_test_exact
    fig11_stats_pairwise_test: str = "sign_test_exact"  # signed_rank_permutation | sign_test_exact
    fig11_stats_multiple_test_correction: bool = True
    fig11_stats_glmm_emmeans_scale: str = "link"  # link | response
    fig21_stats_enabled: bool = True
    fig21_stats_permutations: int = 4000
    fig21_stats_seed: int = 1
    fig21_stats_alpha: float = 0.05
    fig21_stats_show_pairwise: bool = True
    # Figure 21 sensitivity can be negative, so it always uses paired permutation/sign-test stats.
    fig21_stats_pairwise_test: str = "sign_test_exact"  # signed_rank_permutation | sign_test_exact | mean_difference_permutation
    fig21_stats_multiple_test_correction: bool = True
    # Exclude communities from figure generation by community_id pattern.
    # Pattern matching supports shell-style wildcards (fnmatch), and plain-text
    # patterns are also treated as substring filters.
    exclude_community_name_patterns: list[str] = field(default_factory=list)
    # Additional figure-specific exclusions by community_id pattern.
    # Effective exclusions are global + per-figure.
    exclude_community_name_patterns_fig01: list[str] = field(default_factory=list)
    exclude_community_name_patterns_fig02: list[str] = field(default_factory=list)
    exclude_community_name_patterns_fig03: list[str] = field(default_factory=list)
    exclude_community_name_patterns_fig04: list[str] = field(default_factory=list)
    exclude_community_name_patterns_fig05: list[str] = field(default_factory=list)
    exclude_community_name_patterns_fig06: list[str] = field(default_factory=list)
    exclude_community_name_patterns_fig07: list[str] = field(default_factory=list)
    exclude_community_name_patterns_fig08: list[str] = field(default_factory=list)
    exclude_community_name_patterns_fig09: list[str] = field(default_factory=list)
    exclude_community_name_patterns_fig10: list[str] = field(default_factory=list)
    exclude_community_name_patterns_fig11: list[str] = field(default_factory=list)
    exclude_community_name_patterns_fig12: list[str] = field(default_factory=list)
    exclude_community_name_patterns_fig13: list[str] = field(default_factory=list)
    exclude_community_name_patterns_fig14: list[str] = field(default_factory=list)
    exclude_community_name_patterns_fig15: list[str] = field(default_factory=list)
    exclude_community_name_patterns_fig16: list[str] = field(default_factory=list)
    exclude_community_name_patterns_fig17: list[str] = field(default_factory=list)
    exclude_community_name_patterns_fig18: list[str] = field(default_factory=list)
    exclude_community_name_patterns_fig19: list[str] = field(default_factory=list)
    exclude_community_name_patterns_fig20: list[str] = field(default_factory=list)
    exclude_community_name_patterns_fig21: list[str] = field(default_factory=list)


@dataclass
class Config:
    paths: PathsConfig
    discovery: DiscoveryConfig = field(default_factory=DiscoveryConfig)
    normalization: NormalizationConfig = field(default_factory=NormalizationConfig)
    mapping: MappingConfig = field(default_factory=MappingConfig)
    tree: TreeConfig = field(default_factory=TreeConfig)
    distances: DistancesConfig = field(default_factory=DistancesConfig)
    figures: FiguresConfig = field(default_factory=FiguresConfig)


def _normalize_figures_data(raw: Dict[str, Any]) -> Dict[str, Any]:
    """
    Accept both legacy flat figure config and nested subsection style:

      figures:
        fig01:
          enabled: true
          grid_rows: 5
          grid_cols: 6
          exclude_community_name_patterns: []
        fig02:
          enabled: true
          exclude_community_name_patterns: []
        fig03:
          by_sample: true
          by_metric: true
          grid_rows: 5
          grid_cols: 6
          by_metric_grid_rows: 5
          by_metric_grid_cols: 6
          exclude_community_name_patterns: []
        fig04:
          enabled: true
          title: "Species Richness Recovery"
          exclude_community_name_patterns: []
        fig05:
          enabled: true
          title: "Precision-Recall (Species Presence)"
          exclude_community_name_patterns: []
        fig06:
          enabled: true
          title: "Rank-Abundance Curves (Species)"
          grid_rows: 5
          grid_cols: 6
          exclude_community_name_patterns: []
        fig07:
          enabled: true
          title: "Alpha Diversity Error Plots"
          plot_type: boxplot
          exclude_community_name_patterns: []
        fig08:
          enabled: true
          title: "Beta-Diversity Concordance"
          plot_type: boxplot
          exclude_community_name_patterns: []
        fig09:
          enabled: true
          title: "Split/Merge Diagnostics"
          split_index_denominator: mapped_expected_features
          merge_denominator: mapped_inferred_features
          plot_type: boxplot
          split_index_mode: mapping
          graph_min_effective_pident_denoising: 100.0
          graph_min_effective_pident_tic_tac: 98.7
          graph_min_effective_pident_otu: 97.0
          min_query_coverage: 0.0
          no_hit_min_query_coverage: null
          no_hit_min_effective_pident_by_method:
            asv: 100.0
            zotu: 100.0
            otu: 97.0
            tic: 98.7
            tac: 98.7
          graph_edge_threshold_by_metric:
            split_index:
              asv: 100.0
              zotu: 100.0
              otu: 97.0
              tic: 98.7
              tac: 98.7
          exclude_community_name_patterns: []
        fig10:
          enabled: true
          title: "Species Confusion Heatmaps"
          exclude_community_name_patterns: []
        fig11:
          enabled: true
          title: "Distance Distribution Violin Plots"
          plot_type: violin
          exclude_community_name_patterns: []
        fig12:
          enabled: true
          title: "Unmapped Feature Diagnostics"
          font_size: 10
          unmapped_best_hit_min_query_coverage: 0.0
          no_hit_min_effective_pident_by_method:
            asv: 100.0
            zotu: 100.0
            otu: 97.0
            tic: 98.7
            tac: 98.7
          exclude_community_name_patterns: []
        fig13:
          enabled: true
          title: "Split Capture by Clustering"
          font_size: 10
          exclude_community_name_patterns: []
        fig14:
          enabled: true
          title: "Merge Resolution Audit"
          font_size: 10
          exclude_community_name_patterns: []
        fig15:
          enabled: true
          title: "Read-Weighted Mapping Outcome"
          font_size: 10
          no_hit_min_effective_pident_by_method:
            asv: 97.0
            zotu: 97.0
            otu: 97.0
            tic: 97.0
            tac: 97.0
          exclude_community_name_patterns: []
        fig16:
          enabled: true
          title: "Species Driver Plot by Community"
          font_size: 10
          species_driver_top_n: 25
          species_driver_sample_aggregation: mean
          species_driver_title_prefix: "Species Driver Plot"
          species_driver_xlabel: "Relative Abundance Deviation (Inferred - Expected)"
          species_driver_ylabel_template: "Top {n} Species by |Deviation|"
          xtick_rotation: 45
          ytick_rotation: 45
          exclude_community_name_patterns: []
        # Global marker size for overlaid jitter points on box/violin plots
        box_violin_dot_size: 22
    """
    fig = dict(raw or {})

    def _optional_float(value):
        if value is None:
            return None
        return float(value)

    s1 = fig.get("fig01")
    if isinstance(s1, dict):
        if "enabled" in s1:
            fig["fig01"] = bool(s1["enabled"])
        else:
            fig["fig01"] = True
        if "grid_rows" in s1:
            fig["grid_rows"] = s1["grid_rows"]
        if "grid_cols" in s1:
            fig["grid_cols"] = s1["grid_cols"]
        if "exclude_community_name_patterns" in s1:
            fig["exclude_community_name_patterns_fig01"] = s1["exclude_community_name_patterns"]
        if "title" in s1:
            fig["fig01_title"] = s1["title"]
        if "title_with_expected_features" in s1:
            fig["fig01_title_with_expected_features"] = s1["title_with_expected_features"]

    s2 = fig.get("fig02")
    if isinstance(s2, dict):
        if "enabled" in s2:
            fig["fig02"] = bool(s2["enabled"])
        else:
            fig["fig02"] = True
        if "exclude_community_name_patterns" in s2:
            fig["exclude_community_name_patterns_fig02"] = s2["exclude_community_name_patterns"]
        if "title" in s2:
            fig["fig02_title"] = s2["title"]

    s3 = fig.pop("fig03", None)
    if isinstance(s3, dict):
        if "enabled" in s3:
            en = bool(s3["enabled"])
            fig["fig03_by_sample"] = en
            fig["fig03_by_metric"] = en
        if "by_sample" in s3:
            fig["fig03_by_sample"] = bool(s3["by_sample"])
        if "by_metric" in s3:
            fig["fig03_by_metric"] = bool(s3["by_metric"])
        if "grid_rows" in s3:
            fig["grid_rows"] = s3["grid_rows"]
        if "grid_cols" in s3:
            fig["grid_cols"] = s3["grid_cols"]
        if "by_metric_grid_rows" in s3:
            fig["fig03_by_metric_grid_rows"] = s3["by_metric_grid_rows"]
        if "by_metric_grid_cols" in s3:
            fig["fig03_by_metric_grid_cols"] = s3["by_metric_grid_cols"]
        if "exclude_community_name_patterns" in s3:
            fig["exclude_community_name_patterns_fig03"] = s3["exclude_community_name_patterns"]
        if "title_by_sample" in s3:
            fig["fig03_by_sample_title"] = s3["title_by_sample"]
        if "title_by_metric" in s3:
            fig["fig03_by_metric_title"] = s3["title_by_metric"]

    s4 = fig.get("fig04")
    if isinstance(s4, dict):
        if "enabled" in s4:
            fig["fig04"] = bool(s4["enabled"])
        else:
            fig["fig04"] = True
        if "title" in s4:
            fig["fig04_title"] = s4["title"]
        if "exclude_community_name_patterns" in s4:
            fig["exclude_community_name_patterns_fig04"] = s4["exclude_community_name_patterns"]

    s5 = fig.get("fig05")
    if isinstance(s5, dict):
        if "enabled" in s5:
            fig["fig05"] = bool(s5["enabled"])
        else:
            fig["fig05"] = True
        if "title" in s5:
            fig["fig05_title"] = s5["title"]
        if "exclude_community_name_patterns" in s5:
            fig["exclude_community_name_patterns_fig05"] = s5["exclude_community_name_patterns"]

    s6 = fig.get("fig06")
    if isinstance(s6, dict):
        if "enabled" in s6:
            fig["fig06"] = bool(s6["enabled"])
        else:
            fig["fig06"] = True
        if "title" in s6:
            fig["fig06_title"] = s6["title"]
        if "grid_rows" in s6:
            fig["fig06_grid_rows"] = s6["grid_rows"]
        if "grid_cols" in s6:
            fig["fig06_grid_cols"] = s6["grid_cols"]
        if "exclude_community_name_patterns" in s6:
            fig["exclude_community_name_patterns_fig06"] = s6["exclude_community_name_patterns"]

    s7 = fig.get("fig07")
    if isinstance(s7, dict):
        if "enabled" in s7:
            fig["fig07"] = bool(s7["enabled"])
        else:
            fig["fig07"] = True
        if "title" in s7:
            fig["fig07_title"] = s7["title"]
        if "exclude_community_name_patterns" in s7:
            fig["exclude_community_name_patterns_fig07"] = s7["exclude_community_name_patterns"]
        if "stats_permutations" in s7:
            fig["fig07_stats_permutations"] = int(s7["stats_permutations"])
        if "stats_seed" in s7:
            fig["fig07_stats_seed"] = int(s7["stats_seed"])
        if "show_pvalue_label" in s7:
            fig["fig07_show_pvalue_label"] = bool(s7["show_pvalue_label"])
        if "plot_type" in s7:
            fig["fig07_plot_type"] = str(s7["plot_type"])
        if "target_level_by_method" in s7 and isinstance(s7["target_level_by_method"], dict):
            fig["fig07_target_level_by_method"] = s7["target_level_by_method"]
        if "stats_model" in s7:
            fig["fig07_stats_model"] = str(s7["stats_model"])
        if "stats_pairwise_test" in s7:
            fig["fig07_stats_pairwise_test"] = str(s7["stats_pairwise_test"])
        if "stats_glmm_emmeans_scale" in s7:
            fig["fig07_stats_glmm_emmeans_scale"] = str(s7["stats_glmm_emmeans_scale"])
        if "stats_alpha" in s7:
            fig["fig07_stats_alpha"] = float(s7["stats_alpha"])
        if "stats_multiple_test_correction" in s7:
            fig["fig07_stats_multiple_test_correction"] = bool(s7["stats_multiple_test_correction"])
        if "stats_show_pairwise" in s7:
            fig["fig07_stats_show_pairwise"] = bool(s7["stats_show_pairwise"])
        if "include_richness_effective_variant" in s7:
            fig["fig07_include_richness_effective_variant"] = bool(s7["include_richness_effective_variant"])
        if "write_long_tsv" in s7:
            fig["fig07_write_long_tsv"] = bool(s7["write_long_tsv"])

    s8 = fig.get("fig08")
    if isinstance(s8, dict):
        if "enabled" in s8:
            fig["fig08"] = bool(s8["enabled"])
        else:
            fig["fig08"] = True
        if "title" in s8:
            fig["fig08_title"] = s8["title"]
        if "exclude_community_name_patterns" in s8:
            fig["exclude_community_name_patterns_fig08"] = s8["exclude_community_name_patterns"]
        if "plot_type" in s8:
            fig["fig08_plot_type"] = str(s8["plot_type"])

    s9 = fig.get("fig09")
    if isinstance(s9, dict):
        if "enabled" in s9:
            fig["fig09"] = bool(s9["enabled"])
        else:
            fig["fig09"] = True
        if "title" in s9:
            fig["fig09_title"] = s9["title"]
        if "exclude_community_name_patterns" in s9:
            fig["exclude_community_name_patterns_fig09"] = s9["exclude_community_name_patterns"]
        if "plot_type" in s9:
            fig["fig09_plot_type"] = str(s9["plot_type"])
        if "font_size" in s9:
            fig["fig09_font_size"] = float(s9["font_size"])
        if "layout_w_pad" in s9:
            fig["fig09_layout_w_pad"] = float(s9["layout_w_pad"])
        if "layout_left" in s9:
            fig["fig09_layout_left"] = float(s9["layout_left"])
        if "layout_bottom" in s9:
            fig["fig09_layout_bottom"] = float(s9["layout_bottom"])
        if "layout_right" in s9:
            fig["fig09_layout_right"] = float(s9["layout_right"])
        if "layout_top" in s9:
            fig["fig09_layout_top"] = float(s9["layout_top"])
        if "save_bbox_tight" in s9:
            fig["fig09_save_bbox_tight"] = bool(s9["save_bbox_tight"])
        if "save_pad_inches" in s9:
            fig["fig09_save_pad_inches"] = float(s9["save_pad_inches"])
        if "write_individual_panels" in s9:
            fig["fig09_write_individual_panels"] = bool(s9["write_individual_panels"])
        if "split_index_denominator" in s9:
            fig["fig09_split_index_denominator"] = str(s9["split_index_denominator"])
        if "merge_denominator" in s9:
            fig["fig09_merge_denominator"] = str(s9["merge_denominator"])
        if "one_to_one_denominator" in s9:
            fig["fig09_one_to_one_denominator"] = str(s9["one_to_one_denominator"])
        if "split_index_mode" in s9:
            fig["fig09_split_index_mode"] = str(s9["split_index_mode"])
        if "graph_target_level_by_metric" in s9 and isinstance(s9["graph_target_level_by_metric"], dict):
            fig["fig09_graph_target_level_by_metric"] = s9["graph_target_level_by_metric"]
        if "panel_target_level_by_metric" in s9 and isinstance(s9["panel_target_level_by_metric"], dict):
            fig["fig09_panel_target_level_by_metric"] = s9["panel_target_level_by_metric"]
        if "graph_edge_threshold_by_metric" in s9 and isinstance(s9["graph_edge_threshold_by_metric"], dict):
            fig["fig09_graph_edge_threshold_by_metric"] = s9["graph_edge_threshold_by_metric"]
        if "panel_c_enabled" in s9:
            fig["fig09_panel_c_enabled"] = bool(s9["panel_c_enabled"])
        if "panel_titles" in s9 and isinstance(s9["panel_titles"], dict):
            fig["fig09_panel_titles"] = s9["panel_titles"]
        if "graph_labels_by_metric" in s9 and isinstance(s9["graph_labels_by_metric"], dict):
            fig["fig09_graph_labels_by_metric"] = s9["graph_labels_by_metric"]
        s9_panel_d = s9.get("panel_d")
        if isinstance(s9_panel_d, dict):
            if "title" in s9_panel_d:
                fig["fig09_panel_d_title"] = None if s9_panel_d["title"] is None else str(s9_panel_d["title"])
            if "xlabel" in s9_panel_d:
                fig["fig09_panel_d_xlabel"] = None if s9_panel_d["xlabel"] is None else str(s9_panel_d["xlabel"])
            if "x_label" in s9_panel_d:
                fig["fig09_panel_d_xlabel"] = None if s9_panel_d["x_label"] is None else str(s9_panel_d["x_label"])
            if "ylabel" in s9_panel_d:
                fig["fig09_panel_d_ylabel"] = None if s9_panel_d["ylabel"] is None else str(s9_panel_d["ylabel"])
            if "y_label" in s9_panel_d:
                fig["fig09_panel_d_ylabel"] = None if s9_panel_d["y_label"] is None else str(s9_panel_d["y_label"])
            if "methods" in s9_panel_d and isinstance(s9_panel_d["methods"], list):
                fig["fig09_panel_d_methods"] = [str(v) for v in s9_panel_d["methods"]]
            if "plot_type" in s9_panel_d:
                fig["fig09_panel_d_plot_type"] = str(s9_panel_d["plot_type"])
            if "font_size" in s9_panel_d:
                v = s9_panel_d["font_size"]
                fig["fig09_panel_d_font_size"] = None if v is None else float(v)
            if "tick_size" in s9_panel_d:
                v = s9_panel_d["tick_size"]
                fig["fig09_panel_d_tick_size"] = None if v is None else float(v)
            if "xtick_rotation" in s9_panel_d:
                fig["fig09_panel_d_xtick_rotation"] = float(s9_panel_d["xtick_rotation"])
            if "x_tick_rotation" in s9_panel_d:
                fig["fig09_panel_d_xtick_rotation"] = float(s9_panel_d["x_tick_rotation"])
            if "ytick_rotation" in s9_panel_d:
                fig["fig09_panel_d_ytick_rotation"] = float(s9_panel_d["ytick_rotation"])
            if "y_tick_rotation" in s9_panel_d:
                fig["fig09_panel_d_ytick_rotation"] = float(s9_panel_d["y_tick_rotation"])
            if "point_size" in s9_panel_d:
                v = s9_panel_d["point_size"]
                fig["fig09_panel_d_point_size"] = None if v is None else float(v)
            if "dot_size" in s9_panel_d:
                v = s9_panel_d["dot_size"]
                fig["fig09_panel_d_point_size"] = None if v is None else float(v)
            if "point_opacity" in s9_panel_d:
                fig["fig09_panel_d_point_opacity"] = float(s9_panel_d["point_opacity"])
            if "dot_opacity" in s9_panel_d:
                fig["fig09_panel_d_point_opacity"] = float(s9_panel_d["dot_opacity"])
            if "point_jitter" in s9_panel_d:
                fig["fig09_panel_d_point_jitter"] = float(s9_panel_d["point_jitter"])
            if "dot_jitter" in s9_panel_d:
                fig["fig09_panel_d_point_jitter"] = float(s9_panel_d["dot_jitter"])
            if "show_zero_line" in s9_panel_d:
                fig["fig09_panel_d_show_zero_line"] = bool(s9_panel_d["show_zero_line"])
            if "show_reference_note" in s9_panel_d:
                fig["fig09_panel_d_show_reference_note"] = bool(s9_panel_d["show_reference_note"])
            if "reference_note" in s9_panel_d:
                fig["fig09_panel_d_reference_note"] = str(s9_panel_d["reference_note"] or "")
            if "show_excess_removed_annotation" in s9_panel_d:
                fig["fig09_panel_d_show_excess_removed_annotation"] = bool(
                    s9_panel_d["show_excess_removed_annotation"]
                )
            if "y_limits" in s9_panel_d and isinstance(s9_panel_d["y_limits"], list):
                fig["fig09_panel_d_y_limits"] = [float(v) for v in s9_panel_d["y_limits"][:2]]
        s9_transition = s9.get("transition_heatmap")
        if isinstance(s9_transition, dict):
            if "title" in s9_transition:
                fig["fig09_transition_heatmap_title"] = (
                    None if s9_transition["title"] is None else str(s9_transition["title"])
                )
            if "xlabel" in s9_transition:
                fig["fig09_transition_heatmap_xlabel"] = (
                    None if s9_transition["xlabel"] is None else str(s9_transition["xlabel"])
                )
            if "x_label" in s9_transition:
                fig["fig09_transition_heatmap_xlabel"] = (
                    None if s9_transition["x_label"] is None else str(s9_transition["x_label"])
                )
            if "ylabel" in s9_transition:
                fig["fig09_transition_heatmap_ylabel"] = (
                    None if s9_transition["ylabel"] is None else str(s9_transition["ylabel"])
                )
            if "y_label" in s9_transition:
                fig["fig09_transition_heatmap_ylabel"] = (
                    None if s9_transition["y_label"] is None else str(s9_transition["y_label"])
                )
            if "note" in s9_transition:
                fig["fig09_transition_heatmap_note"] = str(s9_transition["note"] or "")
            if "cmap" in s9_transition:
                fig["fig09_transition_heatmap_cmap"] = str(s9_transition["cmap"])
            if "show_counts" in s9_transition:
                fig["fig09_transition_heatmap_show_counts"] = bool(s9_transition["show_counts"])
            if "show_percentages" in s9_transition:
                fig["fig09_transition_heatmap_show_percentages"] = bool(s9_transition["show_percentages"])
            if "show_colorbar" in s9_transition:
                fig["fig09_transition_heatmap_show_colorbar"] = bool(s9_transition["show_colorbar"])
            if "x_tick_rotation" in s9_transition:
                fig["fig09_transition_heatmap_x_tick_rotation"] = float(s9_transition["x_tick_rotation"])
            if "row_label_rotation" in s9_transition:
                fig["fig09_transition_heatmap_row_label_rotation"] = float(s9_transition["row_label_rotation"])
            if "method_label_x_offset" in s9_transition:
                fig["fig09_transition_heatmap_method_label_x_offset"] = float(
                    s9_transition["method_label_x_offset"]
                )
            if "colorbar_pad" in s9_transition:
                fig["fig09_transition_heatmap_colorbar_pad"] = float(s9_transition["colorbar_pad"])
            if "tic_label" in s9_transition:
                fig["fig09_transition_heatmap_tic_label"] = str(s9_transition["tic_label"])
            if "tac_label" in s9_transition:
                fig["fig09_transition_heatmap_tac_label"] = str(s9_transition["tac_label"])
            if "no_match_label" in s9_transition:
                fig["fig09_transition_heatmap_no_match_label"] = str(s9_transition["no_match_label"])
            if "non_bacterial_label" in s9_transition:
                fig["fig09_transition_heatmap_non_bacterial_label"] = str(
                    s9_transition["non_bacterial_label"]
                )
            if "colorbar_label" in s9_transition:
                fig["fig09_transition_heatmap_colorbar_label"] = str(s9_transition["colorbar_label"])
            if "exact_identity" in s9_transition:
                fig["fig09_transition_heatmap_exact_identity"] = float(s9_transition["exact_identity"])
            if "near_exact_min_effective_pident" in s9_transition:
                v = s9_transition["near_exact_min_effective_pident"]
                fig["fig09_transition_heatmap_near_exact_min_effective_pident"] = (
                    None if v is None else float(v)
                )
        s9_homogeneity = s9.get("homogeneity")
        if isinstance(s9_homogeneity, dict):
            if "enabled" in s9_homogeneity:
                fig["fig09_homogeneity_enabled"] = bool(s9_homogeneity["enabled"])
            if "title" in s9_homogeneity:
                fig["fig09_homogeneity_title"] = str(s9_homogeneity["title"])
            if "font_size" in s9_homogeneity:
                v = s9_homogeneity["font_size"]
                fig["fig09_homogeneity_font_size"] = None if v is None else float(v)
            if "tick_size" in s9_homogeneity:
                v = s9_homogeneity["tick_size"]
                fig["fig09_homogeneity_tick_size"] = None if v is None else float(v)
            if "figsize" in s9_homogeneity and isinstance(s9_homogeneity["figsize"], list):
                fig["fig09_homogeneity_figsize"] = [float(v) for v in s9_homogeneity["figsize"][:2]]
            if "count_ylabel" in s9_homogeneity:
                fig["fig09_homogeneity_count_ylabel"] = str(s9_homogeneity["count_ylabel"])
            if "effect_ylabel" in s9_homogeneity:
                fig["fig09_homogeneity_effect_ylabel"] = str(s9_homogeneity["effect_ylabel"])
            if "count_label" in s9_homogeneity:
                fig["fig09_homogeneity_count_label"] = str(s9_homogeneity["count_label"])
            if "effect_label" in s9_homogeneity:
                fig["fig09_homogeneity_effect_label"] = str(s9_homogeneity["effect_label"])
            if "count_legend_label" in s9_homogeneity:
                fig["fig09_homogeneity_count_legend_label"] = str(s9_homogeneity["count_legend_label"])
            if "single_species_label" in s9_homogeneity:
                fig["fig09_homogeneity_single_species_label"] = str(s9_homogeneity["single_species_label"])
            if "multi_species_label" in s9_homogeneity:
                fig["fig09_homogeneity_multi_species_label"] = str(s9_homogeneity["multi_species_label"])
            if "single_species_color" in s9_homogeneity:
                fig["fig09_homogeneity_single_species_color"] = str(s9_homogeneity["single_species_color"])
            if "multi_species_color" in s9_homogeneity:
                fig["fig09_homogeneity_multi_species_color"] = str(s9_homogeneity["multi_species_color"])
            if "effect_hatch" in s9_homogeneity:
                fig["fig09_homogeneity_effect_hatch"] = str(s9_homogeneity["effect_hatch"])
            if "effect_alpha" in s9_homogeneity:
                fig["fig09_homogeneity_effect_alpha"] = float(s9_homogeneity["effect_alpha"])
            if "effect_axis_max" in s9_homogeneity:
                fig["fig09_homogeneity_effect_axis_max"] = float(s9_homogeneity["effect_axis_max"])
            if "count_axis_padding" in s9_homogeneity:
                fig["fig09_homogeneity_count_axis_padding"] = float(s9_homogeneity["count_axis_padding"])
            if "show_value_labels" in s9_homogeneity:
                fig["fig09_homogeneity_show_value_labels"] = bool(s9_homogeneity["show_value_labels"])
            if "show_bar_type_labels" in s9_homogeneity:
                fig["fig09_homogeneity_show_bar_type_labels"] = bool(s9_homogeneity["show_bar_type_labels"])
            if "show_legends" in s9_homogeneity:
                fig["fig09_homogeneity_show_legends"] = bool(s9_homogeneity["show_legends"])
            titles = s9_homogeneity.get("panel_titles")
            if isinstance(titles, dict):
                if "pure_exact" in titles:
                    fig["fig09_homogeneity_pure_exact_title"] = str(titles["pure_exact"])
                if "pure_near" in titles:
                    fig["fig09_homogeneity_pure_near_title"] = str(titles["pure_near"])
                if "pure_no_match" in titles:
                    fig["fig09_homogeneity_pure_no_match_title"] = str(titles["pure_no_match"])
                if "exact_near" in titles:
                    fig["fig09_homogeneity_exact_near_title"] = str(titles["exact_near"])
            if "save_bbox_tight" in s9_homogeneity:
                fig["fig09_homogeneity_save_bbox_tight"] = bool(s9_homogeneity["save_bbox_tight"])
            if "save_pad_inches" in s9_homogeneity:
                fig["fig09_homogeneity_save_pad_inches"] = float(s9_homogeneity["save_pad_inches"])
        if "graph_min_effective_pident_denoising" in s9:
            fig["fig09_graph_min_effective_pident_denoising"] = float(s9["graph_min_effective_pident_denoising"])
        elif "graph_min_identity_denoising" in s9:
            # Backward-compatible alias
            fig["fig09_graph_min_effective_pident_denoising"] = float(s9["graph_min_identity_denoising"])
        if "graph_min_effective_pident_tic_tac" in s9:
            fig["fig09_graph_min_effective_pident_tic_tac"] = float(s9["graph_min_effective_pident_tic_tac"])
        elif "graph_min_identity_tic_tac" in s9:
            # Backward-compatible alias
            fig["fig09_graph_min_effective_pident_tic_tac"] = float(s9["graph_min_identity_tic_tac"])
        if "graph_min_effective_pident_otu" in s9:
            fig["fig09_graph_min_effective_pident_otu"] = float(s9["graph_min_effective_pident_otu"])
        elif "graph_min_identity_otu" in s9:
            # Backward-compatible alias
            fig["fig09_graph_min_effective_pident_otu"] = float(s9["graph_min_identity_otu"])
        if "min_query_coverage" in s9:
            fig["fig09_min_query_coverage"] = float(s9["min_query_coverage"])
        if "no_hit_min_query_coverage" in s9:
            v = s9["no_hit_min_query_coverage"]
            fig["fig09_no_hit_min_query_coverage"] = None if v is None else float(v)
        if "no_hit_min_effective_pident_by_method" in s9 and isinstance(
            s9["no_hit_min_effective_pident_by_method"], dict
        ):
            fig["fig09_no_hit_min_effective_pident_by_method"] = s9[
                "no_hit_min_effective_pident_by_method"
            ]
        if "unmapped_best_hit_min_query_coverage" in s9:
            v = s9["unmapped_best_hit_min_query_coverage"]
            vv = None if v is None else float(v)
            # Backward compatibility: legacy fig09 key controls Figure 09
            # no-hit/unmapped query-coverage gate.
            fig["fig09_no_hit_min_query_coverage"] = vv
            # Keep historical behavior for Figure 12 as well.
            fig["fig12_unmapped_best_hit_min_query_coverage"] = vv
        s9_stats = s9.get("stats")
        if isinstance(s9_stats, dict):
            if "enabled" in s9_stats:
                fig["fig09_stats_enabled"] = bool(s9_stats["enabled"])
            if "permutations" in s9_stats:
                fig["fig09_stats_permutations"] = int(s9_stats["permutations"])
            if "seed" in s9_stats:
                fig["fig09_stats_seed"] = int(s9_stats["seed"])
            if "alpha" in s9_stats:
                fig["fig09_stats_alpha"] = float(s9_stats["alpha"])
            if "show_pairwise" in s9_stats:
                fig["fig09_stats_show_pairwise"] = bool(s9_stats["show_pairwise"])
            if "model" in s9_stats:
                fig["fig09_stats_model"] = str(s9_stats["model"])
            if "pairwise_test" in s9_stats:
                fig["fig09_stats_pairwise_test"] = str(s9_stats["pairwise_test"])
            if "multiple_test_correction" in s9_stats:
                fig["fig09_stats_multiple_test_correction"] = bool(s9_stats["multiple_test_correction"])
            if "glmm_binomial_formula" in s9_stats:
                fig["fig09_stats_glmm_binomial_formula"] = str(s9_stats["glmm_binomial_formula"])
            if "glmm_binomial_null_formula" in s9_stats:
                fig["fig09_stats_glmm_binomial_null_formula"] = str(s9_stats["glmm_binomial_null_formula"])
            if "glmm_count_formula" in s9_stats:
                fig["fig09_stats_glmm_count_formula"] = str(s9_stats["glmm_count_formula"])
            if "glmm_count_null_formula" in s9_stats:
                fig["fig09_stats_glmm_count_null_formula"] = str(s9_stats["glmm_count_null_formula"])
            if "glmm_count_family" in s9_stats:
                fig["fig09_stats_glmm_count_family"] = str(s9_stats["glmm_count_family"])
            if "glmm_emmeans_scale" in s9_stats:
                fig["fig09_stats_glmm_emmeans_scale"] = str(s9_stats["glmm_emmeans_scale"])
        if "stats_enabled" in s9:
            fig["fig09_stats_enabled"] = bool(s9["stats_enabled"])
        if "stats_permutations" in s9:
            fig["fig09_stats_permutations"] = int(s9["stats_permutations"])
        if "stats_seed" in s9:
            fig["fig09_stats_seed"] = int(s9["stats_seed"])
        if "stats_alpha" in s9:
            fig["fig09_stats_alpha"] = float(s9["stats_alpha"])
        if "stats_show_pairwise" in s9:
            fig["fig09_stats_show_pairwise"] = bool(s9["stats_show_pairwise"])
        if "stats_model" in s9:
            fig["fig09_stats_model"] = str(s9["stats_model"])
        if "stats_pairwise_test" in s9:
            fig["fig09_stats_pairwise_test"] = str(s9["stats_pairwise_test"])
        if "stats_multiple_test_correction" in s9:
            fig["fig09_stats_multiple_test_correction"] = bool(s9["stats_multiple_test_correction"])
        if "stats_glmm_binomial_formula" in s9:
            fig["fig09_stats_glmm_binomial_formula"] = str(s9["stats_glmm_binomial_formula"])
        if "stats_glmm_binomial_null_formula" in s9:
            fig["fig09_stats_glmm_binomial_null_formula"] = str(s9["stats_glmm_binomial_null_formula"])
        if "stats_glmm_count_formula" in s9:
            fig["fig09_stats_glmm_count_formula"] = str(s9["stats_glmm_count_formula"])
        if "stats_glmm_count_null_formula" in s9:
            fig["fig09_stats_glmm_count_null_formula"] = str(s9["stats_glmm_count_null_formula"])
        if "stats_glmm_count_family" in s9:
            fig["fig09_stats_glmm_count_family"] = str(s9["stats_glmm_count_family"])
        if "stats_glmm_emmeans_scale" in s9:
            fig["fig09_stats_glmm_emmeans_scale"] = str(s9["stats_glmm_emmeans_scale"])

    s10 = fig.get("fig10")
    if isinstance(s10, dict):
        if "enabled" in s10:
            fig["fig10"] = bool(s10["enabled"])
        else:
            fig["fig10"] = True
        if "title" in s10:
            fig["fig10_title"] = s10["title"]
        if "top_n_species" in s10:
            fig["fig10_top_n_species"] = int(s10["top_n_species"])
        if "exclude_community_name_patterns" in s10:
            fig["exclude_community_name_patterns_fig10"] = s10["exclude_community_name_patterns"]

    s11 = fig.get("fig11")
    if isinstance(s11, dict):
        if "enabled" in s11:
            fig["fig11"] = bool(s11["enabled"])
        else:
            fig["fig11"] = True
        if "title" in s11:
            fig["fig11_title"] = s11["title"]
        if "metric" in s11 and str(s11["metric"]).strip():
            fig["fig11_metrics"] = [str(s11["metric"])]
        elif "metrics" in s11 and isinstance(s11["metrics"], list):
            fig["fig11_metrics"] = [str(v) for v in s11["metrics"] if str(v).strip()]
        if "panel_title_by_metric" in s11 and isinstance(s11["panel_title_by_metric"], dict):
            fig["fig11_panel_title_by_metric"] = s11["panel_title_by_metric"]
        if "combined_2x2_title" in s11:
            fig["fig11_combined_2x2_title"] = s11["combined_2x2_title"]
        if "combined_3x2_title" in s11:
            fig["fig11_combined_3x2_title"] = s11["combined_3x2_title"]
        if "exclude_community_name_patterns" in s11:
            fig["exclude_community_name_patterns_fig11"] = s11["exclude_community_name_patterns"]
        if "plot_type" in s11:
            fig["fig11_plot_type"] = str(s11["plot_type"])
        if "font_size" in s11:
            fig["fig11_font_size"] = float(s11["font_size"])
        if "projection_mapping_mode" in s11:
            fig["fig11_projection_mapping_mode"] = str(s11["projection_mapping_mode"])
        if "target_level_by_method" in s11 and isinstance(s11["target_level_by_method"], dict):
            fig["fig11_target_level_by_method"] = s11["target_level_by_method"]
        if "min_effective_pident_by_method" in s11 and isinstance(s11["min_effective_pident_by_method"], dict):
            fig["fig11_min_effective_pident_by_method"] = s11["min_effective_pident_by_method"]
        s11_stats = s11.get("stats")
        if isinstance(s11_stats, dict):
            if "enabled" in s11_stats:
                fig["fig11_stats_enabled"] = bool(s11_stats["enabled"])
            if "permutations" in s11_stats:
                fig["fig11_stats_permutations"] = int(s11_stats["permutations"])
            if "seed" in s11_stats:
                fig["fig11_stats_seed"] = int(s11_stats["seed"])
            if "alpha" in s11_stats:
                fig["fig11_stats_alpha"] = float(s11_stats["alpha"])
            if "show_pairwise" in s11_stats:
                fig["fig11_stats_show_pairwise"] = bool(s11_stats["show_pairwise"])
            if "model" in s11_stats:
                fig["fig11_stats_model"] = str(s11_stats["model"])
            if "pairwise_test" in s11_stats:
                fig["fig11_stats_pairwise_test"] = str(s11_stats["pairwise_test"])
            if "multiple_test_correction" in s11_stats:
                fig["fig11_stats_multiple_test_correction"] = bool(s11_stats["multiple_test_correction"])
            if "glmm_emmeans_scale" in s11_stats:
                fig["fig11_stats_glmm_emmeans_scale"] = str(s11_stats["glmm_emmeans_scale"])
        if "stats_enabled" in s11:
            fig["fig11_stats_enabled"] = bool(s11["stats_enabled"])
        if "stats_permutations" in s11:
            fig["fig11_stats_permutations"] = int(s11["stats_permutations"])
        if "stats_seed" in s11:
            fig["fig11_stats_seed"] = int(s11["stats_seed"])
        if "stats_alpha" in s11:
            fig["fig11_stats_alpha"] = float(s11["stats_alpha"])
        if "stats_show_pairwise" in s11:
            fig["fig11_stats_show_pairwise"] = bool(s11["stats_show_pairwise"])
        if "stats_model" in s11:
            fig["fig11_stats_model"] = str(s11["stats_model"])
        if "stats_pairwise_test" in s11:
            fig["fig11_stats_pairwise_test"] = str(s11["stats_pairwise_test"])
        if "stats_multiple_test_correction" in s11:
            fig["fig11_stats_multiple_test_correction"] = bool(s11["stats_multiple_test_correction"])
        if "stats_glmm_emmeans_scale" in s11:
            fig["fig11_stats_glmm_emmeans_scale"] = str(s11["stats_glmm_emmeans_scale"])

    s12 = fig.get("fig12")
    if isinstance(s12, dict):
        if "enabled" in s12:
            fig["fig12"] = bool(s12["enabled"])
        else:
            fig["fig12"] = True
        if "title" in s12:
            fig["fig12_title"] = s12["title"]
        if "font_size" in s12:
            fig["fig12_font_size"] = float(s12["font_size"])
        if "unmapped_best_hit_min_query_coverage" in s12:
            v = s12["unmapped_best_hit_min_query_coverage"]
            fig["fig12_unmapped_best_hit_min_query_coverage"] = None if v is None else float(v)
        if "no_hit_min_effective_pident_by_method" in s12 and isinstance(
            s12["no_hit_min_effective_pident_by_method"], dict
        ):
            fig["fig12_no_hit_min_effective_pident_by_method"] = s12[
                "no_hit_min_effective_pident_by_method"
            ]
        if "exclude_community_name_patterns" in s12:
            fig["exclude_community_name_patterns_fig12"] = s12["exclude_community_name_patterns"]

    s13 = fig.get("fig13")
    if isinstance(s13, dict):
        if "enabled" in s13:
            fig["fig13"] = bool(s13["enabled"])
        else:
            fig["fig13"] = True
        if "title" in s13:
            fig["fig13_title"] = s13["title"]
        if "font_size" in s13:
            fig["fig13_font_size"] = float(s13["font_size"])
        if "exclude_community_name_patterns" in s13:
            fig["exclude_community_name_patterns_fig13"] = s13["exclude_community_name_patterns"]

    s14 = fig.get("fig14")
    if isinstance(s14, dict):
        if "enabled" in s14:
            fig["fig14"] = bool(s14["enabled"])
        else:
            fig["fig14"] = True
        if "title" in s14:
            fig["fig14_title"] = s14["title"]
        if "font_size" in s14:
            fig["fig14_font_size"] = float(s14["font_size"])
        if "show_missing_similarity_data_category" in s14:
            fig["fig14_show_missing_similarity_data_category"] = bool(
                s14["show_missing_similarity_data_category"]
            )
        elif "show_missing_similarity_category" in s14:
            fig["fig14_show_missing_similarity_data_category"] = bool(
                s14["show_missing_similarity_category"]
            )
        if "graph_min_effective_pident_denoising" in s14:
            fig["fig14_graph_min_effective_pident_denoising"] = float(s14["graph_min_effective_pident_denoising"])
        elif "graph_min_identity_denoising" in s14:
            fig["fig14_graph_min_effective_pident_denoising"] = float(s14["graph_min_identity_denoising"])
        if "graph_min_effective_pident_tic_tac" in s14:
            fig["fig14_graph_min_effective_pident_tic_tac"] = float(s14["graph_min_effective_pident_tic_tac"])
        elif "graph_min_identity_tic_tac" in s14:
            fig["fig14_graph_min_effective_pident_tic_tac"] = float(s14["graph_min_identity_tic_tac"])
        if "graph_min_effective_pident_otu" in s14:
            fig["fig14_graph_min_effective_pident_otu"] = float(s14["graph_min_effective_pident_otu"])
        elif "graph_min_identity_otu" in s14:
            fig["fig14_graph_min_effective_pident_otu"] = float(s14["graph_min_identity_otu"])
        if "edge_threshold_by_method" in s14 and isinstance(s14["edge_threshold_by_method"], dict):
            fig["fig14_edge_threshold_by_method"] = s14["edge_threshold_by_method"]
        elif "graph_edge_threshold_by_method" in s14 and isinstance(s14["graph_edge_threshold_by_method"], dict):
            fig["fig14_edge_threshold_by_method"] = s14["graph_edge_threshold_by_method"]
        if "min_query_coverage" in s14:
            v = s14["min_query_coverage"]
            fig["fig14_min_query_coverage"] = None if v is None else float(v)
        if "min_query_coverage_by_method" in s14 and isinstance(s14["min_query_coverage_by_method"], dict):
            fig["fig14_min_query_coverage_by_method"] = s14["min_query_coverage_by_method"]
        elif "graph_min_query_coverage_by_method" in s14 and isinstance(s14["graph_min_query_coverage_by_method"], dict):
            fig["fig14_min_query_coverage_by_method"] = s14["graph_min_query_coverage_by_method"]
        # Legacy Figure-14 location for species-driver settings (moved to fig16).
        if "species_driver_top_n" in s14:
            fig["fig16_species_driver_top_n"] = int(s14["species_driver_top_n"])
        elif "driver_top_n" in s14:
            fig["fig16_species_driver_top_n"] = int(s14["driver_top_n"])
        if "species_driver_sample_aggregation" in s14:
            fig["fig16_species_driver_sample_aggregation"] = str(s14["species_driver_sample_aggregation"])
        elif "driver_sample_aggregation" in s14:
            fig["fig16_species_driver_sample_aggregation"] = str(s14["driver_sample_aggregation"])
        if "species_driver_title_prefix" in s14:
            fig["fig16_species_driver_title_prefix"] = str(s14["species_driver_title_prefix"])
        elif "driver_title_prefix" in s14:
            fig["fig16_species_driver_title_prefix"] = str(s14["driver_title_prefix"])
        if "species_driver_xlabel" in s14:
            fig["fig16_species_driver_xlabel"] = str(s14["species_driver_xlabel"])
        elif "driver_xlabel" in s14:
            fig["fig16_species_driver_xlabel"] = str(s14["driver_xlabel"])
        if "species_driver_ylabel_template" in s14:
            fig["fig16_species_driver_ylabel_template"] = str(s14["species_driver_ylabel_template"])
        elif "driver_ylabel_template" in s14:
            fig["fig16_species_driver_ylabel_template"] = str(s14["driver_ylabel_template"])
        if "xtick_rotation" in s14:
            fig["fig16_xtick_rotation"] = float(s14["xtick_rotation"])
        if "ytick_rotation" in s14:
            fig["fig16_ytick_rotation"] = float(s14["ytick_rotation"])
        if "projection_abundance_mode" in s14:
            fig["fig16_projection_abundance_mode"] = str(s14["projection_abundance_mode"])
        elif "driver_projection_abundance_mode" in s14:
            fig["fig16_projection_abundance_mode"] = str(s14["driver_projection_abundance_mode"])
        if "exclude_community_name_patterns" in s14:
            fig["exclude_community_name_patterns_fig14"] = s14["exclude_community_name_patterns"]

    s15 = fig.get("fig15")
    if isinstance(s15, dict):
        if "enabled" in s15:
            fig["fig15"] = bool(s15["enabled"])
        else:
            fig["fig15"] = True
        if "title" in s15:
            fig["fig15_title"] = s15["title"]
        if "font_size" in s15:
            fig["fig15_font_size"] = float(s15["font_size"])
        if "no_hit_min_effective_pident_by_method" in s15 and isinstance(
            s15["no_hit_min_effective_pident_by_method"], dict
        ):
            fig["fig15_no_hit_min_effective_pident_by_method"] = s15[
                "no_hit_min_effective_pident_by_method"
            ]
        if "exclude_community_name_patterns" in s15:
            fig["exclude_community_name_patterns_fig15"] = s15["exclude_community_name_patterns"]

    s16 = fig.get("fig16")
    if isinstance(s16, dict):
        if "enabled" in s16:
            fig["fig16"] = bool(s16["enabled"])
        else:
            fig["fig16"] = True
        if "title" in s16:
            fig["fig16_title"] = s16["title"]
        if "font_size" in s16:
            fig["fig16_font_size"] = float(s16["font_size"])
        if "methods" in s16:
            raw_methods = s16["methods"]
            if isinstance(raw_methods, str):
                fig["fig16_methods"] = [
                    str(v).strip()
                    for v in re.split(r"[,;]", raw_methods)
                    if str(v).strip()
                ]
            elif isinstance(raw_methods, list):
                fig["fig16_methods"] = [str(v).strip() for v in raw_methods if str(v).strip()]
        if "species_driver_top_n" in s16:
            fig["fig16_species_driver_top_n"] = int(s16["species_driver_top_n"])
        elif "driver_top_n" in s16:
            fig["fig16_species_driver_top_n"] = int(s16["driver_top_n"])
        if "species_driver_sample_aggregation" in s16:
            fig["fig16_species_driver_sample_aggregation"] = str(s16["species_driver_sample_aggregation"])
        elif "driver_sample_aggregation" in s16:
            fig["fig16_species_driver_sample_aggregation"] = str(s16["driver_sample_aggregation"])
        if "species_driver_title_prefix" in s16:
            fig["fig16_species_driver_title_prefix"] = str(s16["species_driver_title_prefix"])
        elif "driver_title_prefix" in s16:
            fig["fig16_species_driver_title_prefix"] = str(s16["driver_title_prefix"])
        if "species_driver_xlabel" in s16:
            fig["fig16_species_driver_xlabel"] = str(s16["species_driver_xlabel"])
        elif "driver_xlabel" in s16:
            fig["fig16_species_driver_xlabel"] = str(s16["driver_xlabel"])
        if "species_driver_ylabel_template" in s16:
            fig["fig16_species_driver_ylabel_template"] = str(s16["species_driver_ylabel_template"])
        elif "driver_ylabel_template" in s16:
            fig["fig16_species_driver_ylabel_template"] = str(s16["driver_ylabel_template"])
        if "xtick_rotation" in s16:
            fig["fig16_xtick_rotation"] = float(s16["xtick_rotation"])
        if "ytick_rotation" in s16:
            fig["fig16_ytick_rotation"] = float(s16["ytick_rotation"])
        if "species_driver_xlim" in s16:
            fig["fig16_species_driver_xlim"] = s16["species_driver_xlim"]
        elif "driver_xlim" in s16:
            fig["fig16_species_driver_xlim"] = s16["driver_xlim"]
        if "projection_abundance_mode" in s16:
            fig["fig16_projection_abundance_mode"] = str(s16["projection_abundance_mode"])
        elif "driver_projection_abundance_mode" in s16:
            fig["fig16_projection_abundance_mode"] = str(s16["driver_projection_abundance_mode"])
        if "exclude_community_name_patterns" in s16:
            fig["exclude_community_name_patterns_fig16"] = s16["exclude_community_name_patterns"]

    s17 = fig.get("fig17")
    if isinstance(s17, dict):
        if "enabled" in s17:
            fig["fig17"] = bool(s17["enabled"])
        else:
            fig["fig17"] = True
        if "title" in s17:
            fig["fig17_title"] = s17["title"]
        if "font_size" in s17:
            fig["fig17_font_size"] = float(s17["font_size"])
        if "dot_size" in s17:
            fig["fig17_dot_size"] = float(s17["dot_size"])
        if "opacity" in s17:
            fig["fig17_opacity"] = float(s17["opacity"])
        elif "density_fill_opacity" in s17:
            fig["fig17_opacity"] = float(s17["density_fill_opacity"])
        if "show_upper_panel" in s17:
            fig["fig17_show_upper_panel"] = bool(s17["show_upper_panel"])
        if "method_view" in s17:
            fig["fig17_method_view"] = str(s17["method_view"])
        elif "collapse_method_groups" in s17:
            fig["fig17_method_view"] = "families" if bool(s17["collapse_method_groups"]) else "methods"
        if "lower_panel_y" in s17:
            fig["fig17_lower_panel_y"] = str(s17["lower_panel_y"])
        if "lower_panel" in s17:
            fig["fig17_lower_panel"] = str(s17["lower_panel"])
        if "grid_mode" in s17:
            fig["fig17_grid_mode"] = str(s17["grid_mode"])
        if "lower_panel_height_ratio" in s17:
            fig["fig17_lower_panel_height_ratio"] = float(s17["lower_panel_height_ratio"])
        if "density_bandwidth" in s17:
            fig["fig17_density_bandwidth"] = (
                None if s17["density_bandwidth"] is None else float(s17["density_bandwidth"])
            )
        if "density_points" in s17:
            fig["fig17_density_points"] = int(s17["density_points"])
        if "density_scale" in s17:
            fig["fig17_density_scale"] = float(s17["density_scale"])
        if "x_tick_interval" in s17:
            fig["fig17_xtick_interval"] = float(s17["x_tick_interval"])
        elif "xtick_interval" in s17:
            fig["fig17_xtick_interval"] = float(s17["xtick_interval"])
        if "metric_title_pad" in s17:
            fig["fig17_metric_title_pad"] = float(s17["metric_title_pad"])
        if "metric" in s17 and s17["metric"] not in (None, ""):
            fig["fig17_metrics"] = [str(s17["metric"])]
        elif "metrics" in s17 and isinstance(s17["metrics"], list):
            fig["fig17_metrics"] = [str(v) for v in s17["metrics"] if str(v).strip()]
        stats17 = s17.get("stats")
        if isinstance(stats17, dict):
            if "enabled" in stats17:
                fig["fig17_stats_enabled"] = bool(stats17["enabled"])
            if "alpha" in stats17:
                fig["fig17_stats_alpha"] = float(stats17["alpha"])
            else:
                fig.setdefault("fig17_stats_alpha", float(fig.get("fig09_stats_alpha", 0.05)))
            if "pairwise_test" in stats17:
                fig["fig17_stats_pairwise_test"] = str(stats17["pairwise_test"])
            else:
                fig.setdefault(
                    "fig17_stats_pairwise_test",
                    str(fig.get("fig09_stats_pairwise_test", "signed_rank_permutation")),
                )
            if "multiple_test_correction" in stats17:
                fig["fig17_stats_multiple_test_correction"] = bool(stats17["multiple_test_correction"])
            else:
                fig.setdefault(
                    "fig17_stats_multiple_test_correction",
                    bool(fig.get("fig09_stats_multiple_test_correction", True)),
                )
            if "family_aggregation" in stats17:
                fig["fig17_stats_family_aggregation"] = str(stats17["family_aggregation"])
            if "annotation_mode" in stats17:
                fig["fig17_stats_annotation_mode"] = str(stats17["annotation_mode"])
            if "show_legend" in stats17:
                fig["fig17_stats_show_legend"] = bool(stats17["show_legend"])
        if "exclude_community_name_patterns" in s17:
            fig["exclude_community_name_patterns_fig17"] = s17["exclude_community_name_patterns"]

    s18 = fig.get("fig18")
    if isinstance(s18, dict):
        if "enabled" in s18:
            fig["fig18"] = bool(s18["enabled"])
        else:
            fig["fig18"] = True
        if "title" in s18:
            fig["fig18_title"] = s18["title"]
        if "font_size" in s18:
            fig["fig18_font_size"] = float(s18["font_size"])
        if "plot_type" in s18:
            fig["fig18_plot_type"] = str(s18["plot_type"])
        if "exclude_community_name_patterns" in s18:
            fig["exclude_community_name_patterns_fig18"] = s18["exclude_community_name_patterns"]

    s19 = fig.get("fig19")
    if isinstance(s19, dict):
        if "enabled" in s19:
            fig["fig19"] = bool(s19["enabled"])
        else:
            fig["fig19"] = True
        if "title" in s19:
            fig["fig19_title"] = s19["title"]
        if "font_size" in s19:
            fig["fig19_font_size"] = float(s19["font_size"])
        if "point_size" in s19:
            fig["fig19_point_size"] = float(s19["point_size"])
        if "metrics" in s19 and isinstance(s19["metrics"], list):
            fig["fig19_metrics"] = [str(v) for v in s19["metrics"] if str(v).strip()]
        if "profile_distance" in s19:
            fig["fig19_profile_distance"] = str(s19["profile_distance"])
        if "standardize_profiles" in s19:
            fig["fig19_standardize_profiles"] = bool(s19["standardize_profiles"])
        if "show_umap" in s19:
            fig["fig19_show_umap"] = bool(s19["show_umap"])
        if "umap_n_neighbors" in s19:
            fig["fig19_umap_n_neighbors"] = int(s19["umap_n_neighbors"])
        if "umap_min_dist" in s19:
            fig["fig19_umap_min_dist"] = float(s19["umap_min_dist"])
        if "umap_random_state" in s19:
            fig["fig19_umap_random_state"] = int(s19["umap_random_state"])
        if "exclude_community_name_patterns" in s19:
            fig["exclude_community_name_patterns_fig19"] = s19["exclude_community_name_patterns"]

    s20 = fig.get("fig20")
    if isinstance(s20, dict):
        if "enabled" in s20:
            fig["fig20"] = bool(s20["enabled"])
        else:
            fig["fig20"] = True
        if "title" in s20:
            fig["fig20_title"] = s20["title"]
        if "font_size" in s20:
            fig["fig20_font_size"] = float(s20["font_size"])
        if "methods" in s20:
            raw_methods = s20["methods"]
            if isinstance(raw_methods, str):
                fig["fig20_methods"] = [
                    str(v).strip()
                    for v in re.split(r"[,;]", raw_methods)
                    if str(v).strip()
                ]
            elif isinstance(raw_methods, list):
                fig["fig20_methods"] = [str(v).strip() for v in raw_methods if str(v).strip()]
        if "sample_aggregation" in s20:
            fig["fig20_sample_aggregation"] = str(s20["sample_aggregation"])
        if "exclude_community_name_patterns" in s20:
            fig["exclude_community_name_patterns_fig20"] = s20["exclude_community_name_patterns"]

    s21 = fig.get("fig21")
    if isinstance(s21, dict):
        # Support panel-specific nested config under figures.fig21 while staying
        # backward compatible with the existing flat fig21 keys.
        s21_norm = dict(s21)
        p_a = s21.get("panel_a")
        p_b = s21.get("panel_b")
        p_c = s21.get("panel_c")
        p_d = s21.get("panel_d")
        p_e = s21.get("panel_e")
        p_f = s21.get("panel_f")
        if isinstance(p_a, dict):
            if "enabled" in p_a:
                s21_norm["show_panel_a"] = p_a["enabled"]
            if "title" in p_a:
                s21_norm["panel_a_title"] = p_a["title"]
            if "title_a" in p_a:
                s21_norm["panel_a_a_title"] = p_a["title_a"]
            if "title_b" in p_a:
                s21_norm["panel_a_b_title"] = p_a["title_b"]
            if "font_size" in p_a:
                s21_norm["panel_a_font_size"] = p_a["font_size"]
            if "xtick_rotation" in p_a:
                s21_norm["panel_a_xtick_rotation"] = p_a["xtick_rotation"]
            if "ytick_rotation" in p_a:
                s21_norm["panel_a_ytick_rotation"] = p_a["ytick_rotation"]
            if "show_planned_comparison_annotations" in p_a:
                s21_norm["panel_a_show_planned_comparison_annotations"] = p_a["show_planned_comparison_annotations"]
            if "x_axis_label" in p_a:
                s21_norm["panel_a_x_axis_label"] = p_a["x_axis_label"]
            if "y_axis_label" in p_a:
                s21_norm["panel_a_y_axis_label"] = p_a["y_axis_label"]
        if isinstance(p_b, dict):
            if "enabled" in p_b:
                s21_norm["show_panel_b"] = p_b["enabled"]
            if "title" in p_b:
                s21_norm["panel_b_title"] = p_b["title"]
            if "font_size" in p_b:
                s21_norm["panel_b_font_size"] = p_b["font_size"]
            if "xtick_rotation" in p_b:
                s21_norm["panel_b_xtick_rotation"] = p_b["xtick_rotation"]
            if "ytick_rotation" in p_b:
                s21_norm["panel_b_ytick_rotation"] = p_b["ytick_rotation"]
            if "sample_v4_multiple_label" in p_b:
                s21_norm["panel_b_sample_v4_multiple_label"] = p_b["sample_v4_multiple_label"]
            if "sample_v4_single_label" in p_b:
                s21_norm["panel_b_sample_v4_single_label"] = p_b["sample_v4_single_label"]
            if "dot_size" in p_b:
                s21_norm["panel_b_dot_size"] = p_b["dot_size"]
            if "dot_opacity" in p_b:
                s21_norm["panel_b_dot_opacity"] = p_b["dot_opacity"]
            if "dot_jitter" in p_b:
                s21_norm["panel_b_dot_jitter"] = p_b["dot_jitter"]
            if "x_axis_label" in p_b:
                s21_norm["panel_b_x_axis_label"] = p_b["x_axis_label"]
            if "y_axis_label" in p_b:
                s21_norm["panel_b_y_axis_label"] = p_b["y_axis_label"]
            if "pair_groups_to_show" in p_b and isinstance(p_b["pair_groups_to_show"], list):
                s21_norm["panel_b_pair_groups_to_show"] = p_b["pair_groups_to_show"]
            if "show_within_contrast_stats" in p_b:
                s21_norm["panel_b_show_within_contrast_stats"] = p_b["show_within_contrast_stats"]
            if "method_groups" in p_b and isinstance(p_b["method_groups"], dict):
                s21_norm["panel_b_method_groups"] = p_b["method_groups"]
            if "group_labels" in p_b and isinstance(p_b["group_labels"], dict):
                s21_norm["panel_b_group_labels"] = p_b["group_labels"]
        if isinstance(p_c, dict):
            if "enabled" in p_c:
                s21_norm["show_panel_c"] = p_c["enabled"]
            if "title" in p_c:
                s21_norm["panel_c_title"] = p_c["title"]
            if "font_size" in p_c:
                s21_norm["panel_c_font_size"] = p_c["font_size"]
            if "endpoint_matrix_gap_cells" in p_c:
                s21_norm["endpoint_matrix_gap_cells"] = p_c["endpoint_matrix_gap_cells"]
            if "endpoint_matrix_xtick_rotation" in p_c:
                s21_norm["endpoint_matrix_xtick_rotation"] = p_c["endpoint_matrix_xtick_rotation"]
            if "endpoint_matrix_ytick_rotation" in p_c:
                s21_norm["endpoint_matrix_ytick_rotation"] = p_c["endpoint_matrix_ytick_rotation"]
            if "title_a" in p_c:
                s21_norm["panel_c_a_title"] = p_c["title_a"]
            if "title_b" in p_c:
                s21_norm["panel_c_b_title"] = p_c["title_b"]
            if "a_xtick_rotation" in p_c:
                s21_norm["panel_c_a_xtick_rotation"] = p_c["a_xtick_rotation"]
            if "a_ytick_rotation" in p_c:
                s21_norm["panel_c_a_ytick_rotation"] = p_c["a_ytick_rotation"]
            if "b_xtick_rotation" in p_c:
                s21_norm["panel_c_b_xtick_rotation"] = p_c["b_xtick_rotation"]
            if "b_ytick_rotation" in p_c:
                s21_norm["panel_c_b_ytick_rotation"] = p_c["b_ytick_rotation"]
            if "x_axis_label" in p_c:
                s21_norm["panel_c_x_axis_label"] = p_c["x_axis_label"]
            if "y_axis_label" in p_c:
                s21_norm["panel_c_y_axis_label"] = p_c["y_axis_label"]
        if isinstance(p_d, dict):
            if "enabled" in p_d:
                s21_norm["show_panel_d"] = p_d["enabled"]
            if "title" in p_d:
                s21_norm["panel_d_title"] = p_d["title"]
            if "font_size" in p_d:
                s21_norm["panel_d_font_size"] = p_d["font_size"]
            if "dot_size" in p_d:
                s21_norm["panel_d_dot_size"] = p_d["dot_size"]
            if "dot_opacity" in p_d:
                s21_norm["panel_d_dot_opacity"] = p_d["dot_opacity"]
            if "dot_jitter" in p_d:
                s21_norm["panel_d_dot_jitter"] = p_d["dot_jitter"]
            if "x_axis_label" in p_d:
                s21_norm["panel_d_x_axis_label"] = p_d["x_axis_label"]
            if "y_axis_label" in p_d:
                s21_norm["panel_d_y_axis_label"] = p_d["y_axis_label"]
            if "methods" in p_d:
                s21_norm["panel_d_methods"] = p_d["methods"]
            p_d_stats = p_d.get("stats")
            if isinstance(p_d_stats, dict) and "show_annotations" in p_d_stats:
                s21_norm["panel_d_stats_show_annotations"] = p_d_stats["show_annotations"]
        if isinstance(p_e, dict):
            if "enabled" in p_e:
                s21_norm["show_panel_e"] = p_e["enabled"]
            if "title" in p_e:
                s21_norm["panel_e_title"] = p_e["title"]
            if "font_size" in p_e:
                s21_norm["panel_e_font_size"] = p_e["font_size"]
            if "dot_size" in p_e:
                s21_norm["panel_e_dot_size"] = p_e["dot_size"]
            if "dot_opacity" in p_e:
                s21_norm["panel_e_dot_opacity"] = p_e["dot_opacity"]
            if "dot_jitter" in p_e:
                s21_norm["panel_e_dot_jitter"] = p_e["dot_jitter"]
            if "grid_label_font_size" in p_e:
                s21_norm["panel_e_grid_label_font_size"] = p_e["grid_label_font_size"]
            if "grid_tick_font_size" in p_e:
                s21_norm["panel_e_grid_tick_font_size"] = p_e["grid_tick_font_size"]
            if "grid_text_font_size" in p_e:
                s21_norm["panel_e_grid_text_font_size"] = p_e["grid_text_font_size"]
            if "grid_legend_font_size" in p_e:
                s21_norm["panel_e_grid_legend_font_size"] = p_e["grid_legend_font_size"]
            if "abs_density_label" in p_e:
                s21_norm["panel_e_abs_density_label"] = p_e["abs_density_label"]
            if "summary_density_label" in p_e:
                s21_norm["panel_e_summary_density_label"] = p_e["summary_density_label"]
            if "show_summary_density" in p_e:
                s21_norm["panel_e_show_summary_density"] = p_e["show_summary_density"]
            if "mean_label" in p_e:
                s21_norm["panel_e_mean_label"] = p_e["mean_label"]
            if "sample_v4_multiple_label" in p_e:
                s21_norm["panel_e_sample_v4_multiple_label"] = p_e["sample_v4_multiple_label"]
            if "sample_v4_single_label" in p_e:
                s21_norm["panel_e_sample_v4_single_label"] = p_e["sample_v4_single_label"]
            if "x_axis_label" in p_e:
                s21_norm["panel_e_x_axis_label"] = p_e["x_axis_label"]
            if "y_axis_label" in p_e:
                s21_norm["panel_e_y_axis_label"] = p_e["y_axis_label"]
        if isinstance(p_f, dict):
            if "enabled" in p_f:
                s21_norm["show_panel_f"] = p_f["enabled"]
            if "title" in p_f:
                s21_norm["panel_f_title"] = p_f["title"]
            if "font_size" in p_f:
                s21_norm["panel_f_font_size"] = p_f["font_size"]
            if "dot_size" in p_f:
                s21_norm["panel_f_dot_size"] = p_f["dot_size"]
            if "dot_opacity" in p_f:
                s21_norm["panel_f_dot_opacity"] = p_f["dot_opacity"]
            if "dot_jitter" in p_f:
                s21_norm["panel_f_dot_jitter"] = p_f["dot_jitter"]
            if "x_axis_label" in p_f:
                s21_norm["panel_f_x_axis_label"] = p_f["x_axis_label"]
            if "y_axis_label" in p_f:
                s21_norm["panel_f_y_axis_label"] = p_f["y_axis_label"]
            if "x_tick_rotation" in p_f:
                s21_norm["panel_f_xtick_rotation"] = p_f["x_tick_rotation"]
            if "methods" in p_f:
                s21_norm["panel_f_methods"] = p_f["methods"]
            if "summary_statistic" in p_f:
                s21_norm["panel_f_summary_statistic"] = p_f["summary_statistic"]
            p_f_stats = p_f.get("stats")
            if isinstance(p_f_stats, dict):
                if "pairwise_test" in p_f_stats:
                    s21_norm["panel_f_stats_pairwise_test"] = p_f_stats["pairwise_test"]
                if "alternative" in p_f_stats:
                    s21_norm["panel_f_stats_alternative"] = p_f_stats["alternative"]
                if "practical_difference_delta" in p_f_stats:
                    s21_norm["panel_f_stats_practical_difference_delta"] = p_f_stats["practical_difference_delta"]
                if "bootstrap_resamples" in p_f_stats:
                    s21_norm["panel_f_stats_bootstrap_resamples"] = p_f_stats["bootstrap_resamples"]
                if "ci_level" in p_f_stats:
                    s21_norm["panel_f_stats_ci_level"] = p_f_stats["ci_level"]
                if "show_annotations" in p_f_stats:
                    s21_norm["panel_f_stats_show_annotations"] = p_f_stats["show_annotations"]
                if "show_non_significant_annotations" in p_f_stats:
                    s21_norm["panel_f_stats_show_non_significant_annotations"] = p_f_stats[
                        "show_non_significant_annotations"
                    ]
                if "planned_comparisons" in p_f_stats:
                    s21_norm["panel_f_stats_planned_comparisons"] = p_f_stats["planned_comparisons"]
        s21 = s21_norm

        if "enabled" in s21:
            fig["fig21"] = bool(s21["enabled"])
        else:
            fig["fig21"] = True
        if "title" in s21:
            fig["fig21_title"] = s21["title"]
        if "panel_a_title" in s21:
            fig["fig21_panel_a_title"] = s21["panel_a_title"]
        if "panel_a_a_title" in s21:
            fig["fig21_panel_a_a_title"] = s21["panel_a_a_title"]
        if "panel_a_b_title" in s21:
            fig["fig21_panel_a_b_title"] = s21["panel_a_b_title"]
        if "panel_b_title" in s21:
            fig["fig21_panel_b_title"] = s21["panel_b_title"]
        if "panel_c_title" in s21:
            fig["fig21_panel_c_title"] = s21["panel_c_title"]
        if "panel_c_a_title" in s21:
            fig["fig21_panel_c_a_title"] = s21["panel_c_a_title"]
        if "panel_c_b_title" in s21:
            fig["fig21_panel_c_b_title"] = s21["panel_c_b_title"]
        if "panel_d_title" in s21:
            fig["fig21_panel_d_title"] = s21["panel_d_title"]
        if "panel_e_title" in s21:
            fig["fig21_panel_e_title"] = s21["panel_e_title"]
        if "panel_f_title" in s21:
            fig["fig21_panel_f_title"] = s21["panel_f_title"]
        if "font_size" in s21:
            fig["fig21_font_size"] = float(s21["font_size"])
        if "panel_a_font_size" in s21:
            fig["fig21_panel_a_font_size"] = float(s21["panel_a_font_size"])
        if "panel_b_font_size" in s21:
            fig["fig21_panel_b_font_size"] = float(s21["panel_b_font_size"])
        if "panel_c_font_size" in s21:
            fig["fig21_panel_c_font_size"] = float(s21["panel_c_font_size"])
        if "panel_d_font_size" in s21:
            fig["fig21_panel_d_font_size"] = float(s21["panel_d_font_size"])
        if "panel_e_font_size" in s21:
            fig["fig21_panel_e_font_size"] = float(s21["panel_e_font_size"])
        if "panel_f_font_size" in s21:
            fig["fig21_panel_f_font_size"] = float(s21["panel_f_font_size"])
        if "panel_e_grid_label_font_size" in s21:
            fig["fig21_panel_e_grid_label_font_size"] = _optional_float(s21["panel_e_grid_label_font_size"])
        if "panel_e_grid_tick_font_size" in s21:
            fig["fig21_panel_e_grid_tick_font_size"] = _optional_float(s21["panel_e_grid_tick_font_size"])
        if "panel_e_grid_text_font_size" in s21:
            fig["fig21_panel_e_grid_text_font_size"] = _optional_float(s21["panel_e_grid_text_font_size"])
        if "panel_e_grid_legend_font_size" in s21:
            fig["fig21_panel_e_grid_legend_font_size"] = _optional_float(s21["panel_e_grid_legend_font_size"])
        if "panel_e_abs_density_label" in s21:
            fig["fig21_panel_e_abs_density_label"] = str(s21["panel_e_abs_density_label"])
        if "panel_e_summary_density_label" in s21:
            fig["fig21_panel_e_summary_density_label"] = str(s21["panel_e_summary_density_label"])
        if "panel_e_show_summary_density" in s21:
            fig["fig21_panel_e_show_summary_density"] = bool(s21["panel_e_show_summary_density"])
        if "panel_e_mean_label" in s21:
            fig["fig21_panel_e_mean_label"] = str(s21["panel_e_mean_label"])
        if "panel_e_sample_v4_multiple_label" in s21:
            fig["fig21_panel_e_sample_v4_multiple_label"] = str(s21["panel_e_sample_v4_multiple_label"])
        if "panel_e_sample_v4_single_label" in s21:
            fig["fig21_panel_e_sample_v4_single_label"] = str(s21["panel_e_sample_v4_single_label"])
        if "panel_a_x_axis_label" in s21:
            fig["fig21_panel_a_x_axis_label"] = str(s21["panel_a_x_axis_label"])
        if "panel_a_y_axis_label" in s21:
            fig["fig21_panel_a_y_axis_label"] = str(s21["panel_a_y_axis_label"])
        if "panel_b_x_axis_label" in s21:
            fig["fig21_panel_b_x_axis_label"] = str(s21["panel_b_x_axis_label"])
        if "panel_b_y_axis_label" in s21:
            fig["fig21_panel_b_y_axis_label"] = str(s21["panel_b_y_axis_label"])
        if "panel_c_x_axis_label" in s21:
            fig["fig21_panel_c_x_axis_label"] = str(s21["panel_c_x_axis_label"])
        if "panel_c_y_axis_label" in s21:
            fig["fig21_panel_c_y_axis_label"] = str(s21["panel_c_y_axis_label"])
        if "panel_d_x_axis_label" in s21:
            fig["fig21_panel_d_x_axis_label"] = str(s21["panel_d_x_axis_label"])
        if "panel_d_y_axis_label" in s21:
            fig["fig21_panel_d_y_axis_label"] = str(s21["panel_d_y_axis_label"])
        if "panel_e_x_axis_label" in s21:
            fig["fig21_panel_e_x_axis_label"] = str(s21["panel_e_x_axis_label"])
        if "panel_e_y_axis_label" in s21:
            fig["fig21_panel_e_y_axis_label"] = str(s21["panel_e_y_axis_label"])
        if "panel_f_x_axis_label" in s21:
            fig["fig21_panel_f_x_axis_label"] = str(s21["panel_f_x_axis_label"])
        if "panel_f_y_axis_label" in s21:
            fig["fig21_panel_f_y_axis_label"] = str(s21["panel_f_y_axis_label"])
        if "panel_f_xtick_rotation" in s21:
            fig["fig21_panel_f_xtick_rotation"] = float(s21["panel_f_xtick_rotation"])
        if "panel_b_pair_groups_to_show" in s21 and isinstance(s21["panel_b_pair_groups_to_show"], list):
            fig["fig21_panel_b_pair_groups_to_show"] = [str(v) for v in s21["panel_b_pair_groups_to_show"] if str(v).strip()]
        if "panel_b_show_within_contrast_stats" in s21:
            fig["fig21_panel_b_show_within_contrast_stats"] = bool(s21["panel_b_show_within_contrast_stats"])
        if "panel_b_method_groups" in s21 and isinstance(s21["panel_b_method_groups"], dict):
            fig["fig21_panel_b_method_groups"] = {
                str(group): [str(method) for method in methods]
                for group, methods in s21["panel_b_method_groups"].items()
                if isinstance(methods, list)
            }
        if "panel_b_group_labels" in s21 and isinstance(s21["panel_b_group_labels"], dict):
            fig["fig21_panel_b_group_labels"] = {
                str(group): str(label)
                for group, label in s21["panel_b_group_labels"].items()
                if str(label).strip()
            }
        if "panel_b_sample_v4_multiple_label" in s21:
            fig["fig21_panel_b_sample_v4_multiple_label"] = str(s21["panel_b_sample_v4_multiple_label"])
        if "panel_b_sample_v4_single_label" in s21:
            fig["fig21_panel_b_sample_v4_single_label"] = str(s21["panel_b_sample_v4_single_label"])
        if "panel_b_dot_size" in s21:
            fig["fig21_panel_b_dot_size"] = _optional_float(s21["panel_b_dot_size"])
        if "panel_b_dot_opacity" in s21:
            fig["fig21_panel_b_dot_opacity"] = float(s21["panel_b_dot_opacity"])
        if "panel_b_dot_jitter" in s21:
            fig["fig21_panel_b_dot_jitter"] = float(s21["panel_b_dot_jitter"])
        if "panel_d_dot_size" in s21:
            fig["fig21_panel_d_dot_size"] = _optional_float(s21["panel_d_dot_size"])
        if "panel_d_dot_opacity" in s21:
            fig["fig21_panel_d_dot_opacity"] = float(s21["panel_d_dot_opacity"])
        if "panel_d_dot_jitter" in s21:
            fig["fig21_panel_d_dot_jitter"] = float(s21["panel_d_dot_jitter"])
        if "panel_d_methods" in s21:
            raw_methods = s21["panel_d_methods"]
            fig["fig21_panel_d_methods"] = [str(v) for v in raw_methods] if isinstance(raw_methods, list) else []
        if "panel_d_stats_show_annotations" in s21:
            fig["fig21_panel_d_stats_show_annotations"] = bool(s21["panel_d_stats_show_annotations"])
        if "panel_e_dot_size" in s21:
            fig["fig21_panel_e_dot_size"] = _optional_float(s21["panel_e_dot_size"])
        if "panel_e_dot_opacity" in s21:
            fig["fig21_panel_e_dot_opacity"] = float(s21["panel_e_dot_opacity"])
        if "panel_e_dot_jitter" in s21:
            fig["fig21_panel_e_dot_jitter"] = float(s21["panel_e_dot_jitter"])
        if "panel_f_dot_size" in s21:
            fig["fig21_panel_f_dot_size"] = _optional_float(s21["panel_f_dot_size"])
        if "panel_f_dot_opacity" in s21:
            fig["fig21_panel_f_dot_opacity"] = float(s21["panel_f_dot_opacity"])
        if "panel_f_dot_jitter" in s21:
            fig["fig21_panel_f_dot_jitter"] = float(s21["panel_f_dot_jitter"])
        if "panel_f_methods" in s21:
            raw_methods = s21["panel_f_methods"]
            fig["fig21_panel_f_methods"] = [str(v) for v in raw_methods] if isinstance(raw_methods, list) else []
        if "panel_f_summary_statistic" in s21:
            fig["fig21_panel_f_summary_statistic"] = str(s21["panel_f_summary_statistic"])
        if "panel_f_stats_pairwise_test" in s21:
            fig["fig21_panel_f_stats_pairwise_test"] = str(s21["panel_f_stats_pairwise_test"])
        if "panel_f_stats_alternative" in s21:
            fig["fig21_panel_f_stats_alternative"] = str(s21["panel_f_stats_alternative"])
        if "panel_f_stats_practical_difference_delta" in s21:
            fig["fig21_panel_f_stats_practical_difference_delta"] = float(
                s21["panel_f_stats_practical_difference_delta"]
            )
        if "panel_f_stats_bootstrap_resamples" in s21:
            fig["fig21_panel_f_stats_bootstrap_resamples"] = int(s21["panel_f_stats_bootstrap_resamples"])
        if "panel_f_stats_ci_level" in s21:
            fig["fig21_panel_f_stats_ci_level"] = float(s21["panel_f_stats_ci_level"])
        if "panel_f_stats_show_annotations" in s21:
            fig["fig21_panel_f_stats_show_annotations"] = bool(s21["panel_f_stats_show_annotations"])
        if "panel_f_stats_show_non_significant_annotations" in s21:
            fig["fig21_panel_f_stats_show_non_significant_annotations"] = bool(
                s21["panel_f_stats_show_non_significant_annotations"]
            )
        if "panel_f_stats_planned_comparisons" in s21:
            raw_pairs = s21["panel_f_stats_planned_comparisons"]
            fig["fig21_panel_f_stats_planned_comparisons"] = [
                [str(method) for method in pair]
                for pair in raw_pairs
                if isinstance(pair, list) and len(pair) == 2
            ] if isinstance(raw_pairs, list) else []
        if "panel_a_xtick_rotation" in s21:
            fig["fig21_panel_a_xtick_rotation"] = float(s21["panel_a_xtick_rotation"])
        if "panel_a_ytick_rotation" in s21:
            fig["fig21_panel_a_ytick_rotation"] = float(s21["panel_a_ytick_rotation"])
        if "panel_a_show_planned_comparison_annotations" in s21:
            fig["fig21_panel_a_show_planned_comparison_annotations"] = bool(
                s21["panel_a_show_planned_comparison_annotations"]
            )
        if "panel_b_xtick_rotation" in s21:
            fig["fig21_panel_b_xtick_rotation"] = float(s21["panel_b_xtick_rotation"])
        if "panel_b_ytick_rotation" in s21:
            fig["fig21_panel_b_ytick_rotation"] = float(s21["panel_b_ytick_rotation"])
        if "endpoint_matrix_gap_cells" in s21:
            fig["fig21_endpoint_matrix_gap_cells"] = int(s21["endpoint_matrix_gap_cells"])
        if "endpoint_matrix_xtick_rotation" in s21:
            fig["fig21_endpoint_matrix_xtick_rotation"] = float(s21["endpoint_matrix_xtick_rotation"])
        if "endpoint_matrix_ytick_rotation" in s21:
            fig["fig21_endpoint_matrix_ytick_rotation"] = float(s21["endpoint_matrix_ytick_rotation"])
        if "panel_c_a_xtick_rotation" in s21:
            v = s21["panel_c_a_xtick_rotation"]
            fig["fig21_panel_c_a_xtick_rotation"] = None if v is None else float(v)
        if "panel_c_a_ytick_rotation" in s21:
            v = s21["panel_c_a_ytick_rotation"]
            fig["fig21_panel_c_a_ytick_rotation"] = None if v is None else float(v)
        if "panel_c_b_xtick_rotation" in s21:
            v = s21["panel_c_b_xtick_rotation"]
            fig["fig21_panel_c_b_xtick_rotation"] = None if v is None else float(v)
        if "panel_c_b_ytick_rotation" in s21:
            v = s21["panel_c_b_ytick_rotation"]
            fig["fig21_panel_c_b_ytick_rotation"] = None if v is None else float(v)
        if "plot_type" in s21:
            fig["fig21_plot_type"] = str(s21["plot_type"])
        if "violin_bw_method" in s21:
            fig["fig21_violin_bw_method"] = s21["violin_bw_method"]
        if "point_size" in s21:
            fig["fig21_point_size"] = float(s21["point_size"])
        if "metrics" in s21 and isinstance(s21["metrics"], list):
            fig["fig21_metrics"] = [str(v) for v in s21["metrics"] if str(v).strip()]
        if "agnostic_metrics" in s21 and isinstance(s21["agnostic_metrics"], list):
            fig["fig21_agnostic_metrics"] = [str(v) for v in s21["agnostic_metrics"] if str(v).strip()]
        if "aware_metrics" in s21 and isinstance(s21["aware_metrics"], list):
            fig["fig21_aware_metrics"] = [str(v) for v in s21["aware_metrics"] if str(v).strip()]
        if "profile_distance" in s21:
            fig["fig21_profile_distance"] = str(s21["profile_distance"])
        if "standardize_profiles" in s21:
            fig["fig21_standardize_profiles"] = bool(s21["standardize_profiles"])
        if "show_pair_group_distance_contrast" in s21:
            fig["fig21_show_pair_group_distance_contrast"] = bool(s21["show_pair_group_distance_contrast"])
        if "pair_group_summary_stat" in s21:
            fig["fig21_pair_group_summary_stat"] = str(s21["pair_group_summary_stat"])
        if "show_sensitivity_by_method" in s21:
            fig["fig21_show_sensitivity_by_method"] = bool(s21["show_sensitivity_by_method"])
            if "fig22" not in fig:
                fig["fig22"] = bool(s21["show_sensitivity_by_method"])
        if "show_sample_method_heatmap" in s21:
            fig["fig21_show_sample_method_heatmap"] = bool(s21["show_sample_method_heatmap"])
            if "fig23" not in fig:
                fig["fig23"] = bool(s21["show_sample_method_heatmap"])
        if "show_per_sample_method_metric_pcoa" in s21:
            fig["fig21_show_per_sample_method_metric_pcoa"] = bool(s21["show_per_sample_method_metric_pcoa"])
            if "fig24" not in fig:
                fig["fig24"] = bool(s21["show_per_sample_method_metric_pcoa"])
        if "show_method_metric_profile_pcoa" in s21:
            fig["fig21_show_method_metric_profile_pcoa"] = bool(s21["show_method_metric_profile_pcoa"])
            if "fig25" not in fig:
                fig["fig25"] = bool(s21["show_method_metric_profile_pcoa"])
        if "show_pairwise_difference_pcoa" in s21:
            fig["fig21_show_pairwise_difference_pcoa"] = bool(s21["show_pairwise_difference_pcoa"])
            if "fig26" not in fig:
                fig["fig26"] = bool(s21["show_pairwise_difference_pcoa"])
        if "show_panel_e_raincloud_grid" in s21:
            fig["fig21_show_panel_e_raincloud_grid"] = bool(s21["show_panel_e_raincloud_grid"])
        if "show_panel_a" in s21:
            fig["fig21_show_panel_a"] = bool(s21["show_panel_a"])
        if "show_panel_b" in s21:
            fig["fig21_show_panel_b"] = bool(s21["show_panel_b"])
        if "show_panel_c" in s21:
            fig["fig21_show_panel_c"] = bool(s21["show_panel_c"])
        if "show_panel_d" in s21:
            fig["fig21_show_panel_d"] = bool(s21["show_panel_d"])
        if "show_panel_e" in s21:
            fig["fig21_show_panel_e"] = bool(s21["show_panel_e"])
        if "show_panel_f" in s21:
            fig["fig21_show_panel_f"] = bool(s21["show_panel_f"])
        if "projection_mapping_mode" in s21:
            fig["fig21_projection_mapping_mode"] = str(s21["projection_mapping_mode"])
        if "target_level_by_method" in s21 and isinstance(s21["target_level_by_method"], dict):
            fig["fig21_target_level_by_method"] = s21["target_level_by_method"]
        if "min_effective_pident_by_method" in s21 and isinstance(s21["min_effective_pident_by_method"], dict):
            fig["fig21_min_effective_pident_by_method"] = s21["min_effective_pident_by_method"]
        s21_stats = s21.get("stats")
        if isinstance(s21_stats, dict):
            if "enabled" in s21_stats:
                fig["fig21_stats_enabled"] = bool(s21_stats["enabled"])
            if "permutations" in s21_stats:
                fig["fig21_stats_permutations"] = int(s21_stats["permutations"])
            if "seed" in s21_stats:
                fig["fig21_stats_seed"] = int(s21_stats["seed"])
            if "alpha" in s21_stats:
                fig["fig21_stats_alpha"] = float(s21_stats["alpha"])
            if "show_pairwise" in s21_stats:
                fig["fig21_stats_show_pairwise"] = bool(s21_stats["show_pairwise"])
            if "pairwise_test" in s21_stats:
                fig["fig21_stats_pairwise_test"] = str(s21_stats["pairwise_test"])
            if "multiple_test_correction" in s21_stats:
                fig["fig21_stats_multiple_test_correction"] = bool(s21_stats["multiple_test_correction"])

    s22 = fig.get("fig22")
    if isinstance(s22, dict):
        if "enabled" in s22:
            fig["fig22"] = bool(s22["enabled"])
        else:
            fig["fig22"] = True
        if "title" in s22:
            fig["fig22_title"] = s22["title"]
        if "font_size" in s22:
            fig["fig22_font_size"] = float(s22["font_size"])

    s23 = fig.get("fig23")
    if isinstance(s23, dict):
        if "enabled" in s23:
            fig["fig23"] = bool(s23["enabled"])
        else:
            fig["fig23"] = True
        if "title" in s23:
            fig["fig23_title"] = s23["title"]
        if "font_size" in s23:
            fig["fig23_font_size"] = float(s23["font_size"])

    s24 = fig.get("fig24")
    if isinstance(s24, dict):
        if "enabled" in s24:
            fig["fig24"] = bool(s24["enabled"])
        else:
            fig["fig24"] = True
        if "title" in s24:
            fig["fig24_title"] = s24["title"]
        if "font_size" in s24:
            fig["fig24_font_size"] = float(s24["font_size"])

    s25 = fig.get("fig25")
    if isinstance(s25, dict):
        if "enabled" in s25:
            fig["fig25"] = bool(s25["enabled"])
        else:
            fig["fig25"] = True
        if "title" in s25:
            fig["fig25_title"] = s25["title"]
        if "font_size" in s25:
            fig["fig25_font_size"] = float(s25["font_size"])

    s26 = fig.get("fig26")
    if isinstance(s26, dict):
        if "enabled" in s26:
            fig["fig26"] = bool(s26["enabled"])
        else:
            fig["fig26"] = True
        if "title" in s26:
            fig["fig26_title"] = s26["title"]
        if "font_size" in s26:
            fig["fig26_font_size"] = float(s26["font_size"])
        if "exclude_community_name_patterns" in s21:
            fig["exclude_community_name_patterns_fig21"] = s21["exclude_community_name_patterns"]

    # Backward-compatible flat aliases for Figure 09 effective-identity keys.
    if "fig09_graph_min_effective_pident_denoising" not in fig and "fig09_graph_min_identity_denoising" in fig:
        fig["fig09_graph_min_effective_pident_denoising"] = fig["fig09_graph_min_identity_denoising"]
    if "fig09_graph_min_effective_pident_tic_tac" not in fig and "fig09_graph_min_identity_tic_tac" in fig:
        fig["fig09_graph_min_effective_pident_tic_tac"] = fig["fig09_graph_min_identity_tic_tac"]
    if "fig09_graph_min_effective_pident_otu" not in fig and "fig09_graph_min_identity_otu" in fig:
        fig["fig09_graph_min_effective_pident_otu"] = fig["fig09_graph_min_identity_otu"]
    fig.pop("fig09_graph_min_identity_denoising", None)
    fig.pop("fig09_graph_min_identity_tic_tac", None)
    fig.pop("fig09_graph_min_identity_otu", None)
    fig.pop("fig14_graph_min_identity_denoising", None)
    fig.pop("fig14_graph_min_identity_tic_tac", None)
    fig.pop("fig14_graph_min_identity_otu", None)
    # Deprecated / removed settings (ignored if present in legacy configs).
    fig.pop("fig09_min_alignment_length", None)
    fig.pop("fig09_unmapped_best_hit_min_alignment_length", None)
    fig.pop("fig09_unmapped_best_hit_min_query_coverage", None)

    return fig


def _normalize_mapping_data(raw: Dict[str, Any]) -> Dict[str, Any]:
    """
    Normalize mapping config and keep backward compatibility with legacy names:
      - min_identity -> min_effective_pident
      - min_identity_tic_tac -> min_effective_pident_tic_tac
      - min_identity_clustering -> min_effective_pident_clustering
      - blast_record_min_identity -> blast_record_min_raw_pident
      - blast_record_min_effective_pident -> blast_record_min_raw_pident
    Supports global per-method mapping target level override:
      - target_level_by_method (preferred)
      - mapping_target_level_by_method (legacy alias)
    """
    m = dict(raw or {})

    if "min_effective_pident" not in m and "min_identity" in m:
        m["min_effective_pident"] = m["min_identity"]
    if "min_effective_pident_tic_tac" not in m and "min_identity_tic_tac" in m:
        m["min_effective_pident_tic_tac"] = m["min_identity_tic_tac"]
    if "min_effective_pident_clustering" not in m and "min_identity_clustering" in m:
        m["min_effective_pident_clustering"] = m["min_identity_clustering"]
    if "blast_record_min_raw_pident" not in m:
        if "blast_record_min_effective_pident" in m:
            m["blast_record_min_raw_pident"] = m["blast_record_min_effective_pident"]
        elif "blast_record_min_identity" in m:
            m["blast_record_min_raw_pident"] = m["blast_record_min_identity"]

    # Drop legacy aliases so dataclass construction is clean.
    m.pop("min_identity", None)
    m.pop("min_identity_tic_tac", None)
    m.pop("min_identity_clustering", None)
    # Deprecated module-wide method-specific mapping identity knobs (moved to fig11 scope).
    m.pop("min_effective_pident_by_method", None)
    m.pop("min_identity_by_method", None)
    m.pop("mapping_target_level_by_method", None)
    m.pop("blast_record_min_identity", None)
    m.pop("blast_record_min_effective_pident", None)
    # Deprecated / removed settings (ignored if present in legacy configs).
    m.pop("min_alignment_length", None)
    return m


def load_config(path: Path, overrides: Optional[Dict[str, Any]] = None) -> Config:
    """
    Load YAML config and allow a small set of CLI overrides:
      - root
      - outdir
    """
    data = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
    overrides = overrides or {}

    paths = (data.get("paths", {}) or {}).copy()
    for k in ("root", "outdir"):
        if k in overrides and overrides[k] is not None:
            paths[k] = overrides[k]

    if "root" not in paths or "outdir" not in paths:
        raise ValueError("Config must define paths.root and paths.outdir")

    figures_data = _normalize_figures_data(data.get("figures", {}) or {})
    mapping_data = _normalize_mapping_data(data.get("mapping", {}) or {})

    return Config(
        paths=PathsConfig(
            root=Path(paths["root"]).expanduser(),
            outdir=Path(paths["outdir"]).expanduser(),
        ),
        discovery=DiscoveryConfig(**(data.get("discovery", {}) or {})),
        normalization=NormalizationConfig(**(data.get("normalization", {}) or {})),
        mapping=MappingConfig(**mapping_data),
        tree=TreeConfig(**(data.get("tree", {}) or {})),
        distances=DistancesConfig(**(data.get("distances", {}) or {})),
        figures=FiguresConfig(**figures_data),
    )
    if "target_level_by_method" not in m and "mapping_target_level_by_method" in m:
        m["target_level_by_method"] = m["mapping_target_level_by_method"]
