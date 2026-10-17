from __future__ import annotations

from ..log import get_logger
from ..utils.labels import configure_community_display_mapping_from_manifest
from . import (
    build_tree_centroid,
    build_tree_full,
    compute_distances,
    discover,
    map_features,
    plot_fig01,
    plot_fig02,
    plot_fig03,
    plot_fig04,
    plot_fig05,
    plot_fig06,
    plot_fig07,
    plot_fig08,
    plot_fig09,
    plot_fig10,
    plot_fig11,
    plot_fig12,
    plot_fig13,
    plot_fig14,
    plot_fig15,
    plot_fig16,
    plot_fig17,
    plot_fig18,
    plot_fig19,
    plot_fig20,
    plot_fig21,
    prepare_expected,
    prepare_methods,
    project_abundances,
)


def run(cfg) -> int:
    """
    Run the full benchmarking pipeline in the intended order.

    Steps:
      1) discover
      2) prepare-expected
      3) prepare-methods
      4) build-tree-full (optional)
      5) map-features
      6) build-tree-centroid (optional)
      7) project-abundances
      8) compute-distances
      9) plot figures (01/02/03/04/05/06/07/08/09/10/11/12/13/14/15/16/17/18/19/20/21/22/23/24/25/26 depending on config)
    """
    log = get_logger()
    log.info("Running tacticbench pipeline...")

    discover.run(cfg)
    configure_community_display_mapping_from_manifest(
        cfg.paths.outdir,
        figure_exclude_patterns=list(getattr(cfg.figures, "exclude_community_name_patterns", []) or []),
    )
    prepare_expected.run(cfg)
    prepare_methods.run(cfg)

    if getattr(cfg.tree, "build_full_tree", True):
        build_tree_full.run(cfg)

    map_features.run(cfg)

    if getattr(cfg.tree, "build_centroid_tree", True):
        build_tree_centroid.run(cfg)

    project_abundances.run(cfg)
    compute_distances.run(cfg)

    if getattr(cfg.figures, "fig01", True):
        plot_fig01.run(cfg)
    if getattr(cfg.figures, "fig02", True):
        plot_fig02.run(cfg)
    if getattr(cfg.figures, "fig03_by_sample", True) or getattr(cfg.figures, "fig03_by_metric", True):
        plot_fig03.run(cfg)
    if getattr(cfg.figures, "fig04", True):
        plot_fig04.run(cfg)
    if getattr(cfg.figures, "fig05", True):
        plot_fig05.run(cfg)
    if getattr(cfg.figures, "fig06", True):
        plot_fig06.run(cfg)
    if getattr(cfg.figures, "fig07", True):
        plot_fig07.run(cfg)
    if getattr(cfg.figures, "fig08", True):
        plot_fig08.run(cfg)
    if getattr(cfg.figures, "fig09", True):
        plot_fig09.run(cfg)
    if getattr(cfg.figures, "fig10", True):
        plot_fig10.run(cfg)
    if getattr(cfg.figures, "fig11", True):
        plot_fig11.run(cfg)
    if getattr(cfg.figures, "fig12", True):
        plot_fig12.run(cfg)
    if getattr(cfg.figures, "fig13", True):
        plot_fig13.run(cfg)
    if getattr(cfg.figures, "fig14", True):
        plot_fig14.run(cfg)
    if getattr(cfg.figures, "fig15", True):
        plot_fig15.run(cfg)
    if getattr(cfg.figures, "fig16", True):
        plot_fig16.run(cfg)
    if getattr(cfg.figures, "fig17", True):
        plot_fig17.run(cfg)
    if getattr(cfg.figures, "fig18", True):
        plot_fig18.run(cfg)
    if getattr(cfg.figures, "fig19", True):
        plot_fig19.run(cfg)
    if getattr(cfg.figures, "fig20", True):
        plot_fig20.run(cfg)
    if any(
        bool(getattr(cfg.figures, k, False))
        for k in ("fig21", "fig22", "fig23", "fig24", "fig25", "fig26")
    ):
        plot_fig21.run(cfg)

    log.info("Done.")
    return 0
