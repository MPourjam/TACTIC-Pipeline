from __future__ import annotations

import argparse
from dataclasses import asdict
from pathlib import Path
from typing import Optional

import yaml

from .config import load_config
from .log import get_logger
from .utils.labels import configure_community_display_mapping_from_manifest
from .steps import (
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
    plot_all,
    prepare_expected,
    prepare_methods,
    project_abundances,
    run_all,
)


def _add_common(sp: argparse.ArgumentParser) -> None:
    sp.add_argument("--config", required=True, type=Path, help="Path to YAML config.")
    sp.add_argument("--root", type=Path, default=None, help="Override config.paths.root")
    sp.add_argument("--outdir", type=Path, default=None, help="Override config.paths.outdir")


def _to_yamlable(obj):
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, dict):
        return {k: _to_yamlable(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_to_yamlable(v) for v in obj]
    return obj


def _write_config_snapshots(outdir: Path, config_path: Path, cfg) -> None:
    """
    Persist config provenance in each analysis result directory:
      - config_input.yaml: exact copy of CLI-provided config file
      - config_effective.yaml: resolved config actually used (after CLI overrides)
    Files are always overwritten so they stay current with each rerun.
    """
    input_dst = outdir / "config_input.yaml"
    effective_dst = outdir / "config_effective.yaml"

    input_text = config_path.read_text(encoding="utf-8")
    input_dst.write_text(input_text, encoding="utf-8")

    eff = _to_yamlable(asdict(cfg))
    effective_dst.write_text(yaml.safe_dump(eff, sort_keys=False), encoding="utf-8")


def _cfg_from_args(args) -> "Config":
    overrides = {}
    if args.root is not None:
        overrides["root"] = str(args.root)
    if args.outdir is not None:
        overrides["outdir"] = str(args.outdir)

    cfg = load_config(args.config, overrides=overrides)
    cfg.paths.outdir.mkdir(parents=True, exist_ok=True)
    _write_config_snapshots(cfg.paths.outdir, args.config, cfg)
    get_logger(cfg.paths.outdir / "tacticbench.log")
    configure_community_display_mapping_from_manifest(
        cfg.paths.outdir,
        figure_exclude_patterns=list(getattr(cfg.figures, "exclude_community_name_patterns", []) or []),
    )
    return cfg


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(prog="tacticbench", description="TACTIC benchmarking CLI")
    sub = ap.add_subparsers(dest="cmd", required=True)

    commands = [
        ("discover", discover.run),
        ("prepare-expected", prepare_expected.run),
        ("prepare-methods", prepare_methods.run),
        ("build-tree-full", build_tree_full.run),
        ("map-features", map_features.run),
        ("build-tree-centroid", build_tree_centroid.run),
        ("project-abundances", project_abundances.run),
        ("compute-distances", compute_distances.run),
        ("plot-fig01", plot_fig01.run),
        ("plot-fig02", plot_fig02.run),
        ("plot-fig03", plot_fig03.run),
        ("plot-fig04", plot_fig04.run),
        ("plot-fig05", plot_fig05.run),
        ("plot-fig06", plot_fig06.run),
        ("plot-fig07", plot_fig07.run),
        ("plot-fig08", plot_fig08.run),
        ("plot-fig09", plot_fig09.run),
        ("plot-fig10", plot_fig10.run),
        ("plot-fig11", plot_fig11.run),
        ("plot-fig12", plot_fig12.run),
        ("plot-fig13", plot_fig13.run),
        ("plot-fig14", plot_fig14.run),
        ("plot-fig15", plot_fig15.run),
        ("plot-fig16", plot_fig16.run),
        ("plot-fig17", plot_fig17.run),
        ("plot-fig18", plot_fig18.run),
        ("plot-fig19", plot_fig19.run),
        ("plot-fig20", plot_fig20.run),
        ("plot-fig21", plot_fig21.run),
        ("plot-all", plot_all.run),
        ("run-all", run_all.run),
    ]

    for name, fn in commands:
        sp = sub.add_parser(name)
        _add_common(sp)
        sp.set_defaults(_fn=fn)

    args = ap.parse_args(argv)
    cfg = _cfg_from_args(args)
    return int(args._fn(cfg))


if __name__ == "__main__":
    raise SystemExit(main())
