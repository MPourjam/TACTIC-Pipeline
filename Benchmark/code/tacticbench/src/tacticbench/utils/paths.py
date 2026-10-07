from __future__ import annotations

from pathlib import Path


def qc_dir(outdir: Path) -> Path:
    return outdir / "qc"


def figures_root(outdir: Path) -> Path:
    return outdir / "figures"


def comm_dir(outdir: Path, community_id: str) -> Path:
    return outdir / "communities" / community_id


def expected_dir(outdir: Path, community_id: str) -> Path:
    return comm_dir(outdir, community_id) / "expected"


def method_dir(outdir: Path, community_id: str, method: str) -> Path:
    return comm_dir(outdir, community_id) / "methods" / method


def distances_dir(outdir: Path, community_id: str) -> Path:
    return comm_dir(outdir, community_id) / "distances"
