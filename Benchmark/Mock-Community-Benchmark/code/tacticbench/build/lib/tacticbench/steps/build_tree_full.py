from __future__ import annotations

from pathlib import Path
import subprocess
import shlex

import pandas as pd

from ..log import get_logger
from ..utils.paths import expected_dir, qc_dir


def run(cfg) -> int:
    """
    Build the full expected-tip tree per community from expected sequences.

    Inputs:
      outdir/communities/<cid>/expected/expected_sequences.fasta

    Outputs:
      outdir/communities/<cid>/expected/expected_tree_full.nwk
      outdir/communities/<cid>/expected/expected_sequences.aln.fasta

    Uses external tools:
      - MAFFT
      - FastTree (or fasttree)

    If tools are missing or fail, the pipeline continues and UniFrac/gUniFrac become NaN.
    """
    log = get_logger()
    outdir = Path(cfg.paths.outdir)
    manifest = pd.read_csv(outdir / "manifest.tsv", sep="\t")

    skipped = set()
    sp = qc_dir(outdir) / "skipped_communities.tsv"
    if sp.exists():
        try:
            skipped = set(pd.read_csv(sp, sep="\t")["community_id"].astype(str).tolist())
        except Exception:
            skipped = set()

    mafft = cfg.tree.mafft_exe
    fasttree = cfg.tree.fasttree_exe

    for _, r in manifest.iterrows():
        cid = str(r["community_id"])
        if cid in skipped:
            continue

        ed = expected_dir(outdir, cid)
        fasta = ed / "expected_sequences.fasta"
        if not fasta.exists():
            continue

        aln = ed / "expected_sequences.aln.fasta"
        tree = ed / "expected_tree_full.nwk"
        if tree.exists() and tree.stat().st_size > 0:
            continue

        try:
            # MAFFT alignment
            with open(aln, "w", encoding="utf-8") as out:
                subprocess.run(
                    [mafft, "--auto", str(fasta)],
                    check=True,
                    stdout=out,
                    stderr=subprocess.PIPE,
                    text=True,
                )

            # FastTree
            with open(tree, "w", encoding="utf-8") as out:
                subprocess.run(
                    [fasttree, "-nt", str(aln)],
                    check=True,
                    stdout=out,
                    stderr=subprocess.PIPE,
                    text=True,
                )

            log.info(f"[{cid}] built full tree: {tree}")
        except FileNotFoundError:
            log.info(f"[{cid}] missing MAFFT/FastTree; skipping full tree.")
        except subprocess.CalledProcessError as e:
            stderr = (e.stderr or "").strip()
            cmd = " ".join(shlex.quote(str(x)) for x in (e.cmd or []))
            if stderr:
                log.info(f"[{cid}] full tree build failed; cmd={cmd}; stderr={stderr}")
            else:
                log.info(f"[{cid}] full tree build failed; cmd={cmd}; no stderr captured.")

    return 0
