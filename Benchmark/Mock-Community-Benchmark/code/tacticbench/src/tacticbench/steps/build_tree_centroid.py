from __future__ import annotations

from pathlib import Path
import subprocess

import pandas as pd

from ..log import get_logger
from ..utils.io import load_fasta, write_fasta
from ..utils.paths import expected_dir, qc_dir


def _seq_distance(a: str, b: str) -> float:
    """
    Normalized mismatch distance between two unaligned sequences in [0, 1].
    """
    la = len(a)
    lb = len(b)
    if la == 0 and lb == 0:
        return 0.0
    n = min(la, lb)
    mism = 0
    if n > 0:
        mism = sum(1 for i in range(n) if a[i] != b[i])
    mism += abs(la - lb)
    den = max(la, lb)
    return float(mism / den) if den > 0 else 0.0


def _pick_medoid_tip(tips: list[str], seqs: dict[str, str]) -> str:
    """
    Pick the sequence medoid tip: minimum average pairwise distance to all tips.
    """
    valid = [t for t in tips if t in seqs]
    if not valid:
        return ""
    if len(valid) == 1:
        return valid[0]

    best_tip = valid[0]
    best_score = float("inf")
    for t in valid:
        ds = [_seq_distance(seqs[t], seqs[u]) for u in valid if u != t]
        score = float(sum(ds) / len(ds)) if ds else 0.0
        if score < best_score:
            best_score = score
            best_tip = t
    return best_tip


def run(cfg) -> int:
    """
    Build centroid species tree per community for clustering-based methods (TIC/TAC/OTU).

    Rule for selecting one expected-tip sequence per species:
      - Use the expected-sequence medoid for each species (tip with minimum average
        pairwise sequence distance to other tips of that species).

    Inputs:
      - expected sequences: outdir/communities/<cid>/expected/expected_sequences.fasta
      - tip_to_species:     outdir/communities/<cid>/expected/tip_to_species.tsv
      - species list:       outdir/communities/<cid>/expected/taxonomy_collapsed_species.tsv

    Outputs:
      - centroid fasta: outdir/communities/<cid>/expected/centroid_species_sequences.fasta
      - alignment:     outdir/communities/<cid>/expected/centroid_species_sequences.aln.fasta
      - tree:          outdir/communities/<cid>/expected/centroid_species_tree.nwk

    Uses external tools:
      - MAFFT
      - FastTree

    If tools are missing or fail, the pipeline continues (UniFrac/gUniFrac will be NaN for clustering methods).
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
        seq_p = ed / "expected_sequences.fasta"
        tipmap_p = ed / "tip_to_species.tsv"
        sp_p = ed / "taxonomy_collapsed_species.tsv"
        if not seq_p.exists() or not tipmap_p.exists() or not sp_p.exists():
            continue

        centroid_fa = ed / "centroid_species_sequences.fasta"
        aln = ed / "centroid_species_sequences.aln.fasta"
        tree = ed / "centroid_species_tree.nwk"
        if tree.exists() and tree.stat().st_size > 0:
            continue

        species = pd.read_csv(sp_p, sep="\t", index_col=0).index.astype(str).tolist()
        seqs = load_fasta(seq_p)

        tipmap = pd.read_csv(tipmap_p, sep="\t")
        tip_to_sp = dict(zip(tipmap["tip"].astype(str), tipmap["species"].astype(str)))

        # species -> list of tips
        sp_to_tips: dict[str, list[str]] = {}
        for tip, spx in tip_to_sp.items():
            sp_to_tips.setdefault(spx, []).append(tip)

        centroid_seqs: dict[str, str] = {}
        for spx in species:
            tips = sp_to_tips.get(spx, [])
            chosen_tip = _pick_medoid_tip(tips, seqs)
            if chosen_tip and chosen_tip in seqs:
                centroid_seqs[spx] = seqs[chosen_tip]

        if not centroid_seqs:
            log.info(f"[{cid}] no centroid sequences created; skipping centroid tree.")
            continue

        write_fasta(centroid_seqs, centroid_fa)

        try:
            with open(aln, "w", encoding="utf-8") as out:
                subprocess.run([mafft, "--auto", str(centroid_fa)], check=True, stdout=out)

            with open(tree, "w", encoding="utf-8") as out:
                subprocess.run([fasttree, "-nt", str(aln)], check=True, stdout=out)

            log.info(f"[{cid}] built centroid species tree: {tree}")
        except FileNotFoundError:
            log.info(f"[{cid}] missing MAFFT/FastTree; skipping centroid tree.")
        except subprocess.CalledProcessError:
            log.info(f"[{cid}] centroid tree build failed; skipping.")

    return 0
