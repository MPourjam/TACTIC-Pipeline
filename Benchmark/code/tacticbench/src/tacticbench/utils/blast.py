from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
import subprocess
import random


@dataclass(frozen=True)
class BlastTools:
    makeblastdb: str
    blastn: str


def make_db(tools: BlastTools, fasta: Path, db_prefix: Path) -> bool:
    """
    Create a BLAST nucleotide DB from fasta.
    Returns True on success, False otherwise.
    """
    try:
        subprocess.run(
            [tools.makeblastdb, "-in", str(fasta), "-dbtype", "nucl", "-out", str(db_prefix)],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        return True
    except Exception:
        return False


def run_blast(
    tools: BlastTools,
    query_fasta: Path,
    db_prefix: Path,
    out_tsv: Path,
    min_ident: float,
    max_targets: int = 200,
) -> bool:
    """
    Run blastn query vs db_prefix and write outfmt 6 TSV:
      qseqid sseqid pident length bitscore evalue qlen qcovhsp qstart qend
    Returns True on success, False otherwise.
    """
    try:
        subprocess.run(
            [
                tools.blastn,
                "-query",
                str(query_fasta),
                "-db",
                str(db_prefix),
                "-outfmt",
                "6 qseqid sseqid pident length bitscore evalue qlen qcovhsp qstart qend",
                "-max_target_seqs",
                str(int(max_targets)),
                "-perc_identity",
                str(float(min_ident)),
                "-out",
                str(out_tsv),
            ],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        return True
    except Exception:
        return False


def read_hits(path: Path) -> pd.DataFrame:
    """
    Read BLAST outfmt6 with columns:
      qseqid sseqid pident length bitscore evalue qlen qcovhsp qstart qend
    Strict mode: requires coordinate fields (qstart/qend) so effective identity
    is always computed from exact query coverage.
    """
    raw = pd.read_csv(
        path,
        sep="\t",
        header=None,
    )
    if raw.empty:
        return pd.DataFrame(
            columns=[
                "qseqid",
                "sseqid",
                "pident",
                "length",
                "bitscore",
                "evalue",
                "qlen",
                "qcovhsp",
                "qstart",
                "qend",
                "qcov_exact",
            ]
        )

    if raw.shape[1] >= 10:
        raw = raw.iloc[:, :10].copy()
        raw.columns = ["qseqid", "sseqid", "pident", "length", "bitscore", "evalue", "qlen", "qcovhsp", "qstart", "qend"]
    else:
        raise ValueError(
            f"BLAST hits file {path} is missing required qstart/qend fields "
            "(expected >=10 outfmt columns). Re-run `map-features` to regenerate blast_hits.tsv."
        )

    raw["qseqid"] = raw["qseqid"].astype(str)
    raw["sseqid"] = raw["sseqid"].astype(str)
    for c in ("pident", "length", "bitscore", "evalue", "qlen", "qcovhsp", "qstart", "qend"):
        raw[c] = pd.to_numeric(raw[c], errors="coerce")

    # Exact query coverage fraction based on query coordinates.
    # qcov_exact is stored as fraction in [0,1].
    with np.errstate(divide="ignore", invalid="ignore"):
        qspan = (raw["qend"] - raw["qstart"]).abs() + 1.0
        qcov_exact = qspan / raw["qlen"].replace(0.0, np.nan)
    raw["qcov_exact"] = pd.to_numeric(qcov_exact, errors="coerce").clip(lower=0.0, upper=1.0)

    # Keep qcovhsp synced with exact coordinate-based coverage when missing.
    with np.errstate(divide="ignore", invalid="ignore"):
        est_qcov = 100.0 * raw["qcov_exact"]
    raw["qcovhsp"] = raw["qcovhsp"].fillna(est_qcov)
    _ensure_effective_identity(raw)

    return raw


def _ensure_effective_identity(hits: pd.DataFrame) -> pd.DataFrame:
    """
    Ensure an `effective_identity` column exists:
      effective_identity = pident * qcov_exact
    where qcov_exact = (abs(qend - qstart) + 1) / qlen.
    This function intentionally has no legacy fallback path: qstart/qend/qlen
    are required so malformed or legacy files fail fast.
    """
    if hits.empty:
        if "effective_identity" not in hits.columns:
            hits["effective_identity"] = np.nan
        return hits

    required = {"pident", "qstart", "qend", "qlen"}
    missing = [c for c in sorted(required) if c not in hits.columns]
    if missing:
        raise ValueError(f"Missing required BLAST columns for effective_identity: {', '.join(missing)}")

    pid = pd.to_numeric(hits["pident"], errors="coerce")
    qstart = pd.to_numeric(hits["qstart"], errors="coerce")
    qend = pd.to_numeric(hits["qend"], errors="coerce")
    qlen = pd.to_numeric(hits["qlen"], errors="coerce")

    invalid = pid.isna() | qstart.isna() | qend.isna() | qlen.isna() | (qlen <= 0.0)
    if invalid.any():
        n_bad = int(invalid.sum())
        raise ValueError(
            f"Found {n_bad} BLAST rows with invalid pident/qstart/qend/qlen values; "
            "cannot compute strict effective_identity."
        )

    with np.errstate(divide="ignore", invalid="ignore"):
        qspan = (qend - qstart).abs() + 1.0
        qcov_exact = qspan / qlen
    qcov_exact = qcov_exact.clip(lower=0.0, upper=1.0)
    hits["qcov_exact"] = qcov_exact

    eff = (pid * qcov_exact).clip(lower=0.0, upper=100.0)
    hits["effective_identity"] = eff

    if "qcovhsp" not in hits.columns:
        hits["qcovhsp"] = np.nan
    qcov_from_exact = 100.0 * qcov_exact
    hits["qcovhsp"] = pd.to_numeric(hits["qcovhsp"], errors="coerce").fillna(qcov_from_exact)
    hits["qcovhsp"] = hits["qcovhsp"].clip(lower=0.0, upper=100.0)
    return hits


def filter_hits(
    hits: pd.DataFrame,
    min_effective_pident: float | None = None,
    min_pident: float | None = None,
    min_qcovhsp: float = 0.0,
    strict_pident_gt: bool = False,
) -> pd.DataFrame:
    """
    Apply hit-quality filters in a consistent way.
    - min_effective_pident: threshold on effective identity
      (effective_identity = pident * qcov_exact).
    - min_pident: legacy alias for min_effective_pident.
    - min_qcovhsp: minimum query coverage (percent or fraction in [0,1]),
      evaluated strictly on coordinate-derived qcov_exact.
    """
    if hits.empty:
        return hits.copy()

    out = hits.copy()
    for c in ("pident", "length", "qcovhsp", "qcov_exact", "qlen"):
        if c in out.columns:
            out[c] = pd.to_numeric(out[c], errors="coerce")

    # Ensure strict coordinate-based effective identity is present.
    _ensure_effective_identity(out)

    if min_effective_pident is None and min_pident is not None:
        min_effective_pident = min_pident
    if min_effective_pident is not None:
        mp = float(min_effective_pident)
        metric = pd.to_numeric(out.get("effective_identity", out.get("pident", np.nan)), errors="coerce")
        if strict_pident_gt:
            out = out[metric > mp].copy()
        else:
            out = out[metric >= mp].copy()

    qcov = float(min_qcovhsp or 0.0)
    if qcov > 1.0:
        qcov /= 100.0
    if qcov > 0.0:
        if "qcov_exact" not in out.columns:
            raise ValueError("Missing required BLAST column qcov_exact for coverage filtering.")
        cov = pd.to_numeric(out["qcov_exact"], errors="coerce")
        if cov.isna().any():
            n_bad = int(cov.isna().sum())
            raise ValueError(f"Found {n_bad} BLAST rows with invalid qcov_exact; cannot apply strict coverage filter.")
        out = out[cov >= qcov].copy()

    return out


def _score_tuple(rec: dict) -> tuple:
    """
    Higher is better.
    Order: effective_identity, pident, bitscore, -evalue, length
    """
    eff = rec.get("effective_identity", rec.get("pident", np.nan))
    return (
        float(eff),
        float(rec["pident"]),
        float(rec["bitscore"]),
        -float(rec["evalue"]),
        float(rec["length"]),
    )


def stable_one_to_one(
    candidates: Dict[str, List[dict]],
    rng: random.Random,
) -> Dict[str, dict]:
    """
    Deterministic-ish (seeded) stable greedy assignment enforcing 1-to-1 targets.

    candidates:
      q -> list of candidate dicts, already sorted best->worst.
      Each candidate dict MUST include:
        - "target": str (tip id OR species id)
        - "effective_identity" (optional; falls back to pident)
        - "pident", "bitscore", "evalue", "length"

    Algorithm:
      - each query proposes its best remaining candidate
      - if target free: assign
      - if target taken: keep the query with higher score; loser advances to next candidate and re-queues
      - ties are broken by rng (seeded)
    """
    next_idx = {q: 0 for q in candidates}
    assigned_target: Dict[str, tuple[str, dict]] = {}  # target -> (q, rec)
    assigned_q: Dict[str, dict] = {}

    queue = list(candidates.keys())
    while queue:
        q = queue.pop(0)
        opts = candidates.get(q, [])
        while next_idx[q] < len(opts):
            rec = opts[next_idx[q]]
            tgt = rec["target"]

            if tgt not in assigned_target:
                assigned_target[tgt] = (q, rec)
                assigned_q[q] = rec
                break

            q2, rec2 = assigned_target[tgt]
            s1, s2 = _score_tuple(rec), _score_tuple(rec2)

            if s1 > s2 or (s1 == s2 and rng.random() < 0.5):
                # q takes tgt, q2 gets kicked and must try its next option
                assigned_target[tgt] = (q, rec)
                assigned_q[q] = rec

                if q2 in assigned_q:
                    del assigned_q[q2]
                next_idx[q2] += 1
                queue.append(q2)
                break

            # q loses: try next candidate
            next_idx[q] += 1

    return assigned_q
