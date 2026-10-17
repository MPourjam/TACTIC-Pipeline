from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Dict, List, Set, Tuple


@dataclass
class SampleMapResult:
    mapping: Dict[str, str]                 # expected_sample -> inferred_sample
    missing: List[str]                      # expected samples not found in inferred
    ambiguous: List[Tuple[str, List[str]]]  # expected sample -> list of possible inferred matches


_GENERIC_SAMPLE_TOKENS = {
    "mock",
    "sample",
    "trimmed",
    "untrimmed",
    "fastq",
    "gz",
    "r1",
    "r2",
    "l001",
}


def _normalize_name(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(s).lower())


def _sample_tokens(s: str) -> Set[str]:
    """
    Extract informative sample-name tokens (lowercase).
    Keeps alphanumeric chunks such as even1/staggered2/extreme.
    """
    cleaned = re.sub(r"[^A-Za-z0-9]+", " ", str(s).lower())
    toks = set(re.findall(r"[a-z]+[0-9]*", cleaned))
    return {t for t in toks if t and t not in _GENERIC_SAMPLE_TOKENS and len(t) > 1}


def _fallback_candidates(expected_name: str, inferred: List[str], used: Set[str]) -> List[Tuple[Tuple[int, int, int, int], str]]:
    """
    Build fallback candidates with score:
      (has_norm_containment, overlap_count, max_overlap_token_len, -length_delta)
    Higher is better.
    """
    out: List[Tuple[Tuple[int, int, int, int], str]] = []
    e_norm = _normalize_name(expected_name)
    e_toks = _sample_tokens(expected_name)
    for inf in inferred:
        if inf in used:
            continue
        i_norm = _normalize_name(inf)
        i_toks = _sample_tokens(inf)
        overlap = e_toks & i_toks
        has_containment = int(bool(e_norm) and bool(i_norm) and (e_norm in i_norm or i_norm in e_norm))
        overlap_count = len(overlap)
        max_overlap_len = max((len(t) for t in overlap), default=0)
        if not has_containment and overlap_count == 0:
            continue
        score = (
            has_containment,
            overlap_count,
            max_overlap_len,
            -abs(len(i_norm) - len(e_norm)),
        )
        out.append((score, inf))
    out.sort(key=lambda x: (x[0], -len(x[1]), x[1]), reverse=True)
    return out


def map_by_substring(expected: List[str], inferred: List[str]) -> SampleMapResult:
    """
    Map expected sample names to inferred sample names using substring rule:
      expected_name must be a substring of inferred_name.

    If multiple inferred names match, pick the "closest" by:
      1) minimal extra length (len(inferred)-len(expected))
      2) then shorter inferred name
      3) then lexicographic

    We record ambiguities when multiple inferred samples tie on the scoring.
    """
    mapping: Dict[str, str] = {}
    missing: List[str] = []
    ambiguous: List[Tuple[str, List[str]]] = []
    unresolved: List[str] = []

    # 1) Preserve historical behavior first.
    for e in expected:
        hits = [i for i in inferred if e in i]
        if not hits:
            unresolved.append(e)
            continue

        hits = sorted(hits, key=lambda x: (len(x) - len(e), len(x), x))
        best = hits[0]

        best_key = (len(best) - len(e), len(best))
        tied = [h for h in hits if (len(h) - len(e), len(h)) == best_key]
        if len(tied) > 1:
            ambiguous.append((e, tied))

        mapping[e] = best

    # 2) Robust fallback for naming-variant communities (token/normalized matching).
    used = set(mapping.values())
    still_unresolved: List[str] = []
    for e in unresolved:
        cands = _fallback_candidates(e, inferred, used=used)
        if not cands:
            still_unresolved.append(e)
            continue

        best_score, best_inf = cands[0]
        ties = [inf for score, inf in cands if score == best_score]
        if len(ties) > 1:
            ambiguous.append((e, ties))

        mapping[e] = best_inf
        used.add(best_inf)

    # 3) Last-resort deterministic 1:1 fill for a single unresolved pair.
    if len(still_unresolved) == 1:
        remaining = [i for i in inferred if i not in used]
        if len(remaining) == 1:
            mapping[still_unresolved[0]] = remaining[0]
            used.add(remaining[0])
            still_unresolved = []

    missing.extend(still_unresolved)
    return SampleMapResult(mapping=mapping, missing=missing, ambiguous=ambiguous)
