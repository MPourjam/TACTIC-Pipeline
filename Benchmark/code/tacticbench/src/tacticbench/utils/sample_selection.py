from __future__ import annotations

from dataclasses import dataclass
import re

from .labels import display_community_id, display_sample_id


@dataclass(frozen=True)
class SampleSelection:
    selected_original: list[str]
    label_by_original: dict[str, str]
    profile_by_original: dict[str, str]


def _sample_profile(sample_id: str) -> str:
    """
    Coarse composition profile inferred from sample naming tokens.
    """
    tokens = [t for t in re.split(r"[^A-Za-z0-9]+", str(sample_id).lower()) if t]
    for tok in tokens:
        if tok.startswith("staggered"):
            return "staggered"
        if tok.startswith("even"):
            return "even"
    return "other"


def _dedupe_labels(labels: list[str]) -> list[str]:
    """
    Ensure labels are unique while preserving order.
    """
    seen: dict[str, int] = {}
    out: list[str] = []
    for raw in labels:
        base = str(raw).strip() or "Sample"
        n = seen.get(base, 0) + 1
        seen[base] = n
        out.append(base if n == 1 else f"{base}_{n}")
    return out


def select_and_label_samples(
    community_id: str,
    sample_columns: list[str],
    keep_one_sample_per_community: bool,
) -> SampleSelection:
    """
    Select representative expected sample columns and assign downstream labels.

    Selection behavior when keep_one_sample_per_community=True:
      - If both Even and Staggered profiles exist, keep one of each (first by order).
      - Otherwise keep one sample (first by order).

    Labeling behavior:
      - If one sample selected: label as display community name (e.g., "Mock 8").
      - If multiple samples selected: label each as compact display sample id
        (e.g., "Even1", "Staggered1", "SRR123...").
    """
    ordered = [str(c) for c in sample_columns]
    if not ordered:
        return SampleSelection(selected_original=[], label_by_original={}, profile_by_original={})

    profile_by_original = {s: _sample_profile(s) for s in ordered}
    selected = ordered

    if keep_one_sample_per_community and len(ordered) > 1:
        even = [s for s in ordered if profile_by_original.get(s) == "even"]
        staggered = [s for s in ordered if profile_by_original.get(s) == "staggered"]
        if even and staggered:
            selected = [even[0], staggered[0]]
        else:
            selected = [ordered[0]]

    label_by_original: dict[str, str] = {}
    if len(selected) == 1:
        label_by_original[selected[0]] = display_community_id(str(community_id))
    else:
        compact_labels = _dedupe_labels([display_sample_id(s) for s in selected])
        label_by_original = dict(zip(selected, compact_labels))

    selected_profile = {s: profile_by_original.get(s, "other") for s in selected}
    return SampleSelection(
        selected_original=selected,
        label_by_original=label_by_original,
        profile_by_original=selected_profile,
    )

