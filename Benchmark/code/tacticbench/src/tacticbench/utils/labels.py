from __future__ import annotations

import csv
from fnmatch import fnmatch
from pathlib import Path
import re

METHOD_DISPLAY = {
    "ASV": "DADA2",
    "zOTU": "UNOISE3",
    "OTU": "UPARSE",
    "TIC": "TIC",
    "TAC": "TAC",
}

COMMUNITY_DISPLAY = {
    # Common publication-friendly aliases
    "ZIEL_1": "ZIEL I",
    "ZIEL_2": "ZIEL II",
    "mock-ZIEL_1": "ZIEL I",
    "mock-ZIEL_2": "ZIEL II",
    "mock-ZIEL1": "ZIEL I",
    "mock-ZIEL2": "ZIEL II",
    "mock-ZYMO": "ZYMO",
}

SAMPLE_DISPLAY = {
    # Common non-SRA mock labels
    "mock.1": "Sample1",
    "mock-1": "Sample1",
    "mock_1": "Sample1",
}

_COMMUNITY_DYNAMIC_DISPLAY: dict[str, str] = {}


def display_method(method: str) -> str:
    return METHOD_DISPLAY.get(str(method), str(method))


def display_methods(methods: list[str]) -> list[str]:
    return [display_method(m) for m in methods]


def _community_sort_key(cid: str) -> tuple[int, int | str]:
    s = str(cid)
    m = re.match(r"^mock-(\d+)$", s)
    if m:
        return (0, int(m.group(1)))
    return (1, s.lower())


def _is_excluded(cid: str, patterns: list[str]) -> bool:
    for pat in patterns:
        p = str(pat).strip()
        if not p:
            continue
        if fnmatch(cid, p) or (p in cid):
            return True
    return False


def build_community_display_mapping(community_ids: list[str], figure_exclude_patterns: list[str] | None = None) -> dict[str, str]:
    """
    Build display labels with continuous Mock numbering for communities
    included in figures after applying global figure exclusion patterns.
    """
    pats = [str(p) for p in (figure_exclude_patterns or [])]
    ordered = sorted([str(c) for c in community_ids], key=_community_sort_key)
    included = {cid for cid in ordered if not _is_excluded(cid, pats)}

    out: dict[str, str] = {}
    mock_counter = 1

    for cid in ordered:
        if cid not in included:
            continue
        if cid in COMMUNITY_DISPLAY:
            out[cid] = COMMUNITY_DISPLAY[cid]
            continue
        m = re.match(r"^mock-(\d+)$", cid)
        if m:
            out[cid] = f"Mock {mock_counter}"
            mock_counter += 1
        else:
            out[cid] = cid

    # Give excluded communities a stable fallback label without consuming
    # contiguous figure numbering.
    for cid in ordered:
        if cid in out:
            continue
        if cid in COMMUNITY_DISPLAY:
            out[cid] = COMMUNITY_DISPLAY[cid]
            continue
        m = re.match(r"^mock-(\d+)$", cid)
        if m:
            out[cid] = f"Mock {int(m.group(1))}"
        else:
            out[cid] = cid

    return out


def set_community_display_mapping(mapping: dict[str, str] | None) -> None:
    global _COMMUNITY_DYNAMIC_DISPLAY
    _COMMUNITY_DYNAMIC_DISPLAY = dict(mapping or {})


def configure_community_display_mapping_from_manifest(
    outdir: Path, figure_exclude_patterns: list[str] | None = None
) -> dict[str, str]:
    """
    Load community ids from outdir/manifest.tsv, build continuous display mapping
    using global figure exclusions, and set it as active runtime mapping.
    """
    manifest = Path(outdir) / "manifest.tsv"
    if not manifest.exists():
        set_community_display_mapping({})
        return {}

    ids: list[str] = []
    try:
        with manifest.open("r", encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                cid = str(row.get("community_id", "")).strip()
                if cid:
                    ids.append(cid)
    except Exception:
        set_community_display_mapping({})
        return {}

    mapping = build_community_display_mapping(ids, figure_exclude_patterns=figure_exclude_patterns)
    set_community_display_mapping(mapping)
    try:
        write_community_display_mapping_table(Path(outdir) / "community_name_mapping.tsv", mapping)
    except Exception:
        pass
    return mapping


def write_community_display_mapping_table(path: Path, mapping: dict[str, str]) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["original_community_id", "display_community_name"])
        for cid in sorted(mapping.keys(), key=_community_sort_key):
            w.writerow([cid, mapping[cid]])


def display_community_id(community_id: str) -> str:
    s = str(community_id)
    if s in _COMMUNITY_DYNAMIC_DISPLAY:
        return _COMMUNITY_DYNAMIC_DISPLAY[s]
    if s in COMMUNITY_DISPLAY:
        return COMMUNITY_DISPLAY[s]
    m = re.match(r"^mock-(\d+)$", s)
    if m:
        return f"Mock {int(m.group(1))}"
    return s


def display_sample_id(sample_id: str) -> str:
    s = str(sample_id).strip()
    if s in SAMPLE_DISPLAY:
        return SAMPLE_DISPLAY[s]

    # Keep already-formatted mock labels stable.
    m = re.match(r"^mock\s+(\d+)$", s, flags=re.IGNORECASE)
    if m:
        return f"Mock {int(m.group(1))}"

    # Keep SRA-style accessions unchanged
    if re.match(r"^(?:SRR|ERR|DRR)\d+$", s):
        return s

    # Compact common mock patterns
    m = re.search(r"(Even\d+|Staggered\d+|Extreme\.?\d+)$", s, flags=re.IGNORECASE)
    if m:
        token = m.group(1).replace(".", "")
        if token.lower().startswith("even"):
            return f"Even{token[4:]}"
        if token.lower().startswith("staggered"):
            return f"Staggered{token[9:]}"
        if token.lower().startswith("extreme"):
            return f"Extreme{token[7:]}"

    m = re.match(r"^mock[._-]?(\d+)$", s, flags=re.IGNORECASE)
    if m:
        return f"Sample{int(m.group(1))}"

    # Generic clean fallback
    out = re.sub(r"[^A-Za-z0-9]+", "_", s).strip("_")
    return out or s


def display_sample_key(sample_key: str) -> str:
    """
    sample_key is typically community_id|sample_id. Keep sample_id as-is, but
    display the community_id in a publication-friendly form.
    """
    s = str(sample_key)
    if "|" not in s:
        return display_community_id(s)
    cid, sid = s.split("|", 1)
    disp_c = display_community_id(cid)

    # If sample id is effectively the community label, omit the sample suffix.
    if str(sid).strip() in {"", str(cid).strip(), disp_c}:
        return disp_c
    disp_s = display_sample_id(sid)
    if disp_s == disp_c:
        return disp_c
    return f"{disp_c}|{disp_s}"
