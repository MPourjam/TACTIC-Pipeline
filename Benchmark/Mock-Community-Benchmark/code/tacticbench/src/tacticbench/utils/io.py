from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd
from Bio import SeqIO


FASTA_EXTS = {".fa", ".fasta", ".fna", ".fas"}
TAB_EXTS = {".tab", ".tsv", ".txt"}


def read_expected_taxonomy(path: Path) -> pd.DataFrame:
    """
    Read expected composition table (source/taxonomy.tsv).

    Assumptions:
      - TSV with header
      - first column contains expected feature/taxon labels
      - remaining columns are samples (mock community name / SRA accession / etc.)
    Returns:
      - DataFrame indexed by Taxon labels, columns = samples, values=float (0 if missing)
    """
    df = pd.read_csv(path, sep="\t", header=0, dtype={0: str})
    if df.shape[1] < 2:
        raise ValueError(f"Expected >=2 columns in taxonomy.tsv: {path}")
    df = df.rename(columns={df.columns[0]: "Taxon"})
    df["Taxon"] = df["Taxon"].astype(str)
    df = df.set_index("Taxon")
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return df


def canonical_species_name(label: str) -> str:
    """
    Canonicalize species labels as: Genus_species
      - genus token: Capitalized
      - species token: lowercase
      - separator: underscore
    """
    s = str(label).strip()
    s = re.sub(r"[^A-Za-z0-9]+", "_", s)
    parts = [p for p in s.split("_") if p]
    if len(parts) >= 2:
        genus = parts[0].capitalize()
        species = parts[1].lower()
        return f"{genus}_{species}"
    if len(parts) == 1:
        return parts[0].capitalize()
    return s


def genus_species(label: str) -> str:
    """
    Collapse label to genus_species and canonicalize case/format.
    """
    return canonical_species_name(label)


def collapse_to_species(df: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse expected taxonomy table to genus_species by summing rows.
    """
    tmp = df.copy()
    tmp["__gs__"] = [genus_species(i) for i in tmp.index]
    out = tmp.groupby("__gs__", sort=False).sum(numeric_only=True)
    out.index.name = "genus_species"
    return out


def read_dada2_asv_counts(path: Path) -> pd.DataFrame:
    """
    Read DADA2 ASV_counts.tsv:
      - rows = samples
      - columns = ASVs
    Returns:
      - features x samples (ASVs as rows, expected-sample-names as columns after mapping)
    """
    df = pd.read_csv(path, sep="\t", header=0, dtype={0: str})
    df = df.rename(columns={df.columns[0]: "Sample"})
    df["Sample"] = df["Sample"].astype(str)
    df = df.set_index("Sample")
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return df.T


def read_tactic_table(path: Path) -> pd.DataFrame:
    """
    Read TACTIC tables:
      - first column: Feature ID
      - remaining columns: samples
      - last column may be 'Taxonomy' (case-insensitive) and should be dropped
    Returns:
      - features x samples
    """
    df = pd.read_csv(path, sep="\t", header=0, dtype=str).fillna("")
    if df.shape[1] < 2:
        raise ValueError(f"Unexpected table shape (<2 columns): {path}")

    # drop last column if it's taxonomy
    last = str(df.columns[-1]).strip().lower()
    if last == "taxonomy":
        df = df.iloc[:, :-1]

    df = df.rename(columns={df.columns[0]: "Feature"})
    df["Feature"] = df["Feature"].astype(str)
    df = df.set_index("Feature")

    df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return df


def fasta_ids(path: Path) -> List[str]:
    """
    Return FASTA record IDs (rec.id) in order.
    """
    return [rec.id for rec in SeqIO.parse(str(path), "fasta")]


def load_fasta(path: Path) -> Dict[str, str]:
    """
    Load FASTA into {id: sequence} (uppercase).
    """
    out: Dict[str, str] = {}
    for rec in SeqIO.parse(str(path), "fasta"):
        out[rec.id] = str(rec.seq).upper()
    return out


def write_fasta(seqs: Dict[str, str], path: Path) -> None:
    """
    Write {id: sequence} to FASTA.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        for k, v in seqs.items():
            f.write(f">{k}\n{v}\n")


def load_fasta_records(path: Path) -> List[Tuple[str, str]]:
    """
    Load FASTA preserving record order and duplicate IDs.
    Returns list of (id, sequence_upper).
    """
    out: List[Tuple[str, str]] = []
    for rec in SeqIO.parse(str(path), "fasta"):
        out.append((str(rec.id), str(rec.seq).upper()))
    return out


def write_fasta_records(records: List[Tuple[str, str]], path: Path) -> None:
    """
    Write FASTA records from a list of (id, sequence), preserving order.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        for rid, seq in records:
            f.write(f">{rid}\n{seq}\n")


def uniquify_fasta_ids(
    records: List[Tuple[str, str]],
    copy_suffix: str = "__Copy",
) -> Tuple[List[Tuple[str, str]], Dict[str, List[str]]]:
    """
    Make FASTA IDs unique by suffixing repeated IDs as:
      <id>__Copy1, <id>__Copy2, ...

    Returns:
      - unique records list (same order/sequence content)
      - map original_id -> assigned IDs in appearance order
    """
    used: set[str] = set()
    seen: Dict[str, int] = {}
    assigned_by_original: Dict[str, List[str]] = {}
    out: List[Tuple[str, str]] = []

    for rid, seq in records:
        base = str(rid)
        n = seen.get(base, 0)

        if n == 0 and base not in used:
            new_id = base
        else:
            k = max(1, n)
            new_id = f"{base}{copy_suffix}{k}"
            while new_id in used:
                k += 1
                new_id = f"{base}{copy_suffix}{k}"

        seen[base] = n + 1
        used.add(new_id)
        assigned_by_original.setdefault(base, []).append(new_id)
        out.append((new_id, seq))

    return out, assigned_by_original


def dedupe_fasta_by_sequence(
    records: List[Tuple[str, str]],
) -> Tuple[List[Tuple[str, str]], List[Dict[str, object]]]:
    """
    Deduplicate FASTA records by exact sequence string (uppercase).

    Keeps the first occurrence of each sequence and drops later duplicates.

    Returns:
      - deduplicated records list (order preserved by first occurrence)
      - list of collapsed records with keys:
          dropped_id, kept_id, sequence_length
    """
    seen: Dict[str, str] = {}
    out: List[Tuple[str, str]] = []
    collapsed: List[Dict[str, object]] = []

    for rid, seq in records:
        s = str(seq).upper()
        kept_id = seen.get(s)
        if kept_id is None:
            seen[s] = str(rid)
            out.append((str(rid), s))
            continue

        collapsed.append(
            {
                "dropped_id": str(rid),
                "kept_id": kept_id,
                "sequence_length": int(len(s)),
            }
        )

    return out, collapsed


def _looks_like_fasta(text: str) -> bool:
    return text.lstrip().startswith(">")


_DNA_RE = re.compile(r"^[ACGTUN\-\.]+$", re.IGNORECASE)


def _looks_like_dna(seq: str) -> bool:
    s = seq.strip()
    return bool(s) and bool(_DNA_RE.match(s))


def ensure_fasta(seq_path: Path, out_fasta: Path) -> Path:
    """
    Normalize method sequence file into FASTA at out_fasta.

    Supports:
      - FASTA (.fa/.fasta/.fna/.fas) -> copied
      - TAB/TSV (id<TAB>sequence[<TAB>...]) -> converted
      - If file content starts with '>' it's treated as FASTA regardless of extension.
    """
    out_fasta.parent.mkdir(parents=True, exist_ok=True)

    txt = seq_path.read_text(encoding="utf-8", errors="replace")
    if _looks_like_fasta(txt) or seq_path.suffix.lower() in FASTA_EXTS:
        out_fasta.write_text(txt, encoding="utf-8")
        return out_fasta

    # tab-like: id \t seq
    seqs: Dict[str, str] = {}
    for ln in txt.splitlines():
        if not ln.strip():
            continue
        parts = ln.split("\t")
        if len(parts) < 2:
            continue
        fid = parts[0].strip()
        seq = parts[1].strip().upper()
        # skip header-ish lines
        if fid.lower() in {"id", "feature", "otu", "zotu", "asv"} and not _looks_like_dna(seq):
            continue
        if not _looks_like_dna(seq):
            continue
        if fid:
            seqs[fid] = seq

    if not seqs:
        raise ValueError(f"Could not parse sequences as FASTA or TAB: {seq_path}")

    write_fasta(seqs, out_fasta)
    return out_fasta


def tip_to_species_map(expected_fasta: Path, expected_taxonomy_raw: pd.DataFrame) -> Dict[str, str]:
    """
    Map expected tip IDs (from expected-sequences.fasta) to genus_species.

    Notes:
      - expected-sequences.fasta may contain >= rows than taxonomy.tsv.
      - We try to map a tip id to a taxonomy row label:
          1) exact match to a taxonomy row label
          2) substring match (pick the longest matching taxonomy label)
          3) fallback to deriving genus_species from the tip id itself
    """
    tips = fasta_ids(expected_fasta)
    labels = list(expected_taxonomy_raw.index.astype(str))

    label_to_sp = {lbl: genus_species(lbl) for lbl in labels}

    def best_match(tid: str) -> str | None:
        if tid in label_to_sp:
            return tid
        hits = [lbl for lbl in labels if (lbl in tid) or (tid in lbl)]
        if not hits:
            return None
        # prefer longest label match (most specific)
        return sorted(hits, key=len, reverse=True)[0]

    out: Dict[str, str] = {}
    for tid in tips:
        m = best_match(tid)
        if m is not None:
            out[tid] = label_to_sp[m]
        else:
            out[tid] = genus_species(tid)
    return out


def species_to_tips(tip_to_sp: Dict[str, str]) -> Dict[str, List[str]]:
    """
    Invert tip->species map into species->tips.
    """
   
