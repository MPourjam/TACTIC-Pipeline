import csv
import os
import re
from pathlib import Path, PurePath
import tempfile
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple, Union


_KINGDOM_ATI_RE = re.compile(r"^[A-Za-z]+ati(?:--other)?$", re.IGNORECASE)
_SILVA_PATH_VERSION_RE = re.compile(
    r"(?:^|[/\\_.-])silva[/\\_.-]+v?(?P<ver>\d+(?:\.\d+)*)(?:$|[/\\_.-])",
    re.IGNORECASE,
)
_EXPLICIT_LEVEL_SUFFIX_RE = re.compile(r"-[78]Levels\.(?:tab|fasta|html)$", re.IGNORECASE)
_FASTA_SEQ_ID_RE = re.compile(r"^>(?P<seq_id>[^\s;]+)", re.IGNORECASE)
_FASTA_TAX_RE = re.compile(r"(?P<prefix>^(?:>[^\n]*?[\s;])tax=)(?P<tax>[^\s]*)", re.IGNORECASE)


def _parse_version_tuple(version_str: str) -> Optional[Tuple[int, ...]]:
    cleaned = str(version_str).strip().lstrip("vV")
    if not cleaned or not re.fullmatch(r"\d+(?:\.\d+)*", cleaned):
        return None
    return tuple(int(part) for part in cleaned.split("."))


def detect_silva_version(dbs_loc: str) -> Optional[Tuple[int, ...]]:
    """
    Infer the SILVA database version tuple (e.g. (144,) or (138, 2)) from a
    database directory or file path.
    """
    if not dbs_loc:
        return None
    raw_str = str(dbs_loc).strip()
    match = _SILVA_PATH_VERSION_RE.search(raw_str)
    if match:
        parsed = _parse_version_tuple(match.group("ver"))
        if parsed is not None:
            return parsed

    path_obj = Path(PurePath(raw_str))
    parts = list(path_obj.parts)
    for idx, part in enumerate(parts):
        if part.lower() == "silva" and idx + 1 < len(parts):
            parsed = _parse_version_tuple(parts[idx + 1])
            if parsed is not None:
                return parsed
        parsed_part = _parse_version_tuple(part)
        if parsed_part is not None and idx > 0 and "silva" in parts[idx - 1].lower():
            return parsed_part

    # Check if directory is a parent containing a single SILVA/<version> subfolder
    try:
        resolved = Path(raw_str)
        if resolved.is_file():
            resolved = resolved.parent
        silva_sub = resolved / "SILVA"
        if silva_sub.is_dir():
            candidates = []
            for child in silva_sub.iterdir():
                if child.is_dir():
                    ver_t = _parse_version_tuple(child.name)
                    if ver_t is not None and (child / "SILVA.arb").exists():
                        candidates.append(ver_t)
            if candidates:
                return max(candidates)
    except OSError:
        pass

    return None


def is_silva_144_or_newer(
    dbs_loc: str = "",
    sample_taxonomies: Optional[Iterable[str]] = None,
) -> bool:
    """
    Return True if the SILVA database in use is version >= 144 (which includes a
    Kingdom rank ending in '-ati' between Domain and Phylum, and a source
    identifier at the 8th level). Falls back to inspecting sample SINA taxonomy
    strings when the version cannot be inferred from the path.
    """
    ver = detect_silva_version(dbs_loc)
    if ver is not None:
        return ver >= (144,)

    if sample_taxonomies is not None:
        for raw_tax in sample_taxonomies:
            if not raw_tax:
                continue
            ranks = [r.strip() for r in str(raw_tax).strip().rstrip(";").split(";") if r.strip()]
            if len(ranks) >= 2 and _KINGDOM_ATI_RE.match(ranks[1]):
                return True
            if len(ranks) >= 7:
                return True

    return False


def format_taxonomy_levels(levels: Sequence[str], num_levels: int) -> str:
    """
    Format a sequence of taxonomic ranks to exactly `num_levels` levels where
    every level ends with a semicolon ';' (so an N-level taxonomy contains N ';').
    """
    cleaned = [str(lvl).strip() for lvl in levels]
    padded = (cleaned + [""] * num_levels)[:num_levels]
    return ";".join(padded) + ";"


def parse_sina_taxonomy(
    raw_tax: str,
    is_v144_plus: bool = False,
) -> Tuple[str, str, str]:
    """
    Normalize a raw SINA `lca_tax_slv` taxonomy string into the classic 7-level
    format (`Domain;Phylum;Class;Order;Family;Genus;Species;`, with 7 semicolons
    and an empty 7th Species slot) so downstream tools (`ticlust`, etc.) run on
    a uniform 7-level hierarchy.

    For SILVA >= 144 (`is_v144_plus=True`):
    - Ignores the 8th rank (source identifier, e.g. 'metagenome').
    - Temporarily removes the Kingdom rank (2nd rank, e.g. 'Bacillati',
      'Pseudomonadati--other', 'Methanobacteriati') and returns it so it can be
      restored in post-processing.

    Returns `(classic_7level_tax, domain, kingdom)`.
    """
    cleaned = str(raw_tax or "").strip().rstrip(";").strip()
    if not cleaned:
        return format_taxonomy_levels([], 7), "", ""

    ranks = [r.strip() for r in cleaned.split(";")]
    if not is_v144_plus:
        classic_ranks = ranks[:6]
        domain = classic_ranks[0] if classic_ranks else ""
        return format_taxonomy_levels(classic_ranks, 7), domain, ""

    # SILVA >= 144: up to 8 ranks (Domain, Kingdom, Phylum, Class, Order, Family, Genus, Source)
    ranks = ranks[:7]
    domain = ranks[0] if len(ranks) >= 1 else ""
    if len(ranks) >= 2:
        kingdom = ranks[1]
        classic_ranks = [domain] + ranks[2:7]
    else:
        kingdom = ""
        classic_ranks = ranks[:1]

    return format_taxonomy_levels(classic_ranks, 7), domain, kingdom


def build_kingdom_lookup_from_sina_csv(
    sina_csv_path: str,
) -> Tuple[Dict[str, str], Dict[Tuple[str, str], str]]:
    """
    Parse a SINA CSV output file (`test_ZOTUs-Seqs.csv`, `test_OTUs-Seqs.csv`,
    or `aligned_<id>.csv`) and build:
    1. `seq_kingdom_map`: `{seq_id: kingdom}`
    2. `phylum_kingdom_map`: `{(domain, phylum): kingdom}`
    """
    seq_kingdom_map: Dict[str, str] = {}
    phylum_kingdom_map: Dict[Tuple[str, str], str] = {}
    csv_path = Path(PurePath(sina_csv_path))
    if not csv_path.is_file():
        return seq_kingdom_map, phylum_kingdom_map

    with open(csv_path, "r", encoding="utf-8") as csv_fio:
        reader = csv.reader(csv_fio)
        header_skipped = False
        for tokens in reader:
            if not tokens or not any(t.strip() for t in tokens):
                continue
            first_val = tokens[0].strip()
            if first_val.startswith("#"):
                continue
            if not header_skipped:
                header_skipped = True
                continue
            if len(tokens) <= 6:
                continue
            seq_id = first_val
            raw_tax = tokens[6].strip().rstrip(";").strip()
            if not raw_tax:
                continue
            ranks = [r.strip() for r in raw_tax.split(";")][:7]
            if len(ranks) >= 2 and ranks[1]:
                domain = ranks[0]
                kingdom = ranks[1]
                seq_kingdom_map[seq_id] = kingdom
                if len(ranks) >= 3 and ranks[2]:
                    phylum_kingdom_map[(domain, ranks[2])] = kingdom

    return seq_kingdom_map, phylum_kingdom_map


def build_cluster_id_kingdom_map(
    cluster_map_file: str,
    zotu_kingdom_map: Dict[str, str],
) -> Dict[str, str]:
    """
    Build a mapping from cluster ID (`OTU<N>` in TAC or `SOTU<N>` in TIC) to its
    Kingdom using the cluster-to-ZOTU map file (`Map-ZOTUs-OTUs-TAC.tab` or
    `Map-SOTU-ZOTU.tab`) and `zotu_kingdom_map`.
    """
    cluster_kingdom_map: Dict[str, str] = dict(zotu_kingdom_map)
    map_path = Path(PurePath(cluster_map_file))
    if not map_path.is_file():
        return cluster_kingdom_map

    with open(map_path, "r", encoding="utf-8") as map_fio:
        for raw_line in map_fio:
            line = raw_line.strip()
            if not line or line.lower().startswith("sotu\t") or line.startswith("#"):
                continue
            parts = [p.strip() for p in line.split("\t")]
            if len(parts) < 2:
                continue
            first, second = parts[0], parts[1]
            # In Map-ZOTUs-OTUs-TAC.tab: Zotu1 \t OTU1
            # In Map-SOTU-ZOTU.tab: SOTU1 \t Zotu1
            if first.lower().startswith("zotu"):
                zotu_id, cluster_id = first, second
            else:
                cluster_id, zotu_id = first, second
            if cluster_id not in cluster_kingdom_map and zotu_id in zotu_kingdom_map:
                cluster_kingdom_map[cluster_id] = zotu_kingdom_map[zotu_id]

    return cluster_kingdom_map


def ensure_trailing_semicolon_on_7level_tax(tax_str: str) -> str:
    """
    Ensure a 7-level taxonomy string has 7 levels and ends with a semicolon
    (7 semicolons total).
    """
    raw = str(tax_str or "").strip()
    if not raw:
        return format_taxonomy_levels([], 7)
    levels = [lvl.strip() for lvl in raw.split(";")]
    return format_taxonomy_levels(levels[:7], 7)


def ensure_7level_semicolons_in_table(table_path: str, sep: str = "\t") -> None:
    """
    Normalize the Taxonomy column of an existing 7-level table file in-place so
    each of the 7 levels ends with ';' (7 semicolons total).
    """
    path = Path(PurePath(table_path)).resolve()
    if not path.is_file():
        return
    with open(path, "r", encoding="utf-8") as fio:
        lines = [line.rstrip("\r\n") for line in fio]
    if not lines:
        return

    header = lines[0]
    header_cols = [c.strip().lower() for c in header.split(sep)]
    tax_indices = [i for i, col in enumerate(header_cols) if col == "taxonomy"]
    if not tax_indices:
        return
    tax_idx = tax_indices[-1]

    out_lines = [header]
    for line in lines[1:]:
        if not line.strip() or line.startswith("#"):
            out_lines.append(line)
            continue
        cols = line.split(sep)
        if tax_idx < len(cols):
            cols[tax_idx] = ensure_trailing_semicolon_on_7level_tax(cols[tax_idx])
        out_lines.append(sep.join(cols))

    with open(path, "w", encoding="utf-8") as fio:
        fio.write("\n".join(out_lines) + "\n")


def ensure_7level_semicolons_in_fasta(fasta_path: str) -> None:
    """
    Normalize `tax=...` headers in an existing 7-level FASTA file in-place so
    each of the 7 levels ends with ';' (7 semicolons total).
    """
    path = Path(PurePath(fasta_path)).resolve()
    if not path.is_file():
        return
    with open(path, "r", encoding="utf-8") as fio:
        lines = [line.rstrip("\r\n") for line in fio]

    out_lines = []
    for line in lines:
        if line.startswith(">"):
            match = _FASTA_TAX_RE.search(line)
            if match:
                norm_tax = ensure_trailing_semicolon_on_7level_tax(match.group("tax"))
                line = (
                    line[:match.start("tax")]
                    + norm_tax
                    + line[match.end("tax"):]
                )
        out_lines.append(line)

    with open(path, "w", encoding="utf-8") as fio:
        if out_lines:
            fio.write("\n".join(out_lines) + "\n")


def restore_kingdom_in_taxonomy(
    tax_str: str,
    seq_id: str = "",
    seq_kingdom_map: Optional[Dict[str, str]] = None,
    phylum_kingdom_map: Optional[Dict[Tuple[str, str], str]] = None,
    target_levels: int = 8,
) -> str:
    """
    Re-insert the SILVA >= 144 Kingdom rank between Domain (index 0) and Phylum
    (index 1) of a classic 7-level taxonomy string
    (`Domain;Phylum;Class;Order;Family;Genus;Species;`),
    formatting the output to `target_levels` (7 for ZOTU/OTU, 8 for TAC/TIC) with
    each level terminated by ';'.
    """
    seq_kingdom_map = seq_kingdom_map or {}
    phylum_kingdom_map = phylum_kingdom_map or {}

    raw = str(tax_str or "").strip()
    if not raw:
        return format_taxonomy_levels([], target_levels)

    classic_levels = ([lvl.strip() for lvl in raw.split(";")] + [""] * 7)[:7]
    uses_na_placeholders = any(lvl.startswith("NA-") for lvl in classic_levels)

    clean_seq_id = str(seq_id or "").strip().lstrip(">")
    if classic_levels[0] == "NA-Kingdom":
        domain = "NA-Domain"
        kingdom = seq_kingdom_map.get(clean_seq_id, "") or "NA-Kingdom"
    else:
        domain = classic_levels[0]
        phylum = classic_levels[1]
        kingdom = seq_kingdom_map.get(clean_seq_id, "")
        if not kingdom and domain and phylum and not phylum.startswith("NA"):
            kingdom = phylum_kingdom_map.get((domain, phylum), "")
        if not kingdom and uses_na_placeholders:
            if not domain:
                domain = "NA-Domain"
            kingdom = "NA-Kingdom"

    enriched_levels = [domain, kingdom] + classic_levels[1:7]
    return format_taxonomy_levels(enriched_levels[:target_levels], target_levels)


def write_kingdom_enriched_table(
    input_table: str,
    output_table: str,
    seq_kingdom_map: Optional[Dict[str, str]] = None,
    phylum_kingdom_map: Optional[Dict[Tuple[str, str], str]] = None,
    target_levels: int = 8,
    sep: str = "\t",
) -> None:
    """
    Write a copy of `input_table` to `output_table` with the Kingdom rank
    restored in the Taxonomy column (`target_levels` levels, each ending in ';').
    """
    in_path = Path(PurePath(input_table)).resolve()
    out_path = Path(PurePath(output_table)).resolve()
    with open(in_path, "r", encoding="utf-8") as in_fio:
        lines = [line.rstrip("\r\n") for line in in_fio]
    if not lines:
        return

    header = lines[0]
    header_cols = [c.strip().lower() for c in header.split(sep)]
    tax_indices = [i for i, col in enumerate(header_cols) if col == "taxonomy"]
    tax_idx = tax_indices[-1] if tax_indices else len(header_cols) - 1

    out_lines = [header]
    for line in lines[1:]:
        if not line.strip() or line.startswith("#"):
            out_lines.append(line)
            continue
        cols = line.split(sep)
        seq_id = cols[0].strip() if cols else ""
        if 0 <= tax_idx < len(cols):
            cols[tax_idx] = restore_kingdom_in_taxonomy(
                cols[tax_idx],
                seq_id=seq_id,
                seq_kingdom_map=seq_kingdom_map,
                phylum_kingdom_map=phylum_kingdom_map,
                target_levels=target_levels,
            )
        out_lines.append(sep.join(cols))

    with open(out_path, "w", encoding="utf-8") as out_fio:
        out_fio.write("\n".join(out_lines) + "\n")


def write_kingdom_enriched_fasta(
    input_fasta: str,
    output_fasta: str,
    seq_kingdom_map: Optional[Dict[str, str]] = None,
    phylum_kingdom_map: Optional[Dict[Tuple[str, str], str]] = None,
    target_levels: int = 8,
) -> None:
    """
    Write a copy of `input_fasta` to `output_fasta` with the Kingdom rank
    restored in `tax=...` headers (`target_levels` levels, each ending in ';').
    """
    in_path = Path(PurePath(input_fasta)).resolve()
    out_path = Path(PurePath(output_fasta)).resolve()
    with open(in_path, "r", encoding="utf-8") as in_fio:
        lines = [line.rstrip("\r\n") for line in in_fio]

    out_lines = []
    for line in lines:
        if line.startswith(">"):
            id_match = _FASTA_SEQ_ID_RE.match(line)
            seq_id = id_match.group("seq_id") if id_match else ""
            tax_match = _FASTA_TAX_RE.search(line)
            if tax_match:
                enriched_tax = restore_kingdom_in_taxonomy(
                    tax_match.group("tax"),
                    seq_id=seq_id,
                    seq_kingdom_map=seq_kingdom_map,
                    phylum_kingdom_map=phylum_kingdom_map,
                    target_levels=target_levels,
                )
                line = (
                    line[:tax_match.start("tax")]
                    + enriched_tax
                    + line[tax_match.end("tax"):]
                )
        out_lines.append(line)

    with open(out_path, "w", encoding="utf-8") as out_fio:
        if out_lines:
            out_fio.write("\n".join(out_lines) + "\n")


def _with_level_suffix(file_path: Path, target_levels: int) -> Path:
    suffix_tag = f"-{target_levels}Levels"
    if file_path.stem.endswith(suffix_tag):
        return file_path
    return file_path.with_name(f"{file_path.stem}{suffix_tag}{file_path.suffix}")


def generate_kingdom_enriched_pipeline_outputs(
    pipeline_dir: str,
    target_levels: int = 8,
    seq_kingdom_map: Optional[Dict[str, str]] = None,
    phylum_kingdom_map: Optional[Dict[Tuple[str, str], str]] = None,
    krona_tool: str = "",
    create_krona_fn: Optional[Callable[[str, str, str], str]] = None,
) -> List[str]:
    """
    For a finished pipeline directory (`ZOTU_Pipeline`, `OTU_Pipeline`,
    `TAC_Pipeline`, or `TIC_Pipeline`), generate the Kingdom-enriched version
    (`*-8Levels.*`) of all output tables, FASTA files, and Krona HTML charts
    that contain taxonomy.
    """
    dir_path = Path(PurePath(pipeline_dir)).resolve()
    if not dir_path.is_dir():
        return []

    created_files: List[str] = []
    primary_table_for_krona: Optional[Path] = None
    is_tac_or_tic_dir = dir_path.name.lower() in {"tac_pipeline", "tic_pipeline"}
    staged_zotu_inputs = {"zotus-seqs.fasta", "zotus-table.tab"}

    for item in sorted(dir_path.iterdir()):
        if not item.is_file() or _EXPLICIT_LEVEL_SUFFIX_RE.search(item.name):
            continue
        if is_tac_or_tic_dir and item.name.lower() in staged_zotu_inputs:
            continue

        # Tables with a Taxonomy column (excluding Map-* files)
        if item.suffix.lower() == ".tab" and not item.name.startswith("Map-"):
            with open(item, "r", encoding="utf-8") as fio:
                first_line = fio.readline().strip().lower()
            if "taxonomy" in [c.strip() for c in first_line.split("\t")]:
                out_tab = _with_level_suffix(item, target_levels)
                write_kingdom_enriched_table(
                    str(item),
                    str(out_tab),
                    seq_kingdom_map=seq_kingdom_map,
                    phylum_kingdom_map=phylum_kingdom_map,
                    target_levels=target_levels,
                )
                created_files.append(str(out_tab))
                if re.match(r"^[ZS]?OTUs-Table(?:-TAC|-TIC)?\.tab$", item.name, re.IGNORECASE):
                    primary_table_for_krona = out_tab

        # FASTA files containing tax= headers (or TIC invalid-taxonomy output)
        elif (
            item.suffix.lower() == ".fasta"
            and not item.name.startswith(("test_", "shortened_algn_", "filtered_by_abundance_"))
        ):
            is_tic_invalid_fasta = item.name.lower() == "invalid-tax-seqs-tic.fasta"
            has_tax = False
            with open(item, "r", encoding="utf-8") as fio:
                for _ in range(20):
                    line = fio.readline()
                    if not line:
                        break
                    if line.startswith(">") and "tax=" in line.lower():
                        has_tax = True
                        break
            if has_tax or is_tic_invalid_fasta:
                out_fasta = _with_level_suffix(item, target_levels)
                write_kingdom_enriched_fasta(
                    str(item),
                    str(out_fasta),
                    seq_kingdom_map=seq_kingdom_map,
                    phylum_kingdom_map=phylum_kingdom_map,
                    target_levels=target_levels,
                )
                created_files.append(str(out_fasta))

    has_krona_prereqs = (
        bool(create_krona_fn and krona_tool)
        and primary_table_for_krona is not None
        and primary_table_for_krona.is_file()
    )
    if has_krona_prereqs:
        for html_item in sorted(dir_path.glob("*Krona.html")):
            if _EXPLICIT_LEVEL_SUFFIX_RE.search(html_item.name):
                continue
            out_html = _with_level_suffix(html_item, target_levels)
            create_krona_fn(krona_tool, str(primary_table_for_krona), out_html.name)
            created_files.append(str(out_html))

    return created_files


def write_taxonomy_levels_readme(output_dir: str) -> str:
    """
    Write README_Taxonomy_Levels.txt into `output_dir` explaining the difference
    between standard 7-level outputs and the SILVA >= 144 *-8Levels files.
    """
    dest_dir = Path(PurePath(output_dir)).resolve()
    dest_dir.mkdir(parents=True, exist_ok=True)
    readme_path = dest_dir / "README_Taxonomy_Levels.txt"
    content = (
        "===============================================================================\n"
        "TACTIC Pipeline — Taxonomy Output Formats & SILVA Version Compatibility\n"
        "===============================================================================\n\n"
        "You are using a SILVA reference database version >= 144 (or a database featuring\n"
        "the updated bacterial kingdom nomenclature). In SILVA release 144 onward, a\n"
        "formal Kingdom rank ending in '-ati' (e.g. Bacillati, Pseudomonadati, Fusobacteriati)\n"
        "was introduced between Domain and Phylum.\n\n"
        "To guarantee seamless compatibility with both traditional bioinformatics tools\n"
        "and updated taxonomy standards, TACTIC generates parallel sets of output files:\n\n"
        "1. Standard Outputs (e.g., SOTUs-Table-TIC.tab, OTUs-Table-TAC.tab,\n"
        "   ZOTUs-Table-OnlyBacteria-FullTaxonomy.tab, OTUs-Table.tab):\n"
        "   - Format: Domain;Phylum;Class;Order;Family;Genus;Species; (7 levels, 7 ';')\n"
        "   - Note: In pipelines without species-level taxonomy imputation (ZOTU, OTU, TAC),\n"
        "     unassigned lower ranks such as Species remain empty between semicolons.\n"
        "   - Purpose: Direct compatibility with legacy tools, QIIME, PICRUSt2, and R\n"
        "     packages that strictly expect traditional 7-level taxonomies without Kingdom.\n\n"
        "2. Enriched '*-8Levels.*' Output Files (generated across TIC, TAC, ZOTU, OTU, and -iz):\n"
        "   - Format: Domain;Kingdom;Phylum;Class;Order;Family;Genus;Species; (8 levels, 8 ';')\n"
        "   - Includes the bacterial Kingdom rank between Domain and Phylum down to Species\n"
        "     (with unassigned ranks left empty between semicolons).\n\n"
        "Which file should I use?\n"
        "------------------------\n"
        "• Use Standard files if your downstream scripts (e.g. PICRUSt2, legacy phyloseq\n"
        "  scripts) require classic 7-level taxonomies.\n"
        "• Use *-8Levels files if you want the most up-to-date SILVA 144+\n"
        "  taxonomic nomenclature including bacterial kingdoms (e.g. for Rhea, Namco, or Krona).\n"
        "===============================================================================\n"
    )
    with open(readme_path, "w", encoding="utf-8") as fio:
        fio.write(content)
    return str(readme_path)


def get_pipeline_readme_content(pipeline_type: str) -> str:
    """
    Return simple-language text documentation for the requested pipeline output
    directory ('otu', 'zotu', 'tac', 'tic', 'picrust2', 'shared', or 'root'/'analysis').
    """
    norm = str(pipeline_type).strip().lower().replace("-", "_")
    if norm in {"otu", "otu_pipeline", "de_novo"}:
        return (
            "===============================================================================\n"
            "TACTIC Pipeline — OTU Pipeline Output Files\n"
            "===============================================================================\n\n"
            "What is the OTU Pipeline?\n"
            "-------------------------\n"
            "The OTU (Operational Taxonomic Unit) pipeline clusters closely related 16S rRNA\n"
            "gene sequences together at 97% sequence similarity. Sequences that are 97%\n"
            "identical or higher are grouped into a single representative unit (an OTU).\n"
            "This is the classic, widely established approach for studying microbial\n"
            "communities.\n\n"
            "Below is an explanation of each output file in this folder and how to use it.\n\n\n"
            "1. Primary Abundance Tables\n"
            "---------------------------\n"
            "• OTUs-Table.tab\n"
            "  - What it is: A tab-separated text table containing the count (abundance) of\n"
            "    each 97% OTU across all your samples. Rows represent OTU clusters, columns\n"
            "    represent your sample IDs, and the final column shows the assigned\n"
            "    taxonomy (Domain down to Species).\n"
            "  - What it is used for: This is the main data table for your downstream\n"
            "    analysis. You can open it in Microsoft Excel, Google Sheets, R (using\n"
            "    packages like phyloseq, vegan, or DESeq2), Python, or statistical software\n"
            "    to compare microbial composition, calculate diversity metrics, and identify\n"
            "    differentially abundant bacteria between sample groups.\n\n"
            "• SampleMinCount_Normalized-OTUs-Table.tab\n"
            "  - What it is: The OTU abundance table normalized by subsampling (rarefying)\n"
            "    every sample to match the sample with the lowest total read count.\n"
            "  - What it is used for: Fair comparison between samples. Because sequencing\n"
            "    depth often varies between samples, this normalized table helps remove\n"
            "    library-size bias before computing alpha and beta diversity.\n\n"
            "• Spike_Normalized-OTUs-Table.tab (present only if spike-ins were used)\n"
            "  - What it is: The OTU abundance table normalized against the internal spike-in\n"
            "    DNA standards added to your samples before sequencing.\n"
            "  - What it is used for: Absolute quantification. While standard tables provide\n"
            "    relative percentages, this table allows you to estimate absolute microbial\n"
            "    loads (e.g., cell counts or copy numbers per gram of sample).\n\n\n"
            "2. Representative Sequences\n"
            "---------------------------\n"
            "• OTUs-Seqs.fasta\n"
            "  - What it is: A standard FASTA file containing the representative (centroid)\n"
            "    DNA sequence for each 97% OTU cluster. The header line for each sequence\n"
            "    contains the OTU identifier and its assigned taxonomy.\n"
            "  - What it is used for: Inspecting the actual DNA sequences of your OTUs. Use\n"
            "    this file for sequence alignments, BLAST searches against NCBI, designing\n"
            "    targeted PCR primers, or constructing custom phylogenetic trees.\n\n\n"
            "3. Phylogenetic Trees\n"
            "---------------------\n"
            "• OTUs-Tree-nj.tre\n"
            "  - What it is: A phylogenetic tree of the OTU representative sequences saved\n"
            "    in standard Newick format (.tre), constructed using the Neighbor-Joining\n"
            "    algorithm.\n"
            "  - What it is used for: Phylogeny-aware microbiome analyses, such as calculating\n"
            "    UniFrac distances (weighted or unweighted) in QIIME 2 or phyloseq. You can\n"
            "    also visualize this tree in interactive viewers like FigTree, iTOL, or\n"
            "    Microreact.\n\n\n"
            "4. Interactive Visualizations\n"
            "-----------------------------\n"
            "• OTUs-Krona.html\n"
            "  - What it is: An interactive, visual taxonomy chart (Krona plot) that opens in\n"
            "    any web browser.\n"
            "  - What it is used for: Exploring the bacterial composition of your samples\n"
            "    visually without any coding. Simply double-click the file to open it in\n"
            "    Chrome, Firefox, Safari, or Edge. You can click on the colored concentric\n"
            "    rings to zoom in and out through Domain, Phylum, Class, Order, Family,\n"
            "    Genus, and Species.\n\n\n"
            "5. SILVA 144+ Dual Outputs (*-8Levels.*)\n"
            "----------------------------------------\n"
            "(Generated when a reference database with bacterial kingdoms is used,\n"
            "e.g. SILVA >= 144)\n\n"
            "• OTUs-Table-8Levels.tab\n"
            "• OTUs-Seqs-8Levels.fasta\n"
            "• OTUs-Krona-8Levels.html\n"
            "• SampleMinCount_Normalized-OTUs-Table-8Levels.tab\n"
            "• Spike_Normalized-OTUs-Table-8Levels.tab (if spike-ins were used)\n"
            "  - What they are: The exact same OTU results as above, but with 8 taxonomic\n"
            "    ranks that explicitly include the bacterial Kingdom rank (e.g., Bacillati,\n"
            "    Pseudomonadati, Fusobacteriati) between Domain and Phylum.\n"
            "  - What they are used for: Use the *-8Levels files if you need the latest SILVA\n"
            "    144+ kingdom classification (e.g., for Namco, Rhea, or Krona). Use the\n"
            "    standard (7-level) files above if your downstream software (like PICRUSt2\n"
            "    or older R scripts) expects classic 7-level taxonomies.\n"
            "===============================================================================\n"
        )
    elif norm in {"zotu", "zotu_pipeline", "asv"}:
        return (
            "===============================================================================\n"
            "TACTIC Pipeline — ZOTU Pipeline Output Files\n"
            "===============================================================================\n\n"
            "What is the ZOTU Pipeline?\n"
            "--------------------------\n"
            "The ZOTU (Zero-radius Operational Taxonomic Unit) pipeline performs denoising\n"
            "on your amplicon sequencing data to identify true, exact biological sequences\n"
            "(also widely known as Amplicon Sequence Variants, or ASVs).\n\n"
            "Unlike traditional 97% OTUs, ZOTUs are resolved down to 100% sequence identity.\n"
            "Sequences that differ by even a single nucleotide are kept separate as distinct\n"
            "biological variants, eliminating artificial clustering and preserving maximum\n"
            "biological resolution.\n\n"
            "Below is an explanation of each output file in this folder and how to use it.\n\n\n"
            "1. Primary Abundance Tables\n"
            "---------------------------\n"
            "• ZOTUs-Table.tab\n"
            "  - What it is: A tab-separated text table showing the exact number of times\n"
            "    each ZOTU sequence was observed in each sample. Rows represent unique\n"
            "    ZOTUs (ASVs), columns represent your samples, and the final column provides\n"
            "    the assigned taxonomic classification (Domain down to Species).\n"
            "  - What it is used for: This is your high-resolution abundance table. Open it\n"
            "    in Microsoft Excel, R, Python, or statistical software to analyze community\n"
            "    structure, detect strain-level differences between groups, and perform\n"
            "    differential abundance testing.\n\n"
            "• SampleMinCount_Normalized-ZOTUs-Table.tab\n"
            "  - What it is: A depth-normalized version of the ZOTU abundance table where\n"
            "    all samples are subsampled (rarefied) to match the sequencing depth of the\n"
            "    sample with the fewest reads.\n"
            "  - What it is used for: Fair, unbiased statistical comparisons of diversity\n"
            "    (such as Shannon diversity or Bray-Curtis dissimilarity) without being\n"
            "    distorted by differences in sequencing depth between samples.\n\n"
            "• Spike_Normalized-ZOTUs-Table.tab (present only if spike-ins were used)\n"
            "  - What it is: An abundance table normalized against the known concentration of\n"
            "    internal DNA spike-in standards added to your samples.\n"
            "  - What it is used for: Absolute quantification. This allows you to estimate\n"
            "    absolute microbial abundance (e.g., bacterial cells or copies per unit of\n"
            "    sample) rather than relying only on relative percentages.\n\n\n"
            "2. Representative Sequences\n"
            "---------------------------\n"
            "• ZOTUs-Seqs.fasta\n"
            "  - What it is: A standard FASTA file containing the exact, full-length DNA\n"
            "    sequence for every identified ZOTU. Each sequence header includes the ZOTU\n"
            "    identifier and its taxonomic assignment.\n"
            "  - What it is used for: Sequence-level investigations. You can use these\n"
            "    sequences to query public databases (NCBI BLAST), align against specific\n"
            "    reference genomes, construct phylogenetic trees, or design targeted primers.\n\n\n"
            "3. Phylogenetic Trees\n"
            "---------------------\n"
            "• ZOTUs-Tree-nj.tre\n"
            "  - What it is: A phylogenetic tree of all ZOTU sequences in standard Newick\n"
            "    format (.tre), built using the Neighbor-Joining method.\n"
            "  - What it is used for: Phylogeny-aware analyses, such as calculating weighted\n"
            "    and unweighted UniFrac distances in QIIME 2, phyloseq, or R packages, and\n"
            "    visualizing evolutionary relationships in tree viewers like FigTree or iTOL.\n\n\n"
            "4. Interactive Visualizations\n"
            "-----------------------------\n"
            "• ZOTUs-Krona.html\n"
            "  - What it is: An interactive, visual taxonomy chart (Krona plot) that opens in\n"
            "    any web browser.\n"
            "  - What it is used for: Visually browsing the taxonomic composition of your\n"
            "    ZOTUs. Double-click the file to open it in Chrome, Firefox, Safari, or Edge.\n"
            "    Click through the concentric rings to easily explore proportions from\n"
            "    Domain down to Species.\n\n\n"
            "5. SILVA 144+ Dual Outputs (*-8Levels.*)\n"
            "----------------------------------------\n"
            "(Generated when a reference database with bacterial kingdoms is used,\n"
            "e.g. SILVA >= 144)\n\n"
            "• ZOTUs-Table-8Levels.tab\n"
            "• ZOTUs-Seqs-8Levels.fasta\n"
            "• ZOTUs-Krona-8Levels.html\n"
            "• SampleMinCount_Normalized-ZOTUs-Table-8Levels.tab\n"
            "• Spike_Normalized-ZOTUs-Table-8Levels.tab (if spike-ins were used)\n"
            "  - What they are: The exact same ZOTU results as above, but with 8 taxonomic\n"
            "    ranks that include the bacterial Kingdom rank (e.g., Bacillati,\n"
            "    Pseudomonadati, Fusobacteriati) between Domain and Phylum.\n"
            "  - What they are used for: Use the *-8Levels files if you need the latest SILVA\n"
            "    144+ kingdom classification. Use the standard (7-level) files above if\n"
            "    downstream software (like PICRUSt2 or older R scripts) strictly requires\n"
            "    classic 7-level taxonomies.\n"
            "===============================================================================\n"
        )
    elif norm in {"tac", "tac_pipeline"}:
        return (
            "===============================================================================\n"
            "TACTIC Pipeline — TAC Pipeline Output Files\n"
            "===============================================================================\n\n"
            "What is the TAC Pipeline?\n"
            "-------------------------\n"
            "The TAC (Taxonomy-Agnostic Clustering) pipeline takes the clean, denoised\n"
            "ZOTUs (Amplicon Sequence Variants) and clusters them at 98.7% sequence\n"
            "similarity — the standard threshold for bacterial species — without using any\n"
            "prior taxonomy labels.\n\n"
            "This approach gives you the best of both worlds: it eliminates sequencing noise\n"
            "and PCR errors via denoising first, and then groups closely related variants\n"
            "into biologically sensible, species-level operational units.\n\n"
            "Below is an explanation of each output file in this folder and how to use it.\n\n\n"
            "1. Primary Abundance Tables\n"
            "---------------------------\n"
            "• OTUs-Table-TAC.tab\n"
            "  - What it is: A tab-separated text table containing the abundance count of\n"
            "    each TAC OTU cluster in each sample. Rows are TAC OTUs (clustered at 98.7%\n"
            "    similarity), columns are your sample IDs, and the last column shows the\n"
            "    assigned taxonomy (Domain down to Species).\n"
            "  - What it is used for: This is the primary abundance table for TAC mode.\n"
            "    Open it in Microsoft Excel, R, Python, or statistical software to analyze\n"
            "    microbiome diversity and compare microbial communities between sample groups.\n\n"
            "• SampleMinCount_Normalized-OTUs-Table-TAC.tab\n"
            "  - What it is: A depth-normalized version of the TAC abundance table where\n"
            "    all samples are subsampled (rarefied) to match the sequencing depth of the\n"
            "    sample with the fewest reads.\n"
            "  - What it is used for: Unbiased statistical comparisons of alpha and beta\n"
            "    diversity across samples with different library sizes.\n\n"
            "• Spike_Normalized-OTUs-Table-TAC.tab (present only if spike-ins were used)\n"
            "  - What it is: An abundance table normalized against internal DNA spike-in\n"
            "    standards added to your samples before sequencing.\n"
            "  - What it is used for: Estimating absolute microbial abundance (cell or copy\n"
            "    numbers per unit of sample) instead of relative percentages.\n\n\n"
            "2. Sequence Files and Traceability Mappings\n"
            "-------------------------------------------\n"
            "• OTUs-Seqs-TAC.fasta\n"
            "  - What it is: A FASTA file containing the representative (centroid) DNA\n"
            "    sequence for each TAC OTU cluster, with taxonomy in the headers.\n"
            "  - What it is used for: Sequence searches (BLAST), alignments, and phylogenetic\n"
            "    analyses of the clustered representative sequences.\n\n"
            "• Map-ZOTUs-OTUs-TAC.tab\n"
            "  - What it is: A two-column mapping table that connects each original ZOTU\n"
            "    identifier to the TAC OTU cluster it was assigned to.\n"
            "  - What it is used for: Full traceability. Use this file if you need to know\n"
            "    which exact sequence variants (ZOTUs) were merged into a given TAC OTU.\n\n\n"
            "3. Phylogenetic Trees\n"
            "---------------------\n"
            "• OTUs-Tree-nj-TAC.tre\n"
            "  - What it is: A phylogenetic tree of the TAC OTU centroid sequences in\n"
            "    standard Newick format (.tre), built using Neighbor-Joining.\n"
            "  - What it is used for: Tree-based diversity analyses (such as UniFrac) in\n"
            "    QIIME 2 or phyloseq, and visualization in tree viewers like FigTree or iTOL.\n\n\n"
            "4. Interactive Visualizations\n"
            "-----------------------------\n"
            "• OTUs-Krona.html\n"
            "  - What it is: An interactive Krona plot chart of taxonomic proportions.\n"
            "  - What it is used for: Double-click to open in any web browser and explore\n"
            "    the taxonomic breakdown of your TAC OTUs visually.\n\n\n"
            "5. SILVA 144+ Dual Outputs (*-8Levels.*)\n"
            "----------------------------------------\n"
            "(Generated when a reference database with bacterial kingdoms is used,\n"
            "e.g. SILVA >= 144)\n\n"
            "• OTUs-Table-TAC-8Levels.tab\n"
            "• OTUs-Seqs-TAC-8Levels.fasta\n"
            "• OTUs-Krona-8Levels.html\n"
            "• SampleMinCount_Normalized-OTUs-Table-TAC-8Levels.tab\n"
            "• Spike_Normalized-OTUs-Table-TAC-8Levels.tab (if spike-ins were used)\n"
            "  - What they are: The exact same TAC results as above, but with 8 taxonomic\n"
            "    ranks including the bacterial Kingdom rank (e.g., Bacillati, Pseudomonadati).\n"
            "  - What they are used for: Use the *-8Levels files for the latest SILVA 144+\n"
            "    taxonomy, or use the standard files for downstream tools expecting 7 levels.\n"
            "===============================================================================\n"
        )
    elif norm in {"tic", "tic_pipeline"}:
        return (
            "===============================================================================\n"
            "TACTIC Pipeline — TIC Pipeline Output Files\n"
            "===============================================================================\n\n"
            "What is the TIC Pipeline?\n"
            "-------------------------\n"
            "The TIC (Taxonomy-Informed Clustering) pipeline groups amplicon sequence\n"
            "variants (ZOTUs) into biologically meaningful taxonomic units by using BOTH\n"
            "sequence similarity thresholds AND taxonomic classifications simultaneously.\n\n"
            "Sequences sharing common taxonomy are clustered hierarchically at three\n"
            "biologically validated similarity thresholds:\n"
            "  • Family level  : 90% sequence similarity (default unless changed)  -> FOTUs\n"
            "  • Genus level   : 95% sequence similarity (default unless changed)  -> GOTUs\n"
            "  • Species level : 98.7% sequence similarity (default unless changed) -> SOTUs (Species-level OTUs)\n\n"
            "This prevents unrelated organisms that happen to share sequence similarity\n"
            "from being incorrectly merged, while accurately grouping true species variants.\n\n"
            "Below is an explanation of each output file in this folder and how to use it.\n\n\n"
            "1. Primary Abundance Tables\n"
            "---------------------------\n"
            "• SOTUs-Table-TIC.tab\n"
            "  - What it is: The primary abundance table for TIC mode. Each row represents a\n"
            "    Species-level OTU (SOTU), columns represent your samples, and the final\n"
            "    column provides the complete imputed taxonomy down to species level.\n"
            "  - What it is used for: This is the main data table for TIC analysis. Open it\n"
            "    in Microsoft Excel, R, Python, or statistical packages for community\n"
            "    composition analysis, diversity calculations, and biomarker discovery.\n\n"
            "• SampleMinCount_Normalized-SOTUs-Table-TIC.tab\n"
            "  - What it is: The SOTU abundance table normalized by subsampling (rarefying)\n"
            "    every sample to the lowest sample sequencing depth.\n"
            "  - What it is used for: Unbiased alpha and beta diversity comparisons without\n"
            "    sequencing depth artifacts.\n\n"
            "• Spike_Normalized-SOTUs-Table-TIC.tab (present only if spike-ins were used)\n"
            "  - What it is: The SOTU abundance table normalized using internal DNA spike-in\n"
            "    standards.\n"
            "  - What it is used for: Absolute microbial abundance quantification (e.g.,\n"
            "    cell or copy numbers per gram/ml of sample).\n\n"
            "• ZOTU-Table-OnlyBacteria-FullTaxonomy.tab\n"
            "  - What it is: A high-resolution ZOTU table containing only validated bacterial\n"
            "    sequences with full taxonomic annotations.\n"
            "  - What it is used for: Recommended for standard bacterial diversity and functional\n"
            "    analyses (e.g., PICRUSt2). It excludes non-target amplicon noise, host\n"
            "    mitochondrial/chloroplast contamination, and archaea.\n\n"
            "• ZOTU-Table-All-FullTaxonomy.tab (present if non-bacterial sequences were detected)\n"
            "  - What it is: A full ZOTU table containing all sequence variants, including\n"
            "    bacteria, archaea, and unassigned/non-bacterial sequences.\n"
            "  - What it is used for: Quality control and broad-domain studies. Use this table\n"
            "    if you want to quantify non-bacterial organisms (such as Archaea) or check what\n"
            "    fraction of sequencing reads were non-bacterial or off-target.\n"
            "    (Note: If no non-bacterial sequences are detected in your dataset, this table\n"
            "    is identical to ZOTU-Table-OnlyBacteria-FullTaxonomy.tab and is automatically\n"
            "    omitted to avoid redundant duplicate files).\n\n\n"
            "2. Representative Sequences & Quality Control\n"
            "---------------------------------------------\n"
            "• SOTUs-Seqs-TIC.fasta\n"
            "  - What it is: A standard FASTA file containing the representative (centroid)\n"
            "    DNA sequence for each Species-level OTU (SOTU), with taxonomy in the headers.\n"
            "  - What it is used for: Sequence alignments, database queries (BLAST), and\n"
            "    verifying representative sequences.\n\n"
            "• Invalid-Tax-Seqs-TIC.fasta\n"
            "  - What it is: A FASTA file containing sequences that failed taxonomic\n"
            "    validation (e.g., non-bacterial reads, mitochondrial/chloroplast DNA, or\n"
            "    sequences that could not be assigned valid 16S taxonomy).\n"
            "  - What it is used for: Quality control. Inspect this file to see which reads\n"
            "    were filtered out and ensure your primers did not amplify off-target DNA.\n\n\n"
            "3. Hierarchical Mapping Tables\n"
            "------------------------------\n"
            "• Map-SOTU-ZOTU.tab\n"
            "  - What it is: A mapping table showing which exact ZOTUs were clustered into\n"
            "    each SOTU (Species-level OTU).\n"
            "  - What it is used for: Investigating within-species sequence diversity and\n"
            "    tracking which sequence variants belong to which species cluster.\n\n"
            "• Map-GOTU-SOTU.tab\n"
            "  - What it is: A mapping table showing how SOTUs cluster into Genus-level\n"
            "    OTUs (GOTUs, 95% similarity).\n"
            "  - What it is used for: Exploring genus-level groupings and relationships.\n\n"
            "• Map-FOTU-GOTU.tab\n"
            "  - What it is: A mapping table showing how Genus-level OTUs cluster into\n"
            "    Family-level OTUs (FOTUs, 90% similarity).\n"
            "  - What it is used for: Exploring family-level groupings and relationships.\n\n\n"
            "4. Phylogenetic Trees\n"
            "---------------------\n"
            "• SOTUs-Tree-nj-TIC.tre\n"
            "  - What it is: A phylogenetic tree of SOTU centroids in Newick format (.tre),\n"
            "    built using Neighbor-Joining.\n"
            "  - What it is used for: Phylogeny-aware diversity calculations (UniFrac) in\n"
            "    QIIME 2 or phyloseq, and visualization in tree viewers like FigTree or iTOL.\n\n\n"
            "5. Interactive Visualizations\n"
            "-----------------------------\n"
            "• SOTUs-Krona.html\n"
            "  - What it is: An interactive Krona taxonomy chart that opens in your browser.\n"
            "  - What it is used for: Double-click to view the taxonomic hierarchy of your\n"
            "    SOTUs in Chrome, Firefox, Safari, or Edge.\n\n\n"
            "6. SILVA 144+ Dual Outputs (*-8Levels.*)\n"
            "----------------------------------------\n"
            "(Generated when a reference database with bacterial kingdoms is used,\n"
            "e.g. SILVA >= 144)\n\n"
            "• SOTUs-Table-TIC-8Levels.tab\n"
            "• SOTUs-Seqs-TIC-8Levels.fasta\n"
            "• SOTUs-Krona-8Levels.html\n"
            "• Invalid-Tax-Seqs-TIC-8Levels.fasta\n"
            "• SampleMinCount_Normalized-SOTUs-Table-TIC-8Levels.tab\n"
            "• Spike_Normalized-SOTUs-Table-TIC-8Levels.tab (if spike-ins were used)\n"
            "  - What they are: The exact same TIC results as above, but with 8 taxonomic\n"
            "    ranks including the bacterial Kingdom rank (e.g., Bacillati, Pseudomonadati).\n"
            "  - What they are used for: Use the *-8Levels files for the latest SILVA 144+\n"
            "    kingdom taxonomy, or use the standard files for tools requiring 7 levels.\n"
            "===============================================================================\n"
        )
    elif norm in {"picrust2", "picrust2_inputs", "picrust"}:
        return (
            "===============================================================================\n"
            "TACTIC Pipeline — PICRUSt2 Input Files\n"
            "===============================================================================\n\n"
            "What is PICRUSt2?\n"
            "-----------------\n"
            "PICRUSt2 (Phylogenetic Investigation of Communities by Reconstruction of\n"
            "Unobserved States) is a bioinformatics tool used to predict the functional\n"
            "metabolic potential (enzymes, gene families, and metabolic pathways) of a\n"
            "microbial community from 16S rRNA sequencing data.\n\n"
            "To run PICRUSt2, you need two clean input files:\n"
            "  1. An abundance table of sequence counts across samples.\n"
            "  2. A matching FASTA file of the representative DNA sequences.\n\n"
            "The files in this folder have been specially pre-formatted to be 100%\n"
            "compatible with the official PICRUSt2 pipeline (picrust2_pipeline.py).\n\n\n"
            "Files in this Folder\n"
            "--------------------\n"
            "• <Pipeline>__*OTUs-Table*-For-PICRUSt2.tab\n"
            "  - What it is: Clean abundance tables in tab-separated format where the\n"
            "    taxonomy column has been removed and sequence identifiers match the FASTA.\n"
            "  - What it is used for: Provide this file to the `--input` (-i) argument of\n"
            "    PICRUSt2:\n"
            "      picrust2_pipeline.py -s <Seqs.fasta> -i <Table.tab> -o picrust2_out ...\n\n"
            "• <Pipeline>__*OTUs-Seqs*-For-PICRUSt2.fasta\n"
            "  - What it is: Clean FASTA files containing only the sequence identifier and\n"
            "    the DNA sequence (all spaces, descriptions, and taxonomy headers removed).\n"
            "  - What it is used for: Provide this file to the `--study_fasta` (-s) argument\n"
            "    of PICRUSt2 for sequence alignment and phylogenetic placement.\n"
            "===============================================================================\n"
        )
    elif norm in {"shared", "shared_pipeline"}:
        return (
            "===============================================================================\n"
            "TACTIC Pipeline — Shared Intermediate Sequence Files\n"
            "===============================================================================\n\n"
            "What is the Shared Folder?\n"
            "--------------------------\n"
            "This folder contains pooled and prepared sequence files that serve as the common\n"
            "starting input for the downstream clustering pipelines (OTU, ZOTU, TAC, TIC).\n\n\n"
            "Files in this Folder\n"
            "--------------------\n"
            "• sorted.fasta\n"
            "  - What it is: The combined unique DNA sequences from all samples, sorted by\n"
            "    their abundance (sequence count) across the entire dataset.\n"
            "  - What it is used for: Input to clustering algorithms (e.g., USEARCH cluster_otus\n"
            "    and unoise3), which require sequences sorted from highest to lowest abundance.\n\n"
            "• derep.fasta\n"
            "  - What it is: Dereplicated sequences containing size annotations indicating\n"
            "    how many times each sequence was observed.\n"
            "  - What it is used for: High-speed sequence dereplication and table building.\n\n"
            "• relabel.tab\n"
            "  - What it is: A mapping table connecting internal sequence labels to original\n"
            "    sample identifiers and read IDs.\n"
            "  - What it is used for: Traceability back to the original FASTQ read headers.\n"
            "===============================================================================\n"
        )
    else:  # root / analysis
        return (
            "===============================================================================\n"
            "TACTIC Pipeline — Analysis Results Directory Overview\n"
            "===============================================================================\n\n"
            "Welcome to your TACTIC analysis results! This directory contains all output\n"
            "files, abundance tables, representative sequences, phylogenetic trees, and\n"
            "interactive visualizations generated from your 16S rRNA sequencing data.\n\n"
            "Below is an overview of the folder organization and what each file contains.\n\n\n"
            "1. Pipeline Output Subdirectories\n"
            "---------------------------------\n"
            "Depending on the analysis mode(s) you selected (--analysis-mode), one or more\n"
            "of the following pipeline folders are present:\n\n"
            "• OTU_Pipeline/\n"
            "  - Results from 97% sequence similarity clustering (classic OTUs).\n"
            "  - See OTU_Pipeline/README.txt for descriptions of every file in that folder.\n\n"
            "• ZOTU_Pipeline/\n"
            "  - Results from exact sequence variant denoising (ZOTUs / ASVs).\n"
            "  - See ZOTU_Pipeline/README.txt for descriptions of every file in that folder.\n\n"
            "• TAC_Pipeline/\n"
            "  - Results from Taxonomy-Agnostic Clustering at 98.7% similarity.\n"
            "  - See TAC_Pipeline/README.txt for descriptions of every file in that folder.\n\n"
            "• TIC_Pipeline/\n"
            "  - Results from Taxonomy-Informed Clustering into Species-level OTUs (SOTUs).\n"
            "  - See TIC_Pipeline/README.txt for descriptions of every file in that folder.\n\n"
            "• PICRUSt2-Inputs/\n"
            "  - Pre-formatted abundance tables and FASTA files ready for functional\n"
            "    metabolic pathway prediction with PICRUSt2.\n"
            "  - See PICRUSt2-Inputs/README.txt for instructions on running PICRUSt2.\n\n"
            "• Shared/\n"
            "  - Pooled, dereplicated, and abundance-sorted sequences used as input for\n"
            "    the clustering algorithms.\n"
            "  - See Shared/README.txt for details.\n\n\n"
            "2. Root Directory Files\n"
            "-----------------------\n"
            "• README.txt\n"
            "  - This overview file explaining the results folder structure and files.\n\n"
            "• README_Taxonomy_Levels.txt (present when SILVA >= 144 is used)\n"
            "  - Detailed guide explaining the difference between standard (7-level) outputs\n"
            "    and SILVA 144+ kingdom-enriched (8-level) outputs, and which one to use.\n\n"
            "• Pipeline_log.txt\n"
            "  - Complete log of all preprocessing steps: quality trimming, pair merging,\n"
            "    read filtering, spike-in removal, and read count statistics.\n\n"
            "• Analysis_log.txt\n"
            "  - Complete log of all clustering, taxonomic classification, normalization,\n"
            "    and tree-building steps.\n\n"
            "• spike_mapping_file.csv (if a mapping file was provided)\n"
            "  - A summary table listing all processed samples, total weights, spike amounts,\n"
            "    and input paths.\n\n\n"
            "Quick Guide: Which file should I open first?\n"
            "-------------------------------------------\n"
            "1. For visual exploration:\n"
            "   Open any *.html file (e.g., OTUs-Krona.html or SOTUs-Krona.html) directly in\n"
            "   your web browser to explore your samples' taxonomy interactively.\n\n"
            "2. For statistical analysis in Excel / R / Python:\n"
            "   Open the primary *.tab abundance table inside your chosen pipeline folder\n"
            "   (e.g., OTUs-Table.tab, ZOTUs-Table.tab, OTUs-Table-TAC.tab, or\n"
            "   SOTUs-Table-TIC.tab).\n\n"
            "3. For functional prediction:\n"
            "   Use the matching files inside PICRUSt2-Inputs/ with PICRUSt2.\n"
            "===============================================================================\n"
        )


def write_pipeline_readme(
        output_dir: str,
        pipeline_type: str,
        filename: str = "README.txt") -> str:
    """
    Write a simple-language README.txt file into `output_dir` explaining what information
    each file contains and for what it is used.
    """
    dest_dir = Path(PurePath(output_dir)).resolve()
    dest_dir.mkdir(parents=True, exist_ok=True)
    readme_path = dest_dir / filename
    content = get_pipeline_readme_content(pipeline_type)
    with open(readme_path, "w", encoding="utf-8") as fio:
        fio.write(content)
    return str(readme_path)


def write_all_pipeline_readmes(
        analysis_dir: str,
        active_pipelines: Optional[Sequence[str]] = None,
        filename: str = "README.txt") -> List[str]:
    """
    Write README.txt documentation into the root analysis directory and each pipeline
    subdirectory (OTU_Pipeline, ZOTU_Pipeline, TAC_Pipeline, TIC_Pipeline, PICRUSt2-Inputs, Shared)
    that exists or was active.
    """
    root_dir = Path(PurePath(analysis_dir)).resolve()
    root_dir.mkdir(parents=True, exist_ok=True)
    created: List[str] = []

    # Write root README
    created.append(write_pipeline_readme(str(root_dir), "root", filename=filename))

    active_set = {
        str(p).strip().lower().replace("-", "_")
        for p in active_pipelines
    } if active_pipelines else set()

    sub_dirs = [
        ("otu", root_dir / "OTU_Pipeline"),
        ("zotu", root_dir / "ZOTU_Pipeline"),
        ("tac", root_dir / "TAC_Pipeline"),
        ("tic", root_dir / "TIC_Pipeline"),
        ("picrust2", root_dir / "PICRUSt2-Inputs"),
        ("shared", root_dir / "Shared"),
    ]

    for p_type, p_dir in sub_dirs:
        is_active = p_type in active_set
        has_content = False
        if p_dir.is_dir():
            has_content = any(
                item.is_file() and item.name != filename
                for item in p_dir.iterdir()
            )
        if is_active or has_content:
            p_dir.mkdir(parents=True, exist_ok=True)
            created.append(write_pipeline_readme(str(p_dir), p_type, filename=filename))

    return created


def onelinefasta(fastafilepath):
    proper_filepath = Path(fastafilepath).resolve()
    dirpath = proper_filepath.parent
    with tempfile.NamedTemporaryFile(
            delete=False,
            suffix=str(proper_filepath.suffix),
            dir=str(dirpath)) as temp_file:
        newfile_temp_name = temp_file.name
    with proper_filepath.open('r') as f, open(newfile_temp_name, "w+") as onelinefa:
        line = f.readline()
        sequence = ""
        while line:
            if line[0] == '>':
                if sequence:
                    onelinefa.write(sequence + "\n")
                    sequence = ""
                onelinefa.write(line)
            elif line == '\n' or line == '\r\n':
                pass
            else:
                sequence += line.strip()
            line = f.readline()
        if sequence:
            onelinefa.write(sequence + "\n")
    try:
        os.replace(newfile_temp_name, proper_filepath)
    finally:
        if Path(newfile_temp_name).exists():
            Path(newfile_temp_name).unlink()


def has_fasta_records(fasta_path: Union[str, Path]) -> bool:
    """
    Return True if fasta_path exists, is a regular file, and contains at least
    one sequence record starting with '>'.
    """
    p = Path(fasta_path).resolve()
    if not p.is_file() or p.stat().st_size == 0:
        return False
    with open(p, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                return True
    return False


def keep_seq_id_only(input_fasta: str, output_fasta: str) -> None:
    """
    Removes taxonomy information from fasta headers.
    :param input_fasta: Path to the input fasta file.
    :param output_fasta: Path to the output fasta file without taxonomy info.
    """
    tax_regex = re.compile(r"^>(?P<seq_id>[^\s]+)", re.IGNORECASE)

    with open(input_fasta, 'r') as infile, open(output_fasta, 'w+') as outfile:
        for line in infile:
            if line.startswith('>'):
                cleaned_line = tax_regex.match(line).group(0)
                outfile.write(cleaned_line + '\n')
            else:
                outfile.write(line)


def remove_taxonomy_from_table(input_table: str, output_table: str, sep: str = '\t') -> None:
    """
    Removes taxonomy information from the last column of a tab-separated table.
    :param input_table: Path to the input table file.
    :param output_table: Path to the output table file without taxonomy info.
    """
    with open(input_table, 'r') as infile, open(output_table, 'w+') as outfile:
        header = infile.readline().strip().lower()
        header_l = header.split(sep)
        tax_index = [i for i, col in enumerate(header_l) if col == 'taxonomy']
        if tax_index:
            # Taking everything except the column with taxonomy
            header = sep.join([col for i, col in enumerate(header_l) if i != tax_index[0]])
            header = header.strip()
        outfile.write(header + '\n')
        for line in infile:
            line = line.strip()
            line_l = line.split(sep)
            if tax_index:
                line = sep.join([col for i, col in enumerate(line_l) if i != tax_index[0]])
                line = line.strip()
            outfile.write(line + '\n')


def prepare_picrust2_input(analysis_dir: str) -> None:
    """
    It finds all OTU tables and sequnces fasta files and deletes taxonomy
    from them and only keeps the sequence id in the fastas.
    Any file following the naming pattern '.*[ZS]?OTUs?-Table.*.tab' will be
    considered as an OTU table.
    Any file following the naming pattern '.*[ZS]?OTUs?-Seqs.*.fasta' will be
    considered as a sequences fasta file.
    """
    analysis_dir = Path(PurePath(analysis_dir)).absolute()
    otu_table_pattern = re.compile(r"(?P<prefix>.*[ZS]?OTUs?-Table.*)\.tab$", re.IGNORECASE)
    otu_fasta_pattern = re.compile(r"(?P<prefix>.*[ZS]?OTUs?-Seqs.*)\.fasta$", re.IGNORECASE)
    picrust2_inputs_dir = analysis_dir / "PICRUSt2-Inputs"
    picrust2_inputs_dir.mkdir(exist_ok=True)

    for item in analysis_dir.rglob("*"):
        if item.is_file():
            if picrust2_inputs_dir in item.parents:
                continue
            if _EXPLICIT_LEVEL_SUFFIX_RE.search(item.name):
                continue
            match_table = otu_table_pattern.match(item.name)
            match_fasta = otu_fasta_pattern.match(item.name)

            if match_table:
                prefix = match_table.group("prefix")
                otu_table_file = item
                rel_parent = item.parent.relative_to(analysis_dir)
                rel_parent_str = str(rel_parent).replace(os.sep, "__")
                rel_parent_str = rel_parent_str.strip(".")
                if rel_parent_str:
                    prefix = f"{rel_parent_str}__{prefix}"
                picrust2_otu_table_file = picrust2_inputs_dir / f"{prefix}-For-PICRUSt2.tab"
                remove_taxonomy_from_table(otu_table_file, picrust2_otu_table_file)

            elif match_fasta:
                prefix = match_fasta.group("prefix")
                otu_fasta_file = item
                rel_parent = item.parent.relative_to(analysis_dir)
                rel_parent_str = str(rel_parent).replace(os.sep, "__")
                rel_parent_str = rel_parent_str.strip(".")
                if rel_parent_str:
                    prefix = f"{rel_parent_str}__{prefix}"
                picrust2_otu_fasta_file = picrust2_inputs_dir / f"{prefix}-For-PICRUSt2.fasta"
                keep_seq_id_only(otu_fasta_file, picrust2_otu_fasta_file)
                onelinefasta(picrust2_otu_fasta_file)

    if any(picrust2_inputs_dir.glob("*-For-PICRUSt2.*")):
        write_pipeline_readme(str(picrust2_inputs_dir), "picrust2")
