import csv
import os
import re
from pathlib import Path, PurePath
import tempfile
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple


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
    (index 1) of a classic 7-level taxonomy string (`Domain;Phylum;Class;Order;Family;Genus;Species;`),
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
    target_levels: int,
    seq_kingdom_map: Optional[Dict[str, str]] = None,
    phylum_kingdom_map: Optional[Dict[Tuple[str, str], str]] = None,
    krona_tool: str = "",
    create_krona_fn: Optional[Callable[[str, str, str], str]] = None,
) -> List[str]:
    """
    For a finished pipeline directory (`ZOTU_Pipeline`, `OTU_Pipeline`,
    `TAC_Pipeline`, or `TIC_Pipeline`), generate the Kingdom-enriched version
    (`*-7Levels.*` for ZOTU/OTU or `*-8Levels.*` for TAC/TIC) of all output
    tables, FASTA files, and Krona HTML charts that contain taxonomy.
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
        elif item.suffix.lower() == ".fasta" and not item.name.startswith(("test_", "shortened_algn_", "filtered_by_abundance_")):
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

    if create_krona_fn and krona_tool and primary_table_for_krona and primary_table_for_krona.is_file():
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
    between standard outputs and the SILVA >= 144 *-7Levels and *-8Levels files.
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
        "1. Standard Outputs (e.g., SOTUs-Table-TIC.tab, OTUs-Table-TAC.tab, ZOTUs-Table.tab):\n"
        "   - Format: Domain;Phylum;Class;Order;Family;Genus;Species/SOTU; (7 levels, 7 ';')\n"
        "   - Purpose: Direct compatibility with legacy tools, QIIME, PICRUSt2, and R\n"
        "     packages that strictly expect traditional 7-level taxonomies without Kingdom.\n\n"
        "2. Enriched Output Files:\n"
        "   - '*-8Levels.*' (generated for TAC & TIC pipelines):\n"
        "       Format: Domain;Kingdom;Phylum;Class;Order;Family;Genus;Species/SOTU; (8 levels, 8 ';')\n"
        "       Includes the bacterial kingdom rank between Domain and Phylum.\n"
        "   - '*-7Levels.*' (generated for ZOTU & OTU pipelines, and -iz individual outputs):\n"
        "       Format: Domain;Kingdom;Phylum;Class;Order;Family;Genus; (7 levels, 7 ';')\n"
        "       Reflects the full taxonomic resolution available without species/SOTU clustering.\n\n"
        "Which file should I use?\n"
        "------------------------\n"
        "• Use Standard files if your downstream scripts (e.g. PICRUSt2, legacy phyloseq\n"
        "  scripts) require classic 7-level taxonomies.\n"
        "• Use *-8Levels / *-7Levels files if you want the most up-to-date SILVA 144+\n"
        "  taxonomic nomenclature including bacterial kingdoms (e.g. for Rhea, Namco, or Krona).\n"
        "===============================================================================\n"
    )
    with open(readme_path, "w", encoding="utf-8") as fio:
        fio.write(content)
    return str(readme_path)


def onelinefasta(fastafilepath):
    proper_filepath = Path(fastafilepath).resolve()
    dirpath = proper_filepath.parent
    filename = proper_filepath.name
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
    Any file following the naming pattern '.*[ZS]?OTUs?-Table.*.tab' will be considered as an OTU table
    Any file following the naming pattern '.*[ZS]?OTUs?-Seqs.*.fasta' will be considered as a sequences fasta file
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

