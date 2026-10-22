import datetime
from email import header
from collections import Counter
import csv

from .seqkit import *
from .usearch import *
from pathlib import Path
import multiprocessing
import sys
import logging
import shutil
from itertools import zip_longest
from typing import Iterator, List, Optional, Tuple
from matplotlib.ticker import AutoMinorLocator

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator
except Exception:
    plt = None
    AutoMinorLocator = None

try:
    import numpy as np
except Exception:
    np = None

PLACEHOLDER_TAXA = {
    "uncultured",
    "unclassified",
    "metagenome",
    "candidate",
    "environmental",
    "bacterium",
    "bacteria",
    "sp.",
    "sp",
    "strain",
    "norank",
    "unknown",
}

# Configure logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.FileHandler('silva_algn_cut_and_count.log', mode='w'),
        logging.StreamHandler(sys.stdout)
    ]
)
# Set stdout handler to INFO level only
logging.getLogger().handlers[1].setLevel(logging.INFO)
logger = logging.getLogger(__name__)


SILVA_ALGN_PATHS = [
    # Complete SILVA DB - 0
    "./DBs/SILVA/SILVA_Align/SILVA_138.2_SSURef_tax_silva_full_align_trunc.fasta",
    # NR 99 SILVA DB - 1
    "./DBs/SILVA/SILVA_Align/SILVA_138.2_SSURef_NR99_tax_silva_full_align_trunc.fasta",
    # NR 99 Full Known Bacteria Taxonomy - 2
    "./DBs/SILVA/SILVA_Align/SILVA_138.2_SSURef_NR99_tax_silva_full_align_trunc_full_known_taxonomy.fasta",
    # LTP - 3
    "./DBs/LTP/LTP_10_2024_SILVA_align_DNA_OnlyBacteria_all7.fasta",
    # GG2 - 4
    "./DBs/GG2/GG2_FullSpeciesTax_SILVA_align.fasta",
]


def _iter_fasta_records(fasta_path: Path) -> Iterator[Tuple[str, str]]:
    """Yield FASTA records as (header_without_>, sequence)."""
    header = None
    seq_chunks: list[str] = []
    # Some reference FASTAs contain non-UTF8 bytes in headers/comments.
    # Use replacement decoding so parsing does not fail mid-run.
    with open(fasta_path, "r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)
    if header is not None:
        yield header, "".join(seq_chunks)


def _verify_gapless_version(aligned_fasta: Path, nogap_fasta: Path) -> tuple[bool, str]:
    """
    Validate that `nogap_fasta` is a true gapless rendering of `aligned_fasta`.
    Conditions:
      - same record count and order
      - same headers
      - no '-' or '.' in nogap sequences
      - seq_aligned with '-' and '.' removed equals seq_nogap
    """
    rec_idx = 0
    for rec_aln, rec_nogap in zip_longest(
        _iter_fasta_records(aligned_fasta),
        _iter_fasta_records(nogap_fasta),
        fillvalue=None,
    ):
        rec_idx += 1
        if rec_aln is None or rec_nogap is None:
            return False, f"record count mismatch at record {rec_idx}"
        h_aln, s_aln = rec_aln
        h_nogap, s_nogap = rec_nogap
        if h_aln != h_nogap:
            return False, f"header mismatch at record {rec_idx}: '{h_aln}' != '{h_nogap}'"
        if ("-" in s_nogap) or ("." in s_nogap):
            return False, f"nogap sequence still contains gap chars at record {rec_idx} ('{h_nogap}')"
        if s_aln.replace("-", "").replace(".", "") != s_nogap:
            return False, f"sequence mismatch after ungap at record {rec_idx} ('{h_nogap}')"
    return True, "ok"


def _nogap_path_for(fasta_path: Path) -> Path:
    suffix = fasta_path.suffix if fasta_path.suffix else ".fasta"
    return fasta_path.with_name(f"{fasta_path.stem}_nogap{suffix}")


def _generate_gapless_fasta(input_fasta: Path, output_fasta: Path) -> Path:
    logger.info(f"Generating gapless FASTA: {input_fasta} -> {output_fasta}")
    ret_code, _ = seqkit_seq(
        input_fasta=input_fasta,
        output_fasta=output_fasta,
        args_list=["--remove-gaps", "-j", "10", "-w", "0"],
    )
    if ret_code != 0:
        raise RuntimeError(f"Failed to generate gapless FASTA from {input_fasta}")
    return output_fasta


def _ensure_verified_gapless(aligned_fasta: Path) -> Path:
    """
    Return a verified gapless FASTA path next to aligned_fasta, named *_nogap.fasta.
    If missing, generate it. If existing but invalid, regenerate and validate again.
    """
    nogap_fasta = _nogap_path_for(aligned_fasta)
    if not nogap_fasta.exists():
        _generate_gapless_fasta(aligned_fasta, nogap_fasta)
    else:
        ok, msg = _verify_gapless_version(aligned_fasta, nogap_fasta)
        if not ok:
            logger.warning(f"Existing nogap FASTA failed verification ({msg}); regenerating: {nogap_fasta}")
            _generate_gapless_fasta(aligned_fasta, nogap_fasta)

    ok, msg = _verify_gapless_version(aligned_fasta, nogap_fasta)
    if not ok:
        raise RuntimeError(f"Gapless FASTA verification failed for {nogap_fasta}: {msg}")
    logger.info(f"Verified gapless FASTA: {nogap_fasta}")
    return nogap_fasta


def _normalize_cluster_ids(clust_id: float, clust_ids: Optional[list[float]]) -> list[float]:
    ids = [float(clust_id)] if clust_ids is None else [float(x) for x in clust_ids]
    out: list[float] = []
    seen: set[float] = set()
    for cid in ids:
        if cid <= 0.0 or cid > 1.0:
            raise ValueError(f"Invalid cluster id '{cid}'. Expected 0 < id <= 1.")
        cid_r = round(cid, 6)
        if cid_r in seen:
            continue
        seen.add(cid_r)
        out.append(cid_r)
    if not out:
        raise ValueError("No valid cluster IDs provided.")
    return out


def _cluster_id_label(cluster_id: float) -> str:
    return str(cluster_id).replace(".", "p")


def _to_int(value, default: int = 0) -> int:
    try:
        value_str = str(value).replace(",", "").strip()
        if not value_str:
            return default
        return int(float(value_str))
    except Exception:
        return default


def _to_float(value, default: float = 0.0) -> float:
    try:
        value_str = str(value).replace(",", "").strip()
        if not value_str:
            return default
        return float(value_str)
    except Exception:
        return default


def _extract_taxonomy_from_header(header: str, taxlevel: str = "species", return_full: bool = False) -> str:
    parts = str(header).strip().split(maxsplit=1)
    if len(parts) < 2:
        return ""
    taxonomy = parts[1].strip()
    levels = [lvl.strip().lower() for lvl in taxonomy.split(";") if lvl.strip()]
    if not levels:
        return ""
    if taxlevel == "full":
        return ";".join(levels)
    if taxlevel == "species":
        if len(levels) != 7:
            return ""
        species_tax, _ = _validate_species_taxonomy(header)
        if not species_tax:
            return ""
        levels[6] = species_tax
        # Returns more unique taxa when using the full taxonomy string,
        #  but it may be more noisy and less interpretable than just species-level taxonomy.
        # Depending on the use case, one may prefer one over the other.
        if return_full:
            return ";".join(levels)
        else:
            return levels[6]
    return ""


def _normalize_species_taxonomy(header: str) -> str:
    species_tax, _ = _validate_species_taxonomy(header)
    return species_tax


def _strict_taxonomy_from_header(header: str) -> tuple[str, str]:
    """
    Return taxonomy in strict 7-level form where the last level is exactly
    'Genus species' (extra tokens after species are removed).
    """
    parts = str(header).strip().split(maxsplit=1)
    if len(parts) < 2:
        return "", "missing taxonomy"

    taxonomy = parts[1].strip()
    levels_raw = [lvl.strip() for lvl in taxonomy.split(";") if lvl.strip()]
    if len(levels_raw) != 7:
        return "", f"expected 7 taxonomy levels, found {len(levels_raw)}"

    genus = levels_raw[5]
    species_field = levels_raw[6]
    species_tokens = [tok.strip(".") for tok in species_field.split() if tok.strip()]
    if len(species_tokens) < 2:
        return "", "species rank is not in 'Genus species' format"

    genus_species = species_tokens[:2]
    genus_species_lower = [tok.lower() for tok in genus_species]
    if any(tok in PLACEHOLDER_TAXA for tok in genus_species_lower):
        return "", "species rank contains placeholder taxonomy"

    # Normalize the last level to exactly two tokens.
    # Keep genus from level-6 to enforce '...;Genus;Genus species'.
    levels_raw[6] = f"{genus} {species_tokens[1]}"
    return ";".join(levels_raw), ""


def _validate_species_taxonomy(header: str) -> tuple[str, str]:
    strict_taxonomy, reason = _strict_taxonomy_from_header(header)
    if not strict_taxonomy:
        return "", reason
    levels = [lvl.strip().lower() for lvl in strict_taxonomy.split(";") if lvl.strip()]
    return levels[6], ""


def _filter_fasta_by_valid_species(
    input_fasta: Path,
    output_fasta: Path,
    reuse: bool,
    discarded_file: Optional[Path] = None,
) -> Path:
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    if reuse and output_fasta.exists():
        logger.info(f"Reusing filtered taxonomy FASTA: {output_fasta}")
        if discarded_file and discarded_file.exists():
            logger.info(f"Reusing discarded taxonomy report: {discarded_file}")
        return output_fasta

    kept = 0
    discarded = 0
    if discarded_file:
        discarded_file.parent.mkdir(parents=True, exist_ok=True)
        with open(output_fasta, "w") as out, open(discarded_file, "w") as discarded_out:
            discarded_out.write("header\treason\n")
            for header, seq in _iter_fasta_records(input_fasta):
                species_tax, reason = _validate_species_taxonomy(header)
                if not species_tax:
                    discarded += 1
                    discarded_out.write(f"{header}\t{reason}\n")
                    continue
                seq_id = str(header).split(maxsplit=1)[0]
                strict_taxonomy, _ = _strict_taxonomy_from_header(header)
                out.write(f">{seq_id} {strict_taxonomy}\n")
                out.write(f"{seq}\n")
                kept += 1
    else:
        with open(output_fasta, "w") as out:
            for header, seq in _iter_fasta_records(input_fasta):
                species_tax, _ = _validate_species_taxonomy(header)
                if not species_tax:
                    discarded += 1
                    continue
                seq_id = str(header).split(maxsplit=1)[0]
                strict_taxonomy, _ = _strict_taxonomy_from_header(header)
                out.write(f">{seq_id} {strict_taxonomy}\n")
                out.write(f"{seq}\n")
                kept += 1
    if kept == 0:
        raise RuntimeError(f"No sequences with valid species taxonomy found in {input_fasta}")
    logger.info(f"Filtered {kept} sequences with valid species taxonomy into {output_fasta}")
    if discarded_file:
        logger.info(f"Wrote {discarded} discarded records to {discarded_file}")
    return output_fasta


def _count_unique_taxonomies(fasta_path: Path) -> int:
    taxa = set()
    for header, _ in _iter_fasta_records(fasta_path):
        strict_taxonomy, _ = _strict_taxonomy_from_header(header)
        if strict_taxonomy:
            taxa.add(strict_taxonomy.strip().lower())
    return len(taxa)


def _load_stats_tsv(stats_path: Path) -> dict:
    if not stats_path.exists():
        return {}
    with open(stats_path, "r") as stats_file:
        lines = [line.strip() for line in stats_file if line.strip()]
    if len(lines) < 2:
        return {}
    headers = lines[0].split()
    values = lines[1].split()
    return dict(zip(headers, values))


def _strict_taxonomy_key_from_header(header: str) -> str:
    strict_taxonomy, _ = _strict_taxonomy_from_header(header)
    if not strict_taxonomy:
        return ""
    return strict_taxonomy.strip().lower()


def _count_full_taxonomy_labels_from_header_list(header_list_path: Path) -> Counter:
    """
    Count strict full-taxonomy labels from a header list.
    Labels are normalized to strict 7-level taxonomy and lowercased.
    """
    counts: Counter[str] = Counter()
    if not header_list_path.exists():
        logger.warning(f"Header list not found: {header_list_path}")
        return counts

    with open(header_list_path, "r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                line = line[1:]
            tax_key = _strict_taxonomy_key_from_header(line)
            if not tax_key:
                continue
            if "bacteria" not in tax_key:
                continue
            counts[tax_key] += 1
    return counts


def _count_full_taxonomy_labels_from_fasta(fasta_path: Path) -> Counter:
    """Count strict full-taxonomy labels directly from FASTA headers."""
    counts: Counter[str] = Counter()
    if not fasta_path.exists():
        logger.warning(f"FASTA not found for taxonomy counting: {fasta_path}")
        return counts
    for header, _ in _iter_fasta_records(fasta_path):
        tax_key = _strict_taxonomy_key_from_header(header)
        if not tax_key:
            continue
        if "bacteria" not in tax_key:
            continue
        counts[tax_key] += 1
    return counts


def _plot_species_frequency_histogram(
    counts: Counter[str],
    output_dir: Path,
    bin_size: int = 100,
) -> Optional[Path]:
    """
    Plot histogram of all species-label frequencies with an overlaid smoothed line.
    """
    if plt is None:
        logger.warning("matplotlib is unavailable; skipping species frequency histogram.")
        return None
    if bin_size <= 0:
        logger.warning(f"Invalid histogram bin size ({bin_size}); falling back to 100.")
        bin_size = 100

    freq_values = [int(v) for v in counts.values() if int(v) > 0]
    if not freq_values:
        logger.warning("No positive species frequencies found; skipping histogram.")
        return None

    max_freq = max(freq_values)
    n_bins = max(1, (max_freq // bin_size) + 1)
    bin_edges = [i * bin_size for i in range(n_bins + 1)]
    hist_counts = [0] * n_bins
    for val in freq_values:
        idx = min(int(val // bin_size), n_bins - 1)
        hist_counts[idx] += 1
    bin_centers = [edge + (bin_size / 2.0) for edge in bin_edges[:-1]]

    smooth_counts = [float(v) for v in hist_counts]
    smooth_window = min(7, len(smooth_counts))
    if smooth_window % 2 == 0:
        smooth_window -= 1
    if smooth_window >= 3:
        pad = smooth_window // 2
        padded = [smooth_counts[0]] * pad + smooth_counts + [smooth_counts[-1]] * pad
        kernel = [1.0 / smooth_window] * smooth_window
        smoothed: list[float] = []
        for i in range(len(smooth_counts)):
            smoothed.append(
                sum(padded[i + k] * kernel[k] for k in range(smooth_window))
            )
        smooth_counts = smoothed

    fig, ax = plt.subplots(figsize=(14, 8), dpi=220)
    ax.bar(
        bin_centers,
        hist_counts,
        width=bin_size * 0.92,
        color="#F0E442",
        edgecolor="#555555",
        linewidth=0.5,
        alpha=0.90,
        label=f"Histogram (bin={bin_size})",
    )
    ax.plot(
        bin_centers,
        smooth_counts,
        color="#0072B2",
        linewidth=2.8,
        label="Smoothed trend",
        zorder=5,
    )
    ax.scatter(
        bin_centers,
        smooth_counts,
        color="#0072B2",
        s=14,
        alpha=0.75,
        zorder=6,
    )

    ax.set_title("Histogram of Species Frequencies", fontsize=18)
    ax.set_xlabel("Frequency per species taxonomy label", fontsize=14)
    ax.set_ylabel("Number of species taxonomy labels", fontsize=14)
    ax.tick_params(axis="both", labelsize=12)
    ax.set_xlim(left=0)
    ax.grid(True, axis="y", linestyle=":", alpha=0.35)
    ax.legend(loc="upper right", frameon=False, fontsize=12)
    fig.subplots_adjust(left=0.10, right=0.98, top=0.92, bottom=0.12)

    out_png = output_dir / "gg2_species_frequency_histogram.png"
    out_svg = output_dir / "gg2_species_frequency_histogram.svg"
    fig.savefig(out_png, bbox_inches="tight", pad_inches=0.15)
    fig.savefig(out_svg, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
    logger.info(f"Species frequency histogram saved to: {out_png}")
    logger.info(f"Species frequency histogram saved to: {out_svg}")
    return out_png


def _collect_v1v9_unique_sequences_by_taxonomy(v1v9_nogap_fasta: Path) -> dict[str, set[str]]:
    """Map strict taxonomy -> unique V1-V9 ungapped sequences."""
    by_taxonomy: dict[str, set[str]] = {}
    if not v1v9_nogap_fasta.exists():
        logger.warning(f"V1-V9 nogap FASTA not found: {v1v9_nogap_fasta}")
        return by_taxonomy
    for header, seq in _iter_fasta_records(v1v9_nogap_fasta):
        tax_key = _strict_taxonomy_key_from_header(header)
        if not tax_key:
            continue
        if "bacteria" not in tax_key:
            continue
        by_taxonomy.setdefault(tax_key, set()).add(seq)
    return by_taxonomy


def _cluster_unique_sequences_for_taxonomy(
    unique_sequences: list[str],
    clust_ids: list[float],
    output_dir: Path,
    rank_idx: int,
    reuse: bool,
) -> dict[str, int]:
    out: dict[str, int] = {}
    cluster_cols = [f"Number of Clusters at {cid}" for cid in clust_ids]
    if not unique_sequences:
        for col in cluster_cols:
            out[col] = 0
        return out

    taxon_dir = output_dir / "species_cluster_work" / f"top_{rank_idx:02d}"
    taxon_dir.mkdir(parents=True, exist_ok=True)
    unique_fasta = taxon_dir / "unique_v1-v9_nogap.fasta"
    if (not reuse) or (not unique_fasta.exists()):
        with open(unique_fasta, "w") as out_fa:
            for seq_idx, sequence in enumerate(unique_sequences, start=1):
                out_fa.write(f">top{rank_idx:02d}_seq{seq_idx}\n{sequence}\n")
    cluster_counts = _cluster_counts_for_row(str(unique_fasta), clust_ids)
    for col in cluster_cols:
        out[col] = _to_int(cluster_counts.get(col, 0), 0)
    return out


def plot_species_label_frequency(
    header_list_path: Optional[str | Path],
    output_dir: str | Path,
    top_n: int = 25,
    clust_ids: Optional[list[float]] = None,
    v1v9_nogap_fasta: Optional[str | Path] = None,
    reuse_intermediates: bool = False,
    cluster_colors: Optional[dict] = None,
    unique_bar_color: str = "#CC79A7",
) -> Optional[Path]:
    """
    Draw grouped bars for top-N strict taxonomy labels:
    - frequency in input headers
    - unique V1-V9 sequence count
    - cluster count per threshold on unique V1-V9 sequences
    Also draw an average-frequency horizontal line.
    """
    if plt is None:
        logger.warning("matplotlib is unavailable; skipping species-label frequency plot.")
        return None

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if clust_ids is None:
        _, clust_ids = _read_stats_rows_from_tsv(output_dir / "silva_amplicon_overall_stats.tsv")
    clust_ids = clust_ids or []

    counts: Counter[str]
    if header_list_path is not None and Path(header_list_path).exists():
        counts = _count_full_taxonomy_labels_from_header_list(Path(header_list_path))
    else:
        fallback_fasta = output_dir / "silva_full_length_valid_taxonomy_aligned.fasta"
        counts = _count_full_taxonomy_labels_from_fasta(fallback_fasta)
    if not counts:
        logger.warning("No strict taxonomy labels were counted for top-N frequency plot.")
        return None
    _ = _plot_species_frequency_histogram(
        counts=counts,
        output_dir=output_dir,
        bin_size=100,
    )

    if v1v9_nogap_fasta is None:
        v1v9_nogap_fasta = output_dir / "silva_v1-v9_alignment_1250-1650_filtered_nogap.fasta"
    v1v9_nogap_fasta = Path(v1v9_nogap_fasta)
    v1v9_unique_by_tax = _collect_v1v9_unique_sequences_by_taxonomy(v1v9_nogap_fasta)

    top_labels = counts.most_common(top_n)
    labels: list[str] = []
    freq_vals: list[int] = []
    summary_tsv = output_dir / "gg2_top25_species_label_frequency_with_clusters.tsv"
    cluster_cols = [f"Number of Clusters at {cid}" for cid in clust_ids]
    cached_rows: dict[str, dict] = {}
    if reuse_intermediates and summary_tsv.exists():
        try:
            with open(summary_tsv, "r", encoding="utf-8", errors="replace") as in_tsv:
                reader = csv.DictReader(in_tsv, delimiter="\t")
                required = {"taxonomy_label", "frequency", "unique_v1v9_sequences", *cluster_cols}
                if reader.fieldnames and required.issubset(set(reader.fieldnames)):
                    for row in reader:
                        key = str(row.get("taxonomy_label", "")).strip().lower()
                        if key:
                            cached_rows[key] = row
                    logger.info(f"Reusing species frequency summary table: {summary_tsv}")
        except Exception as e:
            logger.warning(f"Failed reading cached species summary {summary_tsv}: {e}")

    unique_vals: list[int] = []
    cluster_series: dict[float, list[int]] = {cid: [] for cid in clust_ids}
    summary_rows: list[dict] = []

    for rank_idx, (tax_label, freq) in enumerate(top_labels, start=1):
        display_taxon = tax_label.split(";")[-1] if ";" in tax_label else tax_label
        genus, separator, remainder = display_taxon.partition(" ")
        if genus:
            display_taxon = genus[:1].upper() + genus[1:]
            if separator:
                display_taxon += separator + remainder
        labels.append(display_taxon)
        freq_vals.append(int(freq))
        row = {
            "taxonomy_label": tax_label,
            "frequency": int(freq),
            "unique_v1v9_sequences": 0,
        }
        cache_key = tax_label.strip().lower()
        cached = cached_rows.get(cache_key)
        if cached:
            row["unique_v1v9_sequences"] = _to_int(cached.get("unique_v1v9_sequences", 0), 0)
            for cid in clust_ids:
                col = f"Number of Clusters at {cid}"
                row[col] = _to_int(cached.get(col, 0), 0)
        else:
            unique_seqs = sorted(v1v9_unique_by_tax.get(tax_label, set()))
            row["unique_v1v9_sequences"] = int(len(unique_seqs))
            cluster_counts = _cluster_unique_sequences_for_taxonomy(
                unique_sequences=unique_seqs,
                clust_ids=clust_ids,
                output_dir=output_dir,
                rank_idx=rank_idx,
                reuse=reuse_intermediates,
            )
            for cid in clust_ids:
                col = f"Number of Clusters at {cid}"
                row[col] = _to_int(cluster_counts.get(col, 0), 0)

        unique_vals.append(_to_int(row["unique_v1v9_sequences"], 0))
        for cid in clust_ids:
            col = f"Number of Clusters at {cid}"
            count_val = _to_int(row.get(col, 0), 0)
            cluster_series[cid].append(count_val)
            row[col] = count_val
        summary_rows.append(row)

    # Add a 26th grouped column that treats top-N taxa as one pooled entity.
    combined_tax_key = f"__top_{top_n}_combined__"
    combined_display_label = f"Top {top_n} combined"
    combined_freq = int(sum(int(freq) for _, freq in top_labels))

    combined_row = {
        "taxonomy_label": combined_tax_key,
        "frequency": combined_freq,
        "unique_v1v9_sequences": 0,
    }
    combined_cached = cached_rows.get(combined_tax_key.lower())
    if combined_cached:
        combined_row["unique_v1v9_sequences"] = _to_int(
            combined_cached.get("unique_v1v9_sequences", 0), 0
        )
        for cid in clust_ids:
            col = f"Number of Clusters at {cid}"
            combined_row[col] = _to_int(combined_cached.get(col, 0), 0)
    else:
        combined_unique_set: set[str] = set()
        for tax_label, _freq in top_labels:
            combined_unique_set.update(v1v9_unique_by_tax.get(tax_label, set()))
        combined_unique_seqs = sorted(combined_unique_set)
        combined_row["unique_v1v9_sequences"] = int(len(combined_unique_seqs))
        combined_cluster_counts = _cluster_unique_sequences_for_taxonomy(
            unique_sequences=combined_unique_seqs,
            clust_ids=clust_ids,
            output_dir=output_dir,
            rank_idx=top_n + 1,
            reuse=reuse_intermediates,
        )
        for cid in clust_ids:
            col = f"Number of Clusters at {cid}"
            combined_row[col] = _to_int(combined_cluster_counts.get(col, 0), 0)

    # Put combined bars at the left-most position in the plot.
    labels.insert(0, combined_display_label)
    freq_vals.insert(0, combined_freq)
    unique_vals.insert(0, _to_int(combined_row["unique_v1v9_sequences"], 0))
    for cid in clust_ids:
        col = f"Number of Clusters at {cid}"
        count_val = _to_int(combined_row.get(col, 0), 0)
        cluster_series[cid].insert(0, count_val)
        combined_row[col] = count_val
    summary_rows.insert(0, combined_row)

    summary_headers = [
        "taxonomy_label",
        "frequency",
        "unique_v1v9_sequences",
        *[f"Number of Clusters at {cid}" for cid in clust_ids],
    ]
    with open(summary_tsv, "w") as out_tsv:
        out_tsv.write("\t".join(summary_headers) + "\n")
        for row in summary_rows:
            out_tsv.write("\t".join(str(row.get(h, "")) for h in summary_headers) + "\n")

    y_step = 4.80
    y_positions = [float(i) * y_step for i in range(len(labels))]
    series_specs: list[tuple[str, list[int], str]] = []
    # Same color-blind-safe palette family as plot_fig01.py
    series_specs.append(("Unique V1-V9 sequences", unique_vals, unique_bar_color))
    cluster_palette = ["#56B4E9", "#009E73", "#E69F00", "#D55E00", "#999999"]
    for idx, cid in enumerate(clust_ids):
        chosen_color = cluster_palette[idx % len(cluster_palette)]
        if cluster_colors:
            if cid in cluster_colors:
                chosen_color = str(cluster_colors[cid])
            elif str(cid) in cluster_colors:
                chosen_color = str(cluster_colors[str(cid)])
        series_specs.append((f"Clusters at {cid}", cluster_series.get(cid, []), chosen_color))

    n_series = max(1, len(series_specs))
    group_width = 4.80
    bar_width = group_width / n_series
    all_max_vals = [max(vals) for _, vals, _ in series_specs if vals]
    x_max = float(max(all_max_vals)) if all_max_vals else 0.0
    # Breakpoint from the average unique-sequence count across all labels.
    avg_unique = (float(sum(unique_vals)) / float(len(unique_vals))) if unique_vals else 1.0
    x_break = max(1.0, avg_unique)
    use_broken_axis = x_max > (x_break * 1.01)
    # Let the figure grow with label count for readability.
    n_labels = max(1, len(labels))
    fig_height = max(34.0, float(n_labels) * 2.4)
    if use_broken_axis:
        fig, (ax_left, ax_right) = plt.subplots(
            1,
            2,
            sharey=True,
            figsize=(34, fig_height),
            dpi=220,
            gridspec_kw={"width_ratios": [3.0, 1.1], "wspace": 0.06},
        )
    else:
        fig, ax_left = plt.subplots(figsize=(30, fig_height), dpi=220)
        ax_right = None
    # Requested readability bump: +50% font sizes for axis/tick/bar/legend labels.
    axis_label_fs = 46
    tick_label_fs = 40
    bar_value_label_fs = 44
    legend_fs = 36
    label_offset = max(x_max * 0.012, 1.0) if x_max > 0 else 1.0
    label_offset_low = max(x_break * 0.012, 2.0)
    label_offset_high = max((max(x_max, x_break) - x_break) * 0.012, 4.0)

    for series_idx, (series_label, series_vals, color) in enumerate(series_specs):
        offset = (series_idx - (n_series - 1) / 2.0) * bar_width
        y_shifted = [y + offset for y in y_positions]
        bars_left = ax_left.barh(
            y_shifted,
            series_vals,
            height=bar_width * 1.00,
            color=color,
            alpha=0.92,
            label=series_label,
        )
        bars_right = None
        if use_broken_axis and ax_right is not None:
            overflow_vals = [max(float(v) - x_break, 0.0) for v in series_vals]
            bars_right = ax_right.barh(
                y_shifted,
                overflow_vals,
                left=x_break,
                height=bar_width * 1.00,
                color=color,
                alpha=0.92,
                label="_nolegend_",
            )

        for idx, (bar, v) in enumerate(zip(bars_left, series_vals)):
            target_ax = ax_left
            target_bar = bar
            this_offset = label_offset
            if use_broken_axis and ax_right is not None and bars_right is not None:
                if float(v) > (x_break - label_offset_low):
                    target_ax = ax_right
                    target_bar = bars_right[idx]
                    this_offset = label_offset_high
                else:
                    this_offset = label_offset_low
            label = str(_to_int(v, 0))
            value = float(v)
            y_center = target_bar.get_y() + target_bar.get_height() / 2.0
            # Keep labels readable: place inside longer bars and outside only for short bars.
            if value > (this_offset * 2.5):
                text_x = value - (this_offset * 0.55)
                text_ha = "right"
                text_bbox = {
                    "boxstyle": "round,pad=0.12",
                    "facecolor": "white",
                    "alpha": 0.55,
                    "linewidth": 0.0,
                }
            else:
                text_x = value + (this_offset * 0.55)
                text_ha = "left"
                text_bbox = None
            if target_ax is ax_right:
                # Keep labels inside visible right panel when logarithmic scaling is active.
                min_visible = x_break * 1.01
                if text_x <= min_visible:
                    text_x = max(value * 1.01, min_visible)
                    text_ha = "left"
                    text_bbox = None
            else:
                # Prevent labels from crossing the y-axis boundary on very short bars.
                left_safe = x_break * 0.06
                if text_ha == "right" and text_x < left_safe:
                    text_x = max(value + (this_offset * 0.55), left_safe)
                    text_ha = "left"
                    text_bbox = None
            target_ax.text(
                text_x,
                y_center,
                label,
                ha=text_ha,
                va="center",
                fontsize=bar_value_label_fs,
                rotation=0,
                clip_on=False,
                bbox=text_bbox,
            )

    ax_left.set_xlabel("Feature count", fontsize=axis_label_fs)
    ax_left.set_ylabel("Top taxonomy labels", fontsize=axis_label_fs)
    if use_broken_axis and ax_right is not None:
        fig.suptitle(f"Top {top_n} Frequent Taxonomy Labels + V1-V9 Unique/Cluster Counts", fontsize=38)
    else:
        ax_left.set_title(f"Top {top_n} Frequent Taxonomy Labels + V1-V9 Unique/Cluster Counts", fontsize=38)

    if x_max > 0:
        if use_broken_axis and ax_right is not None:
            # Right panel uses logarithmic scale to compress high counts.
            upper_max = max(x_max * 1.06, x_break * 1.5)
            ax_left.set_xlim(0, x_break)
            ax_right.set_xlim(x_break, upper_max)
            ax_right.set_xscale("log")
        else:
            ax_left.set_xlim(0, x_max * 1.10)

    ax_left.set_yticks(y_positions)
    tick_texts = ax_left.set_yticklabels(labels, fontsize=tick_label_fs)
    for tick_text, label in zip(tick_texts, labels):
        if label != combined_display_label:
            tick_text.set_fontstyle("italic")
    ax_left.tick_params(axis="x", labelsize=tick_label_fs)
    ax_left.margins(y=0.02)
    ax_left.grid(True, axis="x", linestyle=":", alpha=0.35)
    ax_left.invert_yaxis()
    if use_broken_axis and ax_right is not None:
        ax_right.grid(True, axis="x", linestyle=":", alpha=0.35)
        ax_right.tick_params(axis="x", labelsize=tick_label_fs)
        ax_right.tick_params(axis="y", left=False, labelleft=False)
        ax_left.spines["right"].set_visible(False)
        ax_right.spines["left"].set_visible(False)

        d = 0.008
        kwargs = dict(transform=ax_left.transAxes, color="k", clip_on=False, linewidth=1.2)
        ax_left.plot((1 - d, 1 + d), (-d, +d), **kwargs)
        ax_left.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)
        kwargs.update(transform=ax_right.transAxes)
        ax_right.plot((-d, +d), (-d, +d), **kwargs)
        ax_right.plot((-d, +d), (1 - d, 1 + d), **kwargs)

    handles, legend_labels = ax_left.get_legend_handles_labels()
    if use_broken_axis and ax_right is not None:
        ax_right.legend(
            handles,
            legend_labels,
            loc="center right",
            bbox_to_anchor=(0.98, 0.5),
            frameon=False,
            fontsize=legend_fs,
            ncol=1,
        )
    else:
        ax_left.legend(
            handles,
            legend_labels,
            loc="center right",
            bbox_to_anchor=(0.98, 0.5),
            frameon=False,
            fontsize=legend_fs,
            ncol=1,
        )
    if use_broken_axis and ax_right is not None:
        fig.subplots_adjust(left=0.24, right=0.985, top=0.92, bottom=0.07, wspace=0.06)
    else:
        fig.subplots_adjust(left=0.24, right=0.985, top=0.93, bottom=0.07)

    out_png = output_dir / "gg2_top25_species_label_frequency.png"
    out_svg = output_dir / "gg2_top25_species_label_frequency.svg"
    fig.savefig(out_png, bbox_inches="tight", pad_inches=0.20)
    fig.savefig(out_svg, bbox_inches="tight", pad_inches=0.20)
    plt.close(fig)
    logger.info(f"Frequency + cluster summary saved to: {summary_tsv}")
    logger.info(f"Species-label frequency plot saved to: {out_png}")
    logger.info(f"Species-label frequency plot saved to: {out_svg}")
    return out_png


def plot_species_label_frequency_from_existing_run(
    output_dir: str | Path,
    header_list_path: Optional[str | Path] = None,
    top_n: int = 25,
    reuse_intermediates: bool = True,
    clust_ids_override: Optional[list[float]] = None,
    cluster_colors: Optional[dict] = None,
    unique_bar_color: str = "#CC79A7",
) -> Optional[Path]:
    """
    Recompute the taxonomy-frequency + V1-V9 clustering bar plot from an existing
    run directory without rerunning the full pipeline.
    """
    output_dir = Path(output_dir)
    stats_tsv = output_dir / "silva_amplicon_overall_stats.tsv"
    if not stats_tsv.exists():
        run_dirs = sorted(
            [p for p in output_dir.glob("start*_end*_step*") if (p / "silva_amplicon_overall_stats.tsv").exists()]
        )
        if len(run_dirs) == 1:
            output_dir = run_dirs[0]
            stats_tsv = output_dir / "silva_amplicon_overall_stats.tsv"
            logger.info(f"Using nested run directory for plotting: {output_dir}")
    _, clust_ids = _read_stats_rows_from_tsv(stats_tsv)
    if clust_ids_override:
        clust_ids = [float(x) for x in clust_ids_override]
    if not clust_ids:
        logger.warning(f"No clustering thresholds found in: {stats_tsv}")
    v1v9_nogap = output_dir / "silva_v1-v9_alignment_1250-1650_filtered_nogap.fasta"
    if header_list_path is None:
        header_list_path = Path(SILVA_ALGN_PATHS[4]).with_name(
            f"{Path(SILVA_ALGN_PATHS[4]).stem}_headers.list"
        )
    return plot_species_label_frequency(
        header_list_path=header_list_path,
        output_dir=output_dir,
        top_n=top_n,
        clust_ids=clust_ids,
        v1v9_nogap_fasta=v1v9_nogap,
        reuse_intermediates=reuse_intermediates,
        cluster_colors=cluster_colors,
        unique_bar_color=unique_bar_color,
    )


def _seqkit_stats_cached(
    fasta_file: Path,
    stats_path: Path,
    reuse: bool,
) -> tuple[int, dict]:
    stats_path.parent.mkdir(parents=True, exist_ok=True)
    if reuse and stats_path.exists():
        logger.info(f"Reusing existing stats: {stats_path}")
        return 0, _load_stats_tsv(stats_path)

    _r, _, stats_dict = seqkit_stats(fasta_file=fasta_file, output_file=stats_path)
    if _r != 0:
        logger.warning(f"Error getting stats for {fasta_file}")
    return _r, stats_dict


def _seqkit_rmdup_cached(
    input_fasta: Path,
    output_fasta: Path,
    reuse: bool,
) -> tuple[int, Path]:
    if reuse and output_fasta.exists():
        logger.info(f"Reusing existing deduplicated FASTA: {output_fasta}")
        return 0, output_fasta
    return seqkit_rmdup(
        input_fasta=input_fasta,
        output_fasta=output_fasta,
        args_list=["-s", "-j", "10"],
    )


def _seqkit_seq_cached(
    input_fasta: Path,
    output_fasta: Path,
    args_list: List[str],
    reuse: bool,
) -> tuple[int, Path]:
    if reuse and output_fasta.exists():
        logger.info(f"Reusing existing seqkit output: {output_fasta}")
        return 0, output_fasta
    return seqkit_seq(
        input_fasta=input_fasta,
        output_fasta=output_fasta,
        args_list=args_list,
    )


# 1. Cut SILVA alignemt at specific positions
def cut_silva_alignment(
    silva_algn_fasta: str,
    start_pos: int,
    end_pos: int,
    output_fasta: str
) -> None:
    """Cut SILVA alignment from start_pos to end_pos and save to output_fasta."""
    logger.info(f"Cutting SILVA alignment: {silva_algn_fasta} -> {output_fasta} (start: {start_pos}, end: {end_pos})")
    args_list = ["-r", f"{start_pos}:{end_pos}", "-j", "1"]

    ret_code, out_fasta = seqkit_subseq(
        input_fasta=Path(silva_algn_fasta),
        output_fasta=Path(output_fasta),
        args_list=args_list
    )
    if ret_code != 0:
        logger.error(f"Error cutting SILVA alignment: {silva_algn_fasta}")
        sys.exit(1)
    logger.info(f"Cut SILVA alignment saved to: {out_fasta}")

# 2. Remove sequences shorter and longer than specific lengths (e.g., 1200-1600 bp)
# The length shoule be calculated without gaps. Write the filtered sequences to a 
# new FASTA file
def filter_fasta_by_length(
    input_fasta: str,
    min_length: int,
    max_length: int,
    output_fasta: str,
    args_list: List[str] = []
) -> None:
    """Filter sequences in input_fasta by length (without gaps) and save to output_fasta."""
    logger.info(f"Filtering FASTA by length: {input_fasta} -> {output_fasta} (min: {min_length}, max: {max_length})")
    args_list = [f"-m", f"{min_length}", f"-M", f"{max_length}"] + args_list
    # NOTE: --remove-gaps is required to calculate length without gaps
    ret_code, out_fasta = seqkit_seq(
        input_fasta=Path(input_fasta),
        output_fasta=Path(output_fasta),
        args_list=args_list
    )
    if ret_code != 0:
        logger.error(f"Error filtering FASTA by length: {input_fasta}")
        sys.exit(1)
    logger.info(f"Filtered FASTA saved to: {out_fasta}")
 
# Define the worker function for parallel execution
def process_amplicon(length, start_pos, silva_algn_fasta, output_dir, reuse_intermediates: bool):
    try:
        logger.info(f"Processing amplicon length: {length}")
        cut_end_pos = start_pos + length - 1
        outdir_cut = output_dir / f"amplicon_{start_pos}_{cut_end_pos}"
        outdir_cut.mkdir(parents=True, exist_ok=True)
        cut_fasta = outdir_cut / f"silva_align_amplicon_{start_pos}_{cut_end_pos}.fasta"
        nogap_fasta = outdir_cut / f"silva_align_amplicon_{start_pos}_{cut_end_pos}_nogap.fasta"
        unique_fasta = outdir_cut / f"silva_align_amplicon_{start_pos}_{cut_end_pos}_unique.fasta"
        unique_nogap_fasta = outdir_cut / f"silva_align_amplicon_{start_pos}_{cut_end_pos}_nogap_unique.fasta"
        ret_codes = []
        # Cut the alignment
        if not (reuse_intermediates and cut_fasta.exists()):
            cut_silva_alignment(
                silva_algn_fasta=silva_algn_fasta,
                start_pos=start_pos,
                end_pos=cut_end_pos,
                output_fasta=str(cut_fasta)
            )
        else:
            logger.info(f"Reusing existing cut FASTA: {cut_fasta}")

        # Aligned (gapped) cut stats
        _r, cut_stats_dict = _seqkit_stats_cached(
            fasta_file=cut_fasta,
            stats_path=outdir_cut / f"silva_align_amplicon_{start_pos}_{cut_end_pos}_stats.tsv",
            reuse=reuse_intermediates,
        )
        ret_codes.append(_r)
        # Generate a unique sequences FASTA file
        ret_code, _ = _seqkit_rmdup_cached(
            input_fasta=cut_fasta,
            output_fasta=unique_fasta,
            reuse=reuse_intermediates,
        )
        if ret_code != 0:
            logger.warning(f"Error generating unique sequences: {cut_fasta}")
            # sys.exit(1)
        ret_codes.append(ret_code)
        # Count unique sequences
        _r, unique_stats_dict = _seqkit_stats_cached(
            fasta_file=unique_fasta,
            stats_path=outdir_cut / f"silva_align_amplicon_{start_pos}_{cut_end_pos}_unique_stats.tsv",
            reuse=reuse_intermediates,
        )
        if _r != 0:
            logger.warning(f"Error getting unique stats: {unique_fasta}")
        ret_codes.append(_r)
        # Remove gaps from the cut alignment
        args_list = [
            "--remove-gaps",
            "-j", "10",
            "-w", "0"
        ]
        ret_code, _ = _seqkit_seq_cached(
            input_fasta=cut_fasta,
            output_fasta=nogap_fasta,
            args_list=args_list,
            reuse=reuse_intermediates,
        )
        if ret_code != 0:
            logger.warning(f"Error removing gaps: {cut_fasta}")
            # sys.exit(1)
        ret_codes.append(ret_code)
        # Stats of nogap fasta
        _r, stats_dict = _seqkit_stats_cached(
            fasta_file=nogap_fasta,
            stats_path=outdir_cut / f"silva_align_amplicon_{start_pos}_{cut_end_pos}_nogap_stats.tsv",
            reuse=reuse_intermediates,
        )
        if _r != 0:
            logger.warning(f"Error getting nogap stats: {nogap_fasta}")
        # Count unique sequences after collapsing gaps
        ret_code, _ = _seqkit_rmdup_cached(
            input_fasta=nogap_fasta,
            output_fasta=unique_nogap_fasta,
            reuse=reuse_intermediates,
        )
        if ret_code != 0:
            logger.warning(f"Error generating unique nogap sequences: {nogap_fasta}")
        ret_codes.append(ret_code)
        _r, unique_nogap_stats_dict = _seqkit_stats_cached(
            fasta_file=unique_nogap_fasta,
            stats_path=outdir_cut / f"silva_align_amplicon_{start_pos}_{cut_end_pos}_unique_nogap_stats.tsv",
            reuse=reuse_intermediates,
        )
        if _r != 0:
            logger.warning(f"Error getting unique nogap stats: {unique_nogap_fasta}")
        ret_codes.append(_r)
        # Update stats dictionary with unique sequence count
        stats_dict["n_unique_align_seq"] = unique_stats_dict.get("num_seqs", 0)
        stats_dict["n_unique_nogap_seq"] = unique_nogap_stats_dict.get("num_seqs", 0)
        stats_dict["n_silva_species"] = _count_unique_taxonomies(cut_fasta)
        stats_dict["fragment_length"] = int(length)
        stats_dict["target_fragment_length"] = _to_float(stats_dict.get("avg_len", 0.0), 0.0)
        stats_dict["input_num_seqs"] = _to_int(cut_stats_dict.get("num_seqs", 0), 0)
        stats_dict["input_min_len"] = _to_int(cut_stats_dict.get("min_len", 0), 0)
        stats_dict["input_max_len"] = _to_int(cut_stats_dict.get("max_len", 0), 0)
        stats_dict["input_avg_len"] = _to_float(cut_stats_dict.get("avg_len", 0.0), 0.0)
        stats_dict["fasta_file"] = str(cut_fasta.name)
        stats_dict["nogap_fasta"] = str(nogap_fasta)
        stats_dict["_sort_key"] = int(length)

        logger.info(f"Amplicon FASTA generated: {nogap_fasta}")

        return stats_dict
    except Exception as e:
        logger.error(f"Error processing amplicon length {length}: {str(e)}")
        return {}


def usearch_cluster(
    input_fasta: str,
    cluster_id: float,
) -> dict:
    """Count ZOTUs in input_fasta using USEARCH."""
    logger.info(f"Counting ZOTUs in FASTA: {input_fasta}")
    # usearch_uclust already performs a length sort; avoid double sorting to
    # prevent input/output collisions on sorted_length.fasta in parallel runs.
    uc_file, centroids_fasta = usearch_uclust(
        input_fasta=input_fasta,
        cluster_id=cluster_id
    )
    renamed_stats_dicts = {}
    # Count number of ZOTUs
    try:
        __, _, stats_dict = seqkit_stats(centroids_fasta)
        renamed_stats_dicts = {
            f"zotus_{key}": value for key, value in stats_dict.items()
        }
    except Exception as e:
        logger.warning(f"Error getting ZOTU stats via seqkit: {centroids_fasta}: {str(e)}")
        # Fallback when seqkit is unavailable: count centroid records directly.
        try:
            centroid_count = sum(1 for _h, _s in _iter_fasta_records(Path(centroids_fasta)))
            renamed_stats_dicts = {
                "zotus_num_seqs": centroid_count,
                "zotus_source": "fasta_header_count_fallback",
            }
        except Exception as fallback_error:
            logger.error(
                f"Fallback centroid counting also failed for {centroids_fasta}: {fallback_error}"
            )
    return renamed_stats_dicts


def _build_alignment_stats_row(
    aligned_fasta: Path,
    nogap_fasta: Path,
    output_dir: Path,
    prefix: str,
    fragment_length: int,
    sort_key: int,
    reuse: bool,
) -> dict:
    """Build the stats row for a given alignment/nogap pair."""
    unique_aligned = output_dir / f"{prefix}_unique_aligned.fasta"
    unique_nogap = output_dir / f"{prefix}_unique_nogap.fasta"

    aligned_stats_file = output_dir / f"{prefix}_aligned_stats.tsv"
    unique_aligned_stats_file = output_dir / f"{prefix}_unique_aligned_stats.tsv"
    nogap_stats_file = output_dir / f"{prefix}_nogap_stats.tsv"
    unique_nogap_stats_file = output_dir / f"{prefix}_unique_nogap_stats.tsv"

    _r, aligned_stats = _seqkit_stats_cached(
        fasta_file=aligned_fasta,
        stats_path=aligned_stats_file,
        reuse=reuse,
    )
    if _r != 0:
        logger.warning(f"Error getting aligned stats: {aligned_fasta}")

    ret_code, _ = _seqkit_rmdup_cached(
        input_fasta=aligned_fasta,
        output_fasta=unique_aligned,
        reuse=reuse,
    )
    if ret_code != 0:
        logger.warning(f"Error generating unique aligned FASTA: {aligned_fasta}")
    _r, unique_aligned_stats = _seqkit_stats_cached(
        fasta_file=unique_aligned,
        stats_path=unique_aligned_stats_file,
        reuse=reuse,
    )
    if _r != 0:
        logger.warning(f"Error getting unique aligned stats: {unique_aligned}")

    _r, nogap_stats = _seqkit_stats_cached(
        fasta_file=nogap_fasta,
        stats_path=nogap_stats_file,
        reuse=reuse,
    )
    if _r != 0:
        logger.warning(f"Error getting nogap stats: {nogap_fasta}")

    ret_code, _ = _seqkit_rmdup_cached(
        input_fasta=nogap_fasta,
        output_fasta=unique_nogap,
        reuse=reuse,
    )
    if ret_code != 0:
        logger.warning(f"Error generating unique nogap FASTA: {nogap_fasta}")
    _r, unique_nogap_stats = _seqkit_stats_cached(
        fasta_file=unique_nogap,
        stats_path=unique_nogap_stats_file,
        reuse=reuse,
    )
    if _r != 0:
        logger.warning(f"Error getting unique nogap stats: {unique_nogap}")

    return {
        "fasta_file": str(aligned_fasta.name),
        "nogap_fasta": str(nogap_fasta),
        "input_num_seqs": _to_int(aligned_stats.get("num_seqs", 0), 0),
        "input_min_len": _to_int(aligned_stats.get("min_len", 0), 0),
        "input_max_len": _to_int(aligned_stats.get("max_len", 0), 0),
        "input_avg_len": _to_float(aligned_stats.get("avg_len", 0.0), 0.0),
        "fragment_length": int(fragment_length),
        "target_fragment_length": _to_float(nogap_stats.get("avg_len", 0.0), 0.0),
        "n_silva_species": _count_unique_taxonomies(aligned_fasta),
        "n_unique_align_seq": _to_int(unique_aligned_stats.get("num_seqs", 0), 0),
        "n_unique_nogap_seq": _to_int(unique_nogap_stats.get("num_seqs", 0), 0),
        "_sort_key": int(sort_key),
    }


def _cluster_counts_for_row(nogap_fasta: str, clust_ids: list[float]) -> dict:
    """
    Compute cluster counts for all requested IDs on one FASTA.
    Done sequentially per FASTA to avoid file overwrite races in usearch_uclust outputs.
    """
    out: dict[str, int] = {}
    for cid in clust_ids:
        stats = usearch_cluster(nogap_fasta, cid)
        col = f"Number of Clusters at {cid}"
        out[col] = _to_int(stats.get("zotus_num_seqs", 0), 0)
    return out


def _stage_for_clustering(input_fasta: Path, output_dir: Path, label: str) -> Path:
    """
    Copy clustering input into output_dir to keep intermediate USEARCH files
    out of reference/source directories.
    """
    stage_dir = output_dir / "cluster_work" / label
    stage_dir.mkdir(parents=True, exist_ok=True)
    staged = stage_dir / input_fasta.name
    if input_fasta.resolve() != staged.resolve():
        shutil.copy2(input_fasta, staged)
    return staged


def _plot_cluster_features_vs_length(
    cut_rows: list[dict],
    full_row: dict,
    clust_ids: list[float],
    output_dir: Path,
    trimmed_row: Optional[dict] = None,
) -> Optional[Path]:
    """
    Plot feature counts against rounded average ungapped amplicon length.
    Solid lines: clusters per ID across window cuts.
    Dashed horizontal lines: expected full-length clusters per ID.
    Black horizontal line: expected species from input SILVA taxonomies.
    """
    if plt is None:
        logger.warning("matplotlib is unavailable; skipping feature plot generation.")
        return None
    if not cut_rows:
        logger.warning("No cut rows available; skipping feature plot generation.")
        return None

    # Okabe-Ito palette (color-blind-friendly), excluding black reserved for species guide line.
    palette = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442"]
    ordered_cut_rows = sorted(cut_rows, key=lambda r: int(r.get("_sort_key", 0)))
    x_vals = [int(round(_to_float(r.get("target_fragment_length", 0.0), 0.0))) for r in ordered_cut_rows]

    def _smooth_xy(x_series: list[int], y_series: list[int]) -> tuple[list[float], list[float]]:
        if np is None or len(x_series) <= 2:
            return x_series, y_series
        dense_x = np.linspace(x_series[0], x_series[-1], len(x_series) * 5)
        dense_y = np.interp(dense_x, x_series, y_series)
        smooth_window = min(7, len(dense_y))
        if smooth_window > 1:
            kernel = np.ones(smooth_window) / smooth_window
            pad_left = smooth_window // 2
            pad_right = smooth_window - 1 - pad_left
            padded = np.pad(dense_y, (pad_left, pad_right), mode="edge")
            smooth_y = np.convolve(padded, kernel, mode="valid")
        else:
            smooth_y = dense_y
        return dense_x, smooth_y

    fig, ax = plt.subplots(figsize=(11, 6), dpi=200)
    unique_nogap_vals = [_to_int(r.get("n_unique_nogap_seq", 0), 0) for r in ordered_cut_rows]
    red_x, red_y = _smooth_xy(x_vals, unique_nogap_vals)
    ax.plot(
        red_x,
        red_y,
        color="red",
        linewidth=2.2,
        linestyle="-",
        label="Unique sequences",
    )
    ax.scatter(
        x_vals,
        unique_nogap_vals,
        color="red",
        edgecolors="white",
        linewidths=0.5,
        s=28,
        zorder=3,
    )
    if trimmed_row is not None:
        ax.axhline(
            y=_to_int(trimmed_row.get("n_unique_nogap_seq", 0), 0),
            color="red",
            linestyle="--",
            linewidth=1.8,
            label="Unique full sequence (V1-V9)",
        )

    for idx, cid in enumerate(clust_ids):
        col = f"Number of Clusters at {cid}"
        color = palette[idx % len(palette)]
        y_vals = [_to_int(r.get(col, 0), 0) for r in ordered_cut_rows]
        plot_x, plot_y = _smooth_xy(x_vals, y_vals)
        ax.plot(
            plot_x,
            plot_y,
            color=color,
            linewidth=2.0,
            label=f"Sequence clusters at {cid}",
        )
        ax.scatter(
            x_vals,
            y_vals,
            color=color,
            edgecolors="white",
            linewidths=0.5,
            s=32,
            zorder=3,
        )
        ref_row = trimmed_row if trimmed_row is not None else full_row
        ref_label = "Full sequence (V1-V9) clusters" if trimmed_row is not None else "Full-length clusters"
        ax.axhline(
            y=_to_int(ref_row.get(col, 0), 0),
            color=color,
            linestyle="--",
            linewidth=1.8,
            label=f"{ref_label} at {cid}",
        )

    ax.axhline(
        y=_to_int(full_row.get("n_silva_species", 0), 0),
        color="black",
        linestyle="-",
        linewidth=2.0,
        label="Expected species",
    )

    ax.set_xlabel("Average sequence length (bp)", fontsize=20)
    ax.set_ylabel("Feature count", fontsize=20)
    ax.set_title("Feature Counts vs. Sequence Length", fontsize=24)
    ax.locator_params(axis="x", nbins=14)
    ax.locator_params(axis="y", nbins=12)
    ax.grid(True, linestyle=":", alpha=0.35)
    if AutoMinorLocator is not None:
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.grid(True, which="minor", linestyle=":", alpha=0.2)
    ax.legend(loc="best", frameon=False, fontsize=12)
    fig.tight_layout()

    output_plot = output_dir / "silva_cluster_features_vs_length.png"
    output_plot_svg = output_dir / "silva_cluster_features_vs_length.svg"
    fig.savefig(output_plot)
    fig.savefig(output_plot_svg)
    plt.close(fig)
    logger.info(f"Feature plot saved to: {output_plot}")
    logger.info(f"Feature plot saved to: {output_plot_svg}")
    return output_plot


def _read_stats_rows_from_tsv(stats_tsv: Path) -> tuple[list[dict], list[float]]:
    stats_tsv = Path(stats_tsv)
    if not stats_tsv.exists():
        logger.warning(f"Stats TSV not found: {stats_tsv}")
        return [], []
    with stats_tsv.open() as fh:
        header_line = fh.readline().strip()
        if not header_line:
            return [], []
        headers = header_line.split("\t")
        cluster_ids = []
        for h in headers:
            if h.startswith("Number of Clusters at "):
                cid_str = h.replace("Number of Clusters at ", "")
                try:
                    cluster_ids.append(float(cid_str))
                except ValueError:
                    logger.warning(f"Unable to parse cluster column '{h}'")
        rows = []
        for line in fh:
            parts = line.strip().split("\t")
            if not parts:
                continue
            row = {}
            for idx, header in enumerate(headers):
                row[header] = parts[idx] if idx < len(parts) else ""
            rows.append(row)
    return rows, cluster_ids


def _tsv_row_to_internal(row: dict) -> dict:
    """Map one silva_amplicon_overall_stats.tsv row to internal plotting keys."""
    out = {
        "fasta_file": row.get("FASTA File", ""),
        "fragment_length": _to_int(row.get("Fragment Length", 0), 0),
        "target_fragment_length": _to_float(row.get("Target Fragment Length", 0.0), 0.0),
        "n_silva_species": _to_int(row.get("Number of SILVA Species", 0), 0),
        "n_unique_nogap_seq": _to_int(row.get("Number of Unique Sequences", 0), 0),
    }
    out["_sort_key"] = out["fragment_length"]
    for k, v in row.items():
        if k.startswith("Number of Clusters at "):
            out[k] = _to_int(v, 0)
    return out


def plot_from_stats_tsv(
    stats_tsv: Path,
    output_dir: Optional[Path] = None,
) -> Optional[Path]:
    """
    Re-create the cluster/features figure directly from a saved TSV.
    """
    stats_tsv = Path(stats_tsv)
    output_dir = Path(output_dir or stats_tsv.parent)
    rows, clust_ids = _read_stats_rows_from_tsv(stats_tsv)
    if not rows or not clust_ids:
        logger.warning("Not enough data in stats TSV to draw the plot.")
        return None
    rows_internal = [_tsv_row_to_internal(r) for r in rows]
    full_row = max(rows_internal, key=lambda r: _to_int(r.get("fragment_length", 0), 0))
    trimmed_row = next(
        (
            r for r in rows_internal
            if str(r.get("fasta_file", "")).startswith("silva_v1-v9_alignment_1250-1650_filtered")
        ),
        None,
    )
    cut_rows = [
        r for r in rows_internal
        if r is not full_row
        and r is not trimmed_row
        and str(r.get("fasta_file", "")).startswith("silva_align_amplicon_")
    ]
    if not cut_rows:
        cut_rows = [r for r in rows_internal if r is not full_row and r is not trimmed_row]
    return _plot_cluster_features_vs_length(
        cut_rows=cut_rows,
        full_row=full_row,
        clust_ids=clust_ids,
        output_dir=output_dir,
        trimmed_row=trimmed_row,
    )

# 3. Get amplicons of various length increased by 1000 bp by each step from cut SILVA
# alignment and store them in a directory
def generate_amplicon_fastas(
    silva_algn_fasta: str,
    start_pos: int,
    end_pos: int,
    step_size: int,
    min_amplicon_length: int,
    max_amplicon_length: int,
    output_dir: str,
    clust_id: float = 0.97,
    clust_ids: Optional[list[float]] = None,
    reuse_intermediates: bool = False,
) -> None:
    """
    Generate amplicon FASTA files of various lengths from cut SILVA alignment and write
    one summary row per FASTA (full-length first, then incremental cuts).
    `reuse_intermediates` instructs the pipeline to skip re-running operations when the
    expected intermediate outputs are already present.
    """
    total_cores = multiprocessing.cpu_count()
    num_workers = max(1, int(total_cores * 0.8))
    clust_ids_norm = _normalize_cluster_ids(clust_id=clust_id, clust_ids=clust_ids)

    today_date = datetime.datetime.now().strftime("%Y%m%d")
    clust_tag = "-".join(_cluster_id_label(cid) for cid in clust_ids_norm)
    output_dir = Path(f"{output_dir}{today_date}_clust{clust_tag}")
    output_dir = output_dir / f"start{start_pos}_end{end_pos}_step{step_size}"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Full-length row: ensure and validate *_nogap.fasta beside original input.
    silva_input_source = Path(silva_algn_fasta).resolve()
    silva_algn_input = silva_input_source
    silva_algn_input_nogap = _ensure_verified_gapless(silva_algn_input)
    silva_filtered_aligned = _filter_fasta_by_valid_species(
        input_fasta=silva_algn_input,
        output_fasta=output_dir / "silva_full_length_valid_taxonomy_aligned.fasta",
        reuse=reuse_intermediates,
        discarded_file=output_dir / "silva_full_length_discarded_invalid_taxonomy_aligned.tsv",
    )
    silva_filtered_nogap = _filter_fasta_by_valid_species(
        input_fasta=silva_algn_input_nogap,
        output_fasta=output_dir / "silva_full_length_valid_taxonomy_nogap.fasta",
        reuse=reuse_intermediates,
        discarded_file=output_dir / "silva_full_length_discarded_invalid_taxonomy_nogap.tsv",
    )
    silva_algn_input = silva_filtered_aligned
    silva_algn_input_nogap = silva_filtered_nogap
    full_row = _build_alignment_stats_row(
        aligned_fasta=silva_algn_input,
        nogap_fasta=silva_algn_input_nogap,
        output_dir=output_dir,
        prefix="silva_full_length",
        fragment_length=50000,
        sort_key=-1,
        reuse=reuse_intermediates,
    )
    # Requirement: run full-length clustering on verified nogap SILVA before any cutting.
    full_cluster_input = _stage_for_clustering(
        input_fasta=Path(str(full_row["nogap_fasta"])),
        output_dir=output_dir,
        label="full_length",
    )
    full_row.update(_cluster_counts_for_row(str(full_cluster_input), clust_ids_norm))

    logger.info(f"Generating V1-V9 alignment: {output_dir / 'silva_v1-v9_alignment.fasta'}")
    silva_alg_cut = output_dir / "silva_v1-v9_alignment.fasta"
    _ = seqkit_subseq(
        input_fasta=silva_algn_input,
        output_fasta=silva_alg_cut,
        args_list=["-r", f"{start_pos}:{end_pos}"],
        output_md5sum="",
    )

    filtered_nogap_silva_algn = output_dir / "silva_v1-v9_alignment_1250-1650_filtered_nogap.fasta"
    if reuse_intermediates and filtered_nogap_silva_algn.exists():
        logger.info(f"Reusing filtered nogap alignment: {filtered_nogap_silva_algn}")
    else:
        filter_fasta_by_length(
            input_fasta=silva_alg_cut,
            min_length=min_amplicon_length,
            max_length=max_amplicon_length,
            output_fasta=filtered_nogap_silva_algn,
            args_list=["--remove-gaps"],
        )

    silva_algn_len_filtered_fasta = output_dir / "silva_v1-v9_alignment_1250-1650_filtered.fasta"
    if reuse_intermediates and silva_algn_len_filtered_fasta.exists():
        logger.info(f"Reusing filtered aligned FASTA: {silva_algn_len_filtered_fasta}")
    else:
        seqkit_common(
            input_fasta1=silva_alg_cut,
            input_fasta2=filtered_nogap_silva_algn,
            output_fasta=silva_algn_len_filtered_fasta,
            args_list=["-n", "-i", "-j", str(num_workers)],
        )
    silva_algn_fasta = silva_algn_len_filtered_fasta
    logger.info(f"Using filtered SILVA alignment for amplicon generation: {silva_algn_fasta}")

    trimmed_fragment_length = end_pos - start_pos + 1

    trimmed_row = _build_alignment_stats_row(
        aligned_fasta=silva_algn_len_filtered_fasta,
        nogap_fasta=filtered_nogap_silva_algn,
        output_dir=output_dir,
        prefix="silva_v1-v9_filtered",
        fragment_length=trimmed_fragment_length,
        sort_key=0,
        reuse=reuse_intermediates,
    )
    trimmed_cluster_input = _stage_for_clustering(
        input_fasta=Path(trimmed_row["nogap_fasta"]),
        output_dir=output_dir,
        label="v1-v9_full",
    )
    trimmed_row.update(
        _cluster_counts_for_row(str(trimmed_cluster_input), clust_ids_norm)
    )

    lengths = list(range(step_size, end_pos - start_pos + 1, step_size))
    logger.info(f"Starting parallel amplicon processing with {num_workers} workers...")

    cut_rows: list[dict] = []
    try:
        with multiprocessing.Pool(processes=num_workers) as pool:
            this_args = [(length, 1, silva_algn_fasta, output_dir, reuse_intermediates) for length in lengths]
            results = pool.starmap(process_amplicon, this_args)
        cut_rows = [r for r in results if r]
    except Exception as e:
        logger.error(f"Error in parallel processing: {str(e)}")
        raise

    if not cut_rows:
        logger.warning("No cut amplicon rows were produced; writing full-length row only.")
    # Compute clustering counts for cut rows (sequential per FASTA to avoid overwrite races).
    try:
        if cut_rows:
            with multiprocessing.Pool(processes=num_workers) as pool:
                cluster_args = [(str(r["nogap_fasta"]), clust_ids_norm) for r in cut_rows]
                cluster_results = pool.starmap(_cluster_counts_for_row, cluster_args)
            for i, cluster_counts in enumerate(cluster_results):
                cut_rows[i].update(cluster_counts)
    except Exception as e:
        logger.error(f"Error in USEARCH clustering: {str(e)}")
        raise
    stats_rows = [full_row, trimmed_row] + cut_rows

    # Final table schema and ordering
    cluster_cols = [f"Number of Clusters at {cid}" for cid in clust_ids_norm]
    ordered_rows = sorted(stats_rows, key=lambda r: int(r.get("_sort_key", 0)))
    output_rows = []
    for r in ordered_rows:
        output_rows.append(
            {
                "FASTA File": r.get("fasta_file", ""),
                "Number of Sequences": r.get("input_num_seqs", 0),
                "Minimum Length": r.get("input_min_len", 0),
                "Maximum Length": r.get("input_max_len", 0),
                "Average Length": r.get("input_avg_len", 0.0),
                "Fragment Length": r.get("fragment_length", 0),
                "Target Fragment Length": r.get("target_fragment_length", 0.0),
                "Number of SILVA Species": r.get("n_silva_species", 0),
                "Number of Unique Aligned Sequences": r.get("n_unique_align_seq", 0),
                "Number of Unique Sequences": r.get("n_unique_nogap_seq", 0),
                **{c: r.get(c, 0) for c in cluster_cols},
            }
        )

    stats_tsv = output_dir / "silva_amplicon_overall_stats.tsv"
    headers = [
        "FASTA File",
        "Number of Sequences",
        "Minimum Length",
        "Maximum Length",
        "Average Length",
        "Fragment Length",
        "Target Fragment Length",
        "Number of SILVA Species",
        "Number of Unique Aligned Sequences",
        "Number of Unique Sequences",
        *cluster_cols,
    ]
    with open(stats_tsv, "w") as stats_file:
        stats_file.write("\t".join(headers) + "\n")
        for row in output_rows:
            stats_file.write("\t".join(str(row[h]) for h in headers) + "\n")
    logger.info(f"Overall statistics saved to: {stats_tsv}")
    _ = _plot_cluster_features_vs_length(
        cut_rows=cut_rows,
        full_row=full_row,
        clust_ids=clust_ids_norm,
        output_dir=output_dir,
        trimmed_row=trimmed_row,
    )
    header_list_candidate = silva_input_source.with_name(f"{silva_input_source.stem}_headers.list")
    header_list_for_plot = header_list_candidate if header_list_candidate.exists() else None
    _ = plot_species_label_frequency(
        header_list_path=header_list_for_plot,
        output_dir=output_dir,
        top_n=25,
        clust_ids=clust_ids_norm,
        v1v9_nogap_fasta=filtered_nogap_silva_algn,
        reuse_intermediates=reuse_intermediates,
    )
