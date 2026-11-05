"""
vsearch commands for the TACTIC pipeline steps.

TACTIC runs USEARCH 11 by default. With

    SequenceToolArgs:
      tool: vsearch

in the arguments YAML file, the steps below run vsearch instead. Each function returns the
command(s) of one step as argument lists for ``system_sub``. vsearch stops with an error when a
command gets an option it does not support, so every command gets only its own options.

Differences from USEARCH 11.0.667 (checked with vsearch 2.29.3 on the test data):
- Read merging: ``-fastq_pctid`` has no vsearch counterpart; ``--fastq_maxdiffpct`` is set to
  100 - pctid. vsearch rejects many more pairs whose reverse read is of poor quality ("alignment
  score too low, or score drop too high"; there is no option for this). Reads left after quality
  filtering, vsearch vs USEARCH: 29,476 vs 32,395 (-9 %) for test sample truncSRR13005876, but
  76 vs 217 (-65 %) for the first 1,000 pairs of truncSRR13005987, whose reverse reads are poor.
- zOTUs: USEARCH ``-unoise3`` is two steps in vsearch (``--cluster_unoise``, then
  ``--uchime3_denovo``). On the test data both gave the same zOTUs and the same zOTU table.
- De novo OTUs: vsearch has no UPARSE (``-cluster_otus``). It is replaced by clustering at 97 %
  (``--cluster_size``) followed by de novo chimera removal (``--uchime_denovo``), so the OTUs
  differ slightly.
- TAC clustering (``--cluster_smallmem``) and the ``.uc`` output are the same as with USEARCH.
"""
from pathlib import Path

USEARCH = "usearch"
VSEARCH = "vsearch"
SUPPORTED_SEQ_TOOLS = (USEARCH, VSEARCH)

# vsearch refuses FASTQ quality values above 41 by default; newer Illumina data can exceed that.
FASTQ_QMAX = "93"
# vsearch cuts sequence labels at the first space; USEARCH keeps them. TACTIC headers can carry the
# taxonomy after a space (e.g. "Zotu5 tax=...;Incertae Sedis;"), so steps that see such headers keep
# the full label.
KEEP_FULL_LABELS = "--notrunclabels"
# The spike references are long (several kb). Checking only the best few candidate spikes per
# read halves the run time without missing spike reads (same hits as the vsearch default on the
# test data).
SPIKE_SEARCH_MAXREJECTS = "4"


def tool_description(vsearch_bin) -> str:
    """'vsearch (<path>)' when vsearch is used, otherwise 'usearch'. For the logs."""
    return f"{VSEARCH} ({vsearch_bin})" if vsearch_bin else USEARCH


def _fmt(value) -> str:
    """Formats numbers without a trailing '.0' (50.0 -> '50', 0.987 -> '0.987')."""
    return f"{value:g}" if isinstance(value, float) else str(value)


def merge_pairs(vsearch, forward_fastq, reverse_fastq, out_fastq,
                maxdiffs, pctid, minmergelen, maxmergelen, threads=1) -> list:
    return [[
        str(vsearch), "--fastq_mergepairs", str(forward_fastq),
        "--reverse", str(reverse_fastq),
        "--fastq_maxdiffs", _fmt(maxdiffs),
        "--fastq_maxdiffpct", _fmt(100 - float(pctid)),
        "--fastq_minmergelen", _fmt(minmergelen),
        "--fastq_maxmergelen", _fmt(maxmergelen),
        # USEARCH merges staggered pairs (reverse read reaching past the forward read) by default
        "--fastq_allowmergestagger",
        "--fastq_qmax", FASTQ_QMAX,
        "--threads", _fmt(threads),
        "--fastqout", str(out_fastq),
    ]]


def strip_ends(vsearch, in_fastq, out_fastq, stripleft, stripright=None) -> list:
    """Primer removal by length (USEARCH -fastx_truncate -stripleft/-stripright)."""
    cmd = [str(vsearch), "--fastx_filter", str(in_fastq), "--fastq_stripleft", _fmt(stripleft)]
    if stripright is not None:
        cmd += ["--fastq_stripright", _fmt(stripright)]
    cmd += ["--fastq_qmax", FASTQ_QMAX, "--fastqout", str(out_fastq)]
    return [cmd]


def quality_filter(vsearch, in_fastq, out_fasta, maxee_rate, minlen=None, trunclen=None) -> list:
    cmd = [str(vsearch), "--fastq_filter", str(in_fastq), "--fastq_maxee_rate", _fmt(maxee_rate)]
    if minlen is not None:
        cmd += ["--fastq_minlen", _fmt(minlen)]
    if trunclen is not None:
        cmd += ["--fastq_trunclen", _fmt(trunclen)]
    cmd += ["--fastq_qmax", FASTQ_QMAX, "--fastaout", str(out_fasta)]
    return [cmd]


def search_spikes(vsearch, query_fasta, spikes_fasta, identity, userout, notmatched,
                  threads=1) -> list:
    return [[
        str(vsearch), "--usearch_global", str(query_fasta),
        "--db", str(spikes_fasta),
        "--id", _fmt(identity),
        "--strand", "both",
        "--maxrejects", SPIKE_SEARCH_MAXREJECTS,
        "--userout", str(userout),
        "--userfields", "query+target",
        "--notmatched", str(notmatched),
        "--threads", _fmt(threads),
    ]]


def dereplicate(vsearch, in_fasta, out_fasta, relabel=None) -> list:
    cmd = [
        str(vsearch), "--fastx_uniques", str(in_fasta),
        "--sizein", "--sizeout", "--strand", "both",
    ]
    if relabel:
        cmd += ["--relabel", str(relabel)]
    cmd += ["--fastaout", str(out_fasta)]
    return [cmd]


def sort_by_size(vsearch, in_fasta, out_fasta) -> list:
    return [[
        str(vsearch), "--sortbysize", str(in_fasta), KEEP_FULL_LABELS, "--output", str(out_fasta)
    ]]


def sort_by_length(vsearch, in_fasta, out_fasta) -> list:
    return [[
        str(vsearch), "--sortbylength", str(in_fasta), KEEP_FULL_LABELS, "--output", str(out_fasta)
    ]]


def relabel_with_sample(vsearch, in_fastx, sample_id, out_fasta, out_fastq=None) -> list:
    """
    Labels reads as 'sample=<sample_id>;<n>;size=<s>' (USEARCH -fastx_relabel -prefix
    'sample=<id>;' -keep_annots), so that --otutabout puts them into the column of their sample.
    """
    cmd = [
        str(vsearch), "--fastx_filter", str(in_fastx),
        "--relabel", f"sample={sample_id};",
        "--sizein", "--sizeout",
        "--fastaout", str(out_fasta),
    ]
    if out_fastq:
        cmd += ["--fastq_qmax", FASTQ_QMAX, "--fastqout", str(out_fastq)]
    return [cmd]


def denoise(vsearch, sorted_fasta, zotus_out, minsize, centroids_tmp, threads=1) -> list:
    """
    USEARCH -unoise3: denoising, then chimera removal. zOTUs are named Zotu1, Zotu2, ... without
    size annotation, as USEARCH names them.
    """
    return [
        [
            str(vsearch), "--cluster_unoise", str(sorted_fasta),
            "--minsize", _fmt(minsize),
            "--sizein", "--sizeout",
            "--strand", "both",
            "--threads", _fmt(threads),
            "--centroids", str(centroids_tmp),
        ],
        [
            str(vsearch), "--uchime3_denovo", str(centroids_tmp),
            "--sizein",
            "--relabel", "Zotu",
            "--nonchimeras", str(zotus_out),
        ],
    ]


def cluster_otus(vsearch, sorted_fasta, otus_out, minsize, abundant_tmp, centroids_tmp,
                 identity=0.97, threads=1) -> list:
    """
    Replacement for USEARCH -cluster_otus (UPARSE): keep sequences seen at least `minsize` times,
    cluster them at 97 %, then remove chimeric OTUs. OTUs are named OTU1, OTU2, ...
    """
    return [
        [
            str(vsearch), "--sortbysize", str(sorted_fasta),
            "--minsize", _fmt(minsize),
            "--sizein", "--sizeout",
            "--output", str(abundant_tmp),
        ],
        [
            str(vsearch), "--cluster_size", str(abundant_tmp),
            "--id", _fmt(identity),
            "--sizein", "--sizeout",
            "--strand", "both",
            "--threads", _fmt(threads),
            "--centroids", str(centroids_tmp),
        ],
        [
            str(vsearch), "--uchime_denovo", str(centroids_tmp),
            "--sizein",
            "--relabel", "OTU",
            "--nonchimeras", str(otus_out),
        ],
    ]


def otu_table(vsearch, reads_fasta, db_fasta, out_table, identity, threads=1) -> list:
    """USEARCH -otutab: maps reads (with size annotations) to the OTUs/zOTUs in db_fasta."""
    return [[
        str(vsearch), "--usearch_global", str(reads_fasta),
        "--db", str(db_fasta),
        "--id", _fmt(identity),
        "--sizein",
        "--strand", "both",
        "--threads", _fmt(threads),
        "--otutabout", str(out_table),
    ]]


def cluster_smallmem(vsearch, sorted_by_length_fasta, centroids, uc_file, identity) -> list:
    """TAC clustering. The input must be sorted by decreasing length."""
    return [[
        str(vsearch), "--cluster_smallmem", str(sorted_by_length_fasta),
        "--id", _fmt(identity),
        "--strand", "both",
        KEEP_FULL_LABELS,
        "--centroids", str(centroids),
        "--uc", str(uc_file),
    ]]


def make_udb(vsearch, in_fasta, out_udb) -> list:
    return [[str(vsearch), "--makeudb_usearch", str(in_fasta), "--output", str(out_udb)]]


def merge_otutab_columns(table_path) -> None:
    """
    Sums all sample columns of an OTU table into one column. Used for the per-sample zOTU table:
    vsearch takes the sample name from the read labels, which can split one sample into several
    columns. The column is named after the first sample column.
    """
    table_path = Path(table_path)
    lines = table_path.read_text().splitlines()
    if not lines:
        return
    header = lines[0].split("\t")
    if len(header) <= 2:
        return
    out = [f"{header[0]}\t{header[1]}"]
    for line in lines[1:]:
        if not line.strip():
            continue
        fields = line.split("\t")
        total = sum(float(v) for v in fields[1:])
        out.append(f"{fields[0]}\t{total:g}")
    table_path.write_text("\n".join(out) + "\n")


def select_fasta_records(in_fasta, labels_file, out_fasta) -> int:
    """
    Writes the records of in_fasta whose label is listed in labels_file (one label per line), like
    USEARCH -fastx_getseqs -labels. Returns the number of records written.
    """
    wanted = {line.strip() for line in Path(labels_file).read_text().splitlines() if line.strip()}
    written = 0
    keep = False
    with open(in_fasta) as in_h, open(out_fasta, "w") as out_h:
        for line in in_h:
            if line.startswith(">"):
                keep = line[1:].strip().split()[0] in wanted if line[1:].strip() else False
                written += keep
            if keep:
                out_h.write(line)
    return written
