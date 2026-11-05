import gzip
import subprocess
from pathlib import Path

import pytest

from tactic_pipeline import processing_helper
from tactic_pipeline import vsearch_commands as vs

DATA_DIR = Path(__file__).parent / "data" / "truncated"
SPIKES_FASTA = Path(__file__).resolve().parents[2] / "spikesidx" / "CFM_spikes.fasta"
VSEARCH = processing_helper.default_vsearch_bin()

needs_vsearch = pytest.mark.skipif(not VSEARCH, reason="vsearch binary not found")


def run(cmds, cwd):
    for cmd in cmds:
        result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
        assert result.returncode == 0, f"{' '.join(cmd)}\n{result.stderr}"


def head_fastq(gz_path: Path, out_path: Path, n_reads: int) -> Path:
    with gzip.open(gz_path, "rt") as in_h, open(out_path, "w") as out_h:
        for i, line in enumerate(in_h):
            if i >= 4 * n_reads:
                break
            out_h.write(line)
    return out_path


def fasta_labels(path: Path) -> list:
    return [line[1:].strip() for line in path.read_text().splitlines() if line.startswith(">")]


@needs_vsearch
def test_vsearch_commands_run_through_the_pipeline_steps(tmp_path):
    """Every vsearch command TACTIC builds is accepted by vsearch and gives the expected output."""
    samples = {
        "S1": "truncSRR13005876_S1_L001",
        "S2": "truncSRR13005987_S2_L001",
    }
    relabeled = []
    for sample_id, prefix in samples.items():
        sdir = tmp_path / sample_id
        sdir.mkdir()
        fwd = head_fastq(DATA_DIR / f"{prefix}_R1_001.fastq.gz", sdir / "R1.fastq", 1500)
        rev = head_fastq(DATA_DIR / f"{prefix}_R2_001.fastq.gz", sdir / "R2.fastq", 1500)
        # Per-sample processing
        run(vs.merge_pairs(VSEARCH, fwd, rev, "merged.fastq", 50, 50, 200, 600), sdir)
        run(vs.strip_ends(VSEARCH, "merged.fastq", "filtered1.fastq", 17, 21), sdir)
        run(vs.quality_filter(
            VSEARCH, "filtered1.fastq", "filtered2.fasta", 0.002, minlen=200
        ), sdir)
        assert fasta_labels(sdir / "filtered2.fasta")
        run(vs.dereplicate(VSEARCH, "filtered2.fasta", "derep.fasta"), sdir)
        run(vs.sort_by_size(VSEARCH, "derep.fasta", "sorted.fasta"), sdir)
        assert ";size=" in fasta_labels(sdir / "sorted.fasta")[0]
        # Analysis: label reads with their sample
        run(vs.relabel_with_sample(VSEARCH, "sorted.fasta", sample_id, "relabeled.fasta"), sdir)
        assert fasta_labels(sdir / "relabeled.fasta")[0].startswith(f"sample={sample_id};1;size=")
        relabeled.append(sdir / "relabeled.fasta")

    shared = tmp_path / "Shared"
    shared.mkdir()
    gathered = shared / "sorted.fasta"
    gathered.write_text("".join(p.read_text() for p in relabeled))
    run(vs.dereplicate(VSEARCH, gathered, "derep.fasta", relabel="Uniq"), shared)
    assert fasta_labels(shared / "derep.fasta")[0].startswith("Uniq1;size=")
    run(vs.sort_by_size(VSEARCH, "derep.fasta", "sorted_size.fasta"), shared)

    # zOTUs: named like USEARCH names them (Zotu1, no size annotation)
    run(vs.denoise(VSEARCH, "sorted_size.fasta", "zotus.fasta", 2, "unoise.tmp", threads=2), shared)
    zotu_labels = fasta_labels(shared / "zotus.fasta")
    assert zotu_labels and zotu_labels[0] == "Zotu1"
    run(vs.otu_table(VSEARCH, gathered, "zotus.fasta", "ZOTUs-Table.tab", 0.99, threads=2), shared)
    header = (shared / "ZOTUs-Table.tab").read_text().splitlines()[0].split("\t")
    assert header == ["#OTU ID", "S1", "S2"]

    # De novo OTUs
    run(vs.cluster_otus(
        VSEARCH, "sorted_size.fasta", "otus1.fa", 8, "abund.tmp", "cent.tmp", threads=2
    ), shared)
    assert fasta_labels(shared / "otus1.fa")[0] == "OTU1"

    # TAC clustering keeps labels that contain spaces
    taxed = shared / "taxed_zotus.fasta"
    lines = (shared / "zotus.fasta").read_text().splitlines()
    taxed.write_text("\n".join(
        f"{line} tax=Bacteria;Bacillota;Incertae Sedis;" if line.startswith(">") else line
        for line in lines
    ) + "\n")
    run(vs.sort_by_length(VSEARCH, taxed, "sorted_length.fasta"), shared)
    run(vs.cluster_smallmem(
        VSEARCH, "sorted_length.fasta", "centroids.fasta", "clusters.uc", 0.987
    ), shared)
    assert fasta_labels(shared / "centroids.fasta")[0].endswith("Incertae Sedis;")
    uc_lines = (shared / "clusters.uc").read_text().splitlines()
    assert any(line.startswith("S\t") for line in uc_lines)

    run(vs.make_udb(VSEARCH, "zotus.fasta", "zotus.udb"), shared)
    assert (shared / "zotus.udb").stat().st_size > 0


@needs_vsearch
@pytest.mark.skipif(not SPIKES_FASTA.is_file(), reason="spike references not found")
def test_vsearch_spike_search_writes_hits_and_unmatched_reads(tmp_path):
    reads = tmp_path / "filtered2.fasta"
    spike_seq = "".join(
        line.strip() for line in SPIKES_FASTA.read_text().splitlines()[1:40]
        if not line.startswith(">")
    )[:400]
    reads.write_text(f">spike_read\n{spike_seq}\n>other_read\n{'ACGT' * 100}\n")
    run(vs.search_spikes(VSEARCH, reads, SPIKES_FASTA, 0.9, "hits.tsv", "nospike.fasta"), tmp_path)
    hits = (tmp_path / "hits.tsv").read_text().splitlines()
    assert [h.split("\t")[0] for h in hits] == ["spike_read"]
    assert fasta_labels(tmp_path / "nospike.fasta") == ["other_read"]


def test_merge_pairs_converts_pctid_to_maxdiffpct():
    cmd = vs.merge_pairs("vsearch", "f.fq", "r.fq", "m.fq", 50, 90, 200, 600)[0]
    assert cmd[cmd.index("--fastq_maxdiffpct") + 1] == "10"
    assert "--fastq_pctid" not in cmd


def test_merge_otutab_columns_sums_samples(tmp_path):
    table = tmp_path / "zotu_table.txt"
    table.write_text("#OTU ID\tM03455\tM0345\nZotu1\t3\t2\nZotu2\t0\t4\n")
    vs.merge_otutab_columns(table)
    assert table.read_text() == "#OTU ID\tM03455\nZotu1\t5\nZotu2\t4\n"


def test_merge_otutab_columns_keeps_single_column(tmp_path):
    table = tmp_path / "zotu_table.txt"
    table.write_text("#OTU ID\tS1\nZotu1\t3\n")
    vs.merge_otutab_columns(table)
    assert table.read_text() == "#OTU ID\tS1\nZotu1\t3\n"


def test_select_fasta_records(tmp_path):
    fasta = tmp_path / "in.fa"
    fasta.write_text(">Zotu1\nACGT\n>Zotu2\nGGGG\nCCCC\n>Zotu10\nTTTT\n")
    labels = tmp_path / "labels.txt"
    labels.write_text("Zotu2\nZotu10\n")
    out = tmp_path / "out.fa"
    assert vs.select_fasta_records(fasta, labels, out) == 2
    assert out.read_text() == ">Zotu2\nGGGG\nCCCC\n>Zotu10\nTTTT\n"


def test_sequence_tool_defaults_to_usearch_and_is_part_of_preprocessing_args(tmp_path):
    default_args = processing_helper.TACTICArgsParser(config_yaml="").preproc_args
    assert processing_helper.seq_tool_name(default_args) == "usearch"

    yml = tmp_path / "args.yml"
    yml.write_text("SequenceToolArgs:\n  tool: VSearch\n")
    vsearch_args = processing_helper.TACTICArgsParser(config_yaml=str(yml)).preproc_args
    assert processing_helper.seq_tool_name(vsearch_args) == "vsearch"
    # A different tool makes the preprocessing arguments different, so samples get reprocessed
    assert vsearch_args != default_args

    # Argument files written before this option existed count as usearch
    old_yml = tmp_path / "old_args.yml"
    old_yml.write_text("MergePairsArgs:\n  fastq_maxdiffs: 50\n")
    old_args = processing_helper.TACTICArgsParser(config_yaml=str(old_yml)).preproc_args
    assert old_args == default_args


@needs_vsearch
def test_vsearch_binary_check():
    assert processing_helper.Vsearch().check_or_get_bin() == (0, VSEARCH)


def test_vsearch_binary_check_missing_file(tmp_path):
    assert processing_helper.Vsearch(str(tmp_path / "no_vsearch")).check_or_get_bin() == (161, "")
