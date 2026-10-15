import gzip
import logging
import os
import subprocess
from pathlib import Path

import pytest

from tactic_pipeline import processing_helper


def test_gunzip_arb_file_streams_to_an_atomic_destination(tmp_path, monkeypatch):
    content = (b"ACGT" * (256 * 1024)) + b"end"
    gzip_path = tmp_path / "SILVA.arb.gz"
    with gzip.open(gzip_path, "wb") as compressed_file:
        compressed_file.write(content)

    real_gzip_open = processing_helper.gzip.open
    read_sizes = []

    class TrackingReader:
        def __init__(self, wrapped):
            self.wrapped = wrapped

        def __enter__(self):
            self.wrapped.__enter__()
            return self

        def __exit__(self, *args):
            return self.wrapped.__exit__(*args)

        def read(self, size=-1):
            read_sizes.append(size)
            assert size > 0, "the complete decompressed database was read at once"
            return self.wrapped.read(size)

    def tracking_gzip_open(*args, **kwargs):
        return TrackingReader(real_gzip_open(*args, **kwargs))

    monkeypatch.setattr(processing_helper.gzip, "open", tracking_gzip_open)

    arb_path = Path(processing_helper.gunzip_arb_file(
        str(gzip_path),
        keep=True,
        chunk_size=64 * 1024
    ))

    assert arb_path.read_bytes() == content
    assert gzip_path.exists()
    assert len(read_sizes) > 1
    assert not list(tmp_path.glob(".SILVA.arb.*.tmp"))


def test_index_silva_database_uses_one_thread_and_reuses_index(
        tmp_path,
        monkeypatch):
    gzip_path = tmp_path / "SILVA.arb.gz"
    with gzip.open(gzip_path, "wb") as compressed_file:
        compressed_file.write(b"test ARB database")

    calls = []

    def fake_system_sub(command, **kwargs):
        calls.append(command)
        arb_path = Path(command[command.index("--db") + 1])
        arb_path.with_suffix(".sidx").write_bytes(b"test SINA index")
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(processing_helper, "system_sub", fake_system_sub)

    index_path = Path(processing_helper.index_silva_database(
        str(gzip_path),
        sina_bin="sina-test",
        logger_obj=logging.getLogger("test.silva")
    ))

    assert calls == [[
        "sina-test",
        "--db", str(tmp_path / "SILVA.arb"),
        "--in", str(tmp_path / "fake1.fasta"),
        "--out", str(tmp_path / "fake2.fasta"),
        "--threads", "1"
    ]]
    assert index_path.read_bytes() == b"test SINA index"
    assert index_path.with_suffix(".sidx.success").exists()
    assert not (tmp_path / "fake1.fasta").exists()
    assert not (tmp_path / "fake2.fasta").exists()

    # A second run must neither re-extract the ARB nor rebuild its fresh index.
    arb_mtime = (tmp_path / "SILVA.arb").stat().st_mtime_ns
    assert processing_helper.index_silva_database(
        str(gzip_path),
        sina_bin="sina-test",
        logger_obj=logging.getLogger("test.silva")
    ) == str(index_path)
    assert (tmp_path / "SILVA.arb").stat().st_mtime_ns == arb_mtime
    assert len(calls) == 1


def test_index_silva_database_removes_partial_index_after_oom(
        tmp_path,
        monkeypatch):
    arb_path = tmp_path / "SILVA.arb"
    arb_path.write_bytes(b"test ARB database")
    index_path = arb_path.with_suffix(".sidx")
    success_path = index_path.with_suffix(".sidx.success")
    success_path.write_text("stale marker", encoding="utf-8")

    def fail_with_partial_index(command, **kwargs):
        index_path.write_bytes(b"partial")
        raise MemoryError("out of memory")

    monkeypatch.setattr(
        processing_helper,
        "system_sub",
        fail_with_partial_index
    )

    with pytest.raises(MemoryError, match="out of memory"):
        processing_helper.index_silva_database(
            str(arb_path),
            logger_obj=logging.getLogger("test.silva")
        )

    assert not index_path.exists()
    assert not success_path.exists()
    assert not (tmp_path / "fake1.fasta").exists()
    assert not (tmp_path / "fake2.fasta").exists()


def test_find_files_and_dirs_owned_by_root(tmp_path, monkeypatch):
    sub_dir = tmp_path / "databases"
    sub_dir.mkdir()
    file1 = sub_dir / "SILVA.arb"
    file1.write_bytes(b"data")
    file2 = sub_dir / "user_file.txt"
    file2.write_bytes(b"user")

    real_stat = Path.stat

    def fake_stat(self, *args, **kwargs):
        res = real_stat(self, *args, **kwargs)
        uid = 0 if self.name in ("SILVA.arb", "databases") else 1000
        return os.stat_result((
            res.st_mode, res.st_ino, res.st_dev, res.st_nlink,
            uid, res.st_gid, res.st_size,
            int(res.st_atime), int(res.st_mtime), int(res.st_ctime)
        ))

    monkeypatch.setattr(Path, "stat", fake_stat)

    root_items = processing_helper.find_files_and_dirs_owned_by_root(sub_dir)
    assert str(sub_dir) in root_items
    assert str(file1) in root_items
    assert str(file2) not in root_items

    non_existent = tmp_path / "does_not_exist"
    assert processing_helper.find_files_and_dirs_owned_by_root(non_existent) == []


# NOTE: do not import run_tactic in unit tests. Importing it opens
# /base/inputs/Pipeline_log.txt in 'w' mode (truncating it) and registers
# signal/atexit handlers. Keep testable logic in processing_helper instead.
def test_correct_files_modes(tmp_path, monkeypatch):
    target_dir = tmp_path / "test_db"
    target_dir.mkdir()
    test_file = target_dir / "test.txt"
    test_file.write_text("hello")

    chmod_calls = []

    def fake_chmod(self, mode):
        chmod_calls.append((str(self), mode))

    monkeypatch.setattr(Path, "chmod", fake_chmod)
    monkeypatch.setattr(
        processing_helper,
        "find_files_and_dirs_owned_by_root",
        lambda d: [str(target_dir), str(test_file)] if str(d) == str(target_dir) else []
    )

    # duplicated and missing directories must not cause repeated or failing chmods
    mode_result = processing_helper.correct_files_modes(
        [target_dir, str(target_dir), tmp_path / "missing", "", None],
        mode=0o777
    )
    assert mode_result == 0o777
    assert sorted(chmod_calls) == sorted([(str(target_dir), 0o777), (str(test_file), 0o777)])


def test_pair_seq_files_deterministic_sorting(tmp_path):
    dir_z = tmp_path / "z_folder"
    dir_a = tmp_path / "a_folder"
    dir_m = tmp_path / "m_folder"
    for d in [dir_z, dir_a, dir_m]:
        d.mkdir()
        (d / "sample_R1.fastq").write_text("@r1\nACGT\n+\nIIII\n")
        (d / "sample_R2.fastq").write_text("@r2\nTGCA\n+\nIIII\n")

    pairs = processing_helper.pair_seq_files(str(tmp_path))
    assert len(pairs) == 3
    expected_order = [
        (str(dir_a / "sample_R1.fastq"), str(dir_a / "sample_R2.fastq")),
        (str(dir_m / "sample_R1.fastq"), str(dir_m / "sample_R2.fastq")),
        (str(dir_z / "sample_R1.fastq"), str(dir_z / "sample_R2.fastq")),
    ]
    assert pairs == expected_order


def test_is_fastq_extension():
    assert processing_helper.is_fastq_extension("sample.fastq") is True
    assert processing_helper.is_fastq_extension("sample.fq") is True
    assert processing_helper.is_fastq_extension("sample_R1_001.fastq.gz") is True
    assert processing_helper.is_fastq_extension("SAMPLE_R2.FQ.GZ") is True
    assert processing_helper.is_fastq_extension("sample.fasta") is False
    assert processing_helper.is_fastq_extension("report.txt") is False
    assert processing_helper.is_fastq_extension("archive.tar.gz") is False


def test_pair_seq_files_skips_non_fastq_and_processed(tmp_path, monkeypatch):
    run_dir = tmp_path / "run_dir"
    run_dir.mkdir()
    proc_dir = run_dir / "sample_01__processed"
    proc_dir.mkdir()

    # Valid FASTQ pair
    (run_dir / "sample_01_R1.fastq").write_text("@seq1\nACGT\n+\nIIII\n")
    (run_dir / "sample_01_R2.fastq").write_text("@seq1\nTGCA\n+\nIIII\n")

    # Non-FASTQ files that should be completely skipped without opening
    (run_dir / "notes.txt").write_text("not a fastq")
    (run_dir / "image.png").write_text("binary data")
    (proc_dir / "processed_R1.fastq").write_text("@proc1\nACGT\n+\nIIII\n")
    (proc_dir / "processed_R2.fastq").write_text("@proc1\nTGCA\n+\nIIII\n")

    opened_files = []
    real_is_seq_file = processing_helper.is_seq_file

    def tracking_is_seq_file(path, *args, **kwargs):
        opened_files.append(str(path))
        return real_is_seq_file(path, *args, **kwargs)

    monkeypatch.setattr(processing_helper, "is_seq_file", tracking_is_seq_file)

    pairs = processing_helper.pair_seq_files(str(run_dir))
    assert len(pairs) == 1
    assert pairs[0] == (str(run_dir / "sample_01_R1.fastq"), str(run_dir / "sample_01_R2.fastq"))
    # notes.txt, image.png, and processed files must never be passed to is_seq_file
    assert str(run_dir / "notes.txt") not in opened_files
    assert str(run_dir / "image.png") not in opened_files
    assert str(proc_dir / "processed_R1.fastq") not in opened_files


def test_find_fastq_pair_for_sample_with_and_without_parent(tmp_path):
    folder_a = tmp_path / "folder_a"
    folder_b = tmp_path / "folder_b"
    folder_a.mkdir()
    folder_b.mkdir()

    (folder_a / "sample_A_R1_001.fastq.gz").write_bytes(gzip.compress(b"@sA\nACGT\n+\nIIII\n"))
    (folder_a / "sample_A_R2_001.fastq.gz").write_bytes(gzip.compress(b"@sA\nTGCA\n+\nIIII\n"))
    (folder_b / "sample_B_R1_001.fastq.gz").write_bytes(gzip.compress(b"@sB\nACGT\n+\nIIII\n"))
    (folder_b / "sample_B_R2_001.fastq.gz").write_bytes(gzip.compress(b"@sB\nTGCA\n+\nIIII\n"))

    # Test with explicit parent_path
    pair_a = processing_helper.find_fastq_pair_for_sample(
        sample_id="sample_A",
        parent_path="folder_a",
        base_fastq_dir=tmp_path
    )
    assert pair_a == (
        str(folder_a / "sample_A_R1_001.fastq.gz"),
        str(folder_a / "sample_A_R2_001.fastq.gz")
    )

    # Test without parent_path (subfolder recursive discovery)
    pair_b = processing_helper.find_fastq_pair_for_sample(
        sample_id="sample_B",
        parent_path="",
        base_fastq_dir=tmp_path
    )
    assert pair_b == (
        str(folder_b / "sample_B_R1_001.fastq.gz"),
        str(folder_b / "sample_B_R2_001.fastq.gz")
    )


def test_find_fastq_pair_exact_sample_id_differentiation(tmp_path):
    (tmp_path / "Sample1_R1_001.fastq").write_text("@s1\nACGT\n+\nIIII\n")
    (tmp_path / "Sample1_R2_001.fastq").write_text("@s1\nTGCA\n+\nIIII\n")
    (tmp_path / "Sample10_R1_001.fastq").write_text("@s10\nACGT\n+\nIIII\n")
    (tmp_path / "Sample10_R2_001.fastq").write_text("@s10\nTGCA\n+\nIIII\n")

    # Sample1 must match Sample1 and not Sample10
    pair_1 = processing_helper.find_fastq_pair_for_sample(
        sample_id="Sample1",
        base_fastq_dir=tmp_path
    )
    assert pair_1 == (
        str(tmp_path / "Sample1_R1_001.fastq"),
        str(tmp_path / "Sample1_R2_001.fastq")
    )

    # Sample10 must match Sample10
    pair_10 = processing_helper.find_fastq_pair_for_sample(
        sample_id="Sample10",
        base_fastq_dir=tmp_path
    )
    assert pair_10 == (
        str(tmp_path / "Sample10_R1_001.fastq"),
        str(tmp_path / "Sample10_R2_001.fastq")
    )


def test_find_fastq_pair_claimed_files(tmp_path):
    (tmp_path / "sample_dup_R1.fastq").write_text("@s\nACGT\n+\nIIII\n")
    (tmp_path / "sample_dup_R2.fastq").write_text("@s\nTGCA\n+\nIIII\n")

    claimed = set()
    pair_1 = processing_helper.find_fastq_pair_for_sample(
        sample_id="sample_dup",
        base_fastq_dir=tmp_path,
        claimed_files=claimed
    )
    assert pair_1[0] == str(tmp_path / "sample_dup_R1.fastq")
    assert pair_1[1] == str(tmp_path / "sample_dup_R2.fastq")
    claimed.add(pair_1[0])
    claimed.add(pair_1[1])

    # Second call with claimed files must return empty
    pair_2 = processing_helper.find_fastq_pair_for_sample(
        sample_id="sample_dup",
        base_fastq_dir=tmp_path,
        claimed_files=claimed
    )
    assert pair_2 == ("", "")


def test_parse_mapping_file_targeted_without_seq_file_pairs(tmp_path, monkeypatch):
    import run_tactic

    folder = tmp_path / "fastqs"
    folder.mkdir()
    (folder / "SampleA_R1_001.fastq").write_text("@sA\nACGT\n+\nIIII\n")
    (folder / "SampleA_R2_001.fastq").write_text("@sA\nTGCA\n+\nIIII\n")
    (folder / "SampleB_R1_001.fastq").write_text("@sB\nACGT\n+\nIIII\n")
    (folder / "SampleB_R2_001.fastq").write_text("@sB\nTGCA\n+\nIIII\n")

    map_content = "#SampleID\ttotal_weight_in_g\tspike_amount\tparent_path\nSampleA\t1.5\t5.0\tfastqs\n"
    map_file = tmp_path / "mapping_file.csv"
    map_file.write_text(map_content)

    monkeypatch.setattr(run_tactic, "FASTQ_DIR", tmp_path)

    # Calling parse_mapping_file with ONLY the mapping_file_path (no files_tups)
    mapping_dict = run_tactic.parse_mapping_file(mapping_file_path=str(map_file))

    assert "SampleA" in mapping_dict
    map_entry = mapping_dict["SampleA"]
    map_line, file_pair = map_entry
    assert map_line.SampleID == "SampleA"
    assert map_line.total_weight_in_g == 1.5
    assert map_line.spike_amount == "5.0"
    assert file_pair == (
        str(folder / "SampleA_R1_001.fastq"),
        str(folder / "SampleA_R2_001.fastq")
    )


def test_find_fastq_pair_warns_when_reverse_pair_missing(tmp_path, caplog):
    (tmp_path / "sample_unpaired_R1_001.fastq").write_text("@s1\nACGT\n+\nIIII\n")
    logger = logging.getLogger("test_unpaired")
    logger.setLevel(logging.WARNING)

    with caplog.at_level(logging.WARNING):
        pair = processing_helper.find_fastq_pair_for_sample(
            sample_id="sample_unpaired",
            base_fastq_dir=tmp_path,
            logger_obj=logger
        )

    assert pair == (str(tmp_path / "sample_unpaired_R1_001.fastq"), "")
    assert "Could not find matching reverse FASTQ pair for forward file 'sample_unpaired_R1_001.fastq'" in caplog.text
    assert "Treating as single-end." in caplog.text


def test_find_fastq_pair_warns_when_parent_path_missing(tmp_path, caplog):
    logger = logging.getLogger("test_missing_dir")
    logger.setLevel(logging.WARNING)

    with caplog.at_level(logging.WARNING):
        pair = processing_helper.find_fastq_pair_for_sample(
            sample_id="sample_missing",
            parent_path="non_existent_folder",
            base_fastq_dir=tmp_path,
            logger_obj=logger
        )

    assert pair == ("", "")
    assert "Search directory for parent_path 'non_existent_folder' does not exist" in caplog.text


def test_find_fastq_pair_warns_when_only_reverse_read_present(tmp_path, caplog):
    (tmp_path / "sample_only_rev_R2_001.fastq").write_text("@s2\nTGCA\n+\nIIII\n")
    logger = logging.getLogger("test_only_rev")
    logger.setLevel(logging.WARNING)

    with caplog.at_level(logging.WARNING):
        pair = processing_helper.find_fastq_pair_for_sample(
            sample_id="sample_only_rev",
            base_fastq_dir=tmp_path,
            logger_obj=logger
        )

    assert pair == ("", "")
    assert "Found reverse FASTQ file(s) for sample 'sample_only_rev'" in caplog.text
    assert "A forward read is required." in caplog.text


def test_tac_and_tic_default_args():
    tic_defaults = processing_helper.ComplexTICArgs.default_args
    tac_defaults = processing_helper.TACClusterArgs.default_args

    assert tic_defaults["abund_limit"] == 0.0
    assert tic_defaults["sample_wise_correction"] is False
    assert tac_defaults["abund_limit"] == 0.0
    assert tac_defaults["sample_wise_correction"] is False


def test_tactic_args_parser_reflects_yaml_and_defaults_when_omitted(tmp_path):
    # 1. Parse official TACTICPipeline_args.yml
    official_yaml = Path("/src/IMNGS2Pipeline/TACTICPipeline_args.yml")
    parser = processing_helper.TACTICArgsParser(str(official_yaml))

    assert parser.analysis_args.complex_tic.abund_limit == 0.0
    assert parser.analysis_args.complex_tic.sample_wise_correction is False
    assert parser.analysis_args.tac_cluster.abund_limit == 0.0
    assert parser.analysis_args.tac_cluster.sample_wise_correction is False

    # 2. Parse YAML where ComplexTICArgs and TACClusterArgs are completely omitted
    sparse_yaml = tmp_path / "sparse.yml"
    sparse_yaml.write_text("BuildZOTUTableArgs:\n  match_id: 0.99\n")
    sparse_parser = processing_helper.TACTICArgsParser(str(sparse_yaml))

    assert sparse_parser.analysis_args.complex_tic.abund_limit == 0.0
    assert sparse_parser.analysis_args.complex_tic.sample_wise_correction is False
    assert sparse_parser.analysis_args.tac_cluster.abund_limit == 0.0
    assert sparse_parser.analysis_args.tac_cluster.sample_wise_correction is False

    # 3. Parse YAML where sections exist but abund_limit / sample_wise_correction are omitted
    partial_yaml = tmp_path / "partial.yml"
    partial_yaml.write_text("ComplexTICArgs:\n  family_sim: 0.85\nTACClusterArgs:\n  cluster_thr: 0.97\n")
    partial_parser = processing_helper.TACTICArgsParser(str(partial_yaml))

    assert partial_parser.analysis_args.complex_tic.family_sim == 0.85
    assert partial_parser.analysis_args.complex_tic.abund_limit == 0.0
    assert partial_parser.analysis_args.complex_tic.sample_wise_correction is False
    assert partial_parser.analysis_args.tac_cluster.cluster_thr == 0.97
    assert partial_parser.analysis_args.tac_cluster.abund_limit == 0.0
    assert partial_parser.analysis_args.tac_cluster.sample_wise_correction is False


def test_args_parser_util_numeric_coercion(tmp_path):
    # Test that integer 0 in YAML properly coerces to float 0.0, and explicit floats/bools work
    test_yaml = tmp_path / "types.yml"
    test_yaml.write_text(
        "ComplexTICArgs:\n"
        "  abund_limit: 0\n"
        "  sample_wise_correction: true\n"
        "TACClusterArgs:\n"
        "  abund_limit: 0.005\n"
        "  sample_wise_correction: false\n"
    )
    parser = processing_helper.TACTICArgsParser(str(test_yaml))

    assert isinstance(parser.analysis_args.complex_tic.abund_limit, float)
    assert parser.analysis_args.complex_tic.abund_limit == 0.0
    assert parser.analysis_args.complex_tic.sample_wise_correction is True

    assert isinstance(parser.analysis_args.tac_cluster.abund_limit, float)
    assert parser.analysis_args.tac_cluster.abund_limit == 0.005
    assert parser.analysis_args.tac_cluster.sample_wise_correction is False


def test_get_base_name_various_formats():
    # 1. Hyphenated R1 / R2 (user screenshot format)
    assert processing_helper.get_base_name("T01-Ad4-V2-R1.fastq") == "T01-Ad4-V2"
    assert processing_helper.get_base_name("T01-Ad4-V2-R2.fastq") == "T01-Ad4-V2"

    # 2. Standard Illumina format with Lane and Sample numbers
    raw_name = processing_helper.get_base_name("537-522-3751-FU-1_S388_L001_R1_001.fastq.gz")
    assert raw_name == "537-522-3751-FU-1_S388_L001"
    clean_name = processing_helper.get_base_name("537-522-3751-FU-1_S388_L001_R1_001.fastq.gz", strip_illumina=True)
    assert clean_name == "537-522-3751-FU-1"

    # 3. Case insensitivity (lowercase r1 / r2)
    assert processing_helper.get_base_name("sample_r1_001.fastq.gz") == "sample"
    assert processing_helper.get_base_name("sample_r2_001.fastq.gz") == "sample"

    # 4. Word indicators: forward / reverse, read1 / read2
    assert processing_helper.get_base_name("sample_forward.fastq") == "sample"
    assert processing_helper.get_base_name("sample_reverse.fastq") == "sample"
    assert processing_helper.get_base_name("sample_read1.fq.gz") == "sample"
    assert processing_helper.get_base_name("sample_read2.fq.gz") == "sample"

    # 5. Numbered 1 / 2 indicators
    assert processing_helper.get_base_name("sample_1.fastq.gz") == "sample"
    assert processing_helper.get_base_name("sample_2.fastq.gz") == "sample"
    assert processing_helper.get_base_name("sample.1.fq") == "sample"
    assert processing_helper.get_base_name("sample-1.fq.gz") == "sample"

    # 6. Rightmost match preservation (embedded R1 / 1 in sample name)
    assert processing_helper.get_base_name("Patient_R1_rep1_R1_001.fastq.gz") == "Patient_R1_rep1"
    assert processing_helper.get_base_name("Cohort_1_pool_1.fastq.gz") == "Cohort_1_pool"

    # 7. Single-end reads without forward / reverse indicators
    assert processing_helper.get_base_name("SingleSample.fastq.gz") == "SingleSample"
    assert processing_helper.get_base_name("my_reads.fq") == "my_reads"


def test_find_reverse_file_hyphenated_and_various_formats(tmp_path):
    # 1. Hyphenated R1 / R2
    fwd_hyphen = tmp_path / "T01-Ad4-V2-R1.fastq"
    rev_hyphen = tmp_path / "T01-Ad4-V2-R2.fastq"
    fwd_hyphen.write_text("@s\nACGT\n+\nIIII\n")
    rev_hyphen.write_text("@s\nTGCA\n+\nIIII\n")
    assert processing_helper.find_reverse_file(fwd_hyphen) == rev_hyphen

    # 2. Standard Illumina
    fwd_illum = tmp_path / "537-522-3751-FU-1_S388_L001_R1_001.fastq.gz"
    rev_illum = tmp_path / "537-522-3751-FU-1_S388_L001_R2_001.fastq.gz"
    fwd_illum.write_text("dummy")
    rev_illum.write_text("dummy")
    assert processing_helper.find_reverse_file(fwd_illum) == rev_illum

    # 3. Lowercase r1 / r2
    fwd_lower = tmp_path / "sample_r1_001.fastq.gz"
    rev_lower = tmp_path / "sample_r2_001.fastq.gz"
    fwd_lower.write_text("dummy")
    rev_lower.write_text("dummy")
    assert processing_helper.find_reverse_file(fwd_lower) == rev_lower

    # 4. Dot-separated 1 / 2
    fwd_dot = tmp_path / "sample.1.fq"
    rev_dot = tmp_path / "sample.2.fq"
    fwd_dot.write_text("dummy")
    rev_dot.write_text("dummy")
    assert processing_helper.find_reverse_file(fwd_dot) == rev_dot


def test_pair_seq_files_hyphenated_and_single_end(tmp_path):
    run_dir = tmp_path / "reads"
    run_dir.mkdir()

    # Hyphenated pair
    fwd_h = run_dir / "T01-Ad4-V2-R1.fastq"
    rev_h = run_dir / "T01-Ad4-V2-R2.fastq"
    fwd_h.write_text("@s1\nACGT\n+\nIIII\n")
    rev_h.write_text("@s1\nTGCA\n+\nIIII\n")

    # Single-end read
    se = run_dir / "Sample_SE.fastq"
    se.write_text("@s2\nAAAA\n+\nIIII\n")

    pairs = processing_helper.pair_seq_files(str(run_dir))
    assert len(pairs) == 2
    assert (str(fwd_h), str(rev_h)) in pairs
    assert (str(se), "") in pairs


def test_find_fastq_pair_for_sample_illumina_clean_id(tmp_path):
    # Setup files: FU-1 and FU-10 to verify no prefix collision occurs
    (tmp_path / "537-522-3751-FU-1_S388_L001_R1_001.fastq").write_text("@fu1\nACGT\n+\nIIII\n")
    (tmp_path / "537-522-3751-FU-1_S388_L001_R2_001.fastq").write_text("@fu1\nTGCA\n+\nIIII\n")
    (tmp_path / "537-522-3751-FU-10_S389_L001_R1_001.fastq").write_text("@fu10\nACGT\n+\nIIII\n")
    (tmp_path / "537-522-3751-FU-10_S389_L001_R2_001.fastq").write_text("@fu10\nTGCA\n+\nIIII\n")

    # Sample ID in mapping file is clean without Illumina lane suffixes
    pair_fu1 = processing_helper.find_fastq_pair_for_sample(
        sample_id="537-522-3751-FU-1",
        base_fastq_dir=tmp_path
    )
    assert pair_fu1 == (
        str(tmp_path / "537-522-3751-FU-1_S388_L001_R1_001.fastq"),
        str(tmp_path / "537-522-3751-FU-1_S388_L001_R2_001.fastq")
    )

    pair_fu10 = processing_helper.find_fastq_pair_for_sample(
        sample_id="537-522-3751-FU-10",
        base_fastq_dir=tmp_path
    )
    assert pair_fu10 == (
        str(tmp_path / "537-522-3751-FU-10_S389_L001_R1_001.fastq"),
        str(tmp_path / "537-522-3751-FU-10_S389_L001_R2_001.fastq")
    )


def test_find_fastq_pair_for_sample_hyphenated(tmp_path):
    (tmp_path / "T01-Ad4-V2-R1.fastq").write_text("@t1\nACGT\n+\nIIII\n")
    (tmp_path / "T01-Ad4-V2-R2.fastq").write_text("@t1\nTGCA\n+\nIIII\n")

    pair = processing_helper.find_fastq_pair_for_sample(
        sample_id="T01-Ad4-V2",
        base_fastq_dir=tmp_path
    )
    assert pair == (
        str(tmp_path / "T01-Ad4-V2-R1.fastq"),
        str(tmp_path / "T01-Ad4-V2-R2.fastq")
    )


def test_file_discovery_skips_analysis_directories(tmp_path):
    # 1. Normal sample directory with valid FASTQ files
    valid_dir = tmp_path / "fastqs"
    valid_dir.mkdir()
    (valid_dir / "sample_A_R1.fastq").write_text("@sA\nACGT\n+\nIIII\n")
    (valid_dir / "sample_A_R2.fastq").write_text("@sA\nTGCA\n+\nIIII\n")

    # Sample whose filename happens to start with Analysis_
    (valid_dir / "Analysis_named_sample_R1.fastq").write_text("@sAn\nACGT\n+\nIIII\n")
    (valid_dir / "Analysis_named_sample_R2.fastq").write_text("@sAn\nTGCA\n+\nIIII\n")

    # 2. Output directory starting with Analysis_
    analysis_dir = tmp_path / "Analysis_20261015_120000" / "Shared" / "relabeling"
    analysis_dir.mkdir(parents=True)
    (analysis_dir / "sample_A_relabeled.fastq").write_text("@sA\nACGT\n+\nIIII\n")
    (analysis_dir / "sample_old_R1.fastq").write_text("@old\nACGT\n+\nIIII\n")
    (analysis_dir / "sample_old_R2.fastq").write_text("@old\nTGCA\n+\nIIII\n")

    # 3. Processed directory
    proc_dir = tmp_path / "sample_A__processed" / "20261015"
    proc_dir.mkdir(parents=True)
    (proc_dir / "sample_A_R1_trimmed.fastq").write_text("@sA\nACGT\n+\nIIII\n")

    # Verify pair_seq_files skips both Analysis_ and __processed directories
    pairs = processing_helper.pair_seq_files(str(tmp_path))
    found_forws = [p[0] for p in pairs]
    assert str(valid_dir / "sample_A_R1.fastq") in found_forws
    assert str(valid_dir / "Analysis_named_sample_R1.fastq") in found_forws
    assert str(analysis_dir / "sample_old_R1.fastq") not in found_forws
    assert str(analysis_dir / "sample_A_relabeled.fastq") not in found_forws
    assert str(proc_dir / "sample_A_R1_trimmed.fastq") not in found_forws

    # Verify find_fastq_pair_for_sample skips Analysis_ directories
    # Sample "sample_old" exists ONLY inside Analysis_ directory, so it must return ("", "")
    pair_old = processing_helper.find_fastq_pair_for_sample(
        sample_id="sample_old",
        base_fastq_dir=tmp_path
    )
    assert pair_old == ("", "")

    # Sample "sample_A" must match the valid directory files, not the Analysis_ or __processed files
    pair_a = processing_helper.find_fastq_pair_for_sample(
        sample_id="sample_A",
        base_fastq_dir=tmp_path
    )
    assert pair_a == (
        str(valid_dir / "sample_A_R1.fastq"),
        str(valid_dir / "sample_A_R2.fastq")
    )






