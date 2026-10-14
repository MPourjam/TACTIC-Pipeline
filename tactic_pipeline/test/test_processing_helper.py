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


def test_correct_created_files_modes(tmp_path, monkeypatch):
    import run_tactic
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

    mode_result = run_tactic.correct_created_files_modes(mode=777, target_dirs=[target_dir])
    assert mode_result == 0o777
    assert (str(target_dir), 0o777) in chmod_calls
    assert (str(test_file), 0o777) in chmod_calls
