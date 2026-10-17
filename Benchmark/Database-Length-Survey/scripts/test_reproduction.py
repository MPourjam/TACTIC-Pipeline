"""Small tests for length-statistic definitions and archive safety."""

import gzip
import hashlib
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

from survey_io import expected_matches, scan, summarize, write_tsv
from verify_input_archive import verify


class ReproductionTests(unittest.TestCase):
    def test_legacy_fasta_statistics(self):
        script = Path(__file__).resolve().parents[1] / "code/length_survey/compare_zotu_locations.py"
        spec = importlib.util.spec_from_file_location("legacy", script)
        legacy = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(legacy)
        lengths = [1, 2, 3, 3, 5, 9]
        actual = summarize({n: lengths.count(n) for n in set(lengths)})
        for key, value in legacy.calculate_stats(lengths).items():
            self.assertAlmostEqual(actual[key], value)

    def test_fastq_prefixes_and_empty_reads(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "reads.fastq.gz"
            with gzip.open(path, "wb") as handle:
                handle.write(b"@a\n\n+\n\n" + b"@b\nACG\n+\nIII\n" * 1001)
            result = scan(path, "fastq")
            self.assertEqual(result["stats"]["count"], 1002)
            self.assertEqual(result["samples"]["1000"]["count"], 1000)
            self.assertEqual(result["stats"]["min"], 0)

    def test_empty_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "empty.fastq"
            path.touch()
            self.assertEqual(scan(path, "fastq")["stats"]["count"], 0)

    def test_corrupt_fastq_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.fastq"
            path.write_bytes(b"@x\nAC\n+\nI\n")
            with self.assertRaises(ValueError):
                scan(path, "fastq")

    def test_saved_rounding(self):
        self.assertEqual(expected_matches({"mean": 1.234}, {"mean_2dp": "1.23"}), [])
        self.assertEqual(expected_matches({"mean": 1.234}, {"mean_2dp": "1.24"}), ["mean_2dp"])

    def test_zip_checksums(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.zip"
            data = b"input data"
            with zipfile.ZipFile(path, "w") as archive:
                archive.writestr("payloads/input", data)
                archive.writestr("input_manifest/SHA256SUMS", hashlib.sha256(data).hexdigest() + "  payloads/input\n")
            self.assertEqual(verify(path)["zip_crc_and_sha256"], "passed")
            with zipfile.ZipFile(path, "w") as archive:
                archive.writestr("payloads/input", b"changed")
                archive.writestr("input_manifest/SHA256SUMS", hashlib.sha256(data).hexdigest() + "  payloads/input\n")
            with self.assertRaises(ValueError):
                verify(path)

    def test_portable_runner_zip_and_extracted(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "input_manifest"
            data = gzip.compress(b"@x\nACG\n+\nIII\n")
            payload = root / "payloads/reads.fastq.gz"
            payload.parent.mkdir()
            payload.write_bytes(data)
            write_tsv(manifest / "inputs.tsv", [{"input_id": "test", "kind": "fastq",
                "data_path": "payloads/reads.fastq.gz", "sha256": hashlib.sha256(data).hexdigest()}])
            write_tsv(manifest / "observations.tsv", [{"observation_id": "o", "script": "example",
                "input_id": "test", "mock": "m", "max_records": "1000",
                "expected_json": json.dumps({"count": 1, "mean": 3})}])
            write_tsv(manifest / "findings.tsv", [], ["severity", "input", "detail"])
            archive_path = root / "data.zip"
            with zipfile.ZipFile(archive_path, "w") as archive:
                archive.write(payload, "payloads/reads.fastq.gz")
            runner = Path(__file__).parent / "run_length_survey.py"
            for name, extra in [("extracted", []), ("streaming", ["--archive", str(archive_path)])]:
                command = [sys.executable, str(runner), "--data-dir", str(root),
                           "--output-dir", str(root / name), *extra]
                subprocess.run(command, check=True, capture_output=True)
                summary = json.loads((root / name / "validation.json").read_text())
                self.assertEqual(summary["mismatches"], 0)
                self.assertEqual(summary["observations"], 1)


if __name__ == "__main__":
    unittest.main()
