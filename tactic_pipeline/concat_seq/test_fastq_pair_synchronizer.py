import unittest
import os
import tempfile
import shutil
import sys
from unittest.mock import patch
from io import StringIO

# Import the functions and main from the script
# Ensure fastq_pair_synchronizer.py is in the same directory or PYTHONPATH
try:
    from fastq_pair_synchronizer import clean_header, match_and_pair_fastq, index_fastq_file, main
except ImportError:
    # Fallback if running directly without the module installed/path set
    # This assumes the script is named exactly fastq_pair_synchronizer.py
    import importlib.util
    spec = importlib.util.spec_from_file_location("fastq_pair_synchronizer", "fastq_pair_synchronizer.py")
    if spec and spec.loader:
        fastq_pair_synchronizer = importlib.util.module_from_spec(spec)
        sys.modules["fastq_pair_synchronizer"] = fastq_pair_synchronizer
        spec.loader.exec_module(fastq_pair_synchronizer)
        from fastq_pair_synchronizer import clean_header, match_and_pair_fastq, index_fastq_file, main

class TestHeaderCleaning(unittest.TestCase):
    """
    Tests the clean_header function to ensure it correctly extracts
    unique read IDs from various FASTQ header formats.
    """

    def test_basic_header(self):
        self.assertEqual(clean_header("@SEQ_ID"), "SEQ_ID")
        self.assertEqual(clean_header("SEQ_ID"), "SEQ_ID")

    def test_casava_1_8_format(self):
        # Casava 1.8+ format: @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
        header = "@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG"
        expected = "EAS139:136:FC706VJ:2:2104:15343:197393"
        self.assertEqual(clean_header(header), expected)

    def test_old_illumina_suffix(self):
        # Older format with /1 or /2 suffix
        self.assertEqual(clean_header("@HWUSI-EAS100R:6:73:941:1973#0/1"), "HWUSI-EAS100R:6:73:941:1973#0")
        self.assertEqual(clean_header("@HWUSI-EAS100R:6:73:941:1973#0/2"), "HWUSI-EAS100R:6:73:941:1973#0")

    def test_whitespace_variations(self):
        # Test headers with spaces
        self.assertEqual(clean_header("@SEQ_ID 1:N:0:1"), "SEQ_ID")
        self.assertEqual(clean_header("@SEQ_ID 2:N:0:1"), "SEQ_ID")

class TestFastqSync(unittest.TestCase):
    """
    Integration tests for the file synchronization logic.
    Creates temporary FASTQ files to test pairing and orphan removal.
    """

    def setUp(self):
        # Create a temporary directory for test files
        self.test_dir = tempfile.mkdtemp()
        
        # Define file paths
        self.r1_path = os.path.join(self.test_dir, "test_R1.fastq")
        self.r2_path = os.path.join(self.test_dir, "test_R2.fastq")
        self.out_r1_path = os.path.join(self.test_dir, "out_R1.fastq")
        self.out_r2_path = os.path.join(self.test_dir, "out_R2.fastq")

    def tearDown(self):
        # Clean up the temporary directory after tests
        shutil.rmtree(self.test_dir)

    def create_fastq_file(self, path, records):
        """Helper to create a dummy FASTQ file from a list of tuples (header, seq)."""
        with open(path, 'w') as f:
            for header, seq in records:
                f.write(f"{header}\n{seq}\n+\nIIII\n")

    def test_match_and_pair_logic(self):
        """Test the core logic function directly."""
        # Scenario:
        # Read A: Present in both (Matched)
        # Read B: Present in both (Matched)
        # Read C: Only in R1 (Orphan R1 - should be discarded)
        # Read D: Only in R2 (Orphan R2 - should be discarded)
        
        r1_records = [
            ("@ReadA/1", "AAAA"),
            ("@ReadC/1", "CCCC"), # Orphan
            ("@ReadB/1", "BBBB"),
        ]
        
        r2_records = [
            ("@ReadD/2", "DDDD"), # Orphan
            ("@ReadB/2", "bbbb"),
            ("@ReadA/2", "aaaa"),
        ]

        self.create_fastq_file(self.r1_path, r1_records)
        self.create_fastq_file(self.r2_path, r2_records)

        # Run the synchronization
        # We capture stdout to avoid cluttering test output with progress bars
        with patch('sys.stdout', new=StringIO()):
            match_and_pair_fastq(
                self.r1_path, 
                self.r2_path, 
                self.out_r1_path, 
                self.out_r2_path
            )

        # Verify R1 Output (Expect ReadA and ReadB)
        with open(self.out_r1_path, 'r') as f:
            content_r1 = f.read()
            self.assertIn("@ReadA/1", content_r1)
            self.assertIn("@ReadB/1", content_r1)
            self.assertNotIn("@ReadC/1", content_r1) # Orphan C removed

        # Verify R2 Output (Expect ReadA and ReadB)
        with open(self.out_r2_path, 'r') as f:
            content_r2 = f.read()
            self.assertIn("@ReadA/2", content_r2)
            self.assertIn("@ReadB/2", content_r2)
            self.assertNotIn("@ReadD/2", content_r2) # Orphan D removed

    def test_cli_integration(self):
        """Test the main function (CLI) using mocked sys.argv."""
        r1_records = [("@ReadA/1", "AAAA")]
        r2_records = [("@ReadA/2", "aaaa")]
        
        self.create_fastq_file(self.r1_path, r1_records)
        self.create_fastq_file(self.r2_path, r2_records)

        # Mock command line arguments
        test_args = [
            "fastq_pair_synchronizer.py",
            "--r1", self.r1_path,
            "--r2", self.r2_path,
            "--out-r1", self.out_r1_path,
            "--out-r2", self.out_r2_path
        ]

        with patch.object(sys, 'argv', test_args):
            with patch('sys.stdout', new=StringIO()) as fake_out:
                main()
                
                # Check if output says something about success or summary
                output = fake_out.getvalue()
                self.assertIn("Synchronization Summary", output)

        # Verify output files were actually created
        self.assertTrue(os.path.exists(self.out_r1_path))
        self.assertTrue(os.path.exists(self.out_r2_path))

    def test_cli_missing_file_error(self):
        """Test that the CLI handles missing input files gracefully."""
        # Don't create the files
        test_args = [
            "fastq_pair_synchronizer.py",
            "--r1", "missing_R1.fastq",
            "--r2", "missing_R2.fastq",
            "--out-r1", self.out_r1_path,
            "--out-r2", self.out_r2_path
        ]

        with patch.object(sys, 'argv', test_args):
            with patch('sys.stdout', new=StringIO()):
                # The main function calls sys.exit(1) on missing files
                with self.assertRaises(SystemExit) as cm:
                    main()
                self.assertEqual(cm.exception.code, 1)

if __name__ == '__main__':
    unittest.main()