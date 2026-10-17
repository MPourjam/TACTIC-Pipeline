#!/usr/bin/env python3
"""
Analyze read lengths in untrimmed_fastq directories for all mocks.
Reports the read length for each R1/R2 pair file.
"""

from __future__ import annotations

import csv
import gzip
from collections import Counter
from pathlib import Path

BASE_DIR = Path(
    "/media/adminuser/MSD-Storage/Mohsen/TACTIC_Benchmark/Mocks/Benchmark_git/All_Analyses/mocks/raw_data"
)
OUTPUT_DIR = Path(
    "/media/adminuser/MSD-Storage/Mohsen/TACTIC_Benchmark/length_survey"
)


def read_lengths(fastq_path: Path, max_reads: int = 1000) -> list[int]:
    """Extract read lengths from a FASTQ or FASTQ.GZ file."""
    opener = gzip.open if fastq_path.suffix == ".gz" else open
    lengths: list[int] = []

    with opener(fastq_path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line_number % 4 == 2:
                lengths.append(len(line.rstrip("\n\r")))
                if len(lengths) >= max_reads:
                    break

    return lengths


def main() -> None:
    untrimmed_dirs = sorted(
        path for path in BASE_DIR.glob("mock-*/untrimmed_fastq") if path.is_dir()
    )

    rows: list[dict[str, object]] = []

    print("=" * 100)
    print("UNTRIMMED FASTQ READ LENGTH SURVEY - ALL MOCKS")
    print("=" * 100)
    print()

    for untrimmed_dir in untrimmed_dirs:
        mock_name = untrimmed_dir.parent.name

        # Find all "read 1" files (both _R1_ and _1 patterns)
        r1_files = sorted(
            path
            for path in untrimmed_dir.iterdir()
            if path.is_file()
            and (
                ("_R1_" in path.name or "_R1." in path.name or path.name.endswith("_1.fastq") or path.name.endswith("_1.fastq.gz"))
            )
            and (path.name.endswith(".fastq") or path.name.endswith(".fastq.gz"))
        )

        if not r1_files:
            print(f"⚠ {mock_name}: No read 1 files found in untrimmed_fastq")
            continue

        print(f"\n{mock_name}")
        print("-" * 100)

        for r1_path in r1_files:
            # Corresponding R2
            if "_R1_" in r1_path.name:
                r2_name = r1_path.name.replace("_R1_", "_R2_")
            elif "_R1." in r1_path.name:
                r2_name = r1_path.name.replace("_R1.", "_R2.")
            elif r1_path.name.endswith("_1.fastq.gz"):
                r2_name = r1_path.name.replace("_1.fastq.gz", "_2.fastq.gz")
            elif r1_path.name.endswith("_1.fastq"):
                r2_name = r1_path.name.replace("_1.fastq", "_2.fastq")
            else:
                continue
            r2_path = untrimmed_dir / r2_name

            if not r2_path.exists():
                print(f"  ✗ {r1_path.name}: No matching R2 file found")
                continue

            # Read lengths
            r1_lengths = read_lengths(r1_path)
            r2_lengths = read_lengths(r2_path)

            if not r1_lengths or not r2_lengths:
                print(f"  ✗ {r1_path.name}/{r2_name}: Empty or unreadable")
                continue

            # Get unique lengths
            r1_unique = sorted(set(r1_lengths))
            r2_unique = sorted(set(r2_lengths))

            # Format output
            r1_str = (
                f"{r1_unique[0]} bp" if len(r1_unique) == 1 else f"{r1_unique[0]}-{r1_unique[-1]} bp"
            )
            r2_str = (
                f"{r2_unique[0]} bp" if len(r2_unique) == 1 else f"{r2_unique[0]}-{r2_unique[-1]} bp"
            )

            print(f"  {r1_path.name:<50} R1: {r1_str:<15} ({len(r1_lengths):>6,} reads sampled)")
            print(f"  {r2_name:<50} R2: {r2_str:<15} ({len(r2_lengths):>6,} reads sampled)")

            rows.append(
                {
                    "mock": mock_name,
                    "r1_file": r1_path.name,
                    "r1_length_min": r1_unique[0],
                    "r1_length_max": r1_unique[-1],
                    "r1_reads_sampled": len(r1_lengths),
                    "r2_file": r2_name,
                    "r2_length_min": r2_unique[0],
                    "r2_length_max": r2_unique[-1],
                    "r2_reads_sampled": len(r2_lengths),
                }
            )

    # Write CSV
    output_csv = OUTPUT_DIR / "untrimmed_fastq_read_length_summary.csv"
    with output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "mock",
                "r1_file",
                "r1_length_min",
                "r1_length_max",
                "r1_reads_sampled",
                "r2_file",
                "r2_length_min",
                "r2_length_max",
                "r2_reads_sampled",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print()
    print("=" * 100)
    print(f"✓ Report saved to: {output_csv}")
    print("=" * 100)


if __name__ == "__main__":
    main()
