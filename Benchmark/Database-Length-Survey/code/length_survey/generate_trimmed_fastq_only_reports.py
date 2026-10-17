#!/usr/bin/env python3
"""
Generate one report per FASTQ file directly inside each mock's trimmed_fastq directory.

Only files in:
  .../mocks/raw_data/mock-*/trimmed_fastq/
are included. Nested FASTQs elsewhere are intentionally ignored.
"""

from __future__ import annotations

import csv
import gzip
import statistics
from collections import Counter
from pathlib import Path


BASE_DIR = Path(
    "/media/adminuser/MSD-Storage/Mohsen/TACTIC_Benchmark/Mocks/Benchmark_git/All_Analyses/mocks/raw_data"
)
OUTPUT_DIR = Path(
    "/media/adminuser/MSD-Storage/Mohsen/TACTIC_Benchmark/length_survey/mockrobiota_trimmed_fastq_reports"
)


def read_lengths(fastq_path: Path) -> list[int]:
    opener = gzip.open if fastq_path.suffix == ".gz" else open
    lengths: list[int] = []

    with opener(fastq_path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line_number % 4 == 2:
                lengths.append(len(line.rstrip("\n\r")))

    return lengths


def summarize_lengths(lengths: list[int]) -> dict[str, object]:
    if not lengths:
        return {
            "count": 0,
            "min": None,
            "max": None,
            "mean": None,
            "median": None,
            "stdev": None,
            "unique_lengths": [],
            "distribution": [],
        }

    counter = Counter(lengths)
    return {
        "count": len(lengths),
        "min": min(lengths),
        "max": max(lengths),
        "mean": statistics.mean(lengths),
        "median": statistics.median(lengths),
        "stdev": statistics.stdev(lengths) if len(lengths) > 1 else 0.0,
        "unique_lengths": sorted(counter),
        "distribution": sorted(counter.items()),
    }


def mock_name_from_path(fastq_path: Path) -> str:
    return fastq_path.relative_to(BASE_DIR).parts[0]


def report_text(fastq_path: Path, stats: dict[str, object]) -> str:
    unique_lengths = stats["unique_lengths"]
    distribution = stats["distribution"]

    lines = [
        "FASTQ FILE REPORT",
        f"{'=' * 80}",
        f"File: {fastq_path}",
        f"Mock: {mock_name_from_path(fastq_path)}",
        f"Directory: trimmed_fastq",
        f"Total reads: {stats['count']}",
    ]

    if stats["count"] == 0:
        lines.append("Status: EMPTY OR UNREADABLE")
        return "\n".join(lines) + "\n"

    lines.extend(
        [
            f"Length range: {stats['min']}-{stats['max']} bp",
            f"Mean length: {stats['mean']:.2f} bp",
            f"Median length: {stats['median']:.2f} bp",
            f"Std dev: {stats['stdev']:.2f}",
            f"Unique lengths: {len(unique_lengths)}",
            "",
            "Length distribution:",
        ]
    )

    for length, count in distribution:
        pct = (count / stats["count"]) * 100
        lines.append(f"  {length} bp: {count} reads ({pct:.2f}%)")

    lines.append("")
    if len(unique_lengths) == 1:
        lines.append(f"Consistency: CONSISTENT ({unique_lengths[0]} bp only)")
    else:
        lines.append(f"Consistency: VARIABLE ({len(unique_lengths)} lengths observed)")

    return "\n".join(lines) + "\n"


def main() -> None:
    trimmed_dirs = sorted(
        path for path in BASE_DIR.glob("mock-*/trimmed_fastq") if path.is_dir()
    )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    reports_dir = OUTPUT_DIR / "per_file_reports"
    reports_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []

    for trimmed_dir in trimmed_dirs:
        mock_name = trimmed_dir.parent.name
        fastq_files = sorted(
            path for path in trimmed_dir.iterdir()
            if path.is_file()
            and (
                path.name.endswith(".fastq")
                or path.name.endswith(".fastq.gz")
                or path.name.endswith(".fq")
                or path.name.endswith(".fq.gz")
            )
        )

        for fastq_path in fastq_files:
            lengths = read_lengths(fastq_path)
            stats = summarize_lengths(lengths)

            report_name = fastq_path.relative_to(BASE_DIR).as_posix().replace("/", "__")
            report_path = reports_dir / f"{report_name}.txt"
            report_path.write_text(report_text(fastq_path, stats))

            rows.append(
                {
                    "mock": mock_name,
                    "file": str(fastq_path),
                    "reads": stats["count"],
                    "min_length": stats["min"],
                    "max_length": stats["max"],
                    "mean_length": f"{stats['mean']:.2f}" if stats["mean"] is not None else "",
                    "median_length": f"{stats['median']:.2f}" if stats["median"] is not None else "",
                    "stdev": f"{stats['stdev']:.2f}" if stats["stdev"] is not None else "",
                    "unique_lengths": len(stats["unique_lengths"]),
                    "report": str(report_path),
                }
            )

    index_path = OUTPUT_DIR / "trimmed_fastq_file_report_index.csv"
    with index_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "mock",
                "file",
                "reads",
                "min_length",
                "max_length",
                "mean_length",
                "median_length",
                "stdev",
                "unique_lengths",
                "report",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    summary_path = OUTPUT_DIR / "trimmed_fastq_file_report_summary.txt"
    total_reads = sum(row["reads"] for row in rows)
    variable_files = sum(1 for row in rows if row["unique_lengths"] > 1)
    summary_path.write_text(
        "\n".join(
            [
                "TRIMMED FASTQ FILE REPORT SUMMARY",
                f"{'=' * 80}",
                f"Input directory: {BASE_DIR}",
                f"Trimmed_fastq directories found: {len(trimmed_dirs)}",
                f"FASTQ files reported: {len(rows)}",
                f"Total reads parsed: {total_reads}",
                f"Files with variable lengths: {variable_files}",
                f"Per-file reports: {reports_dir}",
                f"Index CSV: {index_path}",
                "",
                "Scope constraint:",
                "  Only files directly inside each mock's trimmed_fastq folder were processed.",
            ]
        )
        + "\n"
    )

    print(summary_path)
    print(index_path)
    print(reports_dir)


if __name__ == "__main__":
    main()