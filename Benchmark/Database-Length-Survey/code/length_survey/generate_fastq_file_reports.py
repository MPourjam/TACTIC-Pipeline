#!/usr/bin/env python3
"""
Generate a per-fastq report for every FASTQ file found under the mockrobiota data tree.

The script writes one text report per input file plus a CSV index with the
main statistics for quick comparison.
"""

from __future__ import annotations

import csv
import gzip
import statistics
from collections import Counter
from pathlib import Path


BASE_DIR = Path(
    "/media/adminuser/MSD-Storage/Mohsen/TACTIC_Benchmark/Mocks/Picked DO NOT ENTER ME DO NOT ENTER ME DO NOT ENTER ME/mockrobiota/data"
)
OUTPUT_DIR = Path(
    "/media/adminuser/MSD-Storage/Mohsen/TACTIC_Benchmark/length_survey/mockrobiota_fastq_reports"
)


def read_lengths(fastq_path: Path) -> list[int]:
    """Return all read lengths from a FASTQ or FASTQ.GZ file."""

    opener = gzip.open if fastq_path.suffix == ".gz" else open
    lengths: list[int] = []

    with opener(fastq_path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line_number % 4 == 2:
                lengths.append(len(line.rstrip("\n\r")))

    return lengths


def summarize_lengths(lengths: list[int]) -> dict[str, object]:
    """Compute length statistics for a FASTQ file."""

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

    length_counter = Counter(lengths)
    return {
        "count": len(lengths),
        "min": min(lengths),
        "max": max(lengths),
        "mean": statistics.mean(lengths),
        "median": statistics.median(lengths),
        "stdev": statistics.stdev(lengths) if len(lengths) > 1 else 0.0,
        "unique_lengths": sorted(length_counter),
        "distribution": sorted(length_counter.items()),
    }


def relative_group(fastq_path: Path) -> str:
    """Return the immediate mock directory the file belongs to."""

    try:
        return fastq_path.relative_to(BASE_DIR).parts[0]
    except ValueError:
        return fastq_path.parent.name


def report_text(fastq_path: Path, stats: dict[str, object]) -> str:
    """Format a standalone report for one FASTQ file."""

    unique_lengths = stats["unique_lengths"]
    distribution = stats["distribution"]

    lines = [
        "FASTQ FILE REPORT",
        f"{'=' * 80}",
        f"File: {fastq_path}",
        f"Mock group: {relative_group(fastq_path)}",
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

    if len(unique_lengths) == 1:
        lines.append("")
        lines.append(f"Consistency: CONSISTENT ({unique_lengths[0]} bp only)")
    else:
        lines.append("")
        lines.append(f"Consistency: VARIABLE ({len(unique_lengths)} lengths observed)")

    return "\n".join(lines) + "\n"


def main() -> None:
    fastq_files = sorted(
        path for path in BASE_DIR.rglob("*")
        if path.is_file() and (
            path.name.endswith(".fastq")
            or path.name.endswith(".fastq.gz")
            or path.name.endswith(".fq")
            or path.name.endswith(".fq.gz")
        )
    )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    reports_dir = OUTPUT_DIR / "per_file_reports"
    reports_dir.mkdir(parents=True, exist_ok=True)

    index_rows: list[dict[str, object]] = []

    for fastq_path in fastq_files:
        lengths = read_lengths(fastq_path)
        stats = summarize_lengths(lengths)

        report_name = fastq_path.relative_to(BASE_DIR).as_posix().replace("/", "__")
        report_path = reports_dir / f"{report_name}.txt"
        report_path.write_text(report_text(fastq_path, stats))

        index_rows.append(
            {
                "file": str(fastq_path),
                "mock_group": relative_group(fastq_path),
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

    index_path = OUTPUT_DIR / "fastq_file_report_index.csv"
    with index_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "file",
                "mock_group",
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
        writer.writerows(index_rows)

    summary_path = OUTPUT_DIR / "fastq_file_report_summary.txt"
    total_reads = sum(row["reads"] for row in index_rows)
    variable_files = sum(1 for row in index_rows if row["unique_lengths"] > 1)

    summary_path.write_text(
        "\n".join(
            [
                "FASTQ FILE REPORT SUMMARY",
                f"{'=' * 80}",
                f"Input directory: {BASE_DIR}",
                f"Total FASTQ files: {len(index_rows)}",
                f"Total reads parsed: {total_reads}",
                f"Files with variable lengths: {variable_files}",
                f"Per-file reports: {reports_dir}",
                f"Index CSV: {index_path}",
                "",
                "Report naming convention:",
                "  <relative-path-with-slashes-replaced-by-double-underscores>.txt",
            ]
        )
        + "\n"
    )

    print(summary_path)
    print(index_path)
    print(reports_dir)


if __name__ == "__main__":
    main()