"""Standard-library input readers shared by packaging and reproduction."""

import csv
import gzip
import hashlib
import json
import math
from collections import Counter
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(4 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(path, rows, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def summarize(counter):
    n = sum(counter.values())
    if not n:
        return {"count": 0, "min": None, "max": None, "mean": None,
                "median": None, "stdev": None, "unique_lengths": 0, "q1": None, "q3": None}
    total = sum(k * v for k, v in counter.items())
    mean = total / n

    def at(index):
        cumulative = 0
        for length, count in sorted(counter.items()):
            cumulative += count
            if cumulative > index:
                return length
        raise IndexError(index)

    variance = sum(v * (k - mean) ** 2 for k, v in counter.items()) / (n - 1) if n > 1 else 0
    return {
        "count": n, "min": min(counter), "max": max(counter), "mean": mean,
        "median": (at((n - 1) // 2) + at(n // 2)) / 2,
        "stdev": math.sqrt(variance), "unique_lengths": len(counter),
        "q1": at(n // 4) if n > 3 else min(counter),
        "q3": at(3 * n // 4) if n > 3 else max(counter),
    }


def scan(path, kind, limit=0, open_file=None):
    path = Path(path)
    opener = open_file or (gzip.open if path.suffix == ".gz" else open)
    counter = Counter()
    sampled = {1000: Counter(), 10000: Counter()}
    with opener(path, "rb") as handle:
        if kind == "fastq":
            count = 0
            while header := handle.readline():
                sequence, plus, quality = handle.readline(), handle.readline(), handle.readline()
                if not header.startswith(b"@") or not plus.startswith(b"+") or not quality:
                    raise ValueError(f"Invalid FASTQ record {count + 1}: {path}")
                sequence = sequence.rstrip(b"\r\n")
                if len(sequence) != len(quality.rstrip(b"\r\n")):
                    raise ValueError(f"Sequence/quality length mismatch: {path}, record {count + 1}")
                count += 1
                length = len(sequence)
                counter[length] += 1
                for size, subset in sampled.items():
                    if count <= size:
                        subset[length] += 1
                if limit and count >= limit:
                    break
        else:
            length = 0
            seen_header = False
            for line in handle:
                line = line.strip()
                if line.startswith(b">"):
                    if length:
                        counter[length] += 1
                    length = 0
                    seen_header = True
                else:
                    if line and not seen_header:
                        raise ValueError(f"FASTA sequence before header: {path}")
                    length += len(line)
            if length:
                counter[length] += 1
    return {
        "stats": summarize(counter),
        "distribution": dict(sorted(counter.items())),
        "samples": {str(size): summarize(c) for size, c in sampled.items() if c},
    }


def expected_matches(actual, expected):
    differences = []
    for key, value in expected.items():
        if value is None:
            if actual.get(key) is not None:
                differences.append(key)
        elif key == "distribution":
            if actual.get(key) != value:
                differences.append(key)
        elif key.endswith("_2dp"):
            if actual.get(key[:-4]) is None or f"{actual[key[:-4]]:.2f}" != str(value):
                differences.append(key)
        elif actual.get(key) is None or not math.isclose(float(actual[key]), float(value), abs_tol=1e-9, rel_tol=1e-9):
            differences.append(key)
    return differences


def worker(path, kind, known_sha256=None):
    return {"sha256": known_sha256 or sha256(path), **scan(path, kind)}


if __name__ == "__main__":
    import sys
    print(json.dumps(worker(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None)))
