#!/usr/bin/env python3
"""Verify every ZIP member against the SHA-256 manifest without extracting it."""

import argparse
import hashlib
import json
import time
import zipfile
from pathlib import Path


def verify(path):
    started = time.monotonic()
    checked_bytes = 0
    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()
        if len(names) != len(set(names)):
            raise ValueError("Duplicate ZIP member names")
        expected = {}
        for line in archive.read("input_manifest/SHA256SUMS").decode().splitlines():
            digest, name = line.split("  ", 1)
            expected[name] = digest
        if set(names) != set(expected) | {"input_manifest/SHA256SUMS"}:
            raise ValueError("ZIP inventory differs from SHA256SUMS")
        for index, (name, digest) in enumerate(expected.items(), 1):
            member = Path(name)
            if member.is_absolute() or ".." in member.parts:
                raise ValueError(f"Unsafe member path: {name}")
            actual = hashlib.sha256()
            with archive.open(name) as handle:
                while chunk := handle.read(4 * 1024 * 1024):
                    actual.update(chunk)
                    checked_bytes += len(chunk)
            if actual.hexdigest() != digest:
                raise ValueError(f"SHA-256 mismatch: {name}")
            if index % 200 == 0:
                print(f"Verified {index}/{len(expected)} members", flush=True)
    return {
        "archive": str(Path(path).resolve()),
        "verified_payload_and_manifest_files": len(expected),
        "zip_members": len(names),
        "verified_uncompressed_bytes": checked_bytes,
        "zip_crc_and_sha256": "passed",
        "elapsed_seconds": round(time.monotonic() - started, 2),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    args = parser.parse_args()
    print(json.dumps(verify(args.archive), indent=2))


if __name__ == "__main__":
    main()
