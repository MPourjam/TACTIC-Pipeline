#!/usr/bin/env python3
"""Freeze the archived benchmark inputs and recover their local read provenance."""

import argparse
import concurrent.futures
import csv
import gzip
import hashlib
import json
import re
import subprocess
import sys
import zipfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from verify_input_archive import verify


FASTQ_SUFFIXES = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
CHUNK = 4 * 1024 * 1024
SRA_MATE = re.compile(rb"^([SED]RR[0-9]+\.[0-9]+)\.[123]$")
READ_ID_VERSION = "sra_spot_and_mate_v2"


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields=None):
    if fields is None:
        fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(CHUNK):
            digest.update(chunk)
    return digest.hexdigest()


def is_fastq(path):
    return path.is_file() and path.name.lower().endswith(FASTQ_SUFFIXES)


def read_identity(header):
    token = header.strip().split(None, 1)[0][1:]
    token = re.sub(rb"/[123]$", b"", token)
    return SRA_MATE.sub(rb"\1", token)


def read_ids_digest(path):
    digest = hashlib.sha256()
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    with opener(path, "rb") as handle:
        while header := handle.readline():
            for _ in range(3):
                handle.readline()
            digest.update(read_identity(header) + b"\n")
    return digest.hexdigest()


def fastq_stats(path):
    content = hashlib.sha256()
    ids = hashlib.sha256()
    count = total_bases = errors = 0
    minimum = None
    maximum = 0
    first = last = ""
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    with opener(path, "rb") as handle:
        while header := handle.readline():
            sequence = handle.readline()
            separator = handle.readline()
            quality = handle.readline()
            for line in (header, sequence, separator, quality):
                content.update(line)
            seq = sequence.rstrip(b"\r\n")
            qual = quality.rstrip(b"\r\n")
            if not header.startswith(b"@") or not separator.startswith(b"+") or len(seq) != len(qual) or not quality:
                errors += 1
            identity = read_identity(header)
            ids.update(identity + b"\n")
            last = header.rstrip().decode(errors="replace")
            if not first:
                first = last
            length = len(seq)
            count += 1
            total_bases += length
            minimum = length if minimum is None else min(minimum, length)
            maximum = max(maximum, length)
    return {
        "read_count": count,
        "total_bases": total_bases,
        "minimum_read_length": minimum if minimum is not None else "",
        "maximum_read_length": maximum,
        "invalid_fastq_records": errors,
        "uncompressed_sha256": content.hexdigest(),
        "ordered_read_ids_sha256": ids.hexdigest(),
        "read_id_normalization": READ_ID_VERSION,
        "first_header": first,
        "last_header": last,
        "gzip_crc": "passed" if path.name.lower().endswith(".gz") else "not_gzip",
    }


def identify_read(path, cid):
    match = re.search(r"^(.*?)_R([12])(?:_|\.|$)", path.name)
    if match:
        return match.group(1), "R" + match.group(2)
    for word, mate in (("forward", "R1"), ("reverse", "R2"), ("index", "I1"), ("barcodes", "I1")):
        if word in path.name.lower():
            return cid + "_pooled", mate
    return path.name, "unknown"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--source-root", type=Path)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--cache", type=Path, default=Path("/tmp/tacticbench_input_archive_cache.json"))
    parser.add_argument("--inventory-only", action="store_true")
    parser.add_argument("--fastq-worker", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--read-ids-only", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.fastq_worker is not None:
        result = {"ordered_read_ids_sha256": read_ids_digest(args.fastq_worker)} if args.read_ids_only else fastq_stats(args.fastq_worker)
        print(json.dumps(result))
        return
    bundle = args.bundle.resolve()
    analysis_root = (args.source_root or bundle.parent.parent).resolve()
    output = bundle / "data"
    reports = output / "input_manifest"
    reports.mkdir(exist_ok=True)
    archived = bundle / "results/archive_run"
    original = read_tsv(bundle / "metadata/original_run_manifest.tsv")
    source_links = {r["community_id"]: r for r in read_tsv(bundle / "metadata/FASTQ_and_publication_links.tsv")}
    source_root = Path(original[0]["community_path"]).parent
    if args.source_root:
        source_root = analysis_root / "mocks/cleaned_up_202605"
    payload = {}
    fastqs = []
    benchmark = []
    provenance = []
    issues = []
    cache = json.loads(args.cache.read_text()) if args.cache.exists() else {}

    def cached_hash(path):
        path = path.resolve()
        stat = path.stat()
        key = f"{path}|{stat.st_size}|{stat.st_mtime_ns}"
        rec = cache.setdefault(key, {})
        if "sha256" not in rec:
            rec["sha256"] = sha256(path)
        return rec["sha256"]

    def add(path, name, role, cid=""):
        if not path.is_file():
            raise FileNotFoundError(path)
        if name in payload:
            if cached_hash(payload[name]["source"]) != cached_hash(path):
                alternate = "input_manifest/workspace_snapshots/" + name
                if alternate not in payload:
                    payload[alternate] = {"source": path, "role": role + "_workspace_differs_from_bundle", "community_id": cid}
                    issues.append({"severity": "note", "community_id": cid, "issue": "workspace_snapshot_differs", "detail": f"Preserved bundle {name}; current workspace copy is {alternate}"})
            return name
        payload[name] = {"source": path, "role": role, "community_id": cid}
        return name

    def community_add(path, comm, role):
        return add(path, "communities/" + comm.name + "/" + path.relative_to(comm).as_posix(), role, comm.name)

    for path in sorted((output / "communities").rglob("*")):
        if path.is_file() and "__pycache__" not in path.parts:
            add(path, path.relative_to(output).as_posix(), "existing_benchmark_input_snapshot", path.relative_to(output / "communities").parts[0])
    print(f"Retained {len(payload)} files from the existing data bundle", flush=True)

    selected_runs = {}
    for row in original:
        cid = row["community_id"]
        comm = source_root / cid
        old_comm = Path(row["community_path"])
        analysis = comm / Path(row["TIC_table"]).relative_to(old_comm).parent.parent
        selected_runs[cid] = analysis
        for role, value in row.items():
            if role in {"community_id", "community_path"}:
                continue
            source = comm / Path(value).relative_to(old_comm)
            name = "communities/" + cid + "/" + source.relative_to(comm).as_posix()
            bundled = payload.get(name, {}).get("source")
            if bundled is None:
                raise ValueError(f"Missing archived direct input: {name}")
            equal = cached_hash(source) == cached_hash(bundled)
            if not equal:
                raise ValueError(f"Direct input changed since the bundle was created: {source}")
            payload[name]["role"] = "direct_benchmark_" + role
            benchmark.append({"community_id": cid, "input_role": role, "original_source_path": value,
                              "bundled_source_path": str(bundled), "archive_path": name,
                              "sha256": cached_hash(bundled), "original_and_bundle_identical": equal})

        for path in comm.iterdir():
            if path.is_file() and path.suffix.lower() in {".md", ".tsv", ".txt", ".yaml", ".yml", ".csv"}:
                community_add(path, comm, "community_metadata_or_log")
        for path in (comm / "source").rglob("*"):
            if path.is_file() and path.suffix.lower() in {".fasta", ".fa", ".fna", ".gb", ".gbk", ".txt", ".tsv", ".csv", ".log", ".py", ".json"}:
                community_add(path, comm, "expected_reference_provenance")
        for dirname in ("trimmed_fastq", "untrimmed_fastq", "raw_data", "trimmed_but_spiked", "emp-paired-end-sequences", "dada2_results/work", "dada2_results/trimmed"):
            folder = comm / dirname
            if not folder.is_dir():
                continue
            paths = folder.rglob("*") if dirname.startswith("dada2_results/") else folder.iterdir()
            for path in sorted(paths):
                if is_fastq(path):
                    sample, mate = identify_read(path, cid)
                    fastqs.append({"source": path, "community_id": cid, "stage": dirname, "sample_id": sample, "mate": mate})
                elif path.is_file() and (path.suffix.lower() in {".log", ".txt", ".tsv", ".csv", ".yaml", ".yml", ".r", ".sh"} or path.name == "MANIFEST"):
                    community_add(path, comm, "preprocessing_metadata_or_script")
        checkpoint = comm / "dada2_results/dada2_out/seqtab_nochim.rds"
        if checkpoint.exists():
            community_add(checkpoint, comm, "dada2_sequence_table_checkpoint")
        for path in analysis.iterdir():
            if path.is_file() and path.suffix.lower() in {".txt", ".log", ".csv", ".yml", ".yaml"}:
                community_add(path, comm, "selected_tactic_run_provenance")
        for pipeline in ("TIC_Pipeline", "TAC_Pipeline"):
            for path in (analysis / pipeline).iterdir():
                if path.is_file() and (path.name.startswith(("Map-", "ZOTUs-", "SOTUs-Table-TIC-without", "Invalid-Tax", "archaeal", "README")) or path.name in {"otu_clusters.uc", "cluster_centroids.fasta", "TIC-FullTaxonomy-Sequences.fasta"}):
                    if path.suffix.lower() != ".html":
                        community_add(path, comm, "cluster_membership_or_taxonomy_provenance")

    for directory in (bundle / "config", bundle / "environment", bundle / "metadata", bundle / "code"):
        for path in directory.rglob("*"):
            if path.is_file() and "__pycache__" not in path.parts and path.suffix not in {".pyc", ".pyo"}:
                add(path, "input_manifest/context/" + path.relative_to(bundle).as_posix(), "reproduction_context")
    for filename in ("manifest.tsv", "config_input.yaml", "config_effective.yaml", "tacticbench.log", "community_name_mapping.tsv"):
        path = archived / filename
        if path.exists():
            add(path, "input_manifest/context/archive_run/" + filename, "archived_run_record")
    for filename in ("sample_name_mapping.tsv", "expected_sample_selection.tsv", "skipped_communities.tsv"):
        path = archived / "qc" / filename
        if path.exists():
            add(path, "input_manifest/context/archive_run/qc/" + filename, "archived_sample_selection")
    for path in analysis_root.iterdir():
        if path.is_file() and (path.suffix in {".py", ".sh", ".yml", ".yaml"} or path.name in {"run_tactic.log", "Pipeline_log.txt"}):
            if path.name != "qiime2-2024.10_env.yml":
                add(path, "input_manifest/upstream_scripts/" + path.name, "upstream_driver_or_parameter_snapshot")
    for path in source_root.iterdir():
        if path.is_file() and path.suffix in {".sh", ".py"}:
            add(path, "input_manifest/upstream_scripts/cleaned_up_202605/" + path.name, "reference_or_read_preparation_script")
    dbdir = analysis_root / "databases/SILVA/138.2"
    for filename in ("SILVA.arb.gz", "silva-bac-16s-id90.fasta", "silva-arc-16s-id95.fasta", "SILVA.sidx.success"):
        path = dbdir / filename
        if path.exists():
            add(path, "input_manifest/upstream_databases/SILVA/138.2/" + filename, "upstream_silva_reference")

    priority = {"trimmed_fastq": 0, "untrimmed_fastq": 1, "trimmed_but_spiked": 2, "raw_data": 3,
                "dada2_results/trimmed": 4, "dada2_results/work": 5, "emp-paired-end-sequences": 6}
    fastqs.sort(key=lambda r: (priority[r["stage"]], str(r["source"])))
    hash_tasks = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        for rec in fastqs:
            hash_tasks[pool.submit(cached_hash, rec["source"])] = rec
        for index, future in enumerate(concurrent.futures.as_completed(hash_tasks), 1):
            hash_tasks[future]["sha256"] = future.result()
            if index % 60 == 0:
                print(f"Hashed {index}/{len(fastqs)} FASTQ locations", flush=True)
    stored = {}
    for rec in fastqs:
        digest = rec["sha256"]
        if digest in stored:
            rec["archive_path"] = stored[digest]["archive_path"]
            rec["storage"] = "byte_identical_alias"
        else:
            rec["archive_path"] = community_add(rec["source"], source_root / rec["community_id"], "fastq_" + rec["stage"].replace("/", "_"))
            rec["storage"] = "stored"
            stored[digest] = rec
    print(f"FASTQs: {len(fastqs)} original locations, {len(stored)} distinct payloads", flush=True)
    args.cache.write_text(json.dumps(cache))

    def inspect(rec):
        path = rec["source"].resolve()
        stat = path.stat()
        key = f"{path}|{stat.st_size}|{stat.st_mtime_ns}"
        info = cache[key].get("fastq")
        if info is None:
            # Separate interpreters avoid the GIL bottleneck when reading records.
            worker = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--fastq-worker", str(path)],
                                    capture_output=True, text=True, check=True)
            info = json.loads(worker.stdout)
        elif info.get("read_id_normalization") != READ_ID_VERSION:
            info = dict(info)
            first = info["first_header"].split(None, 1)[0][1:].encode()
            if SRA_MATE.fullmatch(first):
                worker = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--fastq-worker", str(path), "--read-ids-only"],
                                        capture_output=True, text=True, check=True)
                info.update(json.loads(worker.stdout))
            info["read_id_normalization"] = READ_ID_VERSION
        return rec["sha256"], key, info

    stats = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        tasks = [pool.submit(inspect, rec) for rec in stored.values()]
        for index, future in enumerate(concurrent.futures.as_completed(tasks), 1):
            digest, key, stat = future.result()
            cache[key]["fastq"] = stat
            stats[digest] = stat
            if index % 10 == 0:
                args.cache.write_text(json.dumps(cache))
                print(f"Validated FASTQ records and gzip CRCs: {index}/{len(stored)}", flush=True)
    args.cache.write_text(json.dumps(cache))
    lookup = {str(rec["source"].absolute()): rec for rec in fastqs}
    fastq_rows = []
    pair_groups = defaultdict(dict)
    for rec in fastqs:
        path = rec["source"]
        info = stats[rec["sha256"]]
        links = source_links.get(rec["community_id"], {})
        urls = links.get("fastq_source", "").split("; ")
        source_url = "; ".join(url for url in urls if rec["sample_id"].startswith("SRR") and rec["sample_id"] in url)
        if not source_url:
            marker = {"R1": "forward", "R2": "reverse", "I1": "index"}.get(rec["mate"], "")
            source_url = "; ".join(url for url in urls if marker and marker in url)
        accession = re.match(r"@([SED]RR[0-9]+)", info["first_header"])
        fastq_rows.append({"community_id": rec["community_id"], "stage": rec["stage"], "sample_id": rec["sample_id"],
                           "mate": rec["mate"], "source_path": str(path), "archive_path": rec["archive_path"],
                           "storage": rec["storage"], "size_bytes": path.stat().st_size, "sha256": rec["sha256"],
                           "raw_data_origin_url": source_url,
                           "sra_accession_in_read_header": accession[1] if accession else "",
                           "observed_sra_record_url": "https://www.ncbi.nlm.nih.gov/sra/?term=" + accession[1] if accession else "",
                           "publication_doi_or_url": links.get("publication_doi_or_url", ""),
                           "url_evidence": "existing_bundle_metadata; remote URLs not revalidated", **info})
        key = (rec["community_id"], rec["stage"], rec["sample_id"])
        # Plain/gz duplicates are retained, but pair checks prefer the gz copy.
        if rec["mate"] not in pair_groups[key] or path.suffix == ".gz":
            pair_groups[key][rec["mate"]] = rec
        if info["invalid_fastq_records"]:
            issues.append({"severity": "error", "community_id": rec["community_id"], "issue": "invalid_fastq_records", "detail": str(path)})
    pair_rows = []
    for (cid, stage, sample), mates in sorted(pair_groups.items()):
        if "R1" not in mates or "R2" not in mates:
            issues.append({"severity": "error", "community_id": cid, "issue": "missing_fastq_mate", "detail": f"{stage}/{sample}"})
            continue
        r1, r2 = mates["R1"], mates["R2"]
        a, b = stats[r1["sha256"]], stats[r2["sha256"]]
        ok = a["read_count"] == b["read_count"] and a["ordered_read_ids_sha256"] == b["ordered_read_ids_sha256"]
        pair_rows.append({"community_id": cid, "stage": stage, "sample_id": sample,
                          "R1_archive_path": r1["archive_path"], "R2_archive_path": r2["archive_path"],
                          "R1_read_count": a["read_count"], "R2_read_count": b["read_count"],
                          "paired_read_ids_in_order": ok,
                          "index_archive_path": mates.get("I1", {}).get("archive_path", "")})
        if not ok:
            issues.append({"severity": "error", "community_id": cid, "issue": "fastq_mate_mismatch", "detail": f"{stage}/{sample}"})

    pipeline_rows = []
    pipeline_lookup = {}
    zero_length_rows = []
    for rec in stored.values():
        if stats[rec["sha256"]]["minimum_read_length"] != 0:
            continue
        opener = gzip.open if rec["source"].suffix == ".gz" else open
        with opener(rec["source"], "rb") as handle:
            while header := handle.readline():
                sequence = handle.readline().rstrip(b"\r\n")
                handle.readline()
                handle.readline()
                if not sequence:
                    zero_length_rows.append({"community_id": rec["community_id"], "sample_id": rec["sample_id"],
                                             "mate": rec["mate"], "archive_path": rec["archive_path"],
                                             "read_id": read_identity(header).decode(errors="replace")})
    for cid, analysis in selected_runs.items():
        comm = source_root / cid
        log = analysis / "Pipeline_log.txt"
        lines = log.read_text(errors="replace").splitlines()
        counts = {}
        bowtie_counts = {}
        processed = set()
        for lineno, line in enumerate(lines, 1):
            match = re.search(r"Sample '([^']+)': (\d+) spike reads removed from (\d+) total reads", line)
            if match:
                counts[match[1]] = (int(match[3]), int(match[2]))
            bowtie = re.search(r"Warning in running command .*bowtie2.* -1 (\S+) -2 ", line)
            if bowtie:
                sample = Path(bowtie[1]).name.split("_R1", 1)[0]
                block = "\n".join(lines[lineno:lineno + 22])
                total = re.search(r"(?m)^(\d+) reads; of these:", block)
                unpaired = re.search(r"(?m)^\s*(\d+) .*were unpaired;", block)
                if total:
                    bowtie_counts[sample] = (int(total[1]), int(unpaired[1]) if unpaired else 0)
                    provenance.append({"community_id": cid, "evidence_kind": "bowtie2_record_counts", "archive_path": "communities/" + cid + "/" + log.relative_to(comm).as_posix(), "line_number": lineno + 1, "text": block.split("overall alignment rate", 1)[0] + "overall alignment rate"})
            if "Gathering sequence files" in line or "samples were collected" in line or match:
                provenance.append({"community_id": cid, "evidence_kind": "tactic_input_or_read_count", "archive_path": "communities/" + cid + "/" + log.relative_to(comm).as_posix(), "line_number": lineno, "text": line})
        for lineno, line in enumerate(lines, 1):
            match = re.search(r"-fastq_mergepairs\s+(\S+)\s+-reverse\s+(\S+)", line)
            if not match:
                continue
            sample = Path(match[1]).name.split("_R1", 1)[0]
            if sample in processed:
                continue
            processed.add(sample)
            provenance.append({"community_id": cid, "evidence_kind": "tactic_merge_command", "archive_path": "communities/" + cid + "/" + log.relative_to(comm).as_posix(), "line_number": lineno, "text": line})
            for logged in match.groups():
                relative = logged.split("/mocks/cleaned_up_202605/" + cid + "/", 1)[-1]
                folder = comm / Path(relative).parent
                if folder.is_dir():
                    for path in folder.iterdir():
                        if path.is_file() and path.suffix in {".txt", ".csv", ".yaml", ".yml", ".log"}:
                            community_add(path, comm, "selected_sample_preprocessing_evidence")
            recovered = []
            for mate in ("R1", "R2"):
                if cid in {"mock-18", "mock-19"}:
                    path = comm / "trimmed_but_spiked" / f"{sample}_{mate}_trimmed_with_spike.fastq.gz"
                else:
                    path = comm / "trimmed_fastq" / f"{sample}_{mate}_trimmed.fastq.gz"
                recovered.append(lookup[str(path.absolute())])
            dada = [lookup[str((comm / "dada2_results/trimmed" / f"{sample}_{mate}_trim.fastq.gz").absolute())] for mate in ("R1", "R2")]
            input_count, removed = counts[sample]
            total_records, unpaired_records = bowtie_counts.get(sample, (input_count, 0))
            observed = stats[recovered[0]["sha256"]]["read_count"]
            row = {"community_id": cid, "upstream_sample_id": sample, "selected_analysis_directory": str(analysis),
                   "pipeline_log_archive_path": "communities/" + cid + "/" + log.relative_to(comm).as_posix(),
                   "merge_command_line_number": lineno, "logged_R1_working_path": match[1], "logged_R2_working_path": match[2],
                   "TACTIC_recovered_R1_source": str(recovered[0]["source"]), "TACTIC_recovered_R2_source": str(recovered[1]["source"]),
                   "TACTIC_R1_archive_path": recovered[0]["archive_path"], "TACTIC_R2_archive_path": recovered[1]["archive_path"],
                   "TACTIC_reported_paired_reads": input_count, "TACTIC_logged_total_input_records": total_records,
                   "TACTIC_logged_unpaired_reads": unpaired_records, "TACTIC_recovered_pair_records": observed,
                   "TACTIC_logged_spike_pairs_removed": removed, "TACTIC_logged_post_spike_pairs": input_count - removed,
                   "TACTIC_input_count_matches_log": observed == total_records,
                   "DADA2_R1_source": str(dada[0]["source"]), "DADA2_R2_source": str(dada[1]["source"]),
                   "DADA2_R1_archive_path": dada[0]["archive_path"], "DADA2_R2_archive_path": dada[1]["archive_path"],
                   "DADA2_input_pairs": stats[dada[0]["sha256"]]["read_count"],
                   "historical_FASTQ_checksum_available": False,
                   "provenance": "selected_analysis_log_and_read_count; renamed_with_spike_copy" if cid in {"mock-18", "mock-19"} else "selected_analysis_log_and_read_count; original_trimmed_name"}
            pipeline_rows.append(row)
            pipeline_lookup[(cid, sample)] = row
            if observed != total_records:
                issues.append({"severity": "error", "community_id": cid, "issue": "tactic_input_count_differs_from_log", "detail": sample})
        for path in (comm / "trimmed_fastq").glob("*.cutadapt.log"):
            for lineno, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
                if line.startswith("Command line parameters:"):
                    provenance.append({"community_id": cid, "evidence_kind": "cutadapt_command", "archive_path": "communities/" + cid + "/" + path.relative_to(comm).as_posix(), "line_number": lineno, "text": line})

    sample_links = []
    for rec in read_tsv(archived / "qc/sample_name_mapping.tsv"):
        key = (rec["community_id"], rec["inferred_sample"])
        run = pipeline_lookup[key]
        prefix = "DADA2" if rec["method"] == "ASV" else "TACTIC"
        sample_links.append({**rec, "upstream_pipeline": prefix, "R1_archive_path": run[prefix + "_R1_archive_path"],
                             "R2_archive_path": run[prefix + "_R2_archive_path"], "selected_analysis_directory": run["selected_analysis_directory"]})
    issues.extend([
        {"severity": "limitation", "community_id": "all", "issue": "historical_fastq_checksums_absent", "detail": "Logs record sample paths and read counts, not historical FASTQ SHA-256 values. Current recovered files are frozen by this archive."},
        {"severity": "limitation", "community_id": "mock-4;mock-5", "issue": "tic_archaeal_corrections", "detail": "The selected benchmark uses TIC files corrected after the original TACTIC run. Preserve these exact tables, FASTA files, and membership maps; an unmodified upstream rerun need not recreate them."},
        {"severity": "limitation", "community_id": "all", "issue": "upstream_container_not_archived", "detail": "Driver snapshot names ghcr.io/mpourjam/tactic-pipeline:0.7.5-alpha. Container digest, image, licensed binaries, and container-internal spike reference FASTA are not frozen here. No upstream pipeline was run."},
        {"severity": "note", "community_id": "mock-18;mock-19", "issue": "different_tactic_and_dada_inputs", "detail": "TACTIC logs record input before internal spike removal. Copies renamed *_with_spike.fastq.gz are retained under trimmed_but_spiked; later DADA2 input copies are spikefree. Both stages are included."},
        {"severity": "note", "community_id": "all", "issue": "legacy_metadata_filenames", "detail": "Some original MANIFEST and driver snapshots contain obsolete names or early exit statements. pipeline_inputs.tsv and sample_input_links.tsv identify the resolved files; scripts are preserved as evidence, not an automatic upstream replay recipe."},
        {"severity": "note", "community_id": "mock-ZIEL_1", "issue": "obsolete_analysis_memberships_relocated", "detail": "Older membership-only Analysis_20260510_115118 and Analysis_20260511_170504 directories were preserved under input_manifest/legacy_analysis_memberships. Only the selected Analysis_20260511_170817 run remains under data/communities, preventing mtime-based discovery from skipping four methods."},
        {"severity": "note", "community_id": "mock-ZIEL_2", "issue": "bowtie2_pair_vs_record_count", "detail": "SRR13005964 retains 147918 FASTQ pair records. The recorded Bowtie2 summary also says 147918 total, consisting of 147917 paired and one unpaired read. The R2 mate of SRR13005964.10542 has zero length; TACTIC's summary reports the paired subset. Original files are preserved."},
        {"severity": "limitation", "community_id": "all", "issue": "dada2_invocation_not_fully_recorded", "detail": "Generated per-community R scripts, pair lists, and driver snapshots are included, but the exact historical per-run positional R arguments and an immutable DADA2 container digest were not retained. Driver scripts contain early exits and are provenance snapshots, not a verified upstream replay recipe."},
    ])
    for filename in ("package_input_data.py", "verify_input_archive.py", "audit_benchmark_inputs.py"):
        path = Path(__file__).parent / filename
        if path.exists():
            add(path, "input_manifest/scripts/" + filename, "archive_build_or_verification_script")
    report_files = {
        "benchmark_inputs.tsv": benchmark,
        "fastq_files.tsv": fastq_rows,
        "fastq_pairs.tsv": pair_rows,
        "pipeline_inputs.tsv": pipeline_rows,
        "sample_input_links.tsv": sample_links,
        "provenance_evidence.tsv": provenance,
        "zero_length_reads.tsv": zero_length_rows,
        "limitations_and_findings.tsv": issues,
    }
    for filename, records in report_files.items():
        write_tsv(reports / filename, records)
    preparation_report = reports / "benchmark_preparation_validation.tsv"
    preparation = read_tsv(preparation_report) if preparation_report.is_file() else []
    if not preparation or any(r["passed"] != "True" for r in preparation):
        raise ValueError("Run scripts/audit_benchmark_inputs.py successfully before packaging")
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "communities": len(original), "samples": len(pipeline_rows), "direct_benchmark_input_files": len(benchmark),
        "direct_benchmark_input_copies_identical": all(r["original_and_bundle_identical"] for r in benchmark),
        "original_fastq_locations": len(fastqs), "distinct_fastq_payloads": len(stored),
        "fastq_pairs_checked_in_original_layouts": len(pair_rows),
        "fastq_pair_checks_passed": all(r["paired_read_ids_in_order"] for r in pair_rows),
        "invalid_fastq_records": sum(r["invalid_fastq_records"] for r in stats.values()),
        "tactic_samples_matching_logged_read_count": sum(r["TACTIC_input_count_matches_log"] for r in pipeline_rows),
        "sample_method_input_links": len(sample_links),
        "benchmark_preparation_artifacts_matching_archive": len(preparation),
        "fastq_stages": dict(Counter(r["stage"] for r in fastqs)),
        "upstream_tactic_or_dada2_executed": False,
        "errors": [r for r in issues if r["severity"] == "error"],
    }
    (reports / "validation.json").write_text(json.dumps(summary, indent=2) + "\n")
    readme = f"""# Benchmark inputs and FASTQ provenance

This archive extends the TACTICBench reproduction bundle for
`TACTICBENCH_RESULTS-TIC-Archaea-mock4-mock-5_SRP`.
It freezes {len(original)} communities, {len(pipeline_rows)} samples, {len(benchmark)} direct benchmark input
files, and {len(stored)} distinct FASTQ payloads recovered from {len(fastqs)} original file locations.
No TACTIC, DADA2, denoising, clustering, or tree pipeline was run to assemble it.

## Extract and reproduce the benchmark

From the reproduction bundle directory:

```bash
python scripts/verify_input_archive.py data/tacticbench_inputs_with_fastq.zip
unzip -n data/tacticbench_inputs_with_fastq.zip -d data
conda activate tacticbench-reproduction
python -m pip install --no-deps -e code/tacticbench
python -m tacticbench.cli run-all --config config/run_config.yaml
```

The ZIP has `communities/` and `input_manifest/` at its root. Extracting into `data/`
restores the benchmark's expected directory layout. Existing files are retained
with `unzip -n`. Results are written to `results/reproduced_run/` by the supplied
run config. `data/communities/` already contains the smaller benchmark inputs;
extract the archive to add the reads and additional provenance.

## What the code consumes

`discover.py` selects `source/taxonomy.tsv`, `source/expected-sequences.v4.fasta`
(with the configured full-reference fallback), DADA2 `ASV_counts.tsv`/`ASVs.fasta`,
and the feature tables/representative sequences from ZOTU, OTU, TIC, and TAC.
The original run's manifest pins the selected `Analysis_*` directories; do not
replace them with a later analysis selected by modification time. Two obsolete
ZIEL-I membership-only directories have been preserved under
`input_manifest/legacy_analysis_memberships/`, leaving only the selected analysis
in its usable community data directory. TIC/TAC
membership maps are also needed by cluster diagnostics. `prepare_expected.py`
and `prepare_methods.py` select samples, reconcile their names, and normalize
the supplied tables. TACTICBench does not read FASTQs itself.

The package retains the exact direct inputs and additional full reference FASTAs,
sample metadata, raw reads, trimmed reads, selected analysis logs and YAML files,
cluster membership and taxonomy evidence, generated DADA2 R scripts and pair lists,
and locally available SILVA 138.2 source databases. Copies of the bundle's code,
configs, environment, and archived run metadata are under `input_manifest/context/`.
Original publication and read URLs remain in the metadata snapshots; publication
PDFs and upstream container binaries are not included.

## Manifest guide

| File | Meaning |
|---|---|
| `files.tsv` | Every physical payload with original source path, archive path, role, byte size, timestamp, and SHA-256. |
| `benchmark_inputs.tsv` | All {len(benchmark)} files explicitly selected by the archived run manifest, with identity checks against the bundle. |
| `fastq_files.tsv` | Every recovered input FASTQ location, storage/alias status, SHA-256 of stored and decompressed bytes, read counts, length ranges, and first/last read headers. |
| `fastq_pairs.tsv` | Mate counts and ordered read-ID checks in each original directory layout. |
| `pipeline_inputs.tsv` | One row per upstream sample: selected run, exact log command and line number, recovered TACTIC and DADA2 mates, and logged spike-removal counts. |
| `sample_input_links.tsv` | Links the archived expected/inferred sample names to mates for all five methods. |
| `provenance_evidence.tsv` | Exact commands and input/read-count evidence with archive paths and log line numbers. |
| `zero_length_reads.tsv` | Retained empty mates, which can explain differences between physical records and a tool's reported paired subset. |
| `limitations_and_findings.tsv` | Concrete differences, provenance qualifications, and requirements for an upstream replay. |
| `benchmark_preparation_validation.tsv` | Independently regenerates only benchmark preparation in a temporary output and compares it to archived tables. |
| `validation.json` | Input inventory, FASTQ checks, and any errors detected while collecting files. |
| `SHA256SUMS` | Hashes of every ZIP member except this checksum list itself. |

Archive paths are relative to `data/` after extraction. Original paths are evidence
from the analysis machine and need not exist on another computer. Byte-identical
FASTQ copies are stored once; alias rows resolve to that stored archive path. This
includes DADA2 work/trimmed copies and renamed raw downloads. The manifest makes
the usable archive path explicit instead of relying on obsolete source MANIFESTs.
The table `files.tsv` lists physical payloads; `fastq_files.tsv` retains all input
locations including aliases. Filtered DADA2 reads and disposable analysis
intermediates are not required benchmark inputs and are omitted.

## Special handling and reproducibility limits

For mocks 18 and 19, the selected TACTIC runs received reads containing spikes,
then removed spike pairs internally. The original pre-removal trimmed reads were
subsequently moved/renamed under `trimmed_but_spiked/`. DADA2's retained input
copies are the later `spikefree` reads. Both are archived. The per-sample log
comparison checks recovered pre-removal read counts against the TACTIC logs;
it does not assume that the currently named spikefree file was the original
TACTIC input. Consult `pipeline_inputs.tsv` before replaying either pipeline.

For SRR13005964, both retained FASTQs contain 147918 records. Its historical
Bowtie2 output reports 147918 input records, comprising 147917 paired records
and one unpaired read. R2 of `SRR13005964.10542` is empty. TACTIC's shorter
summary reports the paired subset, explaining the apparent one-pair difference.
These input records are preserved; the manifests distinguish physical records,
reported paired/unpaired counts, and the subset remaining after spike removal.
Read-ID checks normalize SRA's `.1`/`.2` mate suffix only when it follows a
separate spot number (for example `SRR2990088.1.1` and `SRR2990088.1.2`), and
also handle `/1` and `/2` suffixes. Spot numbers are retained.

The TIC input files for mocks 4 and 5 incorporate the archaeal corrections made
after their original upstream runs. The exact corrected feature tables, FASTAs,
and membership maps are retained. They permit reproducing the benchmark; a fresh
unmodified upstream TACTIC run can differ from them. Available before-correction
tables/maps are retained as evidence, but no complete automated correction recipe
or immutable upstream container image is present.

Historical logs recorded filenames and read counts but no original FASTQ hashes.
The supplied hashes freeze the recovered files as packaged now. Agreement of
names, counts, pair identities, and DADA2 copies supports their provenance, but
does not replace an unavailable historical checksum. Locally available SILVA
inputs are included; container-internal spike reference sequences, the original
Docker image digest, and original binary images are not archived. Therefore the
bundle supports benchmark reproduction from the exact inferred tables, while an
end-to-end regeneration of those tables from reads still requires the recorded
upstream software and the documented corrections. No claim of rerunning or
verifying all figures from raw reads is made.
Generated DADA2 R scripts and pair lists are retained, but their exact historical
positional argument values and immutable container digest were not recorded.
Current driver scripts are evidence and can contain early exits; they are not
presented as a verified replay of every original upstream invocation.

The build and validation scripts are in `input_manifest/scripts/` and in the
bundle's `scripts/` directory. The archive uses ZIP64 and retains gzip files
without recompression. `data/tacticbench_inputs_with_fastq.zip.sha256` checks
the outer archive; `data/input_archive_validation.json` records the full ZIP
CRC and member SHA-256 verification performed after writing it.
"""
    (reports / "README.md").write_text(readme)
    if summary["errors"]:
        print(json.dumps(summary, indent=2))
        raise ValueError("Input audit found errors; reports were retained for inspection")
    print(json.dumps(summary, indent=2), flush=True)
    if args.inventory_only:
        return
    for path in reports.rglob("*"):
        if path.is_file() and path.name not in {"files.tsv", "SHA256SUMS"}:
            add(path, "input_manifest/" + path.relative_to(reports).as_posix(), "generated_manifest_or_provenance")
    physical = []
    for index, (name, info) in enumerate(sorted(payload.items()), 1):
        path = info["source"].resolve()
        physical.append({"archive_path": name, "source_path": str(path), "community_id": info["community_id"],
                         "role": info["role"], "size_bytes": path.stat().st_size,
                         "source_mtime_utc": datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat(),
                         "sha256": cached_hash(path)})
    write_tsv(reports / "files.tsv", physical)
    add(reports / "files.tsv", "input_manifest/files.tsv", "generated_payload_manifest")
    sums = []
    for name, info in sorted(payload.items()):
        sums.append(cached_hash(info["source"]) + "  " + name + "\n")
    (reports / "SHA256SUMS").write_text("".join(sums))
    add(reports / "SHA256SUMS", "input_manifest/SHA256SUMS", "checksum_manifest")
    args.cache.write_text(json.dumps(cache))
    target = output / "tacticbench_inputs_with_fastq.zip"
    partial = target.with_suffix(".zip.partial")
    print(f"Writing {len(payload)} members to {target}", flush=True)
    with zipfile.ZipFile(partial, "w", allowZip64=True, compression=zipfile.ZIP_DEFLATED, compresslevel=1) as archive:
        for index, (name, info) in enumerate(sorted(payload.items()), 1):
            source = info["source"]
            compress = zipfile.ZIP_STORED if source.suffix.lower() == ".gz" else zipfile.ZIP_DEFLATED
            archive.write(source, name, compress_type=compress, compresslevel=1)
            if index % 150 == 0:
                print(f"Archived {index}/{len(payload)} files", flush=True)
    checked = verify(partial)
    partial.replace(target)
    checked["archive"] = str(target)
    checked["archive_size_bytes"] = target.stat().st_size
    checked["archive_sha256"] = sha256(target)
    (output / (target.name + ".sha256")).write_text(checked["archive_sha256"] + "  " + target.name + "\n")
    (output / "input_archive_validation.json").write_text(json.dumps(checked, indent=2) + "\n")
    print(json.dumps(checked, indent=2), flush=True)


if __name__ == "__main__":
    main()
