#!/usr/bin/env python3
"""
Simple in silico amplicon generation for 16S rRNA using cutadapt.
"""

import argparse
import gzip
import shutil
import subprocess
import csv
import re
import sys
from collections import OrderedDict
from itertools import product
from pathlib import Path
from typing import Dict, List, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed


# Default primers definition
DEFAULT_PRIMERS = {
    "forwards": {
        "27F": "AGAGTTTGATCMTGGCTCAG"
    },
    "reverses": {
        "338R": "TGCTGCCTCCCGTAGGAGT",
        "519R": "GWATTACCGCGGCKGCTG",
        "785R": "GACTACHVGGGTATCTAATCC",
        "926R": "CCGTCAATTCMTTTRAGT",
        "1492R": "GGTTACCTTGTTACGACTT"
    }
}

# IUPAC-aware complement table
IUPAC_COMPLEMENT = {
    "A": "T", "C": "G", "G": "C", "T": "A",
    "R": "Y", "Y": "R", "S": "S", "W": "W",
    "K": "M", "M": "K", "B": "V", "D": "H",
    "H": "D", "V": "B", "N": "N",
    "a": "t", "c": "g", "g": "c", "t": "a",
    "r": "y", "y": "r", "s": "s", "w": "w",
    "k": "m", "m": "k", "b": "v", "d": "h",
    "h": "d", "v": "b", "n": "n",
}


def rc(seq: str) -> str:
    """Reverse-complement, preserving IUPAC degenerates."""
    return "".join(IUPAC_COMPLEMENT.get(b, "N") for b in reversed(seq))


def gunzip_if_needed(path: Path) -> Path:
    """If path ends with .gz, decompress next to it and return new path."""
    if path.suffix == ".gz":
        out = path.with_suffix("")
        with gzip.open(path, "rb") as src, open(out, "wb") as dst:
            shutil.copyfileobj(src, dst)
        return out
    return path

def load_primer_table(primer_table_path):
    """Load primer table TSV into a dict."""
    primer_table_path = Path(primer_table_path)
    if not primer_table_path.is_file():
        raise FileNotFoundError(f"Primer table not found: {primer_table_path}")

    primers = {
        "forwards": {},
        "reverses": {}
    }
    with primer_table_path.open() as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            name = row['name'].strip()
            direction = row.get('direction', '').strip()
            direction = 'forwards' if "F" == direction else 'reverses' if "R" == direction else None
            sequence = row.get('sequence', '').strip()
            if not name or not sequence:
                continue
            primers[direction][name] = sequence

    return primers


def primer_position_from_name(name):
    """Extract numeric position from primer name (e.g., '27F_Isi' -> 27)."""
    m = re.search(r'(\d+)', name)
    if not m:
        raise ValueError(f"Could not extract numeric position from primer name: {name}")
    return int(m.group(1))


def generate_forw_reve_set(primers: Dict[str, Dict[str, str]]) -> List[Tuple[Tuple[str, str], Tuple[str, str]]]:
    forwards = list(primers.get("forwards", {}).keys())
    reverses = list(primers.get("reverses", {}).keys())
    if not forwards or not reverses:
        return []
    # make sure that always reverse primers are matched to their preceding forward primers
    paired_primers = []
    for fwd in forwards:
        fwd_pos = primer_position_from_name(fwd)
        fwd_seq = primers["forwards"][fwd]
        for rev in reverses:
            rev_pos = primer_position_from_name(rev)
            rev_seq = rc(primers["reverses"][rev])
            if rev_pos > fwd_pos:
                paired_primers.append(((fwd, fwd_seq), (rev, rev_seq)))

    return paired_primers


def get_db_name_from_path(ref_fasta_path):
    """Extract database name (SILVA or GG2) from the reference FASTA path."""
    path_str = str(ref_fasta_path)
    if 'SILVA' in path_str:
        return 'SILVA'
    elif 'GG2' in path_str:
        return 'GG2'
    else:
        stem = Path(ref_fasta_path).stem
        # Handle double extensions like .fasta.gz
        if stem.endswith('.fasta') or stem.endswith('.fa'):
            stem = Path(stem).stem
        return stem


def run_cutadapt(
    ref_fasta: str,
    out_fasta: str,
    fwd_seq: str,
    rev_seq: str,
    error_rate: float,
    min_overlap: int,
    minimum_length: int,
    maximum_length: int,
    action: str = "trim",
    threads: int = 1,
    cutadapt_bin: str = "cutadapt"
) -> Path:
    """Run cutadapt for a given primer pair."""
    out_fasta = Path(out_fasta).absolute()
    ref_fasta = Path(ref_fasta).absolute()
    untrimmed_fasta = out_fasta.with_suffix('.untrimmed.fasta')
    cutadapt_json = out_fasta.with_suffix('.cutadapt.json')
    trimmed_fasta = out_fasta.with_suffix('.trimmed.fasta')
    out_fasta_dir = out_fasta.parent
    out_fasta_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        cutadapt_bin,
        "-j", str(threads),
        # "-g", fwd_seq,  # Guiding adapter 5'
        # Linked adapters ref: https://cutadapt.readthedocs.io/en/stable/guide.html#linked-adapters
        "-a", f"{fwd_seq}...{rev_seq}",
        "--action", action,
        '--match-read-wildcards',  # Interpret IUPAC codes in reads
        # "--discard-untrimmed",  # Does not work with --untrimmed-output
        "--overlap", str(min_overlap),
        "--error-rate", str(error_rate),
        "--minimum-length", str(minimum_length),
        "--maximum-length", str(maximum_length),
        "--json", str(cutadapt_json),
        "--report", "full",
        "--untrimmed-output", str(untrimmed_fasta),
        "-o", str(trimmed_fasta),
        str(ref_fasta)
    ]
    # write the command to a file
    with open(out_fasta.with_suffix('.cutadapt.cmd').with_suffix('.txt'), 'w') as f:
        f.write(' '.join(cmd) + '\n')

    subprocess.run(cmd, check=True, capture_output=True, text=True)
    return str(trimmed_fasta), str(cutadapt_json), str(untrimmed_fasta)


def calculate_n50(lengths):
    """Calculates N50 metric."""
    if not lengths:
        return 0
    
    sorted_lengths = sorted(lengths, reverse=True)
    total_len = sum(sorted_lengths)
    target = total_len / 2
    
    current_sum = 0
    for l in sorted_lengths:
        current_sum += l
        if current_sum >= target:
            return l
    return 0


def cutadapt_parameter_sweep(
        ref_fastas: list = ["./DBs/SILVA/SILVA_138.2_SSURef_NR99_tax_silva.fasta"],
        error_rates: list = [0.1],
        min_overlaps: list = [3],
        min_amp_lens: list = [100],
        max_amp_lens: list = [2000],
        primers_tsv: str = "",
        cutadapt_threads: int = 1,
        out_dir: str = "./results",
        n_jobs: int = 1):
    """
    This function creates combinations of parameters from the given main parameters for 
    cutadapt and runs run_cutadapt() with the given prameters and stores the resutls in 
    relevant dirs.
    """
    # if primers_tsv is not given then load the default primers
    primers_tsv_path = Path(primers_tsv) if primers_tsv else None
    primers_d = DEFAULT_PRIMERS

    if primers_tsv_path and primers_tsv_path.is_file():
        primers_d = load_primer_table(primers_tsv_path)
    primers_set = generate_forw_reve_set(primers_d)

    # output_dir
    # make combinations main arguments values
    params_combinations = []
    for ref, er, ov, min_l, max_l, pri_set in product(
        ref_fastas,
        error_rates,
        min_overlaps,
        min_amp_lens,
        max_amp_lens,
        primers_set
    ):
        run_id = f"er{er}_ov{ov}_min{min_l}_max{max_l}_{pri_set[0][0]}_{pri_set[1][0]}"
        db_type = get_db_name_from_path(ref)
        this_out_dir = Path(f"{out_dir}/{db_type}/{run_id}")
        this_out_dir.mkdir(parents=True, exist_ok=True)
        output_fasta = this_out_dir / f"{run_id}.fasta"

        param_dict = {
            'ref_fasta': ref,
            'out_fasta': output_fasta,
            'error_rate': er,
            'min_overlap': ov,
            'minimum_length': min_l,
            'maximum_length': max_l,
            'fwd_seq': pri_set[0][1],
            'rev_seq': pri_set[1][1],
            'threads': cutadapt_threads,
        }

        params_combinations.append((run_id, param_dict, (pri_set[0][0], pri_set[1][0])))

    # Run cutadapt parallel
    results = []
    print(f"Running cutadapt for {len(params_combinations)} parameter combinations")
    with ProcessPoolExecutor(max_workers=n_jobs) as executor:
        futures = {}
        for run_id, param_dict, primer_names in params_combinations:
            # Merge other_cutadapt_args with param_dict
            cutadapt_args = {**param_dict}  #, **other_cutadap_args}
            future = executor.submit(run_cutadapt, **cutadapt_args)
            futures[future] = (run_id, param_dict, primer_names)
        
        for future in as_completed(futures):
            run_id, param_dict, primer_names = futures[future]
            try:
                result = future.result()
                results.append((run_id, param_dict, primer_names, result))
                print(f"Completed: {run_id}")
            except Exception as e:
                print(f"Error in {run_id}: {e}", file=sys.stderr)
                results.append((run_id, param_dict, primer_names, None))
    
    # result: (trimmed_fasta, cutadapt_json, untrimmed_fasta)
    return results

