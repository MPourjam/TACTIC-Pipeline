"""
Docstring for Step-wise_Diversity.pcr_simulation.seqkit
This module contains function wrapper around seqkit functionalities.
"""
import sys
import subprocess
import logging
from pathlib import Path
from typing import List

# Configure logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.propagate = True


def check_md5(file_path: str, expected_md5: str) -> bool:
    """Check if the MD5 checksum of the file matches the expected MD5."""
    import hashlib

    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    file_md5 = hash_md5.hexdigest()
    return file_md5 == expected_md5


def my_subprocess_run(
    cmd_list: List[str],
    *args,
    **kwargs) -> subprocess.CompletedProcess:
    """Run a subprocess command and return the CompletedProcess object."""
    logger.info(f"Running command: {' '.join(cmd_list)}")
    try:
        result = subprocess.run(
            cmd_list,
            *args,
            **kwargs
        )
        logger.info(f"Command finished with return code: {result.returncode}")
        # logger.debug(f"Command stdout: {result.stdout}")
        logger.debug(f"Command stderr: {result.stderr}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running command {' '.join(cmd_list)}: {e.stderr}")
        sys.exit(1)

    return result


def seqkit_common(
    input_fasta1: Path,
    input_fasta2: Path,
    output_fasta: Path,
    args_list: List[str] = []) -> None:
    """Common function to run seqkit commands."""
    ret_code = 1
    try:
        res_ = my_subprocess_run(
            [
            "seqkit",
            "common",
            *args_list,
            "-o", str(output_fasta),
            str(input_fasta1),
            str(input_fasta2)
            ],
            check=True
        )
        ret_code = res_.returncode

    except subprocess.CalledProcessError as e:
        logger.error(f"Error running seqkit common: {e.stderr}")
        sys.exit(1)
    return ret_code, output_fasta


def seqkit_stats(
    fasta_file: Path,
    args_list: List[str] = None,
    output_file: Path = None) -> tuple:
    """Run seqkit stats on the given FASTA file and return the output lines."""
    if args_list is None:
        args_list = []

    cmd_list = [
        "seqkit",
        "stats",
        *args_list
    ]

    if output_file is not None:
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        cmd_list.extend(
            ["-o", str(output_file)]
        )
    
    cmd_list.append(str(fasta_file))

    try:
        result = my_subprocess_run(
            cmd_list,
            capture_output=True,
            text=True,
            check=True
        )
        
        if output_file is None:
            content = result.stdout.strip()
        else:
            with open(output_file, 'r') as of:
                content = of.read().strip()

        output_lines = content.split('\n')
        output_dict = {}
        
        if len(output_lines) >= 2:
            headers = output_lines[0].split()
            values = output_lines[1].split()
            output_dict = dict(zip(headers, values))

        return result.returncode, output_file, output_dict

    except subprocess.CalledProcessError as e:
        logger.error(f"Error running seqkit stats: {e.stderr}")
        sys.exit(1)


def seqkit_subseq(
    input_fasta: Path,
    output_fasta: Path,
    args_list: List[str] = [],
    output_md5sum: str = None) -> None:
    """Cut sequences in the FASTA file from start to end positions using seqkit."""
    ret_code = 1
    cmd_list = [
        "seqkit", "subseq",
        *args_list,
        "-o", str(output_fasta),
        str(input_fasta)
    ]
    if output_md5sum and output_fasta.exists():
        # check the md5sum before running
        if check_md5(str(output_fasta), output_md5sum):
            logger.info(f"Output FASTA already exists with matching MD5: {output_fasta}")
            return 0, output_fasta
    
    logger.info(f"Running command: {' '.join(cmd_list)}")

    try:
        res_ = my_subprocess_run(
            cmd_list,
            check=True
        )
        ret_code = res_.returncode
    except subprocess.CalledProcessError as e:
        logger.error(f"Error cutting sequences with seqkit: {e.stderr}")
        sys.exit(1)
    
    if output_md5sum:
        if not check_md5(str(output_fasta), output_md5sum):
            logger.error(f"MD5 checksum does not match for output FASTA: {output_fasta}")
            sys.exit(1)

    return ret_code, output_fasta


def seqkit_seq(
    input_fasta: Path,
    output_fasta: Path,
    args_list: List[str] = []) -> None:
    """Seq command of seqkit."""
    ret_code = 1
    try:
        res_ = my_subprocess_run(
            [
            "seqkit", "seq",
            *args_list,
            "-o", str(output_fasta),
            str(input_fasta)
            ],
            check=True
        )
        ret_code = res_.returncode
    except subprocess.CalledProcessError as e:
        logger.error(f"Error filtering sequences with seqkit: {e.stderr}")
        sys.exit(1)
    
    return ret_code, output_fasta


def seqkit_watch(
    input_fasta: Path,
    output_fasta: Path,
    args_list: List[str] = []) -> None:
    """Watch command of seqkit."""
    ret_code = 1
    try:
        res_ = my_subprocess_run(
            [
            "seqkit", "watch",
            *args_list,
            "-o", str(output_fasta),
            str(input_fasta)
            ],
            check=True
        )
        ret_code = res_.returncode
    except subprocess.CalledProcessError as e:
        logger.error(f"Error watching sequences with seqkit: {e.stderr}")
        sys.exit(1)
    
    return ret_code, output_fasta


def seqkit_amplicon(
    input_fasta: Path,
    output_fasta: Path,
    args_list: List[str] = []) -> None:
    """Amplicon command of seqkit."""
    ret_code = 1
    try:
        res_ = my_subprocess_run(
            [
            "seqkit", "amplicon",
            *args_list,
            "-o", str(output_fasta),
            str(input_fasta)
            ],
            check=True
        )
        ret_code = res_.returncode
    except subprocess.CalledProcessError as e:
        logger.error(f"Error generating amplicons with seqkit: {e.stderr}")
        sys.exit(1)
    
    return ret_code, output_fasta


def seqkit_rmdup(
    input_fasta: Path,
    output_fasta: Path,
    args_list: List[str] = []) -> None:
    """Remove duplicate sequences using seqkit rmdup."""
    ret_code = 1
    try:
        res_ = my_subprocess_run(
            [
            "seqkit", "rmdup",
            *args_list,
            "-o", str(output_fasta),
            str(input_fasta)
            ],
            check=True
        )
        ret_code = res_.returncode
    except subprocess.CalledProcessError as e:
        logger.error(f"Error removing duplicates with seqkit: {e.stderr}")
        sys.exit(1)
    
    return ret_code, output_fasta
