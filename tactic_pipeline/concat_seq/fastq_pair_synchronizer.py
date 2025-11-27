import gzip
import os
import argparse
import sys
from typing import Dict, Tuple, Iterator, List, TextIO

# Define a type for a single FASTQ record (Header, Sequence, Delimiter, Quality)
FastqRecord = Tuple[str, str, str, str]

def clean_header(header: str) -> str:
    """
    Cleans a FASTQ header to get the unique read ID by removing the leading '@'
    and any paired-end suffixes (e.g., /1, /2, :1:1, :2:1).
    """
    header = header.strip().split(" ")[0].lstrip('@')
    # Common paired-end suffixes to remove for matching:
    for suffix in ['/1', '/2', ' 1:N', ' 2:N', ' 1:Y', ' 2:Y']:
        if header.endswith(suffix):
            header = header[:-len(suffix)].strip()
    return header

def read_fastq_record_generator(file_handle: TextIO) -> Iterator[FastqRecord]:
    """
    Generator to read a FASTQ file four lines at a time.
    Yields (header, sequence, delimiter, quality_score).
    """
    while True:
        try:
            # Read and strip four lines for a single record
            header = next(file_handle).strip()
            sequence = next(file_handle).strip()
            delimiter = next(file_handle).strip() 
            quality_score = next(file_handle).strip()
            
            if not header.startswith('@'):
                 # Handle unexpected file end or formatting error
                 continue
                 
            yield header, sequence, delimiter, quality_score
        except StopIteration:
            break

def index_fastq_file(r_file_path: str) -> Dict[str, FastqRecord]:
    """
    Reads a FASTQ file entirely into a dictionary, keyed by the cleaned read ID.
    This creates the lookup table for the matching process.
    """
    print(f"Indexing reads from file: {r_file_path}")
    
    # Use appropriate open function for gzipped or plain text files
    _open = gzip.open if r_file_path.endswith('.gz') else open
    index = {}
    
    try:
        with _open(r_file_path, 'rt') as r_file:
            for i, record in enumerate(read_fastq_record_generator(r_file)):
                header = record[0]
                unique_id = clean_header(header)
                index[unique_id] = record
                
                if (i + 1) % 100000 == 0:
                    print(f"  Indexed {i + 1} reads...", end='\r')
                    
        print(f"\nIndex complete. Found {len(index)} unique read IDs.")
        return index
    except FileNotFoundError:
        print(f"Error: File not found at {r_file_path}")
        return {}


def match_and_pair_fastq(
    r1_file_path: str, 
    r2_file_path: str, 
    output_r1_path: str, 
    output_r2_path: str
):
    """
    Synchronizes two paired-end FASTQ files based on read IDs.
    It discards any reads that do not have a partner in the corresponding file.

    Args:
        r1_file_path: Path to the R1 (Forward) FASTQ file.
        r2_file_path: Path to the R2 (Reverse) FASTQ file.
        output_r1_path: Path for the new, synchronized R1 output file.
        output_r2_path: Path for the new, synchronized R2 output file.
    """
    
    # 1. Index one file (R1 is chosen arbitrarily)
    r1_index = index_fastq_file(r1_file_path)
    
    if not r1_index:
        print("Cannot proceed without a valid R1 index.")
        return

    # 2. Iterate through the second file (R2) and match
    print(f"Starting synchronization using R2 file: {r2_file_path}")

    _open = gzip.open if r2_file_path.endswith('.gz') else open
    
    processed_count = 0
    paired_count = 0
    
    try:
        with _open(r2_file_path, 'rt') as r2_file, \
             open(output_r1_path, 'w') as out_r1, \
             open(output_r2_path, 'w') as out_r2:
            
            r2_generator = read_fastq_record_generator(r2_file)
            
            for r2_record in r2_generator:
                r2_header, r2_seq, r2_delim, r2_qual = r2_record
                r2_unique_id = clean_header(r2_header)
                
                # Look up the partner read in the R1 index
                if r2_unique_id in r1_index:
                    # Match found! Retrieve R1 record and remove from index
                    r1_record = r1_index.pop(r2_unique_id)
                    
                    # Write synchronized R1 record
                    out_r1.write(f"{r1_record[0]}\n{r1_record[1]}\n{r1_record[2]}\n{r1_record[3]}\n")
                    
                    # Write synchronized R2 record
                    out_r2.write(f"{r2_header}\n{r2_seq}\n{r2_delim}\n{r2_qual}\n")
                    
                    paired_count += 1
                
                processed_count += 1
                
                if processed_count % 100000 == 0:
                    print(f"  Processed {processed_count} reads. Paired: {paired_count}...", end='\r')

            # 3. Report results
            orphaned_r1_count = len(r1_index)
            
            print("\n--- Synchronization Summary ---")
            print(f"Total R2 Reads Processed: {processed_count}")
            print(f"Total Pairs Written to Output: {paired_count}")
            print(f"Orphaned R1 Reads (discarded): {orphaned_r1_count}")
            print(f"Orphaned R2 Reads (discarded): {processed_count - paired_count}")
            print(f"Output files created: {output_r1_path}, {output_r2_path}")

    except Exception as e:
        print(f"An unexpected error occurred during processing: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Synchronize paired-end FASTQ files by discarding orphan reads."
    )
    
    parser.add_argument(
        "--r1", required=True, help="Path to the input R1 (Forward) FASTQ file"
    )
    parser.add_argument(
        "--r2", required=True, help="Path to the input R2 (Reverse) FASTQ file"
    )
    parser.add_argument(
        "--out-r1", required=True, help="Path for the synchronized R1 output file"
    )
    parser.add_argument(
        "--out-r2", required=True, help="Path for the synchronized R2 output file"
    )

    args = parser.parse_args()

    # Check input files exist
    if not os.path.exists(args.r1):
        print(f"Error: R1 file not found at {args.r1}")
        sys.exit(1)
    if not os.path.exists(args.r2):
        print(f"Error: R2 file not found at {args.r2}")
        sys.exit(1)

    match_and_pair_fastq(
        r1_file_path=args.r1,
        r2_file_path=args.r2,
        output_r1_path=args.out_r1,
        output_r2_path=args.out_r2
    )

if __name__ == '__main__':
    main()