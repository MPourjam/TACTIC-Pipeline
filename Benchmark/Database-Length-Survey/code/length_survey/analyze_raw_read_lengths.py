#!/usr/bin/env python3
"""
Analyze read lengths in raw fastq files for ZIEL-1, ZIEL-2, and ZYMO mock communities.
"""

import gzip
import os
from pathlib import Path
from collections import Counter, defaultdict
import csv

def get_read_lengths(fastq_file, max_reads=None):
    """Extract read lengths from a fastq file (gzip-compressed or not)."""
    lengths = []
    read_count = 0
    
    try:
        if fastq_file.endswith('.gz'):
            f = gzip.open(fastq_file, 'rt')
        else:
            f = open(fastq_file, 'r')
        
        with f:
            line_num = 0
            for line in f:
                line = line.strip()
                line_num += 1
                
                # FASTQ format: header, sequence, +, quality
                # We want every 2nd line starting from line 2 (the sequence)
                if line_num % 4 == 2:  # Sequence line
                    lengths.append(len(line))
                    read_count += 1
                    if max_reads and read_count >= max_reads:
                        break
    except Exception as e:
        print(f"Error reading {fastq_file}: {e}")
        return None
    
    return lengths

def analyze_mock_community(mock_name, fastq_dir):
    """Analyze all fastq files for a mock community."""
    print(f"\n{'='*80}")
    print(f"Analyzing: {mock_name}")
    print(f"Directory: {fastq_dir}")
    print(f"{'='*80}\n")
    
    # Find all R1 and R2 fastq files
    fastq_files = sorted(Path(fastq_dir).glob('*.fastq.gz'))
    # Filter out processed files
    fastq_files = [f for f in fastq_files if '__processed' not in str(f)]
    
    if not fastq_files:
        print(f"No fastq files found in {fastq_dir}")
        print(f"Attempting to list directory contents:")
        try:
            for item in os.listdir(fastq_dir):
                print(f"  {item}")
        except Exception as e:
            print(f"  Error: {e}")
        return None
    
    results = {}
    
    for fastq_file in fastq_files:
        print(f"\nFile: {fastq_file.name}")
        lengths = get_read_lengths(str(fastq_file), max_reads=10000)
        
        if not lengths:
            print("  Failed to read file")
            continue
        
        # Analyze lengths
        length_counter = Counter(lengths)
        unique_lengths = sorted(length_counter.keys())
        
        print(f"  Total reads analyzed: {len(lengths)}")
        print(f"  Unique lengths: {unique_lengths}")
        print(f"  Length distribution:")
        for length in unique_lengths:
            count = length_counter[length]
            percentage = (count / len(lengths)) * 100
            print(f"    {length} bp: {count} reads ({percentage:.2f}%)")
        
        # Check consistency
        if len(unique_lengths) == 1:
            print(f"  ✓ CONSISTENT: All reads are {unique_lengths[0]} bp")
        else:
            print(f"  ✗ INCONSISTENT: Found {len(unique_lengths)} different lengths")
        
        results[fastq_file.name] = {
            'total_reads': len(lengths),
            'unique_lengths': unique_lengths,
            'length_distribution': dict(length_counter),
            'is_consistent': len(unique_lengths) == 1
        }
    
    return results

def main():
    base_path = "/media/adminuser/MSD-Storage/Mohsen/TACTIC_Benchmark/Mocks/Picked DO NOT ENTER ME DO NOT ENTER ME DO NOT ENTER ME/Isabel"
    
    mocks = {
        'ZIEL_1': 'ZIEL_1',
        'ZIEL_2': 'ZIEL_2',
        'ZYMO': 'ZYMO'
    }
    
    all_results = {}
    
    for mock_name, dir_name in mocks.items():
        fastq_dir = os.path.join(base_path, dir_name, 'fastqs')
        if os.path.isdir(fastq_dir):
            results = analyze_mock_community(mock_name, fastq_dir)
            if results:
                all_results[mock_name] = results
    
    # Generate summary report
    print("\n\n" + "="*80)
    print("SUMMARY REPORT")
    print("="*80 + "\n")
    
    summary_file = "/media/adminuser/MSD-Storage/Mohsen/TACTIC_Benchmark/raw_read_length_analysis.csv"
    
    with open(summary_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Mock_Community', 'Fastq_File', 'Total_Reads_Sampled', 'Unique_Lengths', 
                        'Min_Length', 'Max_Length', 'Consistent', 'Length_Distribution'])
        
        for mock_name in sorted(all_results.keys()):
            print(f"\n{mock_name}:")
            print("-" * 40)
            
            for file_name in sorted(all_results[mock_name].keys()):
                file_info = all_results[mock_name][file_name]
                unique_lengths = file_info['unique_lengths']
                
                if unique_lengths:
                    min_len = min(unique_lengths)
                    max_len = max(unique_lengths)
                else:
                    min_len = max_len = 'N/A'
                
                is_consistent = file_info['is_consistent']
                dist = file_info['length_distribution']
                dist_str = '; '.join([f"{k}bp:{v}reads" for k, v in sorted(dist.items())])
                
                print(f"  {file_name}:")
                print(f"    Total reads: {file_info['total_reads']}")
                print(f"    Lengths: {unique_lengths}")
                print(f"    Consistent: {'Yes' if is_consistent else 'No'}")
                
                writer.writerow([
                    mock_name,
                    file_name,
                    file_info['total_reads'],
                    len(unique_lengths),
                    min_len,
                    max_len,
                    'Yes' if is_consistent else 'No',
                    dist_str
                ])
    
    print(f"\n✓ Summary saved to: {summary_file}")
    
    # Check for consistency across reads and within mock communities
    print("\n" + "="*80)
    print("CROSS-FILE CONSISTENCY CHECK")
    print("="*80 + "\n")
    
    for mock_name in sorted(all_results.keys()):
        print(f"\n{mock_name}:")
        all_lengths_in_mock = set()
        files_consistent = True
        
        for file_name, file_info in all_results[mock_name].items():
            unique_lengths = file_info['unique_lengths']
            all_lengths_in_mock.update(unique_lengths)
            
            if not file_info['is_consistent']:
                files_consistent = False
        
        if len(all_lengths_in_mock) == 1:
            print(f"  ✓ All files in {mock_name} have consistent read length: {list(all_lengths_in_mock)[0]} bp")
        else:
            print(f"  ✗ {mock_name} has variable read lengths: {sorted(all_lengths_in_mock)}")
        
        if not files_consistent:
            print(f"  ⚠ WARNING: Some files contain variable-length reads")

if __name__ == '__main__':
    main()
