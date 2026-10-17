#!/usr/bin/env python3
"""
Comprehensive read length survey for ZIEL_1 fastq files.
Analyzes ALL reads (not sampled) to get full length distribution.
"""

import gzip
from pathlib import Path
from collections import Counter
import statistics

def analyze_fastq_lengths(fastq_file):
    """Extract ALL read lengths from a fastq file."""
    lengths = []
    
    try:
        if str(fastq_file).endswith('.gz'):
            f = gzip.open(fastq_file, 'rt')
        else:
            f = open(fastq_file, 'r')
        
        with f:
            line_num = 0
            for line in f:
                line = line.strip()
                line_num += 1
                
                # FASTQ format: header, sequence, +, quality
                if line_num % 4 == 2:  # Sequence line
                    lengths.append(len(line))
    except Exception as e:
        print(f"Error reading {fastq_file}: {e}")
        return None
    
    return lengths

def main():
    fastq_dir = Path("/media/adminuser/MSD-Storage/Mohsen/TACTIC_Benchmark/Mocks/Picked DO NOT ENTER ME DO NOT ENTER ME DO NOT ENTER ME/Isabel/ZIEL_1/fastqs")
    
    # Find all raw fastq files
    fastq_files = sorted([f for f in fastq_dir.glob('*.fastq.gz') if '__processed' not in str(f)])
    
    print("\n" + "="*80)
    print("COMPREHENSIVE READ LENGTH SURVEY - ZIEL_1 FASTQ FILES")
    print("="*80 + "\n")
    
    all_results = {}
    
    for fastq_file in fastq_files:
        print(f"\nAnalyzing: {fastq_file.name}")
        print("-" * 80)
        
        lengths = analyze_fastq_lengths(fastq_file)
        
        if not lengths:
            print("  Failed to read file")
            continue
        
        # Calculate statistics
        length_counter = Counter(lengths)
        unique_lengths = sorted(length_counter.keys())
        
        stats = {
            'total_reads': len(lengths),
            'unique_lengths': unique_lengths,
            'length_distribution': dict(length_counter),
            'min': min(lengths),
            'max': max(lengths),
            'mean': statistics.mean(lengths),
            'median': statistics.median(lengths),
            'stdev': statistics.stdev(lengths) if len(lengths) > 1 else 0,
        }
        
        all_results[fastq_file.name] = stats
        
        print(f"  Total reads: {stats['total_reads']:,}")
        print(f"  Unique lengths: {unique_lengths}")
        
        if len(unique_lengths) == 1:
            print(f"  ✓ All reads are {unique_lengths[0]} bp")
        else:
            print(f"  ✗ Found {len(unique_lengths)} different lengths:")
            for length in unique_lengths:
                count = length_counter[length]
                percentage = (count / len(lengths)) * 100
                print(f"      {length} bp: {count:,} reads ({percentage:.4f}%)")
        
        print(f"  Min length: {stats['min']} bp")
        print(f"  Max length: {stats['max']} bp")
        print(f"  Mean length: {stats['mean']:.2f} bp")
        print(f"  Median length: {stats['median']} bp")
        print(f"  Std Dev: {stats['stdev']:.4f}")
    
    # Summary
    print("\n" + "="*80)
    print("SUMMARY - ZIEL_1 FASTQ FILES")
    print("="*80 + "\n")
    
    print("File-by-file summary:")
    for fname in sorted(all_results.keys()):
        stats = all_results[fname]
        unique_lengths = stats['unique_lengths']
        if len(unique_lengths) == 1:
            print(f"  {fname:<30} {stats['total_reads']:>10,} reads  │  {unique_lengths[0]} bp")
        else:
            print(f"  {fname:<30} {stats['total_reads']:>10,} reads  │  {stats['min']}-{stats['max']} bp")
    
    # Check consistency across R1/R2 pairs
    print("\n" + "-"*80)
    print("R1/R2 Pairing Consistency Check:")
    print("-"*80)
    
    samples = {}
    for fname in all_results.keys():
        # Extract sample ID and direction
        parts = fname.replace('.fastq.gz', '').split('_')
        sample_id = '_'.join(parts[:-1])
        direction = parts[-1]
        
        if sample_id not in samples:
            samples[sample_id] = {}
        samples[sample_id][direction] = all_results[fname]
    
    for sample_id in sorted(samples.keys()):
        print(f"\n  {sample_id}:")
        r1_stats = samples[sample_id].get('R1')
        r2_stats = samples[sample_id].get('R2')
        
        if r1_stats and r2_stats:
            r1_lengths = r1_stats['unique_lengths']
            r2_lengths = r2_stats['unique_lengths']
            
            if r1_lengths == r2_lengths:
                print(f"    ✓ R1 and R2 match: {r1_lengths[0]} bp")
            else:
                print(f"    ✗ R1 and R2 differ:")
                print(f"      R1: {r1_lengths}")
                print(f"      R2: {r2_lengths}")
    
    # Overall consistency
    print("\n" + "-"*80)
    print("Overall Consistency Across All Files:")
    print("-"*80)
    
    all_unique_lengths = set()
    for fname, stats in all_results.items():
        all_unique_lengths.update(stats['unique_lengths'])
    
    all_unique_lengths = sorted(all_unique_lengths)
    
    if len(all_unique_lengths) == 1:
        print(f"  ✓ FULLY CONSISTENT: All files contain {all_unique_lengths[0]} bp reads")
    else:
        print(f"  ✗ VARIABLE: Files contain mixed lengths: {all_unique_lengths}")

if __name__ == '__main__':
    main()
