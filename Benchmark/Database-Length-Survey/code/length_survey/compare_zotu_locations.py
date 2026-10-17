#!/usr/bin/env python3
"""
Comprehensive ZOTU sequence length survey and comparison.
Compares ZOTU sequences from:
1. /media/adminuser/MSD-Storage/Mohsen/TACTIC_Benchmark/Mocks/Benchmark_git/All_Analyses/mocks (tactic_results)
2. /media/adminuser/MSD-Storage/Mohsen/TACTIC_Benchmark/Mocks/Benchmark_git/All_Analyses/mocks/raw_data (trimmed_fastq/fastqs)
"""

import os
from pathlib import Path
from collections import Counter, defaultdict
import statistics
import csv

def parse_fasta_lengths(fasta_file):
    """Extract sequence lengths from FASTA file."""
    lengths = []
    try:
        with open(fasta_file, 'r') as f:
            current_seq = ""
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_seq:
                        lengths.append(len(current_seq))
                    current_seq = ""
                else:
                    current_seq += line
            if current_seq:
                lengths.append(len(current_seq))
    except Exception as e:
        print(f"Error reading {fasta_file}: {e}")
        return None
    
    return lengths if lengths else None

def calculate_stats(lengths):
    """Calculate statistics for sequence lengths."""
    if not lengths:
        return None
    
    lengths_sorted = sorted(lengths)
    n = len(lengths)
    
    stats = {
        'count': n,
        'min': min(lengths),
        'max': max(lengths),
        'mean': sum(lengths) / n,
        'median': statistics.median(lengths),
        'stdev': statistics.stdev(lengths) if n > 1 else 0,
        'q1': lengths_sorted[n // 4] if n > 3 else lengths_sorted[0],
        'q3': lengths_sorted[3 * n // 4] if n > 3 else lengths_sorted[-1],
    }
    return stats

def main():
    base_dir = Path("/media/adminuser/MSD-Storage/Mohsen/TACTIC_Benchmark/Mocks/Benchmark_git/All_Analyses/mocks")
    
    # Location 1: tactic_results (original analysis)
    tactic_results_dir = base_dir
    
    # Location 2: raw_data (additional analysis)
    raw_data_dir = base_dir / "raw_data"
    
    print("\n" + "="*100)
    print("COMPREHENSIVE ZOTU SEQUENCE LENGTH SURVEY AND COMPARISON")
    print("="*100)
    print(f"\nLocation 1 (tactic_results): {tactic_results_dir}")
    print(f"Location 2 (raw_data): {raw_data_dir}\n")
    
    # Collect ZOTU files from both locations
    zotu_files_tactic = sorted(tactic_results_dir.glob('*/tactic_results/*/ZOTUs-Seqs.fasta'))
    zotu_files_raw_data = sorted(raw_data_dir.glob('*/trimmed_fastq/*/ZOTUs-Seqs.fasta'))
    zotu_files_raw_data += sorted(raw_data_dir.glob('*/fastqs/*/ZOTUs-Seqs.fasta'))
    
    print(f"Found {len(zotu_files_tactic)} ZOTU files in tactic_results")
    print(f"Found {len(zotu_files_raw_data)} ZOTU files in raw_data\n")
    
    results_by_location = {
        'tactic_results': {},
        'raw_data': {}
    }
    
    # Analyze tactic_results
    print("="*100)
    print("ANALYZING TACTIC_RESULTS LOCATION")
    print("="*100 + "\n")
    
    for i, zotu_file in enumerate(zotu_files_tactic, 1):
        lengths = parse_fasta_lengths(zotu_file)
        if not lengths:
            continue
        
        stats = calculate_stats(lengths)
        
        # Extract mock and analysis info
        parts = zotu_file.parts
        mock_idx = list(parts).index(base_dir.name) + 1
        mock_name = parts[mock_idx]
        analysis_dir = parts[-2]
        
        key = f"{mock_name}/{analysis_dir}"
        results_by_location['tactic_results'][key] = stats
        
        if i <= 5:  # Print first 5
            print(f"[{i}] {key}: {stats['count']} ZOTUs, "
                  f"{stats['min']}-{stats['max']} bp, mean {stats['mean']:.1f} bp")
    
    print(f"\n... and {len(zotu_files_tactic) - 5} more files")
    
    # Analyze raw_data
    print("\n" + "="*100)
    print("ANALYZING RAW_DATA LOCATION")
    print("="*100 + "\n")
    
    for i, zotu_file in enumerate(zotu_files_raw_data, 1):
        lengths = parse_fasta_lengths(zotu_file)
        if not lengths:
            continue
        
        stats = calculate_stats(lengths)
        
        # Extract mock and analysis info
        parts = zotu_file.parts
        raw_data_idx = list(parts).index('raw_data')
        mock_name = parts[raw_data_idx + 1]
        analysis_dir = parts[-2]
        
        # Determine source (trimmed_fastq or fastqs)
        source = parts[raw_data_idx + 2]
        
        key = f"{mock_name}/{source}/{analysis_dir}"
        results_by_location['raw_data'][key] = stats
        
        if i <= 5:  # Print first 5
            print(f"[{i}] {key}: {stats['count']} ZOTUs, "
                  f"{stats['min']}-{stats['max']} bp, mean {stats['mean']:.1f} bp")
    
    print(f"\n... and {len(zotu_files_raw_data) - 5} more files")
    
    # Generate comparison reports
    print("\n" + "="*100)
    print("STATISTICAL COMPARISON")
    print("="*100 + "\n")
    
    # Aggregate statistics by location
    all_counts_tactic = [s['count'] for s in results_by_location['tactic_results'].values()]
    all_mins_tactic = [s['min'] for s in results_by_location['tactic_results'].values()]
    all_maxs_tactic = [s['max'] for s in results_by_location['tactic_results'].values()]
    all_means_tactic = [s['mean'] for s in results_by_location['tactic_results'].values()]
    
    all_counts_raw = [s['count'] for s in results_by_location['raw_data'].values()]
    all_mins_raw = [s['min'] for s in results_by_location['raw_data'].values()]
    all_maxs_raw = [s['max'] for s in results_by_location['raw_data'].values()]
    all_means_raw = [s['mean'] for s in results_by_location['raw_data'].values()]
    
    print("TACTIC_RESULTS Location:")
    print(f"  Files analyzed: {len(results_by_location['tactic_results'])}")
    print(f"  Total ZOTUs: {sum(all_counts_tactic):,}")
    print(f"  Mean ZOTUs per file: {statistics.mean(all_counts_tactic):.1f}")
    print(f"  Length range: {min(all_mins_tactic)}-{max(all_maxs_tactic)} bp")
    print(f"  Mean length across files: {statistics.mean(all_means_tactic):.2f} bp")
    print(f"  Std dev of mean lengths: {statistics.stdev(all_means_tactic):.2f} bp" if len(all_means_tactic) > 1 else "  N/A")
    
    print("\nRAW_DATA Location:")
    print(f"  Files analyzed: {len(results_by_location['raw_data'])}")
    print(f"  Total ZOTUs: {sum(all_counts_raw):,}")
    print(f"  Mean ZOTUs per file: {statistics.mean(all_counts_raw):.1f}")
    print(f"  Length range: {min(all_mins_raw)}-{max(all_maxs_raw)} bp")
    print(f"  Mean length across files: {statistics.mean(all_means_raw):.2f} bp")
    print(f"  Std dev of mean lengths: {statistics.stdev(all_means_raw):.2f} bp" if len(all_means_raw) > 1 else "  N/A")
    
    # Save detailed comparison
    output_file = base_dir.parent / "zotu_comparison_tactic_vs_raw.csv"
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Location', 'Mock_Name', 'Analysis_Dir', 'ZOTU_Count', 'Min_Len', 'Max_Len', 
                        'Mean_Len', 'Median_Len', 'Stdev', 'Q1', 'Q3'])
        
        for key, stats in sorted(results_by_location['tactic_results'].items()):
            mock, analysis = key.split('/')
            writer.writerow(['tactic_results', mock, analysis, stats['count'], stats['min'], stats['max'],
                           f"{stats['mean']:.2f}", stats['median'], f"{stats['stdev']:.2f}", 
                           stats['q1'], stats['q3']])
        
        for key, stats in sorted(results_by_location['raw_data'].items()):
            parts = key.split('/')
            if len(parts) == 3:
                mock, source, analysis = parts
            else:
                continue
            writer.writerow(['raw_data', f"{mock}/{source}", analysis, stats['count'], stats['min'], stats['max'],
                           f"{stats['mean']:.2f}", stats['median'], f"{stats['stdev']:.2f}", 
                           stats['q1'], stats['q3']])
    
    print(f"\n✓ Comparison saved to: {output_file}")
    
    # By-mock comparison
    output_file = base_dir.parent / "zotu_comparison_by_mock.csv"
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Mock_Name', 'Location', 'File_Count', 'Total_ZOTUs', 'Mean_ZOTUs', 
                        'Min_ZOTU_Len', 'Max_ZOTU_Len', 'Mean_ZOTU_Len', 'Stdev_ZOTU_Len'])
        
        # Group by mock
        mocks_tactic = defaultdict(list)
        for key, stats in results_by_location['tactic_results'].items():
            mock = key.split('/')[0]
            mocks_tactic[mock].append(stats)
        
        mocks_raw = defaultdict(list)
        for key, stats in results_by_location['raw_data'].items():
            mock = key.split('/')[0]
            mocks_raw[mock].append(stats)
        
        all_mocks = set(mocks_tactic.keys()) | set(mocks_raw.keys())
        
        for mock in sorted(all_mocks):
            if mock in mocks_tactic:
                stats_list = mocks_tactic[mock]
                total_zotus = sum(s['count'] for s in stats_list)
                writer.writerow([
                    mock, 'tactic_results', len(stats_list), total_zotus,
                    f"{statistics.mean([s['count'] for s in stats_list]):.1f}",
                    min([s['min'] for s in stats_list]),
                    max([s['max'] for s in stats_list]),
                    f"{statistics.mean([s['mean'] for s in stats_list]):.2f}",
                    f"{statistics.stdev([s['mean'] for s in stats_list]):.2f}" if len(stats_list) > 1 else "0"
                ])
            
            if mock in mocks_raw:
                stats_list = mocks_raw[mock]
                total_zotus = sum(s['count'] for s in stats_list)
                writer.writerow([
                    mock, 'raw_data', len(stats_list), total_zotus,
                    f"{statistics.mean([s['count'] for s in stats_list]):.1f}",
                    min([s['min'] for s in stats_list]),
                    max([s['max'] for s in stats_list]),
                    f"{statistics.mean([s['mean'] for s in stats_list]):.2f}",
                    f"{statistics.stdev([s['mean'] for s in stats_list]):.2f}" if len(stats_list) > 1 else "0"
                ])
    
    print(f"✓ Mock-wise comparison saved to: {output_file}")
    
    print("\n" + "="*100)
    print("ANALYSIS COMPLETE")
    print("="*100)

if __name__ == '__main__':
    main()
