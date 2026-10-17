#!/usr/bin/env python3
"""
Generate sequence length statistics for ZOTUs-Seqs.fasta files from TACTIC de-novo analyses.
"""

import os
import sys
from pathlib import Path
from collections import defaultdict
import statistics
import csv
from datetime import datetime
import subprocess

def parse_fasta(fasta_file):
    """Parse FASTA file and return list of sequence lengths."""
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
        print(f"Error parsing {fasta_file}: {e}", file=sys.stderr)
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

def extract_info_from_path(path):
    """Extract mock name and timestamp from analysis directory path."""
    path_parts = path.split(os.sep)
    
    # Find the Analysis_TIMESTAMP directory
    analysis_dir = None
    for part in path_parts:
        if part.startswith('Analysis_') and 'am-de-novo' in part:
            analysis_dir = part
            break
    
    if not analysis_dir:
        return None, None, None
    
    # Extract timestamp and parameters
    parts = analysis_dir.split('__')
    timestamp = parts[0].replace('Analysis_', '')  # e.g., "20260103_160944"
    
    # Find mock name in path
    mock_name = None
    for i, part in enumerate(path_parts):
        if part in ['All_Analyses', 'Isabel'] and i + 1 < len(path_parts):
            potential_mock = path_parts[i + 1]
            # Check if it's a mock community name
            if potential_mock not in ['All_Analyses', 'fastqs', 'tactic_results', 'trimmed_fastq', 'untrimmed_fastq']:
                mock_name = potential_mock
                break
    
    return timestamp, mock_name, analysis_dir

def main():
    base_dir = Path('/media/adminuser/MSD-Storage/Mohsen/TACTIC_Benchmark/Mocks')
    
    # Find all de-novo ZOTU files using subprocess for better handling of paths with spaces
    result = subprocess.run(
        ['find', str(base_dir), '-type', 'f', '-path', '*am-de-novo*/ZOTUs-Seqs.fasta'],
        capture_output=True, text=True
    )
    zotu_files = sorted([Path(p) for p in result.stdout.strip().split('\n') if p])
    
    print(f"Found {len(zotu_files)} de-novo ZOTU files\n")
    
    # Organize results by mock
    results_by_mock = defaultdict(list)
    all_results = []
    
    for i, zotu_file in enumerate(zotu_files, 1):
        file_path = str(zotu_file)
        lengths = parse_fasta(file_path)
        
        if lengths is None:
            print(f"[{i}/{len(zotu_files)}] FAILED: {file_path}")
            continue
        
        stats = calculate_stats(lengths)
        if not stats:
            print(f"[{i}/{len(zotu_files)}] EMPTY: {file_path}")
            continue
        
        timestamp, mock_name, analysis_dir = extract_info_from_path(file_path)
        
        result = {
            'timestamp': timestamp,
            'mock_name': mock_name or 'Unknown',
            'analysis_dir': analysis_dir,
            'file_path': file_path,
            'stats': stats
        }
        
        all_results.append(result)
        if mock_name:
            results_by_mock[mock_name].append(result)
        
        print(f"[{i}/{len(zotu_files)}] ✓ {mock_name or 'Unknown':<20} | "
              f"Seqs: {stats['count']:<5} | "
              f"Len: {stats['min']}-{stats['max']} | "
              f"Mean: {stats['mean']:.1f}")
    
    # Generate summary statistics CSV
    print("\n" + "="*100)
    print("Generating summary reports...\n")
    
    # Report 1: Overall statistics by mock
    output_file = base_dir / 'zotu_length_stats_summary.csv'
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Mock_Name', 'Analysis_Count', 'Total_ZOTUs', 'Min_Len', 'Max_Len', 
                        'Mean_Len', 'Median_Len', 'Stdev', 'Q1', 'Q3'])
        
        for mock_name in sorted(results_by_mock.keys()):
            analyses = results_by_mock[mock_name]
            
            # Aggregate stats
            all_lengths = []
            all_counts = []
            for result in analyses:
                all_counts.append(result['stats']['count'])
            
            for result in analyses:
                # Get individual sequence lengths if available - for now just use the stats
                all_lengths.append(result['stats']['mean'])  # Approximation
            
            total_zotus = sum(all_counts)
            mean_count = statistics.mean(all_counts)
            
            writer.writerow([
                mock_name,
                len(analyses),
                total_zotus,
                min(r['stats']['min'] for r in analyses),
                max(r['stats']['max'] for r in analyses),
                sum(r['stats']['count'] * r['stats']['mean'] for r in analyses) / total_zotus,
                statistics.median([r['stats']['median'] for r in analyses]),
                sum((r['stats']['count'] - 1) * r['stats']['stdev']**2 for r in analyses) / (total_zotus - len(analyses)) ** 0.5 if total_zotus > len(analyses) else 0,
                min(r['stats']['q1'] for r in analyses),
                max(r['stats']['q3'] for r in analyses),
            ])
    
    print(f"✓ Saved: {output_file}")
    
    # Report 2: Detailed results
    output_file = base_dir / 'zotu_length_stats_detailed.csv'
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Mock_Name', 'Timestamp', 'Analysis_Dir', 'Sequence_Count', 'Min_Length', 
                        'Max_Length', 'Mean_Length', 'Median_Length', 'Stdev', 'Q1', 'Q3'])
        
        for result in all_results:
            s = result['stats']
            writer.writerow([
                result['mock_name'],
                result['timestamp'],
                result['analysis_dir'],
                s['count'],
                s['min'],
                s['max'],
                f"{s['mean']:.2f}",
                s['median'],
                f"{s['stdev']:.2f}",
                s['q1'],
                s['q3'],
            ])
    
    print(f"✓ Saved: {output_file}")
    
    # Report 3: Text summary by mock
    output_file = base_dir / 'zotu_length_stats_by_mock.txt'
    with open(output_file, 'w') as f:
        f.write(f"ZOTU Sequence Length Statistics Report\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Total Mock Communities: {len(results_by_mock)}\n")
        f.write(f"Total Analyses: {len(all_results)}\n\n")
        
        for mock_name in sorted(results_by_mock.keys()):
            analyses = results_by_mock[mock_name]
            f.write(f"\n{'='*80}\n")
            f.write(f"Mock: {mock_name} ({len(analyses)} analyses)\n")
            f.write(f"{'='*80}\n\n")
            
            for idx, result in enumerate(analyses, 1):
                s = result['stats']
                f.write(f"  Analysis {idx}: {result['timestamp']}\n")
                f.write(f"    Directory: {result['analysis_dir']}\n")
                f.write(f"    Sequence Count: {s['count']}\n")
                f.write(f"    Length Range: {s['min']} - {s['max']} bp\n")
                f.write(f"    Mean Length: {s['mean']:.2f} bp\n")
                f.write(f"    Median Length: {s['median']} bp\n")
                f.write(f"    Std Dev: {s['stdev']:.2f}\n")
                f.write(f"    Quartiles (Q1/Q3): {s['q1']} / {s['q3']}\n\n")
    
    print(f"✓ Saved: {output_file}")
    
    # Report 4: Statistics overview
    if all_results:
        output_file = base_dir / 'zotu_length_stats_overview.txt'
        with open(output_file, 'w') as f:
            f.write(f"ZOTU Length Statistics Overview\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            # Aggregate across all analyses
            all_lengths = []
            for result in all_results:
                s = result['stats']
                # Estimate all lengths for overall statistics
                all_lengths.append((s['count'], s['min'], s['max'], s['mean'], s['median']))
            
            total_analyses = len(all_results)
            total_zotus = sum(l[0] for l in all_lengths)
            overall_min = min(l[1] for l in all_lengths)
            overall_max = max(l[2] for l in all_lengths)
            overall_mean = sum(l[0] * l[3] for l in all_lengths) / total_zotus
            
            f.write(f"Total Analyses: {total_analyses}\n")
            f.write(f"Total ZOTUs: {total_zotus}\n")
            f.write(f"Unique Mock Communities: {len(results_by_mock)}\n\n")
            
            f.write(f"Sequence Length Statistics (Across All Analyses):\n")
            f.write(f"  Minimum: {overall_min} bp\n")
            f.write(f"  Maximum: {overall_max} bp\n")
            f.write(f"  Mean: {overall_mean:.2f} bp\n")
            f.write(f"  Average ZOTUs per Analysis: {total_zotus / total_analyses:.1f}\n\n")
            
            f.write(f"Mock Community Breakdown:\n")
            for mock_name in sorted(results_by_mock.keys()):
                analyses = results_by_mock[mock_name]
                mock_zotus = sum(r['stats']['count'] for r in analyses)
                f.write(f"  {mock_name:<20}: {len(analyses):>3} analyses, {mock_zotus:>6} total ZOTUs\n")
        
        print(f"✓ Saved: {output_file}")
    
    print("\n" + "="*100)
    print(f"Summary: Processed {len(all_results)} analyses successfully")
    print(f"Output files saved to: {base_dir}")

if __name__ == '__main__':
    main()
