#!/usr/bin/env python3
"""
Comprehensive read length survey for trimmed_fastq files from all mock communities in raw_data.
Analyzes read lengths to understand trimming consistency across mocks.
"""

import gzip
from pathlib import Path
from collections import Counter, defaultdict
import statistics
import csv

def analyze_fastq_lengths(fastq_file, max_reads=None):
    """Extract read lengths from a fastq file (full coverage or sampled)."""
    lengths = []
    read_count = 0
    
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
                    read_count += 1
                    if max_reads and read_count >= max_reads:
                        break
    except Exception as e:
        print(f"Error reading {fastq_file}: {e}")
        return None
    
    return lengths

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
    base_dir = Path("/media/adminuser/MSD-Storage/Mohsen/TACTIC_Benchmark/Mocks/Benchmark_git/All_Analyses/mocks/raw_data")
    
    # Find all mock directories with trimmed_fastq
    mock_dirs = sorted([d for d in base_dir.glob('mock-*/trimmed_fastq')])
    
    print("\n" + "="*100)
    print("COMPREHENSIVE TRIMMED FASTQ READ LENGTH SURVEY - ALL MOCK COMMUNITIES")
    print("="*100)
    print(f"\nFound {len(mock_dirs)} mock communities with trimmed_fastq data\n")
    
    all_results_by_mock = defaultdict(dict)
    
    for mock_path in mock_dirs:
        mock_name = mock_path.parent.name
        
        # Find all fastq files in this mock's trimmed_fastq directory
        fastq_files = sorted(
            [f for f in mock_path.iterdir() if f.is_file() and (f.name.endswith('.fastq') or f.name.endswith('.fastq.gz'))]
        )
        
        if not fastq_files:
            print(f"\n⚠ {mock_name}: No fastq files found")
            continue
        
        print(f"\n{'='*100}")
        print(f"Mock: {mock_name} ({len(fastq_files)} files)")
        print(f"{'='*100}")
        
        file_results = {}
        
        for i, fastq_file in enumerate(fastq_files, 1):
            # Sample first 10,000 reads for speed
            lengths = analyze_fastq_lengths(fastq_file, max_reads=10000)
            
            if not lengths:
                print(f"  [{i}/{len(fastq_files)}] ✗ FAILED: {fastq_file.name}")
                continue
            
            stats = calculate_stats(lengths)
            length_counter = Counter(lengths)
            unique_lengths = sorted(length_counter.keys())
            
            file_results[fastq_file.name] = stats
            
            if len(unique_lengths) == 1:
                print(f"  [{i}/{len(fastq_files)}] ✓ {fastq_file.name:<50} {stats['count']:>6,} reads │ {unique_lengths[0]} bp")
            else:
                print(f"  [{i}/{len(fastq_files)}] ⚠ {fastq_file.name:<50} {stats['count']:>6,} reads │ {stats['min']}-{stats['max']} bp")
        
        all_results_by_mock[mock_name] = file_results
    
    # Summary by mock
    print("\n" + "="*100)
    print("SUMMARY BY MOCK COMMUNITY")
    print("="*100 + "\n")
    
    # Generate summary statistics
    summary_data = []
    
    for mock_name in sorted(all_results_by_mock.keys()):
        file_results = all_results_by_mock[mock_name]
        
        if not file_results:
            continue
        
        # Aggregate stats
        all_counts = [s['count'] for s in file_results.values()]
        all_mins = [s['min'] for s in file_results.values()]
        all_maxs = [s['max'] for s in file_results.values()]
        all_means = [s['mean'] for s in file_results.values()]
        all_medians = [s['median'] for s in file_results.values()]
        
        total_reads = sum(all_counts)
        min_across = min(all_mins)
        max_across = max(all_maxs)
        mean_across = statistics.mean(all_means)
        
        # Check consistency
        unique_mins = set(all_mins)
        unique_maxs = set(all_maxs)
        
        if len(unique_mins) == 1 and len(unique_maxs) == 1:
            consistency = "✓ CONSISTENT"
        elif len(unique_mins) <= 2 and len(unique_maxs) <= 2:
            consistency = "⚠ MOSTLY CONSISTENT"
        else:
            consistency = "✗ VARIABLE"
        
        print(f"{mock_name:<20}")
        print(f"  Files: {len(file_results)}")
        print(f"  Total reads (sampled): {total_reads:,}")
        print(f"  Length range: {min_across}-{max_across} bp")
        print(f"  Mean length: {mean_across:.2f} bp")
        print(f"  Consistency: {consistency}")
        if len(unique_mins) > 1:
            print(f"    Min values: {sorted(unique_mins)}")
        if len(unique_maxs) > 1:
            print(f"    Max values: {sorted(unique_maxs)}")
        print()
        
        summary_data.append({
            'mock': mock_name,
            'files': len(file_results),
            'total_reads': total_reads,
            'min_length': min_across,
            'max_length': max_across,
            'mean_length': mean_across,
            'median_length': statistics.mean(all_medians),
            'consistency': consistency.split()[0]  # Just the symbol
        })
    
    # Save summary to CSV
    output_file = base_dir.parent / "trimmed_fastq_length_summary.csv"
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Mock_Community', 'File_Count', 'Total_Reads_Sampled', 'Min_Length', 
                        'Max_Length', 'Mean_Length', 'Median_Length', 'Consistency_Status'])
        
        for data in summary_data:
            writer.writerow([
                data['mock'],
                data['files'],
                data['total_reads'],
                data['min_length'],
                data['max_length'],
                f"{data['mean_length']:.2f}",
                f"{data['median_length']:.2f}",
                data['consistency']
            ])
    
    print(f"✓ Summary saved to: {output_file}\n")
    
    # Detailed file comparison
    output_file = base_dir.parent / "trimmed_fastq_length_detailed.csv"
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Mock_Community', 'Fastq_File', 'Total_Reads_Sampled', 'Min_Length', 
                        'Max_Length', 'Mean_Length', 'Median_Length', 'Stdev'])
        
        for mock_name in sorted(all_results_by_mock.keys()):
            file_results = all_results_by_mock[mock_name]
            for fname in sorted(file_results.keys()):
                stats = file_results[fname]
                writer.writerow([
                    mock_name,
                    fname,
                    stats['count'],
                    stats['min'],
                    stats['max'],
                    f"{stats['mean']:.2f}",
                    stats['median'],
                    f"{stats['stdev']:.2f}"
                ])
    
    print(f"✓ Detailed comparison saved to: {output_file}\n")
    
    # Overall statistics
    print("="*100)
    print("OVERALL STATISTICS ACROSS ALL MOCKS")
    print("="*100 + "\n")
    
    all_means = [d['mean_length'] for d in summary_data]
    all_totals = [d['total_reads'] for d in summary_data]
    all_mins_overall = [d['min_length'] for d in summary_data]
    all_maxs_overall = [d['max_length'] for d in summary_data]
    
    print(f"Total mock communities: {len(summary_data)}")
    print(f"Total reads sampled: {sum(all_totals):,}")
    print(f"Overall sequence length range: {min(all_mins_overall)}-{max(all_maxs_overall)} bp")
    print(f"Mean sequence length (across mocks): {statistics.mean(all_means):.2f} bp")
    print(f"Std dev of mean lengths: {statistics.stdev(all_means):.2f} bp" if len(all_means) > 1 else "N/A")
    
    # Count consistency
    consistent_count = sum(1 for d in summary_data if '✓' in d['consistency'])
    variable_count = sum(1 for d in summary_data if '✗' in d['consistency'])
    mostly_count = sum(1 for d in summary_data if '⚠' in d['consistency'])
    
    print(f"\nConsistency Summary:")
    print(f"  ✓ Fully consistent: {consistent_count} mocks")
    print(f"  ⚠ Mostly consistent: {mostly_count} mocks")
    print(f"  ✗ Variable: {variable_count} mocks")

if __name__ == '__main__':
    main()
