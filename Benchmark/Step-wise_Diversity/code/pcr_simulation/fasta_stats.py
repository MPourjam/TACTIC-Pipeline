import sys
import statistics
from typing import List, Tuple
from pathlib import Path
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt


def get_sequence_lengths(fasta_file) -> Tuple[List[int], int]:
    """
    Parses a FASTA file and returns a list of sequence lengths.
    Handles multi-line FASTA sequences.
    """
    lengths = []
    current_length = 0
    current_sequence = []
    in_sequence = False

    try:
        with open(fasta_file, 'r') as f:
            unique_sequences = set()
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                if line.startswith(">"):
                    if in_sequence:
                        lengths.append(current_length)
                        unique_sequences.add("".join(current_sequence))
                    current_length = 0
                    current_sequence = []
                    in_sequence = True
                else:
                    current_length += len(line)
                    current_sequence.append(line)
            
            # Append the last sequence
            if in_sequence:
                lengths.append(current_length)
                unique_sequences.add("".join(current_sequence))
                
    except FileNotFoundError:
        print(f"Error: File '{fasta_file}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file: {e}", file=sys.stderr)
        sys.exit(1)

    return lengths, len(unique_sequences)


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


def fasta_stats(fasta_file) -> dict:
    """
    Docstring for fasta_stats

    :param fasta_file: Description
    :return: Description
    :rtype: dict
    """
    lengths, n_unique = get_sequence_lengths(fasta_file)
    if not lengths:
        return {}
    # Calculate statistics
    count = len(lengths)
    total_bases = sum(lengths)
    min_len = min(lengths)
    max_len = max(lengths)
    mean_len = statistics.mean(lengths)
    median_len = statistics.median(lengths)
    n50 = calculate_n50(lengths)
    # Standard deviation requires at least two data points
    stdev = statistics.stdev(lengths) if count > 1 else 0.0
    # Construct dictionary for JSON output
    stats_data = {
        "file": fasta_file,
        "total_sequences": count,
        "n_unique_sequences": n_unique,
        "total_bases": total_bases,
        "min_length": min_len,
        "max_length": max_len,
        "mean_length": round(mean_len, 2),
        "median_length": round(median_len, 2),
        "n50": n50,
        "stdev": round(stdev, 2)
    }

    return stats_data


def plot_fasta_length_histogram(
        fasta_path,
        ax=None,
        bins=50,
        title=None,
        smooth_sigma=2) -> str:
    """
    Plot histogram of sequence lengths on given axis with smooth line overlay.
    Saves the plot next to the fasta file with .png suffix.
    """
    
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure  # Get the figure from the provided axes
    
    lengths, _ = get_sequence_lengths(fasta_path)
    
    if not lengths:
        print(f"No sequences found in {fasta_path}")
        ax.set_visible(False)
        return []

    counts, bin_edges, _ = ax.hist(lengths, bins=bins, alpha=0.7, edgecolor='black')
    
    # Plot smooth line over histogram
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    smooth_counts = gaussian_filter1d(counts, sigma=smooth_sigma)
    ax.plot(bin_centers, smooth_counts, 'r-', linewidth=2, label='Smoothed')
    
    ax.set_xlabel("Sequence length (bp)")
    ax.set_ylabel("Count")
    if title is None:
        title = Path(fasta_path).name
    ax.set_title(title, fontsize=10)
    ax.legend()
    
    # Save the plot next to the fasta file
    fasta_path = Path(fasta_path)
    output_path = fasta_path.with_suffix('.length_histogram.png')
    ax.figure.savefig(output_path, dpi=150, bbox_inches='tight')
    fig.close()  # Close THIS specific figure

    return str(output_path.absolute())
