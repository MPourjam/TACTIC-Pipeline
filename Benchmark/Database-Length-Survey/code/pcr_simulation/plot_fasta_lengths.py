import argparse
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from scipy.ndimage import gaussian_filter1d


def read_fasta_lengths(fasta_path):
    """
    Read a FASTA file and return sequence lengths.
    """
    lengths = []
    print(f"Reading FASTA file: {fasta_path}")
    with open(fasta_path) as fh:
        seq_parts = []
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_parts:
                    seq = "".join(seq_parts)
                    lengths.append(len(seq))
                    seq_parts = []
            else:
                seq_parts.append(line)
        if seq_parts:
            seq = "".join(seq_parts)
            lengths.append(len(seq))
    return lengths


def plot_fasta_length_histogram(fasta_path, ax, bins=50, title=None, smooth_sigma=2):
    """
    Plot histogram of sequence lengths on given axis with smooth line overlay.
    """
    lengths = read_fasta_lengths(fasta_path)
    
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
    
    print(f"{fasta_path}: Mean length: {sum(lengths)/len(lengths):.2f} bp")
    return lengths


def find_and_plot(
    dir_path: str,
    bins: int = 50,
    save_path: str = None,
    glob_pattern: str = "*.fasta",
    nrows: int = None,
    ncols: int = None,
    figsize: tuple = None,
    smooth_sigma: float = 2):
    """
    Find FASTA files in a directory and plot their length histograms in a grid.
    
    Parameters
    ----------
    dir_path : str or Path
        Directory to search for FASTA files.
    bins : int, optional
        Number of histogram bins (default: 50).
    save_path : str or Path, optional
        If given, save the figure to this path.
    glob_pattern : str, optional
        Glob pattern to find FASTA files (default: '*.fasta').
    nrows : int, optional
        Number of rows in subplot grid. Auto-calculated if None.
    ncols : int, optional
        Number of columns in subplot grid. Auto-calculated if None.
    figsize : tuple, optional
        Figure size (width, height) in inches.
    smooth_sigma : float, optional
        Sigma for Gaussian smoothing of the histogram line (default: 2).
    """
    dir_path = Path(dir_path)
    fasta_files = sorted(list(dir_path.glob(glob_pattern)))
    
    if not fasta_files:
        print(f"No FASTA files found in {dir_path} with pattern {glob_pattern}")
        return

    n_files = len(fasta_files)
    
    # Auto-calculate grid dimensions if not provided
    if nrows is None and ncols is None:
        ncols = int(np.ceil(np.sqrt(n_files)))
        nrows = int(np.ceil(n_files / ncols))
    elif nrows is None:
        nrows = int(np.ceil(n_files / ncols))
    elif ncols is None:
        ncols = int(np.ceil(n_files / nrows))
    
    if figsize is None:
        figsize = (4 * ncols, 3 * nrows)
    
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    axes = np.atleast_2d(axes).flatten()
    
    for idx, fasta_file in enumerate(fasta_files):
        plot_fasta_length_histogram(
            fasta_path=fasta_file,
            ax=axes[idx],
            bins=bins,
            smooth_sigma=smooth_sigma
        )
    
    # Hide unused subplots
    for idx in range(n_files, len(axes)):
        axes[idx].set_visible(False)
    
    fig.tight_layout()
    
    if save_path is not None:
        print(f"Saving figure to {save_path}")
        fig.savefig(save_path, dpi=300)
    else:
        default_save_path = dir_path / "combined_length_hist.png"
        print(f"Saving figure to {default_save_path}")
        fig.savefig(default_save_path, dpi=300)
    
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot histogram of FASTA sequence lengths."
    )
    parser.add_argument(
        "dir_path",
        type=str,
        help="Path to a directory containing FASTA files.",
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=50,
        help="Number of histogram bins (default: 50).",
    )
    parser.add_argument(
        "--save-path",
        type=str,
        default=None,
        help="Path to save the plot image (default: None).",
    )
    parser.add_argument(
        "--glob-pattern",
        type=str,
        default="*.fasta",
        help="Glob pattern to find FASTA files in a directory (default: '*.fasta').",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Search for FASTA files recursively in subdirectories.",
    )
    parser.add_argument(
        "--nrows",
        type=int,
        default=None,
        help="Number of rows in subplot grid (auto-calculated if not set).",
    )
    parser.add_argument(
        "--ncols",
        type=int,
        default=None,
        help="Number of columns in subplot grid (auto-calculated if not set).",
    )
    parser.add_argument(
        "--figwidth",
        type=float,
        default=None,
        help="Figure width in inches.",
    )
    parser.add_argument(
        "--figheight",
        type=float,
        default=None,
        help="Figure height in inches.",
    )
    parser.add_argument(
        "--smooth-sigma",
        type=float,
        default=2,
        help="Sigma for Gaussian smoothing of histogram line (default: 2).",
    )
    args = parser.parse_args()
    input_path = Path(args.dir_path)
    
    if not input_path.is_dir():
        print(f"Error: {input_path} is not a directory.")
        exit(1)
    
    figsize = None
    if args.figwidth is not None and args.figheight is not None:
        figsize = (args.figwidth, args.figheight)
    
    glob_pattern = args.glob_pattern
    if args.recursive:
        glob_pattern = "**/" + args.glob_pattern
    
    find_and_plot(
        dir_path=input_path,
        bins=args.bins,
        save_path=args.save_path,
        glob_pattern=glob_pattern,
        nrows=args.nrows,
        ncols=args.ncols,
        figsize=figsize,
        smooth_sigma=args.smooth_sigma,
    )
