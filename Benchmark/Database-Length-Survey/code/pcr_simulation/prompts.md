# Prompts
## `main.py`
A CLI and argument parse for this script:
This script takes as input SILVA or GreenGenes (GG2) fasta files of full 16S rRNA sequences and one or more primers sets (both forward and reverse) and other related parameters adjustable for run_cutadapt(). It generates combinations of parameters for run_cutadapt and stores the results in and output directory (adjustable via the CLI or main function of this script). The output directory should look like "{SILVA or GG2}/err0.2_{rest of parameters}/{output of run_cutadapt}"

When running cutadapt is done with the parameters sets (in multi-threaded mode to speed up the process). The script runs the run_stats_fasta() function on trimmed fasta file output of run_cutadapt. and also runs the find_and_plot() function on each fasta file generated for each run_cutadapt()

