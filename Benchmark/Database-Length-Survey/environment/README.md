# Dependencies

The eight scripts in `code/length_survey/` import only Python standard-library
modules: `csv`, `gzip`, `os`, `sys`, `pathlib`, `collections`, `statistics`,
`datetime`, and `subprocess`. The portable packaging, verification, and rerun
utilities also use only the standard library. Python 3.10 is sufficient.

The original `generate_zotu_stats.py` shell-discovery step called Unix `find`.
The portable runner reads a frozen manifest instead and does not require `find`,
BLAST, MAFFT, FastTree, Biopython, pandas, NumPy, Matplotlib, Seaborn or Conda itself
if a suitable Python interpreter is already installed.

`length-survey-reproduction.yaml` therefore lists only Python. The sibling
TACTICBench environment is compatible and needs no additional packages, but would
install unnecessary scientific/mapping dependencies for this survey alone.

This file specifies a minimal viable environment, not a lockfile. No new Conda
environment is installed by the packaging procedure. The portable runner is
tested with an already available Python interpreter, including Python 3.10 from
the TACTICBench reproduction environment.
