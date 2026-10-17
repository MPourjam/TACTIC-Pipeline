# Minimal Reproduction Environment

`tacticbench-reproduction.yaml` targets the bundled TACTICBench module and its
active run config, not the upstream TACTIC or DADA2 pipelines.

| Dependency | Reason |
|---|---|
| Python 3.10 | Compatible interpreter for the bundled code |
| NumPy, pandas, PyYAML | Numerical calculations, input tables and YAML configuration |
| Matplotlib, Seaborn | Figures, including the enabled Figure 21 KDE/boxplot panel |
| Biopython | FASTA and Newick tree parsing |
| BLAST+ | `blastn` / `makeblastdb` reference mapping |
| MAFFT, FastTree | Expected-sequence alignments and reference trees |
| pip, setuptools, wheel | Editable installation of `code/tacticbench` |

The six scientific Python libraries match `code/tacticbench/pyproject.toml`.
Transitive dependencies (for example SciPy used by Seaborn KDE) are resolved by
Conda; they are not unrelated packages copied from a general-purpose environment.
UMAP is disabled in the supplied run config, so `umap-learn` is omitted. Optional
historical exploratory wave-plot scripts may need Floweaver/widgets, but they are
not used by the configured main figures.

All eight scripts in the sibling `length_survey_reproduction` project use only
the Python standard library. This TACTICBench environment can run that project
unchanged, but its own minimal environment requires only Python 3.10. Its portable
runner also avoids the legacy survey script's external `find` command.

This is a minimal viable environment specification, not a complete platform lock
or a promise of bit-identical renderings across package versions. Historical
inputs and outputs are checksummed separately. The available reproduction
environment passed the benchmark preparation audit; a fresh Conda solve with the
added Seaborn dependency has not been installed as part of packaging.
