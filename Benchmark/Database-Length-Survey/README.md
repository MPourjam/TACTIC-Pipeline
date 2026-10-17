# Step-wise 16S Diversity Reproduction

This collection reproduces the database sequence-length and clustering analyses
in `Step-wise_Diversity`. It includes the original working-tree Python scripts,
the exact prepared aligned reference FASTAs used as analysis inputs, archived
numerical results and figures, portable launchers, and a minimal Conda environment.

The analyses cut increasing windows from aligned 16S reference sequences, remove
alignment gaps, count unique sequences and taxonomy labels, and cluster sequences
with USEARCH 11. The results show how sequence length and clustering identity
change the number of distinguishable features. The Greengenes2 analyses also
compare unique sequences and clusters for the most frequent species labels.

## Data and Settings

| Dataset | Prepared input | V1-V9 alignment coordinates (1-based, inclusive) | Length filter after removing gaps | Clustering identities |
|---|---|---|---|---|
| SILVA 138.2 SSURef NR99 | Full known bacterial taxonomy, aligned | 1159-43016 | 1200-1600 nt | 0.987, 0.97 |
| LTP 10_2024 | Bacterial sequences with seven taxonomy levels, aligned | 1150-43026 | 1200-1600 nt | 0.987, 0.97 |
| Greengenes2 2024.09 | Full species taxonomy, aligned to SILVA coordinates | 1150-43026 | 1200-1600 nt | 0.995, 0.987, 0.97 |

All three analyses increase window size in 2000-alignment-column steps, starting
at the first column of the filtered V1-V9 alignment. Each produces 20 window rows,
one full-reference row, and one filtered V1-V9 row. Alignment columns include
gaps; the graph's x coordinates are mean sequence lengths after gap removal.

The executable settings and input/run paths are recorded in `config/runs.json`.
The historical filenames contain `1250-1650`, but the actual length filter is
**1200-1600**, as confirmed by the source entry point and saved length statistics.
The retained headers must pass the original seven-level taxonomy check. The
historical SILVA and LTP summaries count distinct two-word species names and
preserve the original FASTA headers. GG2 uses the later strict taxonomy profile,
which normalizes headers and counts full seven-level taxonomy strings. These
dataset-specific definitions are explicit in `config/runs.json`; the portable
adapter in `scripts/historical_taxonomy.py` preserves the older definition.

**USEARCH must use 19 threads to reproduce the archived cluster counts.** The
original wrapper selected 80% of the 24 available cores. USEARCH 11's parallel
greedy clustering can produce different counts with another thread count, even
for identical sorted sequences. `--threads` controls this calculation and defaults
to 19. Adjust `--workers` to limit concurrent window jobs while keeping 19 threads
per clustering call when reproducing the historical results.

## Contents

| Location | Contents |
|---|---|
| `code/pcr_simulation/` | Snapshot of the original Python scripts and primer tables |
| `config/runs.json` | Portable settings for the three archived analyses |
| `environment/stepwise-diversity-reproduction.yaml` | Required Conda packages only |
| `tools/` | Exact Linux x86-64 USEARCH 11 executable and attribution |
| `data/stepwise_diversity_inputs.zip` | Three prepared aligned FASTAs, the GG2 header list, and manifests |
| `data/input_manifest/files.tsv` | Original and portable paths, sizes, roles, source links, and SHA-256 checksums |
| `results/archive_run/` | Archived summary tables, per-window statistics, and figures |
| `results/validation_checks/` | Recorded checks carried out while preparing this collection |
| `metadata/` | Source commit, working-tree status, file provenance, and earlier workflow files |
| `scripts/` | Portable input verification, analysis, plotting, and packaging commands |

## Reproduce

Commands below are run from this directory. Use Linux x86-64 for the supplied
USEARCH executable. The FASTA archive and small code collection should be kept
together, or their deposits should be linked explicitly. The ZIP download details
and future Zenodo DOI belong in [data/README.md](data/README.md).
After extracting `Step-wise_Diversity_code_and_results.zip`, place the downloaded
`stepwise_diversity_inputs.zip` under this directory's `data/` folder.

Create the minimal environment:

```bash
conda env create -f environment/stepwise-diversity-reproduction.yaml
conda activate stepwise-diversity-reproduction
python scripts/run_stepwise.py check
```

The environment needs Python, NumPy, Matplotlib, and SeqKit. The supplied USEARCH
11 binary provides clustering. The old read/zOTU length-survey environment needs
NumPy, Matplotlib, and SeqKit added to run this project. The TACTICBench environment
also needs SeqKit; use this project's environment for the recorded tool versions.

Regenerate graphs from the saved numerical summaries without downloading or
extracting the large reference archive:

```bash
python scripts/run_stepwise.py plots --outdir results/reproduced_plots
```

This invokes the original feature-versus-length plotting function for all three
databases and renders GG2 top-25 and vertical top-5 species figures from the saved
measurements. The archived SVGs preserve the original presentation, including
later figure variants. The summary-only species renderer uses a simple grouped
bar layout and preserves numerical values, capitalized genus names, and italics.

Verify the ZIP checksum and every FASTA payload, and extract the inputs:

```bash
python scripts/prepare_inputs.py --extract
```

Run a short end-to-end check using 64 real records from each database, all 20
window lengths, and every configured clustering identity. The smoke selection
retains the first 64 records with valid taxonomy, an eligible V1-V9 length, and
at least 32 ungapped bases in the shortest window, so tiny fixtures are nonempty:

```bash
python scripts/run_stepwise.py reproduce --smoke-records 64 --outdir results/smoke_run
```

Recompute all statistics and species summaries from the full input databases:

```bash
python scripts/run_stepwise.py reproduce --threads 19 --workers 3
```

The full run writes fresh intermediates and figures under
`results/reproduced_run/<dataset>/run_<date>_clust.../start.../`. Its
`validation.json` compares all numerical columns in all 22 rows against the
archived tables and returns a nonzero exit code if any value differs. Commas in
formatted counts are ignored; numerical differences are reported individually.
Use a fresh output directory for each rerun. This workflow generates hundreds of
GB of intermediates; allow at least 300 GB of free working space and several GB
of RAM. Completion time depends on storage, CPU, and database size.

To run just one dataset, append `--datasets ltp`, `--datasets silva`, or
`--datasets gg2`. `prepare_inputs.py` accepts the same selection. `--data-dir`
and `--outdir` allow inputs and working outputs on another disk. `--usearch`
selects a different executable explicitly and records its version and checksum.

Two focused commands validate the selected references and reproduce the GG2
species figures without recomputing the entire window-clustering sweep:

```bash
python scripts/audit_reference_populations.py --outdir results/reference_audit
python scripts/reproduce_species_summaries.py --outdir results/reproduced_species
```

The population audit checks the record and species counts before and after the
V1-V9 length filter for all three complete databases, plus the GG2 top-25 species
frequencies. The species command recomputes the individual and pooled top-25 and
top-5 unique sequences and cluster counts at the two thresholds used in those
archived species figures (0.987 and 0.97). Both commands compare against saved
measurements and return a nonzero exit code for discrepancies.

To exercise the launchers from a relocated copy of the small collection, use a
new test directory and point to the extracted reference inputs:

```bash
python scripts/test_reproduction.py --workdir /tmp/stepwise-isolated-test --data-dir data
```

This checks tools, regenerates plots, and runs all three smoke analyses with
`PYTHONPATH` cleared and only the active environment plus system utilities on
`PATH`. Add `--full-datasets ltp` to include a complete LTP numerical rerun.

## Scope and Provenance

The prepared FASTAs are the exact inputs required by the current alignment-cut
workflow. Their original upstream alignment and taxonomy preparation is preserved
through the files and source references, rather than being recreated by a separate
SINA/QIIME2 database-building pipeline. Downloads of current database releases may
differ from these exact prepared inputs; use the checksummed archive for reruns.

The original scripts are copied from the working tree, including local plotting
edits; `metadata/source_snapshot.json` records the commit and dirty status.
The portable launcher sets input/output paths, process limits, and the executable
location. It also turns external command failures into ordinary worker exceptions
and builds the shared SeqKit index before parallel cutting. This prevents workers
from reading a partially written index and producing empty window files. The
original module performs cutting, length
filtering, deduplication, and clustering. The documented historical taxonomy
adapter preserves the species-label definition and original headers used by the
older SILVA/LTP runs; the bundled source snapshot remains unmodified.

The older primer-sweep/Nextflow materials are retained as historical context under
`metadata/legacy_workflow/`. The validated commands above use the alignment-cut
workflow. Running the old helper branches can require additional tools such as
Cutadapt and SciPy; these are not dependencies of the analyses reproduced here.
Large generated FASTAs, temporary clusters, caches, and duplicate database exports
are regenerated by the workflow and are excluded from the deposit.

Database origins: [SILVA](https://www.arb-silva.de/),
[LTP](https://bioinfo.uib.es/ltp/), and
[Greengenes2 2024.09](https://ftp.microbio.me/greengenes_release/2024.09/).
The exact USEARCH binary is available in the author's
[public-domain archive](https://github.com/rcedgar/usearch_old_binaries/).
See [tools/README.md](tools/README.md) for its version and checksum.

For the completed test results and any numerical differences from the historical
outputs, see [results/validation_checks/README.md](results/validation_checks/README.md).
