# Benchmark inputs and FASTQ provenance

This archive extends the TACTICBench reproduction bundle for
`TACTICBENCH_RESULTS-TIC-Archaea-mock4-mock-5_SRP`.
It freezes 16 communities, 27 samples, 192 direct benchmark input
files, and 119 distinct FASTQ payloads recovered from 249 original file locations.
No TACTIC, DADA2, denoising, clustering, or tree pipeline was run to assemble it.

## Extract and reproduce the benchmark

From the reproduction bundle directory:

```bash
python scripts/verify_input_archive.py data/tacticbench_inputs_with_fastq.zip
unzip -n data/tacticbench_inputs_with_fastq.zip -d data
conda activate tacticbench-reproduction
python -m pip install --no-deps -e code/tacticbench
python -m tacticbench.cli run-all --config config/run_config.yaml
```

The ZIP has `communities/` and `input_manifest/` at its root. Extracting into `data/`
restores the benchmark's expected directory layout. Existing files are retained
with `unzip -n`. Results are written to `results/reproduced_run/` by the supplied
run config. `data/communities/` already contains the smaller benchmark inputs;
extract the archive to add the reads and additional provenance.

## What the code consumes

`discover.py` selects `source/taxonomy.tsv`, `source/expected-sequences.v4.fasta`
(with the configured full-reference fallback), DADA2 `ASV_counts.tsv`/`ASVs.fasta`,
and the feature tables/representative sequences from ZOTU, OTU, TIC, and TAC.
The original run's manifest pins the selected `Analysis_*` directories; do not
replace them with a later analysis selected by modification time. Two obsolete
ZIEL-I membership-only directories have been preserved under
`input_manifest/legacy_analysis_memberships/`, leaving only the selected analysis
in its usable community data directory. TIC/TAC
membership maps are also needed by cluster diagnostics. `prepare_expected.py`
and `prepare_methods.py` select samples, reconcile their names, and normalize
the supplied tables. TACTICBench does not read FASTQs itself.

The package retains the exact direct inputs and additional full reference FASTAs,
sample metadata, raw reads, trimmed reads, selected analysis logs and YAML files,
cluster membership and taxonomy evidence, generated DADA2 R scripts and pair lists,
and locally available SILVA 138.2 source databases. Copies of the bundle's code,
configs, environment, and archived run metadata are under `input_manifest/context/`.
Original publication and read URLs remain in the metadata snapshots; publication
PDFs and upstream container binaries are not included.

## Manifest guide

| File | Meaning |
|---|---|
| `files.tsv` | Every physical payload with original source path, archive path, role, byte size, timestamp, and SHA-256. |
| `benchmark_inputs.tsv` | All 192 files explicitly selected by the archived run manifest, with identity checks against the bundle. |
| `fastq_files.tsv` | Every recovered input FASTQ location, storage/alias status, SHA-256 of stored and decompressed bytes, read counts, length ranges, and first/last read headers. |
| `fastq_pairs.tsv` | Mate counts and ordered read-ID checks in each original directory layout. |
| `pipeline_inputs.tsv` | One row per upstream sample: selected run, exact log command and line number, recovered TACTIC and DADA2 mates, and logged spike-removal counts. |
| `sample_input_links.tsv` | Links the archived expected/inferred sample names to mates for all five methods. |
| `provenance_evidence.tsv` | Exact commands and input/read-count evidence with archive paths and log line numbers. |
| `zero_length_reads.tsv` | Retained empty mates, which can explain differences between physical records and a tool's reported paired subset. |
| `limitations_and_findings.tsv` | Concrete differences, provenance qualifications, and requirements for an upstream replay. |
| `benchmark_preparation_validation.tsv` | Independently regenerates only benchmark preparation in a temporary output and compares it to archived tables. |
| `validation.json` | Input inventory, FASTQ checks, and any errors detected while collecting files. |
| `SHA256SUMS` | Hashes of every ZIP member except this checksum list itself. |

Archive paths are relative to `data/` after extraction. Original paths are evidence
from the analysis machine and need not exist on another computer. Byte-identical
FASTQ copies are stored once; alias rows resolve to that stored archive path. This
includes DADA2 work/trimmed copies and renamed raw downloads. The manifest makes
the usable archive path explicit instead of relying on obsolete source MANIFESTs.
The table `files.tsv` lists physical payloads; `fastq_files.tsv` retains all input
locations including aliases. Filtered DADA2 reads and disposable analysis
intermediates are not required benchmark inputs and are omitted.

## Special handling and reproducibility limits

For mocks 18 and 19, the selected TACTIC runs received reads containing spikes,
then removed spike pairs internally. The original pre-removal trimmed reads were
subsequently moved/renamed under `trimmed_but_spiked/`. DADA2's retained input
copies are the later `spikefree` reads. Both are archived. The per-sample log
comparison checks recovered pre-removal read counts against the TACTIC logs;
it does not assume that the currently named spikefree file was the original
TACTIC input. Consult `pipeline_inputs.tsv` before replaying either pipeline.

For SRR13005964, both retained FASTQs contain 147918 records. Its historical
Bowtie2 output reports 147918 input records, comprising 147917 paired records
and one unpaired read. R2 of `SRR13005964.10542` is empty. TACTIC's shorter
summary reports the paired subset, explaining the apparent one-pair difference.
These input records are preserved; the manifests distinguish physical records,
reported paired/unpaired counts, and the subset remaining after spike removal.
Read-ID checks normalize SRA's `.1`/`.2` mate suffix only when it follows a
separate spot number (for example `SRR2990088.1.1` and `SRR2990088.1.2`), and
also handle `/1` and `/2` suffixes. Spot numbers are retained.

The TIC input files for mocks 4 and 5 incorporate the archaeal corrections made
after their original upstream runs. The exact corrected feature tables, FASTAs,
and membership maps are retained. They permit reproducing the benchmark; a fresh
unmodified upstream TACTIC run can differ from them. Available before-correction
tables/maps are retained as evidence, but no complete automated correction recipe
or immutable upstream container image is present.

Historical logs recorded filenames and read counts but no original FASTQ hashes.
The supplied hashes freeze the recovered files as packaged now. Agreement of
names, counts, pair identities, and DADA2 copies supports their provenance, but
does not replace an unavailable historical checksum. Locally available SILVA
inputs are included; container-internal spike reference sequences, the original
Docker image digest, and original binary images are not archived. Therefore the
bundle supports benchmark reproduction from the exact inferred tables, while an
end-to-end regeneration of those tables from reads still requires the recorded
upstream software and the documented corrections. No claim of rerunning or
verifying all figures from raw reads is made.
Generated DADA2 R scripts and pair lists are retained, but their exact historical
positional argument values and immutable container digest were not recorded.
Current driver scripts are evidence and can contain early exits; they are not
presented as a verified replay of every original upstream invocation.

The build and validation scripts are in `input_manifest/scripts/` and in the
bundle's `scripts/` directory. The archive uses ZIP64 and retains gzip files
without recompression. `data/tacticbench_inputs_with_fastq.zip.sha256` checks
the outer archive; `data/input_archive_validation.json` records the full ZIP
CRC and member SHA-256 verification performed after writing it.
