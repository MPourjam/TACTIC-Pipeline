# Read and zOTU Length Survey

This companion to the TACTICBench reproduction bundle preserves the scripts,
recoverable inputs, and reports from `TACTIC_Benchmark/length_survey`. The survey
measured read lengths before and after processing and sequence lengths in earlier
TACTIC zOTU FASTAs. It is **not** the current 16-community benchmark dataset:
different historical processing stages, extra communities, and repeated analyses
occur in these inputs.

## Contents

| Directory | Contents |
|---|---|
| `code/length_survey/` | All eight original Python scripts, unchanged, including their historical absolute paths |
| `scripts/run_length_survey.py` | Portable runner using explicit recorded input selections and the same per-file length/statistic definitions |
| `scripts/package_inputs.py` | Local-workspace recovery and ZIP builder; not needed for normal reproduction |
| `scripts/verify_input_archive.py` | CRC/SHA-256 verification of every archive member |
| `environment/` | Python-only Conda environment and dependency explanation |
| `data/small_inputs/` | Content-addressed zOTU FASTAs, small enough to accompany the code |
| `data/length_survey_inputs.zip` | Large FASTQ payloads and embedded input manifests; intended for Zenodo |
| `data/input_manifest/` | Original-to-recovered path mappings, input/observation selections, checksums, distributions and historical validation |
| `results/archive_run/` | Original reports, including results originally written outside `length_survey` |
| `results/validation_checks/` | Recorded portable-runner checks on real included data |
| `results/reproduced_run/` | Fresh portable-runner results, separate from the historical reports |

No reference database FASTA is read by any of these eight scripts. The FASTA
inputs are **inferred zOTU sequences**, not SILVA/Greengenes databases. Consequently,
this package does not add unrelated reference databases or upstream denoising
tools. It analyzes the supplied FASTQs and precomputed zOTUs; it does not regenerate
zOTUs from raw reads or rerun the upstream TACTIC pipeline.

## Reproduce

Run these commands from this bundle's root:

```bash
conda env create -f environment/length-survey-reproduction.yaml
conda activate length-survey-reproduction
python scripts/verify_input_archive.py data/length_survey_inputs.zip
unzip -n data/length_survey_inputs.zip -d data
python scripts/run_length_survey.py
```

Alternatively, avoid expanding the large FASTQs:

```bash
python scripts/run_length_survey.py --archive data/length_survey_inputs.zip
```

The sibling TACTICBench environment can run these same commands without adding
packages. No pip installation is required. To repeat one original script's
measurements, use its stem, for example:

```bash
python scripts/run_length_survey.py --script compare_zotu_locations
python scripts/run_length_survey.py --script analyze_untrimmed_fastq_lengths --archive data/length_survey_inputs.zip
```

The runner produces `input_statistics.tsv`, `survey_observations.tsv`, one TSV
per script, exact full length histograms, and `validation.json`. Reads are parsed
as four-line FASTQ records with sequence/quality validation. Each input's full
length distribution is measured once; the sampled analyses use the first 1,000
or 10,000 reads, matching the original scripts. Mean, median, sample standard
deviation, and the original order-statistic quartile definitions are retained.
Identical payloads can serve multiple recorded observations without being stored
twice. Original directory scanning, date stamps, output paths and narrative text
are deliberately not reproduced by the portable runner.

## Historical Coverage and Limits

**The old survey cannot currently be reproduced completely.**
The recovered package supports **541 available file/script observations**.
Of these, 537 have recorded numerical baselines; the remaining four are the
full-read ZIEL-I survey, for which inputs but no saved numerical baseline were recovered.
There are **94 unavailable historical measurements**, corresponding to 49 missing
zOTU FASTAs (some were measured by two scripts) and 29 trimmed-read entries.
Per-script coverage is recorded in `metadata/survey_coverage.tsv`; distinct missing
inputs are summarized in `metadata/missing_input_summary.tsv`.

`data/input_manifest/findings.tsv` explicitly lists the historical measurements
whose source files were not recovered. These include May 6, 2026 zOTU analyses
absent from the current workspace and earlier trimmed-read processing stages
whose surviving namesakes do not match the saved report statistics.
They are not replaced with newer FASTAs,
fabricated from length distributions, or silently omitted from the audit.
Saved reports for those measurements remain under `results/archive_run/`.

Rejected same-name FASTQ candidates are listed separately in
`data/input_manifest/rejected_recovery_candidates.tsv`; they are not silently
treated as the historical inputs. Each original report is used independently:
the sampled trimmed-read report and full-file trimmed-read report were produced
at different processing stages and do not necessarily select the same filenames.

`historical_validation.tsv` compares surviving inputs against the saved per-file
counts, length ranges, rounded mean/median/SD and quartiles where recorded. Exact
historical byte identity cannot be proven without contemporaneous input hashes;
the included SHA-256 values freeze the recovered inputs at packaging time.

`generate_zotu_stats.py` did not save complete input paths in its detailed CSV.
For its recovered rows, the analysis-directory name and all saved length
statistics were matched to a surviving file. This establishes matching recorded
measurements, not a unique historical tree location. The saved groups `mocks`
and `Unknown` are legacy path-parser outputs, not biological community names.

The original broad scans included repeated analysis outputs and intermediate
files. Their aggregate counts must not be interpreted as independent samples or
unique biological features. Original summary prose contains exploratory
interpretations that are preserved as provenance, **not independently validated
scientific conclusions**. The old zOTU pooled summary also uses approximate
median/SD aggregation; the portable runner reproduces and validates per-file
statistics instead of endorsing those pooled formulas.

## Zenodo Deposit

See [data/README.md](data/README.md) for archive details, checksum verification and
editable **Zenodo DOI/download URL placeholders**. Publish the small code/config/
FASTA bundle together with the large ZIP, or cross-link their separate deposits.
Read-source links are in `metadata/input_sources.tsv`; source-specific citation
and redistribution conditions still apply. Publication PDFs are not included.

See `metadata/verification_notes.md` for the completed checks and
`metadata/small_files_manifest.tsv` for checksums of the small reproduction files.
