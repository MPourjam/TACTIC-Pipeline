# Changelog

All notable changes to the TACTIC pipeline will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-10-16

### FASTQ Detection & Pairing
- **Robust Read Direction & Marker Matching**: Read direction is now determined by the rightmost marker in the filename, supporting standard R1/R2 (and r1/r2), read1/read2, forward/reverse, SRA/ENA `_1`/`_2` directly before file extension, and F/R tags.
- **Fix for Hyphenated & Numbered Sample IDs**: Sample names containing digits and hyphens (such as `FU-1_..._R2_001`) are no longer misclassified as forward reads, preventing ghost single-end samples or pairing errors regardless of filesystem file order.
- **Clean IDs in Mapping Files**: You can use clean sample IDs (e.g. `SampleA`) in mapping files; the pipeline accurately matches sequencer filenames like `SampleA_S388_L001_R1_001.fastq.gz` without prefix mix-ups (e.g., `SampleA` vs `SampleA0`).
- **Consistent Reverse Read Validation**: Both directory autodiscovery and mapping-file lookup now validate reverse mates equally. If a reverse file is invalid or corrupted, the forward read runs as single-end with a clear warning.
- **Skips Previous Results**: Re-running inside an existing project folder automatically skips previous `Analysis_*` output folders and `__processed` directories.
- **Single-End & Paired-End Support**: Single-end reads are automatically detected and processed alongside paired-end libraries. Unsupported `.zst` files are ignored.

### Taxonomy & SILVA 144+ Support
- **SILVA 144+ Kingdom Classification**: Fully compatible with SILVA release 144+ databases by detecting and handling the formal bacterial Kingdom rank (ending in `-ati` / `-ati--other`) between Domain and Phylum, while ignoring 8th-level source tags (such as `metagenome`).
- **Dual Taxonomy Outputs**: Simultaneously produces classic 7-level files (`Domain..Genus;` or `Domain..Species;`) for backward compatibility with legacy tools (PICRUSt2, classic R/QIIME scripts) alongside enriched files (`*-8Levels.*` for TAC/TIC, `*-7Levels.*` for ZOTU/OTU) containing the full SILVA 144+ bacterial Kingdom resolution.
- **Strict Trailing Semicolon Normalization**: Every taxonomic rank cleanly terminates with a trailing semicolon across all abundance tables and FASTA headers (7 semicolons for 7 levels, 8 semicolons for 8 levels), preventing rank-shifting issues in downstream parsers.
- **Automatic Output Guidance**: Automatically writes a clear `README_Taxonomy_Levels.txt` guide into each `Analysis_<timestamp>/` results directory and logs an explanatory summary banner in `Analysis_log.txt` detailing output files and tool compatibility.

### Clustering & Output Refinements
- **Clean Full-Taxonomy ZOTU Tables**: When no non-bacterial sequences are detected in TIC mode, the redundant duplicate `ZOTU-Table-All-FullTaxonomy.tab` is automatically dropped, leaving a clean single `ZOTU-Table-OnlyBacteria-FullTaxonomy.tab` (and its `-8Levels` and PICRUSt2 equivalents). When non-bacterial sequences are present, both `OnlyBacteria` and `All` tables are retained and enriched.
- **Guaranteed Non-Bacterial Sequence Logging**: `Invalid-Tax-Seqs-TIC.fasta` (and `Invalid-Tax-Seqs-TIC-8Levels.fasta`) is consistently retained across all runs for auditing and downstream test compatibility.

### Usability & Logging
- **Targeted Lookups**: When a mapping file is supplied, the pipeline only searches for the requested samples, significantly speeding up startup on large datasets.
- **Clearer Warnings**: Provides descriptive warnings in the log for missing reverse pairs, orphan reverse reads, or missing directory paths.
- **Automatic Log Archiving**: `Pipeline_log.txt` is now automatically archived inside each run's `Analysis_<timestamp>/` directory and stamped with the pipeline version.

### Settings & Stability
- **Updated Defaults**: Defaulted `abund_limit` to `0.0` and `sample_wise_correction` to `false` for TAC and TIC clustering modes in `TACTICPipeline_args.yml`, with safe fallbacks when omitted from custom configs.
- **Deterministic Runs**: Sequence pairs and samples are processed in a consistent, reproducible order across different platforms and filesystems.
- **Automatic Permission Fixes**: Reference databases (e.g., SILVA) automatically set appropriate file permissions immediately after download for shared and multi-user environments.

