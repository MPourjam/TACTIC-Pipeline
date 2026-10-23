# Changelog

All notable changes to the TACTIC pipeline will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-10-23

### Output Documentation & Usability
- **Simple-Language Pipeline Output Guides (`Output files explanation.txt`)**: Automatically writes a clear, plain-language `Output files explanation.txt` file into the root `Analysis_<timestamp>/` results directory and each active pipeline folder (`OTU_Pipeline`, `ZOTU_Pipeline`, `TAC_Pipeline`, `TIC_Pipeline`, `PICRUSt2-Inputs`, and `Shared`). Explains what each file is, how to open it (e.g. interactive Krona HTML charts in a web browser, abundance tables in Excel/R), and its intended scientific usage.
- **TIC Full-Taxonomy Guidance**: Clearly documents the role and conditional generation of `ZOTU-Table-All-FullTaxonomy.tab` (retained when non-bacterial sequences are detected for broad-domain studies and contamination auditing) versus `ZOTU-Table-OnlyBacteria-FullTaxonomy.tab` (default for standard bacterial diversity and functional predictions with PICRUSt2).
- **PICRUSt2 Input Guidance**: Added dedicated documentation in `PICRUSt2-Inputs/Output files explanation.txt` explaining pre-formatted inputs and providing exact CLI syntax for running `picrust2_pipeline.py`.

### Taxonomy Formatting Logic
- **Semicolon Separator Format**: Standardized taxonomy formatting to use `;` as a level separator. Standard 7-level outputs contain exactly 6 semicolons (`Domain;Phylum;Class;Order;Family;Genus;Species`), and SILVA 144+ 8-level outputs contain exactly 7 semicolons (`Domain;Kingdom;Phylum;Class;Order;Family;Genus;Species`).

### Security & Dependencies
- **Cloudflare SILVA Downloads Worker**: Upgraded `sharp` dependency to `>=0.35.5` (and updated `undici`) in `cloudflare/silva-downloads-worker/package-lock.json` to resolve security advisories.

### Environment & Database Reliability
- **SINA Binary Shared Library Link Restoration**: Fixed path resolution for the SINA binary and `libsina.so.0`, resolving shared library loading issues in containerized Docker environments.
- **Offline Fallback for Local SILVA Database**: Added automatic fallback to local SILVA database files when remote download endpoints or Cloudflare R2 worker connections are inaccessible or in air-gapped environments.

### Benchmark Data & Documentation
- **Paper Benchmarks & DOI References**: Updated step-wise and length-based 16S rRNA diversity survey benchmark datasets and documentation accompanying the TACTIC paper, including Zenodo DOI references.


### FASTQ Detection & Pairing
- **Robust Read Direction & Marker Matching**: Read direction is now determined by the rightmost marker in the filename, supporting standard R1/R2 (and r1/r2), read1/read2, forward/reverse, SRA/ENA `_1`/`_2` directly before file extension, and F/R tags.
- **Fix for Hyphenated & Numbered Sample IDs**: Sample names containing digits and hyphens (such as `FU-1_..._R2_001`) are no longer misclassified as forward reads, preventing ghost single-end samples or pairing errors regardless of filesystem file order.
- **Clean IDs in Mapping Files**: You can use clean sample IDs (e.g. `SampleA`) in mapping files; the pipeline accurately matches sequencer filenames like `SampleA_S388_L001_R1_001.fastq.gz` without prefix mix-ups (e.g., `SampleA` vs `SampleA0`).
- **Consistent Reverse Read Validation**: Both directory autodiscovery and mapping-file lookup now validate reverse mates equally. If a reverse file is invalid or corrupted, the forward read runs as single-end with a clear warning.
- **Skips Previous Results**: Re-running inside an existing project folder automatically skips previous `Analysis_*` output folders and `__processed` directories.
- **Single-End & Paired-End Support**: Single-end reads are automatically detected and processed alongside paired-end libraries. Unsupported `.zst` files are ignored.

### Taxonomy & SILVA 144+ Support
- **SILVA 144+ Kingdom Classification**: Fully compatible with SILVA release 144+ databases by detecting and handling the formal bacterial Kingdom rank (ending in `-ati` / `-ati--other`) between Domain and Phylum, while ignoring 8th-level source tags (such as `metagenome`).
- **Dual Taxonomy Outputs**: Simultaneously produces classic 7-level files (`Domain..Species;`) for backward compatibility with legacy tools (PICRUSt2, classic R/QIIME scripts) alongside 8-level enriched files (`*-8Levels.*` across TIC, TAC, ZOTU, OTU, and `-iz`) containing the full SILVA 144+ bacterial Kingdom resolution down to Species.
- **Consistent Semicolon Separator Normalization**: Standardized taxonomic level formatting across all abundance tables and FASTA headers using ';' as a separator (6 semicolons for 7 levels, 7 semicolons for 8 levels), preventing rank-shifting issues in downstream parsers.
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

