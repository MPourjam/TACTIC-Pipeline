# Changelog

All notable changes to the TACTIC pipeline will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-10-23

### Taxonomy & SILVA 144+ Support (Dual-Level Resolution)
- **SILVA 144+ Bacterial Kingdom Classification**: Full native compatibility with SILVA release 144+ databases by detecting and handling the formal bacterial Kingdom rank (ending in `-ati` / `-ati--other`, such as *Bacillati*, *Pseudomonadati*) between Domain and Phylum, while cleanly ignoring 8th-level source tags (such as `metagenome`).
- **Dual Taxonomy Outputs**: Simultaneously generates classic 7-level files alongside enriched 8-level files (`*-8Levels.*` across TIC, TAC, ZOTU, OTU, and `-iz`):
  - *Standard Outputs (7 levels)*: Classic `Domain;Phylum;Class;Order;Family;Genus;Species` for downstream tools strictly requiring 7 levels (e.g. PICRUSt2, legacy QIIME / R scripts).
  - *Enriched Outputs (8 levels)*: Full `Domain;Kingdom;Phylum;Class;Order;Family;Genus;Species` for high-resolution profiling, Rhea, Namco, and 8-ring Krona charts.
- **Semicolon Separator Format Standard**: Standardized taxonomy level formatting across all abundance tables and FASTA headers using `;` strictly as a level separator:
  - 7-level outputs contain exactly 6 semicolons (`len(tax.split(";")) == 7`).
  - 8-level outputs contain exactly 7 semicolons (`len(tax.split(";")) == 8`).
  - Correctly preserves multi-word taxon names with spaces (e.g. *Incertae Sedis*, *Candidatus Cryptobacteroides*, *NK4A214 group*).
- **Taxonomy Levels Guide (`README_Taxonomy_Levels.txt`)**: Automatically written into each `Analysis_<timestamp>/` results directory to guide users on when to use standard vs. 8-level files, accompanied by an explanatory summary banner in `Analysis_log.txt`.

### Pre-formatted PICRUSt2 Inputs (`PICRUSt2-Inputs/`)
- **Ready-to-Use Functional Predictions**: Added an automated `PICRUSt2-Inputs/` folder containing pre-formatted abundance tables (taxonomy column removed) and matching single-line FASTA files (pure sequence IDs, all header annotations stripped).
- **Direct CLI Compatibility**: Fully compatible out-of-the-box with `picrust2_pipeline.py -s <Seqs.fasta> -i <Table.tab> -o picrust2_out` without manual reformatting or BIOM conversions. Includes dedicated guidance in `PICRUSt2-Inputs/Output files explanation.txt`.

### Plain-Language Output Documentation (`Output files explanation.txt`)
- **Folder-Specific User Guides**: Automatically writes a clear, plain-language `Output files explanation.txt` guide inside the root `Analysis_<timestamp>/` results folder and within each active pipeline directory (`OTU_Pipeline`, `ZOTU_Pipeline`, `TAC_Pipeline`, `TIC_Pipeline`, `PICRUSt2-Inputs`, and `Shared`).
- **Accessible Guidance**: Written in simple language explaining what each output file contains, how to open it (e.g. interactive Krona HTML charts in a web browser, abundance tables in Excel or R), and its recommended scientific purpose.
- **TIC Full-Taxonomy Clarification**: Clearly explains the purpose of `ZOTU-Table-OnlyBacteria-FullTaxonomy.tab` (default for bacterial diversity and functional predictions) versus `ZOTU-Table-All-FullTaxonomy.tab` (broad-domain studies and contamination auditing).

### Clustering & Analysis Refinements
- **Multi-Analysis Mode Selector**: Supports running multiple clustering pipelines in a single command (e.g., `--analysis-mode OTU TAC TIC`), eliminating redundant preprocessing and sharing denoised sequences.
- **Clean Full-Taxonomy ZOTU Tables**: When no non-bacterial sequences are detected in TIC mode, the redundant duplicate `ZOTU-Table-All-FullTaxonomy.tab` is automatically omitted, leaving a clean single `ZOTU-Table-OnlyBacteria-FullTaxonomy.tab`. When non-bacterial sequences are present, both tables are retained.
- **Guaranteed Non-Bacterial Sequence Logging**: `Invalid-Tax-Seqs-TIC.fasta` (and its `-8Levels` counterpart) is consistently retained across all runs for auditing and downstream test compatibility.

### FASTQ Detection & Pairing
- **Robust Read Direction & Marker Matching**: Read direction is determined by the rightmost marker in the filename, supporting standard `R1`/`R2` (and `r1`/`r2`), `read1`/`read2`, `forward`/`reverse`, SRA/ENA `_1`/`_2` directly before the extension, and `F`/`R` tags.
- **Support for Hyphenated & Numbered Sample IDs**: Sample names containing digits and hyphens (such as `FU-1_..._R2_001`) are accurately classified without false forward pairings.
- **Clean Sample IDs in Mapping Files**: Supports clean sample IDs (e.g. `SampleA`) in mapping files, correctly resolving sequencer filenames like `SampleA_S388_L001_R1_001.fastq.gz` without prefix collisions (e.g., `SampleA` vs `SampleA0`).
- **Targeted Sample Lookups**: When a mapping file is provided, only requested samples are gathered, dramatically speeding up startup on large datasets.
- **Consistent Reverse Read Validation**: Re-running inside an existing project directory automatically skips prior `Analysis_*` and `*__processed` folders. Corrupted reverse mates safely fall back to single-end processing with a clear warning.
- **Single-End & Paired-End Support**: Single-end reads are automatically detected and processed alongside paired-end libraries. Unsupported `.zst` files are ignored.

### Database Reliability & Offline Support
- **Offline Fallback for Local SILVA Database**: Automatically falls back to local SILVA database files when remote download endpoints or Cloudflare R2 worker connections are inaccessible or in air-gapped environments.
- **SINA Binary Shared Library Link Restoration**: Fixed path resolution for the SINA binary and `libsina.so.0` in containerized Docker environments.
- **Automatic Permission Normalization**: Downloaded reference databases automatically set appropriate file permissions immediately after download for shared and multi-user environments.
- **Cloudflare SILVA Downloads Worker**: Upgraded `sharp` dependency to `>=0.35.5` (and updated `undici`) in `cloudflare/silva-downloads-worker/package-lock.json` to resolve security advisories.

### Settings, Logging & Reproducibility
- **Updated Pipeline Defaults**: Defaulted `abund_limit` to `0.0` and `sample_wise_correction` to `false` for TAC and TIC clustering modes in `TACTICPipeline_args.yml`, with safe fallbacks when omitted from custom configs.
- **Versioned Log Archiving**: `Pipeline_log.txt` is stamped with the pipeline version and automatically archived inside each run's `Analysis_<timestamp>/` results directory.
- **Deterministic Runs**: Sequence pairs and samples are processed in a consistent, reproducible order across different platforms and filesystems.
- **Paper Benchmarks & DOI References**: Updated step-wise and length-based 16S rRNA diversity survey benchmark datasets and documentation accompanying the TACTIC paper, including Zenodo DOI references.
