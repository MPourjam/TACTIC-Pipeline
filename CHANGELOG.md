# Changelog

All notable changes to the TACTIC pipeline will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Each version ends with **Changes that affect your results**: what can make the numbers of zOTUs, OTUs or SOTUs,
the counts, the taxonomy or the normalized tables differ from the previous version. A table lists every such change,
whom it affects (**BiotaWiz** users, or **command-line** users who start TACTIC themselves with `docker run`; both
run the same Docker image) and **whether you can reproduce the earlier result with the new version**; the details
follow below the table. Read it before comparing results of different versions. The TACTIC version used is the first
line of `Analysis_log.txt`.

Parameters are named as on the **Parameters** step of BiotaWiz, followed by their name in the arguments file, e.g.
**Minimum relative abundance** (`ClusterZOTUsArgs` › `abund_limit`).

**Effect** grades how much a change can alter the **biological interpretation** of the runs it applies to: richness
and alpha diversity, beta diversity, taxonomic composition, absolute abundance and comparisons between groups. Five
levels:

| Effect | Biological interpretation |
|---|---|
| **1 Very low** | No biological result changes; only technical details of the files differ. |
| **2 Low** | Small numeric differences; diversity, composition and group differences stay the same. |
| **3 Medium** | One kind of result can change for some taxa or samples (e.g. absolute abundances, some taxon names, zOTU richness); the main conclusions usually hold. Recheck that analysis. |
| **4 High** | The conclusions of one kind of analysis can change (e.g. absolute abundances between groups, group statistics, single samples). Redo that analysis with the new version before comparing. |
| **5 Very high** | The biological picture changes: richness, alpha or beta diversity or the taxonomic composition differ substantially. Conclusions drawn from earlier results may not hold. |

The grades are estimates from what each change does. Only the one for **Minimum relative abundance** (1.0.0) is
measured.

## How to reproduce results of an earlier version

1. **With a parameter**, where a change below says so: set the value given there on the **Parameters** step of
   BiotaWiz or in the arguments file (see *a. Finding a parameter* below).
2. **Otherwise, run the earlier version:**
   - BiotaWiz: under **Settings › Environments › TACTIC Pipeline**, set **Update strategy** to *always use version X*.
     The **Base** step of the wizard then shows `(pinned version)`.
   - Command line: use the image of that version, e.g. `ghcr.io/mpourjam/tactic-pipeline:1.0.0`.

   Run it on a copy of the FASTQ folder without the `<sample>__processed/` folders written by newer versions (on the
   command line, `-fp` also works), so that every sample is processed by that version.
3. **Use the same SILVA release.** The taxonomy also depends on the SILVA release, which depends on when you run, not
   only on the TACTIC version: with internet access, TACTIC 1.0.x always uses the newest SILVA release. To keep an
   older release (command line, TACTIC 1.0.0 or newer), keep only that release in the database folder and run without
   internet access (`docker run --network none ...`); TACTIC then uses the local copy. BiotaWiz has no option for
   this.

### a. Finding a parameter

**In BiotaWiz:** go to the **Parameters** step of the wizard (after **Spike-in**). It has two parts:
- **Basic settings** (open): **Threads** and **Primer trimming**.
- **Advanced settings** (closed; click its title to open it): all other parameters, in sections such as **Merging
  paired-end reads**, **Quality filtering: paired-end reads** or **zOTU picking (all samples together)**.

Find the section named in the change, then the field by its label, e.g. **Minimum relative abundance**. The small grey
name under each label is the parameter's name in the arguments file (e.g. `abund_limit`). Some sections only appear
when they apply: **UPARSE (OTU) arguments**, **Taxonomy-informed clustering (TIC)** and **Taxonomy-agnostic clustering
(TAC)** for that analysis mode, and **zOTU picking per sample** when **Create a zOTU table for each sample** is ticked.

**In the arguments file:**
- BiotaWiz writes `TACTICPipeline_args.yml` into the folder chosen on the **Base** step and rewrites it from the
  **Parameters** step for every run, so change BiotaWiz runs on the **Parameters** step, not in the file.
- Command-line runs use the file given with `-y` (relative to the input folder), or `TACTICPipeline_args.yml` in the input
  folder (created from the default file if it is missing). `-tf` writes a template to start from.
- Every `Analysis_<timestamp>/` folder keeps a copy of the arguments file the run used: look there to see which values
  an earlier run had.

Each section name (e.g. `ClusterZOTUsArgs`) starts a block, and its parameters are indented below it:

```yaml
ClusterZOTUsArgs:
  minsize: 2
  abund_limit: 0.0025
  sample_wise_correction: false
  match_id: 0.99
```

For TACTIC versions before 1.0.0, write decimal parameters with a decimal point (`0.0`, not `0`); otherwise they are
ignored (see 1.0.0).

## [Unreleased]

## [1.0.1] - 2026-11-05

Release candidates: 1.0.1-rc.1, 1.0.1-rc.2, 1.0.1-rc.3

### a. Changed
- **Default zOTU filter of 0.25 %**: the shipped `TACTICPipeline_args.yml` sets `abund_limit: 0.0025` under `ClusterZOTUsArgs` (was `0.000`), following [Reitmeier et al. 2021](https://doi.org/10.1038/s43705-021-00033-z): a zOTU is kept only if it reaches 0.25 % of the reads in at least one sample. The built-in default was already `0.0025`.
- **Spike reads are counted after quality filtering**: spike reads used to be counted in the raw read pairs, before merging and quality filtering, while OTU counts come from the merged and filtered reads. Spike normalization was therefore biased by how many reads of each sample survived QC. Spike reads are now counted in, and removed from, the quality-filtered reads (`usearch_global`, 90 % identity, as in Tourlousse et al. 2017), so they go through the same losses as the OTU reads. Per-spike counts are written to `spike_counts.tsv` in each `<sample>__processed/` folder. **Samples processed by an earlier version are processed again** automatically, as they have no `spike_counts.tsv`; spike-normalized values can differ from earlier runs.
- **Output explanations are HTML pages**: `Output files explanation.html` replaces the `.txt` file in every results folder. It opens in any web browser and has a section menu, a filter box, a dark mode and "Open" links to the files in the folder. The `.txt` file is only written if the page cannot be created.
- **Intermediate files are removed**: `test_*`, `shortened_algn_*`, `filtered_by_abundance_*` and `sorted_length.fasta` (the length-sorted clustering input of TAC and TIC) files no longer remain in the pipeline folders, and they are not used for the PICRUSt2 inputs.
- **Output explanation pages tag each file**: in the OTU, ZOTU, TAC and TIC pages, every file is tagged *Downstream analysis* (the results to analyse) or *Complementary* (files to check, trace or explore them). `zotus.fasta` in `ZOTU_Pipeline/` (the unfiltered denoising output) is now explained. The pages say that `ZOTU-Table-OnlyBacteria-FullTaxonomy.tab` and `ZOTU-Table-All-FullTaxonomy.tab` (TIC) are for analyses at ZOTU level with full taxonomy.

### b. Fixed
- **Spike normalization skipped the first sample**: the first sample column of the OTU table was not checked for spike data. All sample columns are now checked.
- **Missing dependency**: `six` is now installed, as `ete3` (tree handling) needs it.
- **Corrected sample weights are now used in spike normalization**: a processed sample is reused in later runs when its settings and spike amount are unchanged. Its stored weight was reused with it, so a weight corrected in the mapping file never reached `Spike_Normalized-*` tables. The weight and spike amount of reused samples are now updated from the current mapping file.
- **`SampleMinCount_Normalized-*` description**: the output explanation pages described these tables as rarefied (randomly subsampled). They are scaled: each sample's counts are multiplied so that its total equals the smallest sample's total. The pages now say so, and recommend rarefying raw counts downstream for richness comparisons.

### c. Changes that affect your results
Compared with 1.0.0.

**Files not affected** (runs without spike-ins, as long as their samples contain no spike reads):
- all files of BiotaWiz runs: their zOTU filter is set by BiotaWiz (see c.e);
- all files of command-line runs whose arguments file sets its own `abund_limit` under `ClusterZOTUsArgs`;
- the files of the OTU mode (`OTUs-Table.tab`, `OTUs-Seqs.fasta`, `OTUs-Tree-nj.tre`, …) of every run.

**Files affected:**
- command-line runs with the shipped arguments file: every file built from the zOTUs, i.e. all files of the ZOTU,
  TAC and TIC modes and the matching files in `PICRUSt2-Inputs/` (see c.e);
- runs with spike-ins, clearly: the `Spike_Normalized-*.tab` tables of every analysis mode, e.g.
  `Spike_Normalized-ZOTUs-Table.tab`, `Spike_Normalized-SOTUs-Table-TIC.tab`;
- runs with spike-ins, slightly, in samples with spike reads: the count tables (`ZOTUs-Table.tab`, `OTUs-Table.tab`,
  `OTUs-Table-TAC.tab`, `SOTUs-Table-TIC.tab`, `SampleMinCount_Normalized-*.tab`) and the sequence files
  (`*-Seqs*.fasta`).

| Change | Effect | Affects | Reproduce the 1.0.0 result with 1.0.1? |
|---|---|---|---|
| `Spike_Normalized-*` values shift between samples | 3 Medium | BiotaWiz and command line, runs with spike-ins | **No**: run 1.0.0 |
| Counts in samples with spike reads can differ slightly | 2 Low | BiotaWiz and command line, samples with spike reads | **No**: run 1.0.0 |
| No `Spike_Normalized-*` tables when the first sample has no valid spike data | 4 High | BiotaWiz and command line | **No**, and not wanted: the old tables were wrong |
| Corrected sample weights are used | 3 Medium | BiotaWiz and command line, reruns of the same FASTQ folder | **Yes**: put the old weight back |
| zOTU filter at 0.25 % by default | **5 Very high** | Command line with the shipped arguments file; BiotaWiz users through BiotaWiz 1.0.3 | **Yes**: set it to `0` |

#### c.a `Spike_Normalized-*` values shift between samples
Spike reads are now counted after merging and quality filtering, like the reads in the tables. Each sample is scaled
by the mean spike reads of all samples divided by its own spike reads, so the overall level stays about the same, but
samples shift against each other: a sample that lost more reads in merging and filtering than the average sample gets
**higher** values than in 1.0.0; one that lost fewer gets **lower** values. Samples processed by an earlier version
are processed again automatically.

**Effect: 3 Medium.** Only absolute abundances (microbial load per gram) change; relative composition, alpha and
beta diversity do not. Load differences between samples or groups can become smaller or larger, most where samples
lost different shares of reads in merging and filtering (e.g. sequencing runs or batches of different read quality).

**Affects:** BiotaWiz and command-line users, runs with a spike-in standard (**Spike-in** step).

**Reproduce the 1.0.0 result with 1.0.1: No.** There is no parameter for the old counting. Run 1.0.0 (see *How to
reproduce results of an earlier version*). The new values are the corrected ones.

#### c.b Counts in samples with spike reads can differ slightly
Spike reads are now removed from the quality-filtered reads (sequence similarity of at least 90 %) instead of from the
raw reads (bowtie2), in every sample. The number of zOTUs and the read counts can change a little in samples that
contain spike reads: the spiked samples, and samples that picked up spike reads from spiked samples sequenced with
them.

**Effect: 2 Low.** Only reads very similar to the synthetic spike-in sequences are handled differently; real taxa
are hardly affected. At most a zOTU made of spike reads appears or disappears.

**Affects:** BiotaWiz and command-line users, samples with spike reads (normally the spiked samples).

**Reproduce the 1.0.0 result with 1.0.1: No.** Run 1.0.0.

#### c.c No `Spike_Normalized-*` tables when the first sample has no valid spike data
1.0.0 did not check the first sample of the table and still wrote the tables. Now no spike-normalized tables are
written in that case, as for any other sample without valid spike data.

**Effect: 4 High** for the runs it applies to: the old `Spike_Normalized-*` tables contained a sample that cannot be
spike-normalized, so comparisons of absolute abundance with that sample were not valid. Relative composition and
diversity are not affected.

**Affects:** BiotaWiz and command-line users, runs whose first sample had no valid spike data.

**Reproduce the 1.0.0 result with 1.0.1: No, and the old tables should not be used:** they included a sample that
cannot be spike-normalized. Correct that sample's spike data in the mapping file, or leave the sample out.

#### c.d Corrected sample weights are used
If you corrected a sample's weight (`total_weight_in_g` in the mapping file) and the sample was reused from an earlier
run, 1.0.0 still used the old weight. Its spike-normalized values now change in proportion (half the weight gives
double the values).

**Effect: 3 Medium** for the corrected samples: their absolute abundance changes by the same factor as the weight
(half the weight, double the load), which can change comparisons of microbial load between groups. Relative
composition and diversity are not affected.

**Affects:** BiotaWiz and command-line users who rerun the same FASTQ folder after correcting a weight.

**Reproduce the 1.0.0 result with 1.0.1: Yes.** Put the old weight back in the mapping file.

#### c.e zOTU filter at 0.25 % by default
The shipped arguments file `TACTICPipeline_args.yml` now sets **Minimum relative abundance** (`ClusterZOTUsArgs` ›
`abund_limit`) to `0.0025` instead of `0.000`: a zOTU is kept only if it reaches 0.25 % of the reads in at least one
sample, as recommended by [Reitmeier et al. 2021](https://doi.org/10.1038/s43705-021-00033-z) before alpha and beta diversity analyses. This gives fewer
zOTUs, and therefore fewer TAC OTUs and SOTUs, than 1.0.0, and the same filter as 0.7.x.

For BiotaWiz runs, BiotaWiz writes this value itself, for every TACTIC version: BiotaWiz 1.0.3 and newer write
`0.0025`, BiotaWiz 1.0.0 to 1.0.2 write `0` (keep all zOTUs). Update BiotaWiz to get the 0.25 % filter. BiotaWiz
before 1.0.0 only runs TACTIC 0.7.x, which already filtered at 0.25 %.

**Effect: 5 Very high.** The reverse of 1.0.0 (i.a): rare zOTUs are removed, so richness and alpha diversity drop and
beta diversity based on presence/absence (Jaccard, unweighted UniFrac) changes. Measured on 24 mouse samples with
**Minimum number of reads per zOTU** = 2: **186 instead of 493 zOTUs** and **141 instead of 348 SOTUs**.
Abundance-weighted measures (Bray-Curtis, weighted UniFrac) and the dominant taxa change much less.

**Affects:** command-line users with the shipped arguments file. Arguments files that set their own `abund_limit`
keep their value. BiotaWiz users only through BiotaWiz 1.0.3 (see above).

**Reproduce the 1.0.0 result with 1.0.1: Yes.** Set **Minimum relative abundance** to `0`:
- BiotaWiz: **Parameters** step › **Advanced settings** › **zOTU picking (all samples together)** › **Minimum relative
  abundance** = `0`.
- Arguments file: `abund_limit: 0.0` under `ClusterZOTUsArgs`.

## [1.0.0] - 2026-10-23

### a. Taxonomy & SILVA 144+ Support (Dual-Level Resolution)
- **SILVA 144+ Bacterial Kingdom Classification**: Full native compatibility with SILVA release 144+ databases by detecting and handling the formal bacterial Kingdom rank (ending in `-ati` / `-ati--other`, such as *Bacillati*, *Pseudomonadati*) between Domain and Phylum, while cleanly ignoring 8th-level source tags (such as `metagenome`).
- **Dual Taxonomy Outputs**: Simultaneously generates classic 7-level files alongside enriched 8-level files (`*-8Levels.*` across TIC, TAC, ZOTU, OTU, and `-iz`):
  - *Standard Outputs (7 levels)*: Classic `Domain;Phylum;Class;Order;Family;Genus;Species` for downstream tools strictly requiring 7 levels (e.g. PICRUSt2, legacy QIIME / R scripts).
  - *Enriched Outputs (8 levels)*: Full `Domain;Kingdom;Phylum;Class;Order;Family;Genus;Species` for high-resolution profiling, Rhea, Namco, and 8-ring Krona charts.
- **Semicolon Separator Format Standard**: Standardized taxonomy level formatting across all abundance tables and FASTA headers using `;` strictly as a level separator:
  - 7-level outputs contain exactly 6 semicolons (`len(tax.split(";")) == 7`).
  - 8-level outputs contain exactly 7 semicolons (`len(tax.split(";")) == 8`).
  - Correctly preserves multi-word taxon names with spaces (e.g. *Incertae Sedis*, *Candidatus Cryptobacteroides*, *NK4A214 group*).
- **Taxonomy Levels Guide (`README_Taxonomy_Levels.txt`)**: Automatically written into each `Analysis_<timestamp>/` results directory to guide users on when to use standard vs. 8-level files, accompanied by an explanatory summary banner in `Analysis_log.txt`.

### b. Pre-formatted PICRUSt2 Inputs (`PICRUSt2-Inputs/`)
- **Ready-to-Use Functional Predictions**: Added an automated `PICRUSt2-Inputs/` folder containing pre-formatted abundance tables (taxonomy column removed) and matching single-line FASTA files (pure sequence IDs, all header annotations stripped).
- **Direct CLI Compatibility**: Fully compatible out-of-the-box with `picrust2_pipeline.py -s <Seqs.fasta> -i <Table.tab> -o picrust2_out` without manual reformatting or BIOM conversions. Includes dedicated guidance in `PICRUSt2-Inputs/Output files explanation.txt`.

### c. Plain-Language Output Documentation (`Output files explanation.txt`)
- **Folder-Specific User Guides**: Automatically writes a clear, plain-language `Output files explanation.txt` guide inside the root `Analysis_<timestamp>/` results folder and within each active pipeline directory (`OTU_Pipeline`, `ZOTU_Pipeline`, `TAC_Pipeline`, `TIC_Pipeline`, `PICRUSt2-Inputs`, and `Shared`).
- **Accessible Guidance**: Written in simple language explaining what each output file contains, how to open it (e.g. interactive Krona HTML charts in a web browser, abundance tables in Excel or R), and its recommended scientific purpose.
- **TIC Full-Taxonomy Clarification**: Clearly explains the purpose of `ZOTU-Table-OnlyBacteria-FullTaxonomy.tab` (default for bacterial diversity and functional predictions) versus `ZOTU-Table-All-FullTaxonomy.tab` (broad-domain studies and contamination auditing).

### d. Clustering & Analysis Refinements
- **Multi-Analysis Mode Selector**: Supports running multiple clustering pipelines in a single command (e.g., `--analysis-mode OTU TAC TIC`), eliminating redundant preprocessing and sharing denoised sequences.
- **Clean Full-Taxonomy ZOTU Tables**: When no non-bacterial sequences are detected in TIC mode, the redundant duplicate `ZOTU-Table-All-FullTaxonomy.tab` is automatically omitted, leaving a clean single `ZOTU-Table-OnlyBacteria-FullTaxonomy.tab`. When non-bacterial sequences are present, both tables are retained.
- **Guaranteed Non-Bacterial Sequence Logging**: `Invalid-Tax-Seqs-TIC.fasta` (and its `-8Levels` counterpart) is consistently retained across all runs for auditing and downstream test compatibility.

### e. FASTQ Detection & Pairing
- **Robust Read Direction & Marker Matching**: Read direction is determined by the rightmost marker in the filename, supporting standard `R1`/`R2` (and `r1`/`r2`), `read1`/`read2`, `forward`/`reverse`, SRA/ENA `_1`/`_2` directly before the extension, and `F`/`R` tags.
- **Support for Hyphenated & Numbered Sample IDs**: Sample names containing digits and hyphens (such as `FU-1_..._R2_001`) are accurately classified without false forward pairings.
- **Clean Sample IDs in Mapping Files**: Supports clean sample IDs (e.g. `SampleA`) in mapping files, correctly resolving sequencer filenames like `SampleA_S388_L001_R1_001.fastq.gz` without prefix collisions (e.g., `SampleA` vs `SampleA0`).
- **Targeted Sample Lookups**: When a mapping file is provided, only requested samples are gathered, dramatically speeding up startup on large datasets.
- **Consistent Reverse Read Validation**: Re-running inside an existing project directory automatically skips prior `Analysis_*` and `*__processed` folders. Corrupted reverse mates safely fall back to single-end processing with a clear warning.
- **Single-End & Paired-End Support**: Single-end reads are automatically detected and processed alongside paired-end libraries. Unsupported `.zst` files are ignored.

### f. Database Reliability & Offline Support
- **Offline Fallback for Local SILVA Database**: Automatically falls back to local SILVA database files when remote download endpoints or Cloudflare R2 worker connections are inaccessible or in air-gapped environments.
- **SINA Binary Shared Library Link Restoration**: Fixed path resolution for the SINA binary and `libsina.so.0` in containerized Docker environments.
- **Automatic Permission Normalization**: Downloaded reference databases automatically set appropriate file permissions immediately after download for shared and multi-user environments.
- **Cloudflare SILVA Downloads Worker**: Upgraded `sharp` dependency to `>=0.35.5` (and updated `undici`) in `cloudflare/silva-downloads-worker/package-lock.json` to resolve security advisories.

### g. Settings, Logging & Reproducibility
- **Updated Pipeline Defaults**: Defaulted `abund_limit` to `0.0` and `sample_wise_correction` to `false` for TAC and TIC clustering modes in `TACTICPipeline_args.yml`, with safe fallbacks when omitted from custom configs.
- **Versioned Log Archiving**: `Pipeline_log.txt` is stamped with the pipeline version and automatically archived inside each run's `Analysis_<timestamp>/` results directory.
- **Deterministic Runs**: Sequence pairs and samples are processed in a consistent, reproducible order across different platforms and filesystems.
- **Paper Benchmarks & DOI References**: Updated step-wise and length-based 16S rRNA diversity survey benchmark datasets and documentation accompanying the TACTIC paper, including Zenodo DOI references.

### h. Fixed
- **Reruns picked up earlier results as samples**: FASTQ files inside `<sample>__processed/` folders of an earlier run could be treated as new samples, with numbered names, again on every rerun. They are now skipped, and read pairs are processed in a fixed order.
- **Missing reverse reads and wrong `parent_path`**: a warning is logged when the reverse read file of a sample is missing or the `parent_path` in the mapping file does not exist.
- **Pipeline version in the logs**: the version is resolved from the Docker image (`GIT_TAG`, `VERSION` file) and written as the first line of `Pipeline_log.txt` and `Analysis_log.txt`.

### i. Changes that affect your results
Compared with 0.7.6.

**Files not affected:** the files of the OTU mode (de novo OTUs in `OTU_Pipeline/`: `OTUs-Table.tab`,
`OTUs-Seqs.fasta`, `OTUs-Tree-nj.tre`, `OTUs-Krona.html` and their `SampleMinCount_Normalized-*` /
`Spike_Normalized-*` tables), as long as the run used SILVA 138 or older, standard FASTQ file names, a FASTQ folder
without earlier results and BiotaWiz or the shipped arguments file.

**Files affected:**
- BiotaWiz runs: every file built from the zOTUs, i.e. all files of the ZOTU (`ZOTUs-*`), TAC (`OTUs-*-TAC*`,
  `Map-ZOTUs-OTUs-TAC.tab`) and TIC (`SOTUs-*`, `Map-*.tab`, `ZOTU-Table-*-FullTaxonomy.tab`,
  `Invalid-Tax-Seqs-TIC.fasta`) modes, and the matching files in `PICRUSt2-Inputs/`;
- runs with SILVA 144 or newer: the taxonomy in every table, sequence file and Krona chart, and the TIC files;
- runs with non-standard FASTQ file names, reruns in a folder with earlier results, and command-line runs without
  primer trimming parameters: all files.

| Change | Effect | Affects | Reproduce the 0.7.6 result with 1.0.0? |
|---|---|---|---|
| More zOTUs, TAC OTUs and SOTUs: **Minimum relative abundance (`abund_limit`)** `0` is now used | **5 Very high** | **All BiotaWiz users**; command line only with whole numbers for decimal parameters | **Yes**: set it to `0.0025` |
| Primer trimming defaults of 17 and 21 bases instead of 5 and 5 | 3 Medium | Command line only, arguments files without these parameters | **Yes**: set `5` and `5` |
| Taxonomy columns with SILVA 144 or newer | **5 Very high** | BiotaWiz and command line | **No**, and not wanted: the old ranks were in the wrong columns |
| Which FASTQ files are paired | 4 High | BiotaWiz and command line, non-standard file names | **No** for wrong pairings; rename the files |
| Reruns no longer pick up earlier results as samples | 4 High | BiotaWiz and command line, reruns | **No**, and not needed: drop the extra samples |

#### i.a More zOTUs, TAC OTUs and SOTUs for runs from BiotaWiz
BiotaWiz sets **Minimum relative abundance** (`ClusterZOTUsArgs` › `abund_limit`) to `0` by default, meaning "keep all
zOTUs", and writes it as the whole number `0`. Up to 0.7.6, a whole number given for a decimal parameter was silently
ignored and the default `0.0025` was used: every zOTU that did not reach more than 0.25 % of the reads in at least one
sample was removed. From 1.0.0 the value is used, so these low-abundance zOTUs are kept. Expect clearly **more
zOTUs**, and therefore more TAC OTUs and SOTUs (which are built from the zOTUs), higher richness and more **rare**
taxa. BiotaWiz 1.0.3 and newer write `0.0025` again, so BiotaWiz users who update get the 0.7.x filter with 1.0.0 as
well (see 1.0.1, c.e).

**Effect: 5 Very high.** Many rare zOTUs are now kept, so richness and alpha diversity rise strongly and beta
diversity based on presence/absence (Jaccard, unweighted UniFrac) changes. Measured on 24 mouse samples with
**Minimum number of reads per zOTU** = 2 (the BiotaWiz default): **493 instead of 186 zOTUs** (2.7 times as many)
and **348 instead of 141 SOTUs** (2.5 times). Abundance-weighted measures (Bray-Curtis, weighted UniFrac) and the
dominant taxa change much less, as each removed zOTU had at most 0.25 % of the reads in every sample.

**Affects:** **all BiotaWiz users**. Command-line users only if their arguments file has whole numbers for decimal
parameters (e.g. `abund_limit: 0`); the shipped `TACTICPipeline_args.yml` (`abund_limit: 0.000`) is not affected.

**Reproduce the 0.7.6 result with 1.0.0: Yes.** Set **Minimum relative abundance (`abund_limit`)** to `0.0025`; 1.0.0 then removes the
same zOTUs as 0.7.6:
- BiotaWiz: **Parameters** step › **Advanced settings** › **zOTU picking (all samples together)** › **Minimum relative
  abundance** = `0.0025`.
- Arguments file: `abund_limit: 0.0025` under `ClusterZOTUsArgs`.

The same applies to the other decimal parameters. If you gave one of them a whole number such as `0` or `1`, 0.7.6
ignored it and used the value below; enter that value to reproduce the 0.7.6 result:

| Parameter in BiotaWiz (**Parameters** step › **Advanced settings**) | Arguments file | Value 0.7.6 used |
|---|---|---|
| Quality filtering: paired-end reads › **Maximum expected errors per base** | `FilterBothSidesArgs` › `fastq_maxee_rate` | `0.002` |
| Quality filtering: single-end reads › **Maximum expected errors per base** | `FilterOneSideArgs` › `fastq_maxee_rate` | `0.002` |
| zOTU picking (all samples together) › **Minimum identity for mapping reads to zOTUs** | `ClusterZOTUsArgs` › `match_id` | `0.99` |
| UPARSE (OTU) arguments › **Minimum relative abundance** | `DeNovoClusterOTUsArgs` › `abund_limit` | `0.0025` |
| Taxonomy-informed clustering (TIC) › **Family**, **Genus**, **Species similarity threshold** | `ComplexTICArgs` › `family_sim`, `genus_sim`, `species_sim` | `0.95`, `0.97`, `0.987` (not the `0.90`, `0.95`, `0.987` that BiotaWiz and the shipped arguments file set) |
| Taxonomy-agnostic clustering (TAC) › **Clustering similarity threshold** | `TACClusterArgs` › `cluster_thr` | `0.987` |
| not in BiotaWiz | `Filter16SArgs` › `e` | `0.1` |
| not in BiotaWiz | `BuildZOTUTableArgs` › `match_id` | `0.99` |

#### i.b Primer trimming defaults of 17 and 21 bases
If an arguments file does not set **Bases to trim from the start of merged reads** (`stripleft`) and **Bases to trim
from the end of merged reads** (`stripright`) under `TrimBothSidesArgs` / `TrimOneSideArgs`, the defaults are now 17
and 21 bases (341F/785R primers) instead of 5 and 5. Sequences are shorter and no longer start with primer remnants,
which can give **fewer zOTUs**.

**Effect: 3 Medium.** zOTUs that differed only in leftover primer bases, which come from the primers and not from
the bacteria, are merged. zOTU richness can drop; the taxonomic composition stays the same.

**Affects:** command-line users with their own arguments file without these parameters. Not BiotaWiz users and not the
shipped arguments file: both set 17 and 21.

**Reproduce the 0.7.6 result with 1.0.0: Yes.** Set `stripleft: 5` and `stripright: 5` under `TrimBothSidesArgs`, and
`stripleft: 5` under `TrimOneSideArgs` in the `TACTICPipeline_args.yml` (Only command-line users).

#### i.c Taxonomy columns with SILVA 144 or newer
SILVA 144 added a kingdom rank (e.g. *Bacillati*) between domain and phylum. 0.7.6 did not know it: the kingdom ended
up in the phylum column and every lower rank one column to the right (the genus in the species column). 1.0.0 puts
each rank in its own column and also writes `*-8Levels.*` files with the kingdom. Because TIC clusters by taxonomy,
the SOTUs (`SOTUs-Table-TIC.tab`, `SOTUs-Seqs-TIC.fasta`) differ as well. With SILVA 138 or older, the taxonomy is the
same as in 0.7.6.

**Effect: 5 Very high.** In 0.7.6 every rank below the domain was in the wrong column (the kingdom as phylum, the
genus as species), so the taxonomic composition at every rank was wrong. TIC also formed different SOTUs, which
changes richness and diversity in TIC mode.

**Affects:** BiotaWiz and command-line users whose runs used **SILVA 144** or newer.

**Reproduce the 0.7.6 result with 1.0.0: No, and it should not be used:** its ranks were in the wrong columns.

#### i.d Which FASTQ files are paired
Read direction is now taken from the last marker in the file name, and bare numbers only count as read numbers as
`_1`/`_2` directly before the extension. Samples that earlier versions paired wrongly, missed or treated as single-end
can now be paired (or the other way round). Check `Pipeline_log.txt` and `Missed_samples.txt` to see that every sample
was found and paired as expected.

**Effect: 4 High** for the samples it applies to: a missing sample changes the group comparisons, and a sample
analysed from single reads instead of merged pairs gets shorter, less accurate sequences and other zOTUs, so it can
look different from the rest of its group. Samples with standard file names do not change.

**Affects:** BiotaWiz and command-line users whose file names are not in a standard pattern (e.g. `_R1_001`/`_R2_001`).

**Reproduce the 0.7.6 result with 1.0.0: No** for wrong pairings. If a sample is now missing or single-end, rename its
files to a supported pattern (listed in the README).

#### i.e Reruns no longer pick up earlier results as samples
FASTQ files inside the `<sample>__processed/` folders of an earlier run could be picked up again as extra samples with
numbered names only if a run was stopped or failed at early stages. They are now skipped.

**Effect: 4 High** for the runs it applies to: copies of samples count as extra replicates, which shifts group
averages and makes group differences look more significant than they are. The composition of each real sample
does not change.

**Affects:** BiotaWiz and command-line users who rerun the same FASTQ folder.

**Reproduce the 0.7.6 result with 1.0.0: No, and not needed:** the extra samples were copies of real samples. Remove
them from the old tables before comparing. Or start from a fresh folder with your raw FASTQs. 

## [0.7.6] - 2026-09-22

### a. Added
- **Several analysis modes in one run**: `--analysis-mode`/`-am` takes several modes, e.g. `-am OTU TAC TIC`. Preprocessing and the denoised ZOTUs are shared between the modes.
- **Prebuilt SILVA indexes**: the SILVA database and its SINA index are downloaded prebuilt and checksum-verified from the TACTIC file server (Cloudflare R2). If that fails, the pipeline falls back to the SILVA FTP server and builds the index itself.

### b. Changed
- **`de-novo` analysis mode**: now an alias of `OTU`, so older BiotaWiz versions that send `de-novo` keep working.
- **Less output**: the list of available SILVA versions is no longer printed.

### c. Fixed
- **Memory use when preparing SILVA**: the ARB file is extracted in chunks and replaced only when complete, an extracted file that is still current is reused, and incomplete indexes are removed after a failure. The number of SINA indexing threads can be set.
- **Tree leaf labels**: labels in the `.tre` files no longer carry extra quotes.

### d. Changes that affect your results
Compared with 0.7.5.

**Files not affected:** all tables (`*.tab`), sequence files (`*.fasta`) and Krona charts (`*-Krona.html`).

**Files affected:** only the trees (`*-Tree-nj*.tre`, e.g. `SOTUs-Tree-nj-TIC.tre`), and only their leaf labels.

| Change | Effect | Affects | Reproduce the 0.7.5 result with 0.7.6? |
|---|---|---|---|
| Tree leaf labels without extra quotes | 1 Very low | BiotaWiz and command line, users of the trees | **Not needed**: same results wherever the labels matched |

#### d.a Tree-based results (UniFrac, phylogenetic diversity)
The leaf labels in the `.tre` files had extra quotes, so tools could fail to match tree leaves to the table rows.
Results computed with the trees can differ now that the labels match. Tree shapes, branch lengths, counts and taxonomy
are not affected.

**Effect: 1 Very low.** Tree shapes and branch lengths did not change, only the labels. UniFrac and phylogenetic
diversity are the same wherever the tool matched the labels; where it did not, the tool could not use the tree at
all, rather than giving other numbers.

**Affects:** BiotaWiz and command-line users who use the trees.

**Reproduce the 0.7.5 result with 0.7.6: Not needed.** Where the old labels matched, the results are the same; where
they did not, the new results are the correct ones.

## [0.7.5] - 2026-07-15

### a. Changed
- **Choosing the newest SILVA release**: decided by the version number; the timestamp of the `.arb.gz` file only breaks ties.
- **Argument and log files kept with the results**: the arguments file used (any `*.yml`) and all `*log.txt` files are copied into, and kept in, the `Analysis_<timestamp>/` folder.
- **README**: explains how to mount a database folder as a Docker volume.

### b. Changes that affect your results
Compared with 0.7.4.

**Files not affected:** the counts and sequences of the ZOTU, OTU and TAC files (`ZOTUs-Table.tab`,
`OTUs-Table.tab`, `OTUs-Table-TAC.tab`, `ZOTUs-Seqs.fasta`, `OTUs-Seqs*.fasta`); only their taxonomy can change.

**Files affected (only if the database folder holds several SILVA releases):** the taxonomy in all tables, sequence
files and Krona charts; the TIC files (`SOTUs-Table-TIC.tab`, `SOTUs-Seqs-TIC.fasta`, `Map-*.tab`), because TIC
clusters by taxonomy; and the trees (`*-Tree-nj*.tre`), which are built from an alignment against SILVA.

| Change | Effect | Affects | Reproduce the 0.7.4 result with newer versions? |
|---|---|---|---|
| Newest SILVA release chosen by version number | 3 Medium | BiotaWiz and command line, database folders with several SILVA releases | **Yes**: keep only the release you want |

#### b.a Taxonomy, if the database folder holds several SILVA releases
The newest release is now chosen by its version number instead of the file date, so another release than before may
be used. Taxon names can then differ, and with them the SOTUs of TIC, which clusters by taxonomy. The counts of
zOTUs, OTUs and TAC OTUs are not affected.

**Effect: 3 Medium.** Taxa that were renamed or reclassified between the two SILVA releases get other names, so the
composition can shift for those taxa, and TIC can form somewhat different SOTUs. Most taxa keep their names;
diversity based on zOTUs, OTUs and TAC OTUs does not change.

**Affects:** BiotaWiz and command-line users whose database folder holds more than one SILVA release.

**Reproduce the 0.7.4 result with newer versions: Yes.** Keep only the release you want in the database folder (see
*How to reproduce results of an earlier version*, step 3).

## [0.7.4] - 2025-12-19

### a. Added
- **Quickstart** section with a ready-to-run command in the README.

### b. Changed
- **`-tf` (template files)**: also writes the arguments file template (`*_TEMPLATE.yml`) into the FASTQ folder, next to the mapping file template.
- **Database folder (`-db`)**: an absolute path is used as given; a relative path is taken relative to the input folder (`-i`). The folder is created if it does not exist.

### c. Fixed
- **Interrupted SILVA indexing**: a signal file (`.sidx.success`) is written after the SINA index was built successfully. If it is missing or invalid, the index is built again instead of using an incomplete one.

### d. Changes that affect your results
Compared with 0.7.3.

**Files not affected:** the counts and sequences of the ZOTU, OTU and TAC files (`ZOTUs-Table.tab`,
`OTUs-Table.tab`, `OTUs-Table-TAC.tab`, `ZOTUs-Seqs.fasta`, `OTUs-Seqs*.fasta`).

**Files affected (only runs whose SILVA indexing was interrupted):** the taxonomy in all tables, sequence files and
Krona charts; the TIC files (`SOTUs-Table-TIC.tab`, `SOTUs-Seqs-TIC.fasta`, `Map-*.tab`); and the trees
(`*-Tree-nj*.tre`).

| Change | Effect | Affects | Reproduce the 0.7.3 result with newer versions? |
|---|---|---|---|
| Incomplete SILVA index is rebuilt | **5 Very high** | BiotaWiz and command line, interrupted database setups | **No**, and not wanted: the old taxonomy was incomplete |

#### d.a Taxonomy after an interrupted database setup
An incomplete SILVA index (e.g. after a crash during indexing) is now rebuilt instead of used. Runs that used an
incomplete index could have missing or wrong taxonomy, and therefore different SOTUs (TIC) and trees.

**Effect: 5 Very high** for the runs it applies to: the taxonomy could be missing or wrong for many sequences, so
the taxonomic composition was wrong, and so were the TIC SOTUs and with them richness and diversity in TIC mode.

**Affects:** BiotaWiz and command-line users whose SILVA indexing was interrupted.

**Reproduce the 0.7.3 result with newer versions: No, and it should not be used:** its taxonomy was incomplete.
