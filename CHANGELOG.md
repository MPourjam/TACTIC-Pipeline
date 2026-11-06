# Changelog

All notable changes to the TACTIC pipeline will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Each version ends with **Changes that affect your results**: what can make the numbers of zOTUs, OTUs or SOTUs,
the counts, the taxonomy or the normalized tables differ from the previous version, graded by **Effect** (1 Very low
to 5 Very high). How to read these sections, what the grades mean and how to reproduce results of an earlier version:
[Comparing results of different versions](https://tum-core-facility-microbiome.github.io/BiotaWiz/docs/getting-started/compare-versions).

A short, plain-language summary of each version is in [CHANGELOG_SIMPLE.md](CHANGELOG_SIMPLE.md).

## [Unreleased]

## [1.0.0] - 2026-11-06

**Re-released on 2026-11-06.** The first build of 1.0.0, released on 2026-10-23, was withdrawn. This build contains
all of its changes plus the fixes tested as 1.0.1-rc.1 to 1.0.1-rc.4. Both builds write `v1.0.0` in the first line of
`Analysis_log.txt`. Runs from before 2026-11-06 used the first build (see j.k). Later runs used it too wherever its
image was still installed: on the command line until `docker pull ghcr.io/mpourjam/tactic-pipeline:1.0.0`, and in
BiotaWiz before 1.0.3. BiotaWiz 1.0.3 and newer replace it automatically.

Release candidates (tested as 1.0.1): 1.0.1-rc.1, 1.0.1-rc.2, 1.0.1-rc.3, 1.0.1-rc.4

### a. Taxonomy & SILVA 144+ Support (Dual-Level Resolution)
- **SILVA 144+ Bacterial Kingdom Classification**: Full native compatibility with SILVA release 144+ databases by detecting and handling the formal bacterial Kingdom rank (ending in `-ati` / `-ati--other`, such as *Bacillati*, *Pseudomonadati*) between Domain and Phylum, while cleanly ignoring 8th-level source tags (such as `metagenome`).
- **Dual Taxonomy Outputs**: With SILVA 144 or newer, generates classic 7-level files alongside enriched 8-level files (`*-8Levels.*` across TIC, TAC, ZOTU, OTU, and `-iz`). With SILVA 138 or older, only the 7-level files are written, in the same format as in 0.7.6:
  - *Standard Outputs (7 levels)*: Classic `Domain;Phylum;Class;Order;Family;Genus;Species` for downstream tools strictly requiring 7 levels (e.g. PICRUSt2, legacy QIIME / R scripts).
  - *Enriched Outputs (8 levels)*: Full `Domain;Kingdom;Phylum;Class;Order;Family;Genus;Species` for high-resolution profiling, Rhea, Namco, and 8-ring Krona charts.
- **Semicolon Separator Format Standard**: Standardized taxonomy level formatting across all abundance tables and FASTA headers using `;` strictly as a level separator:
  - 7-level outputs contain exactly 6 semicolons (`len(tax.split(";")) == 7`).
  - 8-level outputs contain exactly 7 semicolons (`len(tax.split(";")) == 8`).
  - Correctly preserves multi-word taxon names with spaces (e.g. *Incertae Sedis*, *Candidatus Cryptobacteroides*, *NK4A214 group*).
- **Taxonomy Levels Guide (`README_Taxonomy_Levels.txt`)**: With SILVA 144 or newer, written into the `Analysis_<timestamp>/` results directory to guide users on when to use standard vs. 8-level files, accompanied by an explanatory summary banner in `Analysis_log.txt`.

### b. Pre-formatted PICRUSt2 Inputs (`PICRUSt2-Inputs/`)
- **Ready-to-Use Functional Predictions**: Added an automated `PICRUSt2-Inputs/` folder containing pre-formatted abundance tables (taxonomy column removed) and matching single-line FASTA files (pure sequence IDs, all header annotations stripped). Intermediate files are not used for them.
- **Direct CLI Compatibility**: Fully compatible out-of-the-box with `picrust2_pipeline.py -s <Seqs.fasta> -i <Table.tab> -o picrust2_out` without manual reformatting or BIOM conversions. Includes dedicated guidance in `PICRUSt2-Inputs/Output files explanation.html`.

### c. Output explanation pages (`Output files explanation.html`)
- **One page per results folder**: an `Output files explanation.html` page is written in the `Analysis_<timestamp>/` folder and in each pipeline folder (`OTU_Pipeline`, `ZOTU_Pipeline`, `TAC_Pipeline`, `TIC_Pipeline`, `PICRUSt2-Inputs`, `Shared`). It opens in any web browser and has a section menu, a filter box, a dark mode and "Open" links to the files in the folder. A `.txt` file is only written if the page cannot be created.
- **Plain language**: each file is explained with what it contains, how to open it (e.g. Krona charts in a web browser, tables in Excel or R) and what it is for.
- **Files are tagged**: in the OTU, ZOTU, TAC and TIC pages, every file is tagged *Downstream analysis* (the results to analyse) or *Complementary* (files to check, trace or explore them). `zotus.fasta` in `ZOTU_Pipeline/` (the unfiltered denoising output) is explained.
- **TIC full-taxonomy tables**: `ZOTU-Table-OnlyBacteria-FullTaxonomy.tab` and `ZOTU-Table-All-FullTaxonomy.tab` are explained as the tables for analyses at ZOTU level with full taxonomy (bacteria only, or all domains for broad-domain studies and contamination checks).
- **`SampleMinCount_Normalized-*` tables** are described as scaled (each sample's counts multiplied so that its total equals the smallest sample's total), not rarefied, with the advice to rarefy raw counts downstream for richness comparisons.
- **Intermediate files are removed**: `test_*`, `shortened_algn_*`, `filtered_by_abundance_*` and `sorted_length.fasta` (the length-sorted clustering input of TAC and TIC) no longer remain in the pipeline folders.

### d. Clustering & Analysis Refinements
- **Default zOTU filter of 0.25 %**: the shipped `TACTICPipeline_args.yml` sets `abund_limit: 0.0025` under `ClusterZOTUsArgs` (was `0.000`), following [Reitmeier et al. 2021](https://doi.org/10.1038/s43705-021-00033-z): a zOTU is kept only if it has more than 0.25 % of the reads in at least one sample. The built-in default was already `0.0025`.
- **Clean Full-Taxonomy ZOTU Tables**: When no non-bacterial sequences are detected in TIC mode, the redundant duplicate `ZOTU-Table-All-FullTaxonomy.tab` is automatically omitted, leaving a clean single `ZOTU-Table-OnlyBacteria-FullTaxonomy.tab`. When non-bacterial sequences are present, both tables are retained.
- **Non-Bacterial Sequence Logging**: with SILVA 144 or newer, `Invalid-Tax-Seqs-TIC.fasta` also gets a `-8Levels` counterpart.

### e. Spike-in normalization
- **Spike reads are counted after quality filtering**: spike reads used to be counted in the raw read pairs, before merging and quality filtering, while OTU counts come from the merged and filtered reads. Spike normalization was therefore biased by how many reads of each sample survived QC. Spike reads are now counted in, and removed from, the quality-filtered reads (`usearch_global`, 90 % identity, as in Tourlousse et al. 2017), so they go through the same losses as the OTU reads. Per-spike counts are written to `spike_counts.tsv` in each `<sample>__processed/` folder. **Samples processed by an earlier version are processed again** automatically, as they have no `spike_counts.tsv`.

### f. FASTQ Detection & Pairing
- **Robust Read Direction & Marker Matching**: Read direction is determined by the rightmost marker in the filename, supporting standard `R1`/`R2` (and `r1`/`r2`, after `_`, `-` or `.`), `read1`/`read2`, `forward`/`reverse`, SRA/ENA `_1`/`_2` directly before the extension, and `F`/`R` tags directly before the extension or before `_001`. `_F_`/`_R_` elsewhere in the name no longer count as markers.
- **Support for Hyphenated & Numbered Sample IDs**: Sample names containing digits and hyphens (such as `FU-1_..._R2_001`) are accurately classified without false forward pairings.
- **Clean Sample IDs in Mapping Files**: Supports clean sample IDs (e.g. `SampleA`) in mapping files, correctly resolving sequencer filenames like `SampleA_S388_L001_R1_001.fastq.gz` without prefix collisions (e.g., `SampleA` vs `SampleA0`).
- **Targeted Sample Lookups**: When a mapping file is provided, only requested samples are gathered, dramatically speeding up startup on large datasets.
- **Consistent Reverse Read Validation**: Re-running inside an existing project directory automatically skips prior `Analysis_*` and `*__processed` folders. Corrupted reverse mates safely fall back to single-end processing with a clear warning.
- **Single-End Files Without a Read Marker**: without a mapping file, FASTQ files without any read marker are now processed as single-end samples (0.7.6 ignored them; forward files without a reverse file were already processed as single-end). Unsupported `.zst` files are ignored.

### g. Database Reliability & Offline Support
- **Offline Fallback for Local SILVA Database**: Automatically falls back to local SILVA database files when remote download endpoints or Cloudflare R2 worker connections are inaccessible or in air-gapped environments.
- **SINA Binary Shared Library Link Restoration**: Fixed path resolution for the SINA binary and `libsina.so.0` in containerized Docker environments.
- **Automatic Permission Normalization**: Downloaded reference databases automatically set appropriate file permissions immediately after download for shared and multi-user environments.
- **Cloudflare SILVA Downloads Worker**: Upgraded `sharp` dependency to `>=0.35.5` (and updated `undici`) in `cloudflare/silva-downloads-worker/package-lock.json` to resolve security advisories.

### h. Settings, Logging & Reproducibility
- **Settings are applied as written**: a whole number given for a decimal parameter (e.g. `abund_limit: 0`) is now used. Up to 0.7.6 it was silently ignored and the built-in default was used.
- **Shipped arguments file**: `abund_limit: 0.0025` under `ClusterZOTUsArgs` (see d); `abund_limit: 0.0` and `sample_wise_correction: false` under `ComplexTICArgs` and `TACClusterArgs` (not used yet), with safe fallbacks when omitted from custom configs.
- **Versioned Log Archiving**: `Pipeline_log.txt` is stamped with the pipeline version and automatically archived inside each run's `Analysis_<timestamp>/` results directory.
- **Deterministic Runs**: Sequence pairs and samples are processed in a consistent, reproducible order across different platforms and filesystems.
- **Paper Benchmarks & DOI References**: Updated step-wise and length-based 16S rRNA diversity survey benchmark datasets and documentation accompanying the TACTIC paper, including Zenodo DOI references.

### i. Fixed
- **Reruns picked up earlier results as samples**: FASTQ files inside `<sample>__processed/` folders of an earlier run could be treated as new samples, with numbered names, again on every rerun. They are now skipped, and read pairs are processed in a fixed order.
- **Missing reverse reads and wrong `parent_path`**: a warning is logged when the reverse read file of a sample is missing or the `parent_path` in the mapping file does not exist.
- **Pipeline version in the logs**: the version is resolved from the Docker image (`GIT_TAG`, `VERSION` file) and written as the first line of `Pipeline_log.txt` and `Analysis_log.txt`.
- **Missing normalized tables**: if a processed sample had no column in a table (no reads matching a zOTU, e.g. a blank), no `SampleMinCount_Normalized-*` and no `Spike_Normalized-*` table was written for that table. They are now written for the samples in the table.
- **Unchecked first sample in normalization**: the check that every sample of a table is in the spike stats file skipped the first sample. If it was missing, its raw counts were written unchanged into the normalized tables. The run now stops with an error, as it already did for the other samples.
- **Corrected sample weights are now used in spike normalization**: a processed sample is reused in later runs when its settings and spike amount are unchanged. Its stored weight was reused with it, so a weight corrected in the mapping file never reached `Spike_Normalized-*` tables. The weight and spike amount of reused samples are now updated from the current mapping file.
- **Missing dependency**: `six` is now installed, as `ete3` (tree handling) needs it.

### j. Changes that affect your results
Compared with 0.7.6. For runs with the first 1.0.0 build (2026-10-23 to 2026-11-05), see j.k.

**Files not affected:** the values in all files of a run that used
- BiotaWiz 1.0.3 or newer, or an arguments file that sets its own `abund_limit` under `ClusterZOTUsArgs` and the
  primer trimming values (`stripleft`, `stripright`) (on the command line, a FASTQ folder already analysed with
  0.7.6 keeps its arguments file, see j.a),
- no spike-ins, and no spike reads in the samples,
- *SILVA 138 or older*,
- standard FASTQ file names and a FASTQ folder without earlier results.

Even then, a few files are new or no longer written: `ZOTU-Table-All-FullTaxonomy.tab` (TIC) is not written when
there are no non-bacterial sequences, as it would be identical to `ZOTU-Table-OnlyBacteria-FullTaxonomy.tab`; the
intermediate files `test_*`, `shortened_algn_*`, `filtered_by_abundance_*` and `sorted_length.fasta` are removed;
normalized tables that 0.7.6 left out are written (j.f); `PICRUSt2-Inputs/` and the explanation pages are new.

The files of the OTU mode (de novo OTUs in `OTU_Pipeline/`: `OTUs-Table.tab`, `OTUs-Seqs.fasta`, `OTUs-Tree-nj.tre`,
`OTUs-Krona.html` and their `SampleMinCount_Normalized-*` tables) are also not affected by the zOTU filter (j.a, j.b).

**Files affected:**
- command-line runs with the shipped arguments file in a new FASTQ folder, and BiotaWiz 1.0.0 to 1.0.2 runs: every
  file built from the zOTUs, i.e. all files of the ZOTU (`ZOTUs-*`), TAC (`OTUs-*-TAC*`, `Map-ZOTUs-OTUs-TAC.tab`)
  and TIC (`SOTUs-*`, `Map-*.tab`, `ZOTU-Table-*-FullTaxonomy.tab`, `Invalid-Tax-Seqs-TIC.fasta`) modes, and the
  matching files in `PICRUSt2-Inputs/` (j.a, j.b);
- runs with spike-ins, clearly: the `Spike_Normalized-*.tab` tables of every analysis mode; slightly, in samples with
  spike reads: the count tables and sequence files (j.d, j.e, j.g);
- runs in which a sample has no column in a table: the `SampleMinCount_Normalized-*.tab` and `Spike_Normalized-*.tab`
  tables, which 0.7.6 did not write (j.f);
- runs with SILVA 144 or newer: the taxonomy in every table, sequence file and Krona chart, and the TIC files (j.h);
- runs with non-standard FASTQ file names, reruns in a folder with earlier results, and command-line runs without
  primer trimming parameters: all files (j.c, j.i, j.j).

| Change | Effect | Affects | Reproduce the 0.7.6 result with 1.0.0? |
|---|---|---|---|
| j.a Fewer zOTUs, TAC OTUs and SOTUs: zOTU filter at 0.25 % by default | **5 Very high** | Command line with the shipped arguments file in a new FASTQ folder: ZOTU, TAC, TIC files | **Yes**: set it to `0.0` |
| j.b More zOTUs, TAC OTUs and SOTUs: **Minimum relative abundance** `0` is now used | **5 Very high** | BiotaWiz 1.0.0–1.0.2 (command line only with whole numbers for decimal parameters): ZOTU, TAC, TIC files | **Yes**: update BiotaWiz, or set it to `0.0025` |
| j.c Primer trimming defaults of 17 and 21 bases instead of 5 and 5 | 3 Medium | Command line only, own arguments files without these parameters: all files | **Yes**: set `5` and `5` |
| j.d `Spike_Normalized-*` values shift between samples | 3 Medium | BiotaWiz and command line: runs with spike-ins, if you used `Spike_Normalized-*` | **No**: run 0.7.6 |
| j.e Counts in samples with spike reads can differ slightly | 2 Low | BiotaWiz and command line: samples with spike reads, any table | **No**: run 0.7.6 |
| j.f Normalized tables: written when a sample has no column; first sample checked | 1 Very low (4 High in the rare second case) | BiotaWiz and command line: runs in which a sample had no column in a table; runs whose first table sample was missing from the spike stats file | **Not needed**: 0.7.6 wrote no tables, or wrong ones |
| j.g Corrected sample weights are used | 3 Medium | BiotaWiz and command line: reruns after correcting a weight, if you used `Spike_Normalized-*` | **Yes**: put the old weight back |
| j.h Taxonomy columns with SILVA 144 or newer | **5 Very high** | BiotaWiz and command line, runs with SILVA ≥ 144: taxonomy everywhere, TIC files | **No**, and not wanted: the old ranks were in the wrong columns |
| j.i Which FASTQ files are paired | 4 High | BiotaWiz and command line, non-standard file names: all files of those samples | **No** for wrong pairings; rename the files |
| j.j Reruns no longer pick up earlier results as samples | 4 High | BiotaWiz and command line, reruns in the same folder: all tables | **No**, and not needed: drop the extra samples |

#### j.a Fewer zOTUs, TAC OTUs and SOTUs: zOTU filter at 0.25 % by default
The shipped arguments file `TACTICPipeline_args.yml` now sets **Minimum relative abundance** (`ClusterZOTUsArgs` ›
`abund_limit`) to `0.0025` instead of `0.000`: a zOTU is kept only if it has more than 0.25 % of the reads in at
least one sample, as recommended by [Reitmeier et al. 2021](https://doi.org/10.1038/s43705-021-00033-z) before alpha
and beta diversity analyses. 0.7.6 used the `0.000` of the shipped file and kept practically all zOTUs (relative
abundances are rounded to 5 decimals, so a zOTU below 0.0005 % in every sample was still removed). This gives fewer
zOTUs, and therefore fewer TAC OTUs and SOTUs. BiotaWiz runs are not affected by this change: BiotaWiz writes its own
value (see j.b).

TACTIC copies the shipped file into the FASTQ folder only when the folder has no `TACTICPipeline_args.yml` yet. A
folder already analysed with 0.7.6 keeps its copy with `abund_limit: 0.000`, so a new run there still keeps all
zOTUs.

**Effect: 5 Very high.** Rare zOTUs are removed, so richness and alpha diversity drop and beta diversity based on
presence/absence (Jaccard, unweighted UniFrac) changes. Measured on 24 mouse samples with **Minimum number of reads
per zOTU** = 2: **186 instead of 493 zOTUs** and **141 instead of 348 SOTUs**. Abundance-weighted measures
(Bray-Curtis, weighted UniFrac) and the dominant taxa change much less.

**Affects:** you, if all of these apply:
- **Interface:** command line with the shipped `TACTICPipeline_args.yml` of 1.0.0, i.e. a FASTQ folder without an
  arguments file, or a template written with `-tf` by 1.0.0. Arguments files that set their own `abund_limit`,
  including the copy left in a folder by 0.7.6, keep their value.
- **Your run:** ZOTU, TAC or TIC mode. The OTU mode (de novo OTUs) is not affected.
- **Files you used downstream:** the zOTU-based files: `ZOTUs-*` (ZOTU), `OTUs-*-TAC*` and `Map-ZOTUs-OTUs-TAC.tab`
  (TAC), `SOTUs-*`, `Map-*.tab` and `ZOTU-Table-*-FullTaxonomy.tab` (TIC), and `PICRUSt2-Inputs/`, in particular for
  richness, alpha diversity and presence/absence-based beta diversity.

**Reproduce the 0.7.6 result with 1.0.0: Yes.** Set `abund_limit: 0.0` under `ClusterZOTUsArgs`.

#### j.b More zOTUs, TAC OTUs and SOTUs: **Minimum relative abundance** `0` is now used
BiotaWiz 1.0.0 to 1.0.2 set **Minimum relative abundance** (`ClusterZOTUsArgs` › `abund_limit`) to `0`, meaning
"keep all zOTUs", and write it as the whole number `0`. Up to 0.7.6, a whole number given for a decimal parameter was
silently ignored and the default `0.0025` was used: every zOTU that did not have more than 0.25 % of the reads in at
least one sample was removed. From 1.0.0 the value is used, so these low-abundance zOTUs are kept. BiotaWiz 1.0.3 and newer
write `0.0025`, so they give the same zOTUs as with 0.7.6.

**Effect: 5 Very high.** The reverse of j.a: many rare zOTUs are kept, so richness and alpha diversity rise strongly
and beta diversity based on presence/absence changes. Measured on the same 24 mouse samples: **493 instead of 186
zOTUs** (2.7 times as many) and **348 instead of 141 SOTUs** (2.5 times).

**Affects:** you, if all of these apply:
- **Interface:** **BiotaWiz 1.0.0 to 1.0.2**. Command line only if the arguments file has whole numbers for decimal
  parameters (e.g. `abund_limit: 0`).
- **Your run:** ZOTU, TAC or TIC mode. The OTU mode (de novo OTUs) is not affected.
- **Files you used downstream:** the zOTU-based files, as in j.a.

**Reproduce the 0.7.6 result with 1.0.0: Yes.** Update BiotaWiz to 1.0.3 or newer, or set **Minimum relative
abundance** to `0.0025`:
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

#### j.c Primer trimming defaults of 17 and 21 bases
If an arguments file does not set **Bases to trim from the start of merged reads** (`stripleft`) and **Bases to trim
from the end of merged reads** (`stripright`) under `TrimBothSidesArgs` / `TrimOneSideArgs`, the defaults are now 17
and 21 bases (341F/785R primers) instead of 5 and 5. Sequences are shorter and no longer start with primer remnants,
which can give **fewer zOTUs**.

**Effect: 3 Medium.** zOTUs that differed only in leftover primer bases, which come from the primers and not from
the bacteria, are merged. zOTU richness can drop; the taxonomic composition stays the same.

**Affects:** you, if all of these apply:
- **Interface:** command line only. Not BiotaWiz and not the shipped arguments file: both set 17 and 21.
- **Your run:** your own arguments file without `stripleft` / `stripright` under `TrimBothSidesArgs` /
  `TrimOneSideArgs`.
- **Files you used downstream:** all files, in particular zOTU counts and richness.

**Reproduce the 0.7.6 result with 1.0.0: Yes.** Set `stripleft: 5` and `stripright: 5` under `TrimBothSidesArgs`, and
`stripleft: 5` under `TrimOneSideArgs` (command line only).

#### j.d `Spike_Normalized-*` values shift between samples
Spike reads are now counted after merging and quality filtering, like the reads in the tables. Each sample is scaled
by the mean spike reads of all samples divided by its own spike reads, so the overall level stays about the same, but
samples shift against each other: a sample that lost more reads in merging and filtering than the average sample gets
**higher** values than in 0.7.6; one that lost fewer gets **lower** values. Samples processed by an earlier version
are processed again automatically.

**Effect: 3 Medium.** Only absolute abundances (microbial load per gram) change; relative composition, alpha and
beta diversity do not. Load differences between samples or groups can become smaller or larger, most where samples
lost different shares of reads in merging and filtering (e.g. sequencing runs or batches of different read quality).

**Affects:** you, if all of these apply:
- **Interface:** BiotaWiz or command line.
- **Your run:** samples with a spike-in standard (BiotaWiz: **Spike-in** step; command line: `spike_amount` above 0 in
  the mapping file).
- **Files you used downstream:** `Spike_Normalized-*.tab`, e.g. for absolute abundances or microbial load. Relative
  abundances and diversity are not affected.

**Reproduce the 0.7.6 result with 1.0.0: No.** There is no parameter for the old counting. Run 0.7.6 (see [How to
reproduce results of an earlier version](https://tum-core-facility-microbiome.github.io/BiotaWiz/docs/getting-started/compare-versions#reproduce)).
The new values are the corrected ones.

#### j.e Counts in samples with spike reads can differ slightly
Spike reads are now removed from the quality-filtered reads (sequence similarity of at least 90 %) instead of from the
raw reads (bowtie2). As before, this is done in every sample, whatever its spike amount. The number of zOTUs and
the read counts can change a little in samples that contain spike reads: the spiked samples, and samples that picked
up spike reads from spiked samples sequenced with them, also in runs without spike-ins.

**Effect: 2 Low.** Only reads very similar to the synthetic spike-in sequences are handled differently; real taxa
are hardly affected. At most a zOTU made of spike reads appears or disappears.

**Affects:** you, if all of these apply:
- **Interface:** BiotaWiz or command line.
- **Your run:** samples that contain spike reads: the spiked samples, or samples sequenced together with them.
- **Files you used downstream:** any count table or sequence file of those samples (`ZOTUs-Table.tab`,
  `OTUs-Table*.tab`, `SOTUs-Table-TIC.tab`, `*-Seqs*.fasta`).

**Reproduce the 0.7.6 result with 1.0.0: No.** Run 0.7.6.

#### j.f Normalized tables: written when a sample has no column, and the first sample is checked
Two cases:

1. **A sample has no column in a table.** A table only has columns for samples with reads that match a zOTU; a
   sample without such reads (e.g. a blank or negative control, possibly only spike reads) has none. 0.7.6 then
   wrote neither `SampleMinCount_Normalized-*` nor `Spike_Normalized-*` for that table; the log only showed
   "Error while normalizing … on file …". Now both are written for the samples in the table (`Spike_Normalized-*` as
   long as every sample has valid spike data).
2. **The first sample of a table is missing from the spike stats file.** 0.7.6 checked that every sample of a table is
   in the spike stats file, but skipped the first one. If it was missing, its raw counts were written unchanged into
   `SampleMinCount_Normalized-*` and `Spike_Normalized-*`, next to normalized samples. Now the run stops with an error,
   as it already did for any other sample. No normal run is known to lead to this case.

**Effect: 1 Very low** for case 1: tables that were missing are now written; no value that 0.7.6 wrote changes.
**4 High** for case 2, for that one sample: its values in the normalized tables were raw counts, not comparable with
the other samples. The plain count tables were correct.

**Affects:** you, if all of these apply:
- **Interface:** BiotaWiz or command line.
- **Your run:** case 1: a sample with no column in a table, with or without spike-ins (`Analysis_log.txt` of 0.7.6
  contains "Error while normalizing"); case 2: a table whose first sample was missing from the spike stats file.
- **Files you used downstream:** case 1: none (the tables were missing); case 2: `SampleMinCount_Normalized-*.tab` and
  `Spike_Normalized-*.tab`.

**Reproduce the 0.7.6 result with 1.0.0: Not needed.** In case 1, 0.7.6 wrote no tables; in case 2, its values for
the first sample should not be used.

#### j.g Corrected sample weights are used
If you corrected a sample's weight (`total_weight_in_g` in the mapping file) and the sample was reused from an earlier
run, 0.7.6 still used the old weight. Its spike-normalized values now change in proportion (half the weight gives
double the values).

**Effect: 3 Medium** for the corrected samples: their absolute abundance changes by the same factor as the weight,
which can change comparisons of microbial load between groups. Relative composition and diversity are not affected.

**Affects:** you, if all of these apply:
- **Interface:** BiotaWiz or command line.
- **Your run:** you ran the same FASTQ folder again with unchanged settings after correcting a sample's weight
  (`total_weight_in_g` in the mapping file).
- **Files you used downstream:** `Spike_Normalized-*.tab` (the values of the corrected samples).

**Reproduce the 0.7.6 result with 1.0.0: Yes.** Put the old weight back in the mapping file.

#### j.h Taxonomy columns with SILVA 144 or newer
SILVA 144 added a kingdom rank (e.g. *Bacillati*) between domain and phylum. 0.7.6 did not know it: the kingdom ended
up in the phylum column and every lower rank one column to the right (the genus in the species column). 1.0.0 puts
each rank in its own column and also writes `*-8Levels.*` files with the kingdom. Because TIC clusters by taxonomy,
the SOTUs (`SOTUs-Table-TIC.tab`, `SOTUs-Seqs-TIC.fasta`) differ as well. With *SILVA 138 or older*, the taxonomy is
the same as in 0.7.6.

**Effect: 5 Very high.** In 0.7.6 every rank below the domain was in the wrong column (the kingdom as phylum, the
genus as species), so the taxonomic composition at every rank was wrong. TIC also formed different SOTUs, which
changes richness and diversity in TIC mode.

**Affects:** you, if all of these apply:
- **Interface:** BiotaWiz or command line.
- **Your run:** the run used **SILVA 144** or newer. You can see it in the 0.7.6 tables: a kingdom name ending in
  *-ati* (e.g. *Bacillati*) stands where the phylum should be.
- **Files you used downstream:** the taxonomy of any table, sequence file or Krona chart (summaries per phylum,
  family, genus …), and in TIC mode the SOTU files (`SOTUs-Table-TIC.tab`, `SOTUs-Seqs-TIC.fasta`, `Map-*.tab`).

**Reproduce the 0.7.6 result with 1.0.0: No, and it should not be used:** its ranks were in the wrong columns.

#### j.i Which FASTQ files are paired
Read direction is now taken from the last marker in the file name. New markers are recognized (`read1`/`read2`,
`forward`/`reverse`, `R1`/`R2` after `-` or `.`), and `_F_`/`_R_` in the middle of a name no longer count as
markers. Without a mapping file, FASTQ files without any read marker are now processed as single-end samples; 0.7.6
ignored them (e.g. index reads such as `_I1_001.fastq.gz` now show up as samples). Samples that earlier versions
paired wrongly, missed or treated as single-end can now be paired (or the other way round). Check `Pipeline_log.txt`
and `Missed_samples.txt` to see that every sample was found and paired as expected.

**Effect: 4 High** for the samples it applies to: a missing sample changes the group comparisons, and a sample
analysed from single reads instead of merged pairs gets shorter, less accurate sequences and other zOTUs, so it can
look different from the rest of its group. Samples with standard file names do not change.

**Affects:** you, if all of these apply:
- **Interface:** BiotaWiz or command line.
- **Your run:** FASTQ file names without a standard read marker (e.g. `_R1_001`/`_R2_001`, or `_1`/`_2` directly
  before the extension), or, without a mapping file, other FASTQ files in the folder.
- **Files you used downstream:** all files, for the samples concerned (see `Pipeline_log.txt` and
  `Missed_samples.txt`).

**Reproduce the 0.7.6 result with 1.0.0: No** for wrong pairings. If a sample is now missing or single-end, rename its
files to a supported pattern (listed in the README).

#### j.j Reruns no longer pick up earlier results as samples
FASTQ files inside the `<sample>__processed/` folders of an earlier run could be picked up again as extra samples with
numbered names if a run was stopped or failed at an early stage. They are now skipped.

**Effect: 4 High** for the runs it applies to: copies of samples count as extra replicates, which shifts group
averages and makes group differences look more significant than they are. The composition of each real sample
does not change.

**Affects:** you, if all of these apply:
- **Interface:** BiotaWiz or command line.
- **Your run:** you ran TACTIC again in a FASTQ folder that already had `<sample>__processed/` folders. You can see it
  in the old tables: extra sample columns with numbered names.
- **Files you used downstream:** all tables, in particular for group comparisons and statistics.

**Reproduce the 0.7.6 result with 1.0.0: No, and not needed:** the extra samples were copies of real samples. Remove
them from the old tables before comparing, or start from a fresh folder with your raw FASTQs.

#### j.k Runs with the first 1.0.0 build (2026-10-23 to 2026-11-05)
The first build already had j.b, j.c and j.h to j.j. It did not have j.a and j.d to j.g. Compared with the first
build, this build changes:

| Change | Effect | Affects | Reproduce the first-build result with this build? |
|---|---|---|---|
| zOTU filter at 0.25 % by default (j.a) | **5 Very high** | Command line with the shipped arguments file in a new FASTQ folder: ZOTU, TAC, TIC files | **Yes**: set `abund_limit: 0.0` |
| `Spike_Normalized-*` values shift between samples (j.d) | 3 Medium | Runs with spike-ins, if you used `Spike_Normalized-*` | **No** |
| Counts in samples with spike reads can differ slightly (j.e) | 2 Low | Samples with spike reads, any table | **No** |
| Normalized tables written when a sample has no column; first sample checked (j.f) | 1 Very low (4 High in the rare second case) | Runs in which a sample had no column in a table; runs whose first table sample was missing from the spike stats file | **Not needed** |
| Corrected sample weights are used (j.g) | 3 Medium | Reruns after correcting a weight | **Yes**: put the old weight back |

BiotaWiz runs with the first build are affected by the zOTU filter only through the BiotaWiz version (j.b): BiotaWiz
1.0.0 to 1.0.2 kept all zOTUs with either build, BiotaWiz 1.0.3 filters at 0.25 % with either build.

The first build cannot be run as `1.0.0` any more. Find the date of a run in its `Analysis_<timestamp>/` folder name.

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
| Tree leaf labels without extra quotes | 1 Very low | BiotaWiz and command line, if you used the trees (`*-Tree-nj*.tre`) | **Not needed**: same results wherever the labels matched |

#### d.a Tree-based results (UniFrac, phylogenetic diversity)
The leaf labels in the `.tre` files had extra quotes, so tools could fail to match tree leaves to the table rows.
Results computed with the trees can differ now that the labels match. Tree shapes, branch lengths, counts and taxonomy
are not affected.

**Effect: 1 Very low.** Tree shapes and branch lengths did not change, only the labels. UniFrac and phylogenetic
diversity are the same wherever the tool matched the labels; where it did not, the tool could not use the tree at
all, rather than giving other numbers.

**Affects:** you, if all of these apply:
- **Interface:** BiotaWiz or command line.
- **Your run:** any.
- **Files you used downstream:** the trees (`*-Tree-nj*.tre`), for UniFrac or phylogenetic diversity (e.g. in phyloseq or QIIME 2).

**Reproduce the 0.7.5 result with 0.7.6: Not needed.** Where the old labels matched, the results are the same; where
they did not, the new results are the correct ones.

## [0.7.5] - 2026-07-15

### a. Changed
- **Choosing the newest SILVA release**: decided by the version number; the timestamp of the `.arb.gz` file only breaks ties.
- **Argument and log files kept with the results**: the arguments file used is copied into the `Analysis_<timestamp>/` folder, and the `*.yml` and `*log.txt` files in that folder are kept.
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
| Newest SILVA release chosen by version number | 3 Medium | BiotaWiz and command line, database folders with several SILVA releases: taxonomy, TIC files, trees | **Yes**: keep only the release you want |

#### b.a Taxonomy, if the database folder holds several SILVA releases
The newest release is now chosen by its version number instead of the file date, so another release than before may
be used. Taxon names can then differ, and with them the SOTUs of TIC, which clusters by taxonomy. The counts of
zOTUs, OTUs and TAC OTUs are not affected.

**Effect: 3 Medium.** Taxa that were renamed or reclassified between the two SILVA releases get other names, so the
composition can shift for those taxa, and TIC can form somewhat different SOTUs. Most taxa keep their names;
diversity based on zOTUs, OTUs and TAC OTUs does not change.

**Affects:** you, if all of these apply:
- **Interface:** BiotaWiz or command line.
- **Your run:** the database folder held more than one SILVA release.
- **Files you used downstream:** the taxonomy of any table, sequence file or Krona chart; in TIC mode the SOTU files; the trees.

**Reproduce the 0.7.4 result with newer versions: Yes.** Keep only the release you want in the database folder (see
[How to reproduce results of an earlier version](https://tum-core-facility-microbiome.github.io/BiotaWiz/docs/getting-started/compare-versions#reproduce), step 3).

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
| Incomplete SILVA index is rebuilt | **5 Very high** | BiotaWiz and command line, interrupted database setups: taxonomy, TIC files, trees | **No**, and not wanted: the old taxonomy was incomplete |

#### d.a Taxonomy after an interrupted database setup
An incomplete SILVA index (e.g. after a crash during indexing) is now rebuilt instead of used. Runs that used an
incomplete index could have missing or wrong taxonomy, and therefore different SOTUs (TIC) and trees.

**Effect: 5 Very high** for the runs it applies to: the taxonomy could be missing or wrong for many sequences, so
the taxonomic composition was wrong, and so were the TIC SOTUs and with them richness and diversity in TIC mode.

**Affects:** you, if all of these apply:
- **Interface:** BiotaWiz or command line.
- **Your run:** the SILVA indexing was interrupted (e.g. BiotaWiz, Docker or the computer was stopped during the first
  database setup).
- **Files you used downstream:** the taxonomy of any table, sequence file or Krona chart; in TIC mode the SOTU files; the trees.

**Reproduce the 0.7.3 result with newer versions: No, and it should not be used:** its taxonomy was incomplete.
