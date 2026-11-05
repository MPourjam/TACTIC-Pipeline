# What's new in TACTIC: the short version

This page sums up each TACTIC version in plain words. Most updates add features or make results more accurate,
and in most cases you do not need to do anything. Results you already have stay valid for the version that made
them; the version is in the first line of `Analysis_log.txt`.

The full details, with every file and setting, are in the [changelog](CHANGELOG.md). How to compare or reproduce
results of different versions is explained in
[Comparing results of different versions](https://tum-core-facility-microbiome.github.io/BiotaWiz/docs/getting-started/compare-versions).

## Do I need to do anything?

| If you … | then … |
|---|---|
| use BiotaWiz | update BiotaWiz to 1.0.3 or newer. Nothing else is needed. |
| ran TACTIC 1.0.0 with BiotaWiz 1.0.0 to 1.0.2 | your zOTU, TAC and TIC results kept even very rare zOTUs. If you want the usual filter (0.25 %), update BiotaWiz and run the analysis again. This is quick: your processed samples are reused. |
| used spike-in normalized tables (`Spike_Normalized-*`) | run the analysis again with TACTIC 1.0.1 for more accurate absolute abundances. |
| run TACTIC on the command line with your own arguments file | compare it with the new `TACTICPipeline_args.yml`, especially `abund_limit` and the primer trimming values. |
| none of the above | nothing. |

## 1.0.1 (2026-11-05)

- **More accurate spike-in normalization.** Spike-in reads are now counted after the same quality steps as all
  other reads. The `Spike_Normalized-*` values of single samples can shift a little against each other. Samples are
  processed again automatically the first time you run 1.0.1.
- **Corrected sample weights are always used.** If you fix a weight in the mapping file, the next run uses it.
- **zOTU filter of 0.25 % by default.** A zOTU is kept if it makes up at least 0.25 % of the reads in at least one
  sample, as recommended by [Reitmeier et al. 2021](https://doi.org/10.1038/s43705-021-00033-z). This removes most
  sequencing noise. BiotaWiz users get this filter with BiotaWiz 1.0.3; to keep all zOTUs, set
  **Minimum relative abundance** to 0.
- **Clearer results folders.** Each folder has an `Output files explanation.html` page that opens in your browser and
  marks every file as *for your analysis* or *complementary*. Temporary files are cleaned up.
- **Unchanged:** runs without spike-ins give the same numbers as 1.0.0 with the same settings.

## 1.0.0 (2026-10-23)

- **Support for the newest SILVA database (144 and newer).** It adds a kingdom level (e.g. *Bacillati*) to the
  taxonomy. TACTIC now writes the usual 7-level files and, in addition, `*-8Levels.*` files with the kingdom. If you
  used TACTIC 0.7.6 with SILVA 144, please run the analysis again: 0.7.6 did not know the new level, so the taxonomy
  columns were shifted by one.
- **Ready-made PICRUSt2 input** in the `PICRUSt2-Inputs/` folder.
- **More reliable FASTQ detection.** More file name patterns are recognized, and earlier results in the same folder
  are no longer picked up as new samples.
- **Settings are applied as written.** A value such as `0` in the arguments file is now used; earlier versions
  sometimes ignored it. For BiotaWiz 1.0.0 to 1.0.2 this meant that all zOTUs were kept (see the table above).
- **The version is written at the top of the logs**, and the logs and the arguments file are kept with the results.

## 0.7.6 (2026-09-22)

- **Several analysis modes in one run**, e.g. OTU, TAC and TIC together.
- **Faster first setup:** the SILVA database comes ready-indexed from the TACTIC file server.
- **Trees work in more tools:** the labels in the `.tre` files no longer have extra quotes, so tools like phyloseq or
  QIIME 2 match them to the tables.

## 0.7.5 (2026-07-15)

- **The newest SILVA release is chosen by its version number** when the database folder holds several.
- **The arguments file and the logs are kept** in each `Analysis_…` folder.

## 0.7.4 (2025-12-19)

- **An interrupted database setup is repaired automatically** the next time you run TACTIC.
- **Quickstart** in the README, and `-tf` also writes a template of the arguments file.
