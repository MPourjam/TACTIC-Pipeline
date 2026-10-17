# Recovered Historical Settings

The reproduction checks identified two calculation settings that need to be
preserved explicitly, and one indexing race in parallel execution.

## USEARCH Thread Count

The original wrapper calculates `POOL_SIZE = int(cpu_count() * 0.8)`. On the
24-core analysis host this is 19. Changing USEARCH 11's thread count changes its
parallel greedy clustering, even when the sorted sequence input is identical.

For example, the 3184 unique V1-V9 sequences of the most frequent GG2 species,
Faecalibacterium duncaniae, produced 74 clusters at identity 0.987 using two
threads and the archived 60 clusters using 19 threads. The complete top-25 and
top-5 tables match the archived measurements with 19 threads. The pooled top-25
input sequence bodies have the same SHA-256 in the archived intermediates and
the fresh selection:
`a52eaad995f17da72953ca2dd4ffaecb625a9dbf612a24e0fe4b4198d00b232e`.

`config/runs.json` records 19 threads per dataset. The launcher uses these values
unless `--threads` is explicitly supplied. Lower `--workers` to limit concurrent
jobs without changing the clustering thread count.

## Species-Label Definition

The archived SILVA and LTP runs retained the source headers and counted the first
two words in the final species field, ignoring case. These definitions recover
15407 full-reference SILVA species labels and 19322 LTP species labels. The
filtered V1-V9 counts are 15396 and 19195, respectively.

The later GG2 run uses the strict taxonomy normalization in the source snapshot
and full seven-level taxonomy keys, giving 8082 full-reference labels and 8041
filtered V1-V9 labels. Its top-species figures also use normalized full taxonomy
keys to match sequence sets to species labels.

The current source snapshot normalizes all headers and counts full taxonomy
keys, which changes the older SILVA/LTP populations. The small adapter
`scripts/historical_taxonomy.py` preserves the earlier definition for those two
datasets. These profiles are recorded in the run configuration and validation
reports. The original Python source files are retained without modifications.

## Shared FASTA Index

Multiple window workers can ask SeqKit to index the same filtered alignment at
the same time. A worker can then read the index before another worker finishes
writing it. The isolated test exposed this as an empty 8000-column window even
though the shorter 2000-column window contained sequences.

The portable launcher builds that shared index synchronously immediately after
the filtered alignment is written, before the worker pool starts. This changes
execution ordering rather than the sequence selection or clustering parameters.

The small failed diagnostic reports are retained under
`results/validation_checks/diagnostics/` to document why these settings were
recovered. The final verification reports are separate.
