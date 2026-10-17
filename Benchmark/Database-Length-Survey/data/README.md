# Reference Input Archive

Zenodo record DOI: **TODO: add after deposition.**

Zenodo archive download URL: **TODO: add after deposition.**

Upload `stepwise_diversity_inputs.zip` and its `.sha256` checksum file. Link the
deposit to the small code, configuration, environment, original summaries, and
documentation in the parent directory. Both parts are needed to recompute the
analyses from reference sequences.

The companion upload `Step-wise_Diversity_code_and_results.zip` is generated
beside this project directory and contains the small collection and recorded
verification results. It excludes this large reference ZIP and extracted FASTAs.

The ZIP contains the exact prepared SILVA 138.2, LTP 10_2024, and Greengenes2
2024.09 aligned reference FASTAs used by the step-wise 16S diversity analyses,
plus the GG2 taxonomy-header list and an embedded manifest. These inputs are
reference sequences, rather than sample FASTQ reads. The archive omits the large
generated intermediate FASTAs because the supplied workflow regenerates them.

`input_archive.json` records the final ZIP size and SHA-256 checksum.
`input_manifest/files.tsv` records each payload's original path, portable path,
purpose, size, SHA-256, and database source URL. The manifest also exists inside
the ZIP alongside the run configuration.

From the parent directory:

```bash
python scripts/prepare_inputs.py --extract
python scripts/run_stepwise.py reproduce --threads 19 --workers 3
```

`prepare_inputs.py` verifies the ZIP SHA-256 and each uncompressed payload before
finishing extraction. To verify without extracting, omit `--extract`.

Suggested Zenodo description:

> Reference inputs for the step-wise 16S diversity analysis in the TACTIC study.
> The archive contains the prepared aligned SILVA 138.2, Living Tree Project
> 10_2024, and Greengenes2 2024.09 FASTA datasets, the Greengenes2 taxonomy-header
> list, and checksummed manifests. These files support reproduction of sequence
> length, unique-feature, species-label, and clustering analyses using the
> accompanying code and environment.
