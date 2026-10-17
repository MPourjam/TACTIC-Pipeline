# Input Archive and Zenodo Deposit

Zenodo record DOI: **TODO: add the DOI after deposition.**

Zenodo archive download URL: **TODO: add the file download URL after deposition.**

Upload `tacticbench_inputs_with_fastq.zip`, its `.sha256` sidecar, and the small
code/configuration bundle together or link their separate deposits explicitly.
Do not upload only the raw archive: the runnable module and portable config are
provided separately in this reproduction directory.

| Item | Value |
|---|---|
| Archive | `tacticbench_inputs_with_fastq.zip` |
| Bytes | 11,413,439,984 |
| SHA-256 | `8bbfc6606fbe774d48688b1fbd0ebd36504f67408235d7ae3a896f2e5499d565` |
| Scope | 119 distinct FASTQ payloads recovered from 249 original paths, 192 direct benchmark inputs, supporting references and upstream provenance |
| Communities / samples | 16 / 27 |
| Detailed inventory | `input_manifest/README.md`, `files.tsv`, `fastq_files.tsv`, `pipeline_inputs.tsv` |
| Validation | `input_archive_validation.json`; all ZIP members verified by CRC and SHA-256 |

From the bundle root:

```bash
python scripts/verify_input_archive.py data/tacticbench_inputs_with_fastq.zip
unzip -n data/tacticbench_inputs_with_fastq.zip -d data
```

The smaller benchmark inputs are already present under `data/communities/`.
Extraction adds reads and supporting provenance without overwriting these inputs.
The raw archive retains an immutable packaging-time context snapshot; use the
current top-level `environment/tacticbench-reproduction.yaml` when creating the
environment. That environment now explicitly includes the declared Seaborn
dependency required for the enabled Figure 21 density panel.

FASTQs and reference databases may have source-specific redistribution and
citation conditions. Publication PDFs are not redistributed; their DOI/URL links
are in `metadata/FASTQ_and_publication_links.tsv` and the main README.
