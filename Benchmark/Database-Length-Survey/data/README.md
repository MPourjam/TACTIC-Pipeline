# Raw-Input Archive

Zenodo record DOI: **TODO: add the DOI after deposition.**

Zenodo archive download URL: **TODO: add the file download URL after deposition.**

The large file is `length_survey_inputs.zip` (**21,833,283,472 bytes**, approximately
21.8 GB). Its SHA-256 digest is
`d7a1da223d1509649008df655087d12a14d9b2d5d95734186d8315f0e7f3ec1d`.
The size, digest, payload counts and validation result are also recorded in `input_archive_validation.json`;
the adjacent `length_survey_inputs.zip.sha256` is the upload/download checksum.
The ZIP contains 268 unique FASTQ payloads under `payloads/`, plus its input manifest.
Deduplication preserves separate original-path and observation records.

The 62 small required zOTU FASTAs are distributed separately under `small_inputs/`
alongside the code. They are not reference databases. **The ZIP alone is not the
complete reproduction package**: publish or link this directory's code, small
FASTA inputs, environment, manifests and README as well.

From the bundle root:

```bash
python scripts/verify_input_archive.py data/length_survey_inputs.zip
unzip -n data/length_survey_inputs.zip -d data
python scripts/run_length_survey.py
```

Extraction is optional: `run_length_survey.py --archive data/length_survey_inputs.zip`
streams the FASTQs from the ZIP while reading small FASTAs from `data/small_inputs`.

See `input_manifest/findings.tsv` for unavailable historical measurements and
`input_manifest/historical_validation.tsv` for comparisons against saved reports.
Successful archive checksums do not imply every historical survey input survived.
