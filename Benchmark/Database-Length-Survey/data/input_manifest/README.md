# Length-Survey Input Manifest

`inputs.tsv` records each surviving source path, kind, size, SHA-256 and portable
data path. FASTQs live in the raw ZIP under `payloads/`; small zOTU FASTAs live
outside it under `small_inputs/`. Multiple paths may refer to the same deduplicated
payload. `metadata/input_sources.tsv` provides accession/source links when these
can be inferred from recorded paths or the existing benchmark overview.

`observations.tsv` freezes the per-script selections, original path/identifier,
mock/group label, read-sampling limit, report/line evidence, recovery caveat and
expected saved statistics. An observation is **one script measuring one file**,
not one independent community or sample. R1/R2 are separate observations and
repeated analysis outputs remain separate when the original survey included them.

`historical_validation.tsv` compares recovered data with the available report
statistics. `length_distributions.json` contains exact full length distributions
and first-1,000/10,000-read summary statistics from the recovered inputs.

`findings.tsv` lists historical measurements with unavailable source files. A
missing input is not synthesized from a histogram or replaced by another analysis.
The zOTU detailed report lacked full paths; recovered rows identify matching
analysis directory/statistics, not guaranteed original tree position. R1 file
lane-index spelling discrepancies are documented in `recovery_note` and checked
against the saved sampled statistics.

`SHA256SUMS` lists every ZIP member except itself, including the embedded manifests.
The archive verifier rejects duplicate member names, unsafe paths, inventory
differences, checksum mismatches and CRC errors. SHA-256 values are current
recovery-time checksums; no historical checksums were saved by the old survey.
