# Reproduction Checks

Checks were performed on 17 October 2026 in a newly created minimal Linux Conda
environment containing Python 3.11.14, NumPy 2.1.3, Matplotlib 3.9.2, and SeqKit
2.12.0, together with the bundled USEARCH 11.0.667 executable. The resolved Conda
packages are recorded in `environment/tested-linux-64.explicit.txt` at the bundle
root. `validation_summary.json` records the final outcomes.

| Check | Scope | Evidence |
|---|---|---|
| Input archive | ZIP SHA-256, manifest/config consistency, and CRC/SHA-256 of all four payloads | `input_verification.json` |
| Isolated execution | Tools, plots, and smoke runs from a relocated copy with cleared PYTHONPATH and restricted executable PATH | `relocated/test_report.json` |
| Plot regeneration | All three feature-versus-length plots, GG2 top-25 and vertical top-5 plots from saved numerical summaries | `relocated/plots/` |
| End-to-end smoke | 64 eligible real sequences per database, all 20 alignment windows, full/filtered reference rows, and all identities | `relocated/smoke/validation.json` |
| Complete reference populations | All input records; full and filtered record/species counts for SILVA, LTP, GG2; GG2 top-25 species frequencies | `populations/reference_population_validation.json` |
| Complete LTP rerun | All 22 rows and every numerical column compared to the archived table | `full_ltp/validation.json` |
| GG2 species results | Fresh top-25 and top-5 unique sequences and clusters, including pooled groups, compared to archived measurements | `gg2_species/species_summary_validation.json` |

The smoke tests use the original cutting, filtering, deduplication, and clustering
workflow. A smoke pass checks all configured operations and output rows; it does
not compare a 64-record subset to full-database cluster counts.

The complete LTP test uses all 19601 eligible reference records and retains 19452
records after the V1-V9 length filter. The population audit likewise examines all
112334 eligible SILVA and 128098 eligible GG2 reference records. The GG2 species
check recomputes all measurements for the 25 selected species, their pooled
25446 unique sequences, and the pooled 9951 sequences of the top five species.

The full 20-window clustering sweeps for SILVA and GG2 were not rerun during this
validation. Their archived tables and exact FASTA inputs are included, their
complete reference populations were checked, and their entire workflows were
exercised on smoke subsets. Use the full reproduction command to recompute those
complete sweeps.

Validation exposed thread-dependent clustering, historical taxonomy-definition
differences, and a concurrent FASTA-index creation race. The recovered settings
and portable execution fix are described in `metadata/reproduction_findings.md`
at the bundle root. `diagnostics/` retains the initial discrepancy reports; the
final checks use the explicitly recorded historical settings.

`software.json` files record executable paths, versions, worker/thread settings,
and the USEARCH checksum. Absolute paths in test reports identify the local
validation workspace; the portable commands resolve paths from the bundle and
do not require those original locations.
