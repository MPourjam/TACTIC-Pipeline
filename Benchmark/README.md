# TACTIC Reproduction Packages

| Project | Small, runnable bundle | Large data deposit |
|---|---|---|
| Current mock-community benchmark | [tacticbench_TIC_Archaea_mock4_mock5_SRP](tacticbench_TIC_Archaea_mock4_mock5_SRP/README.md) | `data/tacticbench_inputs_with_fastq.zip` inside that bundle |
| Historical read/zOTU length survey | [length_survey_reproduction](length_survey_reproduction/README.md) | `data/length_survey_inputs.zip` inside that bundle |

Each project's `data/README.md` contains placeholders for its future Zenodo DOI
and archive download URL. Keep each large ZIP separate from its small code,
configuration, environment, metadata and reference/feature files when uploading,
and cross-link the deposits so readers obtain both parts.

The local `.gitignore` excludes raw ZIPs, extracted length-survey FASTQ payloads,
temporary archives and fresh rerun outputs. Checksum sidecars and manifests remain
eligible for version control. Upload the raw ZIPs to Zenodo rather than GitHub.

Each bundle has its own minimal Conda environment. TACTICBench needs scientific
Python libraries plus BLAST/MAFFT/FastTree; the length survey requires only Python
and also runs in the TACTICBench environment without additions.

The benchmark input audit passed against the archived run. The historical length
survey has explicitly documented missing earlier source files; see its coverage
notes before describing that package as a complete reproduction of every old
report. No upstream TACTIC/DADA2 pipeline was rerun during packaging.
