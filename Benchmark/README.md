# TACTICBench Reproduction Bundle

This directory packages the inputs, analysis configuration, software, metadata, and archived outputs for the run **`TACTICBENCH_RESULTS-TIC-Archaea-mock4-mock-5_SRP`**. It is intended to let readers inspect the exact analysis and rerun TACTICBench from the ASV/zOTU/OTU/TIC/TAC feature tables used in the run.

## At a Glance

- **Scope:** 16 mock communities and 27 samples.
- **Analysis methods:** DADA2 ASVs, UNOISE3 zOTUs, UPARSE OTUs, TIC, and TAC.
- **Distance calculations:** Jaccard, weighted Jaccard, Bray-Curtis, UniFrac/gUniFrac, and Aitchison as specified in the effective YAML configuration.
- **Phylogenetic mode:** augmented mismatch; the exact mapping thresholds and other settings are recorded in `config/original_config_effective.yaml`.
- **Raw-read storage:** FASTQ files are not duplicated in this bundle. Their download locations and source papers are linked below and in `metadata/FASTQ_and_publication_links.tsv`.

## Directory Contents

| Path | Contents |
|---|---|
| `code/tacticbench/` | TACTICBench Python module, plotting scripts, and package metadata used for this analysis. |
| `config/run_config.yaml` | Relocatable config for rerunning the analysis from this directory. |
| `config/original_config_input.yaml` | Exact input config snapshot from the archived run. |
| `config/original_config_effective.yaml` | Resolved config snapshot used by the archived run; use this as the authoritative record of effective settings. |
| `environment/tacticbench-reproduction.yaml` | Minimal Conda environment for TACTICBench and the external tree/mapping tools. |
| `data/communities/` | Required expected-reference and method-level tables/sequences, sample metadata, and TIC/TAC cluster-membership maps. Raw/intermediate FASTQs are intentionally excluded. |
| `metadata/` | Community overview tables, original run manifest, and links to raw reads/publications. |
| `results/archive_run/` | Archived outputs from the selected run, including figures, QC, distance tables, manifest, and run logs. |
| `results/reproduced_run/` | Destination for a fresh run using `config/run_config.yaml`; created when the commands below are run. |

## Reproduce the Analysis

Run commands from this directory so the relative paths in `config/run_config.yaml` resolve correctly.

```bash
conda env create -f environment/tacticbench-reproduction.yaml
conda activate tacticbench-reproduction
python -m pip install --no-deps -e code/tacticbench
python -m tacticbench.cli run-all --config config/run_config.yaml
```

The run config calls `mafft`, `FastTree`, `blastn`, and `makeblastdb` by executable name. Ensure they are available on `PATH` in the environment. The command regenerates TACTICBench's mapping, distance calculations, QC, and figures from the included feature tables and reference data; it does **not** rerun DADA2, UNOISE3, UPARSE, TIC, or TAC from raw reads.

To rerun only one figure after the analysis tables have been generated, for example:

```bash
python -m tacticbench.cli plot-fig21 --config config/run_config.yaml
```

The archived run outputs are preserved under `results/archive_run/`; reruns write to the separate `results/reproduced_run/` directory.

## Raw FASTQ Data and Source Publications

The mockrobiota FASTQ URLs below point directly to externally hosted compressed reads. For ZIEL and ZYMO, the links point to the corresponding NCBI SRA run records, from which FASTQ files can be downloaded. These links are data provenance, not files included in this package. Cite the [mockrobiota resource](https://doi.org/10.1128/mSystems.00062-16) and the original study where relevant.

| TACTICBench community | FASTQ source | Source publication |
|---|---|---|
| ZIEL I (`mock-ZIEL_1`) | [SRR13005984](https://www.ncbi.nlm.nih.gov/sra/?term=SRR13005984), [SRR13005985](https://www.ncbi.nlm.nih.gov/sra/?term=SRR13005985) | [Abellan-Schneyder et al. 2021](https://doi.org/10.1128/mSphere.01202-20) |
| ZIEL II (`mock-ZIEL_2`) | [SRR13005964](https://www.ncbi.nlm.nih.gov/sra/?term=SRR13005964), [SRR13005965](https://www.ncbi.nlm.nih.gov/sra/?term=SRR13005965), [SRR13005966](https://www.ncbi.nlm.nih.gov/sra/?term=SRR13005966) | [Abellan-Schneyder et al. 2021](https://doi.org/10.1128/mSphere.01202-20) |
| ZYMO (`mock-ZYMO`) | [SRR13006004](https://www.ncbi.nlm.nih.gov/sra/?term=SRR13006004), [SRR13006002](https://www.ncbi.nlm.nih.gov/sra/?term=SRR13006002), [SRR13006003](https://www.ncbi.nlm.nih.gov/sra/?term=SRR13006003) | [Abellan-Schneyder et al. 2021](https://doi.org/10.1128/mSphere.01202-20) |
| Mock 1 (`mock-4`) | [Forward](https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-4/mock-forward-read.fastq.gz), [reverse](https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-4/mock-reverse-read.fastq.gz), [index](https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-4/mock-index-read.fastq.gz) | [Bokulich et al. 2013](https://doi.org/10.1038/nmeth.2276) |
| Mock 2 (`mock-5`) | [Forward](https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-5/mock-forward-read.fastq.gz), [reverse](https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-5/mock-reverse-read.fastq.gz), [index](https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-5/mock-index-read.fastq.gz) | [Bokulich et al. 2015, preprint](https://doi.org/10.7287/peerj.preprints.934v2) |
| Mock 3 (`mock-12`) | [Forward](https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-12/mock-forward-read.fastq.gz), [reverse](https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-12/mock-reverse-read.fastq.gz) | [Callahan et al. 2016](https://doi.org/10.1038/nmeth.3869) |
| Mock 4 (`mock-13`) | [Forward](https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-13/mock-forward-read.fastq.gz), [reverse](https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-13/mock-reverse-read.fastq.gz) | [Kozich et al. 2013](https://doi.org/10.1128/AEM.01043-13) |
| Mock 5 (`mock-14`) | [Forward](https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-14/mock-forward-read.fastq.gz), [reverse](https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-14/mock-reverse-read.fastq.gz) | [Kozich et al. 2013](https://doi.org/10.1128/AEM.01043-13) |
| Mock 6 (`mock-15`) | [Forward](https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-15/mock-forward-read.fastq.gz), [reverse](https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-15/mock-reverse-read.fastq.gz) | [Kozich et al. 2013](https://doi.org/10.1128/AEM.01043-13) |
| Mock 7 (`mock-16`) | [Forward](https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-16/mock-forward-read.fastq.gz), [reverse](https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-16/mock-reverse-read.fastq.gz) | [Schirmer et al. 2015](https://doi.org/10.1093/nar/gku1341) |
| Mock 8 (`mock-18`) | [Forward](https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-18/mock-forward-read.fastq.gz), [reverse](https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-18/mock-reverse-read.fastq.gz) | [Tourlousse et al. 2017](https://doi.org/10.1093/nar/gkw984) |
| Mock 9 (`mock-19`) | [Forward](https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-19/mock-forward-read.fastq.gz), [reverse](https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-19/mock-reverse-read.fastq.gz) | [Tourlousse et al. 2017](https://doi.org/10.1093/nar/gkw984) |
| Mock 10 (`mock-20`) | [Forward](https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-20/mock-forward-read.fastq.gz), [reverse](https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-20/mock-reverse-read.fastq.gz) | [Gohl et al. 2016](https://doi.org/10.1038/nbt.3601) |
| Mock 11 (`mock-21`) | [Forward](https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-21/mock-forward-read.fastq.gz), [reverse](https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-21/mock-reverse-read.fastq.gz) | [Gohl et al. 2016](https://doi.org/10.1038/nbt.3601) |
| Mock 12 (`mock-22`) | [Forward](https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-22/mock-forward-read.fastq.gz), [reverse](https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-22/mock-reverse-read.fastq.gz) | [Gohl et al. 2016](https://doi.org/10.1038/nbt.3601) |
| Mock 13 (`mock-23`) | [Forward](https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-23/mock-forward-read.fastq.gz), [reverse](https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-23/mock-reverse-read.fastq.gz) | [Gohl et al. 2016](https://doi.org/10.1038/nbt.3601) |

The open-access article landing pages above provide publication information; this bundle links to the publications rather than redistributing PDF copies. Please check the publication's license before separately redistributing any downloaded PDF.

## Provenance and Interpretation

- `metadata/mock_communities_paper_table.tsv` summarizes the mock names, origins, primers, sample counts, expected species counts, and composition types.
- `metadata/mock_communities_summary.tsv` and `metadata/mock_metadata_for_paper.tsv` contain additional expected-community and sample-level overview information.
- `metadata/original_run_manifest.tsv` records the inputs selected for the archived run. Its source paths reflect the original analysis workspace; the bundled copies are under `data/communities/`.
- `results/archive_run/config_input.yaml`, `config_effective.yaml`, `manifest.tsv`, and `tacticbench.log` are preserved as the original run record.
- `data/communities/` contains precomputed outputs for the five methods. Recreating those method outputs from FASTQ requires their upstream pipelines and parameters, which are outside this TACTICBench reproduction bundle.

## Citation

Please cite the original publications for the mock datasets used and the TACTICBench/TAC/TIC work associated with this analysis. The selected dataset-specific citations are linked in the table above; the full mock-community overview is in `metadata/mock_communities_paper_table.tsv`.
