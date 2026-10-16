# TACTIC Pipeline

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Docker](https://img.shields.io/badge/docker-ghcr.io%2Fmpourjam%2Ftactic--pipeline-blue)](https://github.com/MPourjam/TACTIC-Pipeline/pkgs/container/tactic-pipeline)

**TACTIC** is an incremental, privacy-aware, and user-friendly pipeline designed for 16S rRNA amplicon sequence analysis. It provides a comprehensive solution for microbial community profiling by integrating four complementary analysis modes:

1. **TIC (Taxonomy-Informed Clustering)**: Novel hybrid clustering utilizing taxonomic boundaries to prevent diversity inflation while preserving true biological variation (*Default*).
2. **TAC (Taxonomy-Agnostic Clustering)**: De novo sequence-identity-based clustering that operates independently of reference taxonomies.
3. **ZOTU (Zero-radius OTUs / Denoised Amplicons)**: High-resolution denoising powered by the UNOISE3 algorithm.
4. **OTU (Operational Taxonomic Units)**: Standard 97% sequence-similarity clustering for backward compatibility.

TACTIC can be run via **Docker** (recommended for all platforms) or directly as a **Python CLI** on Linux x86_64 systems. Results are fully compatible with downstream statistical and visualization suites such as Namco and Rhea.

---

## ⚡ Quick Start (3 Steps)

Run the entire pipeline in 3 simple steps:

### Step 1: Organize Your Data

Place your forward and reverse FASTQ files (e.g. `Sample1_R1.fastq.gz`, `Sample1_R2.fastq.gz`) inside a working directory:

```bash
mkdir -p my_project/fastqs
# Copy or symlink your FASTQ files into my_project/fastqs
```

> 💡 **Supported Paired Read Naming Conventions**:
> - `_R1` / `_R2` or `-R1` / `-R2` (case-insensitive, e.g. `Sample_R1_001.fastq.gz`, `Sample-r1.fq`)
> - `_read1` / `_read2` (case-insensitive, e.g. `Sample_read1.fq.gz`, `Sample_read2.fq.gz`)
> - `_forward` / `_reverse` (case-insensitive, e.g. `Sample_forward.fastq`, `Sample_reverse.fastq`)
> - SRA / ENA numbered reads: `_1` / `_2` directly before file extension (e.g. `SRR123_1.fastq.gz`, `SRR123_2.fastq.gz`)
> - `_F` / `_R` or `-F` / `-R` directly before file extension (e.g. `Sample_F.fastq.gz`, `Sample-r.fq`)
> - `_F_` / `_R_` directly before Illumina machine chunk tail (e.g. `Sample_F_001.fastq.gz`, `Sample_R_001.fastq.gz`)

### Step 2: Run with Docker

```bash
docker run --rm \
  -v "/path/to/my_project:/base/inputs" \
  -v "/path/to/databases:/databases" \
  ghcr.io/mpourjam/tactic-pipeline:latest
```

> 💡 **Tip — Database Volume (`-v ...:/databases`)**:  
> Always mount a persistent local folder to `/databases`. On your first run, TACTIC will automatically download and index the SILVA reference database into this directory. Subsequent runs will reuse it immediately, saving 15–30 minutes per analysis!

### Step 3: View Your Results

All outputs are automatically generated in a timestamped analysis directory inside your project folder:

```text
my_project/
├── fastqs/
├── Analysis_20261015_080000/
│   ├── Pipeline_log.txt              # Complete preprocessing parameters and log
│   ├── Analysis_log.txt              # Detailed clustering and taxonomy log
│   ├── README_Taxonomy_Levels.txt    # Explanation of standard vs 7-level/8-level taxonomy outputs
│   ├── TIC_Pipeline/
│   │   ├── SOTUs-Table-TIC.tab       # Final abundance table for TIC (classic 7 levels)
│   │   ├── SOTUs-Table-TIC-8Levels.tab # Final abundance table with SILVA 144+ Kingdom (8 levels)
│   │   ├── SOTUs-Krona.html          # Interactive Krona taxonomy plot
│   │   └── SOTUs-Krona-8Levels.html  # Interactive Krona plot with bacterial Kingdom ring
│   ├── taxed_ZOTUs.fasta             # Representative sequences with SILVA taxonomy
│   ├── text.krona.html               # Interactive browser-based Krona taxonomy plot
│   ├── R1-per_base_quality.png       # Forward read quality assessment
│   └── reads_report.txt              # Read counts through filtering, merging, and dereplication
```

### 🏷️ SILVA 144+ Dual Taxonomy Outputs (`-7Levels` and `-8Levels`)

In SILVA release 144 onward, a formal **Kingdom** rank ending in `-ati` (e.g. *Bacillati*, *Pseudomonadati*, *Fusobacteriati*) was introduced between Domain and Phylum.

TACTIC automatically detects when SILVA $\ge 144$ is in use and generates parallel versions of all taxonomy-containing files (abundance tables, FASTAs, and Krona HTML charts):

| File Suffix | Analysis Pipelines | Taxonomic Levels | Description & Recommended Use |
| :--- | :--- | :--- | :--- |
| *(None / standard)* | All modes (`TIC`, `TAC`, `ZOTU`, `OTU`) | 7 levels: `Domain;Phylum;Class;Order;Family;Genus;Species/SOTU;` | **Backward compatibility**: Standard 7 ranks (Kingdom temporarily excluded). Use for downstream tools strictly expecting classic 7-level taxonomies (e.g., PICRUSt2, legacy QIIME / R scripts). |
| `*-8Levels.*` | `TIC` & `TAC` Pipelines | 8 levels: `Domain;Kingdom;Phylum;Class;Order;Family;Genus;Species/SOTU;` | **Full SILVA 144+ taxonomy**: Includes the bacterial Kingdom rank (`*-ati`) plus species/SOTU level. Recommended for detailed taxonomic profiling, Rhea, Namco, and 8-ring Krona charts. |
| `*-7Levels.*` | `ZOTU` & `OTU` Pipelines (and `-iz`) | 7 levels: `Domain;Kingdom;Phylum;Class;Order;Family;Genus;` | **Full rank resolution without species**: Includes the bacterial Kingdom rank (`*-ati`) up to genus level. |

> 💡 **Notice**: A descriptive `README_Taxonomy_Levels.txt` file is automatically placed in your `Analysis_*` results folder detailing the level mapping and tool recommendations.

---

## 📦 Installation Options

### Option A: Docker (Recommended)

Docker packages all external bioinformatics binaries (`usearch`, `bowtie2`, `sina`, `krona`) and Python dependencies into a self-contained image.

1. **Pull the latest image**:
   ```bash
   docker pull ghcr.io/mpourjam/tactic-pipeline:latest
   ```

2. Alternatively, download the archived `image.tar.gz` from the [GitHub Releases page](https://github.com/MPourjam/TACTIC-Pipeline/releases) and load it locally:
   ```bash
   docker load -i image.tar.gz
   ```

### Option B: Local Python Installation (Linux x86_64)

If you prefer running without Docker on a native Linux x86_64 machine:

```bash
# 1. Clone the repository
git clone https://github.com/MPourjam/TACTIC-Pipeline.git
cd TACTIC-Pipeline

# 2. Create and activate a virtual environment
python3 -m venv tactic_env
source tactic_env/bin/activate

# 3. Install dependencies
pip install --upgrade pip
pip install -r requirements.txt

# 4. Verify installation
python run_tactic.py --help
```

---

## 📖 How-To Guide & Common Use Cases

### 1. Run Multiple Analysis Modes Simultaneously

You can run one or multiple clustering modes in a single pass. Supported modes are `TIC`, `TAC`, `ZOTU`, and `OTU`:

```bash
docker run --rm \
  -v "/path/to/my_project:/base/inputs" \
  -v "/path/to/databases:/databases" \
  ghcr.io/mpourjam/tactic-pipeline:latest \
  --analysis-mode TIC TAC ZOTU
```

### 2. Narrow Down to a Specific FastQ Subdirectory

If your project directory has multiple sequencing runs or folders, use `--fastq-directory` (specified relative to your input directory):

```bash
docker run --rm \
  -v "/path/to/my_project:/base/inputs" \
  -v "/path/to/databases:/databases" \
  ghcr.io/mpourjam/tactic-pipeline:latest \
  --fastq-directory Sequencing_Run_1/
```

### 3. Customize Pipeline Parameters

To customize primer trimming lengths, error thresholds, or clustering parameters:

1. **Generate template configuration files**:
   ```bash
   docker run --rm \
     -v "/path/to/my_project:/base/inputs" \
     ghcr.io/mpourjam/tactic-pipeline:latest \
     --place-template-files
   ```
   This places `TACTICPipeline_args_TEMPLATE.yml` and `mapping_file_TEMPLATE.csv` into your project directory.

2. **Edit the YAML file**:
   Rename `TACTICPipeline_args_TEMPLATE.yml` to `my_params.yml` and adjust parameters as desired (e.g. primer trimming lengths, maximum expected errors).

3. **Run using your custom configuration**:
   ```bash
   docker run --rm \
     -v "/path/to/my_project:/base/inputs" \
     -v "/path/to/databases:/databases" \
     ghcr.io/mpourjam/tactic-pipeline:latest \
     --yml-file my_params.yml
   ```

### 4. Sample Mapping & Spike-In Absolute Quantification

If your sequencing run includes synthetic spike-in standards for absolute cell abundance quantification, provide a tab-separated mapping file:

1. **Format of `mapping_file.csv`** (tab-delimited):
   ```text
   #SampleID	total_weight_in_g	spike_amount	parent_path
   Sample_01	1.0	nan	
   Sample_02	1.2	50.0	
   Sample_03	0.9	50.0	
   ```
   - `#SampleID`: Unique sample base name matching your FASTQ file prefix.
   - `total_weight_in_g`: Sample biomass / weight (g). Use `nan` or any positive number if unspiked.
   - `spike_amount`: Spike amount added in nanograms (`ng`). If a sample was not spiked, set to `nan`.
   - `parent_path`: Optional subfolder path if duplicate sample names exist across subdirectories.

2. **Execute with mapping file**:
   ```bash
   docker run --rm \
     -v "/path/to/my_project:/base/inputs" \
     -v "/path/to/databases:/databases" \
     ghcr.io/mpourjam/tactic-pipeline:latest \
     --mapping-file mapping_file.csv
   ```

### 5. Multi-Threading for High Performance

TACTIC supports parallel sample preprocessing. Adjust CPU thread allocation with `--threads`:

```bash
docker run --rm \
  -v "/path/to/my_project:/base/inputs" \
  -v "/path/to/databases:/databases" \
  ghcr.io/mpourjam/tactic-pipeline:latest \
  --threads 16
```

---

## 📁 File Structure & Results Interpretation

After the pipeline finishes, your `Analysis_<timestamp>/` results directory contains:

| Output File / Directory | Description |
| :--- | :--- |
| **`Pipeline_log.txt`** | Complete log of the pipeline run, preprocessing parameters, and version details. |
| **`Analysis_log.txt`** | Log of the OTU/ZOTU clustering, taxonomic mapping, and spike normalization steps. |
| **`reads_report.txt`** | Table reporting input, merged, filtered, dereplicated, and retained read counts per sample. |
| **`taxed_ZOTUs.fasta`** | FASTA file of representative amplicon sequences with SILVA taxonomic annotations in headers. |
| **`silva_start_end.txt`** | Estimated 16S rRNA gene alignment coordinates for your sequenced amplicon region. |
| **`text.krona.html`** | Interactive, standalone hierarchical Krona chart viewable in any web browser. |
| **`R1-per_base_quality.png`** | Per-base quality score distribution across forward reads. |
| **`R2-per_base_quality.png`** | Per-base quality score distribution across reverse reads. |
| **`TIC_Pipeline/SOTUs-Table-TIC.tab`** | Final SOTU abundance table generated by Taxonomy-Informed Clustering. |
| **`TAC_Pipeline/OTUs-Table-TAC.tab`** | Final OTU abundance table generated by Taxonomy-Agnostic Clustering. |
| **`ZOTU_Pipeline/ZOTUs-Table.tab`** | Denoised amplicon (ZOTU) table generated by the UNOISE3 algorithm. |
| **`OTU_Pipeline/OTUs-Table.tab`** | Standard 97% sequence-identity clustered OTU table. |

---

## 🛠️ Command-Line Interface (CLI) Reference

```text
usage: run_tactic.py [-h] [-i INPUT_DIRECTORY] [-d FASTQ_DIRECTORY]
                     [-am ANALYSIS_MODE [ANALYSIS_MODE ...]] [-y YML_FILE]
                     [-map MAPPING_FILE] [-stat SPIKE_STAT] [-ut USEARCH_BIN]
                     [-db DB_DIRECTORY] [-spk-ref SPIKES_REFERENCES_DIR] [-sp]
                     [-sa] [-fp] [-tf] [-iz] [-snbf] [-t THREADS]

options:
  -h, --help            Show this help message and exit.
  -i INPUT_DIRECTORY, --input-directory INPUT_DIRECTORY
                        Base directory for all pipeline inputs and outputs.
                        Every needed file and directory should be accessible
                        relative to this directory (Default: current directory).
  -d FASTQ_DIRECTORY, --fastq-directory FASTQ_DIRECTORY
                        Directory containing FASTQ files to be processed;
                        relative to <--input-directory>. (Default: <--input-directory>).
  -am ANALYSIS_MODE [ANALYSIS_MODE ...], --analysis-mode ANALYSIS_MODE [ANALYSIS_MODE ...]
                        Choose one or more analysis modes (space-separated or
                        comma-separated). Supported: TIC, TAC, ZOTU, OTU.
                        (Default: 'TIC').
  -y YML_FILE, --yml-file YML_FILE
                        Path to pipeline configuration YAML file, relative to
                        <--input-directory>.
  -map MAPPING_FILE, --mapping-file MAPPING_FILE
                        Path to mapping file defining sample metadata, sample
                        weight, and spike amount (relative to <--input-directory>).
  -stat SPIKE_STAT, --spike-stat SPIKE_STAT
                        Path to compiled spike statistics file.
  -ut USEARCH_BIN, --usearch-bin USEARCH_BIN
                        Path to USEARCH v11 binary. (Default: bundled 64-bit binary).
  -db DB_DIRECTORY, --db-directory DB_DIRECTORY
                        Path to directory containing reference databases
                        (SILVA, SortMeRNA). (Default: /databases).
  -spk-ref SPIKES_REFERENCES_DIR, --spikes-references-dir SPIKES_REFERENCES_DIR
                        Directory containing spike-in reference FASTA files.
  -sp, --skip-preprocess
                        Skip sample preprocessing and use existing preprocessed data.
  -sa, --skip-analysis  Skip clustering and analysis step (only run preprocessing).
  -fp, --force-preprocess
                        Force re-preprocessing even if existing processed data is detected.
  -tf, --place-template-files
                        Write template YAML and mapping CSV files to <--input-directory>
                        and exit.
  -iz, --individual-zotus
                        Generate individual sample ZOTU tables in preprocessing
                        (slower, more detailed).
  -snbf, --skip-non-bacterial-filter
                        Skip filtering of non-bacterial 16S sequences.
  -t THREADS, --threads THREADS
                        Number of threads to use for parallel processing (Default: 16).
```

---

## ⚙️ Configuration File Details (`TACTICPipeline_args.yml`)

Parameters in `TACTICPipeline_args.yml` allow fine-grained control over read processing and clustering:

### 1. Read Merging (`MergePairsArgs`)
- **`fastq_maxdiffs`**: Maximum alignment mismatches allowed in overlapping region (default: `50`).
- **`fastq_pctid`**: Minimum percentage identity for merging read pairs (default: `50`).
- **`fastq_minmergelen`**: Minimum merged read length in base pairs (default: `200`).
- **`fastq_maxmergelen`**: Maximum merged read length in base pairs (default: `600`).

### 2. Primer Trimming (`TrimBothSidesArgs` / `TrimOneSideArgs`)
- **`stripleft`**: Bases to strip from 5' end, matching forward primer (default: `17`, corresponding to **341F** length).
- **`stripright`**: Bases to strip from 3' end, matching reverse primer (default: `21`, corresponding to **785R** length).

*Reference*: [ProbeBase 341F/785R details](https://probebase.csb.univie.ac.at/pb_report/probe/3844)

### 3. Read Quality Filtering (`FilterBothSidesArgs` / `FilterOneSideArgs`)
- **`fastq_maxee_rate`**: Maximum expected errors per base allowed (default: `0.002`). Reads with expected errors exceeding `(length * maxee_rate)` are discarded.
- **`fastq_truncqual`**: Deprecated in versions >=0.7.3 to avoid premature truncation.

### 4. Clustering Parameters
- **`ClusterZOTUsArgs`** (UNOISE3):
  - `minsize`: Minimum abundance threshold; amplicons with fewer reads are filtered out (default: `2`).
  - `match_id`: Identity threshold for read assignment (default: `0.99`).
  - `abund_limit`: Minimum relative abundance threshold across all samples (default: `0.000`).
- **`ComplexTICArgs`** (Taxonomy-Informed Clustering):
  - `family_sim`: Similarity cutoff for family-level grouping (default: `0.90`).
  - `genus_sim`: Similarity cutoff for genus-level grouping (default: `0.95`).
  - `species_sim`: Similarity cutoff for species-level clustering (default: `0.987`).
- **`TACClusterArgs`** (Taxonomy-Agnostic Clustering):
  - `cluster_thr`: Sequence identity threshold for de novo clustering (default: `0.987`).
- **`DeNovoClusterOTUsArgs`** (Standard OTUs):
  - `minsize`: Minimum amplicon abundance threshold (default: `8`).
  - `abund_limit`: Low-abundance filter threshold (default: `0.0025`).

---

## 📚 References & Citation

If you use TACTIC in your research, please cite:

1. **TACTIC Pipeline**:
   > Pourjam M., et al. *TACTIC: Taxonomy-Agnostic and Taxonomy-Informed Clustering for High-Precision 16S rRNA Amplicon Analysis.* (Manuscript in preparation).
2. **TIC Algorithm**:
   > Lagkouvardos I., et al. (2023). *Taxonomy-informed clustering improves the resolution of 16S rRNA amplicon sequencing data analysis.* **Frontiers in Bioinformatics**, 2:864597. [doi:10.3389/fbinf.2022.864597](https://doi.org/10.3389/fbinf.2022.864597)
3. **UNOISE Algorithm**:
   > Edgar R.C. (2016). *UNOISE2: improved error-correction for Illumina 16S and ITS amplicon sequencing.* **bioRxiv**, 081257. [doi:10.1101/081257](https://doi.org/10.1101/081257)
4. **SINA Classifier & SILVA Database**:
   > Pruesse E., et al. (2012). *SINA: Accurate high-throughput multiple sequence alignment of ribosomal RNA genes.* **Bioinformatics**, 28(14):1823–1829. [doi:10.1093/bioinformatics/bts252](https://doi.org/10.1093/bioinformatics/bts252)
