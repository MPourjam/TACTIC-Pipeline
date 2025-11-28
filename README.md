# TACTIC Pipeline
TACTIC is an incremental, privacy-aware, and user-friendly tool designed for 16S rRNA amplicon sequence analysis, offering a comprehensive and versatile solution to address the challenges of microbial community profiling. It integrates four distinct pipelines—OTU (legacy clustering), zOTU (denoised amplicons), and the novel hybrid approaches, Taxonomy Agnostic Clustering (TAC) and Taxonomy Informed Clustering (TIC)—to provide researchers with flexible options, including conventional methods and advanced strategies that mitigate diversity inflation inherent to denoising. The platform is accessible via a containerized command-line interface, a stand-alone installable software with a graphical user interface, and a browser-based web tool, ensuring accessibility, reproducibility, and scalability for users of all computing environments and skill levels, while its results are compatible with downstream analysis tools like Namco and Rhea.


## Installation
You can run the TACTIC pipeline either by cloning the repository and running the `run_tactic.py` script or via docker image.
### 1. Clone the repository
You can clone this repository as below:
```bash
git clone https://github.com/MPourjam/TACTIC-Pipeline.git
```
or you can alternatively download the repository and unzip it.

Then you need to create a virtual environment and install the required packages in the environment:

```bash
cd TACTIC-Pipeline
python3 -m venv tactic_env
source tactic_env/bin/activate
pip install -r requirements.txt
```

> As the binary files are for linux and CPU architecture of x86, this method only works on __linux__ operating system and
CPU architecture of __x86__.

### 2. Docker (recommended)
There is a docker image for each release available that you can pull and start using. For example:

```bash
docker pull ghcr.io/mpourjam/tactic-pipeline:0.7.2
```

or you can download the exported docker image (`image.tar.gz`) of each version in [the releases page](https://github.com/MPourjam/TACTIC-Pipeline/releases). Then import it to your local docker image repository. 


```bash
docker load -i image.tar.gz
# or
docker import image.tar.gz
```

## Command Line Interface (CLI)
The pipeline is given as single python script (`run_tactic.py`) which you can run by `python3 </path/to/run_tactic.py>`. In order to get the help text on the arguments run `python3 </path/to/run_tactic.py> --help`. Below you find the the CLI arguments output from the help of the pipeline.

```bash
usage: run_tactic.py [-h] [-i INPUT_DIRECTORY] [-d FASTQ_DIRECTORY]
                     [-am {de-novo,TAC,TIC,ZOTU}] [-y YML_FILE]
                     [-map MAPPING_FILE] [-stat SPIKE_STAT]
                     [-ut USEARCH_BIN] [-db DB_DIRECTORY]
                     [-spk-ref SPIKES_REFERENCES_DIR] [-sp] [-sa] [-fp]
                     [-tf] [-iz] [-snbf] [-t THREADS]

options:
  -h, --help            show this help message and exit
  -i INPUT_DIRECTORY, --input-directory INPUT_DIRECTORY
                        Base directory for all pipeline inputs and outputs;
                        every needed file and directory should be findable
                        relative to this directory.
  -d FASTQ_DIRECTORY, --fastq-directory FASTQ_DIRECTORY
                        Directory containing FASTQ files to be processed;
                        should be relative to <--input-directory>. Default
                        is <--input-directory>.
  -am {de-novo,TAC,TIC,ZOTU}, --analysis-mode {de-novo,TAC,TIC,ZOTU}
                        Choose an analysis mode. Default is 'TIC'.
  -y YML_FILE, --yml-file YML_FILE
                        Path to pipeline configuration file in YAML format.
                        This file defines the arguments for the pipeline and
                        should be relative to <--input-directory>.
  -map MAPPING_FILE, --mapping-file MAPPING_FILE
  -stat SPIKE_STAT, --spike-stat SPIKE_STAT
                        The path to a mapping file defining spike count,
                        sample weight and spike amount for each sample.
                        #SampleID SpikeReads spikes_total_weight_in_g
                        spike_amount parent_path. Relative to
                        <--input-directory>
  -ut USEARCH_BIN, --usearch-bin USEARCH_BIN
                        Path to binary of usearch version 11. Default is
                        usearch11.0.667_i86linux64.
  -db DB_DIRECTORY, --db-directory DB_DIRECTORY
                        Path to directory containing required databases
                        (e.g., SILVA, SortMeRNA).
  -spk-ref SPIKES_REFERENCES_DIR, --spikes-references-dir SPIKES_REFERENCES_DIR
                        Directory containing FASTA files of spike-in
                        references (required for spike removal). Should be
                        relative to <--input-directory>. Default is
                        {SPIKESIDX}.
  -sp, --skip-preprocess
                        Should skip preprocessing step
  -sa, --skip-analysis  Should skip analysis step
  -fp, --force-preprocess
                        Force preprocessing samples. It invalidates
                        --skip-preprocess argument.
  -tf, --place-template-files
                        Writes the default argument yaml template file
                        (TACTICPipeline_args.yml) and mapping template file
                        (mapping_file.csv) to <--input-directory>, print
                        help text and exits.
  -iz, --individual-zotus
                        If set, the preprocessing step will produce
                        individual ZOTUs tables for each sample (slower,
                        more detailed output). If not set, only a combined
                        table is produced, which speeds up the pipeline.
  -snbf, --skip-non-bacterial-filter
                        Skip filtering of non-bacterial 16S sequences.
  -t THREADS, --threads THREADS
                        Number of threads to use for parallel processing
                        (default: 16). Recommended range: 1 to 24.
```

## TACTIC Pipeline Configuration Arguments
Anaylsis arguments are defined in a YAML file named `TACTICPipeline_args.yml`. You can modify the arguments to your desired values in the file and give the file to the pipeline to have the arguments taken into account.

### Preprocessing Parameters
These arguemtns affects the pre-processing step. Each (pair) of fastq file goes into pre-processing step. In this pre-processing step:

1. Spike removal
2. Quality check
3. Reads merge (if paired)
4. Trimming
5. Filtering
6. Dereplication

If `-iz / --individual-zotus` flag is set to true then:

7. ZOTUs clustering
8. Filtering for bacterial/archaeal 16S rRNA reads (if opted for by giving `` flag to the pipeline)
9. ZOTU table creation
10. Taxonomy assignment
11. Krona graph generation


#### MergePairsArgs (Adjustable from GUI's expert mode)
The arguments below are set to most lenient argument to merge as many as possible reads. Later in the filtering step the merged reads with bad quality get filtered out.
- **fastq_maxdiffs**: Maximum number of mismatches allowed in the alignment of overlapping region (default: 50)
- **fastq_pctid**: Minimum percentage identity for merging paired reads (default: 50)
- **fastq_minmergelen**: Minimum length of merged reads in base pairs (default: 200)
- **fastq_maxmergelen**: Maximum length of merged reads in base pairs (default: 600)

#### TrimBothSidesArgs (Adjustable from GUI's basic mode)
- **stripleft**: Number of bases to remove from the left side, corresponding to __341F__ primer length (default: 17)
- **stripright**: Number of bases to remove from the right side, corresponding to __785R__ primer length (default: 21)

Reference: [www.probebase.csb.univie.ac.at](https://probebase.csb.univie.ac.at/pb_report/probe/3844)

#### TrimOneSideArgs
- **stripleft**: Number of bases to remove from the left side, corresponding to __341F__ primer length (default: 17)

Reference: [www.probebase.csb.univie.ac.at](https://probebase.csb.univie.ac.at/pb_report/probe/3844)

#### FilterBothSidesArgs (Adjustable from GUI's expert mode)
- **fastq_maxee_rate**: Discard reads with > E expected errors per base.Calculated after any truncation options have been applied. For example, with the fastq_maxee_rate option set to 0.01, then a read of length 100 will be discarded if the expected errors is >1, and a read of length 1,000 will be discarded if the expected errors is >10. (default: 0.002)
- **fastq_truncqual**: Quality score threshold for truncating reads; truncation occurs at first position with Q score ≤ N (default: 10)

> __fastq_truncqual__ is deprecated in versions >=0.7.3 as it could cause drop of some merged reads due to short length after truncation.

Reference: [USEARCH manual for read filtering](https://www.drive5.com/usearch/manual/cmd_fastq_filter.html)

#### FilterOneSideArgs (Adjustable from GUI's expert mode)
- **fastq_maxee_rate**: Discard reads with > E expected errors per base.Calculated after any truncation options have been applied. For example, with the fastq_maxee_rate option set to 0.01, then a read of length 100 will be discarded if the expected errors is >1, and a read of length 1,000 will be discarded if the expected errors is >10. (default: 0.002)
- **fastq_truncqual**: Quality score threshold for truncating reads; truncation occurs at first position with Q score ≤ N (default: 10)

> __fastq_truncqual__ is deprecated in versions >=0.7.3 as it could cause drop of some merged reads due to short length after truncation.

Reference: [USEARCH manual for read filtering](https://www.drive5.com/usearch/manual/cmd_fastq_filter.html)

#### IndividualSampleClusterZOTUsArgs (Adjustable from GUI's basic mode)
- **minsize**: Minimum abundance threshold for merged sequences clustering; sequences with fewer reads are filtered out and not taken to the clustering step (default: 2)

Reference: [USEARCH manual for UNOISE algorithm (ZOTU clustering)](https://www.drive5.com/usearch/manual/cmd_unoise3.html)

#### Filter16SArgs
- **e**: E-value threshold for filtering 16S sequences (default: 0.1)
- **num_alignments**: Number of alignments to report per query (default: 1)

#### BuildZOTUTableArgs
- **match_id**: Minimum sequence identity threshold for ZOTU matching; higher values result in more refined ZOTUs and lower total count of ZOTUs (default: 0.99)

#### AddTaxArgs
- **turn**: Specifies which sequences to add taxonomy information to (default: all)

### Analysis Parameters

#### ClusterZOTUsArgs (Adjustable from GUI's expert mode)
- **minsize**: Minimum abundance threshold; amplicons with fewer reads are filtered out (default: 2)
- **abund_limit**: Abundance limit threshold for filtering low-abundance sequences. ZOTUs with relative abundance below this threshold in all samples will be discarded (default: 0.000)
- **sample_wise_correction**: Enable or disable sample-wise abundance correction (default: false)
- **match_id**: Minimum sequence identity threshold for ZOTU clustering (default: 0.99)

Reference: [USEARCH manual for UNOISE algorithm (ZOTU clustering)](https://www.drive5.com/usearch/manual/cmd_unoise3.html)

#### DeNovoClusterOTUsArgs (Adjustable from GUI's expert mode)
- **minsize**: Minimum abundance threshold for de novo OTU clustering (default: 8)
- **abund_limit**: Abundance limit threshold for filtering low-abundance sequences. ZOTUs with relative abundance below this threshold in all samples will be discarded (default: 0.0025)
- **sample_wise_correction**: Enable or disable sample-wise abundance correction (default: false)

#### ComplexTICArgs (Adjustable from GUI's expert mode)
- **family_sim**: Sequence similarity threshold for family-level taxonomic clustering (default: 0.90)
- **genus_sim**: Sequence similarity threshold for genus-level taxonomic clustering (default: 0.95)
- **species_sim**: Sequence similarity threshold for species-level taxonomic clustering (default: 0.987)
- **abund_limit**: Abundance limit threshold for filtering low-abundance sequences. ZOTUs with relative abundance below this threshold in all samples will be discarded. (not yet implemented) (default: 0.0025)
- **sample_wise_correction**: Enable or disable sample-wise abundance correction. Abudnance values below this threshold will be rounded to __0__ in the S/ZOTUs tables (not yet implemented) (default: false)


References:
1. [TIC paper](https://doi.org/10.3389/fbinf.2022.864597)
2. [Ticlust python package](https://pypi.org/project/ticlust/)
3. [TIC scripts](https://github.com/Lagkouvardos/TIC)


#### TACClusterArgs (Adjustable from GUI's expert mode)
- **abund_limit**: Abundance limit threshold (not yet implemented) (default: 0.0025)
- **sample_wise_correction**: Enable or disable sample-wise abundance correction (not yet implemented) (default: false)
- **cluster_thr**: Clustering threshold for sequence similarity (default: 0.987)


## Tutorial
You can run the pipeline either by directly running the script (`run_tactic.py`) or the docker image.

### Quickstart
You can run the pipeline quickly with _fastq files auto-discovery_ and _default pipeline arguments_. In this mode the pipeline will automatically discover the fastq files (and pair them), creates a mapping file and runs the pipeline with default arguments.

__Running via Script__:
```bash
cd TACTIC-Pipeline
python run_tactic.py --input-dir /path/to/your/project/directory
```
__Running via Docker__:

```bash
docker run --rm -v "/path/to/your/project/directory:/base/inputs" ghcr.io/mpourjam/tactic-pipeline:0.7.3 
```

#### Run with custom pipeline arguments
The pipeline reads the arguments from a [YAML](https://en.wikipedia.org/wiki/YAML) file it expects as input to `-y` or `--yml-file`. To have a template of such YAML file and fill in your desired arguments you can run the pipeline as described in [Quickstart](#quickstart) and with the switch `-tf` or `--place-template-files`. This will put template files of the argument YAML file and the mapping file (described in the next section) in your project directory.

__Running via Script__:
```bash
cd TACTIC-Pipeline
python run_tactic.py --input-dir /path/to/your/project/directory --place-template-files
# or
python run_tactic.py --input-dir /path/to/your/project/directory -y
```
__Running via Docker__:

```bash
docker run --rm -v "/path/to/your/project/directory:/base/inputs" ghcr.io/mpourjam/tactic-pipeline:0.7.3 --place-template-files
# or
docker run --rm -v "/path/to/your/project/directory:/base/inputs" ghcr.io/mpourjam/tactic-pipeline:0.7.3 -y
```

These commands will put two files in your project directory:
1. `mapping_file_TEMPLATE.csv`: This is the template file to give samples in a mapping file with their sample weight and spike amount for spike normalization.
2. `TACTICPipeline_args.yml`: This is a YAML file containing the arguments to the pipeline.

#### Run with custom set of your fastq files

### Inputs
To run the pipeline with default 
### Outputs
### Quickstart
# References
