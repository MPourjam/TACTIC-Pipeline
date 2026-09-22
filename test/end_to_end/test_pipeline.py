import contextlib
import subprocess
from typing import Callable, Tuple, Optional, List
# from contextlib import suppress
import urllib.request as urequest
import zipfile
import pytest
import os
import glob
import tempfile
import re
import csv


global data, IMAGE_NAME
IMAGE_NAME = "tic-pipeline-test"
data = "test/data"
usearch_11_32_url = "https://drive5.com/downloads/usearch11.0.667_i86linux32.gz"


@pytest.fixture(scope="session")
def build_image():
    if os.path.isfile("/.dockerenv"):
        # We are already testing in the docker container
        return
    print("Building the image using file system image")
    # build the image using tar file target
    # docker buildx build -f Dockerfile.TACTICPipeline --output type=tar,dest=image.tar .
    subprocess.run(["docker", "buildx", "build", "-f", "Dockerfile.TACTICPipeline", "--output", "type=tar,dest=image.tar", "."])
    # load with tag
    subprocess.run(["docker", "image", "import", "image.tar", IMAGE_NAME])


class ArgumentSet:
    # Some class-wide arguments
    container_input_dir: str = "/base/inputs"
    yml_file: str = ""  # path to the yml file
    mapping_file: str = ""  # path to the mapping file
    spike_stats_file: str = ""  # path to the spike stats file
    usearch_bin: str = ""  # path to the usearch binary
    db_directory: str = ""  # path to the database directory
    skip_preprocess = False
    skip_analysis = False
    place_template_file = False
    threads = 0
    analysis_mode = ""
    individual_zotus = False
    spikes_references_dir = None

    def __init__(self, input_dir, fastq_dir, **kwargs):
        self.input_dir = input_dir
        self.fastq_dir = fastq_dir
        for key, value in kwargs.items():
            setattr(self, key, value)

    def to_args(self):
        args = []
        if self.input_dir:
            args += ["--input-directory", self.input_dir]
        if self.fastq_dir:
            args += ["--fastq-directory", self.fastq_dir]
        if self.yml_file:
            args += ["--yml-file", self.yml_file]
        if self.mapping_file:
            args += ["--mapping-file", self.mapping_file]
        if self.spike_stats_file:
            args += ["--spike-stat", self.spike_stats_file]
        if self.usearch_bin:
            args += ["--usearch-bin", self.usearch_bin]
        if True:
            args += ["--db-directory", self.db_directory]
        if self.skip_preprocess:
            args += ["--skip-preprocess"]
        if self.skip_analysis:
            args += ["--skip-analysis"]
        if self.place_template_file:
            args = ["--input-directory", self.input_dir, "--fastq-directory", self.fastq_dir]
            args += ["--place-template-files"]
        if self.threads:
            args += ["--threads", str(self.threads)]
        if self.analysis_mode:
            modes = normalize_analysis_modes(self.analysis_mode)
            args += ["--analysis-mode", *modes]
        if self.individual_zotus:
            args += ["-iz"]
        if self.spikes_references_dir:
            args += ["--spikes-references-dir", self.spikes_references_dir]

        return args

    def __repr__(self) -> str:
        return f"ArgumentSet({', '.join([f'{k}={v}' for k, v in self.__dict__.items()])})"

    def __str__(self) -> str:
        return f"ArgumentSet({', '.join([f'{k}={v}' for k, v in self.__dict__.items()])})"


class MappingFile:
    HEADER_ENTRIES = ['#SampleID', 'total_weight_in_g', 'spike_amount', 'parent_path']

    class Entry:
        def __init__(self, sample_id: str, total_weight_in_g: float, spike_amount: float, parent_path: str):
            self.sample_id = sample_id
            self.total_weight_in_g = total_weight_in_g
            self.spike_amount = spike_amount
            self.parent_path = parent_path

        def csv_line(self, sep: str) -> str:
            return sep.join([self.sample_id, str(self.total_weight_in_g), str(self.spike_amount), self.parent_path])

    def __init__(self, entries: List[Entry]):
        self.entries = entries
        self.last_written_file = None

    def write(self, file: str, sep: str = "\t"):
        with open(file, "w") as f:
            f.write(sep.join(self.HEADER_ENTRIES) + "\n")
            for entry in self.entries:
                f.write(entry.csv_line(sep) + "\n")
        self.last_written_file = file


def normalize_analysis_modes(analysis_mode) -> List[str]:
    """
    Normalize one or more analysis modes to a unique list accepted by run_tactic.
    Supports legacy token 'de-novo' by mapping it to 'OTU'.
    """
    if isinstance(analysis_mode, (list, tuple, set)):
        raw_items = [str(item).strip() for item in analysis_mode]
    else:
        raw_items = [str(analysis_mode).strip()] if analysis_mode else ["TIC"]

    modes = []
    for raw_mode in raw_items:
        if not raw_mode:
            continue
        for mode in [m.strip().upper() for m in re.split(r"[,\s]+", raw_mode) if m.strip()]:
            if mode == "DE-NOVO":
                mode = "OTU"
            if mode not in modes:
                modes.append(mode)

    return modes if modes else ["TIC"]


def check_expected_preprocessed_files(preprocessed_dir: str, args: ArgumentSet) -> bool:
    """
    The preprocessed folder should contain the following files:
    ```
    .
    ├── spike_stat_mapping_file.csv  # Must be there
    ├── TACTICPipeline_args.yml  # Must be there
    ├── taxed_ZOTUs.fasta  # Must be there
    ├── silva_start_end.txt  # Must be there
    ├── ZOTUs-table.final.tab
    ├── div_metrics.tab
    ├── reads_report.txt
    ├── R1-per_base_quality.png
    ├── R2-per_base_quality.png
    ├── rank_abundance_plot.png
    ├── rarefaction_curve.png
    ├── text.krona.html
    ├── *.udb
    ├── aligned_*.fasta.tar.gz
    ├── derep.fasta.tar.gz
    ├── other.non16rRNA.fasta
    ```
    """
    return_bool = True
    files_to_be_there = [
        re.compile(r"spike_stat_mapping_file\.csv"),
        re.compile(r"TACTICPipeline_args\.yml"),
        re.compile(r"taxed_ZOTUs\.fasta"),  # also in zip file
        re.compile(r"silva_start_end\.txt"),
        re.compile(r"ZOTUs-table\.final\.tab"),  # also in zip file
        re.compile(r"reads_report\.txt"),
        # re.compile(r"derep\.fasta\.tar\.gz"),  # we do not compress the derep.fasta file anymore
        re.compile(r"derep\.fasta"),
        re.compile(r"R1-per_base_quality\.png"),  # also in zip file
        re.compile(r"R2-per_base_quality\.png"),  # also in zip file
        re.compile(r"aligned_.*\.fasta\.tar\.gz"),
        re.compile(r"text\.krona\.html"),  # also in zip file
        re.compile(r".*\.udb"),
        # NOTE below are some files which could not get produced if some R packages are not installed
        # re.compile(r"div_metrics\.tab"),
        # re.compile(r"other\.non16rRNA\.fasta"),
        # re.compile(r"rank_abundance_plot\.png"),  # also in zip file
        # re.compile(r"rarefaction_curve\.png"),  # also in zip file
    ]
    if not args.individual_zotus:
        files_to_be_there = [
            re.compile(r"spike_stat_mapping_file\.csv"),
            re.compile(r"TACTICPipeline_args\.yml"),
            re.compile(r"R1-per_base_quality\.png"),  # also in zip file
            re.compile(r"R2-per_base_quality\.png"),
            re.compile(r"reads_report\.txt"),
            re.compile(r"sorted\.fasta")
        ]
    existing_files_list = []
    # Checking the format of analysis folder and if it's a zip file we assert that file exist in the zip file
    if preprocessed_dir.endswith(".zip") and os.path.isfile(os.path.join(preprocessed_dir)):
        with zipfile.ZipFile(preprocessed_dir, "r") as z:
            existing_files_list = z.namelist()
        files_to_be_there = [files_to_be_there[ind] for ind in [2, 4, 7, 8, 10]]
    else:
        existing_files_list = os.listdir(preprocessed_dir)
    # use re search to check if patterns in files_to_be_there exist in the existing_files_list
    shoud_be_there = 0
    is_there = 0
    for file_there in files_to_be_there:
        shoud_be_there += 1
        for existing_f in existing_files_list:
            if re.search(file_there, existing_f):
                is_there += 1
                break

    return_bool = shoud_be_there == is_there

    return return_bool


def check_expected_files_exist(analysis_folder: str, analysis_mode: str = "TIC") -> bool:
    """
    Validate the analysis output folder with the new nested structure:
    Shared/, ZOTU_Pipeline/, OTU_Pipeline/, TAC_Pipeline/, TIC_Pipeline/
    and PICRUSt2-Inputs/.
    """
    modes = normalize_analysis_modes(analysis_mode)
    files_to_be_there = [
        re.compile(r"spike_mapping_file\.csv$"),
        re.compile(r"Analysis_log\.txt$"),
        re.compile(r"Shared/sorted\.fasta$"),
        re.compile(r"PICRUSt2-Inputs/.*[ZS]?OTUs?-Table.*-For-PICRUSt2\.tab$"),
        re.compile(r"PICRUSt2-Inputs/.*[ZS]?OTUs?-Seqs.*-For-PICRUSt2\.fasta$")
    ]

    # ZOTU outputs are expected explicitly for ZOTU mode and as TAC/TIC dependencies.
    if any(mode in {"ZOTU", "TAC", "TIC"} for mode in modes):
        files_to_be_there.extend(
            [
                re.compile(r"ZOTU_Pipeline/ZOTUs-Seqs\.fasta$"),
                re.compile(r"ZOTU_Pipeline/ZOTUs-Table\.tab$"),
                re.compile(r"ZOTU_Pipeline/ZOTUs-Krona\.html$"),
                re.compile(r"ZOTU_Pipeline/SampleMinCount_Normalized.*ZOTUs-Table.*\.tab$"),
            ]
        )

    if "OTU" in modes:
        files_to_be_there.extend(
            [
                re.compile(r"OTU_Pipeline/OTUs-Seqs\.fasta$"),
                re.compile(r"OTU_Pipeline/OTUs-Table\.tab$"),
                re.compile(r"OTU_Pipeline/OTUs-Krona\.html$"),
                re.compile(r"OTU_Pipeline/SampleMinCount_Normalized.*OTUs-Table.*\.tab$"),
            ]
        )

    if "TAC" in modes:
        files_to_be_there.extend(
            [
                re.compile(r"TAC_Pipeline/OTUs-Table-TAC\.tab$"),
                re.compile(r"TAC_Pipeline/OTUs-Seqs-TAC\.fasta$"),
                re.compile(r"TAC_Pipeline/OTUs-Krona\.html$"),
                re.compile(r"TAC_Pipeline/Map-ZOTUs-OTUs-TAC\.tab$"),
                re.compile(r"TAC_Pipeline/SampleMinCount_Normalized.*OTUs-Table-TAC.*\.tab$"),
            ]
        )

    if "TIC" in modes:
        files_to_be_there.extend(
            [
                re.compile(r"TIC_Pipeline/Map-FOTU-GOTU\.tab$"),
                re.compile(r"TIC_Pipeline/Map-GOTU-SOTU\.tab$"),
                re.compile(r"TIC_Pipeline/Map-SOTU-ZOTU\.tab$"),
                re.compile(r"TIC_Pipeline/SOTUs-Table-TIC\.tab$"),
                re.compile(r"TIC_Pipeline/SOTUs-Seqs-TIC\.fasta$"),
                re.compile(r"TIC_Pipeline/SOTUs-Krona\.html$"),
                re.compile(r"TIC_Pipeline/Invalid-Tax-Seqs-TIC\.fasta$"),
                re.compile(r"TIC_Pipeline/.*FullTaxonomy\.tab$"),
                re.compile(r"TIC_Pipeline/SampleMinCount_Normalized.*SOTUs-Table-TIC.*\.tab$"),
            ]
        )

    existing_files_list = []
    # Checking the format of analysis folder and if it's a zip file we assert that file exist in the zip file
    # If it's a zip file and the path exists, list all entries inside the zip.
    if analysis_folder.endswith(".zip") and os.path.isfile(analysis_folder):
        with zipfile.ZipFile(analysis_folder, "r") as z:
            existing_files_list = z.namelist()
    else:
        # Walk the directory recursively and collect file paths relative to analysis_folder
        existing_files_list = []
        for root, _, files in os.walk(analysis_folder):
            for fname in files:
                relpath = os.path.relpath(os.path.join(root, fname), analysis_folder)
                existing_files_list.append(relpath)
    missing_patterns = []
    for file_pat in files_to_be_there:
        if not any(re.search(file_pat, existing_f) for existing_f in existing_files_list):
            missing_patterns.append(file_pat.pattern)

    return len(missing_patterns) == 0


def check_otutable_format(otu_table: str, argument_set: ArgumentSet, run_dir: str) -> bool:
    # Find all files named "spike_mapping_file.csv in fastq_dir which has a parent folder ending in "__processed"
    fastq_dir = os.path.join(os.path.abspath(run_dir), argument_set.fastq_dir)
    finished_processes = []
    mapping_file_path = os.path.join(os.path.abspath(run_dir), argument_set.mapping_file)
    if not os.path.isfile(mapping_file_path):
        spike_mapping_files = glob.glob(fastq_dir + "**/*__processed/**/spike_mapping_file.csv", recursive=True)
        for spike_map_f in spike_mapping_files:
            # in case there is a file named taxed_ZOTUs.fasta in the same folder as spike_mapping_file.csv
            if os.path.isfile(os.path.join(os.path.dirname(spike_map_f), "taxed_ZOTUs.fasta")):
                with open(spike_map_f, "r") as f:
                    reader = csv.DictReader(f, delimiter="\t")
                    for row in reader:
                        finished_processes.append(row["#SampleID"])
    else:
        with open(mapping_file_path, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                finished_processes.append(row["#SampleID"])
    # if all samples in finished_processes are in the otu_table header then return true
    with open(otu_table, "r") as f:
        header = f.readline().strip().split("\t")
        for sample in finished_processes:
            if sample not in header:
                return False

    return True


def run_container(setup_file_structure: Callable[[str], Tuple[Optional[str], Optional[str]]]) -> None:
    # set up a temporary directory to run the pipeline
    with contextlib.nullcontext(tempfile.mkdtemp()) as run_dir:
        print("run_dir =", run_dir)

        # copy the data to the temporary directory
        arg_set = setup_file_structure(run_dir)

        # if we are in a container created with image name IMAGE_NAME, the input directory we run run_tactic in /base/run_tactic.py
        # cmd = ["docker", "run", "--rm", "-v", run_dir + ":/base/inputs", IMAGE_NAME, "python", "/base/run_tactic.py"]
        if os.path.isfile("/.dockerenv"):  # then we are in a container
            cmd = ["python", "/base/run_tactic.py"]
            arg_set.input_dir = run_dir
        else:
            cmd = ["docker", "run", "--rm", "-v", run_dir + ":/base/inputs", IMAGE_NAME, "run_tactic"]
        cmd += arg_set.to_args()
        # Write the command to a file
        with open(run_dir + "/cmd.txt", "w") as f:
            f.write(" ".join(cmd) + "\n")
        print("cmd =", cmd)
        fastq_dir = os.path.join(run_dir, arg_set.fastq_dir)
        # Snapshot analysis directories before running the command so we can
        # reliably pick the output folder created by this invocation only.
        pre_existing_analysis_dirs = set(glob.glob(os.path.join(fastq_dir, "Analysis_*")))
        # run docker cmd with the image
        result = subprocess.run(
            cmd,
            capture_output=False, text=True
        )

        file_set_complete = False
        if not arg_set.skip_analysis:
            # Pick analysis folder(s) created by this run. Fall back to
            # timestamp-based selection if needed for backward compatibility.
            post_analysis_dirs = set(glob.glob(os.path.join(fastq_dir, "Analysis_*")))
            created_analysis_dirs = sorted(post_analysis_dirs - pre_existing_analysis_dirs)
            latest_analysis_folder_date = 0
            latest_analysis_folder_time = 0
            latest_analysis_folder = ""
            analysis_candidates = created_analysis_dirs if created_analysis_dirs else list(post_analysis_dirs)
            for analysis_folder in analysis_candidates:
                print("analysis_folder =", analysis_folder)
                _, date, time = str(os.path.basename(analysis_folder).split(".")[0]).split("_")
                date = int(date)
                time = int(time)
                # latest folder has the highest date and time
                if date > latest_analysis_folder_date or (date == latest_analysis_folder_date and time > latest_analysis_folder_time):
                    latest_analysis_folder_date = date
                    latest_analysis_folder_time = time
                    latest_analysis_folder = analysis_folder
            print("latest_analysis_folder =", latest_analysis_folder)
            assert latest_analysis_folder != ""
            selected_modes = normalize_analysis_modes(arg_set.analysis_mode)
            if any(mode in {"ZOTU", "TAC", "TIC"} for mode in selected_modes):
                assert check_otutable_format(
                    os.path.join(latest_analysis_folder, "ZOTU_Pipeline", "ZOTUs-Table.tab"),
                    arg_set,
                    fastq_dir
                ), "Some samples are missing in the ZOTU table"
            if "OTU" in selected_modes:
                assert check_otutable_format(
                    os.path.join(latest_analysis_folder, "OTU_Pipeline", "OTUs-Table.tab"),
                    arg_set,
                    fastq_dir
                ), "Some samples are missing in the OTU table"
            if "TAC" in selected_modes:
                assert check_otutable_format(
                    os.path.join(latest_analysis_folder, "TAC_Pipeline", "OTUs-Table-TAC.tab"),
                    arg_set,
                    fastq_dir
                ), "Some samples are missing in the TAC OTU table"
            if "TIC" in selected_modes:
                assert check_otutable_format(
                    os.path.join(latest_analysis_folder, "TIC_Pipeline", "SOTUs-Table-TIC.tab"),
                    arg_set,
                    fastq_dir
                ), "Some samples are missing in the SOTU table"

            file_set_complete = check_expected_files_exist(latest_analysis_folder, selected_modes)
            assert file_set_complete, f"Some files are missing in the analysis folder for modes: {selected_modes}"
            # Test return code
            assert result.returncode == 0, f"Command failed with return code {result.returncode}"
        elif not arg_set.skip_preprocess and arg_set.skip_analysis:
            file_set_complete_counter = 0
            dirs_checked = 0
            # find all directories with the regex pattern of r'.*__processed/[0-9]{8}_[0-9]{6}' in the run_dir and use re.search
            for preprocessed_folder in glob.glob(fastq_dir + "/**/*__processed/", recursive=True):
                # list directories in preprocessing folder and check if they match the pattern
                for tmstmp_dir in glob.glob(preprocessed_folder + "/*", recursive=True):
                    if re.search(r"[0-9]{8}_[0-9]{6}", tmstmp_dir):
                        dirs_checked += 1
                        print("preprocessed_folder =", tmstmp_dir)
                        file_set_complete_counter += int(check_expected_preprocessed_files(tmstmp_dir, arg_set))
            file_set_complete = bool(dirs_checked > 0 and file_set_complete_counter == dirs_checked)
            assert file_set_complete, "Some files are missing in the preprocessed folder"
            # Test return code
            # List of valid codes here: https://github.com/MPourjam/TACTIC-Pipeline/issues/42
            assert result.returncode == 169, f"Command failed with return code {result.returncode}, expected return code is 169 for skipped analysis"
        elif arg_set.skip_preprocess and arg_set.skip_analysis:
            print("No analysis or preprocessing was done")
            file_set_complete = False
            assert file_set_complete, "No analysis or preprocessing was done"
            assert result.returncode == 169, f"Command failed with return code {result.returncode}, expected return code is 169 for skipped preprocessing and analysis"


def copy_initial_files(run_dir: str) -> None:
    """
    Copy the fastq files from the data directory to the run directory.
    The data directory is expected to be in the same directory as this script.
    """
    # move the contents of the data folder to the run_dir
    for file in glob.glob(data + "/truncated/*"):
        subprocess.run(["cp", file, str(run_dir) + "/"])
    subprocess.run(["cp", data + "/TACTICPipeline_args_test.yml", run_dir])


def copy_spikes_dir(run_dir: str) -> None:
    """
    Copy the spikes directory from the data directory to the run directory.
    The data directory is expected to be in the same directory as this script.
    """
    subprocess.run(["cp", "-r", data + "/custom_spikes", run_dir + "/custom_spikes"])


@pytest.mark.xfail(reason="check_otutable_format() checks complete processing only by existence of spike_mapping_file.csv")
def test_flat_autodiscover(build_image):
    def setup_file_structure(run_dir: str):
        copy_initial_files(run_dir)
        # setting args
        args_set = ArgumentSet("", "")

        return args_set

    run_container(setup_file_structure)


def test_subfolder_autodiscover(build_image):
    def setup_file_structure(run_dir: str):
        # make a subfolder and copy the data to the temporary directory
        subfolder = os.path.join(run_dir, "subfolder")
        os.makedirs(subfolder)
        subprocess.run(["cp", "-r", data, subfolder])
        # setting args
        args_set = ArgumentSet("", "subfolder")
        return args_set

    run_container(setup_file_structure)


def test_flat_mapping(build_image):
    def setup_file_structure(run_dir: str):
        # copy the contents of data to the temporary directory
        copy_initial_files(run_dir)

        # setting args
        args_set = ArgumentSet("", "")
        # create a mapping file
        mapping_file = os.path.join(run_dir, "mapping_file.csv")
        samples = [
            MappingFile.Entry("truncSRR13005876_S1_L001", 1.0, 6, ""),
            MappingFile.Entry("truncSRR13005987_S2_L001", 2.0, 6, ""),
        ]
        mapping = MappingFile(samples)
        mapping.write(mapping_file)
        args_set.mapping_file = "mapping_file.csv"

        return args_set

    run_container(setup_file_structure)


def test_preprocessing(build_image):
    def setup_file_structure(run_dir: str):
        copy_initial_files(run_dir)
        # setting args
        args_set = ArgumentSet("", "")
        # create a mapping file
        mapping_file = os.path.join(run_dir, "mapping_file.csv")
        samples = [
            MappingFile.Entry("truncSRR13005876_S1_L001", 1.0, 6, ""),
            MappingFile.Entry("truncSRR13005987_S2_L001", 2.0, 6, ""),
        ]
        mapping = MappingFile(samples)
        mapping.write(mapping_file)
        args_set.mapping_file = "mapping_file.csv"
        args_set.skip_analysis = True

        return args_set

    run_container(setup_file_structure)


@pytest.mark.xfail(reason="usearch is not a bin file and it should eixt with error code 161")
def test_custom_usearch_bin_arg(build_image):
    def setup_file_structure(run_dir: str):
        copy_initial_files(run_dir)
        # setting args
        args_set = ArgumentSet("", "")
        # create a mapping file
        mapping_file = os.path.join(run_dir, "mapping_file.csv")
        samples = [
            MappingFile.Entry("truncSRR13005876_S1_L001", 1.0, 6, ""),
            MappingFile.Entry("truncSRR13005987_S2_L001", 2.0, 6, ""),
        ]
        mapping = MappingFile(samples)
        mapping.write(mapping_file)
        args_set.mapping_file = "mapping_file.csv"
        # create a fake bin file in run_dir
        bin_file = os.path.join(run_dir, "usearch_custom")
        with open(bin_file, "w") as f:
            f.write("FAKE USEARCH BIN")
        args_set.usearch_bin = "usearch_custom"
        args_set.skip_analysis = True

        return args_set

    run_container(setup_file_structure)


def test_custom_usearch_bin_preprocessing(build_image):
    def setup_file_structure(run_dir: str):
        copy_initial_files(run_dir)
        # setting args
        args_set = ArgumentSet("", "")
        # create a mapping file
        mapping_file = os.path.join(run_dir, "mapping_file.csv")
        samples = [
            MappingFile.Entry("truncSRR13005987_S2_L001", 2.0, 6, ""),
        ]
        mapping = MappingFile(samples)
        mapping.write(mapping_file)
        args_set.mapping_file = "mapping_file.csv"
        # create a fake bin file in run_dir
        # downloading usearch from usearch_11_32_url and passing it as usearch11_custom
        urequest.urlretrieve(usearch_11_32_url, os.path.join(run_dir, "usearch11.0.667_i86linux32.gz"))
        args_set.usearch_bin = "usearch11.0.667_i86linux32.gz"
        args_set.skip_analysis = True

        return args_set

    run_container(setup_file_structure)


def test_custom_usearch_bin_full_analysis(build_image):
    def setup_file_structure(run_dir: str):
        copy_initial_files(run_dir)
        # setting args
        args_set = ArgumentSet("", "")
        # create a mapping file
        mapping_file = os.path.join(run_dir, "mapping_file.csv")
        samples = [
            MappingFile.Entry("truncSRR13005876_S1_L001", 1.0, 6, ""),
            MappingFile.Entry("truncSRR13005987_S2_L001", 2.0, 6, ""),
        ]
        mapping = MappingFile(samples)
        mapping.write(mapping_file)
        args_set.mapping_file = "mapping_file.csv"
        # create a fake bin file in run_dir
        # downloading usearch from usearch_11_32_url and passing it as usearch11_custom
        urequest.urlretrieve(usearch_11_32_url, os.path.join(run_dir, "usearch11.0.667_i86linux32.gz"))
        args_set.usearch_bin = "usearch11.0.667_i86linux32.gz"

        return args_set

    run_container(setup_file_structure)


def test_thread_arg(build_image):
    def setup_file_structure(run_dir: str):
        copy_initial_files(run_dir)
        # setting args
        args_set = ArgumentSet("", "")
        # create a mapping file
        mapping_file = os.path.join(run_dir, "mapping_file.csv")
        samples = [
            MappingFile.Entry("truncSRR13005876_S1_L001", 1.0, 6, ""),
            MappingFile.Entry("truncSRR13005987_S2_L001", 2.0, 6, ""),
        ]
        mapping = MappingFile(samples)
        mapping.write(mapping_file)
        args_set.mapping_file = "mapping_file.csv"
        args_set.threads = 1

        return args_set

    run_container(setup_file_structure)


def test_full_denovo_analysis(build_image):
    def setup_file_structure(run_dir: str):
        copy_initial_files(run_dir)
        # setting args
        args_set = ArgumentSet("", "")
        # create a mapping file
        mapping_file = os.path.join(run_dir, "mapping_file.csv")
        samples = [
            MappingFile.Entry("truncSRR13005876_S1_L001", 1.0, 6, ""),
            MappingFile.Entry("truncSRR13005987_S2_L001", 2.0, 6, ""),
        ]
        mapping = MappingFile(samples)
        mapping.write(mapping_file)
        args_set.mapping_file = "mapping_file.csv"
        args_set.analysis_mode = "OTU"
        args_set.yml_file = "TACTICPipeline_args_test.yml"  # path to the yml file

        return args_set

    run_container(setup_file_structure)


def test_full_multi_analysis_mode_selector(build_image):
    def setup_file_structure(run_dir: str):
        copy_initial_files(run_dir)
        # setting args
        args_set = ArgumentSet("", "")
        # create a mapping file
        mapping_file = os.path.join(run_dir, "mapping_file.csv")
        samples = [
            MappingFile.Entry("truncSRR13005876_S1_L001", 1.0, 6, ""),
            MappingFile.Entry("truncSRR13005987_S2_L001", 2.0, 6, ""),
        ]
        mapping = MappingFile(samples)
        mapping.write(mapping_file)
        args_set.mapping_file = "mapping_file.csv"
        # Multiple-mode selector: OTU + TAC in one run.
        # TAC implies ZOTU dependency outputs should also exist.
        args_set.analysis_mode = ["OTU", "TAC"]
        args_set.yml_file = "TACTICPipeline_args_test.yml"  # path to the yml file

        return args_set

    run_container(setup_file_structure)


def test_full_TAC_analysis(build_image):
    def setup_file_structure(run_dir: str):
        copy_initial_files(run_dir)
        # setting args
        args_set = ArgumentSet("", "")
        # create a mapping file
        mapping_file = os.path.join(run_dir, "mapping_file.csv")
        samples = [
            MappingFile.Entry("truncSRR13005876_S1_L001", 1.0, 6, ""),
            MappingFile.Entry("truncSRR13005987_S2_L001", 2.0, 6, ""),
        ]
        mapping = MappingFile(samples)
        mapping.write(mapping_file)
        args_set.mapping_file = "mapping_file.csv"
        args_set.analysis_mode = "TAC"
        args_set.yml_file = "TACTICPipeline_args_test.yml"  # path to the yml file

        return args_set

    run_container(setup_file_structure)


def test_full_TIC_analysis_no_iz(build_image):
    def setup_file_structure(run_dir: str):
        copy_initial_files(run_dir)
        # setting args
        args_set = ArgumentSet("", "")
        # create a mapping file
        mapping_file = os.path.join(run_dir, "mapping_file.csv")
        samples = [
            MappingFile.Entry("truncSRR13005876_S1_L001", 1.0, 6, ""),
            MappingFile.Entry("truncSRR13005987_S2_L001", 2.0, 6, ""),
        ]
        mapping = MappingFile(samples)
        mapping.write(mapping_file)
        args_set.mapping_file = "mapping_file.csv"
        args_set.analysis_mode = "TIC"
        args_set.yml_file = "TACTICPipeline_args_test.yml"  # path to the yml file

        return args_set

    run_container(setup_file_structure)


def test_full_TIC_analysis_iz(build_image):
    def setup_file_structure(run_dir: str):
        copy_initial_files(run_dir)
        # setting args
        args_set = ArgumentSet("", "")
        # create a mapping file
        mapping_file = os.path.join(run_dir, "mapping_file.csv")
        samples = [
            MappingFile.Entry("truncSRR13005876_S1_L001", 1.0, 6, ""),
            MappingFile.Entry("truncSRR13005987_S2_L001", 2.0, 6, ""),
        ]
        mapping = MappingFile(samples)
        mapping.write(mapping_file)
        args_set.mapping_file = "mapping_file.csv"
        args_set.analysis_mode = "TIC"
        args_set.individual_zotus = True
        args_set.yml_file = "TACTICPipeline_args_test.yml"  # path to the yml file

        return args_set

    run_container(setup_file_structure)


def test_full_TIC_analysis_custom_spike(build_image):
    def setup_file_structure(run_dir: str):
        copy_initial_files(run_dir)
        copy_spikes_dir(run_dir)
        # setting args
        args_set = ArgumentSet("", "")
        # create a mapping file
        mapping_file = os.path.join(run_dir, "mapping_file.csv")
        samples = [
            MappingFile.Entry("truncSRR13005876_S1_L001", 1.0, 6, ""),
            MappingFile.Entry("truncSRR13005987_S2_L001", 2.0, 6, ""),
        ]
        mapping = MappingFile(samples)
        mapping.write(mapping_file)
        args_set.mapping_file = "mapping_file.csv"
        args_set.analysis_mode = "TIC"
        args_set.individual_zotus = False
        args_set.spikes_references_dir = "custom_spikes"
        args_set.yml_file = "TACTICPipeline_args_test.yml"  # path to the yml file

        return args_set

    run_container(setup_file_structure)


def test_full_TIC_analysis_custom_database_path(build_image):
    def setup_file_structure(run_dir: str):
        copy_initial_files(run_dir)
        copy_spikes_dir(run_dir)
        # setting args
        args_set = ArgumentSet("", "")
        # create a mapping file
        mapping_file = os.path.join(run_dir, "mapping_file.csv")
        samples = [
            MappingFile.Entry("truncSRR13005876_S1_L001", 1.0, 6, ""),
            MappingFile.Entry("truncSRR13005987_S2_L001", 2.0, 6, ""),
        ]
        mapping = MappingFile(samples)
        mapping.write(mapping_file)
        args_set.mapping_file = "mapping_file.csv"
        args_set.analysis_mode = "TIC"
        args_set.individual_zotus = False
        # args_set.spikes_references_dir = "custom_spikes"
        args_set.yml_file = "TACTICPipeline_args_test.yml"  # path to the yml file
        args_set.db_directory = "My/custom/database/path"  # If absolute path then not relative to INPUT_DIR
        args_set.skip_analysis = True

        return args_set

    run_container(setup_file_structure)
