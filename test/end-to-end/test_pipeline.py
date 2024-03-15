import contextlib
import subprocess
from typing import Callable, Tuple, Optional, List
# from contextlib import suppress

import zipfile
import pytest
import os
import glob
import tempfile
import re
import shutil


IMAGE_NAME = "tic-pipeline-test"
data = "test/data/truncated"


@pytest.fixture(scope="session")
def build_image():
    if os.path.isfile("/.dockerenv"):
        # We are already testing in the docker container
        return
    print("Building the image using file system image")
    # build the image using tar file target
    # docker buildx build -f Dockerfile.IMNGS2Pipeline --output type=tar,dest=image.tar .
    subprocess.run(["docker", "buildx", "build", "-f", "Dockerfile.IMNGS2Pipeline", "--output", "type=tar,dest=image.tar", "."])
    # load with tag
    subprocess.run(["docker", "image", "import", "image.tar", IMAGE_NAME])


class ArgumentSet:
    # Some class-wide arguments
    container_input_dir: str = "/base/inputs"
    yml_file: str = None  # path to the yml file
    mapping_file: str = None  # path to the mapping file
    spike_stats_file: str = None  # path to the spike stats file
    usearch_bin: str = None  # path to the usearch binary
    db_directory: str = None  # path to the database directory
    skip_preprocess = False
    skip_analysis = False
    place_template_file = False
    threads = 0

    def __init__(self, input_dir, fastq_dir, **kwargs):
        self.input_dir = input_dir
        self.fastq_dir = fastq_dir
        for key, value in kwargs.items():
            setattr(self, key, value)

    def to_args(self):
        args = []
        if self.input_dir:
            args += ["--input-dir", self.input_dir]
        if self.fastq_dir:
            args += ["--fastq-dir", self.fastq_dir]
        if self.yml_file:
            args += ["--yml-file", self.yml_file]
        if self.mapping_file:
            args += ["--mapping-file", self.mapping_file]
        if self.spike_stats_file:
            args += ["--spike-stats-file", self.spike_stats_file]
        if self.usearch_bin:
            args += ["--usearch-bin", self.usearch_bin]
        if self.db_directory:
            args += ["--db-directory", self.db_directory]
        if self.skip_preprocess:
            args += ["--skip-preprocess"]
        if self.skip_analysis:
            args += ["--skip-analysis"]
        if self.place_template_file:
            args = ["--input-dir", self.input_dir, "--fastq-dir", self.fastq_dir]
            args += ["--place-template-file"]
        if self.threads:
            args += ["--threads", str(self.threads)]
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

    def write(self, file: str, sep: str = "\t"):
        with open(file, "w") as f:
            f.write(sep.join(self.HEADER_ENTRIES) + "\n")
            for entry in self.entries:
                f.write(entry.csv_line(sep) + "\n")


def check_expected_preprocessed_files(preprocessed_dir: str) -> bool:
    """
    The preprocessed folder should contain the following files:
    ```
    .
    ├── spike_stat_mapping_file.csv  # Must be there
    ├── IMNGS2Pipeline_args.yml  # Must be there
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
        re.compile(r"IMNGS2Pipeline_args\.yml"),
        re.compile(r"taxed_ZOTUs\.fasta"),  # also in zip file
        re.compile(r"silva_start_end\.txt"),
        re.compile(r"ZOTUs-table\.final\.tab"),  # also in zip file
        re.compile(r"reads_report\.txt"),
        re.compile(r"derep\.fasta\.tar\.gz"),
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


def check_expected_files_exist(analysis_folder: str) -> bool:
    """
    The analysis folder should contain the following files:
    ```
    .
    ├── Analysis_log.txt
    ├── Map-GOTU-FOTU.tab
    ├── Map-SOTU-GOTU.tab
    ├── Map-ZOTU-SOTU.tab
    ├── SOTUs-Seqs.fasta
    ├── SOTUs-Table.tab
    ├── SOTUs-Tree-nj.tre
    ├── spike_mapping_file.csv
    ├── ZOTUs-Seqs.fasta
    ├── ZOTUs-Table.tab
    └── ZOTUs-Tree-nj.tre
    ```
    """
    return_bool = True
    files_to_be_there = [
        re.compile(r"Analysis_log\.txt"),
        re.compile(r"Map-GOTU-FOTU\.tab"),
        re.compile(r"Map-SOTU-GOTU\.tab"),
        re.compile(r"Map-ZOTU-SOTU\.tab"),
        re.compile(r"SOTUs-Seqs\.fasta"),
        re.compile(r"SOTUs-Table\.tab"),
        re.compile(r"SOTUs-Tree-nj\.tre"),
        re.compile(r"spike_mapping_file\.csv"),
        re.compile(r"ZOTUs-Seqs\.fasta"),
        re.compile(r"ZOTUs-Table\.tab"),
        re.compile(r"ZOTUs-Tree-nj\.tre"),
    ]
    existing_files_list = []
    # Checking the format of analysis folder and if it's a zip file we assert that file exist in the zip file
    if analysis_folder.endswith(".zip") and os.path.isfile(os.path.join(analysis_folder)):
        with zipfile.ZipFile(analysis_folder, "r") as z:
            existing_files_list = z.namelist()
    else:
        existing_files_list = os.listdir(analysis_folder)
    # use re search to check if patterns in files_to_be_there exist in the existing_files_list
    for file_there in files_to_be_there:
        is_there = False
        for existing_f in existing_files_list:
            if re.search(file_there, existing_f):
                is_there = True
                break
        return_bool = return_bool and is_there

    return return_bool


def run_container(setup_file_structure: Callable[[str], Tuple[Optional[str], Optional[str]]]) -> None:
    # set up a temporary directory to run the pipeline
    with contextlib.nullcontext(tempfile.mkdtemp()) as run_dir:
        print("run_dir =", run_dir)

        # copy the data to the temporary directory
        arg_set = setup_file_structure(run_dir)

        # if we are in a container created with image name IMAGE_NAME, the input directory we run run_imngs2 in /base/run_imngs2.py
        # cmd = ["docker", "run", "--rm", "-v", run_dir + ":/base/inputs", IMAGE_NAME, "python", "/base/run_imngs2.py"]
        if os.path.isfile("/.dockerenv"):  # then we are in a container
            cmd = ["python", "/base/run_imngs2.py"]
            arg_set.input_dir = run_dir
        else:
            cmd = ["docker", "run", "--rm", "-v", run_dir + ":/base/inputs", IMAGE_NAME, "run_imngs2"]
        cmd += arg_set.to_args()

        # Write the command to a file
        with open(run_dir + "/cmd.txt", "w") as f:
            f.write(" ".join(cmd) + "\n")
        print("cmd =", cmd)
        # run docker cmd with the image
        result = subprocess.run(
            cmd,
            capture_output=False, text=True
        )

        assert result.returncode == 0
        file_set_complete = False
        fastq_dir = os.path.join(run_dir, arg_set.fastq_dir)
        if not arg_set.skip_analysis:
            # get latest analysis folder. those are marked with a timestamp
            latest_analysis_folder_date = 0
            latest_analysis_folder_time = 0
            latest_analysis_folder = ""
            for analysis_folder in glob.glob(fastq_dir + "/Analysis_*"):
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
            file_set_complete = check_expected_files_exist(latest_analysis_folder)
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
                        file_set_complete_counter += int(check_expected_preprocessed_files(tmstmp_dir))
            file_set_complete = bool(dirs_checked > 0 and file_set_complete_counter == dirs_checked)
        assert file_set_complete


def test_flat_autodiscover(build_image):
    def setup_file_structure(run_dir: str):
        # copy the data to the temporary directory
        subprocess.run(["cp", "-r", data, run_dir])
        # move the contents of the data folder to the run_dir
        for file in glob.glob(run_dir + "/" + os.path.basename(data) + "/*"):
            subprocess.run(["mv", file, run_dir])
        # rm the empty folder
        subprocess.run(["rmdir", run_dir + "/" + os.path.basename(data)])
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
        subprocess.run(["cp", "-r", data, run_dir])
        # move the contents of the data folder to the run_dir
        for file in glob.glob(run_dir + "/" + os.path.basename(data) + "/*"):
            subprocess.run(["mv", file, run_dir])
        # rm the empty folder
        shutil.rmtree(run_dir + "/" + os.path.basename(data))

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
        # copy the contents of data to the temporary directory
        subprocess.run(["cp", "-r", data, run_dir])
        # move the contents of the data folder to the run_dir
        for file in glob.glob(run_dir + "/" + os.path.basename(data) + "/*"):
            subprocess.run(["mv", file, run_dir])
        # rm the empty folder
        shutil.rmtree(run_dir + "/" + os.path.basename(data))
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


@pytest.mark.xfail
def test_custom_usearch_bin_arg(build_image):
    def setup_file_structure(run_dir: str):
        # copy the contents of data to the temporary directory
        subprocess.run(["cp", "-r", data, run_dir])
        # move the contents of the data folder to the run_dir
        for file in glob.glob(run_dir + "/" + os.path.basename(data) + "/*"):
            subprocess.run(["mv", file, run_dir])
        # rm the folder
        shutil.rmtree(run_dir + "/" + os.path.basename(data))
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
        # copy the contents of data to the temporary directory
        subprocess.run(["cp", "-r", data, run_dir])
        # move the contents of the data folder to the run_dir
        for file in glob.glob(run_dir + "/" + os.path.basename(data) + "/*"):
            subprocess.run(["mv", file, run_dir])
        # rm the folder
        shutil.rmtree(run_dir + "/" + os.path.basename(data))
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
        custom_usearch_path = os.path.join(run_dir, "usearch11_custom")
        shutil.copy("binaries/usearch_11_64", custom_usearch_path)
        args_set.usearch_bin = "usearch11_custom"
        args_set.skip_analysis = True

        return args_set

    run_container(setup_file_structure)


def test_custom_usearch_bin_full_analysis(build_image):
    def setup_file_structure(run_dir: str):
        # copy the contents of data to the temporary directory
        subprocess.run(["cp", "-r", data, run_dir])
        # move the contents of the data folder to the run_dir
        for file in glob.glob(run_dir + "/" + os.path.basename(data) + "/*"):
            subprocess.run(["mv", file, run_dir])
        # rm the folder
        shutil.rmtree(run_dir + "/" + os.path.basename(data))
        # setting args
        args_set = ArgumentSet("", "")
        # create a mapping file
        mapping_file = os.path.join(run_dir, "mapping_file.csv")
        samples = [
            MappingFile.Entry("truncSRR13005876_S1_L001", 1.0, 6, ""),
        ]
        mapping = MappingFile(samples)
        mapping.write(mapping_file)
        args_set.mapping_file = "mapping_file.csv"
        # create a fake bin file in run_dir
        custom_usearch_path = os.path.join(run_dir, "usearch11_custom")
        shutil.copy("binaries/usearch_11_64", custom_usearch_path)
        args_set.usearch_bin = "usearch11_custom"

        return args_set

    run_container(setup_file_structure)
