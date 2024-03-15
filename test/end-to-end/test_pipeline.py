import contextlib
import subprocess
from typing import Callable, Tuple, Optional
# from contextlib import suppress

import zipfile
import pytest
import os
import glob
import tempfile


IMAGE_NAME = "tic-pipeline-test"
data = "test/data/truncated"


@pytest.fixture(scope="session")
def build_image():
    print("Building the image using file system image")
    # build the image using tar file target
    # docker buildx build -f Dockerfile.IMNGS2Pipeline --output type=tar,dest=image.tar .
    subprocess.run(["docker", "buildx", "build", "-f", "Dockerfile.IMNGS2Pipeline", "--output", "type=tar,dest=image.tar", "."])
    # load with tag
    subprocess.run(["docker", "image", "import", "image.tar", IMAGE_NAME])


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

    def __init__(self, entries: list[Entry]):
        self.entries = entries

    def write(self, file: str, sep: str = "\t"):
        with open(file, "w") as f:
            f.write(sep.join(self.HEADER_ENTRIES) + "\n")
            for entry in self.entries:
                f.write(entry.csv_line(sep) + "\n")


def check_expected_files_exist(analysis_folder: str):
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
        "Analysis_log.txt",
        "Map-GOTU-FOTU.tab",
        "Map-SOTU-GOTU.tab",
        "Map-ZOTU-SOTU.tab",
        "SOTUs-Seqs.fasta",
        "SOTUs-Table.tab",
        "SOTUs-Tree-nj.tre",
        "spike_mapping_file.csv",
        "ZOTUs-Seqs.fasta",
        "ZOTUs-Table.tab",
        "ZOTUs-Tree-nj.tre",
    ]
    for file in files_to_be_there:
        # Checking the format of analysis folder and if it's a zip file we assert that file exist in the zip file
        if analysis_folder.endswith(".zip") and os.path.isfile(os.path.join(analysis_folder, file)):
            with zipfile.ZipFile(analysis_folder, "r") as z:
                try:
                    z.getinfo(file)
                except KeyError:
                    return_bool = False
                break
        else:
            if not os.path.exists(os.path.join(analysis_folder, file)):
                return_bool = False
                break

    return return_bool


def run_container(setup_file_structure: Callable[[str], Tuple[Optional[str], Optional[str]]]) -> None:
    # set up a temporary directory to run the pipeline
    with contextlib.nullcontext(tempfile.mkdtemp()) as run_dir:
        print(f"{run_dir = }")

        # copy the data to the temporary directory
        mapping_file, yaml_file = setup_file_structure(run_dir)

        cmd = ["docker", "run", "--rm", "-v", f"{run_dir}:/base/inputs", IMAGE_NAME, "run_imngs2"]

        if mapping_file:
            cmd += ["--mapping-file", os.path.basename(mapping_file)]
        if yaml_file:
            cmd += ["--yml-file", os.path.basename(yaml_file)]

        # Write the command to a file
        with open(f"{run_dir}/cmd.txt", "w") as f:
            f.write(" ".join(cmd))
        print(f"{cmd = }")
        # run docker cmd with the image
        result = subprocess.run(
            cmd,
            capture_output=False, text=True
        )

        # print the output of the container
        # print(f"{result.stdout = }")

        assert result.returncode == 0

        # get latest analysis folder. those are marked with a timestamp
        latest_analysis_folder_date = 0
        latest_analysis_folder_time = 0
        latest_analysis_folder = ""
        for analysis_folder in glob.glob(f"{run_dir}/Analysis_*"):
            print(f"{analysis_folder = }")
            _, date, time = str(os.path.basename(analysis_folder).split(".")[0]).split("_")
            date = int(date)
            time = int(time)

            # latest folder has the highest date and time
            if date > latest_analysis_folder_date or (date == latest_analysis_folder_date and time > latest_analysis_folder_time):
                latest_analysis_folder_date = date
                latest_analysis_folder_time = time
                latest_analysis_folder = analysis_folder

        print(f"{latest_analysis_folder = }")
        assert latest_analysis_folder != ""
        file_set_complete = check_expected_files_exist(latest_analysis_folder)
        assert file_set_complete


def test_flat_autodiscover(build_image):
    def setup_file_structure(run_dir: str):
        # copy the data to the temporary directory
        subprocess.run(["cp", "-r", data, run_dir])
        # move the contents of the data folder to the run_dir
        for file in glob.glob(f"{run_dir}/{os.path.basename(data)}/*"):
            subprocess.run(["mv", file, run_dir])
        # rm the empty folder
        subprocess.run(["rmdir", f"{run_dir}/{os.path.basename(data)}"])
        return None, None

    run_container(setup_file_structure)


def test_subfolder_autodiscover(build_image):
    def setup_file_structure(run_dir: str):
        # make a subfolder and copy the data to the temporary directory
        subfolder = os.path.join(run_dir, "subfolder")
        os.makedirs(subfolder)
        subprocess.run(["cp", "-r", data, subfolder])
        return None, None

    run_container(setup_file_structure)


def test_flat_mapping(build_image):
    def setup_file_structure(run_dir: str):
        # copy the contents of data to the temporary directory
        subprocess.run(["cp", "-r", data, run_dir])
        # move the contents of the data folder to the run_dir
        for file in glob.glob(f"{run_dir}/{os.path.basename(data)}/*"):
            subprocess.run(["mv", file, run_dir])
        # rm the empty folder
        subprocess.run(["rmdir", f"{run_dir}/{os.path.basename(data)}"])

        # create a mapping file
        mapping_file = os.path.join(run_dir, "mapping_file.csv")
        samples = [
            MappingFile.Entry("truncSRR13005876_S1_L001", 1.0, 6, ""),
            MappingFile.Entry("truncSRR13005987_S2_L001", 2.0, 6, ""),
        ]
        mapping = MappingFile(samples)
        mapping.write(mapping_file)

        return mapping_file, None

    run_container(setup_file_structure)


def test_subfolder_mapping(build_image):
    def setup_file_structure(run_dir: str):
        # copy the contents of data to the temporary directory
        subprocess.run(["cp", "-r", data, run_dir])
        sf_name = os.path.basename(data)

        # create a mapping file
        mapping_file = os.path.join(run_dir, "mapping_file.csv")
        samples = [
            MappingFile.Entry("truncSRR13005876_S1_L001", 1.0, 6, sf_name),
            MappingFile.Entry("truncSRR13005987_S2_L001", 2.0, 6, sf_name),
        ]
        mapping = MappingFile(samples)
        mapping.write(mapping_file)

        return mapping_file, None

    run_container(setup_file_structure)
