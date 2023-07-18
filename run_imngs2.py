import os
import argparse
import shutil
from ENA_Processor.processing_job import main_processing as preprocessing
from ENA_Processor.analyze_sample_job import main as analyze
import ENA_Processor.processing_helper as proc_helper
from multiprocessing import cpu_count, Pool
from sys import version_info
if version_info[0] < 3:
    from pathlib2 import Path, PurePath  # pip2 install pathlib2
else:
    from pathlib import Path, PurePath


global INPUT_DIR, DBS_DIR, POOL_SIZE, PREP_LOG
INPUT_DIR = "/srv/base/inputs/"
DBS_DIR = "/srv/base/databases/"
max_pool = int(cpu_count() * 0.3)
POOL_SIZE = max_pool if max_pool > 0 else 1
max_batch = int(POOL_SIZE * 0.8)
# Preparing the logger
PREP_LOG = proc_helper.gimmelogger()


def split_files_to_directories():
    pass


def run_preprocessing(seq_files_t: tuple, preproc_dir: str, args_yml_path: str):
    """
    It takes a tuple of paths to sequencing files.
    Create directory for basename of files and move
     them there. Then triggers the preprocessing step.
     It returns the path of direcotry in which the
     preprocessd files of sample are.
    """
    try:
        if not hasattr(seq_files_t, "__iter__") or len(seq_files_t) > 2:
            msg = "seq_files_t must be an iterable of max lenght 2"
            # PREP_LOG.error(msg)
            raise ValueError(msg)
        # Creating the preprocessing directory
        seq_files_t = [Path(PurePath(sfi)).absolute() for sfi in seq_files_t]
        # print(seq_files_t)
        new_paths = ["", ""]  # [ForwardNewPath, ReverseNewPath]
        new_dir = seq_files_t[0]
        # print(new_dir)
        # Getting initial stem
        while new_dir.suffixes:
            new_dir = Path(new_dir).absolute().parent.joinpath(new_dir.stem)
        new_dir = Path(PurePath(preproc_dir)).absolute().joinpath(proc_helper.get_base_name(new_dir))
        new_dir.mkdir(parents=True, exist_ok=True)
        for ind, sfi in enumerate(seq_files_t):
            new_path = new_dir.joinpath(sfi.name)
            new_paths.append(new_path)
            # Copying files
            shutil.copy2(str(sfi), str(new_path))
        preprocessing(
            input_dir=new_dir,
            paired="Yes" if len(seq_files_t) == 2 else "No",
            forward_file=new_paths[0],
            reverse_file=new_paths[1],
            input_id=new_dir.name,
            args_file_path=args_yml_path,
        )

    except Exception as exc:
        msg = f"{exc}"
        PREP_LOG.error(msg)
        shutil.rmtree(str(new_dir))
        PREP_LOG.warning("Deleting {}".format(new_dir))
        return ""


def run_analysis():
    pass


def run_imngs2(
        fastq_file_dir: str,
        args_yml_file: str,
        dbs_dir: str,
        only_preproc: bool = False):
    try:
        fastq_file_dir = Path(PurePath(fastq_file_dir))
        args_yml_file = Path(PurePath(args_yml_file))
        dbs_dir = Path(PurePath(dbs_dir))
        preproc_dir = Path(PurePath(INPUT_DIR)).joinpath("Preprocessing")
        # print(preproc_dir)
        preproc_dir.mkdir(parents=True, exist_ok=True)
        assert fastq_file_dir.is_dir(), "fastq_file_dir must be a path to directory"
        assert args_yml_file.is_file(), "args_yml_file must be a path to a file"
        assert dbs_dir.is_dir(), "dbs_dir must be a path to directory"
        # Gathering sequence files
        PREP_LOG.debug("Gathering sequence files in {}".format(str(fastq_file_dir)))
        seq_file_pairs = proc_helper.pair_seq_files(str(fastq_file_dir))
        PREP_LOG.info("{} sampels were collected from {}.".format(str(len(seq_file_pairs)), str(fastq_file_dir)))
        # TODO
        with Pool(POOL_SIZE, maxtasksperchild=1) as pool:
            res_list = [pool.apply_async(run_preprocessing, args=((file_pair), preproc_dir, str(args_yml_file),))
                        for file_pair in seq_file_pairs[:2]]  # TODO REMOVE slicing
            for res in res_list:
                res.wait()
        exit()
        # TODO gathering and running analysis

    except Exception as exc:
        PREP_LOG.error(exc)
        shutil.rmtree(str(preproc_dir))
        PREP_LOG.warning("Deleting {}".format(preproc_dir))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    help_text = "The directory to find recursively all fastq files inside"
    parser.add_argument("-d", "--fastq-directory",
                        type=str,
                        help=help_text,
                        default=".")
    parser.add_argument("-y", "--yml-file",
                        type=str,
                        help="Path to arguments yaml file",
                        default="./IMNGS2Pipeline_args.yml")
    parser.add_argument("-db", "--db-directory",
                        type=str,
                        help="Path to directory containing silva, sortmerna files",
                        default=DBS_DIR)
    parser.add_argument("-op", "--only-preprocess",
                        action="store_true",
                        help="Should run analysis step")
    args = parser.parse_args()
    if args.db_directory != DBS_DIR:
        DBS_DIR = args.db_directory
    run_imngs2(
        args.fastq_directory,
        args.yml_file,
        args.db_directory,
        args.only_preprocess
    )
