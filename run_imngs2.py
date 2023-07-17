import os
from ENA_Processor.processing_job import main_processing as preprocessing
from ENA_Processor.analyze_sample_job import main as analyze
import ENA_Processor.processing_helper as proc_helper
from sys import version_info
if version_info[0] < 3:
    from pathlib2 import Path, PurePath  # pip2 install pathlib2
else:
    from pathlib import Path, PurePath


def gather_files_in_pairs(targ_dir: str):
    """
    It walks through a directory and gathers fastq files in paires
     and returns a list of tuples of fastq files paths paires.
     Single files would have empty strings as pair.
    """
    output_list = []
    targ_dir = Path(PurePath(targ_dir))
    if not targ_dir.is_dir():
        return output_list

    for pd, dirs, fis in os.walk(targ_dir):
        continue
    pass


def split_files_to_directories():
    pass


def run_preprocessing_parallel():
    pass


def run_analysis():
    pass


def run_imngs2(
        fastq_file_dir: str,
        args_yml_file: str,
        dbs_dir: str,
        only_preproc: bool = False):
    fastq_file_dir = Path(PurePath(fastq_file_dir))
    args_yml_file = Path(PurePath(args_yml_file))
    dbs_dir = Path(PurePath(dbs_dir))
    assert fastq_file_dir.is_dir(), "fastq_file_dir must be a path to directory"
    assert args_yml_file.is_file(), "args_yml_file must be a path to a file"
    assert dbs_dir.is_dir(), "dbs_dir must be a path to directory"
    # Preparing the logger
    prep_log = proc_helper.gimmelogger("preparion_logger", fastq_file_dir.joinpath("preparation_logs.txt"))
    # TODO
    gather_files_in_pairs()
    split_files_to_directories()
    run_preprocessing_parallel()
    if not only_preproc:
        run_analysis()


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
                        default="IMNGS2Pipeline_args.yml")
    parser.add_argument("-d", "--db-directory",
                        type=str,
                        help="Path to directory containing silva, sortmerna files",
                        default="./database")
    parser.add_argument("-op", "--only-preprocess",
                        action="store_true",
                        help="Should run analysis step")
    args = parser.parse_args()
    run_imngs2(
        args.fastq_directory,
        args.yml_file,
        args.db_directory,
        args.only_preprocess
    )
