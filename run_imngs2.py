import re
import argparse
import shutil
from imngs2_pipeline.processing_job import main_processing as preprocessing
from imngs2_pipeline.processing_job import calc_spikes, gzip_to_fastq
from imngs2_pipeline.analyze_sample_job import main as main_analysis
from collections import namedtuple
import imngs2_pipeline.processing_helper as proc_helper
from multiprocessing import cpu_count, Pool
from typing import List, Dict
from sys import version_info
if version_info[0] < 3:
    from pathlib2 import Path, PurePath  # pip2 install pathlib2
else:
    from pathlib import Path, PurePath


global INPUT_DIR, DBS_DIR, POOL_SIZE, PREP_LOG, MAPPING_FILE_COLS, MapLineTup
INPUT_DIR = "/srv/base/inputs/"
DBS_DIR = "/srv/base/databases/"
max_pool = int(cpu_count() * 0.3)
POOL_SIZE = max_pool if max_pool > 0 else 1
max_batch = int(POOL_SIZE * 0.8)
# Preparing the logger
PREP_LOG = proc_helper.gimmelogger(
    "run_imngs2",
    log_file=Path(PurePath(INPUT_DIR)).absolute().joinpath("Pipeline_log.txt"),
    only_file=False)
MAPPING_FILE_COLS = ("SampleID", "total_weight_in_g", "amount_spike")
MapLineTup = namedtuple("MapLineTup", [MAPPING_FILE_COLS[0], MAPPING_FILE_COLS[1], MAPPING_FILE_COLS[2]])
global SPIKE_STAT_FILE_NAME, SPIKE_STAT_HEADER
SPIKE_STAT_FILE_NAME = "spike_stat_mapping_file.tsv"
row_fmt = "{}\t{}\t{}\t{}"
SPIKE_STAT_HEADER = row_fmt.format("#SampleID", "SpikeReads", "spikes_total_weight_in_g", "amount_spike")


def parse_mapping_file(mapping_file_path: str, files_tups: list = []):
    """
    files_tups: should be a list of files path tuples. i.e: [(<forward_path>, <reverse_path>), ...]
    """
    mapping_file_path = Path(PurePath(mapping_file_path)).absolute()
    valid_ids = []
    ignored_ids = []  # e.g. entries with NA values
    mapping_lines: Dict[str, List[str]] = {}  # from mapping file, full line for each relevant sample
    index_of_sample_id_col = None
    sample_col_name = MAPPING_FILE_COLS[0]
    index_of_weight_col = None
    weight_col_name = MAPPING_FILE_COLS[1]
    index_of_amounts_col = None
    amounts_col_name = MAPPING_FILE_COLS[2]
    # Preparing default argument for each sample preprocessing
    for file_pair in files_tups:
        sample_id = proc_helper.get_base_name(file_pair[0])
        _row_vals_tup = MapLineTup(sample_id, float("NAN"), 0)
        mapping_lines[sample_id] = (_row_vals_tup, file_pair)
    # If mapping_file exists then we parse it and change the default of sample_weight, spike_mount to actual values.
    if not mapping_file_path.is_file():
        PREP_LOG.warning("Parsing of mapping file failed. Continuing as if there is no spike in the samples.")
        return mapping_lines

    PREP_LOG.warning(f"When mapping file is provided ({mapping_file_path}), only samples in mapping file will get processed!!!")
    # If mapping_file_is there then renew the mapping_lines
    mapping_lines: Dict[str, List[str]] = {}  # from mapping file, full line for each relevant sample

    with open(mapping_file_path, 'r') as mapping_file_h:
        header: str = next(mapping_file_h)

        if not header.startswith('#'):
            PREP_LOG.error('No header in mapping file')
            raise Exception('No header in mapping file')

        columns = header.split('\t')
        for i, col in enumerate(columns):
            colname = col.strip().lstrip('#')

            if colname == sample_col_name:
                index_of_sample_id_col = i
            elif colname == weight_col_name:
                index_of_weight_col = i
            elif colname == amounts_col_name:
                index_of_amounts_col = i

        if index_of_sample_id_col is None:
            msg = f'Column {sample_col_name!r} missing in header of mapping_file: {mapping_file_path}'
            PREP_LOG.error(msg)
            raise Exception(msg)

        if index_of_weight_col is None:
            msg = f'Column {weight_col_name!r} missing in header of mapping_file: {mapping_file_path}'
            PREP_LOG.error(msg)
            raise Exception(msg)

        if index_of_amounts_col is None:
            msg = f'Column {amounts_col_name!r} missing in header of mapping_file: {mapping_file_path}'
            PREP_LOG.error(msg)
            raise Exception(msg)

        for line in mapping_file_h:
            line = line.strip()
            if line.startswith('#'):
                continue

            fields = line.split('\t')
            sample_id = fields[index_of_sample_id_col].strip()
            total_weight_in_g = fields[index_of_weight_col].strip()
            amount = fields[index_of_amounts_col].strip()

            # warn if weight is not a valid number
            try:
                float(total_weight_in_g)
            except ValueError:
                PREP_LOG.warning(
                    f'Weight {total_weight_in_g!r} is not a valid floating point number for {sample_id!r}. Weigth will be processed as NAN.')
                fields[index_of_weight_col] = "NAN"

            # if amount cannot be parsed to float change to 0
            try:
                float(amount)
            except ValueError:
                fields[index_of_amounts_col] = "0"

            valid_ids.append(sample_id)
            map_line = MapLineTup(fields[index_of_sample_id_col], fields[index_of_weight_col], fields[index_of_amounts_col])
            # Mapping file paris to sample_id
            file_pair_tup = ("", "")
            for file_tup in files_tups:
                try:
                    forw_name = file_tup[0]
                except Exception:
                    forw_name = ""
                if sample_id in forw_name:
                    file_pair_tup = file_tup
                    # NOTE We only add valid file names which are given in mapping file
            mapping_lines[sample_id] = (map_line, file_pair_tup)

    return mapping_lines


def remove_spikes(
        files_paths: list,
        sample_id: str,
        sample_weight: float,  # in gram
        spike_amount: int):
    """
    It remove spikes from a sample and adds the number of spikes to a mapping file given.
     the column name should be "SpikeReads". The format of spike_stat_mapping_path must be
     like below:
     #SampleID\tSpikeReads\tspikes_total_weight_in_g\tamount_spike\n'

    Arguments:
        - files_path: list = list of files on which spike removal should be applied
        - sample_id: str = sample identifier to be reported in spike_stat_mapping_path
        - sample_weight: float = sample weight sent for sequencing in gram
        - spike_amount: int = amount of spike in nanogram

    """
    row_fmt = "{}\t{}\t{}\t{}"
    header_rgx = re.compile(SPIKE_STAT_HEADER)
    files_paths_abs = [Path(PurePath(fi_pa)) for fi_pa in files_paths if Path(PurePath(fi_pa)).is_file()]
    spike_stat_mapping_path = files_paths_abs[0].parent.joinpath(SPIKE_STAT_FILE_NAME)
    # NOTE calc_spikes() only works with fastq files. After this step the path to files would change
    # It's important to continue downstream analysis with new files.
    fastq_files_abs_path = gzip_to_fastq(*files_paths_abs)
    # Preparing default output
    out_tup = (0, spike_stat_mapping_path, tuple(fastq_files_abs_path))
    if len(files_paths) != len(fastq_files_abs_path):
        msg = "Some given files for spike removal do not exist!"
        PREP_LOG.error(msg)
        raise ValueError(msg)
    # Writing header
    header_existing = SPIKE_STAT_HEADER + "\n"
    with open(spike_stat_mapping_path, 'w') as stats_h:
        stats_h.write(header_existing)
    # Adding rows to spike stat file.
    header_line_mo = header_rgx.match(header_existing)
    if not header_line_mo:
        msg = f"Existing spike_stat_file: {spike_stat_mapping_path} does not have correct format of : {SPIKE_STAT_HEADER}"
        PREP_LOG.error(msg)
        raise ValueError(msg)

    # postponing file openning to avoid race condition if ran parallel
    real_reads_c, spike_reads_c = calc_spikes(*fastq_files_abs_path, spike_amount=spike_amount)
    out_tup = (spike_reads_c, spike_stat_mapping_path, tuple(fastq_files_abs_path))
    with open(spike_stat_mapping_path, 'a') as stats_h:
        # appending
        stats_h.write(row_fmt.format(sample_id, spike_reads_c, sample_weight, spike_amount) + "\n")

    PREP_LOG.info(f"Spike removal done for {sample_id}.Spike Reads: {spike_reads_c}\tnon-spike reads: {real_reads_c}")

    return out_tup


def run_preprocessing(
        seq_files_t: tuple,
        preproc_dir: str,
        args_yml_path: str,
        sample_id: str = "",  # If it's not provided then the base name of forward file
        sample_weight: float = float("NAN"),  # spike normalizer handles this
        spike_amount: int = 0):
    """
    It takes a tuple of paths to sequencing files.
    Create directory for basename of files and move
     them there. Then triggers the preprocessing step.
     It returns the path of direcotry in which the
     preprocessd files of sample are.
    """
    sample_dir = ""
    try:
        if not hasattr(seq_files_t, "__iter__") or len(seq_files_t) > 2:
            msg = "seq_files_t must be an iterable of max lenght 2"
            # PREP_LOG.error(msg)
            raise ValueError(msg)
        # Creating the preprocessing directory
        seq_files_t = [Path(PurePath(sfi)).absolute() for sfi in seq_files_t]
        new_paths = ["", ""]  # [ForwardNewPath, ReverseNewPath]
        new_dir = seq_files_t[0]
        # Getting initial stem
        while new_dir.suffixes:
            new_dir = Path(new_dir).absolute().parent.joinpath(new_dir.stem)

        sample_base_name = proc_helper.get_base_name(new_dir)
        sample_id = sample_id if sample_id else sample_base_name
        new_dir = Path(PurePath(preproc_dir)).absolute().joinpath(sample_id)
        new_dir.mkdir(parents=True, exist_ok=True)
        for ind, sfi in enumerate(seq_files_t):
            new_path = new_dir.joinpath(sfi.name)
            new_paths[ind] = new_path
            # Copying files
            shutil.copy2(str(sfi), str(new_path))
        # TODO we do spike removal if necessary and add a line to spike_stats file for spike normalization
        new_paths = [el for el in new_paths if bool(el)]
        # removing spikes and decompressing files below
        spike_reads_c, spike_stat_mapping_path, fastq_files_tuple = remove_spikes(
            new_paths,
            sample_id,
            sample_weight,
            spike_amount
        )
        # Running preprocessing
        sample_dir = preprocessing(
            input_dir=new_dir,
            paired="Yes" if len(fastq_files_tuple) == 2 else "No",
            forward_file=fastq_files_tuple[0],
            reverse_file=fastq_files_tuple[1],
            input_id=sample_id,
            args_file_path=args_yml_path,
            spike_amount=0,  # We run it always with 0 as we remove spikes before if there is
        )
    except Exception as exc:
        msg = f"{exc}"
        PREP_LOG.error(msg)
        shutil.rmtree(str(new_dir))
        PREP_LOG.warning("Deleting {}".format(new_dir))

    return sample_dir


def combine_spike_stats_file(*samples_dirs, combined_spike_stat: str = "."):
    """
    It combines samples' spike_stat_file to one file to be input to spike normalization step
    """
    samples_dirs_path = []
    for sam_dir in samples_dirs:
        try:
            sam_dir_path = Path(PurePath(sam_dir)).absolute()
            if not sam_dir_path.is_dir():
                raise ValueError(f"{sam_dir_path} does not exist!")
        except Exception as exc:
            PREP_LOG.warning(f"Invalid sample directory for {sam_dir}. {exc}")
            continue
        samples_dirs_path.append(sam_dir_path)

    combined_spike_stat = Path(PurePath(combined_spike_stat))
    if combined_spike_stat.is_dir():
        combined_spike_stat = samples_dirs_path[0].parent.joinpath(SPIKE_STAT_FILE_NAME)

    sample_stat_files_paths = [sdir.joinpath(SPIKE_STAT_FILE_NAME) for sdir in samples_dirs_path]
    header_rgx = re.compile(SPIKE_STAT_HEADER)
    with open(combined_spike_stat, "w") as combined_stat_file_fio:
        combined_stat_file_fio.write(SPIKE_STAT_HEADER + "\n")
        for ssfp in sample_stat_files_paths:
            with open(ssfp) as sam_stat_file_fio:
                header = sam_stat_file_fio.readline().strip()
                row_values = sam_stat_file_fio.readline().strip()
            header_line_mo = header_rgx.match(header)
            if not header_line_mo:
                msg = f"Existing spike_stat_file: {ssfp} does not have correct format of : {SPIKE_STAT_HEADER}"
                PREP_LOG.error(msg)
                continue
            combined_stat_file_fio.write(row_values + "\n")

    samples_dirs_path = [str(el) for el in samples_dirs_path]
    return samples_dirs_path, combined_spike_stat


def gather_files(*preproc_dirs_path, file_name="taxed_ZOTUs.fasta"):
    # print(preproc_dirs_path)
    sam_dirs = []
    preproc_dirs_path_deepests = []
    for prep_dir in preproc_dirs_path:
        prep_dir = Path(PurePath(prep_dir)).absolute()
        if not prep_dir.is_dir():
            continue
        preproc_dirs_path_deepests.append(str(prep_dir))
    for prep_dir in preproc_dirs_path_deepests:
        prep_dir = Path(PurePath(prep_dir))
        for diri in range(len(preproc_dirs_path_deepests)):
            added_prep_dir = Path(PurePath(preproc_dirs_path_deepests[diri]))
            deeper_path = added_prep_dir if len(added_prep_dir.parts) > len(prep_dir.parts) else prep_dir
            shallower_path = prep_dir if deeper_path == added_prep_dir else added_prep_dir
            if str(shallower_path) in str(deeper_path):
                preproc_dirs_path_deepests[diri] = str(deeper_path)

    preproc_dirs_path = list(set(preproc_dirs_path_deepests))
    size_rgx = re.compile(r"size=([0-9]+);")
    tax_rgx = re.compile(r"tax=([^;]+){0,7};")
    for ds in preproc_dirs_path:
        for fs in Path(PurePath(ds)).rglob(file_name):
            with open(fs) as fs_fio:
                line = fs_fio.readline().strip()
                size_mo = size_rgx.search(line)
                tax_mo = tax_rgx.search(line)
                if not all([line, line.startswith(">"), size_mo, tax_mo]):
                    PREP_LOG.warning(f"{str(fs)} file is not formatted correctly. Proper fasta header format: '>SEQID;size=XXX;tax=KKK,PPP;CCC;OOO")
                    continue
            sam_dirs.append(str(fs))

    return sam_dirs


def run_imngs2(
        fastq_file_dir: str,
        args_yml_file: str,
        dbs_dir: str,
        mapping_file: str = "",
        spike_stat_file: str = "",
        skip_preprocess: bool = False,
        skip_analysis: bool = False):
    fastq_file_dir = Path(PurePath(fastq_file_dir)).absolute()
    args_yml_file = Path(PurePath(args_yml_file)).absolute()
    dbs_dir = Path(PurePath(dbs_dir))
    assert fastq_file_dir.is_dir(), "fastq_file_dir must be a path to directory"
    assert args_yml_file.is_file(), "args_yml_file must be a path to a file"
    assert dbs_dir.is_dir(), "dbs_dir must be a path to directory"
    preproc_dir = fastq_file_dir.joinpath("Preprocessing")
    preproc_dir.mkdir(parents=True, exist_ok=True)
    default_spike_stat_compiled = preproc_dir.joinpath(SPIKE_STAT_FILE_NAME)
    combined_spike_stats_path = default_spike_stat_compiled if not spike_stat_file else Path(PurePath(spike_stat_file))
    reduced_samples_dirs = []
    mapping_file_path = Path(PurePath(str(mapping_file))).absolute() if mapping_file else Path.cwd()
    # Analysis defualt vars
    zotu_file_path, sotu_file_path = ("", "")
    # assert mapping_file_path.is_file(), f"{mapping_file_path} must be a path to an exisitng mapping file"
    if not skip_preprocess:
        try:
            # Gathering sequence files
            PREP_LOG.debug("Gathering sequence files in {}".format(str(fastq_file_dir)))
            seq_file_pairs = proc_helper.pair_seq_files(str(fastq_file_dir))
            PREP_LOG.info("{} sampels were collected from {}.".format(str(len(seq_file_pairs)), str(fastq_file_dir)))
            # If mapping_file exists then we parse it and change the default of sample_weight, spike_mount to actual values.
            mapping_line_tup_dict = parse_mapping_file(mapping_file_path, seq_file_pairs)
            # running preprocessing
            with Pool(POOL_SIZE, maxtasksperchild=1) as pool:
                res_list = [pool.apply_async(run_preprocessing, args=((arg_tup[1]), preproc_dir, str(args_yml_file), arg_tup[0].SampleID, arg_tup[0].total_weight_in_g, arg_tup[0].amount_spike,))
                            for sample_id, arg_tup in mapping_line_tup_dict.items()]
                for res in res_list:
                    res.wait()
            samples_dirs = [el.get() for el in res_list]

        except Exception as exc:
            PREP_LOG.error(f"Preprocessing Failed: {exc}")
            # shutil.rmtree(str(preproc_dir))  # TODO decide if we need to delete Preprocessing direcotory in case of failure.
            PREP_LOG.warning("Deleting {}".format(preproc_dir))
            skip_analysis = True
            samples_dirs = []
    else:
        PREP_LOG.warning(f"Skipping Preprocessing. Processing sammples in directory {preproc_dir}")
        samp_zotu_seq_files = gather_files(*[preproc_dir], file_name="taxed_ZOTUs.fasta")
        samples_dirs = [Path(PurePath(el)).parent for el in samp_zotu_seq_files]

    # Runing analysis
    if not skip_analysis:
        try:
            analysis_dir = fastq_file_dir.joinpath("Analysis")
            analysis_dir.mkdir(parents=True, exist_ok=True)
            # we do the step down because if skip_preprocess is True then tha path to preprocess would be given not all smaple_dirs
            if str(combined_spike_stats_path) != str(default_spike_stat_compiled):
                PREP_LOG.warning(f"Using custom spike_stat file: {combined_spike_stats_path} for spike normalization!! Default spike_stat file {default_spike_stat_compiled} is ignored!")
            reduced_samples_dirs, _ = combine_spike_stats_file(*samples_dirs, combined_spike_stat=combined_spike_stats_path)
            missed_samples = [sam_dir for sam_dir in samples_dirs if str(sam_dir) not in reduced_samples_dirs]
            if missed_samples:
                PREP_LOG.warning(f"Samples in {str(combined_spike_stats_path)} will be sent for analysis. Please Check the file.")
                PREP_LOG.warning(f"Missed Samples are: {missed_samples}")
            samp_zotu_seq_files = gather_files(*reduced_samples_dirs, file_name="taxed_ZOTUs.fasta")
            zotu_file_path, sotu_file_path = main_analysis(samp_zotu_seq_files, analysis_dir, args_yml_file, combined_spike_stats_path, dbs_loc=dbs_dir)
        except Exception as exc:
            PREP_LOG.error(f"Analysis Failed: {exc}")

        ##################
    # TODO Only Normalizing
    if not skip_analysis and mapping_file_path.is_file() and zotu_file_path and sotu_file_path:
        pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    help_text = "The directory to find recursively all fastq files inside"
    parser.add_argument("-d", "--fastq-directory",
                        type=str,
                        help=help_text,
                        default=INPUT_DIR)  # WORKDIR of container is /srv/base/inputs
    parser.add_argument("-y", "--yml-file",
                        type=str,
                        help="Path to arguments yaml file",
                        default=Path(PurePath(INPUT_DIR)).joinpath("IMNGS2Pipeline_args.yml"))
    parser.add_argument("-map", "--mapping-file",
                        type=str,
                        default=Path(PurePath(INPUT_DIR)).joinpath("mapping_file.tsv"),
                        help="The path to a mapping file defining sample weight and spike amount for each sample")
    parser.add_argument("-stat", "--spike-stat",
                        type=str,
                        default="",
                        help=f"The path to a mapping file defining spike count, sample weight and spike amount for each sample.\n{SPIKE_STAT_HEADER}")
    parser.add_argument("-db", "--db-directory",
                        type=str,
                        help="Path to directory containing silva, sortmerna files",
                        default=DBS_DIR)
    parser.add_argument("-sp", "--skip-preprocess",
                        action="store_true",
                        help="Should skip preprocessing step")
    parser.add_argument("-sa", "--skip-analysis",
                        action="store_true",
                        help="Should skip analysis step")
    args = parser.parse_args()
    if args.db_directory != DBS_DIR:
        DBS_DIR = args.db_directory
    run_imngs2(
        args.fastq_directory,
        args.yml_file,
        args.db_directory,
        args.mapping_file,
        args.spike_stat,
        args.skip_preprocess,
        args.skip_analysis,
    )
