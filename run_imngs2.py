#! /usr/local/bin/python
import re
import argparse
import shutil
from os import symlink, chdir
from imngs2_pipeline.processing_job import main_processing as preprocessing
from imngs2_pipeline.processing_job import calc_spikes, gzip_to_fastq
from imngs2_pipeline.analyze_sample_job import main as main_analysis
from collections import namedtuple
import imngs2_pipeline.processing_helper as proc_helper
from multiprocessing import cpu_count, Pool
from typing import List
from sys import version_info
if version_info[0] < 3:
    from pathlib2 import Path, PurePath, PureWindowsPath, PurePosixPath  # pip2 install pathlib2
else:
    from pathlib import Path, PurePath, PureWindowsPath, PurePosixPath


global INPUT_DIR, DBS_DIR, POOL_SIZE, PREP_LOG, MAPPING_FILE_COLS
global FASTQ_DIR, MapLineTup, ARGS_YAML_FILE, PROC_DIR_SUFFIX
INPUT_DIR = Path(PurePath("/base/inputs/")).absolute()
FASTQ_DIR = INPUT_DIR
PROC_DIR_SUFFIX = "__processed"
global DEFAULT_ARG_FILE_NAME
DEFAULT_ARG_FILE_NAME = "IMNGS2Pipeline_args.yml"
ARGS_YAML_FILE = Path(PurePath("/base/" + DEFAULT_ARG_FILE_NAME)).absolute()
MAP_FILE = Path(PurePath("/base/mapping_file_TEMPLATE.tsv")).absolute()
DBS_DIR = "/base/databases/"
max_pool = int(cpu_count() * 0.3)
POOL_SIZE = max_pool if max_pool > 0 else 1
max_batch = int(POOL_SIZE * 0.8)
# Preparing the logger
PREP_LOG = proc_helper.gimmelogger(
    "run_imngs2",
    log_file=INPUT_DIR.joinpath("Pipeline_log.txt"),
    only_file=False)
MAPPING_FILE_COLS = ("SampleID", "total_weight_in_g", "spike_amount", "parent_path")
MapLineTup = namedtuple("MapLineTup", [MAPPING_FILE_COLS[0], MAPPING_FILE_COLS[1], MAPPING_FILE_COLS[2], MAPPING_FILE_COLS[3]])
global SPIKE_STAT_FILE_NAME, SPIKE_STAT_HEADER, PATH_SEP, DEFAULT_MAP_LINE
DEFAULT_MAP_LINE = MapLineTup("", float("NAN"), "0", "")
PATH_SEP = "_-_"
SPIKE_STAT_FILE_NAME = "spike_stat_mapping_file.tsv"
row_fmt = "{}\t{}\t{}\t{}"
SPIKE_STAT_HEADER = row_fmt.format("#SampleID", "SpikeReads", "spikes_total_weight_in_g", "spike_amount")


def is_in_processed_dir(fastq_path: Path, sample_id: str = None) -> bool:
    default = None
    type_test = [
        isinstance(fastq_path, str),
        isinstance(fastq_path, Path),
    ]
    if not any(type_test):
        return default
    default = False
    base_name = proc_helper.get_base_name(str(fastq_path)) if sample_id is None else sample_id
    if f"{str(base_name)}{PROC_DIR_SUFFIX}" in str(fastq_path):
        default = True

    return default


def is_argset_different(existing_arg_file: Path, processed_arg_file: Path = None) -> bool:
    """
    It decides if two set of arguments are the same or not for deciding whether a sample
    should be passed for processing or not.
    """
    is_different = True
    tests = [
        isinstance(existing_arg_file, str) or isinstance(existing_arg_file, Path),
        isinstance(processed_arg_file, Path) or isinstance(processed_arg_file, str) or processed_arg_file is None,
    ]
    if not all(tests):
        raise ValueError("Arguments should be str or Path object")

    existing_arg_file = Path(PurePath(existing_arg_file)).absolute()
    proce_arg_file = Path(PurePath(processed_arg_file)).absolute() if processed_arg_file else None
    if proce_arg_file:
        given_argset_obj = proc_helper.IMNGS2ArgsParser(config_yaml=existing_arg_file)
        proce_argset_obj = proc_helper.IMNGS2ArgsParser(config_yaml=proce_arg_file)

        is_different = not given_argset_obj.preproc_args.__eq__(proce_argset_obj.preproc_args)

    return is_different


def get_processed_argset(processed_dir_path: Path) -> List[Path]:
    out_list = []
    if bool(isinstance(processed_dir_path, str) or isinstance(processed_dir_path, Path)):
        processed_dir_path = Path(PurePath(processed_dir_path)).absolute()
        if not processed_dir_path.is_dir() or not str(processed_dir_path).endswith(PROC_DIR_SUFFIX):
            raise TypeError(f"{processed_dir_path} must be a path to a directory ending with '{PROC_DIR_SUFFIX}'")
        # finding argument files
        out_list = processed_dir_path.rglob(DEFAULT_ARG_FILE_NAME)

    return out_list


def is_processed_dir_healthy(processed_dir_date_version_path: Path) -> bool:
    my_dir = processed_dir_date_version_path
    is_valid_path = [
        isinstance(my_dir, str),
        isinstance(my_dir, Path)
    ]
    if not my_dir.is_dir() or not any(is_valid_path):
        return False
    my_dir = Path(PurePath(my_dir)).absolute()
    should_be_there = [
        DEFAULT_ARG_FILE_NAME,
        "taxed_ZOTUs.fasta",
        SPIKE_STAT_FILE_NAME
    ]
    logic_test = [False for fi in should_be_there]

    for ind in range(len(should_be_there)):
        dir_content = my_dir.iterdir()
        for cont in dir_content:
            if cont.is_file() and str(cont.name) == str(should_be_there[ind]):
                logic_test[ind] = True

    return all(logic_test)


def should_trigger_processing(sample_base_path: Path, given_argset_file: Path) -> (Path, bool):
    this_out = ["", True]
    given_argset_file = Path(PurePath(given_argset_file)).absolute() if isinstance(given_argset_file, str) or isinstance(given_argset_file, Path) else ""
    sample_base_path = Path(PurePath(sample_base_path)).absolute() if isinstance(sample_base_path, str) or isinstance(sample_base_path, Path) else ""
    if not given_argset_file or not sample_base_path:
        raise ValueError(f"Not Valid sample_base_path or not valid argument file. Both should be string or Path object.")

    processed_dir_path = Path(PurePath(str(sample_base_path))).absolute()
    base_name = processed_dir_path.name.replace(PROC_DIR_SUFFIX, "")

    if not str(sample_base_path).endswith(PROC_DIR_SUFFIX):
        processed_dir_path = Path(PurePath(str(sample_base_path) + PROC_DIR_SUFFIX))

    this_out[0] = str(processed_dir_path.joinpath(proc_helper.generate_timestamp()))
    arg_files = processed_dir_path.rglob(DEFAULT_ARG_FILE_NAME)

    for argfile in arg_files:
        if not argfile.is_file():
            continue
        is_already_processed = is_processed_dir_healthy(argfile.parent) and not is_argset_different(argfile, given_argset_file)
        if is_already_processed:
            this_out[0] = str(argfile.parent.absolute())
            this_out[1] = False

    if not this_out[1]:
        PREP_LOG.info(f"Sample '{base_name}' is already processed with given arguemnt set. "
                      f"Skipping processing and using results in {this_out[0]}")

    return tuple(this_out)


def valid_for_analysis(dir_path: Path, given_arg_file: Path) -> bool:
    dir_path = Path(PurePath(dir_path)).absolute() if isinstance(dir_path, str) or isinstance(dir_path, Path) else ""
    given_arg_file = Path(PurePath(given_arg_file)).absolute() if isinstance(given_arg_file, str) and isinstance(given_arg_file, Path) else ""
    if not dir_path and not given_arg_file:
        return False
    is_complete = is_processed_dir_healthy(dir_path)
    old_arg_file = dir_path.joinpath(DEFAULT_ARG_FILE_NAME)
    same_argset = not is_argset_different(given_arg_file, old_arg_file)

    return bool(is_complete and same_argset)


def uniqify_map_lines(dup_ids_maplinetup_list: list) -> list:  # List[(MapLineTup, tuple)]
    sample_ids_count = proc_helper.MyCounter([x[0].SampleID for x in dup_ids_maplinetup_list])
    most_common = proc_helper.MyCounter(dict(sample_ids_count.most_common(1)))
    while most_common.total() > 1:
        most_common_sample_id, most_common_count = most_common.items()
        for ind in range(len(dup_ids_maplinetup_list)):
            map_line_obj = dup_ids_maplinetup_list[ind][0]
            file_pair = dup_ids_maplinetup_list[ind][1]
            if most_common[map_line_obj.SampleID]:
                new_map_line = MapLineTup(
                    str(map_line_obj[3]) + PATH_SEP + str(most_common_sample_id) + "__" + str(most_common_count),
                    map_line_obj[1],
                    map_line_obj[2],
                    map_line_obj[3],
                )
                dup_ids_maplinetup_list[ind] = (new_map_line, file_pair)
                # End Phase of loop
                dup_ids_maplinetup_list = uniqify_map_lines(dup_ids_maplinetup_list)
                sample_ids_count = proc_helper.MyCounter([x[0][0] for x in dup_ids_maplinetup_list])
                most_common = proc_helper.MyCounter(dict(sample_ids_count.most_common(1)))

    return dup_ids_maplinetup_list


def convert_mapping_entries_to_dict(map_line_entries_list: list) -> dict:  # List[(MapLineTup, tuple)]
    # Sorting by sample_id
    # map_line_entries_list = sorted(map_line_entries_list, key= lambda x: x[0][0])
    init_count = len(map_line_entries_list)
    sample_ids_count = proc_helper.MyCounter([x[0][0] for x in map_line_entries_list])
    most_common = proc_helper.MyCounter(dict(sample_ids_count.most_common(1)))

    while most_common.total() > 1:
        map_line_entries_list_filtered = []
        # deleting duplicates from map_line_entries
        sub_list_to_uniqify = []
        for ind in range(len(map_line_entries_list)):
            mapline_tup = map_line_entries_list[ind]
            if most_common[mapline_tup[0].SampleID]:
                sub_list_to_uniqify.append(mapline_tup)
            else:
                map_line_entries_list_filtered.append(map_line_entries_list[ind])
        map_line_entries_list_filtered.extend(uniqify_map_lines(sub_list_to_uniqify))
        map_line_entries_list = map_line_entries_list_filtered
        sample_ids_count = proc_helper.MyCounter([x[0][0] for x in map_line_entries_list])
        most_common = proc_helper.MyCounter(dict(sample_ids_count.most_common(1)))

    mapping_line_dict = {element[0].SampleID: element for element in map_line_entries_list}
    assert (len(mapping_line_dict) == init_count), "Different length for map_line_entries list and output dictionary"
    return mapping_line_dict


def correct_path_for_windows(input_path: str) -> Path:
    logic_test_list = [
        isinstance(input_path, str)
    ]
    if not any(logic_test_list):
        raise TypeError("Not a valid path")
    input_path = str(input_path)
    if "\\\\" in input_path or "\\" in input_path:
        return PureWindowsPath(input_path)
    elif "/" in input_path:
        return PurePosixPath(input_path)

    return PurePath(FASTQ_DIR)


def is_correct_parent(
        fastq_file: Path,
        given_parent_path: str) -> bool:
    """
    Choosing a point in hierarchy of file system, any child path is unique including a
    unique character or name differntiating from other child path. The unique part is not
    in parent path as the selected point have only one parent.
    """
    base_path = FASTQ_DIR
    fastq_file = Path(PurePath(fastq_file))
    pre_checks_bools = [
        not isinstance(given_parent_path, str),
        not isinstance(given_parent_path, Path),
        not isinstance(fastq_file, Path),
    ]
    if not any(pre_checks_bools) and not fastq_file.is_file():
        return False
    given_parent_path = str(FASTQ_DIR) if given_parent_path == "" else given_parent_path
    given_parent_path = correct_path_for_windows(str(given_parent_path))
    fastq_file_parent = fastq_file.parent.relative_to(base_path)
    # We assume that the user gives the full path to parent directory in worst case.
    given_parent_parts = given_parent_path.parts
    ind = -1
    ind_limit = len(given_parent_parts) * -1
    # If given fastq_file_parent is relative and exhausted
    # if fastq_file_parent is given with leading '/'
    # if fastq_file is set in c://
    # if parts of given_parent_parts is exhausted
    while True:
        while_terminator = [
            str(fastq_file_parent) != ".",
            str(fastq_file_parent) != "/",
            not str(fastq_file_parent).endswith(":/"),
            ind >= ind_limit,
        ]

        if not all(while_terminator):
            break
        this_level = fastq_file_parent.name
        if this_level != given_parent_parts[ind]:
            return False
        # Updating while loop parameters
        fastq_file_parent = fastq_file_parent.parent
        ind -= 1

    if str(fastq_file_parent) != ".":
        # The given path is only partially covering the fastq path
        # This is ambiguous
        return False

    return True


def parse_mapping_file(mapping_file_path: str, files_tups: list = []):
    """
    files_tups: should be a list of files path tuples. i.e: [(<forward_path>, <reverse_path>), ...]
    """
    mapping_file_path = Path(PurePath(mapping_file_path)).absolute()
    valid_ids = []
    mapping_lines = []  # List[(MapLineTup, tuple)] -> from mapping file, full line for each relevant sample
    index_of_sample_id_col = None
    sample_col_name = MAPPING_FILE_COLS[0]
    index_of_weight_col = None
    weight_col_name = MAPPING_FILE_COLS[1]
    index_of_amounts_col = None
    amounts_col_name = MAPPING_FILE_COLS[2]
    index_of_parent_path_col = None
    parent_path_col_name = MAPPING_FILE_COLS[3]
    # Preparing default argument for each sample preprocessing

    for file_pair in files_tups:
        sample_id = proc_helper.get_base_name(file_pair[0])
        _row_vals_tup = MapLineTup(
            sample_id,
            DEFAULT_MAP_LINE.total_weight_in_g,
            DEFAULT_MAP_LINE.spike_amount,
            DEFAULT_MAP_LINE.parent_path
        )
        mapping_lines.append((_row_vals_tup, file_pair))
    # If mapping_file exists then we parse it and change the default of sample_weight, spike_mount to actual values.
    mapping_lines = convert_mapping_entries_to_dict(mapping_lines)
    if not mapping_file_path.is_file():
        PREP_LOG.warning("Parsing of mapping file failed. Continuing with fake mapping file.")
        return mapping_lines

    PREP_LOG.warning(f"When mapping file is provided ({mapping_file_path}), only samples in mapping file will get processed!!!")
    # If mapping_file_is there then renew the mapping_lines
    mapping_lines = []  # List[(MapLineTup, file_pair)] -> from mapping file, full line for each relevant sample

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
            elif colname == parent_path_col_name:
                index_of_parent_path_col = i

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

        # injecting parent_path
        index_of_parent_path_col = index_of_parent_path_col if index_of_parent_path_col else len(MAPPING_FILE_COLS)
        # looping over lines
        for line in mapping_file_h:
            line = line.strip().rstrip()
            if line.startswith('#'):
                continue
            # placing default values
            hypo_fields = ["", "NAN", "0", ""]
            # Parsing lines
            fields = line.split('\t')
            try:
                sample_id = fields[index_of_sample_id_col].strip()
                hypo_fields[0] = sample_id
            except Exception:
                raise ValueError(f"Could not parse SampleID from mapping file. Error on line: {line}")
            try:
                total_weight_in_g = fields[index_of_weight_col].strip()
                hypo_fields[1] = total_weight_in_g
            except Exception:
                PREP_LOG.warning(f"Could not parse weight amount from mapping file. Error on line: {line}\n"
                                 f" Check the mapping file format. columns should be"
                                 f" separated by TAB. Continuing wiht default value of NAN")
            try:
                amount = fields[index_of_amounts_col].strip()
                hypo_fields[2] = amount
            except Exception:
                PREP_LOG.warning(f"Could not parse spike amount from mapping file. Error on line: {line}\n"
                                 f" Check the mapping file format. columns should be"
                                 f" separated by TAB. Continuing wiht default value of 0")
            try:
                parent_path = fields[index_of_parent_path_col]
                parent_path = parent_path.strip().rstrip("\\").rstrip("/")
                hypo_fields[3] = parent_path
            except Exception:
                PREP_LOG.debug(f"Could not parse parent path from mapping file. Error on line: {line}\n"
                               f" Check the mapping file format. columns should be"
                               f" separated by TAB. Continuing wiht default value of 0")
            # warn if weight is not a valid number
            try:
                float(hypo_fields[1])
            except ValueError:
                PREP_LOG.warning(
                    f'Weight {total_weight_in_g!r} is not a valid floating point number for {sample_id!r}. Weigth will be processed as NAN.')
                hypo_fields[1] = "NAN"

            # if amount cannot be parsed to float change to 0
            try:
                float(hypo_fields[2])
            except ValueError:
                hypo_fields[2] = "0"

            valid_ids.append(sample_id)
            map_line = MapLineTup(
                hypo_fields[0],
                hypo_fields[1],
                hypo_fields[2],
                hypo_fields[3],
            )
            # Mapping file paris to sample_id
            file_pairs_tup = ("", "")
            # files_tups contain all found fastq files in fastq directory
            for ind in range(len(files_tups)):
                file_tup = files_tups[ind]
                try:
                    forw_file_name = str(file_tup[0]).split("/")[-1]
                except Exception:
                    forw_file_name = ""
                correct_parent = is_correct_parent(file_tup[0], map_line.parent_path)
                if sample_id in forw_file_name and correct_parent:  # sample_id could also be partial path of file_names
                    file_pairs_tup = file_tup
                    del files_tups[ind]
                    break

            mapping_lines.append((map_line, file_pairs_tup))
    mapping_lines = convert_mapping_entries_to_dict(mapping_lines)

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
     #SampleID\tSpikeReads\tspikes_total_weight_in_g\tspike_amount\n'

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

    PREP_LOG.info(f"Spike removal done for {sample_id}. Spike Reads: {spike_reads_c}\tnon-spike reads: {real_reads_c}")

    return out_tup


def run_preprocessing(
        seq_files_t: tuple,
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
    if not any(seq_files_t):
        return sample_dir
    try:
        if not hasattr(seq_files_t, "__iter__") or len(seq_files_t) > 2:
            msg = "seq_files_t must be an iterable of max lenght 2"
            # PREP_LOG.error(msg)
            raise ValueError(msg)
        # Creating the preprocessing directory
        seq_files_t = [Path(PurePath(sfi)).absolute() for sfi in seq_files_t]
        if not any([True for fastq_fi_path in seq_files_t if fastq_fi_path.is_file()]):
            # returning if no provided fastq_file path is actually a file
            return sample_dir
        new_paths = ["", ""]  # [ForwardNewPath, ReverseNewPath]
        file_full_path = Path(PurePath(seq_files_t[0])).absolute()
        # creating processing directory
        sample_base_name = proc_helper.get_base_name(file_full_path)
        # Replacing sample_id if not given with sample_base_name
        sample_id = sample_id if sample_id else sample_base_name
        base_path_dir = file_full_path.parent.joinpath(sample_base_name)  # + PROC_DIR_SUFFIX).joinpath(proc_helper.generate_timestamp())
        sample_dir, shall_continue = should_trigger_processing(base_path_dir, args_yml_path)

        sample_dir = Path(PurePath(sample_dir))
        if not shall_continue:
            return sample_dir

        # TODO:NOTE decide if we need to process the sample or not
        # 1- argumnt check
        # 2- checking if some processed version already exist
        # 3- Choosing the one compliant to current argument set and return success. Otherwise new processing.

        sample_dir.mkdir(parents=True, exist_ok=True)
        # Updating fastq files path
        for ind, sfi in enumerate(seq_files_t):
            new_path = sample_dir.joinpath(sfi.name)
            new_paths[ind] = new_path
            # Copying files
            symlink(str(sfi), str(new_path))  # if sfi is symlink then new_path is symlink to sfi's target
        # We do spike removal if necessary and add a line to spike_stats file for spike normalization
        new_paths = [el for el in new_paths if bool(el)]
        # removing spikes and decompressing files below
        spike_reads_c, spike_stat_mapping_path, fastq_files_tuple = remove_spikes(
            new_paths,
            sample_id,
            sample_weight,
            spike_amount
        )
        # Copying the argument file to sample directory
        sample_arg_file = sample_dir.joinpath(DEFAULT_ARG_FILE_NAME)
        try:
            shutil.copy(args_yml_path, sample_arg_file)
        except shutil.SameFileError:
            pass
        # Running preprocessing
        sample_dir = preprocessing(
            input_dir=sample_dir,
            paired="Yes" if len(fastq_files_tuple) == 2 else "No",
            forward_file=fastq_files_tuple[0],
            reverse_file=fastq_files_tuple[1],
            input_id=sample_id,
            args_file_path=sample_arg_file,
            spike_amount=0,  # We run it always with 0 as we remove spikes before if there is
        )
    except Exception as exc:
        msg = f"{exc}"
        PREP_LOG.error(msg)
        shutil.rmtree(str(sample_dir))
        PREP_LOG.warning("Failed while processing {}, Deleting {}".format(sample_dir.name, sample_dir))

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
        args_yml_file: str = str(ARGS_YAML_FILE),
        dbs_dir: str = str(DBS_DIR),
        mapping_file: str = "",
        spike_stat_file: str = "",
        skip_preprocess: bool = False,
        skip_analysis: bool = False):
    fastq_file_dir = Path(PurePath(fastq_file_dir)).absolute()
    args_yml_file = Path(PurePath(args_yml_file)).absolute()
    dbs_dir = Path(PurePath(dbs_dir))
    assert fastq_file_dir.is_dir(), "fastq_file_dir must be a path to directory"
    assert dbs_dir.is_dir(), "dbs_dir must be a path to directory"
    assert args_yml_file.is_file(), "args_yml_file should be a valid path to YAML file reachable through input_dir"
    # Copying the args file to input_dir to always have the args file.
    # Just to be sure of being in the right place, KEEP THE NEXT LINE
    chdir(fastq_file_dir)
    try:
        shutil.copy(str(args_yml_file), str(fastq_file_dir.joinpath(DEFAULT_ARG_FILE_NAME)))
    except Exception:
        PREP_LOG.debug("Failed to copy given arguments file to {}".format(str(fastq_file_dir.joinpath(DEFAULT_ARG_FILE_NAME).relative_to(INPUT_DIR))))
        pass
    # preproc_dir = fastq_file_dir.joinpath("Preprocessing")
    # preproc_dir.mkdir(parents=True, exist_ok=True)
    default_spike_stat_compiled = fastq_file_dir.joinpath(SPIKE_STAT_FILE_NAME)
    given_spike_stat_file = Path(PurePath(spike_stat_file)).absolute()
    combined_spike_stats_path = default_spike_stat_compiled if not given_spike_stat_file.is_file() else given_spike_stat_file
    reduced_samples_dirs = []
    mapping_file_path = Path(PurePath(str(mapping_file))).absolute()  # it will return a path to current directory if mapping_file = ""
    # Analysis defualt vars
    zotu_file_path, sotu_file_path = ("", "")

    if not skip_preprocess:
        try:
            # Gathering sequence files
            PREP_LOG.debug("Gathering sequence files in {}".format(str(fastq_file_dir)))
            seq_file_pairs = proc_helper.pair_seq_files(str(fastq_file_dir))
            seq_file_pairs = [file_pair_tup for file_pair_tup in seq_file_pairs if bool(not is_in_processed_dir(file_pair_tup[0]) and not is_in_processed_dir(file_pair_tup[1]))]
            PREP_LOG.info("{} sampels were collected from {}.".format(str(len(seq_file_pairs)), str(fastq_file_dir)))
            # Filtering fastq files if they are already in processed directories
            # Sometimes failed processing leaves fastq files in the processing diectories

            mapping_line_tup_dict = parse_mapping_file(mapping_file_path, seq_file_pairs)
            # If mapping_file exists then we parse it and change the default of sample_weight, spike_mount to actual values.
            # running preprocessing
            with Pool(POOL_SIZE, maxtasksperchild=1) as pool:
                # Running the preprocessor function only if parese_mapping_file has return one or two files path for the sample
                res_list = [pool.apply_async(run_preprocessing, args=((arg_tup[1]), str(args_yml_file), arg_tup[0].SampleID, float(arg_tup[0].total_weight_in_g), int(arg_tup[0].spike_amount),))
                            for sample_id, arg_tup in mapping_line_tup_dict.items() if any(arg_tup[0])]
                for res in res_list:
                    res.wait(720)  # After 12 minutes it terminates the thread
            samples_dirs = [el.get() for el in res_list]
            # NOTE result form run_preprocessing could be "" which meand the preprocessing has failed
            samples_dirs = [el for el in samples_dirs if el]
        except Exception as exc:
            PREP_LOG.error(f"Preprocessing Failed: {exc}")
            skip_analysis = True
            samples_dirs = []
    else:
        PREP_LOG.warning(f"Skipping Preprocessing. Processing sammples in directory {fastq_file_dir}")
        samp_zotu_seq_files = gather_files(*[fastq_file_dir], file_name="taxed_ZOTUs.fasta")
        samples_dirs = [Path(PurePath(el)).parent for el in samp_zotu_seq_files]

    # TODO:NOTE the function to decide if the preprocessed samples should go for analysis or not come here.
    # we filter selected samples_dirs if they should be passed for analysis or not
    # 1- Argument set of selected sample should be the same as given one in the argument
    # 2- It should have "taxed_zotu.fasta"
    # Show warning if the sample is filtered out with a reason
    filtered_samples_dir = []
    for ind in range(len(samples_dirs)):
        curr_dir = samples_dirs[ind]
        if valid_for_analysis(curr_dir, args_yml_file):
            filtered_samples_dir.append(curr_dir)
    samples_dirs = filtered_samples_dir

    # Runing analysis
    if not skip_analysis:
        try:
            if not bool(samples_dirs):
                msg = "No samples were completely processed or had same processing "\
                      "argument set as given argument set. SKIPPING Analysisi!"
                raise ValueError(msg)
            analysis_dir = fastq_file_dir.joinpath(f"Analysis_{proc_helper.generate_timestamp()}")
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
            PREP_LOG.error(f"Analysis stopped: {exc}")

    ##################
    # TODO Only Normalizing
    if not skip_analysis and mapping_file_path.is_file() and zotu_file_path and sotu_file_path:
        pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    help_text = "The directory to find recursively all fastq files inside. It should be relative to <--input-directory>. Default is <--input-directory>"
    parser.add_argument("-i", "--input-directory",
                        type=str,
                        help="Every needed file and directory should be findable relative to this directory",
                        default=str(INPUT_DIR))  # WORKDIR of container is /base/inputs
    parser.add_argument("-d", "--fastq-directory",
                        type=str,
                        help=help_text,
                        default=str(INPUT_DIR))  # WORKDIR of container is /base/inputs
    # >BEGIN: Arguments need to be parsed from input and fastq directory
    parser.add_argument("-y", "--yml-file",
                        type=str,
                        help="Path to arguments yaml file. Relative to <--input-directory>",
                        default=None)
    parser.add_argument("-map", "--mapping-file",
                        type=str,
                        default=None,
                        help="The path to a mapping file defining sample weight and spike amount for each sample. Relative to <--input-directory>")
    parser.add_argument("-stat", "--spike-stat",
                        type=str,
                        default=None,
                        help=f"The path to a mapping file defining spike count, sample weight and spike amount for each sample.\n{SPIKE_STAT_HEADER}.\n"
                        "Relative to <--input-directory>")
    # <END
    parser.add_argument("-db", "--db-directory",
                        type=str,
                        help="Path to directory containing silva, sortmerna files. Relative to <--input-directory>",
                        default=DBS_DIR)
    parser.add_argument("-sp", "--skip-preprocess",
                        action="store_true",
                        help="Should skip preprocessing step")
    parser.add_argument("-sa", "--skip-analysis",
                        action="store_true",
                        help="Should skip analysis step")
    parser.add_argument("-tf", "--place-template-files",
                        action="store_true",
                        help=f"Writes the default argument yaml template file ({DEFAULT_ARG_FILE_NAME}) and mapping template file (mapping_file.tsv)"
                        " to <--input-directory>, print help text and exits.")
    args = parser.parse_args()
    # Updating INPUT_DIR
    INPUT_DIR = INPUT_DIR.joinpath(args.input_directory).absolute()
    FASTQ_DIR = INPUT_DIR.joinpath(args.fastq_directory).absolute()  # If they are the same it returns unchanged
    expected_yml_file = INPUT_DIR.joinpath(DEFAULT_ARG_FILE_NAME).absolute()
    expected_mapping_file = FASTQ_DIR.joinpath("mapping_file.tsv").absolute()
    cli_args_file = INPUT_DIR.joinpath(args.yml_file) if args.yml_file else expected_yml_file
    # This will return longest path. mapping file could be anywhere. Difining lower directories as fastq_directory will limit the searched files
    # and then less rows in mapping_file to be found.
    cli_map_file = INPUT_DIR.joinpath(args.mapping_file) if args.mapping_file else ""
    # Exposing default argument files.
    cli_spike_stat_file = INPUT_DIR.joinpath(args.spike_stat) if args.spike_stat else ""
    if not cli_args_file.is_file():
        shutil.copy(
            str(ARGS_YAML_FILE),
            str(cli_args_file)
        )
        PREP_LOG.warning(f"{str(cli_args_file.relative_to(INPUT_DIR))} is not a valid file path. "
                         f"Falling back to default argument set written to {str(cli_args_file.relative_to(INPUT_DIR))}")

    if args.place_template_files:
        shutil.copy(
            str(MAP_FILE),
            str(FASTQ_DIR.joinpath("mapping_file_TEMPLATE.tsv"))
        )
        # log
        PREP_LOG.info(f"Mapping file template files written to {str(FASTQ_DIR.joinpath('mapping_file_TEMPLATE.tsv').relative_to(INPUT_DIR))}")
        exit(0)
    else:
        if args.db_directory != DBS_DIR:
            DBS_DIR = args.db_directory
        # NOTE Every needed file and directory should be in /base/inputs/
        run_imngs2(
            fastq_file_dir=FASTQ_DIR,
            args_yml_file=cli_args_file,
            dbs_dir=INPUT_DIR.joinpath(args.db_directory),
            mapping_file=cli_map_file,
            spike_stat_file=cli_spike_stat_file,
            skip_preprocess=args.skip_preprocess,
            skip_analysis=args.skip_analysis,
        )
