#! /usr/local/bin/python
import re
import signal
import argparse
import atexit
import shutil
import sys
from os import symlink, chdir
from imngs2_pipeline.processing_job import main_processing as preprocessing
from imngs2_pipeline.analyze_sample_job import main as main_analysis
from collections import namedtuple
import imngs2_pipeline.processing_helper as proc_helper
from imngs2_pipeline.processing_helper import gzip_to_fastq, calc_spikes, slice_list
from multiprocessing import cpu_count, Queue, Process
from typing import List, Tuple
if sys.version_info[0] < 3:
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
MAP_FILE = Path(PurePath("/base/mapping_file_TEMPLATE.csv")).absolute()
DBS_DIR = "/base/databases/"
max_pool = int(cpu_count() * 0.7)
POOL_SIZE = max_pool if max_pool > 0 else 1
# Preparing the logger
PREP_LOG = proc_helper.gimmelogger(
    "run_imngs2",
    log_file=INPUT_DIR.joinpath("Pipeline_log.txt"),
    only_file=False)
MAPPING_FILE_COLS = ("SampleID", "total_weight_in_g", "spike_amount", "parent_path")
MapLineTup = namedtuple("MapLineTup", [MAPPING_FILE_COLS[0], MAPPING_FILE_COLS[1], MAPPING_FILE_COLS[2], MAPPING_FILE_COLS[3]])
global SPIKE_STAT_FILE_NAME, SPIKE_STAT_HEADER, PATH_SEP, DEFAULT_MAP_LINE, SPIKE_STAT_FILE_COLS, TAXED_ZOTU_FILE_NAME
DEFAULT_MAP_LINE = MapLineTup("", float("NAN"), "0.0", "")
PATH_SEP = "_-_"
SPIKE_STAT_FILE_NAME = "spike_stat_mapping_file.csv"
SPIKE_STAT_FILE_COLS = ("#SampleID", "SpikeReads", "spikes_total_weight_in_g", "spike_amount", "parent_path")
SPIKE_STAT_HEADER = "\t".join(list(SPIKE_STAT_FILE_COLS))
TAXED_ZOTU_FILE_NAME = "taxed_ZOTUs.fasta"


def handle_system_signals(signum, frame):
    """
    Handles system signals to change the permission of created files and directories by the pipeline.
    """
    _ = correct_created_files_modes()
    # log the event
    PREP_LOG.warning(f"Signal {signum} received.")
    sys.exit(1)


def correct_created_files_modes(mode=755) -> bool:
    """
    It corrects the mode of created files by the pipeline.
    """
    mode = int("0o" + str(mode), 8)
    to_change_mode = proc_helper.find_files_and_dirs_owned_by_root(INPUT_DIR)
    for fi in to_change_mode:
        fi = Path(PurePath(fi)).absolute()
        try:
            fi.chmod(mode)
        except Exception as exc:
            PREP_LOG.warning(f"Could not change mode of {fi}. {exc}")

    return mode


# Register the signal handler
signal.signal(signal.SIGINT, handle_system_signals)
signal.signal(signal.SIGTERM, handle_system_signals)
# Registering exit function
atexit.register(correct_created_files_modes)


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
        TAXED_ZOTU_FILE_NAME,
        SPIKE_STAT_FILE_NAME
    ]
    logic_test = [False for fi in should_be_there]

    for ind in range(len(should_be_there)):
        dir_content = my_dir.iterdir()
        for cont in dir_content:
            if cont.is_file() and str(cont.name) == str(should_be_there[ind]):
                logic_test[ind] = True

    return all(logic_test)


def should_trigger_processing(sample_base_path: Path, given_argset_file: Path) -> Tuple[str, bool]:
    this_out = [None, True]
    given_argset_file = Path(PurePath(given_argset_file)).absolute() if isinstance(given_argset_file, str) or isinstance(given_argset_file, Path) else ""
    sample_base_path = Path(PurePath(sample_base_path)).absolute() if isinstance(sample_base_path, str) or isinstance(sample_base_path, Path) else ""
    if not given_argset_file or not sample_base_path:
        raise ValueError("Not Valid sample_base_path or not valid argument file. Both should be string or Path object.")

    processed_dir_path = Path(PurePath(str(sample_base_path))).absolute()

    if not str(sample_base_path).endswith(PROC_DIR_SUFFIX):
        processed_dir_path = Path(PurePath(str(sample_base_path) + PROC_DIR_SUFFIX))

    this_out[0] = str(processed_dir_path.joinpath(proc_helper.generate_timestamp()))
    arg_files = processed_dir_path.rglob(DEFAULT_ARG_FILE_NAME)
    # TODO:NOTE decide if we need to process the sample or not
    # 1- argument check
    # 2- checking if some processed version already exists
    # 3- Choosing the one compliant to current argument set and return success. Otherwise new processing.
    for argfile in arg_files:
        if not argfile.is_file():
            continue
        is_already_processed = is_processed_dir_healthy(argfile.parent) and not is_argset_different(argfile, given_argset_file)
        if is_already_processed:
            this_out[0] = str(argfile.parent.absolute())
            this_out[1] = False

    if not this_out[1]:
        PREP_LOG.info(f"Sample '{this_out[0]}' is already processed with given argument set. "
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
    while most_common and most_common.total() > 1:
        most_common_sample_id, most_common_count = list(most_common.items())[0]
        for ind in range(len(dup_ids_maplinetup_list)):
            map_line_obj = dup_ids_maplinetup_list[ind][0]
            file_pair = dup_ids_maplinetup_list[ind][1]
            if most_common[map_line_obj.SampleID]:
                new_map_line = MapLineTup(
                    str(most_common_sample_id) + "__" + str(most_common_count),
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

    return PurePath(input_path)


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
    # striping initial / or \ from given_parent_path
    given_parent_path = given_parent_path.strip().lstrip("\\").lstrip("/")
    given_parent_path = correct_path_for_windows(str(given_parent_path))
    given_parent_path = str(FASTQ_DIR) if given_parent_path == "" else Path(PurePath(base_path)).joinpath(given_parent_path)
    # return given_parent_path == fastq_file.parent
    return proc_helper.is_relative_to(fastq_file, given_parent_path)


def select_samples_for_analysis(mapping_line_tup_list: List[Tuple[Tuple[MapLineTup], Tuple[str, str]]], args_yml_path: str) -> list:
    """
    Only processed samples in the given mapping file will be sent for analysis.
    mapping_line_tup_list: list of tuples of mapping file lines and file pairs. e.g:
        [
            (
                MapLineTup(SampleID='truncSRR13005876_S1_L001', total_weight_in_g=1.0, spike_amount=6.0, parent_path='truncated'),
                ('/base/inputs/truncated/truncSRR13005876_S1_L001_R1_001.fastq', '/base/inputs/truncated/truncSRR13005876_S1_L001_R2_001.fastq')
            ),
            ...
        ]
    args_yml_path: path to the argument file
    """
    # find the taxed_zotu.fasta files in the processed samples
    # get the parent of the processed samples
    taxed_files_list = []
    for x in mapping_line_tup_list:
        # NOTE What if the name of processed directory does not match <SampleID> + PROC_DIR_SUFFIX
        # NOTE the sample ID in mapping_line_tup_list could be uniqified before if the base name of fastq files are not unique
        # with __<number> suffix that's why we use get_base_name to get the base name of fastq files
        sample_id = x[0][1]  # MapLineTup.SampleID
        fastq_files_parent_dir = Path(PurePath(x[1][0])).parent.joinpath(f"{str(sample_id)}{PROC_DIR_SUFFIX}")
        gathered_files = []
        for found_file in fastq_files_parent_dir.rglob(f"**/{TAXED_ZOTU_FILE_NAME}"):
            if check_taxed_zotus_fasta(found_file):
                gathered_files.append(found_file)
            else:
                PREP_LOG.warning(f"{str(found_file.relative_to(FASTQ_DIR))} file is not formatted correctly."
                                 " Proper fasta header format: '>SEQID;size=XXX;tax=KKK,PPP;CCC;OOO")

        taxed_files_list.extend(gathered_files)

    # So far we have all taxed_zotu.fasta files of given samples in mapping file
    # Now we need to filter out the samples that are not in the mapping file
    samples_dirs = [Path(PurePath(el)).parent for el in taxed_files_list]
    # TODO:NOTE the function to decide if the preprocessed samples should go for analysis or not come here.
    # we filter selected samples_dirs if they should be passed for analysis or not
    # 1- Argument set of selected sample should be the same as given one in the argument
    # 2- It should have "taxed_zotu.fasta"
    # Show warning if the sample is filtered out with a reason
    filtered_samples_dir = []
    for ind in range(len(samples_dirs)):
        curr_dir = samples_dirs[ind]
        if valid_for_analysis(curr_dir, args_yml_path):
            filtered_samples_dir.append(curr_dir)
        else:
            PREP_LOG.warning(f"{curr_dir.relative_to(FASTQ_DIR)} is not a correctly processed.")
    # If the forcing argument to reprocess samples are used then we might have several copies of same processed set
    # we should get sure that only one processed set from each sample is passed for analysis
    samples_dirs_non_redundant = {str(Path(PurePath(sam_dir))).split(PROC_DIR_SUFFIX)[0]: sam_dir for sam_dir in filtered_samples_dir}
    samples_dirs = list(samples_dirs_non_redundant.values())

    return samples_dirs


def parse_mapping_file(mapping_file_path: str, files_tups: list = []):
    """
    mapping_file_path: should be a path to mapping file
    files_tups: should be a list of files path tuples. i.e: [(<forward_path>, <reverse_path>), ...]
    """
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
        forw_file = Path(PurePath(file_pair[0])).absolute()
        sample_id = proc_helper.get_base_name(forw_file)
        curr_parent = forw_file.parent.relative_to(FASTQ_DIR)
        _row_vals_tup = MapLineTup(
            sample_id,
            DEFAULT_MAP_LINE.total_weight_in_g,
            DEFAULT_MAP_LINE.spike_amount,
            "" if str(curr_parent) == "." else str(curr_parent)
        )
        mapping_lines.append((_row_vals_tup, file_pair))
    # If mapping_file exists then we parse it and change the default of sample_weight, spike_mount to actual values.
    mapping_lines = convert_mapping_entries_to_dict(mapping_lines)
    if not mapping_file_path:
        PREP_LOG.info("No Mapping file provided. Using default values for sample preprocessing.")
        return mapping_lines
    mapping_file_path = Path(PurePath(mapping_file_path)).absolute()
    # if mapping_file_path is not empty string and is not a file then raise an error
    if not mapping_file_path.is_file():
        PREP_LOG.error(f"Mapping file: {mapping_file_path} does not exist!")
        raise FileNotFoundError(f"Mapping file: {mapping_file_path} does not exist!")

    PREP_LOG.info(f"When mapping file is provided ({mapping_file_path}), only samples in mapping file will get processed!!!")
    # If mapping_file_is there then renew the mapping_lines
    mapping_lines = []  # List[(MapLineTup, file_pair)] -> from mapping file, full line for each relevant sample

    with open(mapping_file_path, 'r') as mapping_file_h:
        header: str = next(mapping_file_h)

        if not header.startswith('#'):
            PREP_LOG.error('No header in mapping file')
            raise Exception('No header in mapping file')

        columns = header.strip().split('\t')
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
        index_of_parent_path_col = index_of_parent_path_col if index_of_parent_path_col else len(columns)
        # looping over lines
        for line in mapping_file_h:
            line = line.strip().rstrip()
            if line.startswith('#'):
                continue
            # placing default values
            hypo_fields = [
                "",
                DEFAULT_MAP_LINE.total_weight_in_g,  # float
                DEFAULT_MAP_LINE.spike_amount,  # String float
                DEFAULT_MAP_LINE.parent_path  # String
            ]
            # Parsing lines
            fields = line.split('\t')
            try:
                sample_id = str(fields[index_of_sample_id_col]).strip()
                hypo_fields[0] = sample_id
            except Exception:
                raise ValueError(f"Could not parse SampleID from mapping file. Error on line: {line}")
            try:
                total_weight_in_g = str(fields[index_of_weight_col]).strip()
                # Check if the float value is in german format and if change the format to english
                total_weight_in_g = total_weight_in_g.replace(",", ".")
                hypo_fields[1] = float(total_weight_in_g)
            except Exception:
                raise ValueError(f"Could not parse weight amount from mapping file. Line: {line}\n"
                                 f" Check the mapping file format. Columns should be"
                                 f" separated by TAB. Put NAN for weight if you don't have weight.")
            try:
                amount = str(fields[index_of_amounts_col]).strip()
                # Check if the float value is in german format and if change the format to english
                amount = amount.replace(",", ".")
                hypo_fields[2] = str(float(amount))
            except Exception:
                raise ValueError(f"Could not parse spike amount from mapping file. Line: {line}\n"
                                 f" Check the mapping file format. Columns should be"
                                 f" separated by TAB. Put default value of 0 if you don't have spike amount.")
            try:
                parent_path = str(fields[index_of_parent_path_col])
                parent_path = parent_path.strip().rstrip("\\").rstrip("/")
                hypo_fields[3] = parent_path
            except Exception:
                PREP_LOG.info(f"Could not parse parent path from mapping file. Line: {line}\n"
                              f" Check the mapping file format. Columns should be"
                              f" separated by TAB. Continuing with default value of 0.\n"
                              f" If you are using a mapping_file without 'parent_path' column then ignroe this warning.")

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
                    # update parent path to the correct one
                    map_line = MapLineTup(
                        map_line.SampleID,
                        map_line.total_weight_in_g,
                        map_line.spike_amount,
                        str(Path(PurePath(file_tup[0])).parent.relative_to(FASTQ_DIR))
                    )
                    file_pairs_tup = file_tup
                    del files_tups[ind]
                    break
            found_fastqs = [True if fi else False for fi in file_pairs_tup]
            # NOTE: We should warn the user if we consider a sample being single layout and not paired layout
            if not all(found_fastqs):
                PREP_LOG.warning(f"Could not match fastq  to sample_id: {sample_id} in mapping file.")

            mapping_lines.append((map_line, file_pairs_tup))
    mapping_lines = convert_mapping_entries_to_dict(mapping_lines)

    return mapping_lines


def remove_spikes(
        files_paths: list,
        sample_id: str,
        sample_weight: float,  # in gram
        spike_amount: float):
    """
    It remove spikes from a sample and adds the number of spikes to a mapping file given.
     the column name should be "SpikeReads". The format of spike_stat_mapping_path must be
     like below:
     #SampleID\tSpikeReads\tspikes_total_weight_in_g\tspike_amount\n'

    Arguments:
        - files_path: list = list of files on which spike removal should be applied
        - sample_id: str = sample identifier to be reported in spike_stat_mapping_path
        - sample_weight: float = sample weight sent for sequencing in gram
        - spike_amount: float = amount of spike in nanogram

    """
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

    # putting parent path into the spike_mapping_file
    parent_path = Path(PurePath(files_paths_abs[0])).parent.relative_to(FASTQ_DIR)

    # postponing file opening to avoid race condition if ran parallel
    real_reads_c, spike_reads_c = calc_spikes(*fastq_files_abs_path, spike_amount=spike_amount)
    out_tup = (spike_reads_c, spike_stat_mapping_path, tuple(fastq_files_abs_path))
    with open(spike_stat_mapping_path, 'a') as stats_h:
        # appending\
        new_row = [
            str(sample_id),
            str(spike_reads_c),
            str(sample_weight),
            str(spike_amount),
            str(parent_path) if str(parent_path) != "." else ""
        ]
        stats_h.write("\t".join(new_row) + "\n")

    PREP_LOG.info(f"Spike removal done for {sample_id}. Spike Reads: {spike_reads_c}\tnon-spike reads: {real_reads_c}")

    return out_tup


def run_preprocessing(
        seq_files_t: tuple,
        args_yml_path: str,
        usearch_11_bin: str,
        sample_id: str,
        sample_weight: float = float("NAN"),  # spike normalizer handles this
        spike_amount: float = 0.0):
    """
    It takes a tuple of paths to sequencing files.
    Create directory for basename of files and move
     them there. Then triggers the preprocessing step.
     It returns the path of direcotry in which the
     preprocessd files of sample are.
    """
    sample_dir = None
    if not any(seq_files_t):
        return sample_dir
    try:
        if not hasattr(seq_files_t, "__iter__") or len(seq_files_t) > 2:
            msg = "seq_files_t must be an iterable of max length 2"
            # PREP_LOG.error(msg)
            raise ValueError(msg)
        # Creating the preprocessing directory
        seq_files_t = [Path(PurePath(sfi)).absolute() for sfi in seq_files_t]
        seq_files_t = [sfi for sfi in seq_files_t if sfi.is_file()]
        if not any(seq_files_t):
            # returning if no provided fastq_file path is actually a file
            return sample_dir
        new_paths = ["", ""]  # [ForwardNewPath, ReverseNewPath]
        file_full_path = seq_files_t[0]
        base_path_dir = file_full_path.parent.joinpath(sample_id)  # + PROC_DIR_SUFFIX).joinpath(proc_helper.generate_timestamp())
        sample_dir, shall_continue = should_trigger_processing(base_path_dir, args_yml_path)
        sample_dir = Path(PurePath(sample_dir)) if sample_dir else ""
        if not shall_continue:
            return sample_dir

        sample_dir.mkdir(parents=True, exist_ok=True)
        # Updating fastq files path
        for ind, sfi in enumerate(seq_files_t):
            new_path = sample_dir.joinpath(sfi.name)
            new_paths[ind] = new_path
            # Copying files
            symlink(str(sfi), str(new_path))  # if sfi is symlink then new_path is symlink to sfi's target

        # We do spike removal if necessary and add a line to spike_stats file for spike normalization
        # We only carry valid files in fastq_files_tuple
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
            forward_file=fastq_files_tuple[0],
            reverse_file=fastq_files_tuple[1] if len(fastq_files_tuple) == 2 else "",
            input_id=sample_id,
            args_file_path=sample_arg_file,
            spike_amount=0,  # We run it always with 0 as we remove spikes before if there is
            usearch_11_bin=usearch_11_bin
        )
    except MemoryError as mem_exc:
        err_msg = f"{mem_exc}"
        err_msg += "\n\tIf you are using usearch 32-bit version, consider upgrading to 64-bit version."
        err_msg += "\n\tIf you are using usearch 64-bit then run the programm with lower number of threads."
        PREP_LOG.error(f"{err_msg}")
        # If parent of sample_dir is empty then remove it
        if not list(sample_dir.parent.iterdir()):
            shutil.rmtree(str(sample_dir.parent))
        sample_dir = None
        sys.exit(137)
    except Exception as exc:
        msg = f"{exc}"
        PREP_LOG.error(msg)
        PREP_LOG.warning("Failed while processing {}, Deleting !!!".format(sample_dir.relative_to(FASTQ_DIR)))
        # If parent of sample_dir is empty then remove it
        if not list(sample_dir.parent.iterdir()):
            shutil.rmtree(str(sample_dir.parent))
        sample_dir = None

    return sample_dir


def parallel_preprocessing(args: tuple, process_queue: Queue, res_dict_key: str):
    """
    It's a parallel function to run preprocessing on multiple samples.
    """
    try:
        sample_dir = run_preprocessing(*args)
        process_queue.put((res_dict_key, sample_dir))
    except Exception as exc:
        process_queue.put((res_dict_key, None))
        PREP_LOG.error(f"Failed to run preprocessing for {args[2]}. {exc}")

    return


def combine_spike_stats_file(samples_dirs: list, combined_spike_stat: str):
    """
    It combines samples' spike_stat_file to one file to be input to spike normalization step
    """
    samples_dirs_path = []
    combined_spike_stat = Path(PurePath(combined_spike_stat)).absolute()
    for sam_dir in samples_dirs:
        sam_dir_path = Path(PurePath(sam_dir)).absolute()
        if not sam_dir_path.is_dir():
            raise FileNotFoundError(f"{sam_dir_path} does not exist!")
        if not sam_dir.joinpath(SPIKE_STAT_FILE_NAME).is_file():
            raise FileNotFoundError(f"{sam_dir_path} does not have spike_stat_file: {SPIKE_STAT_FILE_NAME}")
        # cleaned sample dirs
        samples_dirs_path.append(sam_dir_path)
        PREP_LOG.info(f"Found spike_stat_file in {sam_dir}. Listed for Combining...")

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
                msg = f"Existing spike_stat_file: {ssfp} does not have correct format of: {SPIKE_STAT_HEADER}"
                PREP_LOG.error(msg)
                raise ValueError(msg)
            # Updating parent path to include the parent of sample_stat_files_path
            parent_path = ssfp.parent.relative_to(FASTQ_DIR)
            parent_path = "" if str(parent_path) == "." else str(parent_path)
            # Find the corresponding value of each column from ssfp
            row_vals_list = row_values.split("\t")
            new_row_values = [row_vals_list[0], row_vals_list[1], row_vals_list[2], row_vals_list[3], parent_path]
            combined_stat_file_fio.write("\t".join(new_row_values) + "\n")

    return samples_dirs_path, combined_spike_stat


def check_taxed_zotus_fasta(file_path: str) -> bool:
    format_is_correct = True
    size_rgx = re.compile(r"size=([0-9]+);")
    tax_rgx = re.compile(r"tax=([^;]+){0,7};")
    with open(file_path) as fs_fio:
        line = fs_fio.readline().strip()
        size_mo = size_rgx.search(line)
        tax_mo = tax_rgx.search(line)
        if not all([line, line.startswith(">"), size_mo, tax_mo]):
            format_is_correct = False

    return format_is_correct


def run_imngs2(
        fastq_file_dir: str,
        usearch_11_bin: str,
        args_yml_file: str = str(ARGS_YAML_FILE),
        dbs_dir: str = str(DBS_DIR),
        mapping_file: str = "",
        spike_stat_file: str = "",
        skip_preprocess: bool = False,
        skip_analysis: bool = False):
    # NOTE if this function is imported then the default global variables will be used
    global USEARCH_11_BIN, FASTQ_DIR
    FASTQ_DIR = Path(PurePath(fastq_file_dir)).absolute()
    ret_code, usearch_bin_path = proc_helper.Usearch(usearch_11_bin).check_or_get_bin()
    if ret_code != 0:
        PREP_LOG.error(f"Failed to find usearch binary: {usearch_11_bin}")
        sys.exit(ret_code)
    USEARCH_11_BIN = Path(PurePath(usearch_bin_path)).absolute()
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
    except shutil.SameFileError:
        pass
    except Exception:
        PREP_LOG.debug("Failed to copy given arguments file to {}".format(str(fastq_file_dir.joinpath(DEFAULT_ARG_FILE_NAME).relative_to(INPUT_DIR))))
        sys.exit(1)
    # preproc_dir = fastq_file_dir.joinpath("Preprocessing")
    # preproc_dir.mkdir(parents=True, exist_ok=True)
    default_spike_stat_compiled = fastq_file_dir.joinpath(SPIKE_STAT_FILE_NAME)
    given_spike_stat_file = Path(PurePath(spike_stat_file)).absolute()
    combined_spike_stats_path = default_spike_stat_compiled if not given_spike_stat_file.is_file() else given_spike_stat_file
    reduced_samples_dirs = []
    mapping_file_path = Path(PurePath(str(mapping_file))).absolute() if mapping_file else ""  # it will return a path to current directory if mapping_file = ""
    # Analysis default vars
    zotu_file_path, sotu_file_path = ("", "")
    # Preparing mapping file entries
    # Gathering sequence files
    PREP_LOG.debug("Gathering sequence files in {}".format(str(fastq_file_dir)))
    seq_file_pairs = proc_helper.pair_seq_files(str(fastq_file_dir))
    seq_file_pairs = [file_pair_tup for file_pair_tup in seq_file_pairs if bool(not is_in_processed_dir(file_pair_tup[0]) and not is_in_processed_dir(file_pair_tup[1]))]
    PREP_LOG.info("{} samples were collected from {}.".format(str(len(seq_file_pairs)), str(fastq_file_dir)))
    # Filtering fastq files if they are already in processed directories
    # Sometimes failed processing leaves fastq files in the processing directories
    mapping_line_tup_dict = parse_mapping_file(mapping_file_path, seq_file_pairs)
    if not skip_preprocess:
        try:
            # If mapping_file exists then we parse it and change the default of sample_weight, spike_mount to actual values.
            # running preprocessing
            samples_dirs = []
            res_dict = {}
            task_batches = slice_list(list(mapping_line_tup_dict.items()), POOL_SIZE)
            # create a queue for the tasks

            for task_b in task_batches:
                preproc_queue = Queue()
                this_batch = []
                for sample_id, arg_tup in task_b:
                    proc = Process(
                        target=parallel_preprocessing,
                        args=(
                            (
                                arg_tup[1],  # tuple of fastq files
                                str(args_yml_file),  # path to args file
                                USEARCH_11_BIN,
                                arg_tup[0].SampleID,  # sample_id
                                float(arg_tup[0].total_weight_in_g),  # sample_weight
                                float(arg_tup[0].spike_amount),   # spike_amount
                            ),
                            preproc_queue,
                            arg_tup[0].SampleID,  # res_dict_key must be unique to bound process to sample_id
                        )
                    )
                    proc.start()
                    this_batch.append(proc)
                for proc in this_batch:
                    proc.join(timeout=2400)  # 40 minutes of waiting
                    if proc.is_alive():
                        proc.terminate()
                        proc.join()
                while not preproc_queue.empty():
                    sample_id, res = preproc_queue.get()
                    if res:
                        res_dict[sample_id] = res
                    else:
                        PREP_LOG.warning(f"Failed to finalize preprocessing for {sample_id}")
            samples_dirs = list(res_dict.values())
        except Exception as exc:
            PREP_LOG.error(f"Preprocessing Failed: {exc}")
            skip_analysis = True
            samples_dirs = []
            sys.exit(1)
    else:
        PREP_LOG.warning(f"Skipping Preprocessing. Processing samples in directory {fastq_file_dir}")
        samples_dirs = select_samples_for_analysis(list(mapping_line_tup_dict.values()), args_yml_file)
    # Runing analysis
    if not skip_analysis and bool(samples_dirs):
        try:
            if not bool(samples_dirs):
                msg = "No samples were completely processed or had same processing "\
                      "argument set as given argument set. SKIPPING Analysis!"
                raise ValueError(msg)
            analysis_dir = fastq_file_dir.joinpath(f"Analysis_{proc_helper.generate_timestamp()}")
            analysis_dir.mkdir(parents=True, exist_ok=True)
            # we do the step down because if skip_preprocess is True then the path to preprocess would be given not all smaple_dirs
            if str(combined_spike_stats_path) != str(default_spike_stat_compiled):
                PREP_LOG.warning(f"Using custom spike_stat file: {combined_spike_stats_path} for spike normalization!! Default spike_stat file {default_spike_stat_compiled} is ignored!")
            reduced_samples_dirs, _ = combine_spike_stats_file(samples_dirs, combined_spike_stat=combined_spike_stats_path)
            missed_samples = [sam_dir for sam_dir in samples_dirs if sam_dir not in reduced_samples_dirs]
            if missed_samples:
                PREP_LOG.warning(f"Samples in {str(combined_spike_stats_path)} will be sent for analysis. Please Check the file.")
                PREP_LOG.warning(f"Missed Samples are: {missed_samples}")
            zotu_file_path, sotu_file_path = main_analysis(
                analysis_dir=analysis_dir,
                spike_stat_file=combined_spike_stats_path,
                fastqs_dir=str(FASTQ_DIR),
                args_file_path=args_yml_file,
                dbs_loc=dbs_dir,
                threads=POOL_SIZE,
                usearch_11_bin=USEARCH_11_BIN,
            )
        except MemoryError as exc:
            PREP_LOG.error(f"Analysis stopped: {exc}")
            sys.exit(137)
        except Exception as exc:
            PREP_LOG.error(f"Analysis stopped: {exc}")
            sys.exit(1)

    ##################
    # TODO Only Normalizing
    if not skip_analysis and mapping_file_path and zotu_file_path and sotu_file_path:
        pass

    # change mode of files
    correct_created_files_modes()

    sys.exit(0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    help_text = "The directory to find recursively all fastq files inside. It should be relative to <--input-directory>. Default is <--input-directory>"
    parser.add_argument("-i", "--input-directory",
                        type=str,
                        help="Every needed file and directory should be findable relative to this directory",
                        default=".")  # WORKDIR of container is /base/inputs
    parser.add_argument("-d", "--fastq-directory",
                        type=str,
                        help=help_text,
                        default=".")  # WORKDIR of container is /base/inputs
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
    parser.add_argument("-ut", "--usearch-bin",
                        type=str,
                        help="Path to binary of usearch version 11. Default is 11.0.667_i86linux32.",
                        required=True)
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
                        help=f"Writes the default argument yaml template file ({DEFAULT_ARG_FILE_NAME}) and mapping template file (mapping_file.csv)"
                        " to <--input-directory>, print help text and exits.")
    parser.add_argument("-t", "--threads",
                        type=int,
                        help="Number of threads to use for parallel processing",
                        default=POOL_SIZE)
    args = parser.parse_args()
    # Updating INPUT_DIR
    INPUT_DIR = INPUT_DIR.joinpath(args.input_directory).absolute()
    # If they are the same it returns unchanged. If fastq_dir is subpath of input it returns the longest one
    FASTQ_DIR = INPUT_DIR.joinpath(args.fastq_directory).absolute()
    # TODO Later we need to force the user to provide the path to usearch binary. For now we only continue with the default one.
    given_usearch_bin = str(INPUT_DIR.joinpath(args.usearch_bin).absolute())
    # warning the cli users for the given usearch file
    try:
        POOL_SIZE = int(args.threads)
    except ValueError as exc:
        PREP_LOG.warning(f"Could not parse threads argument. Using default value of {POOL_SIZE}. Error: {exc}")
    expected_yml_file = INPUT_DIR.joinpath(DEFAULT_ARG_FILE_NAME).absolute()
    expected_mapping_file = FASTQ_DIR.joinpath("mapping_file.csv").absolute()
    cli_args_file = INPUT_DIR.joinpath(args.yml_file) if args.yml_file else expected_yml_file
    # This will return longest path. mapping file could be anywhere. Defining lower directories as fastq_directory will limit the searched files
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
            str(FASTQ_DIR.joinpath("mapping_file_TEMPLATE.csv"))
        )
        # log
        PREP_LOG.info(f"Mapping file template files written to {str(FASTQ_DIR.joinpath('mapping_file_TEMPLATE.csv').relative_to(INPUT_DIR))}")
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
            usearch_11_bin=given_usearch_bin
        )
