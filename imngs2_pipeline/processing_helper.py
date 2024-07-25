import pickle as pk
import re
import yaml
import logging
import subprocess
import math
import gzip
import shutil
import zipfile
import inspect
import threading
from os import getcwd, makedirs, listdir, access, X_OK
from collections import namedtuple
from datetime import datetime as dt
from collections.abc import MutableMapping
from collections import Counter
import mimetypes as mtypes
from os import path as ospath
from os import remove
from sys import version_info
from copy import deepcopy
if version_info[0] < 3:
    from pathlib2 import Path, PurePath  # pip2 install pathlib2
else:
    from pathlib import Path, PurePath
try:
    # Posix based file locking (Linux, Ubuntu, MacOS, etc.)
    #   Only allows locking on writable files, might cause
    #   strange results for reading.
    import fcntl
    import os

    def lock_file(f):
        if f.writable():
            fcntl.lockf(f, fcntl.LOCK_EX | fcntl.LOCK_NB)

    def unlock_file(f):
        if f.writable():
            fcntl.lockf(f, fcntl.LOCK_UN)
except ModuleNotFoundError:
    # Windows file locking
    import msvcrt
    import os

    def file_size(f):
        return ospath.getsize(ospath.realpath(f.name))

    def lock_file(f):
        msvcrt.locking(f.fileno(), msvcrt.LK_RLCK, file_size(f))

    def unlock_file(f):
        msvcrt.locking(f.fileno(), msvcrt.LK_UNLCK, file_size(f))

# Setting
global bowtie2, SPIKESIDX
BIN_DIR = "/base/binaries/"
bowtie2 = BIN_DIR + "bowtie2/bowtie2"
SPIKESIDX = "/base/spikesidx/spike"

# The Start and end positions are according to reference Escherichia Coli K12 subst. MG1655
os_16S = namedtuple("pos_16S", "start end")
# pos_sil = namedtuple("pos_sil", "start end")
SixteenS_regions_dict = {
    "V1": [
        os_16S(69, 99),
        [1159, 1772]
    ],
    "V2": [
        os_16S(137, 242),
        [2083, 5291]
    ],
    "V3": [
        os_16S(433, 497),
        [9817, 10302]
    ],
    "V4": [
        os_16S(576, 682),
        [15643, 21773]
    ],
    "V5": [
        os_16S(822, 879),
        [25499, 26987]
    ],
    "V6": [
        os_16S(986, 1043),
        [31188, 32827]
    ],
    "V7": [
        os_16S(1117, 1173),
        [35463, 37685]
    ],
    "V8": [
        os_16S(1243, 1294),
        [40315, 40874]
    ],
    "V9": [
        os_16S(1435, 1465),
        [42608, 43016]
    ],
}

forw_file_indicators = [
    re.compile(r".*(_R1_).*"),
    re.compile(r".*(_1)\.(fastq|fq)(\.gz)?"),
    re.compile(r".*(_R1).*"),
    re.compile(r".*(_F_).*"),
    re.compile(r".*(@F).*"),
]

reve_file_indicators = [
    re.compile(r".*(_R2_).*"),
    re.compile(r".*(_2)\.(fastq|fq)(\.gz)?"),
    re.compile(r".*(_R2).*"),
    re.compile(r".*(_R_).*"),
    re.compile(r".*(@R).*"),
]

paired_files_indicators = list(zip(forw_file_indicators, reve_file_indicators))


class ArgsetException(Exception):
    pass


def calc_covered_region(start: int, end: int, regions_dict: dict = SixteenS_regions_dict):
    # Calculating Regions
    regions_coved = "NA"
    covered_regions = {
        k[1:] for k, v in regions_dict.items()
        if (start <= v[1][0] and end >= v[1][1])
    }
    if covered_regions:
        regions_coved = "V{}{}".format(
            min(covered_regions),
            [
                '-V%s' % (max(covered_regions))
                if (len(covered_regions) > 1)
                else ''
            ][0])

    return regions_coved


# overwriting system_sub() to run it with subprocess.run
def loud_subprocess(cmd_args_list: list, shell_bool: bool = False, cap_output: bool = False):
    dev_null_rgx = re.compile(r'(\s+)?([12]?>)(\s+)?(/dev/null|&1|&2)')
    if not isinstance(cmd_args_list, list):
        raise ValueError("cmd_args_list should be a list")
    cmd_string = "\t\t".join(cmd_args_list)
    # remove any null redirector from the given command
    # remove  > /dev/null 2>&1
    cmd_string = re.sub(dev_null_rgx, "", cmd_string)
    cmd_list = cmd_string.split("\t\t")
    # If capture_output_bool is true then do not redirect to /dev/null
    run_output = subprocess.run(
        cmd_list,
        capture_output=cap_output,
        encoding='utf-8',
        shell=shell_bool,
    )

    return run_output, cmd_list


class LogPrint(logging.Logger):
    def debug(msg):
        print(f"\033[95m{msg}\033[0m")

    def info(msg):
        print(f"\033[94m{msg}\033[0m")

    def warning(msg):
        print(f"\033[93m{msg}\033[0m")

    def error(msg):
        print(f"\033[91m{msg}\033[0m")

    def critical(msg):
        print(f"\033[1m\033[91m{msg}\033[0m")


def system_sub(cmd_args_list: list,
               force_log: bool = False,
               shell: bool = False,
               capture_output: bool = False,
               quiet: bool = False,
               logger_obj: logging.Logger = LogPrint("print_logger")):
    # logging
    if force_log:
        msg = f"COMMAND: {' '.join(cmd_args_list)}\n\n"
        logger_obj.info(msg)

    run_output, cmd_list = loud_subprocess(cmd_args_list, shell_bool=shell, cap_output=capture_output)
    if run_output.returncode == 137:  # Process killed due to memory limit
        err_msg = f"Command {' '.join(cmd_list)} exceeded memory limit."
        err_msg += "\n\tIf you are using usearch 32-bit version, consider upgrading to 64-bit version."
        err_msg += "\n\tIf you are using usearch 64-bit then run the programm with lower number of threads."
        raise MemoryError(err_msg)
    elif run_output.returncode != 0 and not quiet:
        raise Exception(f"Error in running command {' '.join(cmd_list)}:\n{run_output.stderr}")
    elif run_output.returncode != 0 and quiet:
        logger_obj.error(f"Error in running command {' '.join(cmd_list)}:\n{run_output.stderr}")
    elif run_output.returncode == 0 and run_output.stderr and not quiet:
        logger_obj.warning(f"Warning in running command {' '.join(cmd_list)}:\n{run_output.stderr}")
    return run_output


def gzip_to_fastq(*files) -> list:
    '''
    It takes files gunzip or zip them and returns the name with proper fastq fuffix.
    '''
    file_path = [Path(PurePath(f)) for f in files if f]
    file_path = [f for f in file_path if f.is_file()]
    if len(files) != len(file_path):
        raise TypeError("Some given arguments are not files.")
    fastq_paths = []
    for f in file_path:
        app, typ = mtypes.guess_type(f)
        file_dir, file_name = f.parent, str(f.name)
        file_name = file_name.replace(".", "_").replace("_gz", "").replace("_bz2", "")
        asciifile_name = file_name.replace(" ", "").replace("_fastq", "") + ".fastq"
        asciifile_name_normalized = ''.join(e for e in asciifile_name if e.isalnum() or e in ['_', '-', '.'])
        asciifile = str(file_dir.joinpath(asciifile_name_normalized))
        # derefrencing f if it's a link
        f_real = f.resolve()
        if 'zip' in str(typ):
            # For gzipped files
            with gzip.open(f_real, 'rb') as gz_file:
                with open(asciifile, 'wb') as ascii_file:
                    shutil.copyfileobj(gz_file, ascii_file)
                # Remove the original gzip file
            os.remove(f)
            fastq_paths.append(asciifile)
        elif 'zip' in str(app):  # For zipped files
            with zipfile.ZipFile(f_real, 'r') as zip_ref:
                zip_ref.extractall(asciifile)
            # Remove the original zip file
            os.remove(f)
            fastq_paths.append(asciifile)
        else:
            try:
                with open(f_real, "r+") as fi:
                    line = fi.readline()
                if len(line) == len(line.encode()):  # If it's ascii
                    try:
                        shutil.copy(f_real, asciifile)
                    except shutil.SameFileError:
                        pass
                    fastq_paths.append(asciifile)
            except Exception as e:
                print("{}\t{}".format(f, e))
                raise e
    if len(file_path) != len(fastq_paths):
        raise FileNotFoundError("Some files could not get converted or were not in utf-8 format!")
    # returns absolute paths
    return fastq_paths


def calc_spikes(*fastq_files, spike_amount: float = 0.0):
    '''
    #NOTE fastq_files MUST NOT be gzipped or zippped.
    spike_amount is in ng
    '''
    fastq_names = [ospath.abspath(f) for f in fastq_files if bool(f)]
    fastqs_abs_paths = [ospath.abspath(f) for f in fastq_names if ospath.isfile(f)]
    # Changing spike_amount to float
    spike_amount = float(spike_amount)
    if math.isclose(spike_amount, 0.0, abs_tol=1e-5) or math.isnan(spike_amount) or spike_amount < 0.0:
        f = fastqs_abs_paths[0]
        with open(f, 'r') as fqfile:
            line_count = 0
            line = fqfile.readline()
            while line:
                line_count += 1
                line = fqfile.readline()
        line_count = line_count // 4  # (total number of reads, spike reads)
        return line_count, 0  # non_spike_reads spike_reads
    else:
        spike_amount = float("{}e-9".format(spike_amount))  # ng to g
        name_regx = r"[a-zA-Z0-9\-]+"
        spike_res_dir = ospath.split(fastqs_abs_paths[0])[0] + "/spike_result/"
        makedirs(spike_res_dir, mode=777, exist_ok=True)
        _, forw_name = ospath.split(fastqs_abs_paths[0])
        fastq_aligned = forw_name[re.search(name_regx, forw_name).start():re.search(name_regx, forw_name).end()]
        fastq_aligned = ospath.join(spike_res_dir, fastq_aligned)
        fastq_unaligned = fastq_aligned + "_unal"

        if len(fastq_names) == 2:
            cmd = [bowtie2, "-x", SPIKESIDX, "-1", fastqs_abs_paths[0], "-2",
                   fastqs_abs_paths[1], "--al-conc", fastq_aligned, "--un-conc",
                   fastq_unaligned]
            call_out, cmd_list = loud_subprocess(cmd, cap_output=True)
            spike_counter = 0
            for fi in [fastq_aligned + ".1", fastq_aligned + ".2"]:
                with open(fi, 'r') as fastq_al:
                    for _ in fastq_al:
                        spike_counter += 1
            spike_counter = int(spike_counter // 4) // 2
            if spike_counter != 0:
                shutil.copy(fastq_unaligned + ".1", fastqs_abs_paths[0])
                shutil.copy(fastq_unaligned + ".2", fastqs_abs_paths[1])

        elif len(fastq_names) == 1:
            cmd = [bowtie2, "-x", SPIKESIDX, "-U", fastqs_abs_paths[0], "--al-conc",
                   fastq_aligned, "--un-conc", fastq_unaligned]
            # system(" ".join(cmd))
            call_out, cmd_list = loud_subprocess(cmd, cap_output=True)
            files_inspike = listdir(spike_res_dir)
            spike_counter = 0
            for fi in [f for f in files_inspike if re.search(fastq_aligned, ospath.abspath(f))]:
                with open(fi, 'r') as fastq_al:
                    for _ in fastq_al:
                        spike_counter += 1
            spike_counter = spike_counter // 4
            # Check if the unaligned file is empty
            if spike_counter != 0:
                shutil.copy(fastq_unaligned + ".1", fastqs_abs_paths[0])

        shutil.rmtree(spike_res_dir, ignore_errors=True)
        unaligned_reads, _ = calc_spikes(*fastqs_abs_paths, spike_amount=0.0)
        # to pass it as reads_number and spike number to seq_met model
        return unaligned_reads, spike_counter


def find_files_and_dirs_owned_by_root(directory):
    root_files_and_dirs = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            filepath = ospath.join(root, file)
            try:
                # Get file owner information
                file_stat = os.stat(filepath)
                file_owner = file_stat.st_uid
                # Check if the owner is root (uid 0)
                if file_owner == 0:
                    root_files_and_dirs.append(filepath)
            except Exception as e:
                root_files_and_dirs.append(filepath)
                argparse_logger.error(f"Error while checking for permissions on {filepath}: {e}")
        for dir_name in dirs:
            dir_path = ospath.join(root, dir_name)
            try:
                # Get directory owner information
                dir_stat = os.stat(dir_path)
                dir_owner = dir_stat.st_uid
                # Check if the owner is root (uid 0)
                if dir_owner == 0:
                    root_files_and_dirs.append(dir_path)
            except Exception as e:
                root_files_and_dirs.append(dir_path)
                argparse_logger.error(f"Error while checking for permissions on {dir_path}: {e}")
    return root_files_and_dirs


def generate_timestamp(thread_safe=False):
    current_time = dt.now()
    timestamp = current_time.strftime("%Y%m%d_%H%M%S")
    if thread_safe:
        # Get the current thread obj hex
        thread_obj = threading.current_thread()
        thread_obj = str(hash(id(thread_obj)))
        timestamp += f"_{thread_obj}"

    return timestamp


def get_base_name(file_path: Path, include_path: Path = None, replace_sep: tuple = ("", "")):
    file_path = Path(PurePath(file_path)).absolute()
    parent_path, file_name = file_path.parent, file_path.name
    if include_path is None or not isinstance(include_path, Path):
        include_path = parent_path

    while file_path.suffixes:
        file_path = Path(PurePath(file_path.parent)).absolute().joinpath(file_path.stem)

    for forw_ind_reg, reve_ind_reg in paired_files_indicators:
        forw_ind_reg_mo = forw_ind_reg.search(file_name)
        reve_ind_reg_mo = reve_ind_reg.search(file_name)
        if forw_ind_reg_mo:
            # split by first group of matched indicator
            file_name = str(file_name.split(str(forw_ind_reg_mo.group(1)))[0])
            break
        elif reve_ind_reg_mo:
            # split by first group of matched indicator
            file_name = str(file_name.split(str(reve_ind_reg_mo.group(1)))[0])
            break
    new_name = parent_path.joinpath(file_name).absolute().relative_to(include_path)
    new_name = str(Path(PurePath(str(new_name).replace(*replace_sep))))
    return new_name


def is_seq_file(file_path, seq_file_format="any"):
    """
    It checks if a file is FASTA or FASTQ sequence file.
    seq_file_format could be 'fasta', 'fastq' or 'any'
    """
    format_opts = [
        "fasta",
        "fastq",
        "any"
    ]
    if seq_file_format not in format_opts:
        raise ValueError("seq_file_format must be 'fasta', 'fastq' or 'any'")
    seq_head_tag_list = [">"] if seq_file_format == "fasta" else ["@"] if seq_file_format == "fastq" else [">", "@"]
    # use realpath to derefrence links to original files
    file_path = Path(PurePath(file_path)).resolve()
    if not file_path.is_file():
        return False
    app, typ = mtypes.guess_type(str(file_path))
    return_bool = []
    for seq_head_tag in seq_head_tag_list:
        try:
            if 'zip' in str(typ):  # For gzipped files
                with gzip.open(file_path, 'rt') as f:
                    return_bool.append(f.readline().startswith(seq_head_tag))
            elif 'zip' in str(app):  # For zipped files
                with zipfile.ZipFile(file_path, 'r') as zip_ref:
                    first_file = zip_ref.namelist()[0]
                    with zip_ref.open(first_file) as f:
                        return_bool.append(f.readline().decode().startswith(seq_head_tag))
            else:  # Text
                try:
                    with open(file_path, "r+") as fi:
                        line = fi.readline()
                        return_bool.append(len(line) == len(line.encode()) and line.startswith(seq_head_tag))
                except Exception as exc:
                    raise exc
        except Exception:
            return_bool.append(False)
    # print(file_path, return_bool)
    return any(return_bool)


def is_forward_file(file_name):
    return any([True for indi_reg in forw_file_indicators if bool(indi_reg.search(file_name))])


def find_reverse_file(file_path):
    """
    It assumes that reverse file exists in the same directory
    as forward file and the naming convention follows the same
    logic. i.e: _R2_ is the reverse of _R1_ and not of @F
    """
    file_path = Path(PurePath(file_path))
    file_name_parts = (file_path.name, file_path.name)
    reverse_name = ""
    rev_file_reg = re.compile(r"")
    stem_path = get_base_name(str(file_path))
    for forw_indic_reg, reve_indic_reg in paired_files_indicators:
        forw_indic_reg_mo = forw_indic_reg.search(file_path.name)
        if bool(forw_indic_reg_mo):
            file_name_parts = file_path.name.split(forw_indic_reg_mo.group(1))
            rev_file_reg = reve_indic_reg
            break
    # Finding most compatibe reverse file
    for caught_file in file_path.parent.glob(f"*{stem_path}*"):
        reve_reg_mo = rev_file_reg.search(caught_file.name)
        if reve_reg_mo and reve_reg_mo.groups():
            reve_file_parts = caught_file.name.split(reve_reg_mo.group(1))
            if reve_file_parts == file_name_parts:
                reverse_name = caught_file.name
                break

    reverse_file_path = file_path.parent.joinpath(reverse_name) if reverse_name else ""
    return reverse_file_path


def pair_seq_files(directory):
    direcotry_path = Path(PurePath(directory))
    paired_files = []
    if not direcotry_path.is_dir():
        return paired_files
    # is_seq_file checks the format of file as well.
    skip_reverses = []
    for file_path in direcotry_path.rglob("*"):
        file_path = Path(PurePath(file_path)).absolute()
        isSeqFile = is_seq_file(file_path, seq_file_format="fastq")
        if isSeqFile and str(file_path) not in skip_reverses:
            if is_forward_file(file_path.name):
                forw_file = str(file_path)
                reve_file = str(find_reverse_file(file_path))
                paired_file = (forw_file, reve_file)
                paired_files.append(paired_file)
                if reve_file:
                    skip_reverses.append(reve_file)

    return paired_files


def is_relative_to(path, *other):
    path = Path(PurePath(path))
    try:
        path.relative_to(*other)
        return True
    except ValueError:
        return False


def gimmelogger(logger_name: str = "", log_file: str = "", only_file: bool = True, propagate: bool = True):
    # finding caller file name and setting logger file path
    caller_frame = inspect.currentframe().f_back
    logger_name = Path(PurePath(caller_frame.f_code.co_filename)).stem if not logger_name else str(logger_name)
    if not log_file:
        log_file_path = Path(PurePath(inspect.getframeinfo(caller_frame).filename)).parent.joinpath(f"{logger_name}_log.txt")
    else:
        log_file_path = Path(PurePath(log_file)).absolute()

    # Setting logger formatter with line number of source of log
    formatter = logging.Formatter('%(asctime)s - %(name)s:%(filename)s:%(lineno)d - %(levelname)s - %(message)s')

    logger = logging.getLogger(logger_name)
    # set the logging level
    logger.setLevel(logging.DEBUG)
    has_stream_handler = [handler for handler in logger.handlers if isinstance(handler, logging.StreamHandler)]
    has_file_handler = [handler for handler in logger.handlers if isinstance(handler, logging.FileHandler) and handler.baseFilename == str(log_file_path)]
    if not only_file:
        ch = logging.StreamHandler() if not has_stream_handler else has_stream_handler[0]  # TODO what if there are several Streamhandler
        ch.setLevel(logging.DEBUG)
        # create a formatter
        # set the formatter to the console handler
        ch.setFormatter(formatter)
        # add the console handler to the logger
        logger.addHandler(ch)
    # Logger
    if log_file:
        if not log_file_path.parent.is_dir():
            log_file_path.parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(log_file_path, mode='w') if not has_file_handler else has_file_handler[0]
        fh.setLevel(logging.INFO)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
    # Propagate option
    logger.propagate = propagate
    return logger


global argparse_logger
argparse_logger = gimmelogger("run_imngs2.ArgumentParser")


def flatten_dict(
        d: MutableMapping,
        parent_key: str = '',
        sep: str = ".") -> MutableMapping:
    items = []
    for k, v in d.items():
        new_key = parent_key + sep + k if parent_key else k
        if isinstance(v, MutableMapping):
            items.extend(flatten_dict(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)


def slice_list(li: list, n: int) -> list:
    return [li[i:i + n] for i in range(0, len(li), n)]


class FileUtil:
    import gzip
    import zipfile
    import os
    from pathlib import Path, PurePath

    @staticmethod
    def is_gzip(file_path):
        file_path = Path(PurePath(file_path)).absolute()
        try:
            with gzip.open(file_path, 'rb') as f:
                f.read(1)
            return True
        except OSError:
            return False

    @staticmethod
    def is_zip(file_path):
        file_path = Path(PurePath(file_path)).absolute()
        try:
            with zipfile.ZipFile(file_path, 'r') as f:
                f.testzip()
            return True
        except zipfile.BadZipFile:
            return False

    @staticmethod
    def gunzip(file_path: str) -> str:
        """
        It gunzips the file and returns the path of the binary file.
        """
        file_path = Path(PurePath(file_path)).absolute()
        with gzip.open(file_path, 'rb') as f:
            file_content = f.read()
        bin_file = Path(PurePath(file_path).parent).joinpath(Path(PurePath(file_path).stem))
        with open(bin_file, 'wb') as f:
            f.write(file_content)
        if FileUtil.is_binary(bin_file):
            return str(bin_file)
        return ""

    @staticmethod
    def is_binary(file_path: str) -> bool:
        file_path = Path(PurePath(file_path)).absolute()
        if not file_path.is_file():
            return False
        # Check if the file is a binary file
        with open(file_path, "rb") as f:
            header = f.read(4)
            if header != b'\x7fELF':
                return False
        return True

    @staticmethod
    def change_mode(file_path: str, mode: int = 0o777):
        file_path = Path(PurePath(file_path)).absolute()
        file_path.chmod(mode)
        return mode

    @staticmethod
    def is_executable(file_path: str) -> bool:
        file_path = Path(PurePath(file_path)).absolute()
        if not file_path.is_file() or not access(file_path, X_OK):
            return False
        return True


class Usearch(FileUtil):
    import subprocess
    from pathlib import Path, PurePath

    desired_version = "usearch v11"

    def __init__(self, file: str, version: str = ""):
        self.file = Path(PurePath(file)).absolute()
        self.binfile = ""
        if not self.file.is_file():
            raise FileNotFoundError(f"File {file} does not exist.")
        self.desired_version = version if version else self.desired_version

    def which_usearch_version(self) -> str:
        """
        It returns the version of the usearch binary file.
        """
        version_cmd = [str(self.binfile), "--version"]
        version_out = subprocess.run(version_cmd, capture_output=True, text=True)
        return version_out.stdout

    def is_version(self, version: str = "") -> bool:
        """
        It checks if the given version is the same as the desired version.
        """
        version = version if version else self.desired_version
        return version in self.which_usearch_version()

    def check_or_get_bin(self) -> (int, str):
        """
        It checks the given file path and if it's a gzip file then it looks for a linux binary file
        made for i86 architecture and returns the path of the binary file.
        In case the file is not a gzip then it checks if the file is binary and made for i86 architecture
        and returns the path of the binary file.
        """
        ret_tup = (1, "")
        if FileUtil.is_binary(self.file):
            # if the file is a binary file
            self.binfile = self.file
        elif FileUtil.is_gzip(self.file):
            # If the file is a gzip file
            self.binfile = FileUtil.gunzip(self.file)
        else:
            # If the file is not a binary file
            # then we check if the binary file is inside
            bin_file = self.file.parent.joinpath(self.file.stem)
            if FileUtil.is_binary(bin_file):
                self.binfile = bin_file
            else:
                # if the binary file is not found
                self.binfile = ""

        if self.binfile:
            # if the binary file is found then we make it executable
            FileUtil.change_mode(self.binfile, 0o777)
            # if the version is not the desired version then we return 162
            ret_tup = (162, "")
            if self.is_version():
                ret_tup = (0, self.binfile)
        else:
            # if the binary file is not found then we return 161
            ret_tup = (161, "")

        return ret_tup

    def __bool__(self):
        return self.binfile != ""

    def __str__(self):
        return str(self.file)

    def __repr__(self):
        return f"Usearch {self.desired_version}" if self.binfile else "Usearch"


class TaskPickle:
    status_keys_num = ["download",
                       "run"]
    status_keys_str = ["msg"]
    status_keys_ess = status_keys_num
    status_keys = status_keys_num + status_keys_str
    args_keys = ["input_dir",
                 "paired",
                 "forward_file",
                 "reverse_file",
                 "input_id",
                 "spike_amount"]
    task_keys = ["args",
                 "status"]
    scode_d = {"Done": 0,
               "Queue": 1,
               "Started": 2,
               "Progress": 3,
               "Error": 4
               }
    # TODO Create the args_dict according to args_keys
    args_dict = {"input_dir": "",
                 "paired": "",
                 "forward_file": "",
                 "reverse_file": "",
                 "input_id": "",
                 "spike_amount": 6.0
                 }

    def check_args_dict(self, args_d=None) -> bool:
        """
        Checks if the arguments dictionary have all the arguments
        necessary for running the main() in processing_job.py
        if args_d is provided then the keys of provided dictionary
        get checked against the class-level args_keys dictionary
        otherwise the class instance's task_dict["args"] gets checked
        against class-level args_keys dictionary.
        """
        args_dict = self.task_dict["args"] if not args_d else args_d
        if not isinstance(args_dict, dict):
            return False
        keys = TaskPickle.args_keys
        has_key_list = [True for k in keys if k not in args_dict]
        if any(has_key_list):
            return False
        return True

    def check_status_dict(self, status_d=None) -> bool:
        """
        Checks if the status dictionary have all the essential keys
        (TaskPickle.status_keys_ess).
        if status_d is provided then the keys of provided dictionary
        get checked against the class-level status_keys_ess dictionary
        otherwise the class instance's task_dict["status"] gets checked
        against class-level status_keys_ess dictionary.
        """
        status_dict = self.task_dict["status"] if not status_d else status_d
        if not isinstance(status_dict, dict):
            return False
        keys = TaskPickle.status_keys_ess
        has_key_list = [True for k in keys if k not in status_dict]
        if any(has_key_list):
            return False
        return True

    def set_task_dict(self, task_d=None) -> bool:
        """
        Checks the structure of task_dict to be like:
        task_dict = {"args": {TaskPicle.args_keys: ...},
                     "status": {TaskPicle.status_keys_ess +
                                TaskPickle.status_keys_str: ...}
                     }
        and also performs the check_args_dict() and check_status_dict()
        If task_d is given then it checks whether the task_d complies with
        the aforementioned structure or not otherwise it checks the instance
        task_dict.
        """
        task_dict = self.task_dict if not task_d else task_d
        if isinstance(task_dict, dict):
            keys = TaskPickle.task_keys
            if any([True for k in keys if k not in task_dict]):
                keys_txt = ", ".join(keys)
                msg = "Task dict must have dictionaries: " + keys_txt
                raise TypeError(msg)
            elif not self.check_args_dict(task_dict["args"]):
                keys_txt = ", ".join(TaskPickle.args_keys)
                msg = "Argument dictionary must have arguments: " + keys_txt
                raise TypeError(msg)
            elif not self.check_status_dict(task_dict["status"]):
                keys_txt = ", ".join(TaskPickle.status_keys)
                msg = "Status dictionary must have keys: " + keys_txt
                raise TypeError(msg)
        else:
            msg = "Given object as task dictionary is not a dictionary."
            raise TypeError()
        self.task_dict = task_dict
        return True

    def __init__(self, filepath=None, task_dict=None) -> None:
        """
        Initialize the TaskPickle instance.
        If filepath is given it checks it to be a pickle file and it contents
        be according to the class task_dict structure usint set_task_dict().
        If task_dict is given it does the content check using set_task_dict() and
        sets the instance's task_dict to given task_dict
        """
        status_dict = {k: TaskPickle.scode_d["Queue"] for k in TaskPickle.status_keys_num}
        for str_k in TaskPickle.status_keys_str:
            status_dict[str_k] = ""
        self.task_dict = {"args": TaskPickle.args_dict, "status": status_dict}
        timestamp = dt.now().strftime('%Y%m%d_%H%M%S')
        file_placeholder = "{}_proc_task_d.pk".format(str(timestamp))
        self.__path = ospath.realpath(filepath) if filepath else ospath.realpath(getcwd()).join(file_placeholder)

        if ospath.isfile(self.__path):
            f_type, f_enc = mtypes.guess_type(str(self.__path))
            if "x-tex-pk" not in f_type:
                raise FileExistsError("File is not a pickle file.")
            # with open(self.__path, "rb") as pkfile:
            _file = self.__open_n_lock(mode="rb")
            file_content = pk.load(_file)
            _file.close()
            self.set_task_dict(file_content)
            self.__file = self.__open_n_lock(mode="w+b")
            self.write()
        else:
            self.__file = self.__open_n_lock(mode="wb")
            self.write()

        if task_dict:
            self.set_task_dict(task_dict)

    def __del__(self):
        # Closing the file will also opens the lock on the file
        try:
            self.__file.close()
            self.__file = False
        except Exception:
            pass

    def close(self):
        self.__del__()

    def write(self) -> None:
        """
        Writes the task_dict of instance to the instance path.
        """
        self.__file.seek(0)  # Setting the position to begining to overwrite
        pk.dump(self.task_dict, self.__file, protocol=pk.HIGHEST_PROTOCOL)
        self.__file.truncate(self.__file.tell())  # To get sure of EOF in correct position

    # @staticmethod
    def __open_n_lock(self, mode: str = "wb"):
        file_path = ospath.realpath(self.__path)
        file_o = open(file_path, mode)
        # Locking the pickle file to prevent other processes changing the file
        lock_file(file_o)
        return file_o

    def args_complete(self) -> bool:
        """
        Performs a set of checks and return boolean value indicating
        the completeness of task_dict["args"]
        """
        arg_d = self.task_dict["args"]
        input_d = arg_d["input_dir"]
        paired = True if str(arg_d["paired"]).lower() == "yes" else False
        is_running = True if self.task_dict["status"]["run"] in [0, 2, 3] else False
        for k, v in arg_d.items():
            if k == "input_dir":
                try:
                    p = ospath.realpath(v)
                except Exception:
                    return False
                if not ospath.exists(p):
                    return False
                if not ospath.isdir(p):
                    return False
            elif k == "forward_file":
                p = ospath.join(ospath.realpath(input_d), str(v))
                # This is for when the processing is close to finish and fastqs are deleted.
                if not ospath.isfile(p) and not is_running:
                    return False
            elif k == "reverse_file":
                p = ospath.join(ospath.realpath(input_d), str(v))
                # This is for when the processing is close to finish and fastqs are deleted.
                if paired and not ospath.isfile(p) and not is_running:
                    return False
            elif k == "paired":
                if v not in ["Yes", "No"]:
                    return False
            elif k == "input_id":
                if not v:
                    return False
            elif k == "spike_amount":
                try:
                    v = float(v)
                except Exception:
                    return False
            else:  # Existence of any other arguments
                return False
        return True

    def __bool__(self) -> bool:
        """
        Magic method for boolean checks.
        """
        if not ospath.isfile(self.__path):
            return False
        elif not self.check_args_dict():
            return False
        elif not self.check_status_dict():
            return False
        elif not self.args_complete():
            return False
        return True

    def __repr__(self) -> str:
        return str(self.__path)

    def __str__(self) -> str:
        return self.__repr__()

    @property
    def path(self):
        return str(self.__path)

    @path.setter
    def path(self, new_path):
        """
        Sets the path attribute of an instance and
        writes the task_dict to the new path.
        """
        new_path = ospath.realpath(new_path)
        old_path = ospath.realpath(self.path)
        if old_path == new_path:
            return
        if not ospath.isfile(new_path):
            raise TypeError("Path of task pk object must point to a file")
        self.__path = str(new_path)
        self.write()
        remove(old_path)


class ArgsParserDunderUtil:

    def __repr__(self):
        return '<%s.%s object at %s>' % (
            self.__class__.__module__,
            self.__class__.__name__,
            hex(id(self))
        )

    def __str__(self):
        dict_print = {}
        for ky, vl in self.__dict__.items():
            dict_print.update({ky: str(vl)})
        return str(dict_print)


class ArgsParserUtil(ArgsParserDunderUtil):
    dict_flatt_sep = "__"

    def __init__(self,
                 args_d: dict = {},
                 parent_key_sep: str = "__"):
        """
        Initializing
        """
        if not isinstance(args_d, dict):
            raise ArgsetException(f"args_d must be a non-empty dictionary. {str(type(args_d))} is given.")
        parsed_args = deepcopy(self.default_args)
        for key, val in self.default_args.items():
            setattr(self, key, val)
            parsed_args[key] = val
            # addTaxUpdating the argument value if it exists in args_d
            not_valid_value = "It's not yet a valid number"
            arg_val_to_put = not_valid_value
            # All args_d is supposed to be flattend with "__" as parent_key seperator
            for ky, vl in args_d.items():
                arg_parser_class_name = ky.split(ArgsParserUtil.dict_flatt_sep)[0:1]
                arg_parser_class_name = arg_parser_class_name[0] if arg_parser_class_name else str(None)
                arg_parser_key = ky.split(ArgsParserUtil.dict_flatt_sep)[1:2]
                arg_parser_key = str(None) if len(arg_parser_key) == 0 else str(arg_parser_key[0])
                """
                NOTE
                If the argument parser part of key matches the class name of current object
                then the argument and it's value is put as an attribute.
                e.g: MergePairsArgs<dict_flatt_sep>fastq_maxdiffs will be only set as an
                attribute to a class parser with the name 'MergePairsArgs'
                """
                if key == arg_parser_key and str(self.__class__.__name__) == arg_parser_class_name:
                    arg_val_to_put = vl
            if str(arg_val_to_put) != not_valid_value and isinstance(arg_val_to_put, type(val)):
                setattr(self, key, arg_val_to_put)
                parsed_args[key] = arg_val_to_put

        setattr(self, "parsed_args", parsed_args)

    def update_attrs(self, args_dict: dict = {}):
        """
        It takes a dictionary and updates the class instance attributes if
        a key of dictionary matches the class instance attributes.
        """
        for key, val in args_dict.items():
            # keys might be preceded by argument parser class name (e.g: MergePairsArgs)
            # TODO Taking key should follow the same logic as __ini__
            key_class = str(key).split(ArgsParserUtil.dict_flatt_sep)[0:1]
            key_class = key_class[0] if key_class else str(None)
            key_ = str(key).split(ArgsParserUtil.dict_flatt_sep)[1:2]
            key_ = key_[0] if key_ and str(self.__class__.__name__) in key_class else str(None)
            if hasattr(self, key) and isinstance(val, type(getattr(self, key, None))):
                setattr(self, key, val)

    @staticmethod
    def parse_yaml(yaml_path: str):
        """
        Parses a yaml to dictionary
        """
        yml_args_d = {}
        yaml_path = Path(PurePath(str(yaml_path))).absolute() if yaml_path else Path(PurePath("."))
        if not yaml_path.is_file():
            return yml_args_d
        with open(yaml_path, "r") as yaml_stream:
            try:
                yml_args_d = yaml.safe_load(yaml_stream)
            except yaml.YAMLError as e:
                argparse_logger.warning(e)
                raise ArgsetException(f"Error while parsing the yaml file: {yaml_path}")

        return yml_args_d

    def __eq__(self, other):
        if not isinstance(self, other.__class__):
            return False
        if not (hasattr(self, "parsed_args") or not hasattr(other, "parsed_args")):
            return False
        eq_tests = {fi: False for fi in self.parsed_args.keys()}
        for ind_key in self.parsed_args.keys():
            el = self.parsed_args.get(ind_key) == other.parsed_args.get(ind_key, "This is not valid at all")
            eq_tests[ind_key] = el
        return all(list(eq_tests.values()))


class SpikeRemovalArgs(ArgsParserUtil):
    default_args = {
        "spike_amount": 0,  # unit: ng
    }


class MergePairsArgs(ArgsParserUtil):
    default_args = {
        # "fastq_mergepairs": [FORWARD_FILE],
        # "reverse": [REVERSE_FILE],
        # "fastqout": merged.fasta,
        "fasq_maxdiffs": 50,
        "fastq_pctid": 50,
        "fastq_minmergelen": 200,
        "fastq_maxmergelen": 600,
    }


class TrimBothSidesArgs(ArgsParserUtil):
    default_args = {
        # "fastq_truncate": "merged.fasta",
        # "fastqout": "filtered1.fasta",
        "stripleft": 5,
        "stripright": 5,
    }


class TrimOneSideArgs(ArgsParserUtil):
    default_args = {
        "stripleft": 5,
    }


class FilterBothSidesArgs(ArgsParserUtil):
    default_args = {
        # "fastq_filter": "filtered1.fasta",
        # "fastaout": "filtered2.fasta",
        "fastq_maxee_rate": 0.002,
    }


class FilterOneSideArgs(ArgsParserUtil):
    default_args = {
        # "fastq_filter": "filtered1.fasta",
        # "fastaout": "filtered2.fasta",
        # "fastq_trunclen": "[MINLEN-0.1MINLEN-5]",
        "fastq_truncqual": 10,
        "fastq_maxee_rate": 0.002,
    }


class DereplicationArgs(ArgsParserUtil):
    default_args = {
        # "fastx_uniques": "filtered2.fasta",
        # "fsataout": "derep.fasta",
        "sizein": True,
        "sizeout": True,
    }


class SortArgs(ArgsParserUtil):
    default_args = {
        # "sortbysize": "derep.fasta",
        # "fastaout": "sorted.fasta",
    }


class ClusterZOTUsArgs(ArgsParserUtil):
    default_args = {
        # "unoise3": sorted.fasta,
        # "zotus": zotus.fasta,
        # "tabbedout": denoising.tab,
        "minsize": 2,
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # print(self.default_args)


class Filter16SArgs(ArgsParserUtil):
    default_args = {
        # "fastx": good-ZOTUs[.fasta],
        # "ref90": "silva-bac-16s-id90.fasta",
        # "ref95": "silva-arc-16s-id95.fasta",
        # "reads": "zotus.fasta",
        "other": "other.non16rRNA",
        "num_alignments": 1,
        "workdir": ".",
        "e": 0.1,
    }


class BuildZOTUTableArgs(ArgsParserUtil):
    default_args = {
        # "otutab": "filtered2.fasta",
        # "zotus": "ZOTUs.fasta",
        # "otutabout": "zotu_table.txt",
        "id": 0.97,
    }


class SelectZOTUsSeqsArgs(ArgsParserUtil):
    default_args = {
        # "fastx_getseqs": "good_ZOTUs.fa",
        # "labels": "filtered_zotu_table_list.txt",
        # "fastaout": "ZOTUs-Seqs.fasta",
    }


class AddTaxArgs(ArgsParserUtil):
    default_args = {
        # "in": "ZOTUs-Seqs.fasta",
        # "threads": 10,
        # "db": "<SINA_ARB>",
        # "out": "<FASTA_FILE>",
        "search": True,
        "meta-fmt": "csv",
        "lca-fields": "tax_slv",
        "turn": "all",
    }


class TrimSidesArgs(ArgsParserUtil):
    default_args = {
        # "fastq_truncate": "analysis.fasta",
        # "fastaout": "filtered1.fasta",
        "stripleft": 0,
        "stripright": 0,
    }


class ComplexTICArgs(ArgsParserUtil):
    default_args = {
        "family_sim": 0.90,
        "genus_sim": 0.95,
        "species_sim": 0.97,
    }


class CreateTableTICArgs(ArgsParserUtil):
    default_args = {
        "abund_limit": 0.0025,
        "sample_wise_correction": True,
    }


class PreprocessingArgsParser(ArgsParserDunderUtil):

    def __init__(
            self,
            config_yaml: str = "",
            config_dict: dict = {},
            *args,
            **kwargs):
        """
        This class reads a yaml file and parse it to a dictionary and finally arguments
        of preprocessing arguments needed for IMNGS2 Pipeline. If additional config_dict
        is given, this dictionary items gets added to the the dictionary created from
        yaml file without updating the yaml dictionary.
        """
        yml_args_dict = ArgsParserUtil.parse_yaml(config_yaml) if config_yaml else {}
        yml_args_dict = flatten_dict(yml_args_dict, sep=ArgsParserUtil.dict_flatt_sep)
        try:
            # NOTE values of the yaml file would overwrite values of config_dict
            config_dict.update(yml_args_dict)
        except Exception as e:
            argparse_logger.warning(e)
            raise ArgsetException("Error in updating the config_dict with yaml file.")

        # Updating attributes of class instance
        self.merge_pairs = MergePairsArgs(config_dict)
        self.trim_both_sides = TrimBothSidesArgs(config_dict)
        self.trim_one_side = TrimOneSideArgs(config_dict)
        self.filter_merged = FilterBothSidesArgs(config_dict)
        self.filter_single_reads = FilterOneSideArgs(config_dict)
        self.dereplication = DereplicationArgs(config_dict)
        self.sort_seq = SortArgs(config_dict)
        self.cluster_zotus = ClusterZOTUsArgs(config_dict)
        self.filter_16S = Filter16SArgs(config_dict)
        self.build_zotus_table = BuildZOTUTableArgs(config_dict)
        self.add_tax = AddTaxArgs(config_dict)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        eq_tests = [
            self.merge_pairs == other.merge_pairs,
            self.trim_both_sides == other.trim_both_sides,
            self.trim_one_side == other.trim_one_side,
            self.filter_merged == other.filter_merged,
            self.filter_single_reads == other.filter_single_reads,
            self.dereplication == other.dereplication,
            self.sort_seq == other.sort_seq,
            self.cluster_zotus == other.cluster_zotus,
            self.filter_16S == other.filter_16S,
            self.build_zotus_table == other.build_zotus_table,
            self.add_tax == other.add_tax,
        ]
        return all(eq_tests)


class AnalysisArgsParser(ArgsParserDunderUtil):

    def __init__(
            self,
            config_yaml: str = "",
            config_dict: dict = {},
            *args,
            **kwargs):
        """
        This class reads a yaml file and parse it to a dictionary and finally arguments
        of analysis needed for IMNGS2 Pipeline. If additional config_dict is given, this
        dictionary items gets added to the the dictionary created from yaml file without
        updating the yaml dictionary.
        """
        yml_args_dict = ArgsParserUtil.parse_yaml(config_yaml) if config_yaml else {}
        yml_args_dict = flatten_dict(yml_args_dict, sep=ArgsParserUtil.dict_flatt_sep)
        try:
            config_dict.update(yml_args_dict)
        except Exception as e:
            argparse_logger.warning(e)

        # TODO add the argument parser classes of Analysis here.
        self.trimsides = TrimSidesArgs(config_dict)
        self.complex_tic = ComplexTICArgs(config_dict)
        self.create_table = CreateTableTICArgs(config_dict)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        eq_tests = [
            self.trimsides == other.trimsides,
            self.complex_tic == other.complex_tic,
            self.create_table == other.create_table,
        ]
        return all(eq_tests)


class IMNGS2ArgsParser(ArgsParserDunderUtil):

    def __init__(
            self,
            config_yaml: str = "",
            config_dict: dict = {},
            *args,
            **kwargs):
        """
        Parses all arguments needed for IMNGS2 Pipeline
        """
        self.preproc_args = PreprocessingArgsParser(config_yaml, config_dict)
        self.analysis_args = AnalysisArgsParser(config_yaml, config_dict)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        prep_eq = self.preproc_args == other.preproc_args
        analysis_eq = self.analysis_args == other.analysis_args
        return prep_eq and analysis_eq


class MyCounter(Counter):

    def total(self):
        return sum([va for ke, va in self.items()])
