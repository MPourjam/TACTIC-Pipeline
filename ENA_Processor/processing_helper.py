import pickle as pk
import yaml
import logging
from os import getcwd
from pathlib import Path, PurePath
from mimetypes import guess_type
from datetime import datetime as dt
from collections.abc import MutableMapping


def gimmelogger(logger_name: str, log_file: Path = None):
    # Logger
    logger_dir = Path(PurePath(getcwd())).name.stem
    log_file_path = logger_dir.joinpath(f"{logger_name}_log.txt")
    if log_file:
        custom_log_file = Path(str(log_file))
        if not custom_log_file.parent.is_dir():
            custom_log_file.parent.mkdir(parents=True, exist_ok=True)
        log_file_path = custom_log_file

    logger = logging.getLogger(logger_name)
    # set the logging level
    logger.setLevel(logging.DEBUG)

    # create a console and file handler
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    fh = logging.FileHandler(log_file_path)
    fh.setLevel(logging.WARNING)
    # create a formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    # set the formatter to the console handler
    ch.setFormatter(formatter)
    fh.setFormatter(formatter)
    # add the console handler to the logger
    logger.addHandler(ch)
    logger.addHandler(fh)

    return logger


global argparse_logger
argparse_logger = gimmelogger("ArgumentParser")


def flatten_dict(
        d: MutableMapping,
        parent_key: str = '',
        sep: str = ".") -> MutableMapping:
        items = []
        for k, v in d.items():
            new_key = parent_key + sep + k if parent_key else k
            if isinstance(v, MutableMapping):
                items.extend(ArgsParserUtil.flatten_dict(v, new_key, sep=sep).items())
            else:
                items.append((new_key, v))
        return dict(items)


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
               "Error": -1
               }
    # TODO Create the args_dict according to args_keys
    args_dict = {"input_dir": "",
                 "paired": "",
                 "forward_file": "",
                 "reverse_file": "",
                 "input_id": "",
                 "spike_amount": 6
                 }

    def check_args_dict(self, args_d=None):
        args_dict = self.task_dict["args"] if not args_d else args_d
        if not isinstance(args_dict, dict):
            return False
        keys = TaskPickle.args_keys
        has_key_list = [True for k in keys if k not in args_dict]
        if any(has_key_list):
            return False
        return True

    def check_status_dict(self, status_d=None):
        status_dict = self.task_dict["status"] if not status_d else status_d
        if not isinstance(status_dict, dict):
            return False
        keys = TaskPickle.status_keys_ess
        has_key_list = [True for k in keys if k not in status_dict]
        if any(has_key_list):
            return False
        return True

    def set_task_dict(self, task_d=None):
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

    def __init__(self, filepath=None, task_dict=None):
        status_dict = {k: TaskPickle.scode_d["Queue"] for k in TaskPickle.status_keys_num}
        for str_k in TaskPickle.status_keys_str:
            status_dict[str_k] = ""
        self.task_dict = {"args": TaskPickle.args_dict, "status": status_dict}
        timestamp = dt.now().strftime('%Y%m%d_%H%M%S')
        file_placeholder = "{}_proc_task_d.pk".format(str(timestamp))
        path = Path(filepath) if filepath else Path(getcwd()).joinpath(file_placeholder)
        self.__path = str(path.resolve())

        if Path(self.__path).exists():
            f_type, f_enc = guess_type(str(self.__path))
            if "x-tex-pk" not in f_type:
                raise FileExistsError("File is not a pickle file.")
            with open(self.__path, "rb") as pkfile:
                file_content = pk.load(pkfile)
            self.set_task_dict(file_content)

        if task_dict:
            self.set_task_dict(task_dict)

    def write(self):
        with open(self.__path, 'wb') as pkfile:
            pk.dump(self.task_dict, pkfile, protocol=pk.HIGHEST_PROTOCOL)

    def args_complete(self):
        arg_d = self.task_dict["args"]
        input_d = arg_d["input_dir"]
        paired = True if str(arg_d["paired"]).lower() == "yes" else False
        is_running = True if self.task_dict["status"]["run"] in [0, 2, 3] else False
        for k, v in arg_d.items():
            if k == "input_dir":
                try:
                    p = Path(v)
                except Exception:
                    return False
                if not p.exists():
                    return False
                if not p.is_dir():
                    return False
            elif k == "forward_file":
                p = Path(input_d)/str(v)
                # This is for when the processing is close to finish and fastqs are deleted.
                if not p.is_file() and not is_running:
                    return False
            elif k == "reverse_file":
                p = Path(input_d)/str(v)
                # This is for when the processing is close to finish and fastqs are deleted.
                if paired and not p.is_file() and not is_running:
                    return False
            elif k == "paired":
                if v not in ["Yes", "No"]:
                    return False
            elif k == "input_id":
                if not v:
                    return False
            elif k == "spike_amount":
                try:
                    v = int(v)
                except Exception:
                    return False
            else:  # Existence of any other arguments
                return False
        return True

    def __bool__(self):
        if not Path(self.__path).is_file():
            return False
        elif not self.check_args_dict():
            return False
        elif not self.check_status_dict():
            return False
        elif not self.args_complete():
            return False
        return True

    def __repr__(self):
        return str(self.__path)

    def __str__(self):
        return self.__repr__()

    @property
    def path(self):
        return str(self.__path)

    @path.setter
    def path(self, new_path):
        new_path = Path(new_path)
        old_path = Path(self.path)
        if old_path == new_path:
            return
        if new_path.is_dir():
            raise TypeError("Path of task pk object must point to a directory")
        self.__path = str(new_path)
        self.write()
        old_path.unlink()


class ArgsParserUtil:
    dict_flatt_sep = "__"

    def __init__(self,
                 args_d: dict = {},
                 parent_key_sep: str = "__"):
        """
        Initializing
        """
        if not isinstance(args_d, dict):
            raise ValueError(f"args_d must be a non-empty dictionary. {str(type(args_d))} is given.")
        # if not args_d:
        #     args_d = self.default_args
        for key, val in self.default_args.items():
            setattr(self, key, val)
            # addTaxUpdating the argument value if it exists in args_d
            arg_val_to_put = None
            # All args_d is supposed to be flattend with "__" as parent_key seperator
            for ky, vl in args_d.items():
                arg_parser_class_name = ky.split(ArgsParserUtil.dict_flatt_sep)[0:1]
                """
                NOTE
                If the argument parser part of key matches the class name of current object
                then the argument and it's value is put as an attribute.
                e.g: MergePairsArgs<dict_flatt_sep>fastq_maxdiffs will be only set as an
                attribute to a class parser with the name 'MergePairsArgs'
                """
                if key in ky and str(self.__class__.__name__) in arg_parser_class_name:
                    arg_val_to_put = vl
            if arg_val_to_put and isinstance(arg_val_to_put, type(val)):
                setattr(self, key, arg_val_to_put)

    def update_attrs(self, args_dict: dict = {}):
        """
        It takes a dictionary and updates the class instance attributes if
        a key of dictionary matches the class instance attributes.
        """
        for key, val in args_dict.items():
            # keys might be preceded by argument parser class name (e.g: MergePairsArgs)
            # TODO Taking key should follow the same logic as __ini__
            if hasattr(self, key) and isinstance(val, type(getattr(self, key, None))):
                setattr(self, key, val)

    @staticmethod
    def parse_yaml(yaml_path: str):
        """
        Parses a yaml to dictionary
        """
        yml_args_d = {}
        yaml_path = Path(yaml_path)
        with open(yaml_path, "r") as yaml_stream:
            try:
                yml_args_d = yaml.safe_load(yaml_stream)
            except yaml.YAMLError as e:
                argparse_logger.warning(e)

        return yml_args_d

    def __repr__(self):
        return '<%s.%s object at %s>' % (
            self.__class__.__module__,
            self.__class__.__name__,
            hex(id(self))
        )

    def __str__(self):
        return self.__repr__()


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
        "fastq_filter": "filtered1.fasta",
        "fastaout": "filtered2.fasta",
        "fastq_maxee_rate": 0.002,
    }


class FilterOneSideArgs(ArgsParserUtil):
    default_args = {
        # "fastq_filter": "filtered1.fasta",
        # "fastaout": "filtered2.fasta",
        # "fastq_trunclen": "[MINLEN-0.1MINLEN-5]",
        "fastq_truncqual": 10,
        "fastq_maxee_rate": 0.002,
        "fastq_maxee_rate": 0.002,
    }


class DereplicationArgs(ArgsParserUtil):
    default_args = {
        # "fastx_uniques": "filtered2.fasta",
        # "fsataout": "derep.fasta",
        "sizein": True,
        "sizeou": True,
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


class Filter16SArgs(ArgsParserUtil):
    default_args = {
        # "fastx": good-ZOTUs[.fasta],
        "ref90": "silva-bac-16s-id90.fasta",
        "ref95": "silva-arc-16s-id95.fasta",
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


class FilterZOTUAbundanceArgs(ArgsParserUtil):
    default_args = {
        # "otutab": "filtered2.fasta",
        # "zotus": "ZOTUs.fasta",
        # "otutabout": "zotu_table.txt",
        "abundance_cutoff": 0.0025,
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


class PreprocessingArgs:

    def __init__(
            self,
            config_yaml: str,
            config_dict: dict = {},
            *args,
            **kwargs):
        """
        This class reads a yaml file and parse it to a dictionary. If additional
        config_dict is given, this dictionary items gets added to the the dictionary
        created from yaml file without updating the yaml dictionary.
        """
        yml_args_dict = ArgsParserUtil.parse_yaml(config_yaml)
        yml_args_dict = flatten_dict(yml_args_dict, ArgsParserUtil.dict_flatt_sep)
        try:
            config_dict.update(yml_args_dict)
        except Exception as e:
            argparse_logger.warning(e)

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
        self.filter_zotu_abundance = FilterZOTUAbundanceArgs(config_dict)
        self.add_tax = AddTaxArgs(config_dict)


class YamlArgs:

    analysis_args_dict = {
        "somekey": "someattr"
    }

    @staticmethod
    def config_file_parser(confile_path: str) -> dict:
        """
        Parses and flatten the config file to a dictionary.
        """
        confile_path = Path(confile_path)

    def __init__(self,
                 preproc_arg_filepath: str = None,
                 analysis_arg_filepath: str = None,
                 args_dict: dict = None):
        """
        Initializing
        ============

        By giving the address to yml files containing arguments and their values
        it parses the arguments and sets them as attributes of classes.
        """
        pass
