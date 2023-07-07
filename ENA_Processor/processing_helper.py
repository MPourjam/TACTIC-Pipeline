import pickle as pk
import yaml
from os import getcwd
from pathlib import Path
from mimetypes import guess_type
from datetime import datetime as dt


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

    def __init__(self,
                 args_d: dict = {}):
        """
        Initializing
        """
        if not isinstance(args_d, dict):
            raise ValueError(f"args_d must be a non-empty dictionary. {str(type(args_d))} is given.")
        if not args_d:
            args_d = self.default_args
        for key, val in args_d.items():
            if key in self.default_args \
                    and isinstance(val, type(self.default_args.get(key, None))):
                setattr(self, key, val)


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


class PreprocessingArgs(ArgsParserUtil):
    prep_args_dict = {
        # merge_pairs
        "fastq_maxdiffs": (int, 50),
        "fastq_pctid": (int, 50),
        "fastq_minmergelen": (int, 200),
        "fastq_maxmergelen": (int, 600),
        # trim_both_sides
        "stripleft": (int, 5),
        "stripright": (int, 5),
        # trim_one_side
        "stripleft": (int, 5),
        # filter_merged:
        "fastq_maxee_rate": (float, 0.002),
        # filter_one_side
        "fastq_truncqual": (int, 10),
        # "fastq_trunclen": "lambda x: mean(x) - 0.1 * mean(x) - 5"
        # dereplication
        "sizein": (bool, True),
        "sizeout": (bool, True),
        # sort_sequences:
        # clusterZOTUs
        "minsize": 2,
        # filter16S
        "e": (float, 0.1),
        "num_alignments": (int, 1),
        # prepare_zotus:
        # build_ZOTU_table:
        "id": (float, 0.97),
        # filter_zotu_abundance
        "abundance_cutoff": (float, 0.0),
        # select_zotu_seqs
        # addTax
        "turn": (str, "all"),
        # create_final_ZOTU_table:
    }


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
