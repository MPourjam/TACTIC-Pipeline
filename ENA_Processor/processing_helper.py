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


class YamlArgs:
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


    def __init__(self,
                 preproc_arg_filepath: str = None,
                 analysis_arg_filepath: str = None,
                 args_dict: dict = None):
        """
        Cool Class :)
        """
        pass
        
