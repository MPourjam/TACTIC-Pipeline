import pickle as pk
from os import getcwd
from pathlib import Path
from mimetypes import guess_type
from datetime import datetime as dt
from collections import OrderedDict as OD


class task_pickle:
    path = getcwd()
    scode_d = {"Done": 0,
               "Queue": 1,
               "Started": 2,
               "Progress": 3,
               "Error": -1
               }
    args_dict = OD({"input_dir": "",
                    "paired": "",
                    "forward_file": "",
                    "reverse_file": "",
                    "input_id": "",
                    "spike_amount": 6
                    })
    status_dict = {"download": scode_d["Queue"], "run": scode_d["Queue"]}
    task_dict = {"args": args_dict, "status": status_dict}

    def check_args_dict(self, args_dict):
        if not isinstance(args_dict, dict):
            return False
        keys = ["input_dir",
                "paired",
                "forward_file",
                "reverse_file",
                "input_id",
                "spike_amount"]
        has_key_list = [True for k in keys if k not in args_dict]
        if any(has_key_list):
            return False
        return True

    def check_status_dict(self, status_dict):
        if not isinstance(status_dict, dict):
            return False
        keys = ["download", "run"]
        has_key_list = [True for k in keys if k not in status_dict]
        if any(has_key_list):
            return False
        return True

    def set_task_dict(self, task_dict):
        if isinstance(task_dict, dict):
            keys = ["args", "status"]
            if any([True for k in keys if k not in task_dict]):
                msg = "Task dict must have dictionaries 'args' and 'status'."
                raise TypeError(msg)
            if not self.check_args_dict(task_dict["args"]):
                msg = "Argument dictionary is not formatted correclty."
                raise TypeError(msg)
            if not self.check_status_dict(task_dict["status"]):
                msg = "Status dictionary is not formatted correclty."
                raise TypeError(msg)
        else:
            msg = "Given object as task dictionary is not a dictionary."
            raise TypeError()
        self.task_dict = task_dict
        return True

    def __init__(self, filepath=None, task_dict=None):
        timestamp = dt.now().strftime('%Y%m%d_%H%M%S')
        path = Path(filepath) if filepath else Path(self.path)
        if path.is_dir():
            path = path.joinpath("{}_proc_task_d.pk".format(str(timestamp)))
        self.path = path.resolve()

        if self.path.exists():
            f_type, f_enc = guess_type(str(self.path))
            if "x-tex-pk" not in f_type:
                raise FileExistsError("File is not a pickle file.")
            with open(self.path, "rb") as pkfile:
                file_content = pk.load(pkfile)
            self.set_task_dict(file_content)

        if task_dict:
            self.set_task_dict(task_dict)

    def write(self):
        with open(self.path, 'wb') as pkfile:
            pk.dump(self.task_dict, pkfile, protocol=pk.HIGHEST_PROTOCOL)

    def args_complete(self):
        arg_d = self.task_dict["args"]
        input_d = arg_d["input_dir"]
        for k, v in arg_d.items():
            if k == "input_dir":
                try:
                    p = Path(v)
                except Exception as e:
                    return False
                if not p.exists():
                    return False
            elif k == "forward_file":
                p = Path(input_d)/str(v)
                if not p.exists():
                    return False
            elif k == "reverse_file":
                p = Path(input_d)/str(v)
                if not p.exists():
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
                except Exception as e:
                    return False
            else:  # Existence of any other arguments
                return False
        return True
