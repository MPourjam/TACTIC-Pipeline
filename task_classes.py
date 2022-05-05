import pickle as pk
from os import getcwd, remove
from os import path as ospath
from mimetypes import guess_type
from datetime import datetime as dt
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
        return os.path.getsize(os.path.realpath(f.name))

    def lock_file(f):
        msvcrt.locking(f.fileno(), msvcrt.LK_RLCK, file_size(f))

    def unlock_file(f):
        msvcrt.locking(f.fileno(), msvcrt.LK_UNLCK, file_size(f))


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
                 "spike_amount": 6
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
            f_type, f_enc = guess_type(str(self.__path))
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
                    v = int(v)
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
