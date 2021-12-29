#! /usr/bin/python3
from processing_job import main_processing
import pickle as pk
from multiprocessing import Process, cpu_count, Pool
import os
import time
import shutil
import re
from datetime import datetime as dt

def till_now(datetime_obj):
    time_delta = datetime_obj - dt.now()
    return time_delta.total_seconds()

def main(pickle_files_path=None):
    max_pool = int(cpu_count() * 0.5)
    pool_size = max_pool if max_pool > 0 else 1
    max_batch = int(pool_size * 0.5)
    batch_size = max_batch if max_batch > 0 else 1
    pickle_files_path = os.path.abspath(pickle_files_path) if pickle_files_path else "/srv/cfm/inputs/"
    if os.path.isdir(pickle_files_path):
        dir_path = pickle_files_path
        pk_files = os.scandir(dir_path)
        file_names = [e.name for e in pk_files if (e.is_file() and e.name.endswith(".pk"))]
        times_dict = {}
        for fn in file_names:
            t_match = re.search("[0-9]{8}_[0-9]{4}", fn)
            if t_match:
                time_delta = till_now(dt.strptime(str(fn[t_match.start():t_match.end()]), '%Y%m%d_%H%M'))
                times_dict[fn] = time_delta
        if times_dict:
            latest_pk_file = sorted(list(times_dict.items()), key=lambda x: x[1], reverse=True)[0][0]
        else:
            raise FileNotFoundError("No pickle file in {}.".format(dir_path))
        main(pickle_files_path=os.path.join("/srv/cfm/inputs/", latest_pk_file))
    elif os.path.isfile(pickle_files_path):
        file_path = pickle_files_path
        try:
            with open(file_path, 'rb') as f:
                dicts_list = pk.load(f)
        except Exception as e:
            return e
        process_args_list = dicts_list[::]
        p_list = []
        for argset in process_args_list:
            report_path = os.path.join(argset['input_dir'], "report.txt")
            if not os.path.isfile(report_path):
                p_list.append(tuple(argset.values()))
        print(p_list)
        with Pool(pool_size, maxtasksperchild=1) as pool:
            l = len(p_list)
            m = l // batch_size
            r = l - (m * batch_size)
            if m > 0:
                for i in range(1, m + 1):
                    print("Multiple of {}".format(batch_size))
                    s = (i - 1) * batch_size
                    e = i * batch_size
                    res_list = [pool.apply_async(main_processing, args=argset) for argset in p_list[s:e]]
                    for res in res_list:
                        res.wait()

                if r > 0:
                    print("Doing the remaining {}.".format(r))
                    s = (m * batch_size) 
                    res_list = [pool.apply_async(main_processing, args=argset) for argset in p_list[s:]]
                    for res in res_list:
                        res.wait()

            elif r > 0:
                print("Less than {} jobs.".format(batch_size))
                res_list = [pool.apply_async(main_processing, args=argset) for argset in p_list]
                for res in res_list:
                        res.wait()
            else:
                print("All samples are processed in {} file.".format(file_path))
                print("Deleting {}.".format(file_path))
                os.remove(file_path)
                
    else:
        msg = "The {} should be either a pickle file or a directory containing the pickle files.".format(pickle_files_path)
        raise FileNotFoundError(msg)
