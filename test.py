#! /usr/bin/python3
from processing_job import main_processing
import pickle as pk
from multiprocessing import Process, cpu_count, Pool
import os
import time

file_path = "/srv/cfm/inputs/proc_dicts_20211227_1442.pk"
pool_size = int(cpu_count() * 0.5)
batch_size = int(pool_size * 0.5)

with open(file_path, 'rb') as f:
    dicts_list = pk.load(f)

process_args_list = dicts_list[::]
p_list = []
for argset in process_args_list:
    report_path = os.path.join(argset['input_dir'], "report.txt")
    if not os.path.isfile(report_path):
        # p_list.append(Process(target=main_processing, args=tuple(argset.values())))
        p_list.append(tuple(argset.values()))

print(p_list)
print(batch_size)
with Pool(pool_size, maxtasksperchild=1) as pool:
    l = len(p_list)
    m = l // batch_size
    # print(m)
    r = l - (m * batch_size)
    print(r)
    if m > 0:
        for i in range(1, m + 1):
            print("Multiple of {}".format(batch_size))
            s = (i - 1) * batch_size
            e = i * batch_size
            res_list = [pool.apply_async(main_processing, args=argset) for argset in p_list[s:e]]
            for res in res_list:
                res.wait()
            # last_res = res_list[batch_size - 1]
            # last_res.wait()
        if r > 0:
            print("Doint the remaining {}.".format(r))
            s = m * batch_size 
            res_list = [pool.apply_async(main_processing, args=argset) for argset in p_list[s:]]
            for res in res_list:
                res.wait()
            # last_res = res_list[batch_size - 1]
            # last_res.wait()
    else:
        print("Less than {} jobs.".format(batch_size))
        res_list = [pool.apply_async(main_processing, args=argset) for argset in p_list]
        for res in res_list:
                res.wait()
