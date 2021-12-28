#! /usr/bin/python3.7
from processing_job import main_processing
import pickle as pk
from multiprocessing import Process
import os

file_path = "/srv/cfm/inputs/proc_dicts_20211227_1442.pk"

with open(file_path, 'rb') as f:
    dicts_list = pk.load(f)

process_args_list = dicts_list[3:]
p_list = []
for argset in process_args_list:
    report_path = os.path.join(argset['input_dir'], "report.txt")
    if not os.path.isfile(report_path):
        p_list.append(Process(target=main_processing, args=tuple(argset.values())))

for p in p_list:
    p.start()
    print("Process {} started.".format(p))

for p in p_list:
    print("Process {} joined.".format(p))
    p.join()
