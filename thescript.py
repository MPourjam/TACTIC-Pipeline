#!/usr/bin/python3.7

import os
import re
import glob
from pprint import pprint
import datetime
from processing_job import main_processing
import shutil
import random
import pickle as pk

os.chdir("/srv/cfm/inputs/")
cwd = os.getcwd()
fastqs = glob.glob("/srv/cfm/inputs/**/*Wasser*.fastq*")
times = {}
for pd in fastqs:
    stat = os.stat(pd)
    time = datetime.timedelta(seconds=stat.st_ctime)
    times.update({os.path.split(pd)[1]: time})

times_sorted = sorted(times.items(), key=lambda x: x[1])
# pprint(fastqs)
# history_origin = times_sorted[0][1] - datetime.timedelta(days=1)

fors = sorted([fq for fq in fastqs if re.search("_R1_", os.path.basename(fq))])
backs = sorted([fq for fq in fastqs if re.search("_R2_", os.path.basename(fq))])

till = min(len(fors), len(backs))

processing_dicts = []

for i in range(0,till):
    # process_dict = dict()
    forw_file_dir, forw_file_name = os.path.split(fors[i])
    back_file_dir, back_file_name = os.path.split(backs[i])
    path_tail = str(fors[i]).split("inputs/")[-1]
    path_tail_dir = path_tail.split(".")[0]
    path_tail_dir = path_tail_dir.split("_R")[0]
    file_suffix = ".".join(path_tail.split(".")[1:])
    path_base = str(fors[i]).split("inputs/")[0]
    output_dir_path = "{}outputs/{}/".format(path_base, path_tail_dir)
    os.makedirs(output_dir_path, exist_ok=True)
    for_file_abspath = os.path.join(output_dir_path, forw_file_name)
    rev_file_abspath = os.path.join(output_dir_path, back_file_name)
    shutil.move(fors[i], for_file_abspath)
    shutil.move(backs[i], rev_file_abspath)
    for_file = os.path.basename(for_file_abspath)
    rev_file = os.path.basename(rev_file_abspath)
    input_id_match = re.search("[0-9]+", forw_file_name)
    if input_id_match:
        input_id = int(forw_file_name[input_id_match.start():input_id_match.end()])
    else:
        input_id = random.randint(0, 99)
    
    process_dict = {"input_dir": output_dir_path,
                    "paired": "Yes",
                    "forward_file": for_file,
                    "reverse_file": rev_file,
                    "input_id": input_id,
                    "spike_amount": 6}
    processing_dicts.append(process_dict)

timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M')

pickle_file = os.path.join("/srv/cfm/inputs/", "proc_dicts_{}.pk".format(timestamp))

with open(pickle_file, 'wb') as dest:
    pk.dump(processing_dicts, dest, protocol=pk.HIGHEST_PROTOCOL)

pprint(processing_dicts)

# main_processing()
