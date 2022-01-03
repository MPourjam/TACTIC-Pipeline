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
from collections import OrderedDict as OD
from task_caller import main as process_multi

os.chdir("/srv/cfm/inputs/")
cwd = os.getcwd()

def proper_length(forw_fastqfile):
    from random import randint
    from statistics import mean, stdev
    from pathlib import Path
    import re
    import time
    start = time.time()
    seq_match = re.compile(r"^[ATCGUN]+$")
    to_assess = [randint(1,1000) for i in range(1, 100)]
    file_path = Path(forw_fastqfile)
    seq_counter = 0
    random_length = []
    with open(file_path, 'r+') as fastqfile:
        line = fastqfile.readline()[:-1]
        if not line.startswith("@"):
            raise TypeError("Faulty fastq file: {}".format(file_path))
        while line:
            if line.startswith("@"):
                line = fastqfile.readline()
            elif line.startswith("+"):
                line = fastqfile.readline()
            elif re.search(seq_match, line):
                seq_counter += 1
                if seq_counter in to_assess:
                    random_length.append(len(line[:-1]))
                line = fastqfile.readline()
            else:
                line = fastqfile.readline()
    len_mean = mean(random_length)
    len_stdev = int(stdev(random_length))
    end = time.time()
    exe_time = end - start
    # print(len(random_length))
    if len_mean < 150 and len_stdev > 1:
        return (False, seq_counter, exe_time)
    else:
        return (True, seq_counter, exe_time)

fastqs = glob.glob("/srv/cfm/inputs/2021_12_27/DNAstab*.fastq*")
print(fastqs)
times = {}
for pd in fastqs:
    stat = os.stat(pd)
    time = datetime.timedelta(seconds=stat.st_ctime)
    times.update({os.path.split(pd)[1]: time})

times_sorted = sorted(times.items(), key=lambda x: x[1])
# pprint(fastqs)
# history_origin = times_sorted[0][1] - datetime.timedelta(days=1)

fors = sorted([fq for fq in fastqs if re.search("R1", os.path.basename(fq))])
backs = sorted([fq for fq in fastqs if re.search("R2", os.path.basename(fq))])

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
    shutil.copy(fors[i], for_file_abspath)
    shutil.copy(backs[i], rev_file_abspath)
    os.system("gunzip -f {}".format(for_file_abspath))
    os.system("gunzip -f {}".format(rev_file_abspath))
    for_file = os.path.basename(for_file_abspath)
    rev_file = os.path.basename(rev_file_abspath)
    for_file = for_file.replace(".gz", "")
    rev_file = rev_file.replace(".gz", "")
    input_id = path_tail_dir.split("/")[-1]
    process_dict = OD()
    process_dict["input_dir"] = output_dir_path
    process_dict["paired"] = "Yes"
    process_dict["forward_file"] = for_file
    process_dict["reverse_file"] = rev_file
    process_dict["input_id"] = input_id
    process_dict["spike_amount"] = 6

    processing_dicts.append(process_dict)

timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M')

pickle_file = os.path.join("/srv/cfm/inputs/", "proc_dicts_{}.pk".format(timestamp))

with open(pickle_file, 'wb') as dest:
    pk.dump(processing_dicts, dest, protocol=pk.HIGHEST_PROTOCOL)

process_multi(pickle_file)
