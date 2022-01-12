#!/usr/bin/env python
import os
import re
import random
import argparse
import datetime
import mimetypes
import pickle as pk
from pathlib import Path
from pprint import pprint as print
from statistics import mean, stdev
from task_classes import task_pickle
from multiprocessing import cpu_count
from collections import OrderedDict as OD
from task_caller import main as process_multi
from multiprocessing.pool import ThreadPool

max_pool = int(cpu_count() * 0.4)
pool_size = max_pool if max_pool > 0 else 1
forw_file_match = re.compile(r"(_R1_|_1\.|_F)")


def proper_length(forw_fastqfile):
    n_first_seqs = 10000
    assess_coeff = 0.5
    min_len = 150
    max_len_std = 1
    sample_range = round(assess_coeff * n_first_seqs)
    to_assess = random.sample(range(1, n_first_seqs), sample_range)
    file_path = Path(forw_fastqfile)
    seq_counter = 0
    random_length = []
    with open(file_path, 'r+') as fastqfile:
        line = fastqfile.readline()[:-1]
        lc = 0
        if not line.startswith("@"):
            raise TypeError("Faulty fastq file: {}".format(file_path))
        while line and seq_counter <= n_first_seqs:
            mod = lc % 4
            if mod == 1:
                seq_counter += 1
                if seq_counter in to_assess:
                    random_length.append(len(line))
            lc += 1
            line = fastqfile.readline()[:-1]
    len_mean = mean(random_length)
    len_stdev = int(stdev(random_length))
    if len_mean < min_len and len_stdev > max_len_std:
        return False
    else:
        return True


def preparation_main(dir_path):
    fastqs_dir = Path(dir_path)
    if fastqs_dir.is_file():
        exit("Given path as directory is actually a file!!!")
    timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M')
    pickle_file = Path(os.getcwd()) / "{}_tasks_to_call.pk".format(timestamp)
    pks = fastqs_dir.rglob("./**/*.pk*")
    parents = OD()
    for pkfile in pks:
        # print(fq)
        try:
            pkf_obj = task_pickle(pkfile)
        except Exception as e:
            continue
        if pkf_obj.task_dict["status"]["download"] != pkf_obj.scode_d["Done"]:
            continue
        p = Path(pkf_obj.task_dict["args"]["input_dir"])
        files = parents.get(p, [])
        fqfiles = p.glob("*.fastq*")
        if not fqfiles:
            continue
        fqfiles = sorted(fqfiles)
        files.extend(fqfiles)
        files = sorted(list(set(files)))
        if len(files) == 2:
            assert(re.search(forw_file_match, str(files[0]))), "Incorrect order of forward and reverse files."
        parents[p] = files

    processing_dicts = []
    with ThreadPool(pool_size) as p:
        for p, fs in parents.items():
            for f in fs:
                full_path = p.joinpath(f)
                assert(full_path.is_file()), "File {} does not exist.".format(full_path)
                typ, enc = mimetypes.guess_type(str(full_path))
                if enc == 'gzip':
                    os.system("gunzip -f {}".format(full_path))
                elif typ == "application/zip":
                    os.system("unzip {}".format(full_path))
            fqfiles = p.glob("*.fastq")
            fqfiles = sorted(list(set(fqfiles)))
            parents[p] = fqfiles

    with ThreadPool(pool_size) as pool:
        proper_length_dict = {p: proper_length(p.joinpath(fs[0])) for p, fs in parents.items() if fs}

    for p, proper in proper_length_dict.items():
        if proper:
            acc_ = p.stem
            pkfile = p.joinpath("{}.pk".format(acc_))
            pkf_obj = task_pickle(str(pkfile.resolve()))
            if len(parents[p]) == 1 and pkf_obj.task_dict["args"]["paired"] == "Yes":
                # TODO demultiplexing
                print("Fastq file in {} need to get split!!!".format(p))
                continue
                # Get the updated fastq files names
            if len(parents[p]) > 2:
                print("More than two fastq files in {}!!!".format(p))
                continue
            pkf_obj.task_dict["args"]["forward_file"] = parents[p][0].name
            pkf_obj.task_dict["args"]["reverse_file"] = ''
            if pkf_obj.task_dict["args"]["paired"] == "Yes":
                pkf_obj.task_dict["args"]["reverse_file"] = parents[p][1].name
            pkf_obj.task_dict["status"]["run"] = pkf_obj.scode_d["Started"]
            pkf_obj.write()
            if pkf_obj.args_complete():
                processing_dicts.append(pkf_obj.task_dict["args"])
            else:
                print("Incomplete task arguments in {}".format(p))

    with open(pickle_file, 'wb') as dest:
        pk.dump(processing_dicts, dest, protocol=pk.HIGHEST_PROTOCOL)
    return pickle_file


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory",
                        help="The directory to find recursively all process pk files inside",
                        default=".")
    parser.add_argument("-nc", "--no-call",
                        help="If call the tasks from pickle files.",
                        default=False, action='store_true')
    args = parser.parse_args()
    fastqs_dir = Path(args.directory).resolve()
    pickle_file_path = preparation_main(fastqs_dir)
    if not args.no_call:
        process_multi(pickle_file_path)
    else:
        print(str(pickle_file_path))
