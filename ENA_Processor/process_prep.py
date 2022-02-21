#!/usr/bin/env python
import re
import random
import argparse
import datetime
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
forw_file_match = re.compile(r"(_R1_|_1\.|_F)?")


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


def preparation_main(dir_path, res_pk=False):
    '''
    if res_pk is false then it returns a lis containing the paths to
    run pickle files.
    '''
    fastqs_dir = Path(dir_path)
    if fastqs_dir.is_file():
        exit("Given path as directory is actually a file!!!")
    timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M')
    pickle_file = Path(__file__).parent / "{}_tasks_to_call.pk".format(timestamp)
    pks = fastqs_dir.rglob("./**/*.pk*")
    parents = OD()
    # Getting only those that are finished with download.
    for pkfile in pks:
        try:
            pkf_obj = task_pickle(pkfile)
        except Exception:
            continue
        if pkf_obj.task_dict["status"]["download"] != pkf_obj.scode_d["Done"]:
            continue
        p = Path(pkf_obj.task_dict["args"]["input_dir"])
        '''
        # There is no need to add forward and reverse as it is already added.
        files = parents.get(p, [])
        fqfiles = p.glob("*.fastq*")
        if not fqfiles:
            continue
        fqfiles = sorted(fqfiles)
        files.extend(fqfiles)
        files = sorted(list(set(files)))
        if len(files) == 2:
            assert(re.search(forw_file_match, str(files[0]))), "Incorrect order of forward and reverse files."
        '''
        parents[p] = pkf_obj

    processing_dicts = []
    with ThreadPool(pool_size) as pool:  # TODO use the pool
        for p, pkfo in parents.items():
            dir_path = Path(pkfo.task_dict["args"]["input_dir"])
            fastqs_files = []
            f_path = dir_path.joinpath(pkfo.task_dict["args"]["forward_file"])
            r_path = dir_path.joinpath(pkfo.task_dict["args"]["reverse_file"])
            if not f_path.is_file():
                pkfo.task_dict["args"]["forward_file"] = ""
                pkfo.task_dict["status"]["download"] = pkfo.scode_d["Error"]
                pkfo.task_dict["status"]["msg"] = "Forward file does not exist."
                pkfo.write()
                continue
            fastqs_files.append(f_path)
            paired = True if pkfo.task_dict["args"]["paired"] == "Yes" else False
            if paired and not r_path.is_file():
                pkfo.task_dict["status"]["download"] = pkfo.scode_d["Error"]
                pkfo.task_dict["status"]["msg"] = "Reverse file does not exist."
                pkfo.write()
                continue
            elif not paired:
                pkfo.task_dict["args"]["reverse_file"] = ""
                r_path = ""
            pkfo.write()
            if r_path:
                fastqs_files.append(r_path)
            # Checking if the fastq file has proper length
            prop_len = [proper_length(fi) for fi in fastqs_files]
            if prop_len and not all(prop_len):
                pkfo.task_dict["status"]["run"] = pkfo.scode_d["Error"]
                pkfo.task_dict["status"]["msg"] = "Not a proper length for fastq file(s)."
                pkfo.write()
                continue
    # Demultiplexing
            acc_ = p.stem
            if len(fastqs_files) == 1 and paired:
                # TODO demultiplexing
                pkfo.task_dict["status"]["download"] = pkfo.scode_d["Error"]
                msg = "Fastq file in {} need to get split!!!".format(p)
                pkfo.task_dict["status"]["msg"] = msg
                continue
                # Get the updated fastq files names
            # TODO Checking if the forward and reverse file are correctly assigned.
            pkfo.task_dict["args"]["forward_file"] = parents[p][0].name
            pkfo.task_dict["args"]["reverse_file"] = ''
            if pkfo.task_dict["args"]["paired"] == "Yes":
                pkfo.task_dict["args"]["reverse_file"] = parents[p][1].name
            pkfo.task_dict["status"]["run"] = pkfo.scode_d["Started"]
            pkfo.write()
            if pkfo.args_complete():
                processing_dicts.append(pkfo.task_dict)
            else:
                print("Incomplete task arguments in {}".format(p))
    if res_pk:
        processing_args_list = [pk_c["args"] for pk_c in processing_dicts]
        with open(pickle_file, 'wb') as dest:
            pk.dump(processing_args_list, dest, protocol=pk.HIGHEST_PROTOCOL)
        return pickle_file
    else:
        return processing_dicts


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
        process_multi(pickle_file_path, res_pk=True)
    else:
        print(str(pickle_file_path))
