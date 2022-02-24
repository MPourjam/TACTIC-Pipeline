#!/usr/bin/env python
import re
import random
import argparse
import datetime
import pickle as pk
from pathlib import Path
# from pprint import pprint as print
from statistics import mean, stdev
from task_classes import task_pickle
from multiprocessing import cpu_count
from task_caller import main as process_multi
from multiprocessing.pool import ThreadPool
from processing_job import gzip_to_fastq

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


def proper_pkarg_o(task_pk_path):
    '''
    It gets a the path to a task_pickle file (file of task_pickle obj)
    and checks the forward and reverse file. Then updated the pickle
    file of the dataset.
    '''
    task_path = Path(task_pk_path)
    if not task_path.is_file():
        raise TypeError("Input should be path to a task_pickle object.")
    try:
        pkfo = task_pickle(str(task_path))
    except Exception as e:
        print(e)
        return False
    if pkfo.task_dict["status"]["download"] != pkfo.scode_d["Done"]:
        return False
    dir_path = Path(pkfo.task_dict["args"]["input_dir"])
    if not dir_path.is_dir():
        pkfo.task_dict["status"]["run"] = pkfo.scode_d["Error"]
        msg = "Input_dir should be path to a directory."
        pkfo.task_dict["status"]["msg"] = msg
        pkfo.write()
        return False
    fastqs_files = []
    f_path = dir_path.joinpath(pkfo.task_dict["args"]["forward_file"])
    r_path = dir_path.joinpath(pkfo.task_dict["args"]["reverse_file"])
    # Checking existence of forward file
    if not f_path.is_file():
        pkfo.task_dict["args"]["forward_file"] = ""
        pkfo.task_dict["status"]["download"] = pkfo.scode_d["Error"]
        pkfo.task_dict["status"]["msg"] = "Forward file does not exist."
        pkfo.write()
        return False
    fastqs_files.append(f_path)
    # Checking existence of reverse file
    paired = True if pkfo.task_dict["args"]["paired"] == "Yes" else False
    if paired and not r_path.is_file():
        pkfo.task_dict["status"]["download"] = pkfo.scode_d["Error"]
        pkfo.task_dict["status"]["msg"] = "Reverse file does not exist."
        pkfo.write()
        return False
    elif not paired:
        pkfo.task_dict["args"]["reverse_file"] = ""
        r_path = ""
    if r_path:
        fastqs_files.append(r_path)
    # Checking if the file can get opened
    for fi_path in fastqs_files:
        if fi_path:
            try:
                fi_cont = open(fi_path)
                fi_cont.close()
            except Exception:
                return False
    # Demultiplexing
    if len(fastqs_files) == 1 and paired:
        # TODO demultiplexing
        pkfo.task_dict["status"]["download"] = pkfo.scode_d["Error"]
        msg = "Fastq file in {} need to get split!!!".format(pkfo)
        pkfo.task_dict["status"]["msg"] = msg
        pkfo.write()
        print(pkfo.task_dict)
        return False
        # Get the updated fastq files names
    # TODO Checking if the forward and reverse file are correctly assigned.
    pkfo.task_dict["args"]["forward_file"] = f_path.name if f_path else ""
    pkfo.task_dict["args"]["reverse_file"] = r_path.name if r_path else ""
    if not pkfo.args_complete():
        print("Incomplete argument(s) for {}".format(pkfo))
        pkfo.task_dict["status"]["run"] = pkfo.scode_d["Error"]
        pkfo.task_dict["status"]["msg"] = "Argument(s) missing."
        pkfo.write()
        return False
    pkfo.task_dict["status"]["run"] = pkfo.scode_d["Queue"]
    # Making sure the download status is correctly saved
    pkfo.task_dict["status"]["download"] = pkfo.scode_d["Done"]
    pkfo.write()
    return pkfo


def decompress(task_pko):
    dir_path = Path(task_pko.task_dict["args"]["input_dir"])
    f_path = dir_path.joinpath(task_pko.task_dict["args"]["forward_file"])
    f_path = gzip_to_fastq(f_path) if f_path else ""
    f_path = Path(f_path[0]).name if f_path else ""
    r_path = dir_path.joinpath(task_pko.task_dict["args"]["reverse_file"])
    r_path = gzip_to_fastq(r_path) if r_path else ""
    r_path = Path(r_path[0]).name if r_path else ""
    task_pko.task_dict["args"]["forward_file"] = f_path
    task_pko.task_dict["args"]["reverse_file"] = r_path
    task_pko.write()
    return task_pko


def preparation_main(dir_path, res_pk=False):
    '''
    if res_pk is false then it returns a lis containing the paths to
    run pickle files.
    '''
    fastqs_dir = Path(dir_path)
    if fastqs_dir.is_file():
        exit("Given path as directory is actually a file!!!")
    timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M')
    pickle_file = Path(__file__).parent/"{}_tasks_to_call.pk".format(timestamp)
    pks = fastqs_dir.rglob("./**/*.pk")
    # Getting only those that are finished with download.
    with ThreadPool(pool_size) as pool:
            pks_list = list(pks)
            pfiles_res_list = pool.map(proper_pkarg_o, pks_list)
            print(pfiles_res_list)
            ppks_list = []
            err_pks = []
            for pk_path, pk_o in zip(pks_list, pfiles_res_list):
                if not pk_o:
                    err_pks.append(str(pk_path))
                    continue
                ppks_list.append(pk_o)
            ppks_forw_files = []
            # TODO Check if waits for the resutls
            ppks_list = pool.map(decompress, ppks_list)
            for pk_o in ppks_list:
                dir_path = Path(pk_o.task_dict["args"]["input_dir"])
                f_path = dir_path.joinpath(pk_o.task_dict["args"]["forward_file"])
                ppks_forw_files.append(str(f_path))
            final_ppks = []
            if ppks_forw_files:
                plen_res_list = pool.map(proper_length, ppks_forw_files)
                for i in range(0, len(plen_res_list)):
                    pk_o = ppks_list[i]
                    if not plen_res_list[i]:
                        err_pks.append(pk_o.path)
                        continue
                    final_ppks.append(pk_o.path)

    if res_pk:
        with open(pickle_file, 'wb') as dest:
            pk.dump(final_ppks, dest, protocol=pk.HIGHEST_PROTOCOL)
        print(pickle_file)
        return pickle_file
    else:
        return (final_ppks, err_pks)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    help_text = "The directory to find recursively all process pk files inside"
    parser.add_argument("-d", "--directory",
                        help=help_text,
                        default=".")
    parser.add_argument("-nc", "--no-call",
                        help="If call the tasks from pickle files.",
                        default=False, action='store_true')
    args = parser.parse_args()
    fastqs_dir = Path(args.directory).resolve()
    if not args.no_call:
        pickle_files_path = preparation_main(fastqs_dir, res_pk=True)
        process_multi(pickle_files_path)
    else:
        pickle_files_path, err_pks = preparation_main(fastqs_dir)
        good_pks = "\t" + "\n\t".join(pickle_files_path)
        inproper_pks = "\t" + "\n\t".join(err_pks)
        msg_prop = "Proper Task Pickles:\n{}".format(good_pks)
        msg_inprop = "Inproper Task Pickles:\n{}".format(inproper_pks)
        msg = str(msg_prop) + "\n" + str(msg_inprop)
        print(msg)
