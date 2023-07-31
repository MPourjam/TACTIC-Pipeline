#!/usr/bin/env python
import re
import random
import argparse
import datetime
import pickle as pk
from time import sleep
from pathlib import Path
# from pprint import pprint as print
from shutil import rmtree
from statistics import mean, stdev
from processing_helper import TaskPickle
from multiprocessing import cpu_count
from task_caller import main as process_multi
from multiprocessing.pool import Pool
from processing_job import gzip_to_fastq

max_pool = int(cpu_count() * 0.4)
pool_size = max_pool if max_pool > 0 else 1
forw_file_match = re.compile(r"(_R1_|_1\.|_F)?")


def proper_length(forw_fastqfile):
    n_first_seqs = 10000
    assess_coeff = 0.1
    min_len = 150
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
    # Upper Bound of length mean of min_len, It's not possible to have varying lenght but only to get assured.
    len_limit = min_len + int(2.33 * len_stdev)
    if len_mean < len_limit:
        return False
    else:
        return True


def proper_pkarg_o(task_pk_path):
    '''
    It gets a the path to a TaskPickle file (file of TaskPickle obj)
    and checks the forward and reverse file. Then updated the pickle
    file of the dataset.
    '''
    task_path = Path(task_pk_path)
    if not task_path.is_file():
        raise TypeError("Input should be path to a TaskPickle object.")
    try:
        pkfo = TaskPickle(str(task_path))
    except Exception as e:
        print(e)
        return False
    dir_path = Path(pkfo.task_dict["args"]["input_dir"])
    down_status = pkfo.task_dict["status"]["download"]
    run_status = pkfo.task_dict["status"]["run"]
    # Excluding finished and failed samples from process list. Keeping the status unchanged
    if run_status == pkfo.scode_d["Done"]:
        pk_name = Path(pkfo.path).name
        print("Processed Successfully: {}".format(str(pk_name)))
        return False
    # Excluding the failed runs to save the message. When run fails the fastq gets deleted and then message becomes {line:113}
    if run_status == pkfo.scode_d["Error"]:
        pk_name = Path(pkfo.path).name
        print("Process has failed: {}".format(str(pk_name)))
        return False
    if down_status != pkfo.scode_d["Done"]:
        print("Download not finished: {}".format(dir_path.name))
        return False
    # Checking existence of the directory
    if not dir_path.is_dir():
        pkfo.task_dict["status"]["run"] = pkfo.scode_d["Error"]
        msg = "Input_dir should be path to a directory."
        pkfo.task_dict["status"]["msg"] = msg
        pkfo.write()
        print("Not a directory: {}".format(dir_path))
        return False
    fastqs_files = []
    f_path = dir_path.joinpath(pkfo.task_dict["args"]["forward_file"])
    r_path = dir_path.joinpath(pkfo.task_dict["args"]["reverse_file"])
    # Checking existence of forward file
    if not f_path.is_file():
        pkfo.task_dict["args"]["forward_file"] = ""
        pkfo.task_dict["status"]["download"] = pkfo.scode_d["Queue"]
        print("File does not exist: {}".format(f_path.name))
    else:
        fastqs_files.append(f_path)
    # Checking existence of reverse file
    paired = True if "yes" == str(pkfo.task_dict["args"]["paired"]).lower() else False
    if paired and not r_path.is_file():  # TODO Are we sure that the download is in progress???
        pkfo.task_dict["status"]["download"] = pkfo.scode_d["Queue"]
        print("File does not exist: {}".format(r_path.name))
    elif not paired:
        pkfo.task_dict["args"]["reverse_file"] = ""
        # If not paired then r_path is a directory path (dir_path)
    if r_path.is_file():
        fastqs_files.append(r_path)
    if not fastqs_files:
        pkfo.task_dict["status"]["msg"] = "No fastq file."
        pkfo.write()
        return False
    # Checking if the file can get opened
    for fi_path in fastqs_files:
        write = 0
        while write < 3:
            if fi_path.is_file():
                try:
                    fi_cont = open(fi_path)
                    fi_cont.close()
                    write = 3
                except Exception as e:
                    if write > 1:
                        print("{}".format(str(e)))
                        return False
                    sleep(3)
            write += 1
    # Demultiplexing
    if len(fastqs_files) == 1 and paired:
        # TODO demultiplexing
        pkfo.task_dict["status"]["download"] = pkfo.scode_d["Error"]
        msg = "Fastq file needs to get split!!!"
        pkfo.task_dict["status"]["msg"] = msg
        pkfo.write()
        fastq_file = pkfo.task_dict["args"]["forward_file"]
        print("Demultiplexing needed for: {}".format(fastq_file))
        return False
        # Get the updated fastq files names
    # TODO Checking if the forward and reverse file are correctly assigned.
    pkfo.task_dict["args"]["forward_file"] = f_path.name if f_path else ""
    pkfo.task_dict["args"]["reverse_file"] = r_path.name if r_path else ""
    if not pkfo.args_complete():
        pkfo.task_dict["status"]["run"] = pkfo.scode_d["Queue"]
        pkfo.task_dict["status"]["msg"] = "Argument(s) missing."
        pkfo.write()
        print("Incomplete argument(s) for {}".format(pkfo))
        return False
    pkfo.task_dict["status"]["run"] = pkfo.scode_d["Queue"]
    # Making sure the download status is correctly saved
    pkfo.task_dict["status"]["download"] = pkfo.scode_d["Done"]
    pkfo.task_dict["status"]["msg"] = ""
    pkfo.write()
    return pkfo


def decompress(task_pko):
    dir_path = Path(task_pko.task_dict["args"]["input_dir"])
    forw_name = task_pko.task_dict["args"]["forward_file"]
    reve_name = task_pko.task_dict["args"]["reverse_file"]
    f_path = dir_path.joinpath(forw_name)
    r_path = dir_path.joinpath(reve_name)
    f_path = gzip_to_fastq(str(f_path)) if f_path.is_file() else ""
    r_path = gzip_to_fastq(str(r_path)) if r_path.is_file() else ""
    f_path = Path(f_path[0]).name if f_path else ""
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
    with Pool(pool_size) as pool:
        pks_list = list(pks)
        pfiles_res_list = pool.map(proper_pkarg_o, pks_list)
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
                    pk_o.task_dict["status"]["run"] = pk_o.scode_d["Error"]
                    pk_o.task_dict["status"]["msg"] = "Improper Sequence Length"
                    pk_o.task_dict["args"]["forward_file"] = ""
                    pk_o.task_dict["args"]["reverse_file"] = ""
                    pk_o.write()
                    pk_name = str(pk_o.task_dict["args"]["input_id"]) + ".pk"
                    input_dir_path = Path(pk_o.task_dict["args"]["input_dir"])
                    for fi in input_dir_path.iterdir():
                        if fi.is_file() and str(fi.name) != pk_name:
                            fi.unlink()
                        elif fi.is_dir():
                            rmtree(str(fi))
                    err_pks.append(pk_o.path)
                    continue
                final_ppks.append(pk_o.path)

    if res_pk:
        with open(pickle_file, 'wb') as dest:
            pk.dump(final_ppks, dest, protocol=pk.HIGHEST_PROTOCOL)
        print("Task pickles were written to: {}".format(pickle_file))
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
        done_pkos = []
        for err_pk in err_pks:
            pko = TaskPickle(err_pk)
            if pko.task_dict["status"]["run"] == pko.scode_d["Done"]:
                done_pkos.append(err_pk)
        for dpkos in done_pkos:
            err_pks.remove(dpkos)
        good_pks = "\t" + "\n\t".join(pickle_files_path)
        done_pks = "\t" + "\n\t".join(done_pkos)
        improper_pks = "\t" + "\n\t".join(err_pks)
        msg_prop = "Proper Task Pickles:\n{}".format(good_pks)
        msg_done = "Finished Task Pickles:\n{}".format(done_pks)
        msg_inprop = "Improper Task Pickles:\n{}".format(improper_pks)
        msg = str(msg_prop) + "\n" + str(msg_done) + "\n" + str(msg_inprop)
        print(msg)
