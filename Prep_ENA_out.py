#!/usr/bin/python3
import os
import re
import argparse
from pprint import pprint
import datetime
import mimetypes
import shutil
from random import randint
import pickle as pk
from collections import OrderedDict as OD
from task_caller import main as process_multi
from pathlib import Path
from statistics import mean, stdev
from multiprocessing import cpu_count
from multiprocessing.pool import ThreadPool

max_pool = int(cpu_count() * 0.4)
pool_size = max_pool if max_pool > 0 else 1
forw_file_match = re.compile("(_R1_|_1\.|_F)")

def proper_length(forw_fastqfile):
    # import re
    # seq_match = re.compile(r'[^ACTGN]')
    to_assess = [randint(1,10000) for i in range(1, 5000)]
    file_path = Path(forw_fastqfile)
    seq_counter = 0
    random_length = []
    with open(file_path, 'r+') as fastqfile:
        line = fastqfile.readline()[:-1]
        if not line.startswith("@"):
            raise TypeError("Faulty fastq file: {}".format(file_path))

        while line and seq_counter < 10000:
            #NOTE This is only to get quality line or seq line and assess the length
            #NOTE (even some quality lines might be discarded if they start with @ or + signs.)
            if line.startswith("@") or line.startswith("+"):
                # prviously_matched = False
                line = fastqfile.readline()[:-1]
            else:
                seq_counter += 1
                if seq_counter in to_assess:
                    length = len(line[:-1])
                    random_length.append(length)
                line = fastqfile.readline()[:-1]
    len_mean = mean(random_length)
    len_stdev = int(stdev(random_length))
    if len_mean < 150 and len_stdev > 1:
        return False
    else:
        return True


def preparation_main(dir_path):
    fastqs_dir = Path(args.directory).resolve()
    if fastqs_dir.is_file():
        exit("Given path as directory is a file!!!")
    timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M')
    pickle_file = fastqs_dir / "{}.pk".format(timestamp)
    fastqs = fastqs_dir.rglob("./**/*.fastq*")
    parents = OD()
    for fq in fastqs:
        # print(fq)
        p = fq.parent
        files = parents.get(p, [])
        files.append(fq.name)
        files = sorted(files)
        if len(files) == 2:
            assert(re.search(forw_file_match, files[0])), "Incorrect order of forward and reverse files."
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
                if typ == "application/zip":
                    os.system("unzip {}".format(full_path))

    with ThreadPool(pool_size) as pool:
        proper_length_dict = {p: proper_length(p.joinpath(fs[0].replace(".gz", ""))) for p, fs in parents.items()}

    # print(proper_length_dict)
    for p, proper in proper_length_dict.items():
        if proper:
            new_dict = OD()
            paired = "Yes" if len(parents[p]) == 2 else "No"
            new_dict["input_dir"] = str(p.resolve())
            new_dict["paired"] = paired
            new_dict["forward_file"] = parents[p][0].replace(".gz", "")
            new_dict["reverse_file"] = parents[p][-1].replace(".gz", "") if paired == "Yes" else ''
            new_dict["input_id"] = p.stem
            new_dict["spike_amount"] = 6
            processing_dicts.append(new_dict)
    with open(pickle_file, 'wb') as dest:
        pk.dump(processing_dicts, dest, protocol=pk.HIGHEST_PROTOCOL)
    return pickle_file

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory",
                        help="The directory to find recursively all fastq filese inside",
                        default=".")
    args = parser.parse_args()
    fastqs_dir = Path(args.directory).resolve()
    pickle_file_path = preparation_main(fastqs_dir)
    print(pickle_file_path)
    process_multi(pickle_file_path)