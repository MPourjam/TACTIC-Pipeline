#! /usr/bin/env python3
import random
import sys
from download_ena import download_main
from process_prep import preparation_main
from task_caller import main as main_caller
from processing_helper import TaskPickle
# from pprint import pprint as print
from pathlib import Path


met_file_all = "/srv/cfm/docs/extracted_metadata.tab"
output_dir_path = "/srv/cfm/outputs/"
ena_out = output_dir_path + "ENA_out/"
one_file = ena_out + "PRJNA400115/9606/SAMN09410203/SRR7582793"
one_file += "/SRR7582793.pk"

argv_1 = str(sys.argv[1]) if len(sys.argv) > 1 else str("00010")
argv_2 = int(sys.argv[2]) if len(sys.argv) > 2 else 0
if len(argv_1) < 5:
    argv_1 = argv_1.zfill(5)
argv_1 = argv_1[:5]
argv_list = [int(i) for i in argv_1]
rand_size = argv_2 if argv_2 else 0  # If 0 then no random sample is taken from the met_fiel_all
seed_val = 250  # If 0 then no seed is set in random sampling
chunk_s = 18 if rand_size == 0 else 0  # If 0 then the whole metadata file gets downladed
num_df = 1 if rand_size == 0 else 2
act = {
    "down": argv_list[0],
    "prep": argv_list[1],
    "proc": argv_list[2],
    "repo": argv_list[3],
    "allp": argv_list[4],
}


def take_random(metadat_path, app_size: int, random_size: int = 10, seed: int = 0):
        """
        Takes the path to metadata file and takes random sample from with random_size
        """
        if random_size == 0:
            return metadat_path
        print(f"Selecting random sample of size {random_size} from {metadat_path}.")
        if seed != 0:
                random.seed(seed)
        met_file = Path(metadat_path)
        met_file_size = app_size
        if met_file.is_dir():
                raise TypeError("met_file must be a tsv file.")
        met_file_o = met_file.parent.joinpath(met_file.stem + f"_rand{str(random_size)}.tab")
        select_ind = random.sample(range(1, met_file_size), k=random_size)
        select_ind = sorted(select_ind)
        print(f"min: {min(select_ind)}\tmax: {max(select_ind)}")
        # exit()
        with open(met_file, 'r') as metfile:
            counter = 0
            with open(met_file_o, "w") as metfile_o:
                # Header
                line = metfile.readline()
                metfile_o.write(line)
                line = metfile.readline()
            with open(met_file_o, "a") as metfile_o:
                # Rest
                for ind in select_ind:
                    while counter != ind:
                        line = metfile.readline()
                        counter += 1
                        # print(counter)
                    metfile_o.write(line)
        print(f"New sample with size of {random_size} saved to: {met_file_o}")
        return str(met_file_o)


met_file = take_random(met_file_all, 1400000, random_size=rand_size, seed=seed_val)
if act["down"]:
    print(download_main(met_file,
          outdir=output_dir_path,
          chunk_size=chunk_s,
          ndf=num_df,
          redownload_faileds=True))
if act["prep"]:
    ppks_path, err_ppks = preparation_main(ena_out)
    err_ms = ""
    ready_ms = ""
    other_ms = ""
    succe_ms = ""
    for ppk in ppks_path:
        pko = TaskPickle(ppk)
        in_dir = pko.task_dict["args"]["input_dir"]
        forw_fi = pko.task_dict["args"]["forward_file"]
        run_st = pko.task_dict["status"]["run"]
        pri_li = [in_dir.ljust(75), forw_fi.ljust(25), run_st]
        pri_li = [str(el) for el in pri_li]
        ready_ms += "\t" + "\t".join(pri_li) + "\n"
    for err_p in err_ppks:
            pko = TaskPickle(err_p)
            in_dir = pko.task_dict["args"]["input_dir"]
            forw_fi = pko.task_dict["args"]["forward_file"]
            msg = pko.task_dict["status"]["msg"]
            run_st = pko.task_dict["status"]["run"]
            msg_li = [in_dir.ljust(75), forw_fi.ljust(25), str(run_st).rjust(2), msg]
            msg_li = [str(el) for el in msg_li]
            if run_st == pko.scode_d["Error"]:
                    err_ms += "\t" + "\t".join(msg_li) + "\n"
            elif run_st == pko.scode_d["Done"]:
                    succe_ms += "\t" + "\t".join(msg_li) + "\n"
            else:
                    other_ms += "\t" + "\t".join(msg_li) + "\n"
    fin_msg = "Ready Tasks:\n"
    fin_msg += ready_ms + "\n"
    fin_msg += "Failed Tasks:\n"
    fin_msg += err_ms + "\n"
    fin_msg += "Successful Tasks:\n"
    fin_msg += succe_ms + "\n"
    fin_msg += "Other Tasks:\n"
    fin_msg += other_ms + "\n"
    print(fin_msg)
if act["proc"]:
    pk_path = preparation_main(ena_out, res_pk=True)
    main_caller(pk_path)

if act["repo"]:
    pks = Path(ena_out).rglob("./**/*.pk")
    pkolist = [TaskPickle(str(tpkf)) for tpkf in pks]
    for pko in pkolist:
        in_dir = pko.task_dict["args"]["input_dir"]
        if str(Path(pko.path).parent) != in_dir:
            print(f"ERRROOOOORRR {pko}")
        for_fi = pko.task_dict["args"]["forward_file"]
        pk_msg = pko.task_dict["status"]["msg"]
        spike = pko.task_dict["args"]["spike_amount"]
        run_st = str(pko.task_dict["status"]["run"]).rjust(2)
        down_st = str(pko.task_dict["status"]["download"]).rjust(2)
        to_write = [in_dir.ljust(75),
                    for_fi.ljust(25),
                    spike, down_st,
                    run_st, pk_msg]
        to_write = [str(e) for e in to_write]
        msg = '\t'.join(to_write)
        print(msg)

if act["allp"]:
    pko = TaskPickle(one_file)
    # pko.task_dict["status"]["run"] = pko.scode_d["Queue"]
    # pko.write()
    # print(bool(pko))
    print(pko.__dict__)
