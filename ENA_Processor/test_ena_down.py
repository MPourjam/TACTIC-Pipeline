#! /usr/bin/env python3
from download_ena import download_main
from process_prep import preparation_main, proper_length, proper_pkarg_o
from task_caller import main as main_caller
from task_classes import TaskPickle
# from pprint import pprint as print
from pathlib import Path
from os import stat

met_file = "/srv/cfm/outputs/extracted_metadata.tab"
ena_out = "/srv/cfm/outputs/ENA_out/"
chunk_s = 10
act = {
       "down": 0,
       "prep": 0,
       "proc": 0,
       "repo": 1,
       "allp": 1,
       }

if act["down"]:
       print(download_main(met_file,
                     outdir="/srv/cfm/outputs/",
                     chunk_size=chunk_s,
                     ndf=1))
if act["prep"]:
       ppks_path, err_ppks = preparation_main(ena_out)
       err_ms = ""
       ready_ms = ""
       other_ms = ""
       for ppk in ppks_path:
              pko = TaskPickle(ppk)
              in_dir = pko.task_dict["args"]["input_dir"]
              forw_fi = pko.task_dict["args"]["forward_file"]
              run_st = pko.task_dict["status"]["run"]
              pri_li = [in_dir, forw_fi, run_st]
              pri_li = [str(el) for el in pri_li]
              ready_ms += "\t" + "\t".join(pri_li) + "\n"
       for err_p in err_ppks:
              pko = TaskPickle(err_p)
              in_dir = pko.task_dict["args"]["input_dir"]
              msg = pko.task_dict["status"]["msg"]
              run_st = pko.task_dict["status"]["run"]
              msg_li = [in_dir, run_st, msg]
              msg_li = [str(el) for el in msg_li]
              if run_st == pko.scode_d["Error"]:
                     err_ms += "\t" + "\t".join(msg_li) + "\n"
              else:
                     other_ms += "\t" + "\t".join(msg_li) + "\n"
       fin_msg = "Ready Tasks:\n"
       fin_msg += ready_ms + "\n"
       fin_msg += "Erroneous Tasks:\n"
       fin_msg += err_ms + "\n"
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
              run_st = pko.task_dict["status"]["run"]
              down_st = pko.task_dict["status"]["download"]
              to_write = [in_dir, for_fi, spike, down_st, run_st, pk_msg]
              to_write = [str(e) for e in to_write]
              msg = '\t'.join(to_write)
              print(msg)

if act["allp"]:
       fq_file = "/srv/cfm/outputs/ENA_out/PRJNA449565/410658/SAMN08911735/SRR7004509"
       fq_file += "/SRR7004509.pk"
       pko = TaskPickle(fq_file)
       # pko.task_dict["status"]["run"] = pko.scode_d["Queue"]
       # pko.write()
       # print(bool(pko))
       print(pko.__dict__)

