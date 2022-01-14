#! /usr/bin/env python
from download_ena import download_main
from process_prep import preparation_main
from task_caller import main as process_multi

f = "/srv/cfm/outputs/EBI_1400t_extract.tab"
down_dir = "/srv/cfm/outputs/"
dest_dir = down_dir + "ENA_out"

download_main(f, outdir=down_dir, chunk_size=3, ndf=1)
tasks_pickle_file = preparation_main(dest_dir)
process_multi(tasks_pickle_file)

