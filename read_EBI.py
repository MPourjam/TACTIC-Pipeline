#! /usr/bin/env python
import random
from getENA_modified import download_main
from Prep_ENA_out import preparation_main
from task_caller import main as process_multi

f = "/srv/cfm/outputs/EBI_1400t_extract.tab"
dest_dir = "/srv/cfm/outputs/ENA_out"

random.seed(140)
to_down = random.sample(range(1400000), 10)

download_main(f, outdir="/srv/cfm/outputs/", chunk_size=10, ndf=1)
process_multi(preparation_main(dest_dir))

