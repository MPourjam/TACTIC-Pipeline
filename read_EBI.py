#! /usr/bin/env python
from random import randint
from getENA_modified import main
from Prep_ENA_out import preparation_main
from task_caller import main as process_multi

f = "/cfm/EBI_1400t_extract.tab"
dest_dir = "/srv/cfm/outputs/ENA_out"

to_down = [randint(1, 100) for i in range(1, 20)]

accs = []
with open(f, "r") as ebi:
    i = 0
    line = ebi.readline()[:-1]
    while i < max(to_down) + 1:
        if i in to_down:
            acc = line.split("\t")[3]
            # print(f"{line} -- {acc} -- {i}\n")
            accs.append(acc)
        line = ebi.readline()[:-1]
        i += 1

main(accs, dest_dir)
process_multi(preparation_main(dest_dir))

