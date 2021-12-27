#! /usr/bin/python3.7
from processing_job import main_processing
import pickle as pk

file_path = "/srv/cfm/inputs/proc_dicts_20211227_1442.pk"

with open(file_path, 'rb') as f:
    dicts_list = pk.load(f)

first = dicts_list[0]

main_processing(**first)