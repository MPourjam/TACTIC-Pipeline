#! /usr/bin/python3.8
from pathlib import Path
import argparse


def get_derep_n(filepath):
    line_n = 0
    with open(filepath) as derep:
        line = derep.readline()
        while line:
            if line.startswith(">"):
                line_n += 1
            line = derep.readline()
    return line_n


parser = argparse.ArgumentParser()
parser.add_argument("-d", "--directory",
                    help="The directory to find all derep.fasta inside",
                    default=".")
parser.add_argument("-r", "--recursive",
                    help="find derep.fasta files recursively or not.",
                    default=False,
                    action="store_true")
parser.add_argument("-o", "--output",
                    help="The path (file or directory) to write the report.")
args = parser.parse_args()

dirpath = Path(args.directory)
if args.recursive:
    derep_files = dirpath.glob("./**/derep.fasta")
else:
    derep_files = dirpath.glob("./derep.fasta")

if args.output:
    outputpath = Path(args.output).resolve()
else:
    outputpath = dirpath.resolve()
if not outputpath.is_file():
    outputpath = outputpath.joinpath("reads_n.csv")
print(outputpath)

with open(outputpath, "w+") as outputfile:
    outputfile.write("Sample,Reads\n")
    for df in derep_files:
        sample_name = df.parent.stem
        read_n = get_derep_n(df)
        outputfile.write("{},{}\n".format(sample_name, read_n))
