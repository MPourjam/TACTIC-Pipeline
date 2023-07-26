#!/usr/bin/env python3

import logging
import argparse as ap

logging.basicConfig(
    level=logging.DEBUG
)

parser = ap.ArgumentParser()

parser.add_argument('cutoff', type=int)

parser.add_argument('readcountstatsfile', type=ap.FileType('r'))
parser.add_argument('pairedfile', type=ap.FileType('r'))
parser.add_argument('mappingfile', type=str)

parser.add_argument('outstatsfiltered', type=ap.FileType('w'))
parser.add_argument('outpairedfiltered', type=ap.FileType('w'))
parser.add_argument('outmappingfilefiltered', type=str)

parser.add_argument('filteredfiles', type=ap.FileType('w'))

args = parser.parse_args()


# filter stats file and identify matching IDs
good_samples = []
number_of_samples = 0
for line in args.readcountstatsfile:
    if line.startswith('#'):
        args.outstatsfiltered.write(line)
        continue

    number_of_samples += 1

    fields = line.split('\t')
    count = int(fields[1])

    if count >= args.cutoff:
        good_samples.append(fields[0])
        args.outstatsfiltered.write(line)
    else:
        args.filteredfiles.write(fields[0] + "\n")

# filter pairing file for good IDs
for line in args.pairedfile:
    if line.startswith('#'):
        args.outpairedfiltered.write(line)
        continue

    if line == "Forward\tReverse\tID\tfasta\n":
        args.outpairedfiltered.write(line)
        continue

    fields = line.split('\t')
    identifier = fields[2]

    if identifier in good_samples:
        args.outpairedfiltered.write(line)

# filter mapping file for good IDs
if args.mappingfile != '0':
    with open(args.mappingfile, 'r') as mapfile, open(args.outmappingfilefiltered, 'w') as mapfilefiltered:
        for line in mapfile:
            if line.startswith('#'):
                mapfilefiltered.write(line)
                continue

            fields = line.split('\t')
            identifier = fields[0]

            if identifier in good_samples:
                mapfilefiltered.write(line)


logging.info(f'Readcount filter (cutoff={args.cutoff}): Kept {100*len(good_samples)/number_of_samples}% ({len(good_samples)}/{number_of_samples}) samples')

args.readcountstatsfile.close()
args.outstatsfiltered.close()
args.pairedfile.close()
args.outpairedfiltered.close()
