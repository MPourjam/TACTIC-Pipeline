#!/usr/bin/env python3

import argparse as ap
import configparser
import subprocess
import sys
from typing import List


def log_to_status_file(msg):
    with open('/usr/local/bin/status.txt', 'w') as stat_f_h:
        stat_f_h.write(msg)


bin = "/usr/local/bin"
executable = "offline-analysis.pl"

parser = ap.ArgumentParser()
parser.add_argument('datadir', type=str, help='Folder with fastq folder')
parser.add_argument('config', type=str, help='Config file for IMNGS')

args = parser.parse_args()
print(args)

config = configparser.ConfigParser()
config.read(args.config)


def exit_on_error(msg):
    print(msg, file=sys.stderr)
    exit(1)


def run_command(command: List[str]):
    log_to_status_file('Starting IMNGS worker')
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    for line in iter(process.stdout.readline, b''):  # replace '' with b'' for Python 3
        print(line.decode('utf-8'))
    process.wait()
    return process.returncode

log_to_status_file('Reading in IMNGS Config')

if 'IMNGS_Settings' not in config.sections():
    exit_on_error("Faulty config file settings.ini")

imngs_settings = config['IMNGS_Settings']

mapfilepath = imngs_settings.get('mapping_file')
if mapfilepath == 0:
    mapfilepath = 0

outpath = imngs_settings.get('outpath')

pipeline_reference = imngs_settings.get('pipeline_reference')

mode = imngs_settings.get('mode')  # otu or zotu
if mode not in ['otu', 'zotu']:
    exit_on_error("unknown mode. otu or zotu")

isPaired = imngs_settings.getint('isPaired')
if isPaired not in [0, 1]:
    exit_on_error("isPaired must be 0 or 1")

twoIndexes = imngs_settings.getint('twoIndexes')
if twoIndexes not in [0, 1]:
    exit_on_error("twoIndexes must be 0 or 1")

runDemux = imngs_settings.getint('runDemux')
if runDemux not in [0, 1]:
    exit_on_error("runDemux must be 0 or 1")

allow_barcode_mismatch = imngs_settings.getint('allow_barcode_mismatch')
if twoIndexes not in [0, 1, 2]:
    exit_on_error("allow_barcode_mismatch must be <=2")

forward_trim = imngs_settings.getint('forward_trim')
reverse_trim = imngs_settings.getint('reverse_trim')
trim_score = imngs_settings.getint('trim_score')

expected_error_rate = imngs_settings.getfloat('expected_error_rate')
if expected_error_rate < 0 or expected_error_rate > 1:
    exit_on_error("expected_error_rate must be between 0 and 1")

minmergelen = imngs_settings.getint('minmergelen')
maxmergelen = imngs_settings.getint('maxmergelen')

abundance = imngs_settings.getfloat('abundance')
if abundance < 0 or abundance > 1:
    exit_on_error("abundance must be between 0 and 1")

maxdiffpct = imngs_settings.getint('maxdiffpct')
if maxdiffpct < 0 or maxdiffpct > 100:
    exit_on_error("maxdiffpct must be between 0 and 100")

lowreadsamplecutoff = imngs_settings.getint('lowreadsamplecutoff')

clean_output = imngs_settings.getint('cleanoutput')
if not (clean_output == 0 or clean_output == 1):
    exit_on_error(f"clean_output must be 0 or 1; got {clean_output!r}")

print("Run settings:")
print(f'  pipeline_reference: {pipeline_reference}')
print(f'  mapping file: {mapfilepath}')
print(f'  output path: {outpath}')
print(f'  isPaired: {isPaired}')
print(f'  twoIndexes: {twoIndexes}')
print(f'  runDemux: {runDemux}')
print(f'  allow_barcode_mismatch: {allow_barcode_mismatch}')
print(f'  forward_trim: {forward_trim}')
print(f'  reverse_trim: {reverse_trim}')
print(f'  trim_score: {trim_score}')
print(f'  expected_error_rate: {expected_error_rate}')
print(f'  minmergelen: {minmergelen}')
print(f'  maxmergelen: {maxmergelen}')
print(f'  abundance: {abundance}')
print(f'  maxdiffpct: {maxdiffpct}')
print(f'  lowreadsamplecutoff: {lowreadsamplecutoff}')
print(f'  clean output directory: {clean_output}')

args_list = [args.datadir, isPaired, twoIndexes, runDemux, allow_barcode_mismatch, forward_trim,
             reverse_trim, trim_score, expected_error_rate, minmergelen, maxmergelen, abundance,
             maxdiffpct, lowreadsamplecutoff, mapfilepath, outpath, clean_output, pipeline_reference]

if mode == 'zotu':
    cmd = [f'{bin}/{executable}', '1']
else:
    cmd = [f'{bin}/{executable}', '0']

for elem in args_list:
    cmd.append(str(elem))

print(f'Running {cmd}')
print("[offline-analysis.pl] >>>")
exit_code = run_command(cmd)

log_to_status_file('IMNGS worker finished')
print('exit_code', exit_code)
exit(exit_code)
