#! /usr/bin/env python
import argparse
from sys import version_info
import pandas as pd
import pickle as pk
from multiprocessing.pool import ThreadPool
from tqdm import tqdm  # pip3 install tqdm
from hashlib import md5
from itertools import repeat
from os import getcwd, cpu_count
from datetime import datetime
# from mimetypes import guess_type
# from datetime import datetime as dt
from task_classes import task_pickle
if version_info[0] < 3:
    from StringIO import StringIO
    from pathlib2 import Path  # pip2 install pathlib2
else:
    from io import StringIO
    from pathlib import Path
try:
    from urllib.request import urlretrieve
    from urllib.request import urlopen
    from urllib.error import (ContentTooShortError, URLError)
except ImportError:
    from urllib import urlretrieve
    from urllib import urlopen
    from urllib.error import (ContentTooShortError, URLError)


global tsvurl
tsvurl = 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession='
tsvurl += '{}&fields=all&result=read_run'
status_dict = {"Done": 0, "Progress": 1, "Error": -1}


def md5sum(filename, blocksize=65536):
    hash = md5()
    with open(filename, 'rb') as f:
        for block in iter(lambda: f.read(blocksize), b''):
            hash.update(block)
    return hash.hexdigest()


def download_fastq(inputdata):
    df, outpath = inputdata
    index, row = df
    codename = row.name
    md5cheked = False
    while not md5cheked:
        listmd5 = []
        for pair, md5pair in zip(row['fastq_ftp'].split(';'),
                                 row['fastq_md5'].split(';')):
            # TODO replace the experiment_accession with tax_id if applicable
            outfile = outpath.joinpath(row['study_accession'])
            outfile = outfile.joinpath(row['experiment_accession'])
            outfile = outfile.joinpath(row['sample_accession'])
            outfile = outfile.joinpath(row['run_accession'])
            outfile.mkdir(exist_ok=True, parents=True)
            # Creating the Args dictionary an storing it as pickle
            paired = "Yes" if row["library_layout"] == "PAIRED" else "No"
            pickle_file_path = outfile.joinpath("{}.pk".format(row['run_accession']))
            pk_obj = task_pickle(pickle_file_path)
            pk_obj.task_dict["args"]["input_dir"] = str(outfile.resolve())
            pk_obj.task_dict["args"]["paired"] = paired
            pk_obj.task_dict["args"]["input_id"] = str(row['run_accession'])
            pk_obj.task_dict["status"]["download"] = pk_obj.scode_d["Started"]
            pk_obj.write()
            outfile = outfile.joinpath(pair.split('/')[-1])
            if outfile.is_file() and md5sum(outfile) == md5pair:
                listmd5.append(1)
            else:
                try:
                    pk_obj.task_dict["status"]["download"] = pk_obj.scode_d["Progress"]
                    pk_obj.write()
                    urlretrieve('ftp://' + pair, outfile)
                except Exception as e:
                    listmd5.append(0)
                    break
                if md5sum(outfile) == md5pair:
                    listmd5.append(1)
                else:
                    listmd5.append(0)
        if all(listmd5):
            pk_obj.task_dict["status"]["download"] = pk_obj.scode_d["Done"]
            pk_obj.write()
            md5cheked = True
    return '[OK] {}'.format(codename)


def urlretrieve_converter(url_path, attmp=0):
    if attmp > 15:
        print('Error', *url_path)
        return False
    try:
        return urlretrieve(*url_path)
    except ContentTooShortError as e:
        print('Retry {}/15'.format(attmp + 1), *url_path)
        urlretrieve_converter(url_path, attmp + 1)
    except URLError as e:
        print('Retry {}/15'.format(attmp + 1), *url_path)
        urlretrieve_converter(url_path, attmp + 1)


def main(acc, outdir=None, threads=None):
    if not isinstance(acc, list):
        exit("acc argument must be a list of accession numbers")
    if not threads or not isinstance(threads, int):
        threads = int(cpu_count() * 0.3)
    threads = threads if threads > 0 else 1
    outputpath = Path(outdir) if outdir else Path(getcwd())
    tmpoutputpath = outputpath/'tmp'
    tmpoutputpath.mkdir(exist_ok=True, parents=True)
    accurl = [tsvurl.format(accid) for accid in acc]
    accout = [tmpoutputpath/(accid+'.tsv') for accid in acc]
    metadata = []
    multipleargs = [(u, a) for (u, a) in zip(accurl, accout) if not a.is_file()]
    with ThreadPool(threads) as p:
        for result in tqdm(p.imap_unordered(urlretrieve_converter, multipleargs),
                           total=len(multipleargs),
                           desc='Downloading metadatas using {} threads'.format(threads),
                           unit='metadatas'):
            metadata.append(result)
    frames = [pd.read_csv(tsv, sep='\t') for tsv in accout]  # , index_col=5
    concat_frames = pd.concat(frames, ignore_index=True)
    genome = []
    with ThreadPool(threads) as p:
        multipleargs = list(zip(concat_frames.iterrows(), repeat(outputpath)))
        for result in tqdm(p.imap_unordered(download_fastq, multipleargs),
                           total=len(multipleargs),
                           desc='Downloading Genomes using {} threads'.format(threads),
                           unit='Genomes'):
            genome.append(result)
    concat_frames.to_csv(tmpoutputpath/'metadata_{}.tsv'.format(datetime.today().strftime('%Y%m%d')),
                         index=False,
                         sep='\t')


if __name__ == '__main__':
    V = '%(prog)s v1.2.5'
    parser = argparse.ArgumentParser(description='Download FASTQ files from ENA ({})'.format(V))
    parser.add_argument('-acc', '--acc', nargs='+')
    parser.add_argument('-o', '--outfiles', type=str,
                        default=getcwd() + '/' + 'ENA_out')
    parser.add_argument('-t', '--threads', type=int,
                        default=int(cpu_count() * 0.5),
                        help='Number of threads to use (default 12)')
    parser.add_argument('-V', '--version', action='version',
                        version=V)
    args = parser.parse_args()
    args.outfiles = Path(args.outfiles)/"ENA_out"
    main([args.acc], args.outfiles, args.threads)
