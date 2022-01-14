#! /usr/bin/env python
import argparse
from sys import version_info
import pandas as pd
from multiprocessing.pool import ThreadPool
from tqdm import tqdm  # pip3 install tqdm
from hashlib import md5
from itertools import repeat
from os import getcwd, cpu_count
from task_classes import task_pickle
if version_info[0] < 3:
    from pathlib2 import Path  # pip2 install pathlib2
else:
    from pathlib import Path
try:
    from urllib.request import urlretrieve
    from urllib.error import (ContentTooShortError, URLError)
except ImportError:
    from urllib import urlretrieve
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
            outfile = outpath.joinpath(row['study_accession'])
            outfile = outfile.joinpath(str(row['tax_id']))
            outfile = outfile.joinpath(row['sample_accession'])
            outfile = outfile.joinpath(row['run_accession'])
            outfile.mkdir(exist_ok=True, parents=True)
            # Creating the Args dictionary an storing it as pickle
            paired = "Yes" if row["library_layout"] == "PAIRED" else "No"
            pk_name = "{}.pk".format(row['run_accession'])
            pickle_file_path = outfile.joinpath(pk_name)
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
                    progress_code = pk_obj.scode_d["Progress"]
                    pk_obj.task_dict["status"]["download"] = progress_code
                    pk_obj.write()
                    urlretrieve('ftp://' + pair, outfile)
                except Exception:
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
    except ContentTooShortError:
        print('Retry {}/15'.format(attmp + 1), *url_path)
        urlretrieve_converter(url_path, attmp + 1)
    except URLError:
        print('Retry {}/15'.format(attmp + 1), *url_path)
        urlretrieve_converter(url_path, attmp + 1)


def download_main(metadata_file,
                  outdir=None,
                  threads=None,
                  chunk_size=0,
                  ndf=0):
    # if not isinstance(acc, list):
        # exit("acc argument must be a list of accession numbers")
    if not threads or not isinstance(threads, int):
        threads = int(cpu_count() * 0.3)
    threads = threads if threads > 0 else 1
    if chunk_size and not isinstance(chunk_size, int):
        try:
            chunk_size = int(chunk_size)
        except Exception:
            exit("Chunk size must be an integer value!")
    if ndf and not isinstance(chunk_size, int):
        try:
            ndf = int(ndf)
        except Exception:
            exit("Number of dataframes to parse must be an integer value!")
    metadata_file = Path(metadata_file)
    outputpath = Path(outdir) if outdir else Path(getcwd())
    outputpath = outputpath/'ENA_out'
    #
    if not metadata_file.exists() or not metadata_file.is_file():
        exit("Metadata file does not exit.")
    with open(metadata_file, 'r+') as metfile:
        line = metfile.readline()[:-1]
        column_names = line.split("\t")
        lc = 0
        while line:
            lc += 1
            line = metfile.readline()[:-1]
    chunk_size = chunk_size if chunk_size else lc
    # Spliting the tsv file to dfs
    m = lc // chunk_size
    rem = lc - (m * chunk_size)
    if ndf < m and ndf != 0:
        m = ndf
        rem = 0
    starts = []
    for r in range(m):
        start = int(r * chunk_size) + 1
        end = chunk_size
        starts.append((start, end))
    if rem > 0:
        start = int(m * chunk_size) + 1
        end = chunk_size
        starts.append((start, end))
    for start, read_n in starts:
        arg_d = {"filepath_or_buffer": metadata_file,
                 "sep": "\t",
                 "skiprows": start,
                 "nrows": read_n,
                 "header": None,
                 "names": column_names
                 }
        concat_frames = pd.read_csv(**arg_d)
        genome = []
        with ThreadPool(threads) as p:
            multipleargs = list(zip(concat_frames.iterrows(),
                                    repeat(outputpath)))
            tqdm_desc = 'Downloading Genomes using {} threads'.format(threads)
            for result in tqdm(p.imap_unordered(download_fastq, multipleargs),
                               total=len(multipleargs),
                               desc=tqdm_desc,
                               unit='Genomes'):
                genome.append(result)


def file_path(string):
    try:
        path = Path(string)
        if not path.is_file():
            raise
        else:
            return path.resolve()
    except Exception:
        raise FileNotFoundError(string)


if __name__ == '__main__':
    V = '%(prog)s v1.2.5'
    args_desc = 'Download FASTQ files from ENA ({})'.format(V)
    parser = argparse.ArgumentParser(description=args_desc)
    tab_help = "Path to the tab-delimited file from ENA metadata"
    parser.add_argument('-tab', '--tab',
                        type=file_path,
                        help=tab_help)
    parser.add_argument('-o', '--outfiles', type=str,
                        default=getcwd() + '/' + 'ENA_out')
    parser.add_argument('-t', '--threads', type=int,
                        default=int(cpu_count() * 0.5),
                        help='Number of threads to use (default 12)')
    parser.add_argument('-csize', '--chunk-size', type=int,
                        default=10,
                        help='Number of rows of each dataframe.'
                             '0 implies the whole rows of tab file')
    parser.add_argument('-ndf', '--number-dfs', type=int,
                        default=0,
                        help='Number of dataframes to parse.'
                             '0 implies all dataframes')
    parser.add_argument('-V', '--version', action='version',
                        version=V)
    args = parser.parse_args()
    args.outfiles = Path(args.outfiles)/"ENA_out"
    download_main(args.tab, args.outfiles, args.threads,
                  args.chunk_size, args.number_dfs)
