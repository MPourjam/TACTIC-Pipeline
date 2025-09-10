import os
import re
import shutil
import tarfile
from ete3 import Tree
from os import chdir
import os.path as ospath
from multiprocessing import cpu_count, Pool
from processing_helper import TACTICArgsParser, gimmelogger, MyCounter
from processing_helper import system_sub as sys_sub
from spike_normalizer import normalize_otu_table
from pathlib import Path, PurePath
from .analysis_helper import onelinefasta
import logging
from ticlust import TICAnalysis

import warnings
# In Python ≤3.11, invalid escape sequences (like "\d" outside raw strings) were allowed but discouraged.
# In Python 3.12, they raise SyntaxWarning by default to encourage cleaner, future-proof code.
warnings.filterwarnings("ignore", category=SyntaxWarning)


BIN_DIR = "/base/binaries/"
USEARCH_11_BIN = BIN_DIR + "usearch11.0.667_i86linux64"
USEARCH_8_bin = BIN_DIR + "usearch8.1"
DB_LOC = "/base/inputs/databases/"
SINA_ARB = DB_LOC + "SILVA_LATEST.arb"
SINA_BIN = BIN_DIR + "sina/sina"
GOLD_REFDB_USEARCH = DB_LOC + "SILVA-bac-16s-90.udb"
# USERS_DIR = '/srv/crc/users/'
ref16RNAdb_1 = DB_LOC + "silva-bac-16s-id90.fasta"
ref16RNAdb_2 = DB_LOC + "silva-arc-16s-id95.fasta"
SORT_ME_RNA_BIN = BIN_DIR + "sortmerna"
USEARCH_TAIL = '> /dev/null 2>&1'
krona_importtext = BIN_DIR + "Krona/KronaTools/scripts/ImportText.pl"
global SPIKE_STAT_HEADER, TAXED_ZOTU_FILE_NAME, DEREP_NEW_LABEL
SPIKE_STAT_FILE_COLS = ("#SampleID", "SpikeReads", "spikes_total_weight_in_g", "spike_amount", "parent_path")
SPIKE_STAT_HEADER = "\t".join(list(SPIKE_STAT_FILE_COLS))
TAXED_ZOTU_FILE_NAME = "taxed_ZOTUs.fasta"
global FILTERED_READS_FILE, TRIMMED_READS_FILE, DEREP_FILE_NAME
TRIMMED_READS_FILE = "filtered1.fastq"
FILTERED_READS_FILE = "filtered2.fasta"
DEREP_FILE_NAME = "sorted.fasta"
max_pool = int(cpu_count() * 0.7)
DEREP_NEW_LABEL = "Uniq"
global POOL_SIZE
POOL_SIZE = max_pool if max_pool > 0 else 1
global ANA_LOG
ANA_LOG = gimmelogger(
    only_file=False,
)


def system_sub(*args, **kwargs):
    return sys_sub(*args, logger_obj=ANA_LOG, **kwargs)


def read_file(filename):
    with open(filename) as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    return content


def relabel_fasta(fasta_file: str, new_label: str):
    work_dir = Path("./relabeling").resolve()
    # make work_dir even if it exists
    work_dir.mkdir(parents=True, exist_ok=True)
    seq_fasta_path = Path(fasta_file).resolve()
    # if tar.gz version of this_file is found, extract it
    tar_gz_this_file = seq_fasta_path.with_suffix(seq_fasta_path.suffix + ".tar.gz")
    # extract tar_gz_this_file to this_file
    if not seq_fasta_path.is_file() and tar_gz_this_file.is_file():
        with tarfile.open(tar_gz_this_file, "r:gz") as tar:
            tar.extractall(path=seq_fasta_path.parent)
        tar_gz_this_file.unlink()
    elif not seq_fasta_path.is_file():
        raise ValueError(f"{str(seq_fasta_path)} must contain path to each samples sequence file!")
    dataset_name = f"{new_label}"
    new_file_fasta = work_dir.joinpath(f"{dataset_name}_relabeled.fasta")
    new_file_fastq = work_dir.joinpath(f"{dataset_name}_relabeled.fastq")
    cmd_to_call_list = [
        USEARCH_11_BIN,
        '-fastx_relabel',
        str(seq_fasta_path),
        '-threads',
        # We ignore POOL_SIZE here as the operation is IO bound and not memory bound
        '1',
        '-prefix',
        f"sample={dataset_name};",
        '-keep_annots',
        '-fastaout',
        str(new_file_fasta),
        '-fastqout',
        str(new_file_fastq)
    ]
    if str(seq_fasta_path).endswith(".fasta"):
        cmd_to_call_list.pop(-1)
        cmd_to_call_list.pop(-1)
    system_sub(cmd_to_call_list, force_log=True)
    onelinefasta(str(new_file_fasta))
    return str(new_file_fasta), str(new_file_fastq)


def append_reads(seq_fasta, sample_id, append_to):
    seq_fasta_path = os.path.abspath(seq_fasta)
    dataset_name = f"{sample_id}"
    new_file = 'new_' + str(os.path.basename(seq_fasta_path).split('.')[0])
    new_file_fasta = new_file + '.fasta'
    new_file_fastq = new_file + '.fastq'
    cmd_to_call_list = [
        USEARCH_11_BIN,
        '-fastx_relabel',
        seq_fasta_path,
        '-threads',
        str(POOL_SIZE),
        '-prefix',
        f"sample={dataset_name};",
        '-keep_annots',
        '-fastaout',
        new_file_fasta,
        '-fastqout',
        new_file_fastq
    ]
    if seq_fasta_path.endswith(".fasta"):
        cmd_to_call_list.pop(-1)
        cmd_to_call_list.pop(-1)
    # cmd = USEARCH_11_BIN + ' -fastx_relabel ' + taxed_ZOTUs_file_path + ' -prefix ' + dataset_name
    # cmd += '. -fastaout new -keep_annots'
    system_sub(cmd_to_call_list, force_log=True)
    with open(new_file_fasta, 'r') as new_file:
        with open(append_to + ".fasta", 'a') as analysis_file:
            analysis_file.write(new_file.read())

    os.remove(new_file_fasta)
    onelinefasta(append_to + ".fasta")
    if ospath.exists(new_file_fastq):
        with open(new_file_fastq, 'r') as new_file:
            with open(append_to + ".fastq", 'a') as analysis_file:
                analysis_file.write(new_file.read())
        os.remove(new_file_fastq)


def gather_samples_files(map_lines: list, file_to_gather: str):
    sample_seq_files_path = []
    for entries in map_lines:
        sample_processed_path = Path(PurePath(entries[4])).joinpath(file_to_gather)
        sample_seq_files_path.append((sample_processed_path, entries[0]))
    with Pool(int(cpu_count() * 0.7)) as relabel_pool:
        result = relabel_pool.starmap_async(
            relabel_fasta,
            [
                (str(seq_file), sam_id)
                for seq_file, sam_id in sample_seq_files_path
            ]
        )
        relabel_pool.close()
        relabel_pool.join()
        relabeled_files = result.get()
    # gathering to a single file
    rlbld_fasta = Path(file_to_gather).resolve().with_suffix(".fasta")
    rlbld_fastq = Path(file_to_gather).resolve().with_suffix(".fastq")
    output_tuple = [rlbld_fasta, rlbld_fastq]
    with open(rlbld_fasta, 'wb') as rlbld_fasta_fio, open(rlbld_fastq, 'wb') as rlbld_fastq_fio:
        for relabeled_fasta, relabeled_fastq in relabeled_files:
            with open(relabeled_fasta, 'rb') as rel_fasta_fio:
                shutil.copyfileobj(rel_fasta_fio, rlbld_fasta_fio)
            if ospath.exists(relabeled_fastq):
                with open(relabeled_fastq, 'rb') as rel_fastq_fio:
                    shutil.copyfileobj(rel_fastq_fio, rlbld_fastq_fio)
    # delete parent directory of relabeled files
    first_parent = Path(relabeled_files[0][0]).parent
    shutil.rmtree(first_parent)
    output_tuple = [
        file_path if ospath.exists(file_path) and os.path.getsize(file_path) > 0 else ""
        for file_path in output_tuple
    ]
    return tuple(output_tuple)


def trim_sides(five_end_trim, three_end_trim):
    # cmd_0 = USEARCH_11_BIN + " -fastx_truncate analysis.fasta -stripright " + str(five_end_trim)
    # cmd_1 = " -stripleft " + str(three_end_trim) + " -fastaout filtered1.fasta "
    # system_sub(cmd_0 + cmd_1 + USEARCH_TAIL)
    cmd_to_call_list = [
        USEARCH_11_BIN,
        '-fastx_truncate',
        'analysis.fasta',
        '-stripright',
        str(five_end_trim),
        '-stripleft',
        str(three_end_trim),
        '-fastaout',
        'trimmed.fasta',
        '-threads',
        str(POOL_SIZE)
    ]
    system_sub(cmd_to_call_list, force_log=True)
    # cmd_2 = 'mv filtered1.fasta analysis.fasta'
    cmd_to_call_list = [
        "mv",
        "filtered1.fasta",
        "analysis.fasta",
    ]
    # system_sub(cmd_to_call_list)
    onelinefasta('trimmed.fasta')


def filter_trimmed_seqs(trimmed_seqs_file: str, maxee_rate: float) -> str:
    # cmd_0 = USEARCH_11_BIN + " -fastq_filter filtered1.fastq -fastq_maxee_rate " + str(ARGS_CLS.filter_merged.fastq_maxee_rate)
    # cmd_1 = " -fastaout filtered2.fasta >/dev/null 2>/dev/null"
    # system_sub(cmd_0 + cmd_1)
    quality_filtered_fasta = "filtered2.fasta"
    system_sub(
        [
            *list(USEARCH_11_BIN.split(" ")),
            "-fastq_filter",
            trimmed_seqs_file,
            "-fastq_maxee_rate",
            maxee_rate,
            "-fastaout",
            quality_filtered_fasta,
            '-threads',
            str(POOL_SIZE)
        ],
        force_log=True
    )
    return quality_filtered_fasta


def dereplication(reads_file, with_taxonomy=True) -> str:
    # cmd_part_0 = USEARCH_11_BIN + " -fastx_uniques analysis.fasta -relabel Zotu"
    # cmd_part_0 += " -fastaout derep.fasta -tabbedout relabel.tab --threads " + str(POOL_SIZE)
    # cmd_part_0 += " > /dev/null 2>/dev/null"
    # system_sub(cmd_part_0)
    dereped_reads_file = "derep.fasta"
    cmd_to_call_list = [
        USEARCH_11_BIN,
        "-fastx_uniques",
        reads_file,
        "-relabel",
        DEREP_NEW_LABEL,  # dereplication relaling
        "-fastaout",
        dereped_reads_file,
        "-tabbedout",
        "relabel.tab",
        "--threads",
        str(POOL_SIZE),
        '--sizein',
        '--sizeout',
        '-strand',
        'both'
    ]
    # dereplicate reads and anotate them by multiplicity
    system_sub(cmd_to_call_list, force_log=True)
    onelinefasta(dereped_reads_file)
    return dereped_reads_file


def sort_seqs(to_sort_fasta: str, by: str = 'size') -> str:

    to_sort_path = os.path.abspath(to_sort_fasta)
    dir_path, file_name = os.path.split(to_sort_path)
    sorted_fasta_file = os.path.join(dir_path, "sorted_" + str(by) + ".fasta")
    if by == 'size':
        cmd_to_call_list = [
            USEARCH_11_BIN,
            "-sortbysize",
            str(to_sort_path),
            "-fastaout",
            str(sorted_fasta_file),
            '-threads',
            str(POOL_SIZE)
        ]
        system_sub(cmd_to_call_list, force_log=True)
    elif by == 'length':
        cmd_to_call_list = [
            USEARCH_11_BIN,
            "-sortbylength",
            str(to_sort_path),
            "-fastaout",
            str(sorted_fasta_file),
            '-threads',
            str(POOL_SIZE)
        ]
        system_sub(cmd_to_call_list, force_log=True)
        # print(" ".join(cmd_to_call_list))
    else:
        raise ValueError("Sorting should be by either 'size' or 'length'")
    onelinefasta(sorted_fasta_file)
    return sorted_fasta_file


# cluster sequences to OTUs
def clusterZOTUs(sorted_uniques_file: str, minsize: int):

    cmd_to_call_list = [
        USEARCH_11_BIN,
        "-unoise3",
        sorted_uniques_file,
        "-minsize",
        str(minsize),
        "-zotus",
        "zotus.fasta",
        "-sizein",
        "-sizeout",
        "-strand",
        "both",
        "-threads",
        str(POOL_SIZE)
    ]
    system_sub(cmd_to_call_list, force_log=True)
    onelinefasta("zotus.fasta")


# Filter non 16S sequences
def filter16S(input_fasta_name: str, output_fasta_name: str, skip: bool = False):

    if skip:
        shutil.copyfile(input_fasta_name, output_fasta_name)
        onelinefasta(output_fasta_name)
    else:
        system_sub(
            [
                SORT_ME_RNA_BIN,
                "--threads",
                str(POOL_SIZE),
                "--ref",
                ref16RNAdb_1,
                "--ref",
                ref16RNAdb_2,
                "--reads",
                input_fasta_name,
                "--fastx",
                "good-" + input_fasta_name,
                "--other",
                str(ARGS_PREC.filter_16S.other),
                "--workdir",
                ".",
                "-e",
                str(ARGS_PREC.filter_16S.e),
                "--num_alignments",
                str(ARGS_PREC.filter_16S.num_alignments)
            ],
            force_log=True
        )
        system_sub(
            [
                "mv",
                "out/aligned." + input_fasta_name.split(".")[-1],
                output_fasta_name
            ]
        )
        system_sub(["rm", "-r", "idx", "kvdb", "out"])
    return output_fasta_name


# cluster sequences to OTUs
def clusterOTUs(
        sorted_uniques_file: str,
        minsize: int = 8):
    # cmd_part_0 = USEARCH_11_BIN + ' -cluster_otus good-ZOTUs-sized.fasta -fulldp -otus otus1.fa -minsize 1 '
    # cmd_part_1 = ' -uparseout z2o.tab'
    # cmd_part_2 = ' > /dev/null 2>/dev/null'
    # clusering of seq in OTUs (clustered)
    # system_sub(cmd_part_0 + cmd_part_1 + cmd_part_2)
    sorted_uniques_file = os.path.abspath(sorted_uniques_file)
    cmd_to_call_list = [
        USEARCH_11_BIN,
        "-cluster_otus",
        # Input is a FASTA file containing quality filtered and globally trimmed reads from
        # a marker gene amplicon sequencing experiment, e.g. 16S or ITS.
        sorted_uniques_file,
        "-fulldp",
        "-otus",
        "otus1.fa",
        "-minsize",
        str(minsize),
        "-uparseout",
        "z2o.tab",
        "-relabel",
        "OTU",
        "-strand",
        "both",
        "-threads",
        str(POOL_SIZE)
    ]
    system_sub(cmd_to_call_list, force_log=True)
    onelinefasta("otus1.fa")


def build_ZOTU_table(raw_seq_file: str, zotu_seq_file: str, match_id: float) -> None:
    # cmd_0 = USEARCH_11_BIN + ' -otutab analysis.fasta -top_hit_only -zotus nochi_ZOTUs.fasta '
    # cmd_1 = ' -otutabout zotu_table.txt -id 0.97'
    # system_sub(cmd_0 + cmd_1 + USEARCH_TAIL)
    raw_seq_file_path = Path(raw_seq_file).absolute()
    zotu_seq_file_path = Path(zotu_seq_file).absolute()
    otu_tab = Path("ZOTUs-Table.tab").absolute()
    cmd_to_call_list = [
        USEARCH_11_BIN,
        "-otutab",
        str(raw_seq_file_path),
        "-top_hit_only",
        "-zotus",  # if the database sequences are already denoised which is the case here
        str(zotu_seq_file_path),  # This file contains ZOTUs sequences (both centroids of OTUs and matched ZOTUs to that OTU and non-chimeric)
        "-otutabout",
        str(otu_tab),
        "-id",
        str(round(float(match_id), 4)),  # ZOTU identity threshold
        "-threads",
        f"{str(POOL_SIZE)}",
        "-strand",
        "both"
    ]

    system_sub(cmd_to_call_list, force_log=True)

    return str(otu_tab)


def build_OTU_table(raw_seq_file: str, otu_seq_file: str) -> None:
    # cmd_0 = USEARCH_11_BIN + ' -otutab analysis.fasta -top_hit_only -zotus nochi_ZOTUs.fasta '
    # cmd_1 = ' -otutabout zotu_table.txt -id 0.97'
    # system_sub(cmd_0 + cmd_1 + USEARCH_TAIL)
    raw_seq_file_path = Path(raw_seq_file).absolute()
    otu_seq_file_path = Path(otu_seq_file).absolute()
    otu_tab = Path("OTUs-Table.tab").absolute()
    cmd_to_call_list = [
        USEARCH_11_BIN,
        "-otutab",
        str(raw_seq_file_path),
        "-top_hit_only",
        "-otus",  # if the database sequences are already denoised which is the case here
        str(otu_seq_file_path),  # This file contains ZOTUs sequences (both centroids of OTUs and matched ZOTUs to that OTU and non-chimeric)
        "-otutabout",
        str(otu_tab),
        "-id",
        "0.97",
        "-threads",
        f"{str(POOL_SIZE)}",
        "-strand",
        "both"
    ]
    # NOTE zotu_table.txt could miss some ZOTUs as the ZOTU's unique sequence in analysis.fasta could already match to
    # another ZOTU before reaching its corresponding original ZOTU sequence at smimilarities higher than 0.97. It
    # happens when the ZOTU is created
    system_sub(cmd_to_call_list, force_log=True)
    return str(otu_tab)


def get_samples_sizes(input_file, with_taxonomy=True):
    contents = read_file(input_file)[1:]
    tot_sizes_list = list()
    num_tokens = len(contents[0].split('\t')) - int(2 if with_taxonomy else 1)
    # Adding 0.001 to each sum to avoid division by zero
    for i in range(0, num_tokens):
        tot_sizes_list.append(0.0000001)
    for line in contents:
        line_tokens = line.split('\t')
        for i in range(len(tot_sizes_list)):
            tot_sizes_list[i] += float(line_tokens[i + 1])
    tot_sizes_list = [round(sz, 6) for sz in tot_sizes_list]
    return tot_sizes_list


def select_zotu_seqs():
    # cmd_0 = USEARCH_11_BIN + ' -fastx_getseqs nochi_ZOTUs.fasta -labels '  # good-ZOTUs.fasta -> nochi_ZOTUs.fasta
    # cmd_1 = 'filtered_zotu_table_list.txt -fastaout ZOTUs-Seqs.fasta '
    # system_sub(cmd_0 + cmd_1 + USEARCH_TAIL)
    cmd_to_call_list = [
        USEARCH_11_BIN,
        "-fastx_getseqs",
        "nochi_ZOTUs.fasta",
        "-labels",
        "ZOTU_map.tab",
        "-fastaout",
        "ZOTUs-Seqs.fasta",
        "-strand",
        "both",
        '-threads',
        str(POOL_SIZE)
    ]
    system_sub(cmd_to_call_list, force_log=True)
    onelinefasta("ZOTUs-Seqs.fasta")
    system_sub(["rm", "filtered_zotu_table_list.txt"])


def addTax_new(zotu_fasta_path):
    '''
    Input file is a fasta file of seqs
    '''
    zotu_fasta_path = os.path.abspath(zotu_fasta_path)
    dirname, filename = os.path.split(zotu_fasta_path)
    os.chdir(dirname)
    # cmd_part_0 = 'sina --in ' + str(zotu_fasta_path) + ' --search --meta-fmt csv '
    # cmd_part_1 = f'--threads {POOL_SIZE} --lca-fields tax_slv '
    # cmd_part_2 = '--db ' + SINA_ARB + ' --out test_{}'.format(filename)
    # cmd_part_3 = ' > /dev/null 2>/dev/null'
    # system_sub(classifier_dir + cmd_part_0 + cmd_part_1 + cmd_part_2 + cmd_part_3)
    cmd_to_call_list = [
        SINA_BIN,
        "--in",
        str(zotu_fasta_path),
        "--search",
        "--meta-fmt",
        "csv",
        "--threads",
        str(POOL_SIZE),
        "--lca-fields",
        "tax_slv",
        "--turn",
        "all",
        "--db",
        SINA_ARB,
        "--out",
        "test_{}".format(filename),
    ]
    system_sub(cmd_to_call_list, force_log=True)

    new_taxed_path = os.path.join(dirname, "taxed_{}".format(filename))
    seqs = read_file(str(zotu_fasta_path))
    zotunames = []
    zotuseqs = []
    for seql in seqs:
        if seql[0] == ">":
            zotunames.append(seql[1:])
        else:
            zotuseqs.append(seql)
    zotu_seq_dict = dict(zip(zotunames, zotuseqs))
    sina_test_out = os.path.join(dirname, 'test_{}.csv'.format("".join(filename.split(".")[:-1])))
    silva_contents_header = read_file(sina_test_out)
    silva_contents = silva_contents_header[1:]
    with open(new_taxed_path, 'w+') as out_file:
        for line in silva_contents:
            line_tokens = line.split(',')
            silva_taxo = line_tokens[6]
            silva_taxo_list = silva_taxo.strip().split(';')
            padded_silva_taxo_list = silva_taxo_list + [""] * (7 - len(silva_taxo_list))
            padded_silva_taxo_list = [str(p).strip() for p in padded_silva_taxo_list]
            full_silva_taxo = ";".join(padded_silva_taxo_list[:7])
            seq_name = line_tokens[0]
            out_file.write('>' + seq_name + ' ' + "tax=" + full_silva_taxo + '\n' + zotu_seq_dict[seq_name] + '\n')
    shutil.move(new_taxed_path, filename)


def addTax_to_table(table_file: str, taxed_seq_file: str, replace_originals: bool = True):
    """
    It finds the taxonomy of each sequence in the taxed_seq_file and adds it to the table_file

    :param table_file: str
        The path to the table file
    :param taxed_seq_file: str
        The path to the taxed sequence file
    :param replace_originals: bool
    :return: none
    """
    # Path handling
    tax_all_empty = ';;;;;;'
    table_file = os.path.abspath(table_file)
    taxed_seq_file = os.path.abspath(taxed_seq_file)
    tax_regex = re.compile(r"tax=(?P<tax>([^;]*;)*[^;]*);*", re.IGNORECASE)
    seq_id_regex = re.compile(r"^>(?P<seq_header>[^;\s]+).*$", re.IGNORECASE)
    seq_headers = {}
    # seq_taxonomies = []
    with open(taxed_seq_file) as taxed_seq_fio:
        seq_line = taxed_seq_fio.readline().strip()
        while seq_line:
            if seq_line[0] == '>':
                _tax = ''
                seq_id = seq_id_regex.search(seq_line).group('seq_header')
                curr_tax = tax_regex.search(seq_line)
                _tax += curr_tax.group('tax') if curr_tax.group('tax') else tax_all_empty
                seq_headers[seq_id] = _tax
            seq_line = taxed_seq_fio.readline().strip()

    taxed_table_file = 'taxed_' + os.path.basename(table_file)
    with open(table_file) as table_fio, open(taxed_table_file, 'w+') as out_fio:
        table_header = table_fio.readline().strip()
        if table_header.endswith('Taxonomy'):
            out_fio.write(table_header + "\n")
        else:
            out_fio.write(table_header + "\t" + "Taxonomy\n")
        table_line = table_fio.readline().strip()
        while table_line:
            if table_line[0] == '#':
                out_fio.write(table_line + "\n")
            else:
                seq_name = str(table_line.split('\t')[0]).strip()
                _tax = seq_headers.get(seq_name, tax_all_empty)
                _tax_list = _tax.strip().split(';')
                # Pad _tax_list to ensure it has at least 7 elements
                # This operation is redundant here as we already ensure that
                # _tax_list has 7 elements in the function addTax_new() but
                # we keep it anyway
                padded_tax_list = _tax_list + [""] * (7 - len(_tax_list))
                padded_tax_list = [str(p).strip() for p in padded_tax_list]
                new_line = table_line + '\t' + ";".join(padded_tax_list[:7]) + '\n'
                out_fio.write(new_line)
            table_line = table_fio.readline().strip()

    if replace_originals:
        shutil.move(taxed_table_file, table_file)


def shorten_sina_algn(sina_algn_file):
    sina_algn_path = os.path.abspath(sina_algn_file)
    sina_algn_dirname, sina_algn_filename = os.path.split(sina_algn_path)
    shortened_file_path = os.path.join(sina_algn_dirname, f'shortened_algn_{sina_algn_filename}')
    os.chdir(sina_algn_dirname)
    # Getting start and end positions of all zotus in sina alignment
    start_pos_set = set()
    end_pos_set = set()
    with open(sina_algn_path) as zotus_sina_algn:
        zotus_sina_algn_line = zotus_sina_algn.readline().strip()
        while zotus_sina_algn_line:
            if zotus_sina_algn_line.startswith('>'):
                # zotu_name = zotu_regex.search(zotus_sina_algn_line).group(1)
                # getting sequence start and end positions
                zotus_sina_algn_line = zotus_sina_algn.readline().strip()
                start_pos = re.search(r"-*[ACTGUN]", zotus_sina_algn_line).end()
                end_pos = re.search(r"[ACTGUN]-*$", zotus_sina_algn_line).start()
                start_pos_set.add(start_pos)
                end_pos_set.add(end_pos)
            zotus_sina_algn_line = zotus_sina_algn.readline().strip()
    # Writing common part of sina alignment for all zotus in a file
    lowest_start = min(start_pos_set)
    highest_end = max(end_pos_set)
    with open(sina_algn_file) as zotus_sina_algn, open(shortened_file_path, 'w+') as sina_algn_common_part:
        zotus_sina_algn_line = zotus_sina_algn.readline().strip()
        while zotus_sina_algn_line:
            if zotus_sina_algn_line.startswith('>'):
                sina_algn_common_part.write(zotus_sina_algn_line + '\n')
                # getting sequence start and end positions
            else:
                common_aln_seq = zotus_sina_algn_line[lowest_start:highest_end + 1]
                sina_algn_common_part.write(common_aln_seq + '\n')
            zotus_sina_algn_line = zotus_sina_algn.readline().strip()

    return shortened_file_path


def plant_tree(algn_file, method: list = ["rapidnj", "FastTree"], output_file_name_root: str = None):
    algn_file_path = os.path.abspath(algn_file)
    algn_dirname, algn_filename = os.path.split(algn_file_path)
    os.chdir(algn_dirname)
    out_file_name_root = str(algn_filename.split(".")[0])
    output_file_name = output_file_name_root if output_file_name_root else out_file_name_root
    output_file_name = os.path.join(algn_dirname, output_file_name)

    # create both trees
    if "FastTree" in method:
        # cmd_part_0 = 'FastTree -quiet -nosupport -gtr -nt filtered_aligned.fasta'
        # cmd_part_1 = ' > ' + sample + '_FastTree.tre'
        # system_sub(cmd_part_0 + cmd_part_1)
        cmd_to_call_list = [
            BIN_DIR + "FastTree",
            "-quiet",
            "-nosupport",
            "-gtr",
            "-nt",
            algn_file_path,
            ">",
            f"{output_file_name}-Tree-FastTree.tre",
        ]

    if "rapidnj" in method:
        # cmd_2_part_0 = BIN_DIR + 'rapidnj filtered_aligned.fasta --input-format fa --cores 6 --alignment-type d '
        # cmd_2_part_1 = '--output-format t --no-negative-length -x ' + sample + '_rapidNJ_tree.tre'
        # system_sub(cmd_2_part_0 + cmd_2_part_1)
        cmd_to_call_list = [
            BIN_DIR + "rapidnj",
            algn_file_path,
            "--input-format",
            "fa",
            "--cores",
            str(POOL_SIZE),
            "--alignment-type",
            "d",
            "--output-format",
            "t",
            "--no-negative-length",
            "-x",
            f"{output_file_name}-Tree-nj.tre",
        ]
    system_sub(cmd_to_call_list)


def create_krona(krona_tool: str, taxed_otu_table: str, output_html_name: str):
    taxed_otu_table_path = os.path.abspath(taxed_otu_table)
    taxed_otu_dirname, taxed_otu_filename = os.path.split(taxed_otu_table_path)
    OTUfinal_line = read_file(taxed_otu_table_path)[1:]
    krona_txt_file = taxed_otu_dirname + "/krona_" + os.path.basename(taxed_otu_filename)
    if output_html_name:
        output_html_name = os.path.join(taxed_otu_dirname, output_html_name)
    else:
        output_html_name = os.path.join(taxed_otu_dirname, "krona.html")
    krona_text_fio = open(krona_txt_file, 'w+')
    krona_dct = {}
    for OTU_taxonomy in OTUfinal_line:
        OTUFinal_list = OTU_taxonomy.split('\t')[1:]
        taxo = '\t'.join(str(OTUFinal_list[-1]).strip().split(";"))
        size = OTUFinal_list[:-1]
        curr_sum = krona_dct.get(taxo, 0) + sum([float(s) for s in size])
        krona_dct[taxo] = curr_sum
    for tax, size in krona_dct.items():
        krona_text_fio.write(str(size) + '\t' + tax + '\n')
    krona_text_fio.close()
    system_sub(
        [
            krona_tool,
            krona_txt_file,
            "-o",
            output_html_name
        ],
        capture_output=True,
        quiet=True
    )
    system_sub(
        [
            "rm",
            krona_txt_file
        ],
        capture_output=True,
        quiet=True
    )

    return output_html_name


def filter_otus_by_abundance(
        abundance_limit: float,
        OTU_table_file_name: str,
        OTU_fasta_file_name: str = "",
        sep: str = '\t',
        with_taxonomy: bool = True,
        sample_wise_corr_flag: bool = True,
        replace_originals: bool = True) -> None:
    # Abundance float fixation
    abundance_limit = round(float(abundance_limit), 5)
    # Managing paths of input files
    # OTU TABLE
    OTU_table_file_path = os.path.abspath(OTU_table_file_name)
    OTU_table_file_parent, OTU_table_file_name = os.path.split(OTU_table_file_path)
    filtered_OTU_table_path = os.path.join(OTU_table_file_parent, 'filtered_by_abundance_' + OTU_table_file_name)

    unf_tot_sizes = get_samples_sizes(OTU_table_file_path, with_taxonomy=with_taxonomy)
    # Filtering OTU Table
    removed_sotus = list()
    remained_sotus = list()
    with open(OTU_table_file_path, "r") as otu_table_fio, open(filtered_OTU_table_path, "w+") as filtered_otu_fio:
        otu_header = otu_table_fio.readline().strip()
        datasets_samples = len(otu_header.split(sep)) - int(2 if with_taxonomy else 1)
        sample_sizes = list()
        for i in range(0, datasets_samples):
            sample_sizes.append(0)
        filtered_otu_fio.write(otu_header + '\n')
        line = otu_table_fio.readline()
        while line:
            line = line.strip()
            line_list = line.split(sep)
            curr_sizes = list()
            sotu_name = str(line_list[0]).strip()
            rem_sotu = True
            curr_abundances = list()
            for j in range(len(sample_sizes)):
                abundance = round(float(line_list[j + 1]) / unf_tot_sizes[j], 5)
                curr_sizes.append(line_list[j + 1])
                curr_abundances.append(abundance)
            # if sample_wise_corr_flag:
            curr_sizes = [round(float(csz), 5) for csz in curr_sizes]
            corrected_curr_sizes = [round(float(csz), 5) for csz in curr_sizes]
            for ii in range(len(curr_abundances)):
                if curr_abundances[ii] <= abundance_limit:
                    corrected_curr_sizes[ii] = str(round(0.0, 1))
                else:
                    rem_sotu = False
            curr_sizes = corrected_curr_sizes if sample_wise_corr_flag else curr_sizes
            if not rem_sotu:
                # sometimes the Table creation causes some Z/OTUs to get lost
                # from the table file while they are present in the fasta file.
                # We create a white list of Z/OTUs that are present in the table
                # file and we filter the fasta file based on this white list
                remained_sotus.append(sotu_name)
                new_corrected_line = [sotu_name] + curr_sizes
                new_corrected_line += line_list[-1:] if with_taxonomy else []
                new_corrected_line = [str(el).strip() for el in new_corrected_line]
                new_line = sep.join(new_corrected_line)
                filtered_otu_fio.write(new_line + '\n')
            else:
                removed_sotus.append(sotu_name)
            # read new line
            line = otu_table_fio.readline()
    if replace_originals:
        shutil.move(filtered_OTU_table_path, OTU_table_file_path)

    # OTU FASTA
    OTU_fasta_file_path = os.path.abspath(OTU_fasta_file_name)
    if OTU_fasta_file_path and os.path.isfile(OTU_fasta_file_path):
        OTU_fasta_file_parent, OTU_fasta_file_name = os.path.split(OTU_fasta_file_path)
        filtered_otu_fasta_path = os.path.join(OTU_fasta_file_parent, 'filtered_by_abundance_' + OTU_fasta_file_name)
        seq_header_reg = re.compile(r"^>([^;\s]+)[;\s]+", re.IGNORECASE)
        with open(OTU_fasta_file_path, 'r') as otu_fasta_fio, open(filtered_otu_fasta_path, "w+") as filtered_otu_fasta_fio:
            derep_line = otu_fasta_fio.readline()
            write_seq = False
            while derep_line:
                if derep_line.startswith(">"):
                    curr_sotu = seq_header_reg.search(derep_line).group(1)
                    if curr_sotu and curr_sotu in remained_sotus:
                        new_line = derep_line
                        filtered_otu_fasta_fio.write(new_line)
                        write_seq = True
                    else:
                        write_seq = False
                else:
                    if write_seq:
                        filtered_otu_fasta_fio.write(derep_line)
                derep_line = otu_fasta_fio.readline()

        if replace_originals:
            shutil.move(filtered_otu_fasta_path, OTU_fasta_file_path)


def cleanup(directory: str, to_keep: list = []) -> bool:
    directory = Path(PurePath(directory)).absolute()
    must_keep = [
        re.compile(r"TACTICPipeline_args\.yml"),
        re.compile(r"spike_mapping_file\.csv"),
        # Pipelne outputs
        re.compile(r"[ZS]?OTUs-Seqs(-TAC|-TIC)?\.fasta"),
        re.compile(r"[ZS]?OTUs-Table(-TAC|-TIC)?\.tab"),
        re.compile(r"[ZS]?OTUs-Tree-nj(-TAC|-TIC)?\.tre"),
        re.compile(r"[ZS]?OTUs-Tree-FastTree(-TAC|-TIC)?\.tre"),
        re.compile(r"[ZS]?OTUs-Krona(-TAC|-TIC)?\.html"),
        re.compile(r"Map-.*\.tab"),
        re.compile(r"SampleMinCount_Normalized-[ZS]?OTUs-Table(-TAC|-TIC)?\.tab"),
        re.compile(r"Spike_Normalized-[ZS]?OTUs-Table(-TAC|-TIC)?\.tab"),
        re.compile(r"Analysis_log\.txt"),
        re.compile(r"Seqs(-TIC|-TAC)?\.fasta"),
        re.compile(r".*FullTaxonomy\.tab"),  # For ZOTU tables with full taxonomy
    ]
    white_list = [re.compile(f"{tk}") for tk in to_keep] + must_keep
    files = list(directory.glob("*"))
    for file in files:
        if file.is_file() and not any([pat.match(file.name) for pat in white_list]):
            file.unlink()
    return True


def uniqify_map_lines(map_lines_dict: dict) -> dict:
    """
    Given a list of map lines, it will return a list of map line tuples with
    """
    keys_list = list([entries[0] for ky, entries in map_lines_dict.items()])
    items_list = list(map_lines_dict.items())
    sample_ids_count = MyCounter(keys_list)
    most_common = MyCounter(dict(sample_ids_count.most_common(1)))
    while most_common and most_common.total() > 1:
        most_common_sample_id, most_common_count = list(most_common.items())[0]
        for ind in range(len(keys_list)):
            sam_id = keys_list[ind]
            entries = items_list[ind][1]
            long_id = items_list[ind][0]
            if most_common[sam_id]:
                # updating the sample id
                sam_id = str(most_common_sample_id) + "__" + str(most_common_count)
                entries[0] = sam_id
            map_lines_dict[long_id] = entries
        # End Phase of loop
        keys_list = list([entries[0] for ky, entries in map_lines_dict.items()])
        sample_ids_count = MyCounter(keys_list)
        most_common = MyCounter(dict(sample_ids_count.most_common(1)))

    return map_lines_dict


def parse_spike_stat_file(spike_stat_file: str, fastq_dir: str, parsed_spike_stat_path: str) -> dict:
    """
    Given a spike_stat file, it will parse it and return a dictionary
    containing the spike count and the original weight of each spike
    """
    spike_stat_file = Path(PurePath(spike_stat_file)).absolute()
    samples = {}
    fastq_dir = Path(PurePath(fastq_dir)).absolute()
    with open(spike_stat_file, 'r') as stats_h:
        for line in stats_h:
            if line.startswith("#"):
                if SPIKE_STAT_HEADER not in line.strip():
                    raise ValueError(f"Header: {line.strip()} does not match {SPIKE_STAT_HEADER}")
                continue
            elif line == "\n":
                continue
            else:
                fields = line.strip().split("\t")
                sample_id = str(fields[0]).strip()
                spike_reads = int(fields[1])
                try:
                    original_total_weight_in_g = round(float(fields[2]), 5)
                except ValueError:
                    original_total_weight_in_g = float("nan")
                amount = round(float(fields[3]), 5)
                try:
                    parent_path = str(fields[4]).strip()
                    # strip initial / or \
                    parent_path = parent_path.lstrip("/").lstrip("\\")
                except IndexError:
                    parent_path = ""
                abs_parent = fastq_dir.joinpath(Path(PurePath(parent_path)))
                full_id = abs_parent.joinpath(sample_id)
                samples[full_id] = (sample_id, spike_reads, original_total_weight_in_g, amount, abs_parent)

    samples = uniqify_map_lines(samples)
    parsed_spike_stat_file = spike_stat_file.parent.joinpath(f"parsed_{spike_stat_file.name}") if not parsed_spike_stat_path else parsed_spike_stat_path
    # write the parsed spike stat file
    with open(parsed_spike_stat_file, 'w') as parsed_stats_h:
        parsed_stats_h.write(f"{SPIKE_STAT_HEADER}\n")
        for sample_id, values in samples.items():
            values = list(values)
            values[4] = values[4].relative_to(fastq_dir)
            parsed_stats_h.write("\t".join([str(el) for el in values]) + "\n")
    ANA_LOG.info(f'Parsed spike stat file written to {str(parsed_spike_stat_file)}')
    return samples, parsed_spike_stat_file


def extract_and_rename_subtree(
        tree_file,
        subset,
        rename_map,
        output_file):
    # Load the tree from a file
    tree = Tree(tree_file)

    # Extract the sub-tree
    subtree = tree.copy()
    # print the tree
    subset = [f"'{sub}'" for sub in subset]
    rename_map = {f"'{k}'": f"'{v}'" for k, v in rename_map.items()}
    subtree.prune(subset, preserve_branch_length=True)

    # Rename the leaves
    for leaf in subtree:
        if leaf.name in rename_map:
            leaf.name = rename_map[leaf.name]

    # Save the sub-tree to a new file
    subtree.write(outfile=output_file)


def parse_uc_file(
        uc_file: str,
        cut_tax: bool = True) -> dict:
    """
    It parsed uc files from -uclust_smallmem. It returns a dictionary with
    the centroid as the key and the members as the values. The centroid is
    always the first element in the list of members

    :param uc_file: str
        The path to the uc file
    :return: dict
    """
    uc_file = Path(uc_file).absolute()
    tax_reg = re.compile(r"\s(?P<tax_tag>tax=)?(?P<tax>([^;]+;)*([^;]+)?;?)", re.IGNORECASE)

    uc_dict = {}
    with open(uc_file, 'r', encoding='utf-8') as uc_h:
        line = uc_h.readline().strip()
        while line:
            if line.startswith("H"):
                line_tokens = line.split("\t")
                query = line_tokens[8]
                target_centroid = line_tokens[9]
                sim_ = round(float(line_tokens[3]), 5)
                # omitting the taxonomic information
                if cut_tax:
                    query = tax_reg.sub("", query).strip()
                    target_centroid = tax_reg.sub("", target_centroid).strip()
                matched_clusters = uc_dict.get(query, [])
                matched_clusters.append((target_centroid, sim_))
                uc_dict[query] = matched_clusters
            elif line.startswith("S"):
                # case for singletons
                line_tokens = line.split("\t")
                centroid = line_tokens[8]
                sim_ = round(float(100), 5)
                # omitting the taxonomic information
                if cut_tax:
                    centroid = tax_reg.sub("", centroid).strip()
                query = centroid
                matched_clusters = uc_dict.get(query, [])
                matched_clusters.append((centroid, sim_))
                uc_dict[query] = matched_clusters
            line = uc_h.readline().strip()
    # some sequences might match more than one more centers
    # we will keep the one with the highest similarity
    cent_members_dict = {}
    for member, cent_matches in uc_dict.items():
        cent_matches = sorted(cent_matches, key=lambda x: x[1], reverse=True)
        cent_id = cent_matches[0][0]
        mems_list = cent_members_dict.get(cent_id, [])
        mems_list.append(member)
        cent_members_dict[cent_id] = mems_list

    return cent_members_dict


def collapse_zotu_table(
        zotu_cluster_map: str,
        zotu_tab_file: str,
        otu_tab_file: str) -> dict:
    """
    It creates an OTU table from the uc file

    :param zotu_cluster_map: dict
        The uc dictionary output of parse_uc_file
    :param zotu_tab_file: str
        The path to the zotu table file
    :return: dict
    """
    zotu_re = re.compile(r"(?P<zotu_id>Zotu\d+)[\s;]?", re.IGNORECASE)
    zotu_tab_path = Path(PurePath(zotu_tab_file)).absolute()
    otu_tab_path = Path(PurePath(otu_tab_file)).absolute()
    otu_zotu_map_path = otu_tab_path.parent.joinpath("Map-ZOTUs-OTUs-TAC.tab")
    # renaming the keys
    zotu_cluster_map = {"OTU" + str(i + 1): v for i, v in enumerate(zotu_cluster_map.values())}

    otus_row_dict: dict = {}
    with open(zotu_tab_path, 'r') as zotus_tab_fio:
        zotus_table_header = zotus_tab_fio.readline().strip()
        taxonomy_mo = re.search(r"\ttaxonomy$", zotus_table_header, re.IGNORECASE)
        with_taxonomy = True if taxonomy_mo else False
        zotus_tab_line = zotus_tab_fio.readline().strip()
        while zotus_tab_line:
            zotu_mo = zotu_re.search(zotus_tab_line)
            curr_zotu = zotu_mo.group("zotu_id")
            line_list = zotus_tab_line.split("\t")
            till = -1 if with_taxonomy else len(line_list)
            curr_counts = [float(count) for count in zotus_tab_line.split("\t")[1:till]]

            for otu_id, cluster in zotu_cluster_map.items():
                this_cluster_counts = [0.0] * len(curr_counts)
                this_otu_row_els = otus_row_dict.get(
                    otu_id,
                    [
                        otu_id,
                        *this_cluster_counts,
                        "",
                    ]
                )

                c_centroid = cluster[0]  # Centroid of the cluster is the first element
                cent_zotu = zotu_re.search(c_centroid).group("zotu_id")
                c_mems = cluster[1:]
                if curr_zotu == cent_zotu:
                    this_otu_row_els[1:-1] = [sum(x) for x in zip(curr_counts, this_otu_row_els[1:-1])]
                    this_otu_row_els[-1] = line_list[-1] if with_taxonomy else ""
                else:
                    for mem in c_mems:
                        mem_zotu = zotu_re.search(mem).group("zotu_id")
                        if curr_zotu == mem_zotu:
                            this_otu_row_els[1:-1] = [sum(x) for x in zip(curr_counts, this_otu_row_els[1:-1])]
                            break
                otus_row_dict[otu_id] = this_otu_row_els
            zotus_tab_line = zotus_tab_fio.readline().strip()

    sorted_rows = dict(sorted(otus_row_dict.items(), key=lambda x: int(x[0][3:])))
    otus_row_dict = sorted_rows
    with open(otu_tab_file, 'w+') as otu_tab_fio:
        otu_tab_fio.write(zotus_table_header + "\n")
        for otu_id, otu_row in otus_row_dict.items():
            otu_tab_fio.write("\t".join([str(el) for el in otu_row]) + "\n")

    sorted_map = dict(sorted(zotu_cluster_map.items(), key=lambda x: int(x[0][3:])))
    zotu_cluster_map = sorted_map
    with open(otu_zotu_map_path, 'w+') as otu_zotu_map_fio:
        for otu_id, zotus in zotu_cluster_map.items():
            for zotu_id in zotus:
                otu_zotu_map_fio.write(zotu_id + "\t" + otu_id + "\n")

    return zotu_cluster_map


def collapse_fasta_file(
        otu_zotu_map: dict,
        input_fasta_file: str,
        output_fasta_file: str) -> None:
    """
    It creates an OTU fasta file from the uc file

    :param zotu_otu_map: dict
        The zotu to otu map
    :param input_fasta_file: str
        The path to the zotu fasta file
    :param output_fasta_file: str
        The path to the otu fasta file
    :return: None
    """
    input_fasta_file_path = Path(PurePath(input_fasta_file)).absolute()
    output_fasta_file_path = Path(PurePath(output_fasta_file)).absolute()
    zotu_re = re.compile(r"(?P<zotu_id>Zotu\d+)[\s;]?", re.IGNORECASE)
    with open(input_fasta_file_path, 'r') as cent_fio, open(output_fasta_file_path, 'w+') as temp_fio:
        cent_line = cent_fio.readline().strip()
        while cent_line:
            if cent_line.startswith(">"):
                zotu_id = zotu_re.search(cent_line).group("zotu_id")
                for otu, zotus in otu_zotu_map.items():
                    if zotu_id in zotus:
                        otu_id = otu
                        cent_line = cent_line.replace(zotu_id, otu_id)
                        temp_fio.write(cent_line + "\n")
                        break
            else:
                temp_fio.write(cent_line + "\n")
            cent_line = cent_fio.readline().strip()


def uclust(
        input_fasta: str,
        cluster_id: float = 0.987) -> (str, str):
    """
    It runs uclust_smallmem on the input fasta file and
     returns the uc file path and the centroids fasta file path

    :param input_fasta: str
        The path to the input fasta file
    :return: str
    """
    input_fasta_path = Path(PurePath(input_fasta)).absolute()
    sorted_seq_file = sort_seqs(str(input_fasta_path), by="length")
    centroids_file = input_fasta_path.parent.joinpath("cluster_centroids.fasta")
    uc_file = input_fasta_path.parent.joinpath("otu_clusters.uc")
    cmd_to_call = [
        USEARCH_11_BIN,
        "-cluster_smallmem",
        str(sorted_seq_file),
        "-centroids",
        str(centroids_file),
        "-uc",
        str(uc_file),
        "-id",
        str(round(float(cluster_id), 4)),
        "-strand",
        "both",
        '-threads',
        str(POOL_SIZE)
    ]
    system_sub(cmd_to_call, force_log=True)
    onelinefasta(centroids_file)
    return str(uc_file), str(centroids_file)


def move_content_to_parent_directory(directory: str):
    """
    Moves all contents of 'directory' to its parent directory and removes 'directory'.
    :param directory: Path to the directory whose contents should be moved.
    """
    dir_path = Path(directory).resolve()

    if not dir_path.is_dir():
        raise ValueError(f"Directory '{dir_path}' does not exist or is not a directory.")

    parent_dir = dir_path.parent

    for item in dir_path.iterdir():
        target_path = parent_dir / item.name

        if target_path.exists():
            raise FileExistsError(f"Target '{target_path}' already exists. Move aborted to prevent overwriting.")

        shutil.move(str(item), str(target_path))

    dir_path.rmdir()  # Remove the empty directory


def zotu_pipeline(
        sorted_uniques_file: str,
        zotu_minsize: int,
        raw_seq_file: str,
        abund_limit: float,
        sample_wise_correction: bool,
        zotu_p_logger: logging.Logger = ANA_LOG,
        match_id: float = 0.99,
        filtered_non_bacterial: bool = True) -> None:
    """
    This function will run the ZOTU pipeline.
    """
    ZOTU_P_LOGGER = zotu_p_logger
    ZOTU_P_LOGGER.info('# Clustering at ZOTU level')
    clusterZOTUs(sorted_uniques_file, zotu_minsize)
    # we have already checked the sequences to be non-human reads.
    # let's not do filter16S step for now
    ZOTU_P_LOGGER.info(f'# Filtering out non-bacterial 16S: {not filtered_non_bacterial}')
    # filtered_non_bacterial means that the sequences does not include non-bacterial reads
    non_human_zotus_file = filter16S(
        input_fasta_name="zotus.fasta",
        output_fasta_name="ZOTUs-Seqs.fasta",
        skip=filtered_non_bacterial
    )

    ZOTU_P_LOGGER.info('# Creating ZOTU table')
    zotu_tab_file = build_ZOTU_table(raw_seq_file, non_human_zotus_file, match_id)
    ZOTU_P_LOGGER.info('# Filtering ZOTUs by abundance cutt-off of %s', abund_limit)
    filter_otus_by_abundance(
        abundance_limit=abund_limit,
        OTU_table_file_name=zotu_tab_file,
        OTU_fasta_file_name=non_human_zotus_file,
        with_taxonomy=False,
        replace_originals=True,
        sample_wise_corr_flag=sample_wise_correction
    )
    ZOTU_P_LOGGER.info('# Adding Taxonomy to ZOTU sequences')
    # taxed_ZOTUs-Seqs.fasta gets created here and then replaces ZOTUs-Seqs.fasta
    addTax_new(non_human_zotus_file)
    ZOTU_P_LOGGER.info('# Adding taxonomy to ZOTU table')
    addTax_to_table(zotu_tab_file, non_human_zotus_file)
    ZOTU_P_LOGGER.info('# Creating ZOTUs trees')
    sina_algn_shortened_file = shorten_sina_algn("test_ZOTUs-Seqs.fasta")
    try:
        plant_tree(
            algn_file=sina_algn_shortened_file,
            method="rapidnj",
            output_file_name_root="ZOTUs"
        )
    except Exception as exc:
        ZOTU_P_LOGGER.warning(f"Tree creation skipped: {exc}")
    ZOTU_P_LOGGER.info('# Adding ZOTUs Krona HTML')
    create_krona(
        krona_importtext,
        zotu_tab_file,
        output_html_name="ZOTUs-Krona.html"
    )
    # ZOTU_P_LOGGER.info('# Starting OTU pipeline')
    return (
        str(Path(PurePath(zotu_tab_file)).absolute()),
        str(Path(PurePath(non_human_zotus_file)).absolute())
    )


def otu_pipeline(
        sorted_uniques_file: str,
        raw_seq_file: str,
        abund_limit: float,
        sample_wise_correction: bool,
        otu_p_logger: logging.Logger = ANA_LOG,
        filtered_non_bacterial: bool = True,
        otus_minsize: int = 8) -> None:
    """
    This function will run the OTU pipeline
    """
    OTU_LOGGER = otu_p_logger
    OTU_LOGGER.info('# Clustering at OTU level')
    clusterOTUs(sorted_uniques_file, minsize=otus_minsize)
    # This filtering step depends on the passed argument to the whole pipeline
    # let's not do filter16S step
    OTU_LOGGER.info(f'# Filtering out non-bacterial 16S: {not filtered_non_bacterial}')
    # filtered_non_bacterial means that the sequences does not include non-bacterial reads
    non_human_otus_file = filter16S(
        input_fasta_name="otus1.fa",
        output_fasta_name="OTUs-Seqs.fasta",
        skip=filtered_non_bacterial
    )
    OTU_LOGGER.info('# Creating OTU table')
    otu_tab_file = build_OTU_table(raw_seq_file, non_human_otus_file)
    OTU_LOGGER.info('# Filtering OTUs by abundance cut-off of %s', abund_limit)
    filter_otus_by_abundance(
        abundance_limit=abund_limit,
        OTU_table_file_name=otu_tab_file,
        OTU_fasta_file_name=non_human_otus_file,
        with_taxonomy=False,
        replace_originals=True,
        sample_wise_corr_flag=sample_wise_correction
    )
    OTU_LOGGER.info('# Adding Taxonomy to OTU sequences')
    addTax_new(non_human_otus_file)
    OTU_LOGGER.info('# Adding taxonomy to OTU table')
    addTax_to_table(otu_tab_file, non_human_otus_file)
    OTU_LOGGER.info('# Creating OTUs trees')
    sina_algn_shortened_file = shorten_sina_algn("test_OTUs-Seqs.fasta")
    create_krona(
        krona_importtext,
        otu_tab_file,
        output_html_name="OTUs-Krona.html"
    )
    try:
        plant_tree(
            sina_algn_shortened_file,
            method="rapidnj",
            output_file_name_root="OTUs"
        )
    except Exception as exc:
        OTU_LOGGER.warning(f"Tree creation skipped: {exc}")
    OTU_LOGGER.info('# Adding OTUs Krona HTML')
    create_krona(
        krona_importtext,
        otu_tab_file,
        output_html_name="OTUs-Krona.html"
    )
    return (
        str(Path(PurePath(otu_tab_file)).absolute()),
        str(Path(PurePath(non_human_otus_file)).absolute())
    )


def tac_pipeline(
        zotus_seq_fasta: str,
        zotus_table_file: str,
        zotus_tree_file: str = None,
        cluster_id: float = 0.987,
        tac_logger: logging.Logger = ANA_LOG) -> None:
    """
    TAC (Taxonomy Agnostic Clustering) pipeline, takes a set of ZOTUs and cluster
     them at 98.7% similarity to create OTUs.

    :param zotus_seq_fasta: str
        The path to the sequences of non-human ZOTUs
    :param zotus_table_file: str
        The path to the ZOTUs table file

    :return: None
    """
    # Clustering ZOTUs at 98.7% similarity
    TAC_LOGGER = tac_logger
    zotus_seq_path = Path(PurePath(zotus_seq_fasta)).absolute()
    zotus_table_path = Path(PurePath(zotus_table_file)).absolute()
    otu_table_path = zotus_table_path.parent.joinpath("OTUs-Table-TAC.tab")
    TAC_LOGGER.info(f"# Clustering ZOTUs at {str(cluster_id)} similarity")
    uc_file, centroids_file = uclust(
        str(zotus_seq_path),
        cluster_id
    )
    uc_dict = parse_uc_file(uc_file)
    # print([zotus for key, zotus in uc_dict.items() if len(zotus) > 1])
    # exit()
    # Creating OTU table
    TAC_LOGGER.info('# Creating OTU table')
    otu_zotu_map = collapse_zotu_table(
        uc_dict,
        str(zotus_table_path),
        str(otu_table_path)
    )
    # renaming zotus in centroids file to otus
    otus_seq_fasta = zotus_seq_path.parent.joinpath("OTUs-Seqs-TAC.fasta")
    TAC_LOGGER.info('# Collapsing ZOTUs sequences to OTUs')
    collapse_fasta_file(
        otu_zotu_map,
        str(centroids_file),
        str(otus_seq_fasta)
    )
    # Updating Tree
    if zotus_tree_file:
        otus_tree_file = Path(PurePath(zotus_tree_file)).absolute()
        otus_tree_file = zotus_seq_path.parent.joinpath("OTUs-Tree-nj-TAC.tre")
        TAC_LOGGER.info('# Updating ZOTUs tree to OTUs tree')
        extract_and_rename_subtree(
            tree_file=zotus_tree_file,
            # centroids are the first in zotus list
            subset=[zotus[0] for otu, zotus in otu_zotu_map.items()],
            rename_map={zotus[0]: otu for otu, zotus in otu_zotu_map.items()},
            output_file=str(otus_tree_file)
        )

    # adding krona
    TAC_LOGGER.info('# Adding OTUs Krona HTML')
    create_krona(
        krona_importtext,
        str(otu_table_path),
        output_html_name="OTUs-Krona.html"
    )

    return (
        str(otu_table_path),
        str(otus_seq_fasta)
    )


def tic_pipeline(
        zotus_seq_fasta: str,
        zotus_table_file: str,
        zotus_tree_file: str = None,
        species_id: float = 0.987,
        genus_id: float = 0.95,  # These thresholds are taken from TIC paper
        family_id: float = 0.90,  # These thresholds are taken from TIC paper
        tic_logger: logging.Logger = ANA_LOG) -> None:
    """
    TIC (Taxonomy Informed Clustering) pipeline, takes a set of ZOTUs and cluster
     them at 98.7% similarity to create SOTUs, at 97% to create GUTUs and at 95%
     to create FOTUs.
    :param zotus_seq_fasta: str
        The path to the sequences of non-human ZOTUs
    :param zotus_table_file: str
        The path to the ZOTUs table file
    :param zotus_tree_file: str
        The path to the ZOTUs tree file
    :param species_id: float
        The similarity threshold to cluster at species level
    :param genus_id: float
        The similarity threshold to cluster at genus level
    :param family_id: float
        The similarity threshold to cluster at family level
    :param tic_logger: logging.Logger
        The logger to use
    :return None
    """
    zotus_seq_path = Path(zotus_seq_fasta).absolute()
    tic_analysis = TICAnalysis(
        taxed_fasta_file_path=zotus_seq_path,
        zotu_table_file=zotus_table_file,
        logger_obj=tic_logger
    )
    tic_logger.info(
        '# Clustering ZOTUs at %s, %s, %s similarity thresholds for '
        'species, genus and family levels respectively',
        species_id,
        genus_id,
        family_id
    )

    tic_analysis.run(
        threads=POOL_SIZE,
        cluster_thresholds_d={
            "species": species_id,
            "genus": genus_id,
            "family": family_id
        }
    )
    sotu_zotu_map = {}
    with open(tic_analysis.sotu_zotu_file_path, 'r') as sotu_zotu_fio:
        sotu_zotu_fio.readline().strip()
        line = sotu_zotu_fio.readline().strip()
        while line:
            sotu, zotu = line.split("\t")
            # map file is sorted and the first zotu is the centroid
            curr_zotu_list = sotu_zotu_map.get(sotu, [])
            curr_zotu_list.append(zotu)
            sotu_zotu_map[sotu] = curr_zotu_list
            line = sotu_zotu_fio.readline().strip()

    if zotus_tree_file:
        sotus_tree_file = Path(PurePath(zotus_tree_file)).absolute()
        sotus_tree_file = zotus_seq_path.parent.joinpath("SOTUs-Tree-nj-TIC.tre")
        tic_logger.info('# Updating ZOTUs tree to SOTUs tree')
        extract_and_rename_subtree(
            tree_file=zotus_tree_file,
            # centroids are the first in zotus list
            subset=[zotus[0] for otu, zotus in sotu_zotu_map.items()],
            rename_map={zotus[0]: otu for otu, zotus in sotu_zotu_map.items()},
            output_file=str(sotus_tree_file)
        )
    tic_logger.info('# Adding Krona HTML')
    create_krona(
        krona_importtext,
        str(tic_analysis.sotu_table_path),
        output_html_name="SOTUs-Krona.html"
    )
    # moving files in tic_analysis.tic_wd to parent directory
    # TODO check if it works. It seems that the file is not moved.
    shutil.move(tic_analysis.non_bact_fasta_path, tic_analysis.tic_wd / "Invalid-Tax-Seqs-TIC.fasta")
    shutil.move(tic_analysis.sotu_fasta_path, tic_analysis.tic_wd / "SOTUs-Seqs-TIC.fasta")
    shutil.move(tic_analysis.sotu_table_path, tic_analysis.tic_wd / "SOTUs-Table-TIC.tab")
    move_content_to_parent_directory(tic_analysis.tic_wd)

    return (
        str(tic_analysis.tic_wd.parent / "SOTUs-Table-TIC.tab"),
        str(tic_analysis.tic_wd.parent / "SOTUs-Seqs-TIC.fasta")
    )


def main(
        analysis_dir: str,
        spike_stat_file: str,
        fastqs_dir: str,
        args_file_path: str = "",
        dbs_loc: str = DB_LOC,
        threads=POOL_SIZE,
        usearch_11_bin: str = USEARCH_11_BIN,
        run_otu_pipeline: bool = False,
        run_tac_pipeline: bool = False,
        run_tic_pipeline: bool = False,
        skip_non_bacterial_filter: bool = False
    ) -> list:
    """
    sample_seq_files_path: a list of file path to
     each samples' sequence file which is going to be combined with other samples passed to analysis.
    analysis_dir: The destination directory to save results
    """
    global USEARCH_11_BIN, POOL_SIZE, ANALYSIS_DIR, ARGS_CLS, ANA_LOG, DB_LOC, ARGS_PREC
    # updating USEARCH_11_BIN
    USEARCH_11_BIN = usearch_11_bin

    spike_stat_file = Path(PurePath(spike_stat_file))
    assert spike_stat_file.is_file(), "spike_stat_file must be a path to a file"

    DB_LOC = Path(PurePath(dbs_loc))
    assert DB_LOC.is_dir(), "DBS_LOC should be path to direcotry containing SILVA database files (arb)"
    ANALYSIS_DIR = Path(PurePath(analysis_dir)).absolute()
    gimmelogger(
        logger_name="run_tactic.analysis",
        log_file=ANALYSIS_DIR.joinpath("Analysis_log.txt"),
        only_file=True,
    )
    ANA_LOG = gimmelogger(
        logger_name="run_tactic.analysis.main",
    )
    # updating POOL_SIZE
    try:
        POOL_SIZE = int(threads)
        ANA_LOG.info(f"Using {POOL_SIZE} threads")
    except ValueError:
        ANA_LOG.warning(f"threads must be an integer. Using default value of {POOL_SIZE}")

    assert ANALYSIS_DIR.is_dir(), "analysis_dir must be a path to a directory"
    ANALYSIS_DIR = str(ANALYSIS_DIR) + "/"
    cleanup(str(ANALYSIS_DIR))
    # copy the args_file_path to the analysis directory
    this_args_file = Path(PurePath(args_file_path)).absolute()
    new_args_file = Path(PurePath(ANALYSIS_DIR + this_args_file.name))
    try:
        shutil.copy2(this_args_file, new_args_file)
    except Exception as exc:
        ANA_LOG.warning(f"Copying args file to {ANALYSIS_DIR} failed: {exc}")

    ARGS_CLS = TACTICArgsParser(config_yaml=args_file_path).analysis_args
    ARGS_PREC = TACTICArgsParser(config_yaml=args_file_path).preproc_args
    # Parsing spike_stat_file
    try:
        parsed_spike_stat_file_path = Path(PurePath(ANALYSIS_DIR + 'spike_mapping_file.csv'))
        map_lines_dict, parsed_spike_stat_file_path = parse_spike_stat_file(
            spike_stat_file,
            fastq_dir=fastqs_dir,
            parsed_spike_stat_path=parsed_spike_stat_file_path)
    except Exception as exc:
        raise ValueError(f"spike_stat_file: {str(exc)}")

    ANA_LOG.info("# Gathering Sequences")
    chdir(ANALYSIS_DIR)
    # gather_samples_files(list(map_lines_dict.values()), TRIMMED_READS_FILE)
    # gather_samples_files(list(map_lines_dict.values()), FILTERED_READS_FILE)
    gather_samples_files(list(map_lines_dict.values()), DEREP_FILE_NAME)
    # trimming the sequences
    # NOTE for now we discard trimming in the analysis
    # trim_sides(ARGS_CLS.trimsides.stripleft, ARGS_CLS.trimsides.stripright)
    ANA_LOG.info('# Dereplication')
    derep_file = dereplication(DEREP_FILE_NAME, with_taxonomy=False)
    # Sorting the sequences
    ANA_LOG.info('# Sorting')
    sorted_file = sort_seqs(derep_file)
    # ##########################################################################################
    # #################################### ZOTU Pipline ########################################
    # ##########################################################################################
    to_return = ["", "", "", ""]
    try:
        ANA_LOG.info('# Starting ZOTU pipeline')
        zotu_tab_file, non_human_zotus_seq_file = zotu_pipeline(
            sorted_file,
            ARGS_CLS.cluster_zotus.minsize,
            # using dereplicated reads with size tag saves time in creation of zotu table
            # NOTE tested against using filtered reads and the results are the same
            DEREP_FILE_NAME,
            ARGS_CLS.cluster_zotus.abund_limit,
            ARGS_CLS.cluster_zotus.sample_wise_correction,
            ANA_LOG,
            ARGS_CLS.cluster_zotus.match_id,
            skip_non_bacterial_filter
        )
        # ANA_LOG.info('# Starting OTU pipeline')
        ANA_LOG.info('# Normalizing ZOTU Table')
        zotu_norm_methods = "No"
        zotu_norm_methods = normalize_otu_table(zotu_tab_file, str(parsed_spike_stat_file_path))
        ANA_LOG.info(f"{zotu_norm_methods} normalization method(s) applied on {zotu_tab_file}")
    except Exception as exc:
        ANA_LOG.warning(f"ZOTU pipeline failed: {exc}")
        raise exc
    else:
        # preserving important files for next steps
        # and cleaning up the rest
        cleanup(
            str(ANALYSIS_DIR),
            [
                Path(PurePath(sorted_file)).name,
                Path(PurePath(DEREP_FILE_NAME)).name,
            ]
        )
        to_return[0] = zotu_tab_file

    # ##########################################################################################
    # #################################### OTU Pipline #########################################
    # ##########################################################################################
    if run_otu_pipeline:
        try:
            ANA_LOG.info('# Starting OTU pipeline')
            otu_tab_file, non_human_otu_seq_file = otu_pipeline(
                sorted_file,
                DEREP_FILE_NAME,  # contains sequences with sample ids and corresponding sizes
                ARGS_CLS.denovo_cluster_otus.abund_limit,
                ARGS_CLS.denovo_cluster_otus.sample_wise_correction,
                ANA_LOG,
                skip_non_bacterial_filter,  # gets defined from the CLI.
                ARGS_CLS.denovo_cluster_otus.minsize
            )
            # Normalizing Tables
            otu_norm_methods = "No"
            ANA_LOG.info('# Normalizing OTU Table')
            otu_norm_methods = normalize_otu_table(otu_tab_file, str(parsed_spike_stat_file_path))
            ANA_LOG.info(f"{otu_norm_methods} normalization method(s) applied on {otu_tab_file}")
        except Exception as exc:
            msg = f"OTU pipline failes: {exc}"
            ANA_LOG.error(msg)
            raise exc
        else:
            # preserving important files for next steps
            # and cleaning up the rest
            cleanup(
                str(ANALYSIS_DIR),
                [
                    non_human_zotus_seq_file,
                    DEREP_FILE_NAME,
                    zotu_tab_file
                ]
            )
            to_return[1] = otu_tab_file

    if run_tac_pipeline:
        try:
            ANA_LOG.info('# Starting TAC pipeline')
            tac_otu_table, tac_otu_seq_file = tac_pipeline(
                zotus_seq_fasta=non_human_zotus_seq_file,
                zotus_table_file=zotu_tab_file,
                zotus_tree_file="ZOTUs-Tree-nj.tre",
                cluster_id=ARGS_CLS.tac_cluster.cluster_thr,
                tac_logger=ANA_LOG
            )
            # Normalizing Tables
            otu_norm_methods = "No"
            ANA_LOG.info('# Normalizing TAC OTU Table')
            otu_norm_methods = normalize_otu_table(tac_otu_table, str(parsed_spike_stat_file_path))
            ANA_LOG.info(f"{otu_norm_methods} normalization method(s) applied on {tac_otu_table}")
        except Exception as exc:
            msg = f"TAC pipeline failed: {exc}"
            ANA_LOG.error(msg)
            raise exc
        else:
            # preserving important files for next steps
            # and cleaning up the rest
            cleanup(
                str(ANALYSIS_DIR),
                [
                    non_human_zotus_seq_file,
                    DEREP_FILE_NAME,
                ]
            )
            to_return[2] = tac_otu_table

    if run_tic_pipeline:
        try:
            ANA_LOG.info('# Starting TIC pipeline')
            tic_otu_table, sotu_fasta_path = tic_pipeline(
                zotus_seq_fasta=non_human_zotus_seq_file,
                zotus_table_file=zotu_tab_file,
                zotus_tree_file="ZOTUs-Tree-nj.tre",
                species_id=ARGS_CLS.complex_tic.species_sim,
                genus_id=ARGS_CLS.complex_tic.genus_sim,
                family_id=ARGS_CLS.complex_tic.family_sim,
                tic_logger=ANA_LOG
            )
            # Normalizing Tables
            otu_norm_methods = "No"
            ANA_LOG.info('# Normalizing SOTU Table')
            otu_norm_methods = normalize_otu_table(tic_otu_table, str(parsed_spike_stat_file_path))
            ANA_LOG.info(f"{otu_norm_methods} normalization method(s) applied on {tic_otu_table}")
        except Exception as exc:
            msg = f"TIC pipeline failed: {exc}"
            ANA_LOG.error(msg)
            raise exc
        else:
            # preserving important files for next steps
            # and cleaning up the rest
            cleanup(
                str(ANALYSIS_DIR),
                [
                    non_human_zotus_seq_file,
                    DEREP_FILE_NAME,
                ]
            )
            to_return[3] = tic_otu_table

    ANA_LOG.info('ANALYSIS DONE')
    # to keep the same return value as the main function
    ANA_LOG.info('# Cleaning up')
    cleanup(str(ANALYSIS_DIR))
    return tuple(to_return)
