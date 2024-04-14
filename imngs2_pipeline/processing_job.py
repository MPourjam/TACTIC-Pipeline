from os import (chdir, mkdir, listdir, remove,
                path, getcwd
                )
from collections import Counter
from statistics import stdev, mean
from .processing_helper import TaskPickle
from .processing_helper import (
    IMNGS2ArgsParser,
    gimmelogger,
    loud_subprocess,
    calc_covered_region,
    calc_spikes)
import re
import random
import shutil
import glob
import logging
import zipfile


BIN_DIR = "/base/binaries/"
USEARCH_8_BIN = BIN_DIR + "usearch8.1"
USEARCH_11_BIN = BIN_DIR + "usearch_11_64"
USEARCH_11_BIN = USEARCH_11_BIN + " -strand both"
SORT_ME_RNA_BIN = BIN_DIR + 'sortmerna'
USEARCH8_1 = USEARCH_8_BIN + " -threads 1"
SINA_BIN = BIN_DIR + 'sina/sina'
DBS_DIR = "/base/databases/"
GOLD_REFDB_USEARCH = DBS_DIR + 'SILVA-bac-16s-90.udb'
SINA_ARB = DBS_DIR + 'SILVA_138.1_SSURef_NR99_12_06_20_opt.arb'
ref16RNAdb_1 = DBS_DIR + "silva-bac-16s-id90.fasta"
ref16RNAdb_2 = DBS_DIR + "silva-arc-16s-id95.fasta"
USEARCH_TAIL = '> /dev/null 2>&1'
bowtie2 = BIN_DIR + "bowtie2/bowtie2"
krona_importtext = BIN_DIR + "Krona/KronaTools/scripts/ImportText.pl"
SPIKESIDX = "/base/spikesidx/spike"
R_processing_stat = "/base/imngs2_pipeline/processing_stats.R"
S_FLAT_LOCATION = '/base/s_flat.txt'


# overwriting system to run it with subprocess.run
def system_sub(cmd_args_list: list, force_log: bool = False, shell: bool = False, capture_output: bool = False, quiet: bool = False):
    run_output, cmd_list = loud_subprocess(cmd_args_list, shell_bool=shell, cap_output=capture_output)
    # logging
    if force_log:
        msg = f"COMMAND: {' '.join(cmd_list)}\n\n"
        PREPPROC_LOG.info(msg)
    if run_output.returncode == 137:  # Process killed due to memory limit
        err_msg = f"Command {' '.join(cmd_list)} exceeded memory limit."
        err_msg += "\n\tIf you are using usearch 32-bit version, consider upgrading to 64-bit version."
        err_msg += "\n\tIf you are using usearch 64-bit then run the programm with lower number of threads."
        raise MemoryError(err_msg)
    elif run_output.returncode != 0 and not quiet:
        raise Exception(f"Error in running command {' '.join(cmd_list)}:\n{run_output.stderr}")
    elif run_output.returncode != 0 and quiet:
        PREPPROC_LOG.error(f"Error in running command {' '.join(cmd_list)}:\n{run_output.stderr}")
    elif run_output.returncode == 0 and run_output.stderr and not quiet:
        PREPPROC_LOG.warning(f"Warning in running command {' '.join(cmd_list)}:\n{run_output.stderr}")
    return run_output


def seqFileStats(seqFileName):
    n_first_seqs = 10000
    assess_coeff = 0.5
    sample_range = round(assess_coeff * n_first_seqs)
    to_assess = random.sample(range(1, n_first_seqs), sample_range)
    seq_counter = 0
    random_length = []
    with open(seqFileName, 'r+') as fastqfile:
        line = fastqfile.readline()[:-1]
        lc = 0
        if not line.startswith("@"):
            raise TypeError("Faulty fastq file: {}".format(seqFileName))
        while line and seq_counter <= n_first_seqs:
            mod = lc % 4
            if mod == 1:
                seq_counter += 1
                if seq_counter in to_assess:
                    random_length.append(len(line))
            lc += 1
            line = fastqfile.readline()[:-1]
    len_mean = mean(random_length)
    len_stdev = int(stdev(random_length))
    return (len_mean, len_stdev)


def read_file(filename):
    with open(filename) as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    return content


def unzip():
    shutil.unpack_archive("upload.zip")


def clean_FastQC(input_file, reverse_file=False):
    zips = sorted(glob.glob('*.zip'))
    for i, zipf in enumerate(zips):
        dir_name = zipf.replace(".zip", "")
        # cmd_0 = 'unzip {} > /dev/null 2>&1'.format(zipf)
        # cmd_1 = 'mv ' + dir_name + '/Images/per_base_quality.png ../R{}-per_base_quality.png'.format(i + 1)
        # cmd_2 = 'rm -r {} {} {} > /dev/null 2>&1'.format(dir_name, "*.html", zipf)
        # system(cmd_0)
        # system(cmd_1)
        # system(cmd_2)
        try:
            shutil.unpack_archive(zipf)
            shutil.move(dir_name + "/Images/per_base_quality.png", "../R{}-per_base_quality.png".format(i + 1))
            shutil.rmtree(dir_name)
            remove(zipf)
        except BaseException:
            pass
    # system("cd ../ && rm -r fastqc_output > /dev/null 2>&1")
    shutil.rmtree("../fastqc_output", ignore_errors=True)


def mymkdir(input_dir):
    try:
        mkdir(input_dir)
    except BaseException:
        pass


def merge_pairs(forward_file, reverse_file):
    # cmd_0 = USEARCH_11_BIN + " -fastq_mergepairs " + forward_file + ' -reverse ' + reverse_file
    # cmd_1 = " -fastq_maxdiffs " + str(ARGS_CLS.merge_pairs.fasq_maxdiffs)
    # cmd_1 += " -fastq_pctid " + str(ARGS_CLS.merge_pairs.fastq_pctid) + " -fastqout merged.fastq"
    # cmd_2 = " -fastq_minmergelen " + str(ARGS_CLS.merge_pairs.fastq_minmergelen)
    # cmd_2 += " -fastq_maxmergelen " + str(ARGS_CLS.merge_pairs.fastq_maxmergelen) + " >/dev/null 2>/dev/null"
    # print(cmd_0 + cmd_1 + cmd_2)
    # system_sub(cmd_0 + cmd_1 + cmd_2)
    cmd_to_call_list = [
        *list(USEARCH_11_BIN.split(" ")),
        "-fastq_mergepairs",
        forward_file,
        "-reverse",
        reverse_file,
        "-fastq_maxdiffs",
        str(ARGS_CLS.merge_pairs.fasq_maxdiffs),
        "-fastq_pctid",
        str(ARGS_CLS.merge_pairs.fastq_pctid),
        "-fastqout",
        "merged.fastq",
        "-fastq_minmergelen",
        str(ARGS_CLS.merge_pairs.fastq_minmergelen),
        "-fastq_maxmergelen",
        str(ARGS_CLS.merge_pairs.fastq_maxmergelen),
    ]
    system_sub(cmd_to_call_list, force_log=True)


def trim_sides():
    # cmd_0 = USEARCH_11_BIN + " -fastx_truncate merged.fastq -stripright " + str(ARGS_CLS.trim_both_sides.stripright)
    # cmd_1 = " -stripleft " + str(ARGS_CLS.trim_both_sides.stripleft) + " -fastqout filtered1.fastq >/dev/null 2>/dev/null"
    # system_sub(cmd_0 + cmd_1)
    system_sub(
        [
            *list(USEARCH_11_BIN.split(" ")),
            "-fastx_truncate",
            "merged.fastq",
            "-stripright",
            str(ARGS_CLS.trim_both_sides.stripright),
            "-stripleft",
            str(ARGS_CLS.trim_both_sides.stripleft),
            "-fastqout",
            "filtered1.fastq"
        ],
        force_log=True
    )


def filter_merged_reads():
    # cmd_0 = USEARCH_11_BIN + " -fastq_filter filtered1.fastq -fastq_maxee_rate " + str(ARGS_CLS.filter_merged.fastq_maxee_rate)
    # cmd_1 = " -fastaout filtered2.fasta >/dev/null 2>/dev/null"
    # system_sub(cmd_0 + cmd_1)
    system_sub(
        [
            *list(USEARCH_11_BIN.split(" ")),
            "-fastq_filter",
            "filtered1.fastq",
            "-fastq_maxee_rate",
            str(ARGS_CLS.filter_merged.fastq_maxee_rate),
            "-fastaout",
            "filtered2.fasta"
        ],
        force_log=True
    )
    system_sub(
        [
            "rm",
            "filtered1.fastq"
        ]
    )


def run_FastQC(forward_file, reverse_file):
    out_dir = 'fastqc_output'
    mymkdir(out_dir)
    # cmd_0 = "fastqc " + forward_file
    # cmd_1 = " --outdir " + out_dir + " --format fastq --threads 1 --quiet"
    # if not reverse_file == '':
    #     cmd_0 = cmd_0 + " " + reverse_file
    try:
        system_sub(
            [
                "fastqc",
                forward_file,
                reverse_file if reverse_file else "",
                "--outdir",
                out_dir,
                "--format",
                "fastq",
                "--threads",
                "1",
                "--quiet"
            ]
        )
        # system(cmd_0 + cmd_1)
    except BaseException:
        raise RuntimeError("FastQC did NOT run")
    chdir(out_dir)
    if not reverse_file == '':
        clean_FastQC(input_file=forward_file, reverse_file=reverse_file)
    else:
        clean_FastQC(input_file=forward_file)


def dereplicate_seqs():
    # SIZE_ARGS = "-sizein " if ARGS_CLS.dereplication.sizein else ""
    # SIZE_ARGS += "-sizeout " if ARGS_CLS.dereplication.sizeout else ""
    # cmd_0 = USEARCH_11_BIN + " -fastx_uniques filtered2.fasta -fastaout derep.fasta "
    # cmd_1 = SIZE_ARGS + USEARCH_TAIL
    # system_sub(cmd_0 + cmd_1)
    system_sub(
        [
            *list(USEARCH_11_BIN.split(" ")),
            "-fastx_uniques",
            "filtered2.fasta",
            "-fastaout",
            "derep.fasta",
            "-sizein",
            "-sizeout"
        ],
        force_log=True
    )
    line_n = 0
    with open('derep.fasta') as derep:
        line = derep.readline()
        while line:
            if line.startswith(">"):
                line_n += 1
            line = derep.readline()
    return line_n


def sort_seqs():
    # cmd_0 = USEARCH_11_BIN + ' -sortbysize derep.fasta -fastaout sorted.fasta'
    # cmd_1 = ' ' + USEARCH_TAIL
    # system_sub(cmd_0 + cmd_1)
    system_sub(
        [
            *list(USEARCH_11_BIN.split(" ")),
            "-sortbysize",
            "derep.fasta",
            "-fastaout",
            "sorted.fasta"
        ],
        capture_output=False,
        force_log=True
    )


def trim_one_side(forward_file):
    # cmd_0 = USEARCH_11_BIN + " -fastx_truncate " + forward_file + " -stripleft " + str(ARGS_CLS.trim_one_side.stripleft)
    # cmd_1 = " -fastqout filtered1.fastq >/dev/null 2>/dev/null"
    # system_sub(cmd_0 + cmd_1)
    system_sub(
        [
            *list(USEARCH_11_BIN.split(" ")),
            "-fastx_truncate",
            forward_file,
            "-stripleft",
            str(ARGS_CLS.trim_one_side.stripleft),
            "-fastqout",
            "filtered1.fastq"
        ],
        force_log=True
    )


def filter_merged_one_side(forward_file):
    curr_mean, curr_sd = seqFileStats(forward_file)
    minLength = curr_mean - int(0.1 * curr_mean) - 5  # remove the primer triming size plus 10% of the mean size
    # cmd_0 = USEARCH_11_BIN + " -fastq_filter filtered1.fastq -fastq_truncqual " + str(ARGS_CLS.filter_single_reads.fastq_truncqual)
    # cmd_1 = " -fastq_maxee_rate " + str(ARGS_CLS.filter_single_reads.fastq_maxee_rate) + " -fastq_trunclen " + str(minLength)
    # cmd_2 = " -fastaout filtered2.fasta >/dev/null 2>/dev/null"
    # system_sub(cmd_0 + cmd_1 + cmd_2)
    system_sub(
        [
            *list(USEARCH_11_BIN.split(" ")),
            "-fastq_filter",
            "filtered1.fastq",
            "-fastq_truncqual",
            str(ARGS_CLS.filter_single_reads.fastq_truncqual),
            "-fastq_maxee_rate",
            str(ARGS_CLS.filter_single_reads.fastq_maxee_rate),
            "-fastq_trunclen",
            str(minLength),
            "-fastaout",
            "filtered2.fasta"
        ],
        force_log=True
    )


# cluster sequences to OTUs
def clusterZOTUs():
    # cmd_part_0 = USEARCH_11_BIN + " -unoise3 sorted.fasta -minsize " + str(ARGS_CLS.cluster_zotus.minsize) + " -zotus zotus.fasta"
    # cmd_part_1 = ' -tabbedout denoising.tab'
    # cmd_part_2 = ' ' + USEARCH_TAIL
    # # clusering of seq in OTUs (clustered)
    # system_sub(cmd_part_0 + cmd_part_1 + cmd_part_2)
    system_sub(
        [
            *list(USEARCH_11_BIN.split(" ")),
            "-unoise3",
            "sorted.fasta",
            "-minsize",
            str(ARGS_CLS.cluster_zotus.minsize),
            "-zotus",
            "zotus.fasta",
            "-tabbedout",
            "denoising.tab"
        ],
        force_log=True
    )


# Filter non 16S sequences
def filter16S():
    # cmd_0 = SORT_ME_RNA_BIN + " --ref " + ref16RNAdb_1 + " --ref " + ref16RNAdb_2 + " --reads "
    # cmd_1 = " zotus.fasta --fastx good-ZOTUs --other " + str(ARGS_CLS.filter_16S.other) + " --workdir ."
    # cmd_2 = " -e " + str(ARGS_CLS.filter_16S.e) + " --num_alignments " + str(ARGS_CLS.filter_16S.num_alignments)
    # cmd_2 += " >/dev/null 2>&1"
    # # run a RNA filtering step
    # # (The program currently do not distinquish between 16S and 18S)
    # system_sub(cmd_0 + cmd_1 + cmd_2)
    # system_sub('mv out/aligned.fasta good_ZOTUs.fa')
    # system('rm -r idx out kvdb')
    system_sub(
        [
            SORT_ME_RNA_BIN,
            "--ref",
            ref16RNAdb_1,
            "--ref",
            ref16RNAdb_2,
            "--reads",
            "zotus.fasta",
            "--fastx",
            "good-ZOTUs",
            "--other",
            str(ARGS_CLS.filter_16S.other),
            "--workdir",
            ".",
            "-e",
            str(ARGS_CLS.filter_16S.e),
            "--num_alignments",
            str(ARGS_CLS.filter_16S.num_alignments)
        ],
        force_log=True
    )
    system_sub(
        [
            "mv",
            "out/aligned.fasta",
            "good_ZOTUs.fa"
        ]
    )


def prepare_zotus():
    contents = read_file('good_ZOTUs.fa')
    out_file = open('ZOTUs.fasta', 'w+')
    for line in contents:
        if line[0] == '>':
            new_line = line + ';size=1;\n'
            out_file.write(new_line)
        else:
            out_file.write(line + '\n')
    out_file.close()


def build_ZOTU_table():
    # cmd_0 = USEARCH_11_BIN + ' -otutab derep.fasta -zotus good_ZOTUs.fa'
    # cmd_1 = " -otutabout zotu_table.txt -id " + str(ARGS_CLS.build_zotus_table.id) + " "
    # # cmd_1 += " -mapout map.txt -notmatched unmapped.fa -dbmatched otus_with_sizes.fa -sizeout "  # for debugging
    # system_sub(cmd_0 + cmd_1 + USEARCH_TAIL)
    system_sub(
        [
            *list(USEARCH_11_BIN.split(" ")),
            "-otutab",
            "derep.fasta",
            "-zotus",
            "good_ZOTUs.fa",
            "-otutabout",
            "zotu_table.txt",
            "-id",
            str(ARGS_CLS.build_zotus_table.id)
        ],
        force_log=True)


def get_sample_size(input_file):
    contents = read_file(input_file)[1:]
    tot_size = 0
    for line in contents:
        curr_size = int(line.split('\t')[1])
        tot_size += curr_size
    return tot_size


def filter_zotu_abundance():
    unf_tot_size = get_sample_size('zotu_table.txt')
    contents = read_file('zotu_table.txt')[1:]
    out_file = open('zotu_table_filtered.txt', 'w+')
    out_file_2 = open('abundant_zotus_table.txt', 'w+')
    for line in contents:
        curr_size = float(line.split('\t')[1])
        if round(float(curr_size / unf_tot_size), 4) >= round(0.0, 4):
            out_file.write(line + '\n')
            out_file_2.write(line.split('\t')[0] + '\n')  # only the zotu ID is written to this file
    out_file.close()
    out_file_2.close()


def select_zotu_seqs():
    # cmd_0 = USEARCH_11_BIN + ' -fastx_getseqs good_ZOTUs.fa -labels abundant_zotus_table.txt'
    # cmd_1 = ' -fastaout ZOTUs-Seqs.fasta '
    # system_sub(cmd_0 + cmd_1 + USEARCH_TAIL)
    system_sub(
        [
            *list(USEARCH_11_BIN.split(" ")),
            "-fastx_getseqs",
            "good_ZOTUs.fa",
            "-labels",
            "abundant_zotus_table.txt",
            "-fastaout",
            "ZOTUs-Seqs.fasta"
        ],
        force_log=True
    )


def addTax(input_id):
    """
    NOTE There is a bug in SINA align by which it automatically cut fasta headers by space and
    take only the first part before space (" ")
    """
    filebasename = 'aligned_' + str(input_id)
    # cmd_part_0 = SINA_BIN + ' --in ZOTUs-Seqs.fasta --search --meta-fmt csv '
    # cmd_part_1 = '--threads 1 --lca-fields tax_slv --turn all '
    # cmd_part_2 = '--db ' + SINA_ARB + ' --out ' + filebasename + '.fasta'
    # cmd_part_3 = ''  # >/dev/null 2>/dev/null'  # for better logging
    # # print(cmd_part_0 + cmd_part_1 + cmd_part_2 + cmd_part_3)
    # system_sub(cmd_part_0 + cmd_part_1 + cmd_part_2 + cmd_part_3)
    system_sub(
        [
            SINA_BIN,
            "--in",
            "ZOTUs-Seqs.fasta",
            "--search",
            "--meta-fmt",
            "csv",
            "--threads",
            "1",
            "--lca-fields",
            "tax_slv",
            "--turn",
            "all",
            "--db",
            SINA_ARB,
            "--out",
            filebasename + ".fasta"
        ],
        force_log=True,
        capture_output=True,
        quiet=True
    )

    out_file = open('classifiedF.txt', 'w+')
    silva_contents_header = read_file(filebasename + '.csv')
    silva_contents = silva_contents_header[1:]
    for line in silva_contents:
        line_tokens = line.split(',')
        silva_taxo = line_tokens[6]
        seq_name = line_tokens[0]
        out_file.write(seq_name + '\t' + silva_taxo + '\n')
    out_file.close()


def create_final_ZOTU_table():
    out_file = open('ZOTUs-table.final.tab', 'w+')
    out_file.write('#ZOTUId\tsize\ttaxonomy\n')
    taxa_file = read_file('classifiedF.txt')
    taxa_dict = {}
    for line in taxa_file:
        OTU, taxonomy = line.split('\t')
        taxa_dict[OTU] = taxonomy
    table_file = read_file('zotu_table_filtered.txt')
    for line in table_file:
        OTU, size = line.split('\t')
        OTU_taxonomy = taxa_dict[OTU]
        write_line = OTU + '\t' + size + '\t' + OTU_taxonomy + '\n'
        out_file.write(write_line)
    out_file.close()


def add_taxonomy_to_fasta():
    out_file = open('taxed_ZOTUs.fasta', 'w+')
    table_file = read_file('ZOTUs-table.final.tab')[1:]
    OTUS = list()
    sizes = list()
    taxonomies = list()
    for line in table_file:
        OTU, size, taxonomy = line.split('\t')
        OTUS.append(OTU)
        sizes.append(size)
        taxonomies.append(taxonomy)
    fasta_file = read_file('good_ZOTUs.fa')
    for i in range(0, len(fasta_file), 2):
        header = fasta_file[i]
        sequence = fasta_file[i + 1]
        OTU_found = header[1:]
        try:
            found_index = OTUS.index(OTU_found)
        except BaseException:
            continue
        size_line = ';size=' + sizes[found_index]
        taxa_line = ';tax=' + taxonomies[found_index]
        new_header = header + size_line + taxa_line
        out_file.write(new_header + '\n')
        out_file.write(sequence + '\n')
    out_file.close()


def addKrona(KRONA_TOOL):
    OTUfinal_line = read_file('ZOTUs-table.final.tab')[1:]
    krona_text = open('krona.txt', 'w+')
    krona_dct = {}
    for OTU_taxonomy in OTUfinal_line:
        OTUFinal_list = OTU_taxonomy.split('\t')[1:]
        taxo = '\t'.join(OTUFinal_list[-1].split(";")[:-1])
        size = OTUFinal_list[:-1]
        if not krona_dct.get(taxo, False):
            krona_dct[taxo] = sum([int(s) for s in size])
        else:
            krona_dct[taxo] = krona_dct[taxo] + sum([int(s) for s in size])
    for tax, size in krona_dct.items():
        krona_text.write(str(size) + '\t' + tax + '\n')
    krona_text.close()
    # system("{} ./krona.txt".format(KRONA_TOOL))
    system_sub(
        [
            KRONA_TOOL,
            "./krona.txt"
        ],
        capture_output=True,
        quiet=True
    )
    # system("rm krona.txt")
    system_sub(
        [
            "rm",
            "krona.txt"
        ],
        capture_output=True,
        quiet=True
    )


def create_zip(input_id):
    files_to_zip = [
        "ZOTUs-table.final.tab",
        "taxed_ZOTUs.fasta",
        "spike_stat_mapping_file.csv",
        "krona.html",
        "silva_start_end.txt",
        "reads_report.txt",
        "IMNGS2Pipeline.log"
    ]

    with zipfile.ZipFile(f"../{input_id}_processed.zip", "w") as zipf:
        for file_tozip in files_to_zip:
            # try to write the file to the zip and if it fails, ignore it
            try:
                zipf.write(file_tozip)
            except FileNotFoundError:
                pass

        for file in glob.glob("*.png"):
            zipf.write(file)
        for file in glob.glob("*.html"):
            zipf.write(file)
    shutil.make_archive(f"aligned_{input_id}.fasta", "gztar", ".", f"aligned_{input_id}.fasta")
    shutil.make_archive("derep.fasta", "gztar", ".", "derep.fasta")


def create_udb(input_id):
    # cmd_0 = USEARCH_11_BIN + ' -makeudb_ublast taxed_ZOTUs.fasta -output ' + str(input_id) + '.udb '
    # system_sub(cmd_0 + USEARCH_TAIL)
    system_sub(
        [
            *list(USEARCH_11_BIN.split(" ")),
            "-makeudb_ublast",
            "taxed_ZOTUs.fasta",
            "-output",
            f"{input_id}.udb"
        ],
        force_log=True
    )


def cleanup(input_id, full_clean=False, dir_path=None):
    """
    NOTE Be careful where you run this command.
    """
    if not dir_path:
        dir_path = getcwd()
    chdir(dir_path)
    if full_clean:
        files = [f for f in listdir(dir_path) if not f.endswith("_logs.txt")]
        for file in files:
            file_path = path.join(dir_path, file)
            if path.isfile(file_path):
                remove(file_path)
            elif path.isdir(file_path):
                shutil.rmtree(str(file_path))
    else:
        deleted_files = [
            'zotu_table.txt',
            'good_ZOTUs.fa',
            'test.csv',
            'filtered1.fastq',
            'filtered2.fasta',
            f'aligned_{input_id}.csv',
            f'aligned_{input_id}.fasta',
            'merged.fastq',
            'sorted.fasta',
            'zotus.fasta',
            'otus1.fa',
            'z2o.tab',
            'mOTUs-Seqs.fasta',
            'ZOTUs-Table.tab',
            'derep.fasta',
            'classifiedF.txt',
            'abundant_zotus_table.txt',
            'denoising.tab',
            'matched_ZOTUS.txt',
            'ZOTUs.fasta',
            'ZOTUs-Seqs.fasta',
            'zotu_table_filtered.txt',
            'nochi-ZOTUs.fasta',
        ]
        deleted_dirs = [
            'kvdb',
            'out',
            'idx',
            'fastqc_output',
            'spike_result'
        ]
        for entry in listdir(dir_path):
            entry_path = path.join(dir_path, entry)
            if (entry in deleted_files or entry.endswith(".fastq")) and path.isfile(entry_path):
                remove(str(entry_path))
            elif entry in deleted_dirs and path.isdir(entry_path):
                shutil.rmtree(str(entry_path))


def update_s_flat(input_id, origin):
    phrase = str(input_id) + ',' + origin + ','
    table_file = read_file('ZOTUs-table.final.tab')[1:]
    sizes = list()
    for line in table_file:
        size = int(line.split('\t')[1])
        sizes.append(size)
    phrase += str(len(sizes)) + ',' + str(sum(sizes)) + '\n'
    print(S_FLAT_LOCATION)
    out_file = open(S_FLAT_LOCATION, 'a+')
    out_file.write(phrase)
    out_file.close()


def find_silva_start_end(fasta_file, max_seq=200):
    '''
    fasta_file = the path to the alignment fasta output file
    max_seq = the nth first sequence to calculate start and end mode accordingly
    '''
    if path.isfile(path.abspath(fasta_file)):
        filedir, filepath = path.split(path.abspath(fasta_file))
        chdir(filedir)
    else:
        raise FileNotFoundError("Fasta file not found!")
    try:
        max_seq = int(max_seq)
        if max_seq < 1:
            raise ValueError("max_seq must be more than 0")
    except Exception:
        max_seq = 200

    start_reg = re.compile("^[-NACGTU]")
    seq_reg_start = re.compile(r'^-*[NACGTU]')
    seq_reg_end = re.compile(r'[NACGTU]-*$')
    # Getting the index of headers
    s_e = []
    with open(filepath) as f:
        # zotu_names_idx = []
        counter = 1
        line = f.readline()
        while line and counter <= max_seq:  # or counter < max_seq:
            if start_reg.match(line):
                seque = line.strip().upper()
                m_start = seq_reg_start.match(seque)
                m_end = seq_reg_end.search(seque)
                start = m_start.end() if m_start else 1
                end = m_end.start() if m_end else 0
                s_e.append((start - 1, end))
                # zotu_names_idx.append(counter)
                counter += 1
            line = f.readline()

    # Finding the modes
    start_list = [s[0] for s in s_e]
    end_list = [s[1] for s in s_e]
    starts_count_dict = Counter(start_list)
    ends_count_dict = Counter(end_list)
    mode_start_tuple = list(starts_count_dict.items())
    mode_end_tuple = list(ends_count_dict.items())
    start_mode = mode_start_tuple[0][0]
    end_mode = start_mode
    for end_tu in mode_end_tuple:
        if end_tu[0] > start_mode:
            end_mode = end_tu[0]
    if start_mode >= end_mode:
        start_mode = 0
        end_mode = 0

    return start_mode, end_mode


def update_task_pko(status, msg, pk_dir=None) -> bool:
    '''
    It takes a run status and a message and updates the TaskPickle of
    current process.
    status would be ["Started", "Progress", "Error", "Done"]
    '''
    try:
        latest_pk = pko  # pko is available in global scope
        st_code = latest_pk.scode_d.get(status, "")
        latest_pk.task_dict["status"]["run"] = st_code
        latest_pk.task_dict["status"]["msg"] = msg
        if status == "Error":
            latest_pk.task_dict["args"]["forward_file"] = ""
            latest_pk.task_dict["args"]["reverse_file"] = ""
        latest_pk.write()
        return True
    except Exception as e:
        latest_pk.task_dict["status"]["msg"] = "Could not update the status: " + str(e)
        # We do not change the run status as lack of logging is better not to stop the process
        return False


class Filehanlder_pk(logging.FileHandler):
    '''
    Regular file handler with option to update pk file
    '''
    def emit(self, record):
        super().emit(record)
        update_task_pko("Progress", record.msg)


def write_reads_report(input_id, **kwargs):
    initial_report = "{}\n".format(str(input_id))
    for i, v in enumerate(kwargs.items()):
        k_v = "{}:{}".format(v[0], v[1])
        initial_report += k_v
        if i != len(kwargs) - 1:
            initial_report += "\t"
    with open('reads_report.txt', 'w+') as read_rep:
        read_rep.write(initial_report + "\n")
    return initial_report


def main_processing(
        input_dir,
        forward_file: str,
        reverse_file: str,
        input_id: str,
        args_file_path: str = "",
        spike_amount: int = 0,
        usearch_11_bin: str = USEARCH_11_BIN):
    global PREPPROC_LOG, USEARCH_11_BIN
    # usearch_11_bin must be a global variable and pointing to a file
    USEARCH_11_BIN = usearch_11_bin + " -strand both"
    logger_file_path = path.join(path.abspath(input_dir), f"{str(input_id)}_logs.txt")
    PREPPROC_LOG = gimmelogger(
        logger_name=f"run_imngs2.preprocessing_{str(input_id)}",
        log_file=logger_file_path,
        only_file=True
    )
    paired = "Yes" if reverse_file else "No"
    forward_file = path.join(input_dir, forward_file) if forward_file else ""
    reverse_file = path.join(input_dir, reverse_file) if reverse_file else ""
    files_paths = [forward_file, reverse_file]
    # We define the related TaskPickle object here and make it avialble globally as we need the
    # TaskPickle file to stay closed while processing
    global pko
    pk_dir = path.abspath(input_dir)
    pk_file = path.join(pk_dir, "{}.pk".format(input_id))
    # Parsing arguments from args_file_path
    # if args_file_path is empty then default values in processing_helper.py
    # are loaded
    global ARGS_CLS
    ARGS_CLS = IMNGS2ArgsParser(config_yaml=args_file_path).preproc_args
    try:
        if not path.exists(pk_file):
            # Creation of TaskPickle instance for the run
            pko = TaskPickle(pk_file)
            pko.task_dict["args"]["input_dir"] = input_dir
            pko.task_dict["args"]["paired"] = paired
            pko.task_dict["args"]["forward_file"] = forward_file
            pko.task_dict["args"]["reverse_file"] = reverse_file
            pko.task_dict["args"]["input_id"] = input_id
            pko.task_dict["args"]["spike_amount"] = spike_amount
            pko.task_dict["status"]["run"] = pko.scode_d["Started"]
            abs_path_forw = path.exists(path.abspath(path.join(input_dir, forward_file)))
            if reverse_file:
                abs_path_reve = path.exists(path.abspath(path.join(input_dir, reverse_file)))
                if path.exists(abs_path_forw) and path.exists(abs_path_reve):
                    pko.task_dict["status"]["download"] = pko.scode_d["Done"]
                else:
                    pko.task_dict["status"]["download"] = pko.scode_d["Error"]
                    pko.write()
                    raise ValueError("Fastq Files Error.")
            else:
                if path.exists(abs_path_forw):
                    pko.task_dict["status"]["download"] = pko.scode_d["Done"]
                else:
                    pko.task_dict["status"]["download"] = pko.scode_d["Error"]
                    pko.write()
                    raise ValueError("Fastq Files Error.")
            # pko.task_dict["status"]["download"] = pko.scode_d["Started"]
            pko.write()
        else:
            pko = TaskPickle(pk_file)
        chdir(input_dir)
        PREPPROC_LOG.info("# Spike Removal: Started")
        real_reads_c, spike_reads_c = calc_spikes(*files_paths, spike_amount=spike_amount)
        PREPPROC_LOG.info("Actual_reads:{}\tSpike_reads:{}".format(str(real_reads_c), str(spike_reads_c)))
        PREPPROC_LOG.info('# FastQC: Started')
        run_FastQC(forward_file, reverse_file)
        chdir(input_dir)  # This is CRUCIAL to be here
        if reverse_file:
            PREPPROC_LOG.info('# Merging Pairs: Started')
            merge_pairs(forward_file, reverse_file)
            PREPPROC_LOG.info('# Trim Sides: Started')
            trim_sides()
            PREPPROC_LOG.info('# Filter Merged Reads: Started')
            filter_merged_reads()
        else:
            PREPPROC_LOG.info('# Trim One Side: Started')
            trim_one_side(forward_file)
            PREPPROC_LOG.info('# Filter Amplicons: Started')
            filter_merged_one_side(forward_file)
        PREPPROC_LOG.info('# Dereplication: Started')
        dereped_read_n = dereplicate_seqs()
        read_report = write_reads_report(input_id,
                                         Dereplicated_reads=dereped_read_n)
        PREPPROC_LOG.debug(read_report)
        PREPPROC_LOG.info('# Sort Sequences: Started')
        sort_seqs()
        PREPPROC_LOG.info('# Cluster ZOTUs: Started')
        clusterZOTUs()
        PREPPROC_LOG.info('# Filter non 16S sequences: Started')
        filter16S()
        # prepare_zotus()  # Adds size=1 to end of zotus header
        PREPPROC_LOG.info('# Build ZOTU Table: Started')
        build_ZOTU_table()
        PREPPROC_LOG.info('# Filter ZOTUs by Abundance: Started')
        filter_zotu_abundance()  # ignore this step because the required abundance is 0
        PREPPROC_LOG.info('# Select ZOTU Sequences: Started')
        select_zotu_seqs()
        PREPPROC_LOG.info('# Add Taxonomy: Started')
        addTax(input_id)
        create_final_ZOTU_table()
        add_taxonomy_to_fasta()
        PREPPROC_LOG.info("# Add Krona Graph: Started")
        addKrona(krona_importtext)
        # system('Rscript {} >/dev/null 2>/dev/null'.format(R_processing_stat))
        system_sub(
            [
                "Rscript",
                R_processing_stat
            ],
            force_log=False,
            capture_output=True,
            quiet=True
        )
        # udb for both similarity queries
        PREPPROC_LOG.info('# Create UDB: Started')
        create_udb(input_id)
        # update_s_flat(input_id, origin)
        start_mode, end_mode = find_silva_start_end('aligned_' + str(input_id) + '.fasta')
        calced_regions = calc_covered_region(start_mode, end_mode)
        # Writing the start and end mode to the a file
        with open('silva_start_end.txt', 'w+') as s_e_file:
            # writing header
            s_e_file.write("SilvaAlignmentStartPos\tSilvaAlignementEndPos\tCoveredRegion\n")
            s_e_file.write(str(start_mode) + '\t' + str(end_mode) + '\t' + str(calced_regions) + '\n')
        PREPPROC_LOG.info('# Zip Up: Started')
        create_zip(input_id)
        PREPPROC_LOG.info('# Clean Up: Started')
        cleanup(input_id)
    except BaseException as e:
        err_msg = str(e).split("]")[-1]  # To exclude possible '[Errno 2]' from the message
        PREPPROC_LOG.error(err_msg)
        update_task_pko("Error", str(err_msg).strip())
        cleanup(input_id, full_clean=True, dir_path=path.abspath(input_dir))
        pko.close()
        raise e

    chdir(input_dir)
    udb_file = path.join(input_dir, '{}.udb'.format(str(input_id)))
    status_code = "Done"
    status_msg = ""
    if not path.isfile(udb_file):
        status_code = "Error"
        status_msg = "Process Failed! no UDB!"
    update_task_pko(status_code, status_msg)
    pko.close()
    # if not path.isfile(udb_file):
    #    raise FileNotFoundError("No UDB created for {}".format(input_dir))

    return input_dir
