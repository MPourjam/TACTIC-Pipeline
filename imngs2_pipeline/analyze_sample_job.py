import os
import re
import shutil
from os import chdir
from multiprocessing import cpu_count
from .TIC.complex_TIC import main_complex_TIC
from .TIC.split_based_on_taxonomy import split_based_on_taxonomy
from .processing_helper import IMNGS2ArgsParser, gimmelogger, MyCounter, loud_subprocess
from .spike_normalizer import normalize_otu_table
from sys import version_info
if version_info[0] < 3:
    from pathlib2 import Path, PurePath  # pip2 install pathlib2
else:
    from pathlib import Path, PurePath


BIN_DIR = "/base/binaries/"
# Here we don't add -strand both to usearch as we have already added in preprocessing pipeline
USEARCH_11_BIN = BIN_DIR + "usearch_11_64"
USEARCH_8_bin = BIN_DIR + "usearch8.1"
DB_LOC = "/base/databases/"
SINA_ARB = DB_LOC + "SILVA_138.1_SSURef_NR99_12_06_20_opt.arb"
GOLD_REFDB_USEARCH = DB_LOC + "SILVA-bac-16s-90.udb"
# USERS_DIR = '/srv/crc/users/'
ref16RNAdb_1 = DB_LOC + "silva-bac-16s-id90.fasta"
ref16RNAdb_2 = DB_LOC + "silva-arc-16s-id95.fasta"
SORT_ME_RNA_BIN = BIN_DIR + "sortmerna"
USEARCH_TAIL = '> /dev/null 2>&1'
global SPIKE_STAT_HEADER, TAXED_ZOTU_FILE_NAME
SPIKE_STAT_FILE_COLS = ("#SampleID", "SpikeReads", "spikes_total_weight_in_g", "spike_amount", "parent_path")
SPIKE_STAT_HEADER = "\t".join(list(SPIKE_STAT_FILE_COLS))
TAXED_ZOTU_FILE_NAME = "taxed_ZOTUs.fasta"
max_pool = int(cpu_count() * 0.7)
global POOL_SIZE
POOL_SIZE = max_pool if max_pool > 0 else 1


def system_sub(cmd_args_list: list, force_log: bool = False, shell: bool = False, capture_output: bool = False, quiet: bool = False):
    run_output, cmd_list = loud_subprocess(cmd_args_list, shell_bool=shell, cap_output=capture_output)
    # logging
    if force_log:
        msg = f"COMMAND: {' '.join(cmd_list)}\n\n"
        # msg += f"STDOUT: {run_output.stdout}\n\n"
        ANA_LOG.info(msg)
    if run_output.stderr and not quiet:
        raise Exception(f"Error in running command {' '.join(cmd_list)}:\n{run_output.stderr}")
    if run_output.stderr and quiet:
        ANA_LOG.warning(f"Warning in running command {' '.join(cmd_list)}:\n{run_output.stderr}")
    return run_output


def read_file(filename):
    with open(filename) as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    return content


def onelinefasta(fastafilepath):
    proper_filepath = os.path.abspath(fastafilepath)
    dirpath, filename = os.path.split(proper_filepath)
    os.chdir(dirpath)
    f = open(proper_filepath, 'r')
    newfile_temp_name = "oneline_{}.fasta".format(str(".".join(filename.split(".")[:-1])))
    assert (not os.path.isfile(newfile_temp_name)), "Destination file {} already exists".format(newfile_temp_name)
    with open(newfile_temp_name, "w+") as onelinefa:
        line = f.readline()
        first = True
        while line:
            if line[0] == '>':
                if first:
                    onelinefa.write(line)
                    first = False
                else:
                    onelinefa.write('\n' + line)
            elif line == '\n':
                pass
            else:
                onelinefa.write(line[:-1])
            line = f.readline()
        onelinefa.write("\n")
    shutil.move(newfile_temp_name, filename)  # handling space in name of file


def append_reads(taxed_ZOTUs_file_path, sample_id, analysis_dir):
    chdir(analysis_dir)
    taxed_ZOTUs_file_path = os.path.abspath(taxed_ZOTUs_file_path)
    dataset_name = f"{sample_id}"
    cmd_to_call_list = [
        USEARCH_11_BIN,
        '-fastx_relabel',
        taxed_ZOTUs_file_path,
        '-prefix',
        f"{dataset_name}.",
        '-fastaout',
        'new',
        '-keep_annots',
    ]
    # cmd = USEARCH_11_BIN + ' -fastx_relabel ' + taxed_ZOTUs_file_path + ' -prefix ' + dataset_name
    # cmd += '. -fastaout new -keep_annots'
    system_sub(cmd_to_call_list, force_log=True)
    with open('new', 'r') as new_file:
        with open('analysis.fasta', 'a') as analysis_file:
            analysis_file.write(new_file.read())
    onelinefasta('analysis.fasta')


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
        'filtered1.fasta'
    ]
    system_sub(cmd_to_call_list, force_log=True)
    # cmd_2 = 'mv filtered1.fasta analysis.fasta'
    cmd_to_call_list = [
        "mv",
        "filtered1.fasta",
        "analysis.fasta",
    ]
    system_sub(cmd_to_call_list)
    onelinefasta('analysis.fasta')


def dereplication():
    # cmd_part_0 = USEARCH_11_BIN + " -fastx_uniques analysis.fasta -relabel Zotu"
    # cmd_part_0 += " -fastaout derep.fasta -tabbedout relabel.tab --threads " + str(POOL_SIZE)
    # cmd_part_0 += " > /dev/null 2>/dev/null"
    # system_sub(cmd_part_0)

    cmd_to_call_list = [
        USEARCH_11_BIN,
        "-fastx_uniques",
        "analysis.fasta",
        "-relabel",
        "Zotu",
        "-fastaout",
        "derep.fasta",
        "-tabbedout",
        "relabel.tab",
        "--threads",
        str(POOL_SIZE),
    ]
    # dereplicate reads and anotate them by multiplicity
    system_sub(cmd_to_call_list, force_log=True)
    onelinefasta("derep.fasta")
    contents = open('relabel.tab', 'r')
    zotus_dict = dict()
    zotus_tax_dict = dict()
    line = contents.readline()
    while line:
        tokens = line.split('\t')
        assigned_zotu = tokens[1]
        if assigned_zotu not in zotus_dict.keys():
            zotus_dict[assigned_zotu] = [tokens[0]]
        else:
            zotus_dict[assigned_zotu] += [tokens[0]]
        # creating taxonomy without last semicolon and without 'unclutured' or 'unknown'
        # removing last empty array element gets created by last ;
        # . tokens[-1][:-1] of relabel.tab which is the tax for that cluster
        tax_arr = str(tokens[-1][:-1].split('tax=')[1]).split(';')[:-1]  # [:-1] to exclude \n at the end
        to_exclude = ['uncultured', 'unknown', 'incertae']
        # NOTE New TIC will also check for removing this unknown levels
        for tax_i in range(len(tax_arr) - 1, -1, -1):
            tax_lev_val = tax_arr[tax_i]
            ex_bool = [True for te in to_exclude if bool(re.search(te, str(tax_lev_val).strip().lower()))]
            if any(ex_bool):
                tax_arr = tax_arr[:tax_i]
        corrected_tax = ";".join(tax_arr)
        zotus_tax_dict[assigned_zotu] = corrected_tax
        line = contents.readline()
    # zotus_dict is line --> {"Zotu1" : ["sample_name1", "sample_name2", ...], ... }
    contents.close()
    samples_list = list()
    sample_name_regex = re.compile(r"(.*)\.([0-9]+);size=([0-9]+);")
    for vals in zotus_dict.values():
        for val in vals:
            curr_sample = sample_name_regex.search(val).group(1)
            if curr_sample not in samples_list:
                samples_list.append(curr_sample)
    # print(zotus_dict)
    zotu_hit_dict = dict()
    for key, values_list in zotus_dict.items():
        # We rely on uniqness of dataset ids by the db
        cols_dict = {str(s): 0 for s in samples_list}
        for hit in values_list:
            sample_hit = sample_name_regex.search(hit).group(1)
            sample_size = int(sample_name_regex.search(hit).group(3))
            cols_dict[str(sample_hit)] += sample_size
        zotu_hit_dict[key] = cols_dict  # {'Zotu1': {"D3": 20, "D4": 0, "D6": 36, ...}, ...}
    zotus_dict.clear()
    # print(zotu_hit_dict)
    #####################
    # Creating ZOTU Map #
    #####################
    out_file = open('ZOTU_map.tab', 'w+')
    curr_line = ''
    for item in samples_list:
        curr_line += '\t' + str(item)
    out_line = '#Zotu' + curr_line + '\t' + '' + '\n'
    out_file.write(out_line)  # This is the header
    for key, value in zotu_hit_dict.items():
        curr_line = ''
        for sample_key in samples_list:
            # print(value[str(sample_key)])
            curr_line += '\t' + str(value[str(sample_key)])
        out_line = key + curr_line + '\t' + zotus_tax_dict[key] + '\n'
        out_file.write(out_line)
    out_file.close()
    # Writing taxonomy to derep.fasta
    with open("derep.fasta", "r") as derep_fio, open("derep_with_tax.fasta", "w+") as derep_tax_fio:
        derep_line = derep_fio.readline()
        while derep_line:
            derep_line = derep_line.strip()
            if derep_line.startswith(">"):
                curr_zotu = str(derep_line[1:]).strip()
                new_line = ">" + curr_zotu + ';tax=' + zotus_tax_dict[curr_zotu] + '\n'
            else:
                new_line = derep_line + "\n"
            derep_tax_fio.write(new_line)
            derep_line = derep_fio.readline()
    cmd_to_call_list = [
        "rm",
        "relabel.tab",
        "analysis.fasta",
    ]
    system_sub(cmd_to_call_list)


def sort_seqs():
    # cmd_0 = USEARCH_11_BIN + ' -sortbysize derep.fasta -fastaout sorted.fasta'
    # cmd_1 = ' ' + USEARCH_TAIL
    # system_sub(cmd_0 + cmd_1)
    cmd_to_call_list = [
        USEARCH_11_BIN,
        "-sortbysize",
        "derep.fasta",
        "-fastaout",
        "sorted.fasta",
    ]
    system_sub(cmd_to_call_list, force_log=True)
    onelinefasta("sorted.fasta")


# cluster sequences to OTUs
def clusterZOTUs():
    # cmd_part_0 = USEARCH_11_BIN + ' -unoise3 sorted.fasta -minsize 4 -zotus zotus.fasta'
    # cmd_part_1 = ' ' + USEARCH_TAIL
    # # clusering of seq in OTUs (clustered)
    # system_sub(cmd_part_0 + cmd_part_1)
    cmd_to_call_list = [
        USEARCH_11_BIN,
        "-unoise3",
        "sorted.fasta",
        "-minsize",
        "4",
        "-zotus",
        "zotus.fasta",
    ]
    system_sub(cmd_to_call_list, force_log=True)
    onelinefasta("zotus.fasta")
    cmd_to_call_list = [
        "rm",
        "sorted.fasta",
        "derep.fasta",
    ]
    system_sub(cmd_to_call_list)


# Filter non 16S sequences
def filter16S():
    # cmd_0 = SORT_ME_RNA_BIN + " --ref " + ref16RNAdb_1 + " --ref " + ref16RNAdb_2 + " --reads "
    # cmd_1 = " zotus.fasta --fastx --aligned good-ZOTUs --other non16SrRNA --workdir . -e 0.1 --num_alignments 1"  # --aligned
    # cmd_2 = " > /dev/null 2>&1"
    # run a RNA filtering step
    # (The program currently do not distinquish between 16S and 18S)
    # system_sub(cmd_0 + cmd_1 + cmd_2)
    # system_sub('mv out/aligned.fasta good-ZOTUs.fasta')
    cmd_to_call_list = [
        SORT_ME_RNA_BIN,
        "--ref",
        ref16RNAdb_1,
        "--ref",
        ref16RNAdb_2,
        "--reads",
        "zotus.fasta",
        "--fastx",
        "--aligned",
        "good-ZOTUs",
        "--other",
        "non16SrRNA",
        "--workdir",
        ".",
        "-e",
        "0.1",
        "--num_alignments",
        "1",
    ]
    system_sub(cmd_to_call_list, force_log=True)
    system_sub(["rm", "-r", "idx", "kvdb"])


def prepare_zotus():
    contents = read_file('good-ZOTUs.fasta')
    out_file = open('good-ZOTUs-sized.fasta', 'w+')
    for line in contents:
        if line[0] == '>':
            new_line = line + ';size=1;\n'
            out_file.write(new_line)
        else:
            out_file.write(line + '\n')
    out_file.close()


# cluster sequences to OTUs
def clusterOTUs():
    # cmd_part_0 = USEARCH_11_BIN + ' -cluster_otus good-ZOTUs-sized.fasta -fulldp -otus otus1.fa -minsize 1 '
    # cmd_part_1 = ' -uparseout z2o.tab'
    # cmd_part_2 = ' > /dev/null 2>/dev/null'
    # clusering of seq in OTUs (clustered)
    # system_sub(cmd_part_0 + cmd_part_1 + cmd_part_2)
    cmd_to_call_list = [
        USEARCH_11_BIN,
        "-cluster_otus",
        "good-ZOTUs-sized.fasta",
        "-fulldp",
        "-otus",
        "otus1.fa",
        "-minsize",
        "1",
        "-uparseout",
        "z2o.tab",
    ]
    system_sub(cmd_to_call_list, force_log=True)
    onelinefasta("otus1.fa")


def remove_size_from_OTUS():
    contents = read_file('otus1.fa')
    out_file = open('mOTUs-Seqs.fasta', 'w+')
    for line in contents:
        if line[0] == '>':
            new_line = line.split(';')[0] + ';\n'
            out_file.write(new_line)
        else:
            out_file.write(line + '\n')
    out_file.close()
    system_sub(["rm", "otus1.fa"])


def assign_zotus_to_otus():
    zotus_map = read_file('z2o.tab')
    out_file = open('ZOTUs-OTUs-map.tab', 'w+')
    out_file_2 = open('matched_ZOTUS.txt', 'w+')
    valid_zotus_dict = dict()
    for line in zotus_map:
        line_tokens = line.split('\t')
        if 'otu' in line_tokens[1]:
            zotu_number = line_tokens[0].split(';')[0]
            valid_zotus_dict[zotu_number] = line_tokens[1]
        elif 'match' == line_tokens[1]:
            zotu_number = line_tokens[0].split(';')[0]
            matched_zotu = line_tokens[2].split(';')[1].split('=')[1]
            matched_otu_of_match = valid_zotus_dict[matched_zotu]
            valid_zotus_dict[zotu_number] = matched_otu_of_match
    for key, value in valid_zotus_dict.items():
        line = key + '\t' + value + '\n'
        out_file.write(line)
    out_file.close()
    for key in valid_zotus_dict.keys():
        out_file_2.write(key + '\n')
    out_file_2.close()
    system_sub(["rm", "z2o.tab"])


def keep_good_ZOTUs():
    # cmd_0 = USEARCH_11_BIN + ' -fastx_getseqs good-ZOTUs.fasta -labels '  # aligned_16S_ZOTUS.fa -> good-ZOTUs.fasta
    # cmd_1 = 'matched_ZOTUS.txt -fastaout nochi_ZOTUs.fasta '
    # system_sub(cmd_0 + cmd_1 + USEARCH_TAIL)
    cmd_to_call_list = [
        USEARCH_11_BIN,
        "-fastx_getseqs",
        "good-ZOTUs.fasta",
        "-labels",
        "matched_ZOTUS.txt",
        "-fastaout",
        "nochi_ZOTUs.fasta",
    ]
    system_sub(cmd_to_call_list, force_log=True)
    onelinefasta("nochi_ZOTUs.fasta")


def build_ZOTU_table():
    # cmd_0 = USEARCH_11_BIN + ' -otutab analysis.fasta -top_hit_only -zotus nochi_ZOTUs.fasta '
    # cmd_1 = ' -otutabout zotu_table.txt -id 0.97'
    # system_sub(cmd_0 + cmd_1 + USEARCH_TAIL)
    cmd_to_call_list = [
        USEARCH_11_BIN,
        "-otutab",
        "analysis.fasta",
        "-top_hit_only",
        "-zotus",
        "nochi_ZOTUs.fasta",
        "-otutabout",
        "zotu_table.txt",
        "-id",
        "0.97"
    ]
    system_sub(cmd_to_call_list, force_log=True)


def get_samples_sizes(input_file):
    contents = read_file(input_file)[1:]
    tot_sizes_list = list()
    num_tokens = len(contents[0].split('\t')) - 2
    # Adding 0.001 to each sum to avoid division by zero
    for i in range(0, num_tokens):
        tot_sizes_list.append(0.001)
    for line in contents:
        line_tokens = line.split('\t')[1:-1]
        for i in range(len(line_tokens)):
            tot_sizes_list[i] += float(line_tokens[i])
    tot_sizes_list = [round(sz, 5) for sz in tot_sizes_list]
    return tot_sizes_list


def filter_zotu_abundance(abund_limit):
    unf_tot_sizes = get_samples_sizes('ZOTU_map.tab')
    contents = read_file('ZOTU_map.tab')
    sample_sizes = list()
    datasets_samples = len(contents[0].split('\t')) - 1
    for i in range(0, datasets_samples):
        sample_sizes.append(0)
    out_file = open('filtered_zotu_table.txt', 'w+')
    out_file.write(contents[0] + '\n')
    # out_file_2 = open('filtered_zotu_table_list.txt', 'w+')
    for line in contents[1:]:
        curr_sizes = line.split('\t')[1:-1]
        curr_abundances = list()
        for j in range(len(curr_sizes)):
            abundance = float(curr_sizes[j]) / unf_tot_sizes[j]
            curr_abundances.append(abundance)
        if any([x >= abund_limit for x in curr_abundances]):
            out_file.write(line + '\n')
    out_file.close()
    system_sub(["mv", "filtered_zotu_table.txt", "ZOTU_map.tab"])
    zotus_tax_dict = dict()
    map_f = open("ZOTU_map.tab", "r")
    map_line = map_f.readline()
    while map_line:
        if map_line.startswith("Zotu"):
            cols = map_line.split("\t")
            zotu_name = cols[0]
            tax = cols[-1][:-1]  # removing "\n"
            zotus_tax_dict[zotu_name] = tax
        map_line = map_f.readline()
    map_f.close()
    fasta_contents = open('derep.fasta', 'r')
    out_file = open('derep_with_tax.fasta', 'w+')
    derep_line = fasta_contents.readline()
    write_seq = True
    while derep_line:
        if derep_line[0] == '>':
            curr_zotu = derep_line[1:-1]  # removing "\n"
            if zotus_tax_dict.get(curr_zotu, False):
                new_line = '>' + curr_zotu + ';tax=' + zotus_tax_dict[curr_zotu] + '\n'
                out_file.write(new_line)
                write_seq = True
            else:
                write_seq = False
        elif write_seq:
            out_file.write(derep_line)
        derep_line = fasta_contents.readline()
    out_file.close()
    fasta_contents.close()


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
        "ZOTUs-Seqs.fasta"
    ]
    system_sub(cmd_to_call_list, force_log=True)
    onelinefasta("ZOTUs-Seqs.fasta")
    system_sub(["rm", "filtered_zotu_table_list.txt"])
    addTax_new("ZOTUs-Seqs.fasta")


def addTax():
    classifier_dir = '/crc/crc/binaries/sina/'
    # cmd_part_0 = 'sina --in ZOTUs-Seqs.fasta --search --meta-fmt csv '
    # cmd_part_1 = f'--threads {POOL_SIZE} --lca-fields tax_slv '
    # cmd_part_2 = '--db ' + SINA_ARB + ' --out test.fasta'
    # cmd_part_3 = ' > /dev/null 2>/dev/null'
    # system_sub(classifier_dir + cmd_part_0 + cmd_part_1 + cmd_part_2 + cmd_part_3)
    cmd_to_call_list = [
        classifier_dir + "sina",
        "--in",
        "ZOTUs-Seqs.fasta",
        "--search",
        "--meta-fmt",
        "csv",
        "--threads",
        str(POOL_SIZE),
        "--lca-fields",
        "tax_slv",
        "--db",
        SINA_ARB,
        "--out",
        "test.fasta"
    ]
    system_sub(cmd_to_call_list, force_log=True)

    out_file = open('classifiedF.txt', 'w+')
    silva_contents_header = read_file('test.csv')
    silva_contents = silva_contents_header[1:]
    for line in silva_contents:
        line_tokens = line.split(',')
        silva_taxo = line_tokens[6]
        seq_name = line_tokens[0]
        out_file.write(seq_name + '\t' + silva_taxo + '\n')
    out_file.close()


def addTax_new(zotu_fasta_path):
    '''
    Input file is a fasta file of seqs
    '''
    zotu_fasta_path = os.path.abspath(zotu_fasta_path)
    dirname, filename = os.path.split(zotu_fasta_path)
    os.chdir(dirname)
    classifier_dir = '/crc/crc/binaries/sina/'
    # cmd_part_0 = 'sina --in ' + str(zotu_fasta_path) + ' --search --meta-fmt csv '
    # cmd_part_1 = f'--threads {POOL_SIZE} --lca-fields tax_slv '
    # cmd_part_2 = '--db ' + SINA_ARB + ' --out test_{}'.format(filename)
    # cmd_part_3 = ' > /dev/null 2>/dev/null'
    # system_sub(classifier_dir + cmd_part_0 + cmd_part_1 + cmd_part_2 + cmd_part_3)
    cmd_to_call_list = [
        classifier_dir + "sina",
        "--in",
        str(zotu_fasta_path),
        "--search",
        "--meta-fmt",
        "csv",
        "--threads",
        str(POOL_SIZE),
        "--lca-fields",
        "tax_slv",
        "--db",
        SINA_ARB,
        "--out",
        "test_{}".format(filename),
    ]
    system_sub(cmd_to_call_list, force_log=True)

    new_taxed_path = os.path.join(dirname, "taxed_{}".format(filename))
    out_file = open(new_taxed_path, 'w+')
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
    for line in silva_contents:
        line_tokens = line.split(',')
        silva_taxo = line_tokens[6]
        seq_name = line_tokens[0]
        out_file.write('>' + seq_name + ';' + "tax=" + silva_taxo + '\n' + zotu_seq_dict[seq_name] + '\n')
    out_file.close()


def create_final_ZOTU_table():
    out_file = open('ZOTUs-table.final.tab', 'w+')
    header = '\t'.join(read_file('zotu_table.txt')[0].split('\t')[1:])
    out_file.write('#ZOTUId\t' + header + '\ttaxonomy\n')
    taxa_file = read_file('classifiedF.txt')
    taxa_dict = {}
    for line in taxa_file:
        ZOTU, taxonomy = line.split('\t')
        taxa_dict[ZOTU] = taxonomy
    table_file = read_file('zotu_table_filtered.txt')
    for line in table_file:
        ZOTU = line.split('\t')[0]
        ZOTU_taxonomy = taxa_dict[ZOTU]
        write_line = line + '\t' + ZOTU_taxonomy + '\n'
        out_file.write(write_line)
    out_file.close()


def add_taxonomy_to_fasta():
    out_file = open('taxed_ZOTUs.fasta', 'w+')
    table_file = read_file('ZOTUs-table.final.tab')[1:]
    OTUS = list()
    taxonomies = list()
    for line in table_file:
        OTU = line.split('\t')[0]
        taxonomy = line.split('\t')[-1]
        OTUS.append(OTU)
        taxonomies.append(taxonomy)
    fasta_file = read_file('ZOTUs-Seqs.fasta')
    for line in fasta_file:
        if line[0] == '>':
            header = line
            OTU_found = header[1:]
            try:
                found_index = OTUS.index(OTU_found)
            except BaseException:
                continue
            taxa_line = ';tax=' + taxonomies[found_index]
            new_header = header + taxa_line
            out_file.write(new_header + '\n')
        else:
            out_file.write(line + '\n')
    out_file.close()


def create_tree(sample):
    sotu_content = read_file(os.path.join(os.getcwd(), "TICOut/sotu_with_taxonomy.fasta"))
    align_content = read_file("test_ZOTUs-Seqs.fasta")
    align_zotus = []
    align_seqs = []
    for ac in align_content:
        if ac[0] == ">":
            zotu_name = ac.split(">")[1]
            align_zotus.append(str(zotu_name))
        else:
            align_seqs.append(str(ac))
    align_dict = dict(zip(align_zotus, align_seqs))
    # print(align_dict.keys())
    with open("filtered_aligned.fasta", "w+") as fas:
        for line in sotu_content:
            if line[0] == ">":
                sotu_name = line.split(";")[-2]
                zotu_name = line.split(";")[0][1:]
                # print(zotu_name)
                seq = align_dict[zotu_name]
                fas.write(">" + str(sotu_name) + "\n")
                fas.write(str(seq) + "\n")
    # create both trees
    # cmd_part_0 = 'FastTree -quiet -nosupport -gtr -nt filtered_aligned.fasta'
    # cmd_part_1 = ' > ' + sample + '_FastTree.tre'
    # system_sub(cmd_part_0 + cmd_part_1)
    cmd_to_call_list = [
        "FastTree",
        "-quiet",
        "-nosupport",
        "-gtr",
        "-nt",
        "filtered_aligned.fasta",
        ">",
        sample + "_FastTree.tre",
    ]
    system_sub(cmd_to_call_list)
    # cmd_2_part_0 = BIN_DIR + 'rapidnj filtered_aligned.fasta --input-format fa --cores 6 --alignment-type d '
    # cmd_2_part_1 = '--output-format t --no-negative-length -x ' + sample + '_rapidNJ_tree.tre'
    # system_sub(cmd_2_part_0 + cmd_2_part_1)
    cmd_to_call_list = [
        BIN_DIR + "rapidnj",
        "filtered_aligned.fasta",
        "--input-format",
        "fa",
        "--cores",
        "6",
        "--alignment-type",
        "d",
        "--output-format",
        "t",
        "--no-negative-length",
        "-x",
        sample + "_rapidNJ_tree.tre",
    ]
    system_sub(cmd_to_call_list)


def cleanup_old():
    to_keep = ["mapping.tab", "OTUs-table.final.tab", "step.png",
               "ZOTUs-table.final.tab", "taxed_ZOTUs-Seqs.fasta", "text.krona.html",
               "denoised_map_sotu.tab", "gotus_to_fotus_map.tab", "krona_plot.html",
               "sotu_sizes.tab", "sotus_to_gotus_map.tab", "sotu_with_taxonomy.fasta",
               "_rapidNJ_tree.tre", "_FastTree.tre"
               ]
    files = os.listdir(os.getcwd())
    to_rem = []
    for f in files:
        abspath_f = os.path.abspath(f)
        if not any([re.search(tk, f) for tk in to_keep]) and os.path.isfile(abspath_f):
            to_rem.append(f)
    cmd_to_call_list = ["rm", "-f"] + to_rem
    system_sub(cmd_to_call_list)


def rename_final_files():
    system_sub("mv zotus_table.tab ZOTUs-table.final.tab".split(" "))
    system_sub("mv sotus_table.tab OTUs-table.final.tab".split(" "))


def create_krona():
    OTUfinal_line = read_file('ZOTUs-table.final.tab')[1:]
    krona_text = open('krona.txt', 'w+')
    for OTU_taxonomy in OTUfinal_line:
        OTUFinal_list = OTU_taxonomy.split('\t')[1:]
        taxo = '\t'.join(OTUFinal_list[-1].split(";")[:-1])
        size = OTUFinal_list[0]
        krona_text.write(size + '\t' + taxo + '\n')
    krona_text.close()
    system_sub("ktImportText krona.txt".split(" "))
    system_sub("rm krona.txt".split(" "))


def create_OTUs_from_ZOTUS_table():
    clean_zotus = read_file('ZOTUs-table.final.tab')
    map_file = read_file('ZOTUs-OTUs-map.tab')
    otus_contents = read_file('mOTUs-Seqs.fasta')
    out_file = open('taxed_OTUs.fasta', 'w+')
    out_file_2 = open('OTUs-table.final.tab', 'w+')
    header = '\t'.join(clean_zotus[0].split('\t')[1:-1])
    out_file_2.write('#OTUId\t' + header + '\ttaxonomy\n')
    zotu_size_dict = dict()
    zotu_taxo_dict = dict()
    zotu_map_otus = dict()
    for line in clean_zotus[1:]:
        zotu_name = line.split('\t')[0]
        taxonomy = line.split('\t')[-1]
        curr_size = '\t'.join(line.split('\t')[1:-1])
        zotu_size_dict[zotu_name] = curr_size
        zotu_taxo_dict[zotu_name] = taxonomy
    for line in map_file:
        zotu_name, otu_name = line.split('\t')
        zotu_map_otus[zotu_name] = otu_name
    curr_keys = list(zotu_map_otus.keys())
    for key in curr_keys:
        if key not in zotu_taxo_dict.keys():
            del zotu_map_otus[key]
    # Writing the new table based on OTUs (summing size of zouts belonging to the same OTU)
    otu_summed_size_dict = dict()
    for line in clean_zotus[1:]:
        zotu_ = line.split('\t')[0]
        otu_ = zotu_map_otus[zotu_]
        sizes_list = zotu_size_dict[zotu_].split("\t")
        curr_size = otu_summed_size_dict.get(otu_, [0 for i in range(0, len(sizes_list))])
        otu_summed_size_dict[otu_] = [str(int(sizes_list[i]) + int(curr_size[i])) for i in range(0, len(curr_size))]

    for line in otus_contents:
        if line[0] == '>':
            zotu_name = line[1:-1]
            if zotu_name in zotu_map_otus.keys():
                # It drops all zotus belonging to an OTU and only takes the centriod
                out_file.write('>' + zotu_map_otus[zotu_name] + ';tax=' + zotu_taxo_dict[zotu_name] + '\n')
                out_line_2 = zotu_map_otus[zotu_name] + '\t'
                out_line_2 += "\t".join(otu_summed_size_dict[zotu_map_otus[zotu_name]]) + '\t'
                out_line_2 += zotu_taxo_dict[zotu_name] + '\n'
                out_file_2.write(out_line_2)
                seq_flag = 1
            else:
                seq_flag = 0
        else:
            if seq_flag:
                out_file.write(line + '\n')
    out_file.close()
    out_file_2.close()
    ###
    testfasta_contents = read_file('test.fasta')
    testfasta_zotus_dict = {l[1:]: testfasta_contents[i + 1] for i, l in enumerate(testfasta_contents) if (l[0] == ">")}
    # otus_contents_zotus = [o[1:] for o in otus_contents if (o[0] == ">")]
    with open('test_otu.fasta', "w+") as test_otu:
        for oc in otus_contents:
            if oc[0] == ">":
                zotu = oc[1:-1]
                if zotu in zotu_map_otus.keys():
                    otu_name = zotu_map_otus[zotu]
                    test_otu.write(">" + otu_name + '\n')
                    write_seq = zotu
                else:
                    write_seq = False
            else:
                if write_seq:
                    test_otu.write(testfasta_zotus_dict[write_seq] + '\n')


def sina_alignment():
    classifier_dir = '/crc/crc/binaries/sina/'
    # cmd_part_0_0 = 'sina --in zotus_with_taxonomy.fasta --search --meta-fmt csv '
    # cmd_part_0_1 = 'sina --in sotus_with_taxonomy.fasta --search --meta-fmt csv '
    # cmd_part_1 = f'--threads {POOL_SIZE} --lca-fields tax_slv '
    # cmd_part_2_0 = '--db ' + SINA_ARB + ' --out test_z.fasta'
    # cmd_part_2_1 = '--db ' + SINA_ARB + ' --out test_s.fasta'
    # cmd_part_3 = ' >/dev/null 2>/dev/null'
    # system_sub(classifier_dir + cmd_part_0_0 + cmd_part_1 + cmd_part_2_0 + cmd_part_3)
    cmd_to_call_list = [
        classifier_dir + "sina",
        "--in",
        "zotus_with_taxonomy.fasta",
        "--search",
        "--meta-fmt",
        "csv",
        "--threads",
        str(POOL_SIZE),
        "--lca-fields",
        "tax_slv",
        "--db",
        SINA_ARB,
        "--out",
        "test_z.fasta",
    ]
    system_sub(cmd_to_call_list, force_log=True)
    # system_sub(classifier_dir + cmd_part_0_1 + cmd_part_1 + cmd_part_2_1 + cmd_part_3)
    cmd_to_call_list = [
        classifier_dir + "sina",
        "--in",
        "sotus_with_taxonomy.fasta",
        "--search",
        "--meta-fmt",
        "csv",
        "--threads",
        str(POOL_SIZE),
        "--lca-fields",
        "tax_slv",
        "--db",
        SINA_ARB,
        "--out",
        "test_s.fasta",
    ]
    system_sub(cmd_to_call_list, force_log=True)


def extract_columns(input_file_name):
    filepath = os.path.abspath(input_file_name)
    dir_, file_ = os.path.split(filepath)
    good_points = list()
    with open(filepath) as fp:
        line = fp.readline()
        while line:
            if not line[0] == '>':
                counter = 0
                for curr_char in line:
                    if curr_char in ['A', 'C', 'G', 'T', 'U']:
                        good_points.append(counter)
                    counter += 1
            line = fp.readline()
    #
    good_points = list(set(good_points))
    #
    out_file = open('for_tree_{}'.format(file_), 'w+')
    with open(filepath) as fp:
        line = fp.readline()
        while line:
            if str(line).startswith('>'):
                header = line.split(' ')[0] + '\n'
                line = fp.readline()
                curr_seq = ''
                for number in good_points:
                    curr_seq += line[number]
                out_file.write(header)
                out_file.write(curr_seq + '\n')
            line = fp.readline()

    out_file.close()


def create_trees():
    # Removing the ; between the taxonomy and the SOTU/ZOTU header file for skipping the
    # the error in Rhea
    for f in ['test_s.fasta', 'test_z.fasta']:
        with open(f, 'r') as test_sfasta:
            with open('corr_header_{}'.format(f), 'w') as new_test_sfasta:
                for ll in test_sfasta:
                    if ll.startswith('>'):
                        lsp = ll.split(';tax=')
                        ll = lsp[0] + ' tax=' + ';'.join(lsp[1].split(';'))
                    new_test_sfasta.write(ll)
        system_sub('mv corr_header_{} {}'.format(f, f))
    for f in ['test_s.fasta', 'test_z.fasta']:
        extract_columns(f)
    # cmd1 = "/crc/crc/binaries/rapidnj for_tree_test_s.fasta -n -o t -x sotu_nj.tre"
    # cmd2 = "/crc/crc/binaries/rapidnj for_tree_test_z.fasta -n -o t -x zotu_nj.tre"
    # cmd3 = "/crc/crc/binaries/FastTree -gtr -gamma -quiet -nt < for_tree_test_s.fasta > sotu_aml.tre"
    # cmd4 = "/crc/crc/binaries/FastTree -gtr -gamma -quiet -nt < for_tree_test_z.fasta > zotu_aml.tre"
    # system_sub(cmd1)
    # system_sub(cmd2)
    # system_sub(cmd3)
    # system_sub(cmd4)
    cmd_to_call_list = [
        BIN_DIR + "rapidnj",
        "for_tree_test_s.fasta",
        "-n",
        "-o",
        "t",
        "-x",
        "sotu_nj.tre",
    ]
    system_sub(cmd_to_call_list)
    cmd_to_call_list = [
        BIN_DIR + "rapidnj",
        "for_tree_test_z.fasta",
        "-n",
        "-o",
        "t",
        "-x",
        "zotu_nj.tre",
    ]
    system_sub(cmd_to_call_list)
    cmd_to_call_list = [
        BIN_DIR + "FastTree",
        "-gtr",
        "-gamma",
        "-quiet",
        "-nt",
        "<",
        "for_tree_test_s.fasta",
        ">",
        "sotu_aml.tre",
    ]
    system_sub(cmd_to_call_list, shell=True)
    cmd_to_call_list = [
        BIN_DIR + "FastTree",
        "-gtr",
        "-gamma",
        "-quiet",
        "-nt",
        "<",
        "for_tree_test_z.fasta",
        ">",
        "zotu_aml.tre",
    ]
    system_sub(cmd_to_call_list, shell=True)


def cleanup(directory: str):
    directory = Path(PurePath(directory)).absolute()
    to_remove = ["new", "derep.fasta", "ZOTU_map.tab", "derep_with_tax.fasta", "zotu_to_sotu_map.tab",
                 "for_tree_test_s.fasta", "for_tree_test_z.fasta", "test_s.fasta", "test_z.fasta",
                 "test_s.csv", "test_z.csv", "log_file.txt"]
    for f in to_remove:
        file_path = directory.joinpath(f)
        if file_path.exists():
            if file_path.is_file():
                file_path.unlink()
            elif file_path.is_dir():
                shutil.rmtree(file_path)


def create_zip():
    system_sub('zip -r ../Analysis.zip .'.split(" "), capture_output=True, shell=True)


def addKrona(KRONA_TOOL):
    OTUfinal_line = read_file('OTUs-table.final.tab')[1:]
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
    cmd_to_call_list = [KRONA_TOOL, "krona.txt"]
    system_sub(cmd_to_call_list)
    system_sub("rm ./krona.txt".split(""))


def fasta_name_with_space(filepath):
    filepath = os.path.abspath(filepath)
    filedir, filename = filepath.split(filepath)
    fout_path = os.path.join(filedir, "placeholder" + filename)
    fopen = open(filepath, 'r')
    fout = open(fout_path, 'w+')
    line = fopen.readline()
    while line:
        if line[0] == ">":
            parts = line.split(';tax=')
            fout.write(parts[0] + ' tax=' + parts[-1])
        else:
            fout.write(line)
        line = fopen.readline()
    fopen.close()
    fout.close()
    cmd_to_call_list = ["mv", fout_path, filepath]
    system_sub(cmd_to_call_list)


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


def main(
        analysis_dir: str,
        spike_stat_file: str,
        fastqs_dir: str,
        args_file_path: str = "",
        dbs_loc: str = DB_LOC,
        threads=POOL_SIZE):
    """
    sample_seq_files_path: a list of file path to
     each samples' sequence file which is going to be combined with other samples passed to analysis.
    analysis_dir: The destination directory to save results
    """
    spike_stat_file = Path(PurePath(spike_stat_file))
    assert spike_stat_file.is_file(), "spike_stat_file must be a path to a file"

    global ANALYSIS_DIR, ARGS_CLS, ANA_LOG, DB_LOC
    # updating POOL_SIZE
    try:
        POOL_SIZE = int(threads)
    except ValueError:
        ANA_LOG.warning(f"threads must be an integer. Using default value of {POOL_SIZE}")
    DB_LOC = Path(PurePath(dbs_loc))
    assert DB_LOC.is_dir(), "DBS_LOC should be path to direcotry containing SILVA database files (arb)"
    ANALYSIS_DIR = Path(PurePath(analysis_dir)).absolute()
    ANA_LOG = gimmelogger(
        logger_name="run_imngs2.analysis",
        log_file=ANALYSIS_DIR.joinpath("Analysis_log.txt"),
        only_file=True
    )
    assert ANALYSIS_DIR.is_dir(), "analysis_dir must be a path to a directory"
    ANALYSIS_DIR = str(ANALYSIS_DIR) + "/"
    ARGS_CLS = IMNGS2ArgsParser(config_yaml=args_file_path).analysis_args
    # Parsing spike_stat_file
    try:
        parsed_spike_stat_file_path = Path(PurePath(ANALYSIS_DIR + 'spike_mapping_file.csv'))
        map_lines_dict, parsed_spike_stat_file_path = parse_spike_stat_file(
            spike_stat_file,
            fastq_dir=fastqs_dir,
            parsed_spike_stat_path=parsed_spike_stat_file_path)
    except Exception as exc:
        raise ValueError(f"spike_stat_file: {str(exc)}")
    sample_seq_files_path = []
    for full_id, entries in map_lines_dict.items():
        sample_seq_files_path.append((Path(PurePath(entries[4])).joinpath(TAXED_ZOTU_FILE_NAME), entries[0]))

    for taxed_path, sam_id in sample_seq_files_path:
        if not taxed_path.is_file():
            raise ValueError(f"{str(taxed_path)} must contain path to each samples sequence file!")

    chdir(ANALYSIS_DIR)
    for sample, sam_id in sample_seq_files_path:
        append_reads(sample, sam_id, ANALYSIS_DIR)
    ANA_LOG.info('READS CONCATENATED')
    chdir(ANALYSIS_DIR)
    trim_sides(ARGS_CLS.trimsides.stripleft, ARGS_CLS.trimsides.stripright)
    dereplication()
    ANA_LOG.info('Dereplication DONE')
    #################
    ## TIC is here ##
    #################
    TIC_Result_DIR = Path(ANALYSIS_DIR).joinpath("TICResult")
    TIC_Output_DIR = Path(ANALYSIS_DIR).joinpath("TICOut")
    dereplicated_fasta_path = Path(ANALYSIS_DIR).joinpath("derep_with_tax.fasta")
    if TIC_Result_DIR.is_dir():
        shutil.rmtree(str(TIC_Result_DIR))
    try:
        split_based_on_taxonomy(str(TIC_Result_DIR), "derep_with_tax.fasta")
    except Exception as exc:
        raise ValueError(f"Split: {str(exc)}")
    try:
        main_complex_TIC(tool=USEARCH_11_BIN, data_dir=str(TIC_Result_DIR), threads=POOL_SIZE)
    except Exception as exc:
        raise ValueError(f"Main: {str(exc)}")
    if os.path.isdir(TIC_Output_DIR):
        shutil.rmtree(TIC_Output_DIR)
    # Abundance filtering
    abund_limit = ARGS_CLS.create_table.abund_limit
    # output files names
    ZOTUs_fasta_name = "ZOTUs-Seqs.fasta"
    SOTUs_fasta_name = "SOTUs-Seqs.fasta"
    ZOTUs_table_name = "ZOTUs-Table.tab"
    SOTUs_table_name = "SOTUs-Table.tab"
    chdir(ANALYSIS_DIR)
    main_create_args = [
        "python3.7",
        "/base/imngs2_pipeline/TIC/create_fasta_and_table.py",  # 0
        str(TIC_Output_DIR),  # OUTPUT_FOLDER
        ZOTUs_fasta_name,  # OUTPUT_ASV_FASTA_WITH_TAXONOMY
        ZOTUs_table_name,  # OUTPUT_ASV_TABLE
        str(TIC_Result_DIR),  # CLUSTERING_DIRECTORY
        str(dereplicated_fasta_path),  # INPUT_FASTA_CLUSTERING  # 5
        "/usr/local/bin/ktImportText",  # KRONA_TOOL
        SOTUs_fasta_name,  # OUTPUT_SOTU_FASTA_WITH_TAXONOMY
        SINA_ARB,  # SILVA_ARB
        BIN_DIR + "sina/sina",  # SINA_EXECUTABLE
        ANALYSIS_DIR,  # USER_FASTQ_FOLDER  # 10
        BIN_DIR + "rapidnj",  # RAPID_NJ_BIN
        SOTUs_table_name,  # SOTUS_TABLE_NAME
        "ZOTU_map.tab",  # INITIAL_ZOTUS_TABLE
        str(abund_limit),  #
        str(ARGS_CLS.create_table.sample_wise_correction)  #
    ]
    try:
        system_sub(main_create_args, capture_output=False, shell=False)
    except Exception as exc:
        raise ValueError(f"Table: {str(exc)}")
    ANA_LOG.info("TIC Done")
    chdir(ANALYSIS_DIR)
    shutil.rmtree(TIC_Result_DIR)
    for file in TIC_Output_DIR.glob("*"):
        shutil.copy(file, ANALYSIS_DIR)
        file.unlink()
    shutil.rmtree(TIC_Output_DIR)
    cleanup(ANALYSIS_DIR)
    ANA_LOG.info('Cleanup DONE')
    # Normalizing Tables
    try:
        zotu_norm_methods = normalize_otu_table(str(Path(PurePath(ANALYSIS_DIR + ZOTUs_table_name))), str(parsed_spike_stat_file_path))
    except Exception as exc:
        ANA_LOG.warning(f"No Normalization applied on {ZOTUs_table_name}: {exc}")
    else:
        ANA_LOG.info(f"{zotu_norm_methods} normalization method(s) applied on {ZOTUs_table_name}")
    try:
        otu_norm_methods = normalize_otu_table(str(Path(PurePath(ANALYSIS_DIR + SOTUs_table_name))), str(parsed_spike_stat_file_path))
    except Exception as exc:
        ANA_LOG.info(f"No Normalization applied on {SOTUs_table_name}: {exc}")
    else:
        ANA_LOG.info(f"{otu_norm_methods} normalization method(s) applied on {SOTUs_table_name}")
    # creating zip file
    shutil.make_archive(
        Path(ANALYSIS_DIR).parent.joinpath("Analysis"),
        'zip',
        root_dir=str(Path(ANALYSIS_DIR).parent),
        base_dir=str(Path(ANALYSIS_DIR).relative_to(Path(ANALYSIS_DIR).parent))
    )
    ANA_LOG.info('Zip file created')
    ANA_LOG.info('ANALYSIS DONE')
    return ANALYSIS_DIR + ZOTUs_table_name, ANALYSIS_DIR + SOTUs_table_name
