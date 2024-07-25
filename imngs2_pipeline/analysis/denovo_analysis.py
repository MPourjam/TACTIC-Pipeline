import os
import re
import shutil
from os import chdir
import os.path as ospath
from multiprocessing import cpu_count
from processing_helper import IMNGS2ArgsParser, gimmelogger, MyCounter
from processing_helper import system_sub as sys_sub
from spike_normalizer import normalize_otu_table
from pathlib import Path, PurePath


BIN_DIR = "/base/binaries/"
# Here we don't add -strand both to usearch as we have already added in preprocessing pipeline
USEARCH_11_BIN = BIN_DIR + "usearch_11_64"
USEARCH_8_bin = BIN_DIR + "usearch8.1"
DB_LOC = "/base/databases/"
SINA_ARB = DB_LOC + "SILVA_138.1_SSURef_NR99_12_06_20_opt.arb"
SINA_BIN = BIN_DIR + "sina/sina"
GOLD_REFDB_USEARCH = DB_LOC + "SILVA-bac-16s-90.udb"
# USERS_DIR = '/srv/crc/users/'
ref16RNAdb_1 = DB_LOC + "silva-bac-16s-id90.fasta"
ref16RNAdb_2 = DB_LOC + "silva-arc-16s-id95.fasta"
SORT_ME_RNA_BIN = BIN_DIR + "sortmerna"
USEARCH_TAIL = '> /dev/null 2>&1'
krona_importtext = BIN_DIR + "Krona/KronaTools/scripts/ImportText.pl"
global SPIKE_STAT_HEADER, TAXED_ZOTU_FILE_NAME, DEREP_READS_FILE_NAME, DEREP_NEW_LABEL
SPIKE_STAT_FILE_COLS = ("#SampleID", "SpikeReads", "spikes_total_weight_in_g", "spike_amount", "parent_path")
SPIKE_STAT_HEADER = "\t".join(list(SPIKE_STAT_FILE_COLS))
TAXED_ZOTU_FILE_NAME = "taxed_ZOTUs.fasta"
DEREP_READS_FILE_NAME = "derep.fasta"
max_pool = int(cpu_count() * 0.7)
DEREP_NEW_LABEL = "Uniq"
global POOL_SIZE
POOL_SIZE = max_pool if max_pool > 0 else 1


def system_sub(*args, **kwargs):
    return sys_sub(*args, logger_obj=ANA_LOG, **kwargs)


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


def append_reads(seq_fasta, sample_id, analysis_dir):
    chdir(analysis_dir)
    seq_fasta_path = os.path.abspath(seq_fasta)
    dataset_name = f"{sample_id}"
    cmd_to_call_list = [
        USEARCH_11_BIN,
        '-fastx_relabel',
        seq_fasta_path,
        '-prefix',
        f"sample={dataset_name};",
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


def dereplication(with_taxonomy=True):
    # cmd_part_0 = USEARCH_11_BIN + " -fastx_uniques analysis.fasta -relabel Zotu"
    # cmd_part_0 += " -fastaout derep.fasta -tabbedout relabel.tab --threads " + str(POOL_SIZE)
    # cmd_part_0 += " > /dev/null 2>/dev/null"
    # system_sub(cmd_part_0)
    sample_name_regex = re.compile(r"sample=(.*);([0-9]+);size=([0-9]+);")

    cmd_to_call_list = [
        USEARCH_11_BIN,
        "-fastx_uniques",
        "analysis.fasta",
        "-relabel",
        DEREP_NEW_LABEL,  # dereplication relaling
        "-fastaout",
        "derep.fasta",
        "-tabbedout",
        "relabel.tab",
        "--threads",
        str(POOL_SIZE),
        '--sizein',
        '--sizeout'
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
        assigned_uniq = str(tokens[1]).strip()
        if assigned_uniq not in zotus_dict.keys():
            zotus_dict[assigned_uniq] = [sample_name_regex.search(str(tokens[0]).strip()).group(0)]
        else:
            zotus_dict[assigned_uniq] += [sample_name_regex.search(str(tokens[0]).strip()).group(0)]
        # creating taxonomy without last semicolon and without 'unclutured' or 'unknown'
        # removing last empty array element gets created by last ;
        # . tokens[-1][:-1] of relabel.tab which is the tax for that cluster
        if with_taxonomy:
            tax_arr = str(tokens[-1][:-1].split('tax=')[1]).split(';')[:-1]  # [:-1] to exclude \n at the end
            to_exclude = ['uncultured', 'unknown', 'incertae']
            # NOTE New TIC will also check for removing this unknown levels
            for tax_i in range(len(tax_arr) - 1, -1, -1):
                tax_lev_val = tax_arr[tax_i]
                ex_bool = [True for te in to_exclude if bool(re.search(te, str(tax_lev_val).strip().lower()))]
                if any(ex_bool):
                    tax_arr = tax_arr[:tax_i]
            corrected_tax = ";".join(tax_arr)
            zotus_tax_dict[assigned_uniq] = corrected_tax
        line = contents.readline()
    # zotus_dict is line --> {"Zotu1" : ["sample_name1;size=XX;", "sample_name2", ...], ... }

    contents.close()
    samples_list = list()

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
    taxonomy_col = '' if not with_taxonomy else '\t' + 'taxonomy'
    out_line = '#Zotu' + curr_line + taxonomy_col + '\n'
    out_file.write(out_line)  # This is the header

    for key, value in zotu_hit_dict.items():
        curr_line = ''
        for sample_key in samples_list:
            # print(value[str(sample_key)])
            curr_line += '\t' + str(value[str(sample_key)])
        if with_taxonomy:
            out_line = key + curr_line + '\t' + zotus_tax_dict[key] + '\n'
        else:
            out_line = key + curr_line + '\n'
        out_file.write(out_line)
    out_file.close()
    # print(list(zotu_hit_dict.keys()))
    # Writing taxonomy to derep.fasta
    if with_taxonomy:
        with open("derep.fasta", "r") as derep_fio, open("derep_with_tax.fasta", "w+") as derep_tax_fio:
            derep_line = derep_fio.readline()
            while derep_line:
                derep_line = derep_line.strip()
                if derep_line.startswith(">"):
                    curr_zotu = str(derep_line[1:]).strip().split(";")[0]
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
    # system_sub(cmd_to_call_list)  # TODO uncomment this line


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
        # "--sizein",
        # "--sizeout",
    ]
    system_sub(cmd_to_call_list, force_log=True)
    onelinefasta("zotus.fasta")
    cmd_to_call_list = [
        "rm",
        # "sorted.fasta",
        "derep.fasta",
    ]
    # system_sub(cmd_to_call_list)


# Filter non 16S sequences
def filter16S():
    # cmd_0 = SORT_ME_RNA_BIN + " --ref " + ref16RNAdb_1 + " --ref " + ref16RNAdb_2 + " --reads "
    # cmd_1 = " zotus.fasta --fastx --aligned good-ZOTUs --other non16SrRNA --workdir . -e 0.1 --num_alignments 1"  # --aligned
    # cmd_2 = " > /dev/null 2>&1"
    # run a RNA filtering step
    # (The program currently do not distinquish between 16S and 18S)
    # system_sub(cmd_0 + cmd_1 + cmd_2)
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
    system_sub(["mv", "out/aligned.fasta", "good-ZOTUs.fasta"])
    system_sub(["rm", "-r", "idx", "kvdb"])


def prepare_zotus():
    contents = read_file('zotus.fasta')
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
        "good-ZOTUs-sized.fasta",  # Although Edgar suggests that unoise3 and cluster_otus should have same input file, but we use
        "-fulldp",
        "-otus",
        "otus1.fa",
        "-minsize",
        "1",  # minsize 1 to keep all OTUs including singletons
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
    # system_sub(["rm", "otus1.fa"]) # TODO uncomment this line


def assign_zotus_to_otus():
    zotus_map = read_file('z2o.tab')
    out_file = open('ZOTUs-OTUs-map.tab', 'w+')
    out_file_2 = open('matched_ZOTUS.txt', 'w+')
    match_regex = re.compile(r";top=(Zotu[0-9]+);")
    zotu_regex = re.compile(r"(Zotu[0-9]+);")
    valid_zotus_dict = dict()
    for line in zotus_map:
        line_tokens = line.split('\t')
        zotu_number = zotu_regex.search(str(line_tokens[0]).strip()).group(1)
        match_type = str(line_tokens[1]).strip()
        if 'otu' in match_type:
            valid_zotus_dict[zotu_number] = str(line_tokens[1]).strip()
        elif 'match' == match_type:
            matched_zotu = match_regex.search(str(line_tokens[2]).strip()).group(1)
            if matched_zotu in valid_zotus_dict.keys():
                # if the zotu matches to a valid OTU
                matched_otu_of_match = valid_zotus_dict[matched_zotu]
                valid_zotus_dict[zotu_number] = matched_otu_of_match

    for key, value in valid_zotus_dict.items():
        line = key + '\t' + value + '\n'
        out_file.write(line)
    out_file.close()
    for key in valid_zotus_dict.keys():
        out_file_2.write(key + '\n')
    out_file_2.close()
    # system_sub(["rm", "z2o.tab"])  # TODO uncomment this line


def keep_good_ZOTUs():
    # cmd_0 = USEARCH_11_BIN + ' -fastx_getseqs good-ZOTUs.fasta -labels '  # aligned_16S_ZOTUS.fa -> good-ZOTUs.fasta
    # cmd_1 = 'matched_ZOTUS.txt -fastaout nochi_ZOTUs.fasta '
    # system_sub(cmd_0 + cmd_1 + USEARCH_TAIL)
    cmd_to_call_list = [
        USEARCH_11_BIN,
        "-fastx_getseqs",
        "zotus.fasta",  # we changed "good-ZOTUs.fasta" to "zotus.fasta" as we already have filtered sequences for human origin reads
        "-labels",
        "matched_ZOTUS.txt",
        "-fastaout",
        "ZOTUs-Seqs.fasta",  # This file contains a list of non-chimeric ZOTUs sequences. Equivalent to nochi_ZOTUs.fasta
    ]
    system_sub(cmd_to_call_list, force_log=True)
    onelinefasta("ZOTUs-Seqs.fasta")
    return "ZOTUs-Seqs.fasta"


def build_ZOTU_table():
    # cmd_0 = USEARCH_11_BIN + ' -otutab analysis.fasta -top_hit_only -zotus nochi_ZOTUs.fasta '
    # cmd_1 = ' -otutabout zotu_table.txt -id 0.97'
    # system_sub(cmd_0 + cmd_1 + USEARCH_TAIL)
    cmd_to_call_list = [
        USEARCH_11_BIN,
        "-otutab",
        "analysis.fasta",
        "-top_hit_only",
        "-zotus",  # if the database sequences are already denoised which is the case here
        "ZOTUs-Seqs.fasta",  # This file contains ZOTUs sequences (both centroids of OTUs and matched ZOTUs to that OTU and non-chimeric)
        "-otutabout",
        "zotu_table.txt",
        "-id",
        "0.97",
        "-threads",
        f"{str(POOL_SIZE)}",
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
    # addTax_new("ZOTUs-Seqs.fasta")


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
    shutil.move(new_taxed_path, filename)


def create_final_ZOTU_table():
    tax_regex = re.compile(r";tax=(.*)$")
    zotu_regex = re.compile(r"(Zotu[0-9]+);")
    zotu_tax_dict = dict()
    with open('ZOTUs-Seqs.fasta') as zotu_seq_fio:
        zotu_seq_line = zotu_seq_fio.readline().strip()
        while zotu_seq_line:
            if zotu_seq_line[0] == '>':
                zotu_name = zotu_regex.search(zotu_seq_line).group(1)
                zotu_tax = tax_regex.search(zotu_seq_line).group(1)
                zotu_tax_dict[zotu_name] = zotu_tax
            zotu_seq_line = zotu_seq_fio.readline().strip()
    out_file = open('ZOTUs-Table.tab', 'w+')
    zotu_table_lines_list = read_file('zotu_table.txt')
    header = '\t'.join(zotu_table_lines_list[0].split('\t')[1:])
    # Writing the header
    out_file.write('#Zotu\t' + header + '\tTaxonomy\n')
    for zotu_table_line in zotu_table_lines_list[1:]:
        curr_zotu = zotu_table_line.split('\t')[0]
        curr_line = zotu_table_line + '\t' + zotu_tax_dict[curr_zotu] + '\n'
        out_file.write(curr_line)
    out_file.close()

    return 'ZOTUs-Table.tab'


def add_taxonomy_to_fasta():
    out_file = open('taxed_ZOTUs.fasta', 'w+')
    table_file = read_file('ZOTUs-table.tab')[1:]
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
                zotus_sina_algn_line = zotus_sina_algn.readline().strip()
                common_aln_seq = zotus_sina_algn_line[lowest_start:highest_end]
                sina_algn_common_part.write(common_aln_seq + '\n')
            zotus_sina_algn_line = zotus_sina_algn.readline().strip()

    return shortened_file_path


def write_otus_sina_algn_file():
    pass


def create_tree(algn_file, method: list = ["rapidnj", "FastTree"], output_file_name_root: str = None):
    algn_file_path = os.path.abspath(algn_file)
    algn_dirname, algn_filename = os.path.split(algn_file_path)
    os.chdir(algn_dirname)
    out_file_name_root = str(algn_filename.split(".")[0])
    output_file_name = output_file_name_root if output_file_name_root else out_file_name_root
    # create both trees
    if "FastTree" in method:
        # cmd_part_0 = 'FastTree -quiet -nosupport -gtr -nt filtered_aligned.fasta'
        # cmd_part_1 = ' > ' + sample + '_FastTree.tre'
        # system_sub(cmd_part_0 + cmd_part_1)
        cmd_to_call_list = [
            "FastTree",
            "-quiet",
            "-nosupport",
            "-gtr",
            "-nt",
            algn_filename,
            ">",
            f"{output_file_name}-Tree-FastTree.tre",
        ]
        system_sub(cmd_to_call_list)

    if "rapidnj" in method:
        # cmd_2_part_0 = BIN_DIR + 'rapidnj filtered_aligned.fasta --input-format fa --cores 6 --alignment-type d '
        # cmd_2_part_1 = '--output-format t --no-negative-length -x ' + sample + '_rapidNJ_tree.tre'
        # system_sub(cmd_2_part_0 + cmd_2_part_1)
        cmd_to_call_list = [
            BIN_DIR + "rapidnj",
            algn_filename,
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


def create_krona(krona_tool):
    OTUfinal_line = read_file('ZOTUs-Table.tab')[1:]
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
    system_sub(
        [
            krona_tool,
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
    system_sub(
        [
            "mv",
            "text.krona.html",
            "krona.html"
        ],
        capture_output=True,
        quiet=True
    )


def create_OTUs_from_ZOTUS_table():
    clean_zotus = read_file('ZOTUs-Table.tab')
    map_file = read_file('ZOTUs-OTUs-map.tab')
    otus_contents = read_file('mOTUs-Seqs.fasta')
    otu_seqs_fio = open('OTUs-Seqs.fasta', 'w+')
    otu_table_fio = open('OTUs-Table.tab', 'w+')
    header = '\t'.join(clean_zotus[0].split('\t')[1:-1])  # ZOTUs-Table.tab already have taxonomy
    otu_table_fio.write('#OTU\t' + header + '\tTaxonomy\n')
    zotu_size_dict = dict()
    zotu_taxo_dict = dict()
    zotu_map_otus = dict()
    zotu_regex = re.compile(r"^>?(Zotu[0-9]+);?")
    for line in clean_zotus[1:]:
        zotu_name = zotu_regex.search(line).group(1)
        taxonomy = line.split('\t')[-1]
        curr_size = '\t'.join(line.split('\t')[1:-1])
        zotu_size_dict[zotu_name] = curr_size
        zotu_taxo_dict[zotu_name] = taxonomy
    for line in map_file:
        zotu_name, otu_name = line.split('\t')
        zotu_map_otus[zotu_name] = otu_name
    # zotus_from_map_file = list(zotu_map_otus.keys())
    # for key in zotus_from_map_file:
    #     if key not in zotu_taxo_dict.keys():
    #         del zotu_map_otus[key]
    # Writing the new table based on OTUs (summing size of zouts belonging to the same OTU)
    otu_summed_size_dict = dict()
    for line in clean_zotus[1:]:
        zotu_ = zotu_regex.search(line).group(1)
        otu_ = zotu_map_otus[zotu_]
        sizes_list = [int(s) for s in zotu_size_dict[zotu_].split("\t")]
        curr_size = otu_summed_size_dict.get(otu_, [0 for i in sizes_list])
        otu_summed_size_dict[otu_] = [a + b for a, b in zip(curr_size, sizes_list)]
    for line in otus_contents:
        if line[0] == '>':
            zotu_name = zotu_regex.search(line).group(1)
            # The if below is useful for cases in which some filteration has happened on the ZOTUs after clusterZOTUs() step
            if zotu_name in zotu_map_otus.keys():
                # It drops all zotus belonging to an OTU and only takes the centriod
                otu_seqs_fio.write('>' + zotu_map_otus[zotu_name] + ';tax=' + zotu_taxo_dict[zotu_name] + '\n')
                out_line_2 = zotu_map_otus[zotu_name] + '\t'
                out_line_2 += "\t".join([str(count) for count in otu_summed_size_dict[zotu_map_otus[zotu_name]]])
                out_line_2 += "\t" + zotu_taxo_dict[zotu_name] + '\n'
                otu_table_fio.write(out_line_2)
                seq_flag = 1
            else:
                seq_flag = 0
        else:
            if seq_flag and re.search(r"^[ACGTU]", line):
                otu_seqs_fio.write(line + '\n')
    otu_seqs_fio.close()
    otu_table_fio.close()

    return 'OTUs-Table.tab', 'OTUs-Seqs.fasta'


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


def create_trees(zotu_algn_file):
    zotu_algn_file = os.path.abspath(zotu_algn_file)
    zotu_algn_dir, zotu_algn_filename = os.path.split(zotu_algn_file)
    os.chdir(zotu_algn_dir)
    # extracint otu entries in zotu_algn_file
    otu_zotu_map_lines = read_file('z2o.tab')
    otu_zotu_map_dict = dict()
    zotu_regex = re.compile(r"^>?(Zotu[0-9]+);?")
    otu_regex = re.compile(r"\s(otu[0-9]+)\s")
    for line in otu_zotu_map_lines:
        zotu, otu_mo = zotu_regex.search(line).group(1), otu_regex.search(line)
        if otu_mo:
            # only centroid zotus
            otu_zotu_map_dict[zotu] = otu_mo.group(1)
    # extracting otu sequences from zotu_algn_file
    with open(zotu_algn_file) as zotu_algn_fio, open('shortened_algn_test_OTUs-Seqs.fasta', 'w+') as otu_algn_fio:
        zotu_algn_line = zotu_algn_fio.readline().strip()
        while zotu_algn_line:
            if zotu_algn_line.startswith('>'):
                zotu_name = zotu_regex.search(zotu_algn_line).group(1)
                otu_ = otu_zotu_map_dict.get(zotu_name, False)
                if otu_:  # only otus assigned to centroid zotus
                    otu_algn_fio.write(">" + str(otu_) + '\n')
                    zotu_algn_line = zotu_algn_fio.readline().strip()  # alignment seq
                    otu_algn_fio.write(zotu_algn_line + '\n')
            zotu_algn_line = zotu_algn_fio.readline().strip()

    # creating trees
    create_tree(zotu_algn_filename, ["rapidnj"], "ZOTUs")
    create_tree("shortened_algn_test_OTUs-Seqs.fasta", ["rapidnj"], "OTUs")


def filter_abundance(abundance_limit,
                     OTU_table_file_name,
                     OTU_fasta_file_name,
                     ZOTU_table_file_name,
                     ZOTU_fasta_file_name,
                     ZOTU_OTU_map_file,
                     sample_wise_corr_flag=True):
    # Managing paths of input files
    # OTU TABLE
    OTU_table_file_path = ospath.abspath(OTU_table_file_name)
    OTU_table_file_parent, OTU_table_file_name = ospath.split(OTU_table_file_path)
    # OTU FASTA
    OTU_fasta_file_path = ospath.abspath(OTU_fasta_file_name)
    OTU_fasta_file_parent, OTU_fasta_file_name = ospath.split(OTU_fasta_file_path)
    # ZOTU TABLE
    ZOTU_table_file_path = ospath.abspath(ZOTU_table_file_name)
    ZOTU_table_file_parent, ZOTU_table_file_name = ospath.split(ZOTU_table_file_path)
    # ZOTU FASTA
    ZOTU_fasta_file_path = ospath.abspath(ZOTU_fasta_file_name)
    ZOTU_fasta_file_parent, ZOTU_fasta_file_name = ospath.split(ZOTU_fasta_file_path)
    # ZOTU OTU MAP
    ZOTU_OTU_map_file_path = ospath.abspath(ZOTU_OTU_map_file)
    ZOTU_OTU_map_file_parent, ZOTU_OTU_map_file_name = ospath.split(ZOTU_OTU_map_file_path)
    zotu_otu_map_dict = dict()
    with open(ZOTU_OTU_map_file_path, "r") as zotu_otu_map_fio:
        zotu_otu_map_fio.readline()
        for line in zotu_otu_map_fio:
            line = line.strip()
            zotu, otu = line.split('\t')
            zotu_otu_map_dict[zotu] = otu

    unf_tot_sizes = get_samples_sizes(OTU_table_file_path)

    zotu_regex = re.compile(r"^>?(Zotu[0-9]+)[;\s]?")

    abundance_limit = round(float(abundance_limit), 5)
    # TIC_LOG.info("SOTUs got filtered for abundance value of {}".format(abundance_limit))
    sample_sizes = list()
    filtered_OTU_table_path = ospath.join(OTU_table_file_parent, 'filtered_' + OTU_table_file_name)
    removed_sotus = []
    with open(OTU_table_file_path, "r") as otu_table_fio, open(filtered_OTU_table_path, "w+") as filtered_otu_fio:
        otu_header = otu_table_fio.readline().strip()
        datasets_samples = len(otu_header.split('\t')) - 2
        for i in range(0, datasets_samples):
            sample_sizes.append(0)
        filtered_otu_fio.write(otu_header + '\n')
        line = otu_table_fio.readline()
        while line:
            line = line.strip()
            line_list = line.split('\t')
            curr_sizes = line_list[1:-1]
            sotu_name = str(line_list[0]).strip()
            rem_sotu = True
            curr_abundances = list()
            for j in range(len(curr_sizes)):
                abundance = round(float(curr_sizes[j]) / unf_tot_sizes[j], 5)
                curr_abundances.append(abundance)
            if sample_wise_corr_flag:
                new_curr_sizes = [str(csz) for csz in curr_sizes]
                for ii in range(len(curr_abundances)):
                    if curr_abundances[ii] <= abundance_limit:
                        new_curr_sizes[ii] = str(round(0.0, 1))
                    else:
                        rem_sotu = False
                if not rem_sotu:
                    new_line = "\t".join([sotu_name, "\t".join(new_curr_sizes), str(line_list[-1]).strip()])

            elif any([x >= abundance_limit for x in curr_abundances]):
                rem_sotu = False
                new_line = line
            else:
                new_line = line

            if not rem_sotu:
                filtered_otu_fio.write(new_line + '\n')
            else:
                removed_sotus.append(sotu_name)
            # read new line
            line = otu_table_fio.readline()

    shutil.move(filtered_OTU_table_path, OTU_table_file_path)

    # Removing the filtered SOTUs from the fasta file
    filtered_otu_fasta_path = ospath.join(OTU_fasta_file_parent, 'filtered_' + OTU_fasta_file_name)
    with open(OTU_fasta_file_path, 'r') as otu_fasta_fio, open(filtered_otu_fasta_path, "w+") as filtered_otu_fasta_fio:
        derep_line = otu_fasta_fio.readline()
        write_seq = True
        while derep_line:
            if derep_line[0] == '>':
                curr_sotu = str(derep_line[1:-1]).split(";")[0]  # removing "\n"
                if curr_sotu not in removed_sotus:
                    new_line = derep_line
                    filtered_otu_fasta_fio.write(new_line)
                    write_seq = True
                else:
                    write_seq = False
            elif write_seq:
                filtered_otu_fasta_fio.write(derep_line)
            derep_line = otu_fasta_fio.readline()

    shutil.move(filtered_otu_fasta_path, OTU_fasta_file_path)

    # Removing related ZOTUs to deleted SOTUs in ZOTU Table
    filtered_zotu_table_path = ospath.join(ZOTU_table_file_parent, "filtered_" + ZOTU_table_file_name)
    with open(ZOTU_table_file_path, "r") as zotu_table_fio, open(filtered_zotu_table_path, "w+") as filtered_zotu_fio:
        # Writing header
        zotu_header = zotu_table_fio.readline()
        filtered_zotu_fio.write(zotu_header)
        # other lines
        next_line = zotu_table_fio.readline()
        while next_line:
            sotu_of_zotu = zotu_otu_map_dict.get(zotu_regex.search(next_line).group(1), None)
            if sotu_of_zotu not in removed_sotus:
                filtered_zotu_fio.write(next_line)
            next_line = zotu_table_fio.readline()

    shutil.move(filtered_zotu_table_path, ZOTU_table_file_path)

    # Removing related ZOTUs to deleted SOTUs in ZOTU Table
    filtered_zotu_fasta_path = ospath.join(ZOTU_fasta_file_parent, "filtered_" + ZOTU_fasta_file_name)
    with open(ZOTU_fasta_file_path, 'r') as zotu_fasta_fio, open(filtered_zotu_fasta_path, "w+") as filtered_zotu_fasta_fio:
        zotu_fasta_seq_line = zotu_fasta_fio.readline()
        write_seq = False
        while zotu_fasta_seq_line:
            if str(zotu_fasta_seq_line).startswith(">"):
                sotu_of_zotu = zotu_otu_map_dict.get(zotu_regex.search(zotu_fasta_seq_line).group(1), None)
                if sotu_of_zotu not in removed_sotus:
                    write_seq = True
                    filtered_zotu_fasta_fio.write(zotu_fasta_seq_line)
                else:
                    write_seq = False
            elif write_seq:
                filtered_zotu_fasta_fio.write(zotu_fasta_seq_line)
            zotu_fasta_seq_line = zotu_fasta_fio.readline()

    shutil.move(filtered_zotu_fasta_path, ZOTU_fasta_file_path)


def cleanup(directory: str):
    directory = Path(PurePath(directory)).absolute()
    to_keep = [
        "spike_mapping_file.csv",
        "ZOTUs-Seqs.fasta",
        "ZOTUs-Table.tab",
        "OTUs-Seqs.fasta",
        "OTUs-Table.tab",
        "ZOTUs-OTUs-map.tab",
        "ZOTUs-Tree-nj.tre",
        "ZOTUs-Tree-FastTree.tre",
        "OTUs-Tree-nj.tre",
        "OTUs-Tree-FastTree.tre",
        "krona.html",
    ]
    files = list(directory.glob("*"))
    for file in files:
        if file.name not in to_keep and file.is_file():
            file.unlink()
    return True


def create_zip(analysis_id):
    analysis_dir = ANALYSIS_DIR.parent.joinpath(analysis_id)
    chdir(analysis_dir)
    shutil.make_archive(analysis_id, 'zip', analysis_dir)
    # cmd_zip = 'zip -r ' + query_id + '.zip ' + query_id
    # system(cmd_zip)


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


def main_de_novo(
        analysis_dir: str,
        spike_stat_file: str,
        fastqs_dir: str,
        args_file_path: str = "",
        dbs_loc: str = DB_LOC,
        threads=POOL_SIZE,
        usearch_11_bin: str = USEARCH_11_BIN):

    """
    sample_seq_files_path: a list of file path to
     each samples' sequence file which is going to be combined with other samples passed to analysis.
    analysis_dir: The destination directory to save results
    """
    global USEARCH_11_BIN, POOL_SIZE, ANALYSIS_DIR, ARGS_CLS, ANA_LOG, DB_LOC
    # updating USEARCH_11_BIN
    USEARCH_11_BIN = usearch_11_bin

    spike_stat_file = Path(PurePath(spike_stat_file))
    assert spike_stat_file.is_file(), "spike_stat_file must be a path to a file"

    # updating POOL_SIZE
    try:
        POOL_SIZE = int(threads)
    except ValueError:
        ANA_LOG.warning(f"threads must be an integer. Using default value of {POOL_SIZE}")
    DB_LOC = Path(PurePath(dbs_loc))
    assert DB_LOC.is_dir(), "DBS_LOC should be path to direcotry containing SILVA database files (arb)"
    ANALYSIS_DIR = Path(PurePath(analysis_dir)).absolute()
    gimmelogger(
        logger_name="run_imngs2.analysis",
        log_file=ANALYSIS_DIR.joinpath("Analysis_log.txt"),
        only_file=True,
    )
    ANA_LOG = gimmelogger(
        logger_name="run_imngs2.analysis.denovo_analysis",
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
        # sample_seq_files_path.append((Path(PurePath(entries[4])).joinpath(DEREP_READS_FILE_NAME), entries[0]))
        sample_seq_files_path.append((Path(PurePath(entries[4])).joinpath(DEREP_READS_FILE_NAME), entries[0]))

    for derep_path, sam_id in sample_seq_files_path:
        if not derep_path.is_file():
            raise ValueError(f"{str(derep_path)} must contain path to each samples sequence file!")

    chdir(ANALYSIS_DIR)
    ANA_LOG.info('### Analysis Started ###')
    ANA_LOG.info('# Gathering Sequences: Started')
    # TODO as derep.fasta is archived in the zip file, it should be extracted first
    for seq_file, sam_id in sample_seq_files_path:
        append_reads(seq_file, sam_id, ANALYSIS_DIR)
    chdir(ANALYSIS_DIR)
    # trimming the sequences
    trim_sides(ARGS_CLS.trimsides.stripleft, ARGS_CLS.trimsides.stripright)
    ANA_LOG.info('# Dereplication: Started')
    dereplication(with_taxonomy=False)
    # Sorting the sequences
    ANA_LOG.info('# Sorting: Started')
    sort_seqs()
    # Clustering the sequences
    ANA_LOG.info('# Clustering at ZOTU level: Started')
    clusterZOTUs()
    # Removing the size from the ZOTUs
    prepare_zotus()
    # Clustering at OTU level
    ANA_LOG.info('# Clustering at OTU level: Started')
    clusterOTUs()
    ANA_LOG.info('# Creating ZOTU tables: Started')
    remove_size_from_OTUS()
    assign_zotus_to_otus()
    ZOTUs_seq_name = keep_good_ZOTUs()
    build_ZOTU_table()
    ANA_LOG.info('# Adding taxonomy to ZOTU sequences: Started')
    addTax_new(ZOTUs_seq_name)  # taxed_ZOTUs-Seqs.fasta gets created here and then replaces ZOTUs-Seqs.fasta
    ZOTUs_table_name = create_final_ZOTU_table()
    ANA_LOG.info('# Creating OTU table: Started')
    OTUs_table_name, OTU_seq_name = create_OTUs_from_ZOTUS_table()
    # Creating trees
    ANA_LOG.info('# Creating trees: Started')
    sina_algn_shortened_file = shorten_sina_algn("test_ZOTUs-Seqs.fasta")
    ANA_LOG.info('# Creating Krona: Started')
    create_krona(krona_importtext)
    ANA_LOG.info('# Filtering OTUs abundance: Started')
    filter_abundance(
        ARGS_CLS.create_table.abund_limit,
        OTUs_table_name,
        OTU_seq_name,
        ZOTUs_table_name,
        ZOTUs_seq_name,
        "ZOTUs-OTUs-map.tab")
    ANA_LOG.info('# Cleaning up: Started')
    try:
        create_trees(sina_algn_shortened_file)  # TODO check if we need to use create_trees
    except Exception as exc:
        ANA_LOG.warning(f"Tree creation skipped: {exc}")
    cleanup(str(ANALYSIS_DIR))
    # Normalizing Tables
    try:
        zotu_norm_methods = normalize_otu_table(str(Path(PurePath(ANALYSIS_DIR + ZOTUs_table_name))), str(parsed_spike_stat_file_path))
    except Exception as exc:
        ANA_LOG.warning(f"No Normalization applied on {ZOTUs_table_name}: {exc}")
    else:
        ANA_LOG.info(f"{zotu_norm_methods} normalization method(s) applied on {ZOTUs_table_name}")
    try:
        otu_norm_methods = normalize_otu_table(str(Path(PurePath(ANALYSIS_DIR + OTUs_table_name))), str(parsed_spike_stat_file_path))
    except Exception as exc:
        ANA_LOG.info(f"No Normalization applied on {OTUs_table_name}: {exc}")
    else:
        ANA_LOG.info(f"{otu_norm_methods} normalization method(s) applied on {OTUs_table_name}")

    ANA_LOG.info('ANALYSIS DONE')
    # to keep the same return value as the main function
    return ANALYSIS_DIR + ZOTUs_table_name, ANALYSIS_DIR + OTUs_table_name
