from os import chdir, system, mkdir, listdir, makedirs
from statistics import stdev, mean
from re import search
import os
import re
import random
import shutil
import glob
import logging
import sys


BIN_DIR = "/cfm/binaries/"
USEARCH_8_BIN = BIN_DIR + "usearch8.1"
USEARCH_11_bin = BIN_DIR + "usearch_11_64"
SORT_ME_RNA_BIN = '/cfm/binaries/sortmerna'
# USEARCH8_1 = 'usearch8.1 -threads 4'
USEARCH8_1 = USEARCH_8_BIN + " -threads 4"
# GOLD_REFDB_USEARCH = '/srv/cfm/databases/silva_700K.udb'
GOLD_REFDB_USEARCH = '/srv/cfm/databases/SILVA-bac-16s-90.udb'
SINA_ARB = '/srv/cfm/databases/SILVA_138_SSURef_NR99_05_01_20_opt.arb'
ref16RNAdb_1 = "/srv/cfm/databases/silva-bac-16s-id90.fasta"
ref16RNAdb_2 = "/srv/cfm/databases/silva-arc-16s-id95.fasta"
USEARCH_TAIL = '> /dev/null 2>&1'
bowtie2 = "/cfm/binaries/bowtie2/bowtie2"
spikeidx = "/cfm/binaries/bowtie2/spikesidx/spike"
krona_importtext = "/cfm/binaries/Krona/KronaTools/scripts/ImportText.pl"


def calc_spikes(*fastq_files, spike_amount):
    fastq_names = [f for f in fastq_files if bool(f)]
    fastqs_abs_paths = [os.path.abspath(os.path.join(os.getcwd(), f)) for f in fastq_names]

    if str(spike_amount) == "0":
        line_count = 0
        for f in fastqs_abs_paths:
            with open(f, 'r') as fqfile:
                for _ in fqfile:
                    line_count += 1
        line_count = line_count // 4  # (total number of reads, spike reads)
        return line_count, 0  # non_spike_reads spike_reads
    else:
        spike_amount = float("{}e-9".format(spike_amount))  # ng to g
        name_regx = r"[a-zA-Z0-9\-]+"
        spike_res_dir = os.path.split(fastqs_abs_paths[0])[0] + "/spike_result/"
        makedirs(spike_res_dir, mode=777, exist_ok=True)
        _, forw_name = os.path.split(fastqs_abs_paths[0])
        fastq_aligned = forw_name[search(name_regx, forw_name).start():search(name_regx, forw_name).end()]
        fastq_aligned = os.path.join(spike_res_dir, fastq_aligned)
        fastq_unaligned = fastq_aligned + "_unal"

        if len(fastq_names) == 2:
            cmd = [bowtie2, "-x", spikeidx, "-1", fastqs_abs_paths[0], "-2",
                   fastqs_abs_paths[1], "--al-conc", fastq_aligned, "--un-conc",
                   fastq_unaligned, USEARCH_TAIL]
            os.system(" ".join(cmd))
            spike_counter = 0
            for fi in [fastq_aligned + ".1", fastq_aligned + ".2"]:
                with open(fi, 'r') as fastq_al:
                    for _ in fastq_al:
                        spike_counter += 1
            spike_counter = spike_counter // 4
            if spike_counter != 0:
                shutil.copy(fastq_unaligned + ".1", fastqs_abs_paths[0])
                shutil.copy(fastq_unaligned + ".2", fastqs_abs_paths[1])

        elif len(fastq_names) == 1:
            cmd = [bowtie2, "-x", spikeidx, "-U", fastqs_abs_paths[0], "--al-conc",
                   fastq_aligned, "--un-conc", fastq_unaligned, USEARCH_TAIL]
            os.system(" ".join(cmd))
            files_inspike = listdir(spike_res_dir)
            spike_counter = 0
            for fi in [f for f in files_inspike if search(fastq_aligned, os.path.abspath(f))]:
                with open(fi, 'r') as fastq_al:
                    for _ in fastq_al:
                        spike_counter += 1
            spike_counter = spike_counter // 4
            # Check if the unaligned file is empty
            if spike_counter != 0:
                shutil.copy(fastq_unaligned + ".1", fastqs_abs_paths[0])

        shutil.rmtree(spike_res_dir, ignore_errors=True)
        unaligned_reads, _ = calc_spikes(*fastqs_abs_paths, spike_amount=0)
        # to pass it as reads_number and spike number to seq_met model
        return unaligned_reads, spike_counter


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
    return(len_mean, len_stdev)


def read_file(filename):
    with open(filename) as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    return content


def unzip():
    system('unzip upload.zip')


def only_keep_dataset_fastqs(*args):
    mkdir('uploadedfastqs')
    system('mv upload.zip ./uploadedfastqs/')
    chdir("./uploadedfastqs/")
    unzip()
    fastqs = glob.glob('*')
    lose = []
    for i, fq in enumerate(fastqs):
        if fq not in list(args):
            lose.append(fq)
    system('rm -rf {}'.format(" ".join(lose)))
    system('zip upload.zip *')
    system('mv upload.zip ../upload.zip')
    chdir("../")
    system('rm -rf uploadedfastqs')
    unzip()


def clean_FastQC(input_file, reverse_file=False):
    part_F = ".".join(input_file.split('/')[-1].split('.')[:-1])
    output = part_F + '_fastqc'
    cmd_0 = 'unzip ' + output + '.zip > /dev/null 2>&1'
    cmd_1 = 'mv ' + output + '/Images/per_base_quality.png ../R1-per_base_quality.png'
    cmd_2 = 'rm -r ' + output + '.zip ' + output + '.html'
    try:
        system(cmd_0)
        system(cmd_1)
        system(cmd_2)
    except BaseException:
        pass

    if reverse_file:
        part_R = ".".join(reverse_file.split('/')[-1].split('.')[:-1])
        output = part_R + '_fastqc'
        cmd_0 = 'unzip ' + output + '.zip > /dev/null 2>&1'
        cmd_1 = 'mv ' + output + '/Images/per_base_quality.png ../R2-per_base_quality.png'
        cmd_2 = 'rm -r ' + output + '.zip ' + output + '.html'
        try:
            system(cmd_0)
            system(cmd_1)
            system(cmd_2)
        except BaseException:
            pass


def mymkdir(input_dir):
    try:
        mkdir(input_dir)
    except BaseException:
        pass


def merge_pairs(forward_file, reverse_file):
    cmd_0 = USEARCH_11_bin + " -fastq_mergepairs " + forward_file + ' -reverse ' + reverse_file + ' -fastq_maxdiffs 50'
    cmd_1 = " -fastq_pctid 20 -fastqout merged.fasta -fastq_trunctail 20"
    cmd_2 = " -fastq_minmergelen 200 -fastq_maxmergelen 600 >/dev/null 2>/dev/null"
    system(cmd_0 + cmd_1 + cmd_2)


def trim_sides():
    cmd_0 = USEARCH_11_bin + " -fastx_truncate merged.fasta -stripright 10"
    cmd_1 = " -stripleft 10 -fastqout filtered1.fasta >/dev/null 2>/dev/null"
    system(cmd_0 + cmd_1)


def filter_merged_reads():
    cmd_0 = USEARCH_11_bin + ' -fastq_filter filtered1.fasta -fastq_maxee_rate 0.005'
    cmd_1 = " -fastaout filtered2.fasta >/dev/null 2>/dev/null"
    system(cmd_0 + cmd_1)


def run_FastQC(forward_file, reverse_file):
    out_dir = 'fastqc_output'
    mymkdir(out_dir)
    cmd_0 = "fastqc " + forward_file
    cmd_1 = " --outdir " + out_dir + " --format fastq --threads 6 --quiet"
    if not reverse_file == '':
        cmd_0 = cmd_0 + " " + reverse_file
    try:
        system(cmd_0 + cmd_1)
    except BaseException:
        raise("FastQC did NOT run")
    chdir(out_dir)
    if not reverse_file == '':
        clean_FastQC(input_file=forward_file, reverse_file=reverse_file)
    else:
        clean_FastQC(input_file=forward_file)


def dereplicate_seqs():
    cmd_0 = USEARCH_11_bin + ' -fastx_uniques filtered2.fasta -fastaout derep.fasta'
    cmd_1 = ' -sizein -sizeout ' + USEARCH_TAIL
    system(cmd_0 + cmd_1)
    line_n = 0
    with open('derep.fasta') as derep:
        line = derep.readline()
        while line:
            if line.startswith(">"):
                line_n += 1
            line = derep.readline()
    return line_n


def sort_seqs():
    cmd_0 = USEARCH_11_bin + ' -sortbysize derep.fasta -fastaout sorted.fasta'
    cmd_1 = ' ' + USEARCH_TAIL
    system(cmd_0 + cmd_1)


def trim_one_side(forward_file):
    cmd_0 = USEARCH_11_bin + " -fastx_truncate " + forward_file + " -stripleft 10"
    cmd_1 = " -fastqout filtered1.fasta >/dev/null 2>/dev/null"
    system(cmd_0 + cmd_1)


def filter_merged_one_side(forward_file):
    curr_mean, curr_sd = seqFileStats(forward_file)
    minLength = curr_mean - int(0.1 * curr_mean) - 5  # remove the primer triming size plus 10% of the mean size
    cmd_0 = USEARCH_8_BIN + ' -fastq_filter filtered1.fasta -fastq_truncqual 20'
    cmd_1 = ' -fastq_maxee_rate 0.005 -fastq_trunclen ' + str(minLength)
    cmd_2 = ' -fastaout filtered2.fasta >/dev/null 2>/dev/null'
    system(cmd_0 + cmd_1 + cmd_2)


# cluster sequences to OTUs
def clusterZOTUs():
    cmd_part_0 = USEARCH_11_bin + ' -unoise3 sorted.fasta -minsize 4 -zotus zotus.fasta'
    cmd_part_1 = ' -tabbedout denoising.tab'
    cmd_part_2 = ' ' + USEARCH_TAIL
    # clusering of seq in OTUs (clustered)
    system(cmd_part_0 + cmd_part_1 + cmd_part_2)


# Filter non 16S sequences
def filter16S():
    cmd_0 = SORT_ME_RNA_BIN + " --ref " + ref16RNAdb_1 + " --ref " + ref16RNAdb_2 + " --reads "
    cmd_1 = " zotus.fasta --fastx good-ZOTUs --other other.non16rRNA  --workdir . -e 0.1"
    cmd_2 = " > /dev/null 2>&1"
    # run a RNA filtering step
    # (The program currently do not distinquish between 16S and 18S)
    system(cmd_0 + cmd_1 + cmd_2)
    system('mv out/aligned.fasta good_ZOTUs.fa')
    # system('rm -r idx out kvdb')


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
    cmd_0 = USEARCH_11_bin + ' -otutab filtered2.fasta -zotus ZOTUs.fasta '
    cmd_1 = ' -otutabout zotu_table.txt -id 0.97'
    system(cmd_0 + cmd_1 + USEARCH_TAIL)


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
    out_file_2 = open('filtered_zotu_table_list.txt', 'w+')
    for line in contents:
        curr_size = float(line.split('\t')[1])
        if (curr_size / unf_tot_size) >= 0:
            out_file.write(line + '\n')
            out_file_2.write(line.split('\t')[0] + '\n')
    out_file.close()
    out_file_2.close()


def select_zotu_seqs():
    cmd_0 = USEARCH_11_bin + ' -fastx_getseqs good_ZOTUs.fa -labels filtered_zotu_table_list.txt '
    cmd_1 = ' -fastaout ZOTUs-Seqs.fasta '
    system(cmd_0 + cmd_1 + USEARCH_TAIL)


def addTax(input_id):
    classifier_dir = '/cfm/binaries/sina/'
    filebasename = 'aligned_' + str(input_id)
    cmd_part_0 = 'sina --in ZOTUs-Seqs.fasta --search --meta-fmt csv '
    cmd_part_1 = '--threads 4 --lca-fields tax_slv '
    cmd_part_2 = '--db ' + SINA_ARB + ' --out ' + filebasename + '.fasta'
    cmd_part_3 = ' >/dev/null 2>/dev/null'
    system(classifier_dir + cmd_part_0 + cmd_part_1 + cmd_part_2 + cmd_part_3)
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
    system("{} ./krona.txt".format(KRONA_TOOL))
    # system("rm ./krona.txt")


def create_zip(input_id):
    cmd_1 = 'zip ../' + str(input_id) + '_processed.zip ZOTUs-table.final.tab '
    cmd_1 += 'taxed_ZOTUs.fasta *.png *.html'
    system(cmd_1)


def create_udb(input_id):
    cmd_0 = USEARCH_11_bin + ' -makeudb_ublast aligned_' + str(input_id) + '.fasta -output ' + str(input_id) + '.udb'
    system(cmd_0 + USEARCH_TAIL)


def cleanup(input_id):
    deleted_files = 'zotu_table.txt good_ZOTUs.fa test.csv filtered1.fasta filtered2.fasta'
    deleted_files += ' aligned_' + str(input_id) + '.csv merged.fasta'
    deleted_files += ' derep.fasta sorted.fasta zotus.fasta otus1.fa z2o.tab mOTUs-Seqs.fasta ZOTUs-Table.tab'
    deleted_files += ' classifiedF.txt filtered_zotu_table_list.txt denoising.tab'
    deleted_files += ' matched_ZOTUS.txt ZOTUs.fasta ZOTUs-Seqs.fasta zotu_table_filtered.txt nochi-ZOTUs.fasta'
    deleted_files += '*.fastq'
    deleted_dirs = ' -r kvdb out idx'
    system('rm ' + deleted_files + " 2> /dev/null")
    system('rm ' + deleted_dirs + " 2> /dev/null")


def update_s_flat(input_id, origin):
    phrase = str(input_id) + ',' + origin + ','
    table_file = read_file('ZOTUs-table.final.tab')[1:]
    sizes = list()
    for line in table_file:
        size = int(line.split('\t')[1])
        sizes.append(size)
    phrase += str(len(sizes)) + ',' + str(sum(sizes)) + '\n'
    out_file_name = '/crc/crc/crcapp/jobs/s_flat.txt'
    out_file = open(out_file_name, 'a+')
    out_file.write(phrase)
    out_file.close()


def find_silva_start_end(fasta_file, max_seq=None, write=None):
    '''
    fasta_file = the path to the alignment fasta output file
    max_seq = the nth first sequence to calculate start and end mode accordingly
    wirte = path to write the new alignment header with start and end of each sequence
    or not, Default won't write anything and just return the tuple of (start_mode, end_mode)
    '''
    if os.path.isfile(os.path.abspath(fasta_file)):
        filedir, filepath = os.path.split(os.path.abspath(fasta_file))
        os.chdir(filedir)
    else:
        raise FileNotFoundError("Fasta file not found!")
    try:
        max_seq = int(max_seq)
    except Exception:
        max_seq = 0

    # Handling output
    if bool(write):
        try:
            outdir, outpath = os.path.split(os.path.abspath(write))
            if outdir != filedir:
                print(f"Writing output to {os.path.abspath(write)}")
            else:
                if outpath == filepath:
                    msg = '''Writing output to the same file!!
                            THE ORIGINAL FILE WOULD BE OVERWRITTEN'''
                    print(msg)
        except Exception as e:
            write = False
            return e

    with open(filepath) as f:
        content = f.readlines()
        silvout_lines = [x.strip() for x in content]
    zotu_names_idx = []
    for i, l in enumerate(silvout_lines):
        if l[0] == ">":
            zotu_names_idx.append(i)
    seq_reg_start = re.compile(r'^-*[ACGTU]')
    seq_reg_end = re.compile(r'[ACGTU]-*$')
    s_e = []
    # Calclating start and end
    # Trying to do it for only the first 20 reads
    try:
        divi = len(zotu_names_idx) // int(max_seq)
        if not divi:
            max_seq = len(zotu_names_idx)
        else:
            max_seq = int(max_seq) + 1
            msg = "Calculating number start and end position mode "
            msg += f"for thefor the first {str(max_seq - 1)} sequences"
            print(msg)
    except Exception:
        max_seq = len(zotu_names_idx)

    for zi in zotu_names_idx[:max_seq]:
        seque = str(silvout_lines[zi + 1])
        m_start = seq_reg_start.match(seque)
        m_end = seq_reg_end.search(seque)
        start = m_start.end()
        end = m_end.start()
        s_e.append((start - 1, end))
    start_set = {s[0] for s in s_e}
    end_set = {s[1] for s in s_e}
    start_dict = dict(zip(start_set, [0 for i in range(len(start_set))]))
    end_dict = dict(zip(end_set, [0 for i in range(len(end_set))]))
    for o in s_e:
        curr_s = start_dict[o[0]]
        curr_e = end_dict[o[1]]
        start_dict[o[0]] = curr_s + 1
        end_dict[o[1]] = curr_e + 1
    # Calculating the mode of start and end positions
    # If mutiple max and min exists.
    list_start_ordered = dict(sorted(start_dict.items(), key=lambda x: x[0], reverse=True))
    list_end_ordered = dict(sorted(end_dict.items(), key=lambda x: x[0], reverse=False))
    start_mode_value = max(list(list_start_ordered.values()))
    end_mode_value = max(list(list_end_ordered.values()))
    start_mode = int(list(list_start_ordered)[list(list_start_ordered.values()).index(start_mode_value)])
    end_mode = int(list(list_end_ordered)[list(list_end_ordered.values()).index(end_mode_value)])

    # Adding start and end position to the end of file
    if write:
        print(f"Writing output to {os.path.abspath(write)}")
        with open(os.path.abspath(write), "w+") as silvout_st_en:
            for i, idx in enumerate(zotu_names_idx):  # + 1 because user finds the bp 1-based
                if i in range(max_seq):
                    start_end = " " + ":".join([str(se + 1) if isinstance(se, int)
                                                else "NotFound" for se in list(s_e[i])])
                else:
                    start_end = ""
                newheader = silvout_lines[idx] + start_end
                silvout_st_en.write(newheader + "\n")
                if i == len(zotu_names_idx) - 1:
                    silvout_st_en.write(silvout_lines[idx + 1])
                else:
                    silvout_st_en.write(silvout_lines[idx + 1] + "\n")

    # print(str(start_mode) +  " " +str(end_mode))
    return start_mode, end_mode


def gimmelogger(name, dir_path):
    dir_path = os.path.abspath(dir_path)
    log = logging.getLogger(str(name))
    log.setLevel(logging.DEBUG)
    log_file = os.path.join(dir_path, "{}_logs.txt".format(str(name)))
    FH = logging.FileHandler(log_file)
    FH.setLevel(logging.INFO)
    log.addHandler(FH)
    SH = logging.StreamHandler(sys.stdout)
    SH.setLevel(logging.DEBUG)
    log.addHandler(SH)
    return log


def main_processing(input_dir, paired, forward_file, reverse_file, input_id, spike_amount=0):
    log = gimmelogger(input_id, input_dir)
    try:
        chdir(input_dir)
        # only_keep_dataset_fastqs(forward_file, reverse_file)
        log.info("Spike removal started.")
        real_reads_c, spike_reads_c = calc_spikes(*[forward_file, reverse_file], spike_amount=spike_amount)
        log.debug("Actual_reads:{}\tSpike_reads:{}\n".format(real_reads_c, spike_reads_c))
        run_FastQC(forward_file, reverse_file)
        log.info('fastQC DONE')
        chdir(input_dir)
        system("rm -r fastqc_output")
        if reverse_file:
            merge_pairs(forward_file, reverse_file)
            log.info('Merging Pairs DONE')
            trim_sides()
            log.info('Trim Sides DONE')
            filter_merged_reads()
            log.info('Filter merged DONE')
        else:
            trim_one_side(forward_file)
            log.info('Trim One Side DONE')
            filter_merged_one_side(forward_file)
            log.info('Filter one DONE')
        dereped_read_n = dereplicate_seqs()
        initial_report = "{}\n".format(str(input_id))
        initial_report += "Actual_reads:{}\tSpike_reads:{}\tDereplicated_reads:{}\n"
        init_rep = initial_report.format(real_reads_c, spike_reads_c, dereped_read_n)
        log.info(init_rep)
        log.info('Dereplication DONE')
        sort_seqs()
        log.info('Sorting DONE')
        clusterZOTUs()
        log.info('Cluster ZOTUs DONE')
        filter16S()
        log.info('Filtering our non 16S ZOTUs')
        prepare_zotus()
        log.info('Prepared Zotus Done')
        build_ZOTU_table()
        log.info('Build ZOTU table')
        filter_zotu_abundance()  # ignore this step because the required abundance is 0
        log.info('ZOTUs Abundance Filtered')
        select_zotu_seqs()
        log.info('Select ZOTUs Filtered')
        addTax(input_id)
        log.info('Sina taxonomy added.')
        create_final_ZOTU_table()
        add_taxonomy_to_fasta()
        addKrona(krona_importtext)
        log.info("Krona graph added.")
        system('Rscript /crc/crc/crcapp/jobs/processing_stats.R >/dev/null 2>/dev/null')
        log.info('Relabing DONE')
        # udb for both similarity queries
        create_udb(input_id)
        log.info('UDB created')
        # update_s_flat(input_id, origin)
        start_mode, end_mode = find_silva_start_end('aligned_' + str(input_id) + '.fasta')
        with open("report.txt", 'w+') as report:
            to_w = "Start_mode\tEnd_mode\tActual_Raw_reads\tSpike_reads\tDereplicated_reads"
            to_w += "\n{}\t{}\t{}\t{}\t{}\n".format(start_mode, end_mode, real_reads_c, spike_reads_c, dereped_read_n)
            log.info(to_w)
            report.write(to_w)
        cleanup(input_id)
        log.info("Cleaned up: {}".format(input_id))
        create_zip(input_id)
        log.info("Zipped!")
        system("rm {}_logs.txt".format(input_id))
    except BaseException:
        raise ValueError
