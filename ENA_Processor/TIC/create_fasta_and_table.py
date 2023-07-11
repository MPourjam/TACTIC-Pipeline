from sys import argv
from os import mkdir, system, getcwd
from os import path as ospath
from glob import glob


def read_file(filename):
    with open(filename) as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    return content


def fasta2dict(fil):
    dic = {}
    cur_scaf = ''
    cur_seq = []
    for line in open(fil):
        if line.startswith(">") and cur_scaf == '':
            cur_scaf = line[1:].rstrip()
        elif line.startswith(">") and cur_scaf != '':
            dic[cur_scaf] = ''.join(cur_seq)
            cur_scaf = line[1:].rstrip()
            cur_seq = []
        else:
            cur_seq.append(line.rstrip())
    dic[cur_scaf] = ''.join(cur_seq)
    new_dict = dict()
    for key, value in dic.items():
        new_key = key.split(';')[0]
        new_dict[new_key] = value
    return new_dict


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


def filter_abundance(abundance_limit,
                     OTU_table_file_name,
                     OTU_fasta_file_name,
                     ZOTU_table_file_name,
                     ZOTU_fasta_file_name,
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
    unf_tot_sizes = get_samples_sizes(OTU_table_file_path)
    # 
    abundance_limit = round(float(abundance_limit), 5)
    print("Filtered for abundance value of {}".format(abundance_limit))
    sample_sizes = list()
    filtered_OTU_table_path = ospath.join(OTU_table_file_parent , 'filtered_' + OTU_table_file_name)
    with open(OTU_table_file_path, "r") as otu_table_fio, open(filtered_OTU_table_path, "w+") as filtered_otu_fio:
        otu_header = otu_table_fio.readline().strip()
        datasets_samples = len(otu_header.split('\t')) - 2
        for i in range(0, datasets_samples):
            sample_sizes.append(0)
        filtered_otu_fio.write(otu_header + '\n')
        removed_sotus = []
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

    system("mv " + filtered_OTU_table_path + " " + OTU_table_file_path)

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

    system("mv " + filtered_otu_fasta_path + " " + OTU_fasta_file_path)

    # Removing related ZOTUs to deleted SOTUs in ZOTU Table
    filtered_zotu_table_path = ospath.join(ZOTU_table_file_parent, "filtered_" + ZOTU_table_file_name)
    with open(ZOTU_table_file_path, "r") as zotu_table_fio, open(filtered_zotu_table_path, "w+") as filtered_zotu_fio:
        # Writing header
        zotu_header = zotu_table_fio.readline()
        filtered_zotu_fio.write(zotu_header)
        # other lines
        next_line = zotu_table_fio.readline()
        while next_line:
            sotu_of_zotu = str(next_line.split(";")[-2]).strip()
            if sotu_of_zotu not in removed_sotus:
                filtered_zotu_fio.write(next_line)
            next_line = zotu_table_fio.readline()

    system("mv " + filtered_zotu_table_path + " " + ZOTU_table_file_path)

    # Removing related ZOTUs to deleted SOTUs in ZOTU Table
    filtered_zotu_fasta_path = ospath.join(ZOTU_fasta_file_parent, "filtered_" + ZOTU_fasta_file_name)
    with open(ZOTU_fasta_file_path, 'r') as zotu_fasta_fio, open(filtered_zotu_fasta_path, "w+") as filtered_zotu_fasta_fio:
        zotu_fasta_seq_line = zotu_fasta_fio.readline()
        write_seq = False
        while zotu_fasta_seq_line:
            if str(zotu_fasta_seq_line).startswith(">"):
                sotu_of_zotu = str(zotu_fasta_seq_line).split(";")[-2]
                if sotu_of_zotu not in removed_sotus:
                    write_seq = True
                    filtered_zotu_fasta_fio.write(zotu_fasta_seq_line)
                else:
                    write_seq = False
            elif write_seq:
                filtered_zotu_fasta_fio.write(zotu_fasta_seq_line)
            zotu_fasta_seq_line = zotu_fasta_fio.readline()

    system("mv " + filtered_zotu_fasta_path + " " + ZOTU_fasta_file_path)


OUTPUT_FOLDER = argv[1]
OUTPUT_ASV_FASTA_WITH_TAXONOMY = argv[2]
OUTPUT_ASV_TABLE = argv[3]
CLUSTERING_DIRECTORY = argv[4]
INPUT_FASTA_CLUSTERING = argv[5]
KRONA_TOOL = argv[6]
OUTPUT_SOTU_FASTA_WITH_TAXONOMY = argv[7]
SILVA_ARB = argv[8]
SINA_EXECUTABLE = argv[9]
USER_FASTQ_FOLDER = argv[10]
RAPID_NJ_BIN = argv[11]
SOTUS_TABLE_NAME = argv[12]
INITIAL_ZOTUS_TABLE = argv[13]
ABUNDANCE_CUTOFF = argv[14]
SAMPLE_WISE_CORR_FLAG = argv[15]

mkdir(OUTPUT_FOLDER)
zotus_seqs_dict = fasta2dict(INPUT_FASTA_CLUSTERING)
all_stats = glob(CLUSTERING_DIRECTORY + '/species_stats/**/*.stats', recursive=True)
output_fasta = open(OUTPUT_FOLDER + '/' + OUTPUT_ASV_FASTA_WITH_TAXONOMY, 'w+')
output_fasta_sotu_centroids = open(OUTPUT_FOLDER + '/' + OUTPUT_SOTU_FASTA_WITH_TAXONOMY, 'w+')
for curr_stats in all_stats:
    stats_contents = read_file(curr_stats)
    for line in stats_contents:
        if not line:
            continue
        zotu_with_taxonomy = line.split('\t')[1]
        clean_zotu_name = zotu_with_taxonomy.split(';')[0]
        output_fasta.write('>' + zotu_with_taxonomy + '\n')
        output_fasta.write(zotus_seqs_dict[clean_zotu_name] + '\n')
        status = line.split('\t')[0]
        if status == 'S':
            # The format of header would be like -> >SOTU1;tax=X1;X2;X3;X4;X5;X6;SOTU1;Zotu2;
            fasta_header_list = zotu_with_taxonomy.split(";")
            sotu_name = fasta_header_list[-2]
            fasta_header_list[0] = sotu_name
            fasta_header_list[-1] = clean_zotu_name
            new_fasta_header = str(";".join(fasta_header_list)) + ";"
            output_fasta_sotu_centroids.write('>' + new_fasta_header + '\n')
            output_fasta_sotu_centroids.write(zotus_seqs_dict[clean_zotu_name] + '\n')
output_fasta_sotu_centroids.close()
output_fasta.close()

zotus_with_taxonomy_contents = read_file(OUTPUT_FOLDER + '/' + OUTPUT_ASV_FASTA_WITH_TAXONOMY)
zotus_taxonomy_dict = dict()
for line in zotus_with_taxonomy_contents:
    if line[0] == '>':
        tokens = line.split('tax=')
        zotu_name = tokens[0].split(';')[0][1:]
        taxonomy = tokens[1]
        zotus_taxonomy_dict[zotu_name] = taxonomy


##############################################
# Updating ZOTUs-Table with new TIC taxonomy #
##############################################

out_tab = open(OUTPUT_FOLDER + '/' + OUTPUT_ASV_TABLE, 'w+')
# INITIAL_ZOTUS_TABLE = ZOTU_map.tab is the combined zotus table of selecetd samples by the user
tab_contents = read_file(ospath.join(getcwd(), INITIAL_ZOTUS_TABLE))
for line in tab_contents:
    if line[0] == '#':
        out_tab.write(line + '\tTaxonomy\n')
    clean_zotu_name = line.split('\t')[0]
    if clean_zotu_name in zotus_taxonomy_dict.keys():
        taxonomy = str(zotus_taxonomy_dict[clean_zotu_name]).strip()
        new_line = "\t".join(line.split("\t")[:-1]) + '\t' + taxonomy + '\n'
        out_tab.write(new_line)
out_tab.close()

"""
taxonomy_counters_dict = dict()
for value in zotus_taxonomy_dict.values():
    if value not in taxonomy_counters_dict.keys():
        taxonomy_counters_dict[value] = 1
    else:
        taxonomy_counters_dict[value] += 1

out_file = open(OUTPUT_FOLDER + '/for_krona.tab', 'w+')
for key, value in taxonomy_counters_dict.items():
    out_file.write(str(value) + '\t' + key.replace(';', '\t') + '\n')
out_file.close()

system(KRONA_TOOL + ' ' + OUTPUT_FOLDER + '/for_krona.tab 2>>' + USER_FASTQ_FOLDER + '/log_file.txt' + ' 1>>' + USER_FASTQ_FOLDER + '/log_file.txt')
system('mv text.krona.html ' + OUTPUT_FOLDER + '/krona_plot.html')
system('rm ' + OUTPUT_FOLDER + '/for_krona.tab')

def create_annotation():
    annot_1 = open(OUTPUT_FOLDER + '/annot.txt', 'w+')
    guide = open(OUTPUT_FOLDER + '/guide.txt', 'w+')
    phyla = list()
    classes = list()
    orders = list()
    families = list()
    for i in taxonomy_counters_dict.keys():
        tokens = i.split(';')
        cut_taxonomy = '.'.join([tokens[0], tokens[1], tokens[2], tokens[3], tokens[4]])
        guide.write(cut_taxonomy + '\n')
        curr_phylum = tokens[0] + '.' + tokens[1]
        if curr_phylum not in phyla:
            phyla.append(curr_phylum)
        curr_class = curr_phylum + '.' + tokens[2]
        if curr_class not in classes:
            classes.append(curr_class)
        curr_order = curr_class + '.' + tokens[3]
        if curr_order not in orders:
            orders.append(curr_order)
        curr_family = curr_order + '.' + tokens[4]
        if curr_family not in families:
            families.append(curr_family)

    annot_1.write('title_font_size\t33\n')
    annot_1.write('total_plotted_degrees\t340\n')
    annot_1.write('annotation_background_alpha\t0.1\n')
    annot_1.write('start_rotation\t270\n')
    annot_1.write('internal_label\t1\tDomain\n')
    annot_1.write('internal_label\t2\tPhyla\n')
    annot_1.write('internal_label\t3\tClasses\n')
    annot_1.write('internal_label\t4\tOrders\n')
    annot_1.write('internal_label\t5\tFamilies\n')
    annot_1.write('internal_labels_rotation\t270\n')
    for i in phyla:
        phrase_1 = i + '\tclade_marker_shape\th\n'
        annot_1.write(phrase_1)
        if 'UNK' in i:
            phrase_1 = i + '\tclade_marker_color\tred\n'
            annot_1.write(phrase_1)
    for i in classes:
        phrase_1 = i + '\tclade_marker_shape\tp\n'
        annot_1.write(phrase_1)
        if 'UNK' in i:
            phrase_1 = i + '\tclade_marker_color\tred\n'
            annot_1.write(phrase_1)
    for i in orders:
        phrase_1 = i + '\tclade_marker_shape\td\n'
        annot_1.write(phrase_1)
        if 'UNK' in i:
            phrase_1 = i + '\tclade_marker_color\tred\n'
            annot_1.write(phrase_1)
    for i in families:
        phrase_1 = i + '\tclade_marker_shape\ts\n'
        annot_1.write(phrase_1)
        if 'FOTU' in i:
            phrase_1 = i + '\tclade_marker_color\tred\n'
            annot_1.write(phrase_1)
    annot_1.close()
    guide.close()

# We do not produce step.png anymore
# create_annotation()
# cmd = "graphlan_annotate --annot " + OUTPUT_FOLDER + "/annot.txt " + OUTPUT_FOLDER + '/guide.txt '
# cmd += OUTPUT_FOLDER + "/guide.xml"
# system(cmd)
# system("graphlan " + OUTPUT_FOLDER + "/guide.xml " + OUTPUT_FOLDER + "/step.png --dpi 300 --size 6.5")
# system("rm " + OUTPUT_FOLDER + "/annot.txt " + OUTPUT_FOLDER + '/guide.txt ' + OUTPUT_FOLDER + "/guide.xml")
"""

########################
# Updating SOTUs-Table #
########################

sotu_zotu_map_dict = dict()
zotus_taxonomy = read_file(OUTPUT_FOLDER + '/' + OUTPUT_ASV_TABLE)[1:]
for line in zotus_taxonomy:
    tokens = line.split('\t')
    taxonomy = tokens[-1]
    sotu = taxonomy.split(';')[-2]
    zotu = tokens[0]
    if sotu in sotu_zotu_map_dict.keys():
        sotu_zotu_map_dict[sotu] = sotu_zotu_map_dict[sotu] + ',' + zotu
    else:
        sotu_zotu_map_dict[sotu] = zotu

taxonomies = list()
for line in zotus_taxonomy:
    tokens = line.split('\t')
    taxonomy = tokens[-1]
    if taxonomy not in taxonomies:
        taxonomies.append(taxonomy)

out_file = open(OUTPUT_FOLDER + '/' + SOTUS_TABLE_NAME, 'w+')
header = read_file(OUTPUT_FOLDER + '/' + OUTPUT_ASV_TABLE)[0]
out_file.write(header.replace("#Zotu", "#SOTU") + '\n')
for taxonomy in taxonomies:
    out_line = taxonomy.split(';')[-2]
    curr_taxonomy_sample_sizes = list()
    for line in zotus_taxonomy:
        tokens = line.split('\t')
        curr_taxonomy = tokens[-1]
        if curr_taxonomy == taxonomy:
            samples_reads = '\t'.join(tokens[1:-1])
            curr_taxonomy_sample_sizes.append(samples_reads)
    samples_num = curr_taxonomy_sample_sizes[0].count('\t') + 1
    for i in range(samples_num):
        athroisma = 0
        for line in curr_taxonomy_sample_sizes:
            token = line.split('\t')[i]
            athroisma += int(token)
        out_line += '\t' + str(athroisma)
    out_line += '\t' + taxonomy + '\n'
    out_file.write(out_line)
out_file.close()


#######################
# Filtering Abundance #
#######################

SAMPLE_WISE_CORR_FLAG = True if SAMPLE_WISE_CORR_FLAG == "True" else False
filter_abundance(
    abundance_limit=ABUNDANCE_CUTOFF,
    OTU_table_file_name=OUTPUT_FOLDER + "/" + SOTUS_TABLE_NAME,
    OTU_fasta_file_name=OUTPUT_FOLDER + "/" + OUTPUT_SOTU_FASTA_WITH_TAXONOMY,
    ZOTU_table_file_name=OUTPUT_FOLDER + "/" + OUTPUT_ASV_TABLE,
    ZOTU_fasta_file_name=OUTPUT_FOLDER + "/" + OUTPUT_ASV_FASTA_WITH_TAXONOMY,
    sample_wise_corr_flag=SAMPLE_WISE_CORR_FLAG
)


######################
# Creating Map Files #
######################


out_file = open(OUTPUT_FOLDER + '/Map-ZOTU-SOTU.tab', 'w+')
out_file.write('ZOTUs\tSOTUs\n')
contents = read_file(OUTPUT_FOLDER + '/' + OUTPUT_ASV_TABLE)[1:]
for line in contents:
    taxonomy = line.split('\t')[-1]
    zotu = line.split('\t')[0]
    sotu = taxonomy.split(';')[-2]
    out_file.write(zotu + '\t' + sotu + '\n')
out_file.close()

out_file = open(OUTPUT_FOLDER + '/Map-SOTU-GOTU.tab', 'w+')
out_file.write('SOTUs\tGOTUs\n')
contents = read_file(OUTPUT_FOLDER + '/' + OUTPUT_ASV_TABLE)[1:]
passed_sotus = list()
for line in contents:
    taxonomy = line.split('\t')[-1]
    sotu = taxonomy.split(';')[-2]
    gotu = taxonomy.split(';')[-3]
    if sotu not in passed_sotus:
        out_file.write(sotu + '\t' + gotu + '\n')
        passed_sotus.append(sotu)
out_file.close()

out_file = open(OUTPUT_FOLDER + '/Map-GOTU-FOTU.tab', 'w+')
out_file.write('GOTUs\tFOTUs\n')
contents = read_file(OUTPUT_FOLDER + '/' + OUTPUT_ASV_TABLE)[1:]
passed_gotus = list()
for line in contents:
    taxonomy = line.split('\t')[-1]
    gotu = taxonomy.split(';')[-3]
    fotu = taxonomy.split(';')[-4]
    if gotu not in passed_gotus:
        out_file.write(gotu + '\t' + fotu + '\n')
        passed_gotus.append(gotu)
out_file.close()

#################
# Tree Creation #
#################


def sina_alignment():
    # Aligning zotus
    cmd = SINA_EXECUTABLE + ' --in=' + OUTPUT_FOLDER + "/" + OUTPUT_ASV_FASTA_WITH_TAXONOMY
    cmd += ' --out=' + OUTPUT_FOLDER + '/test_z.fasta --db='
    cmd += SILVA_ARB + ' --turn all --fasta-write-dna >/dev/null 2>/dev/null'
    system(cmd)
    # Aligning sotus
    cmd = SINA_EXECUTABLE + ' --in=' + OUTPUT_FOLDER + "/" + OUTPUT_SOTU_FASTA_WITH_TAXONOMY
    cmd += ' --out=' + OUTPUT_FOLDER + '/test_s.fasta --db='
    cmd += SILVA_ARB + ' --turn all --fasta-write-dna >/dev/null 2>/dev/null'
    system(cmd)


# from the alignment get every column that has 1 base at least aligned
# so tree creation is correct and goes faster
def extract_columns(input_file_name):
    good_points = list()
    filepath = input_file_name
    with open(filepath) as fp:
        line = fp.readline()
        while line:
            if not line[0] == '>':
                counter = 0
                for curr_char in line:
                    if curr_char in ['A', 'C', 'G', 'T', 'N']:
                        good_points.append(counter)
                    counter += 1
            line = fp.readline()
    #
    good_points = list(set(good_points))
    #
    out_file = open(OUTPUT_FOLDER + '/for_tree.fasta', 'w+')
    with open(filepath) as fp:
        line = fp.readline()
        while line:
            if line[0] == '>':
                header = line.split(';')[0] + '\n'
            else:
                curr_seq = ''
                for number in good_points:
                    curr_seq += line[number]
                out_file.write(header)
                out_file.write(curr_seq + '\n')
            line = fp.readline()
    out_file.close()
    system('mv ' + OUTPUT_FOLDER + '/for_tree.fasta ' + input_file_name)


def create_trees():
    extract_columns(OUTPUT_FOLDER + '/test_s.fasta')
    extract_columns(OUTPUT_FOLDER + '/test_z.fasta')
    rapidnj = RAPID_NJ_BIN + " "
    cmd1 = rapidnj + OUTPUT_FOLDER + "/test_s.fasta -n -o t -x " + OUTPUT_FOLDER + "/SOTUs-Tree-nj.tre"
    cmd2 = rapidnj + OUTPUT_FOLDER + "/test_z.fasta -n -o t -x " + OUTPUT_FOLDER + "/ZOTUs-Tree-nj.tre"
    # cmd3 = "/crc/crc/binaries/FastTree -gtr -gamma -quiet -nt < test_s.fasta > sotu_aml.tre"
    # cmd4 = "/crc/crc/binaries/FastTree -gtr -gamma -quiet -nt < test_z.fasta > zotu_aml.tre"
    system(cmd1)
    system(cmd2)
    rem_quot_cmd = r"sed -i s/\'//g {}"
    system(rem_quot_cmd.format(OUTPUT_FOLDER + "/SOTUs-Tree-nj.tre"))
    system(rem_quot_cmd.format(OUTPUT_FOLDER + "/ZOTUs-Tree-nj.tre"))
    # print(cmd3)
    # print(cmd4)


sina_alignment()
create_trees()
print("tree created")


def fasta_name_with_space(filepath):
    filepath = ospath.abspath(filepath)
    filedir, filename = ospath.split(filepath)
    fout_path = ospath.join(filedir, "placeholder" + filename)
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
    system("mv {} {}".format(fout_path, filepath))


for fasta_file in [OUTPUT_ASV_FASTA_WITH_TAXONOMY, OUTPUT_SOTU_FASTA_WITH_TAXONOMY]:
    fasta_name_with_space(OUTPUT_FOLDER + '/' + fasta_file)


def remove_last_semi(table_file):
    lines = read_file(table_file)
    with open(table_file + "TO_REMOVE", "w+") as without_semi:
        without_semi.write(str(lines[0]).strip() + "\n")
        for line in lines[1:]:
            line_w = str(line).strip()
            line_w = line_w[:-1] if line_w.endswith(";") else line_w
            without_semi.write(line_w + "\n")
    system("mv {} {}".format(table_file + "TO_REMOVE", table_file))


for table_file in [OUTPUT_ASV_TABLE,
                   SOTUS_TABLE_NAME]:
    remove_last_semi(OUTPUT_FOLDER + '/' + table_file)
