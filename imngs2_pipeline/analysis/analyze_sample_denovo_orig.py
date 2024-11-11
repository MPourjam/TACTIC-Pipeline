from os import chdir, system
from crcapp.models import mimarks
from crcapp.models_job import job_anaquery
# from django.db import models as djmodels
# from django.contrib.auth.models import User # For excluding User models while recurFieldValgetter()
# from crcapp.recursiveField import create_mapping_file_n
# import pandas as pd
# import glob
# import re
import os


ORIGINAL_DIR = '/srv/crc/projects/'
BIN_DIR = "/crc/crc/binaries/"
USEARCH_11_bin = BIN_DIR + "usearch11.0.667_i86linux32"
SINA_ARB = '/srv/crc/databases/SILVA_138_SSURef_NR99_05_01_20_opt.arb'
GOLD_REFDB_USEARCH = '/srv/crc/databases/SILVA-bac-16s-90.udb'
USERS_DIR = '/srv/crc/users/'
ref16RNAdb_1 = "/srv/crc/databases/silva-bac-16s-id90.fasta"
ref16RNAdb_2 = "/srv/crc/databases/silva-arc-16s-id95.fasta"
SORT_ME_RNA_BIN = '/crc/crc/binaries/sortmerna'
USEARCH_TAIL = '> /dev/null 2>&1'


def read_file(filename):
    with open(filename) as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    return content


def onelinefasta(fastafilepath):
    proper_filepath = os.path.abspath(fastafilepath)
    dirpath, filename = os.path.split(proper_filepath)
    os.chdir(dirpath)
    with open(proper_filepath) as f:
        content = f.readlines()
        rawcontent = [x.strip() for x in content]
    newfile_temp_name = "oneline_{}.fasta".format(filename.split(".")[0])
    assert (not os.path.isfile(newfile_temp_name)), "Destination file {} already exists".format(newfile_temp_name)
    with open(newfile_temp_name, "w+") as onelinefa:
        for i in range(0, len(rawcontent)):
            if rawcontent[i][0] == ">":
                if i == 0:
                    onelinefa.write(rawcontent[i] + "\n")
                else:
                    onelinefa.write("\n" + rawcontent[i] + "\n")
            else:
                onelinefa.write(rawcontent[i])
    os.system('mv {fbase} {f}'.format(**{"fbase": newfile_temp_name, "f": filename}))


def append_reads(sample, analysis_dir):
    reads = ORIGINAL_DIR + str(sample.sample_f.organism_f.project_f.id) + '/'
    reads += str(sample.sample_f.organism_f.id) + '/'
    reads += str(sample.sample_f.id) + '/16S/' + str(sample.id) + '/taxed_ZOTUs.fasta'
    print(reads)
    system(USEARCH_11_bin + ' -fastx_relabel ' + reads + ' -prefix D' + str(sample.id) + '. -fastaout new -keep_annots')
    system('cat new >> ' + analysis_dir + 'analysis.fasta')
    onelinefasta('{}analysis.fasta'.format(analysis_dir))


def trim_sides(five_end_trim, three_end_trim):
    cmd_0 = USEARCH_11_bin + " -fastx_truncate analysis.fasta -stripright " + str(five_end_trim)
    cmd_1 = " -stripleft " + str(three_end_trim) + " -fastaout filtered1.fasta "
    system(cmd_0 + cmd_1 + USEARCH_TAIL)
    cmd_2 = 'mv filtered1.fasta analysis.fasta'
    system(cmd_2)
    onelinefasta('analysis.fasta')


def dereplication():
    cmd_part_0 = USEARCH_11_bin + ' -fastx_uniques analysis.fasta'
    cmd_part_1 = ' -sizein -sizeout -fastaout derep.fasta'
    cmd_part_2 = ' >/dev/null 2>/dev/null'
    # dereplicate reads and anotate them by multiplicity
    system(cmd_part_0 + cmd_part_1 + cmd_part_2)
    onelinefasta("derep.fasta")


def sort_seqs():
    cmd_0 = USEARCH_11_bin + ' -sortbysize derep.fasta -fastaout sorted.fasta'
    cmd_1 = ' ' + USEARCH_TAIL
    system(cmd_0 + cmd_1)
    onelinefasta("sorted.fasta")


# cluster sequences to OTUs
def clusterZOTUs():
    cmd_part_0 = USEARCH_11_bin + ' -unoise3 sorted.fasta -minsize 4 -zotus zotus.fasta'
    cmd_part_1 = ' ' + USEARCH_TAIL
    # clusering of seq in OTUs (clustered)
    system(cmd_part_0 + cmd_part_1)
    system('rm derep.fasta sorted.fasta')
    onelinefasta("zotus.fasta")


# Filter non 16S sequences
def filter16S():
    cmd_0 = SORT_ME_RNA_BIN + " --ref " + ref16RNAdb_1 + " --ref " + ref16RNAdb_2 + " --reads "
    cmd_1 = " zotus.fasta --fastx --aligned good-ZOTUs --other non16SrRNA --workdir . -e 0.1"  # --aligned
    cmd_2 = " > /dev/null 2>&1"
    # run a RNA filtering step
    # (The program currently do not distinquish between 16S and 18S)
    system(cmd_0 + cmd_1 + cmd_2)
    # system('mv out/aligned.fasta good-ZOTUs.fasta')
    system('rm -r idx kvdb')


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
    cmd_part_0 = USEARCH_11_bin + ' -cluster_otus good-ZOTUs-sized.fasta -fulldp -otus otus1.fa -minsize 1 '
    cmd_part_1 = ' -uparseout z2o.tab'
    cmd_part_2 = ' >/dev/null 2>/dev/null'
    # clusering of seq in OTUs (clustered)
    system(cmd_part_0 + cmd_part_1 + cmd_part_2)
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
    system('rm otus1.fa')


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
    system('rm z2o.tab')


def keep_good_ZOTUs():
    cmd_0 = USEARCH_11_bin + ' -fastx_getseqs good-ZOTUs.fasta -labels matched_ZOTUS.txt '  # aligned_16S_ZOTUS.fa -> good-ZOTUs.fasta
    cmd_1 = ' -fastaout nochi_ZOTUs.fasta '
    system(cmd_0 + cmd_1 + USEARCH_TAIL)
    onelinefasta("nochi_ZOTUs.fasta")


def build_ZOTU_table():
    cmd_0 = USEARCH_11_bin + ' -otutab analysis.fasta -top_hit_only -zotus nochi_ZOTUs.fasta '
    cmd_1 = ' -otutabout zotu_table.txt -id 0.97'
    system(cmd_0 + cmd_1 + USEARCH_TAIL)


def get_samples_sizes(input_file):
    contents = read_file(input_file)[1:]
    tot_sizes_list = list()
    num_tokens = len(contents[0].split('\t')) - 1
    for i in range(0, num_tokens):
        tot_sizes_list.append(0)
    for line in contents:
        line_tokens = line.split('\t')[1:]
        for i in range(len(line_tokens)):
            tot_sizes_list[i] += float(line_tokens[i])
    return tot_sizes_list


def filter_zotu_abundance(abund_limit):
    unf_tot_sizes = get_samples_sizes('zotu_table.txt')
    contents = read_file('zotu_table.txt')
    sample_sizes = list()
    datasets_samples = len(contents[0].split('\t')) - 1
    for i in range(0, datasets_samples):
        sample_sizes.append(0)
    out_file = open('zotu_table_filtered.txt', 'w+')
    out_file_2 = open('filtered_zotu_table_list.txt', 'w+')
    for line in contents[1:]:
        curr_sizes = line.split('\t')[1:]
        abundances_list = list()
        for j in range(len(curr_sizes)):
            abundance = float(curr_sizes[j]) / unf_tot_sizes[j]
            abundances_list.append(abundance)
        if any(x >= abund_limit for x in abundances_list):
            out_file.write(line + '\n')
            out_file_2.write(line.split('\t')[0] + '\n')
        del abundances_list[:]
    out_file.close()
    out_file_2.close()


def select_zotu_seqs():
    cmd_0 = USEARCH_11_bin + ' -fastx_getseqs good-ZOTUs.fasta -labels filtered_zotu_table_list.txt '
    cmd_1 = ' -fastaout ZOTUs-Seqs.fasta '
    system(cmd_0 + cmd_1 + USEARCH_TAIL)
    system('rm filtered_zotu_table_list.txt')
    onelinefasta("ZOTUs-Seqs.fasta")


def addTax():
    classifier_dir = '/crc/crc/binaries/sina/'
    cmd_part_0 = 'sina --in ZOTUs-Seqs.fasta --search --meta-fmt csv '
    cmd_part_1 = '--threads 4 --lca-fields tax_slv '
    cmd_part_2 = '--db ' + SINA_ARB + ' --out test.fasta'
    cmd_part_3 = ' >/dev/null 2>/dev/null'
    system(classifier_dir + cmd_part_0 + cmd_part_1 + cmd_part_2 + cmd_part_3)
    out_file = open('classifiedF.txt', 'w+')
    silva_contents_header = read_file('test.csv')
    silva_contents = silva_contents_header[1:]
    for line in silva_contents:
        line_tokens = line.split(',')
        silva_taxo = line_tokens[6]
        seq_name = line_tokens[0]
        out_file.write(seq_name + '\t' + silva_taxo + '\n')
    out_file.close()


def addTax_nochi(nochi_zotu_path):
    nochi_zotu_path = os.path.abspath(nochi_zotu_path)
    dirname, filename = os.path.split(nochi_zotu_path)
    os.chdir(dirname)
    classifier_dir = '/crc/crc/binaries/sina/'
    cmd_part_0 = 'sina --in ' + str(nochi_zotu_path) + ' --search --meta-fmt csv '
    cmd_part_1 = '--threads 4 --lca-fields tax_slv '
    cmd_part_2 = '--db ' + SINA_ARB + ' --out test_nochi.fasta'
    cmd_part_3 = ' >/dev/null 2>/dev/null'
    system(classifier_dir + cmd_part_0 + cmd_part_1 + cmd_part_2 + cmd_part_3)

    taxed_nochi_path = os.path.join(dirname, "taxed_nochi_ZOTUs.fasta")
    out_file = open(taxed_nochi_path, 'w+')
    nochi_seqs = read_file(str(nochi_zotu_path))
    zotunames = []
    zotuseqs = []
    for ll in nochi_seqs:
        if ll[0] == ">":
            zotunames.append(ll[1:])
        else:
            zotuseqs.append(ll)
    zotu_seq_dict = dict(zip(zotunames, zotuseqs))
    test_nochi_path = os.path.join(dirname, 'test_nochi.csv')
    silva_contents_header = read_file(test_nochi_path)
    silva_contents = silva_contents_header[1:]
    for line in silva_contents:
        line_tokens = line.split(',')
        # print(line_tokens)
        silva_taxo = line_tokens[6]
        seq_name = line_tokens[0]
        out_file.write('>' + seq_name + ';' + "tax=" + silva_taxo + '\n' + zotu_seq_dict[seq_name] + '\n')
    out_file.close()

    # out_file = open('classifiedF.txt', 'w+')
    # silva_contents_header = read_file('test.csv')
    # silva_contents = silva_contents_header[1:]
    # for line in silva_contents:
    #     line_tokens = line.split(',')
    #     silva_taxo = line_tokens[6]
    #     seq_name = line_tokens[0]
    #     out_file.write(seq_name + '\t' + silva_taxo + '\n')
    # out_file.close()


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
    # create both trees
    cmd_part_0 = 'FastTree -quiet -nosupport -gtr -nt test_otu.fasta'
    cmd_part_1 = ' > ' + sample + '_FastTree.tre'
    system(cmd_part_0 + cmd_part_1)
    cmd_2_part_0 = BIN_DIR + 'rapidnj test_otu.fasta --input-format fa --cores 6 --alignment-type d '
    cmd_2_part_1 = '--output-format t --no-negative-length -x ' + sample + '_rapidNJ_tree.tre'
    system(cmd_2_part_0 + cmd_2_part_1)


def cleanup():
    deleted_files = 'test.fasta test_otu.fasta test.csv good-ZOTUs-sized.fasta new  non16SrRNA.fasta '  # analysis.fasta
    deleted_files += 'good-ZOTUs.log '  # zotus.fasta mOTUs-Seqs.fasta classifiedF.txt good-ZOTUs.fasta
    deleted_files += 'matched_ZOTUS.txt ZOTUs-Seqs.fasta zotu_table_filtered.txt '  # zotu_table.txt
    system('rm ' + deleted_files)


def create_zip(user_name, query_id):
    analysis_dir = USERS_DIR + user_name + '/AnalQuery/'
    chdir(analysis_dir)
    cmd_zip = 'zip -r ' + query_id + '.zip ' + query_id
    system(cmd_zip)


def create_krona():
    OTUfinal_line = read_file('ZOTUs-table.final.tab')[1:]
    krona_text = open('krona.txt', 'w+')
    for OTU_taxonomy in OTUfinal_line:
        OTUFinal_list = OTU_taxonomy.split('\t')[1:]
        taxo = '\t'.join(OTUFinal_list[-1].split(";")[:-1])
        size = OTUFinal_list[0]
        krona_text.write(size + '\t' + taxo + '\n')
    krona_text.close()
    system("ktImportText krona.txt")
    system("rm krona.txt")


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
                out_line_2 = zotu_map_otus[zotu_name] + '\t' + "\t".join(otu_summed_size_dict[zotu_map_otus[zotu_name]]) + '\t'
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
    with open('test_otu.fasta', "w+") as test_otu:
        for ll in otus_contents:
            if ll[0] == ">":
                zotu = ll[1:-1]
                if zotu in zotu_map_otus.keys():
                    otu_name = zotu_map_otus[zotu]
                    test_otu.write(">" + otu_name + '\n')
                    write_seq = zotu
                else:
                    write_seq = False
            else:
                if write_seq:
                    test_otu.write(testfasta_zotus_dict[write_seq] + '\n')


def create_mapping_file(selected_samples):

    dataset_fields = ['mimarks_name', 'accession', 'dna_isolation', 'sequencer', 'run_length',
                      'amplification_steps', 'step_1_PCR', 'step_2_PCR']
    dataset_sequencing_metadata_fields = ['primers_used', 'paired_end', 'target_region']
    sample_fields = ['source_f', 'name', 'description', 'age', 'age_unit', 'weight', 'weight_unit',
                     'collection_date', 'collection_time', 'collection_country',
                     'collection_location_gps', 'preservation']
    sample_optional_mouse_fields = [
        'feed_provider',
        'diet_category',
        'animal_facility',
        'housing_hygiene_level',
        'caging',
        'basal_microbiota',
        'biotic_challenge',
        'abiotic_challenge'
    ]
    sample_optional_human_fields = [
        'cancer_related_symptoms',
        'arterial_hypertension',
        'hypercholesterolemia',
        'smoking',
        'alcohol_dependance',
        'physical_activity',
        'regular_medication',
        'regular_medication_categories',
        'antibiotics',
        'probiotics',
        'supplements',
        'bristol_score',
        'tissue_available',
        'tissue_type'
    ]
    organism_fields = ['name', 'species', 'description']
    organism_optional_human_fields = ['place_of_birth', 'medical_history', 'ibd', 'cancer']
    organism_optional_mouse_fields = ['genetic_modification', 'genotype']
    project_fields = ['project_acronym', 'project_description', 'project_name', 'project_accession', 'create_date']
    out_lines = list()
    out_file = open('mapping.tab', 'w+')

    for curr_dataset in selected_samples:
        dataset_id = 'D' + str(curr_dataset.id)
        curr_metadata_list = list()
        for curr_field in dataset_fields:
            curr_metadata_list.append(str(getattr(curr_dataset, curr_field, 'NA')))
        curr_seq_meta = curr_dataset.sequencing_metadata_f
        for curr_field in dataset_sequencing_metadata_fields:
            curr_metadata_list.append(str(getattr(curr_seq_meta, curr_field, 'NA')))
        print('dataset mapping complete')
        #
        curr_sample = curr_dataset.sample_f
        for curr_field in sample_fields:
            curr_metadata_list.append(str(getattr(curr_sample, curr_field, 'NA')))
        curr_sample_optional_human = curr_sample.sample_human_optional_metadata_f
        for curr_field in sample_optional_human_fields:
            curr_metadata_list.append(str(getattr(curr_sample_optional_human, curr_field, 'NA')))
        curr_sample_optional_mouse = curr_sample.sample_mouse_optional_metadata_f
        for curr_field in sample_optional_mouse_fields:
            curr_metadata_list.append(str(getattr(curr_sample_optional_mouse, curr_field, 'NA')))
        print('sample mapping complete')
        #
        curr_organism = curr_sample.organism_f
        for curr_field in organism_fields:
            curr_metadata_list.append(str(getattr(curr_organism, curr_field, 'NA')))
        curr_organism_optional_human = curr_organism.human_metadata_f
        for curr_field in organism_optional_human_fields:
            curr_metadata_list.append(str(getattr(curr_organism_optional_human, curr_field, 'NA')))
        curr_organism_optional_mouse = curr_organism.mouse_metadata_f
        for curr_field in organism_optional_mouse_fields:
            curr_metadata_list.append(str(getattr(curr_organism_optional_mouse, curr_field, 'NA')))
        print('organism mapping complete')
        #
        curr_project = curr_organism.project_f
        for curr_field in project_fields:
            curr_metadata_list.append(str(getattr(curr_project, curr_field, 'NA')))
        curr_line = dataset_id + '\t' + '\t'.join(curr_metadata_list) + '\n'
        out_lines.append(curr_line)
        del curr_metadata_list[:]
    #  change name to sample_name cant do it before because i need it for the getattr function
    dataset_fields = ['mimarks_name', 'accession', 'dna_isolation', 'sequencer', 'run_length',
                      'amplification_steps', 'step_1_PCR', 'step_2_PCR']
    dataset_sequencing_metadata_fields = ['primers_used', 'paired_end', 'target_region']
    sample_fields = ['source_f', 'sample_name', 'description', 'age', 'age_unit', 'weight', 'weight_unit',
                     'collection_date', 'collection_time', 'collection_country',
                     'collection_location_gps', 'preservation']
    organism_fields = ['organism_name', 'species', 'description']
    header = '##DatasetID\t' + '\t'.join(dataset_fields) + '\t' + '\t'.join(dataset_sequencing_metadata_fields) + '\t' + '\t'.join(sample_fields) + '\t'
    header += '\t'.join(sample_optional_human_fields) + '\t' + '\t'.join(sample_optional_mouse_fields) + '\t' + '\t'.join(organism_fields) + '\t'
    header += '\t'.join(organism_optional_human_fields) + '\t' + '\t'.join(organism_optional_mouse_fields) + '\t' + '\t'.join(project_fields) + '\n'
    out_lines = [header] + out_lines
    for line in out_lines:
        out_file.write(line)
    out_file.close()


def main(user_name, analysis_q_id):
    ANALYSIS_DIR = USERS_DIR + user_name + '/AnalQuery/' + analysis_q_id + '/'
    chdir(ANALYSIS_DIR)
    selected_samples_ids = read_file(USERS_DIR + user_name + '/AnalQuery/selected_samples.txt')
    selected_mimarks = mimarks.objects.filter(pk__in=selected_samples_ids)
    print('# Gathering Sequences')
    for sample in selected_mimarks:
        append_reads(sample, ANALYSIS_DIR)
    chdir(ANALYSIS_DIR)
    analysis_obj = job_anaquery.objects.get(id=int(analysis_q_id))
    print('# Trim Sides')
    trim_sides(analysis_obj.five_end_trim, analysis_obj.three_end_trim)
    print('# Dereplication')
    dereplication()
    print('# Sort Sequences')
    sort_seqs()
    print('# Cluster ZOTUs')
    clusterZOTUs()
    print('# Filter 16S')
    filter16S()
    prepare_zotus()
    print('# Cluster OTUs')
    clusterOTUs()
    remove_size_from_OTUS()
    assign_zotus_to_otus()
    print('# Filter Chimeras')
    keep_good_ZOTUs()
    print('# Build ZOTU Table')
    build_ZOTU_table()
    print('# Filter ZOTUs by Abundance')
    filter_zotu_abundance(analysis_obj.abundance)
    select_zotu_seqs()
    # print('Select ZOTUs Filtered')
    addTax()
    print('# Add Taxonomy')
    create_final_ZOTU_table()
    add_taxonomy_to_fasta()
    # print('Taxonomy Added')
    print('# Analysis Cleanup')
    create_OTUs_from_ZOTUS_table()
    # create_mapping_file(selected_mimarks) # It now gets produced in views_job.py analysis_job()
    # print('Mapping File Created')
    print('# Krona Creation')
    create_krona()
    print('# Tree Creation')
    create_tree(analysis_q_id)
    cleanup()
    print('Cleanup DONE')
    create_zip(user_name, analysis_q_id)
    # print('ANALYSIS DONE')
