import argparse
from os import chdir, mkdir, remove
from os.path import isfile, getsize, isdir
from shutil import copyfile
from tqdm import tqdm
import multiprocessing
from glob import glob
from random import shuffle
from processing_helper import gimmelogger
from processing_helper import system_sub as sys_sub


global TIC_LOG
TIC_LOG = gimmelogger(
    logger_name="run_imngs2.analysis.TIC",
    # log_file=ANALYSIS_DIR.joinpath("Analysis_log.txt"),
    only_file=True
)


def system_sub(*args, **kwargs):
    return sys_sub(*args, logger_obj=TIC_LOG, **kwargs)


VSEARCH_TAIL = ' 1>/dev/null 2>/dev/null'
MAIN_DIR = ''


species_number = 0
species_dict = dict()
genera_number = 0
genera_dict = dict()
family_number = 0
family_dict = dict()


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
    return dic


def vsearch_clustering(curr_fasta, similarity_limit):
    start_of_curr_fasta = curr_fasta.split('.')[0]
    centroids_file = start_of_curr_fasta + '.centroids.fasta'
    uc_file = start_of_curr_fasta + '.uc'
    # cmd_0 = USEARCH_BIN_CLUST + curr_fasta + ' -id ' + similarity_limit + ' -strand both'
    # cmd_1 = ' -top_hits_only -centroids ' + centroids_file + ' -uc ' + uc_file + VSEARCH_TAIL
    # TIC_LOG.info(cmd_0 + cmd_1)
    # systm(cmd_0 + cmd_1)
    system_sub(
        [
            *str(USEARCH_BIN_CLUST).strip().split(),
            curr_fasta,
            '-id',
            similarity_limit,
            '-strand',
            'both',
            '-top_hits_only',
            '-centroids',
            centroids_file,
            '-uc',
            uc_file
        ]
    )


def my_mkdir(curr_dir):
    try:
        mkdir(curr_dir)
    except BaseException:
        return


def create_directory(curr_name, stats_or_fasta):
    tokens = curr_name.split('_')
    curr_dir = stats_or_fasta + '/' + tokens[0] + '/'
    my_mkdir(curr_dir)
    for curr_token in tokens[1:]:
        curr_dir += curr_token + '/'
        my_mkdir(curr_dir)
    return(curr_dir)


def update_species_stats(curr_fasta):
    chdir(MAIN_DIR)
    global species_number
    global species_dict
    curr_name = curr_fasta.split('.')[0]
    uc_file = curr_name + '.uc'
    if not isfile(uc_file) or not(getsize(uc_file)):
        return
    uc_content = read_file(uc_file)
    for line in uc_content:
        line_tokens = line.split('\t')
        if (line[0] == 'S'):
            species_number += 1
            curr_dir = create_directory(curr_name, 'species_stats')
            out_file_name = curr_dir + curr_name + '_SOTU' + str(species_number) + '.stats'
            out_file = open(MAIN_DIR + '/' + out_file_name, 'w+')
            genera_original_centroid_header = line_tokens[8]
            centroids_with_SOTU = genera_original_centroid_header + 'SOTU' + str(species_number) + ';'
            out_line = line_tokens[0] + '\t' + centroids_with_SOTU + '\n'
            out_file.write(out_line)
            out_file.close()
            species_dict[genera_original_centroid_header] = species_number
        elif(line[0] == 'H'):
            target_hit_species_number = str(species_dict[line_tokens[9]])
            curr_dir = create_directory(curr_name, 'species_stats')
            species_hit_with_SOTU = line_tokens[8] + 'SOTU' + target_hit_species_number + ';'
            out_file_name = curr_dir + curr_name + '_SOTU' + target_hit_species_number + '.stats'
            out_file = open(MAIN_DIR + '/' + out_file_name, 'a')
            out_line = line_tokens[0] + '\t' + species_hit_with_SOTU + '\n'
            out_file.write(out_line)
            out_file.close()
    remove(uc_file)


def update_species_centroids(curr_fasta):
    global species_dict
    curr_name = curr_fasta.split('.')[0]
    centroids_file = curr_name + '.centroids.fasta'
    if not isfile(centroids_file) or not getsize(centroids_file):
        return
    centroids_seqs_dict = fasta2dict(centroids_file)
    for key in centroids_seqs_dict.keys():
        curr_species_number = str(species_dict[key])
        curr_dir = create_directory(curr_name, 'species_centroids')
        updated_fasta_file_name = curr_dir + curr_name + '_SOTU' + curr_species_number + '.fasta'
        new_header = key.split(';')[0]
        old_taxonomy = key.split('tax=')[1]
        centroid_sequence = centroids_seqs_dict[key]
        new_taxonomy = old_taxonomy + 'SOTU' + curr_species_number + ';\n'
        updated_fasta_file = open(MAIN_DIR + '/' + updated_fasta_file_name, 'w+')
        updated_fasta_file.write('>' + new_header + ';tax=' + new_taxonomy)
        updated_fasta_file.write(centroid_sequence + '\n')
        updated_fasta_file.close()
        # TIC_LOG.info(MAIN_DIR + '/' + updated_fasta_file_name)
    remove(centroids_file)


# OTUS --> Species
def create_species_of_known_genera(known_genera_fastas):
    chdir(MAIN_DIR)
    mkdir('species_stats')
    mkdir('species_centroids')
    for i in tqdm(range(len(known_genera_fastas))):
        curr_genera_fasta = known_genera_fastas[i]
        vsearch_clustering(curr_genera_fasta, species_similarity)
        update_species_stats(curr_genera_fasta)
        update_species_centroids(curr_genera_fasta)
        remove(curr_genera_fasta)


def get_curr_level_fastas(num_underscores):
    chdir(MAIN_DIR)
    all_fastas = glob('*.fasta')
    curr_level_fastas = list()
    for curr_fasta in all_fastas:
        if (curr_fasta.count('_') == num_underscores) and ('base' not in curr_fasta):
            curr_level_fastas.append(curr_fasta)
    # TIC_LOG.info(curr_level_fastas)
    return(curr_level_fastas)


def threaded_gather_species(all_species_centroids, thread_number, curr_level, suffix):
    for curr_centroid in all_species_centroids:
        curr_centroid_taxonomy = curr_centroid.split('.')[0].split('/')[-1]
        asked_taxonomy = curr_centroid_taxonomy.split('_')[:curr_level]
        target_name_of_level = '_'.join(asked_taxonomy)
        gathered_file_name = MAIN_DIR + '/' + target_name_of_level + suffix + '.thread.' + str(thread_number)
        # system('cat ' + curr_centroid + ' >> ' + gathered_file_name)
        with open(gathered_file_name, 'a') as f:
            with open(curr_centroid, 'r') as f_centroid:
                f.write(f_centroid.read())


def gather_species(curr_level, suffix):
    chdir(MAIN_DIR + '/species_centroids')
    all_species_centroids = glob('**/*.fasta', recursive=True)
    shuffle(all_species_centroids)
    splitted_species_centroids_list = split_list(all_species_centroids, int(THREADS))
    processes = [multiprocessing.Process(target=threaded_gather_species,
                 args=(splitted_species_centroids_list[x], x, curr_level, suffix)) for x in range(int(THREADS))]
    # Run processes
    for p in processes:
        p.start()
    # Exit the completed processes
    for p in processes:
        p.join()
    chdir(MAIN_DIR)
    threaded_fastas = glob('*' + suffix + '.thread.*')
    for threaded_fasta in threaded_fastas:
        final_name = threaded_fasta.split('.thread')[0]
        # system('cat ' + threaded_fasta + ' >> ' + final_name)
        with open(final_name, 'a') as f:
            with open(threaded_fasta, 'r') as f_threaded:
                f.write(f_threaded.read())
        remove(threaded_fasta)


def vsearch_blast(curr_fasta, similarity_limit):
    chdir(MAIN_DIR)
    start_of_curr_fasta = curr_fasta.split('.')[0]
    targets_file = start_of_curr_fasta + '.base.fasta'
    b6_file = start_of_curr_fasta + '.b6'
    not_matched_file = start_of_curr_fasta + '.not_matched.fasta'
    # if targets file not exists --> family with no known genus
    # all of those sequence should go immediately to genera clustering joining the not matched
    # so just copy the file to not matched
    if not isfile(targets_file):
        copyfile(curr_fasta, not_matched_file)
    else:
        # cmd_0 = USEARCH_BIN_BLAST + curr_fasta + ' -id ' + similarity_limit + ' -db ' + targets_file
        # cmd_1 = ' -dbmask none -qmask none -strand both -blast6out ' + b6_file
        # cmd_2 = ' -notmatched ' + not_matched_file + VSEARCH_TAIL
        # system(cmd_0 + cmd_1 + cmd_2)
        system_sub(
            [
                *str(USEARCH_BIN_BLAST).strip().split(),
                curr_fasta,
                '-id',
                similarity_limit,
                '-db',
                targets_file,
                '-dbmask',
                'none',
                '-qmask',
                'none',
                '-strand',
                'both',
                '-blast6out',
                b6_file,
                '-notmatched',
                not_matched_file
            ]
        )
        remove(targets_file)
    remove(curr_fasta)


def handle_blast_matches_species(curr_fasta):
    chdir(MAIN_DIR)
    start_of_curr_fasta = curr_fasta.split('.')[0]
    b6_file = start_of_curr_fasta + '.b6'
    if not isfile(b6_file):
        return
    elif not getsize(b6_file):
        pass
    # in order not to repeat the remove command at the end of function
    # and not run anything else but clean the files at the end
    else:
        b6_contents = read_file(b6_file)
        for line in b6_contents:
            # TIC_LOG.info(line)
            query = line.split('\t')[0]
            target = line.split('\t')[1]
            target_taxonomy = target.split('tax=')[1]
            file_location_stats = '/'.join(target_taxonomy.split(';')[:-2])
            target_stats_file = file_location_stats + '/' + '_'.join(target_taxonomy.split(';')[:-1]) + '.stats'
            query_seq_identifier = query.split(';')[0]
            if not isfile(MAIN_DIR + '/species_stats/' + target_stats_file):
                continue
            target_file = open(MAIN_DIR + '/species_stats/' + target_stats_file, 'a')
            target_file.write('H\t' + query_seq_identifier + ';tax=' + target_taxonomy + '\n')
            target_file.close()
    remove(b6_file)


def cluster_not_matched(curr_fasta, similarity_limit):
    chdir(MAIN_DIR)
    start_of_curr_fasta = curr_fasta.split('.')[0]
    not_matched_file = start_of_curr_fasta + '.not_matched.fasta'
    if getsize(not_matched_file):
        vsearch_clustering(not_matched_file, similarity_limit)
    remove(not_matched_file)


def search_unknown_families_in_known_genera(family_level_fastas):
    chdir(MAIN_DIR)
    for i in tqdm(range(len(family_level_fastas))):
        curr_family = family_level_fastas[i]
        vsearch_blast(curr_family, species_similarity)
        handle_blast_matches_species(curr_family)
        cluster_not_matched(curr_family, species_similarity)
        update_species_stats(curr_family)
        update_species_centroids(curr_family)
    unused_bases = glob('*.base.fasta')
    for i in unused_bases:
        remove(i)


# split list of files for threads
def split_list(alist, wanted_parts=1):
    length = len(alist)
    return [alist[i * length // wanted_parts: (i + 1) * length // wanted_parts]
            for i in range(wanted_parts)]


def threaded_gather_known_genera(all_species_centroids, thread_number, cut_limit):
    for curr_centroid in all_species_centroids:
        curr_centroid_taxonomy = curr_centroid.split('.')[0].split('/')[-1]
        if curr_centroid_taxonomy.count('_') == 6:
            gathered_file_name = MAIN_DIR + '/' + '_'.join(curr_centroid_taxonomy.split('_')[:-cut_limit])
            gathered_file_name += '.base.fasta.thread.' + str(thread_number)
            # system('cat ' + curr_centroid + ' >> ' + gathered_file_name)
            with open(gathered_file_name, 'a') as f:
                with open(curr_centroid, 'r') as f_centroid:
                    f.write(f_centroid.read())


def gather_known_genera(cut_limit):
    chdir(MAIN_DIR + '/species_centroids')
    all_species_centroids = glob('**/*.fasta', recursive=True)
    # TIC_LOG.info(all_species_centroids)
    shuffle(all_species_centroids)
    splitted_species_centroids_list = split_list(all_species_centroids, int(THREADS))
    processes = [multiprocessing.Process(target=threaded_gather_known_genera,
                 args=(splitted_species_centroids_list[x], x, cut_limit)) for x in range(int(THREADS))]
    # Run processes
    for p in processes:
        p.start()
    # Exit the completed processes
    for p in processes:
        p.join()
    chdir(MAIN_DIR)
    threaded_fastas = glob('*.base.fasta.thread.*')
    for threaded_fasta in threaded_fastas:
        final_name = threaded_fasta.split('.thread')[0]
        # system('cat ' + threaded_fasta + ' >> ' + final_name)
        with open(final_name, 'a') as f:
            with open(threaded_fasta, 'r') as f_threaded:
                f.write(f_threaded.read())

        remove(threaded_fasta)


def threaded_gather_unknown_genera(all_species_centroids, thread_number, underscores, cut_limit):
    for curr_centroid in all_species_centroids:
        curr_centroid_taxonomy = curr_centroid.split('.')[0].split('/')[-1]
        if curr_centroid_taxonomy.count('_') == underscores:
            gathered_file_name = '_'.join(curr_centroid_taxonomy.split('_')[:-cut_limit]) + '.genera_queries.fasta'
            gathered_file_name += '.thread.' + str(thread_number)
            # system('cat ' + curr_centroid + ' >> ' + MAIN_DIR + '/' + gathered_file_name)
            with open(MAIN_DIR + '/' + gathered_file_name, 'a') as f:
                with open(curr_centroid, 'r') as f_centroid:
                    f.write(f_centroid.read())


def gather_unknown_genera(underscores, cut_limit):
    chdir(MAIN_DIR + '/species_centroids')
    all_species_centroids = glob('**/*.fasta', recursive=True)
    shuffle(all_species_centroids)
    splitted_species_centroids_list = split_list(all_species_centroids, int(THREADS))
    processes = [multiprocessing.Process(target=threaded_gather_unknown_genera,
                 args=(splitted_species_centroids_list[x], x, underscores, cut_limit)) for x in range(int(THREADS))]
    # Run processes
    for p in processes:
        p.start()
    # Exit the completed processes
    for p in processes:
        p.join()
    chdir(MAIN_DIR)
    threaded_fastas = glob('*.genera_queries.fasta.thread.*')
    for threaded_fasta in threaded_fastas:
        final_name = threaded_fasta.split('.thread')[0]
        # system('cat ' + threaded_fasta + ' >> ' + final_name)
        with open(final_name, 'a') as f:
            with open(threaded_fasta, 'r') as f_threaded:
                f.write(f_threaded.read())
        remove(threaded_fasta)


def handle_blast_matches_genera(curr_fasta):
    chdir(MAIN_DIR)
    start_of_curr_fasta = curr_fasta.split('.')[0]
    b6_file = start_of_curr_fasta + '.b6'
    if not isfile(b6_file):
        return
    elif not getsize(b6_file):
        pass
    # in order not to repeat the remove command at the end of function
    # and not run anything else but clean the files at the end
    else:
        b6_contents = read_file(b6_file)
        for line in b6_contents:
            # TIC_LOG.info(line)
            query = line.split('\t')[0]
            target = line.split('\t')[1]
            target_taxonomy = target.split('tax=')[1]
            found_genus = target_taxonomy.split(';')[5]
            initial_taxonomy = query.split('tax=')[1]
            initial_file_location_stats = MAIN_DIR + '/species_stats/' + '/'.join(initial_taxonomy.split(';')[:-2]) + '/'
            initial_stats_file = initial_file_location_stats + initial_taxonomy.replace(';', '_')[:-1] + '.stats'
            # TIC_LOG.info(initial_stats_file)
            new_taxonomy = ';'.join(initial_taxonomy.split(';')[:-2]) + ';'
            new_taxonomy += found_genus + ';' + initial_taxonomy.split(';')[-2]
            new_stats_file_location_stats = MAIN_DIR + '/species_stats/' + '/'.join(new_taxonomy.split(';')[:-1]) + '/'
            new_stats_file_name = new_stats_file_location_stats + new_taxonomy.replace(';', '_') + '.stats'
            # TIC_LOG.info(new_stats_file_name)
            if not isfile(initial_stats_file):
                continue
            old_stats_contents = read_file(initial_stats_file)
            new_stats_file = open(new_stats_file_name, 'w+')
            for line in old_stats_contents:
                line_tokens = line.split('\t')
                sequence_identifier = line_tokens[1].split('tax=')[0]
                new_stats_file.write(line_tokens[0] + '\t' + sequence_identifier + 'tax=' + new_taxonomy + ';\n')
            new_stats_file.close()
            remove(initial_stats_file)
            old_centroids_file_name = MAIN_DIR + '/species_centroids/' + '/'.join(initial_taxonomy.split(';')[:-2]) + '/'
            old_centroids_file_name += initial_taxonomy.replace(';', '_')[:-1] + '.fasta'
            old_centroids_seqs_dict = fasta2dict(old_centroids_file_name)
            new_centroids_file_name = MAIN_DIR + '/species_centroids/' + '/'.join(new_taxonomy.split(';')[:-1]) + '/'
            new_centroids_file_name += new_taxonomy.replace(';', '_') + '.fasta'
            # TIC_LOG.info(old_centroids_file_name)
            # TIC_LOG.info(new_centroids_file_name)
            new_centroids_file = open(new_centroids_file_name, 'w+')
            query_sequence_identifier = query.split('tax=')[0]
            new_centroids_file.write('>' + query_sequence_identifier + 'tax=' + new_taxonomy + ';\n')
            new_centroids_file.write(old_centroids_seqs_dict[query] + '\n')
            new_centroids_file.close()
            remove(old_centroids_file_name)
    remove(b6_file)


def update_genera_stats(curr_fasta):
    chdir(MAIN_DIR)
    global genera_number
    global genera_dict
    curr_name = curr_fasta.split('.')[0]
    uc_file = curr_name + '.uc'
    if not isfile(uc_file) or not(getsize(uc_file)):
        return
    uc_content = read_file(uc_file)
    for line in uc_content:
        line_tokens = line.split('\t')
        if (line[0] == 'S'):
            genera_number += 1
            curr_species_number = line_tokens[8].split(';')[-2][4:]
            full_curr_name = curr_name + '_GOTU' + str(genera_number)
            curr_dir = create_directory(full_curr_name, 'species_stats')
            new_stats_file_name = curr_dir + curr_name + '_GOTU' + str(genera_number)
            new_stats_file_name += '_SOTU' + curr_species_number + '.stats'
            new_stats_file = open(MAIN_DIR + '/' + new_stats_file_name, 'w+')
            initial_dir = curr_name.replace('_', '/') + '/'
            initial_stats_file = MAIN_DIR + '/species_stats/' + initial_dir + curr_name
            initial_stats_file += '_SOTU' + curr_species_number + '.stats'
            initial_stats_contents = read_file(initial_stats_file)
            for stats_line in initial_stats_contents:
                stats_line_tokens = stats_line.split('\t')
                sequence_identifier = stats_line_tokens[1].split('tax=')[0]
                old_taxonomy = stats_line_tokens[1].split('tax=')[1]
                new_taxonomy = ';'.join(old_taxonomy.split(';')[:-2]) + ';GOTU' + str(genera_number)
                new_taxonomy += ';' + old_taxonomy.split(';')[-2] + ';\n'
                new_stats_file.write(stats_line_tokens[0] + '\t' + sequence_identifier + 'tax=' + new_taxonomy)
            new_stats_file.close()
            genera_dict[line_tokens[8]] = genera_number
            remove(initial_stats_file)
        elif(line[0] == 'H'):
            target_hit_genera_num = str(genera_dict[line_tokens[9]])
            query_seq_identifier = line_tokens[8].split('tax=')[0]
            curr_species_number = line_tokens[8].split(';')[-2][4:]
            old_taxonomy = line_tokens[8].split('tax=')[1]
            new_taxonomy = ';'.join(old_taxonomy.split(';')[:-2]) + ';GOTU' + str(target_hit_genera_num)
            new_taxonomy += ';' + old_taxonomy.split(';')[-2] + ';'
            out_file_name = curr_name + '_GOTU' + str(target_hit_genera_num) + '_SOTU' + curr_species_number + '.stats'
            new_taxonomy_path = '_'.join(out_file_name.split('_')[:-1]).replace('_', '/')
            out_file = open(MAIN_DIR + '/species_stats/' + new_taxonomy_path + '/' + out_file_name, 'a')
            species_hit_with_GOTU = query_seq_identifier + 'tax=' + new_taxonomy
            out_line = 'S' + '\t' + species_hit_with_GOTU + '\n'
            out_file.write(out_line)
            genera_dict[line_tokens[8]] = target_hit_genera_num
            old_taxonomy_dir = '/'.join((old_taxonomy.split(';')[:-2]))
            initial_stats_file = MAIN_DIR + '/species_stats/' + old_taxonomy_dir + '/'
            initial_stats_file += old_taxonomy.replace(';', '_')[:-1] + '.stats'
            initial_contents = read_file(initial_stats_file)
            for init_line in initial_contents[1:]:
                header = init_line.split('\t')[1].split('tax=')[0]
                palia_taxonomy = init_line.split('\t')[1].split('tax=')[1]
                palio_species = palia_taxonomy.split(';')[-2][4:]
                nea_grammi = header + 'tax=' + ';'.join(palia_taxonomy.split(';')[:-2]) + ';GOTU'
                nea_grammi += str(target_hit_genera_num) + ';SOTU' + palio_species + ';'
                out_file.write('H\t' + nea_grammi + '\n')
            out_file.close()
            remove(initial_stats_file)
    remove(uc_file)


def update_genera_centroids():
    # to genera_dict exei gia ka8e ena species pou den egine matched se poio neo genero anoikei
    # anoi3e neo fasta me onoma pou na periexei mesa to gotu kai grapse tis nees plirofories
    # 8imisou sto header na baleis to gotu
    # to dict exei 51 keys, sta centroids exeis 50
    # giati to 1 einai to hit se agnwsto centroids opote
    # gia ka8e ena apo ta keys tou dictionary, anoi3e to fasta kai kane tin parapanw diadikasia
    # oti exw skeftei san sanity check gia auta to exw perasei
    for key, value in genera_dict.items():
        # TIC_LOG.info(key, value)
        curr_GOTU = value
        old_taxonomy = key.split('tax=')[1]
        # TIC_LOG.info(old_taxonomy)
        sequence_identifier = key.split('tax=')[0]
        old_centroid_dir = ';'.join(old_taxonomy.split(';')[:-2]).replace(';', '/') + '/'
        old_fasta_file_name = MAIN_DIR + '/species_centroids/' + old_centroid_dir
        old_fasta_file_name += old_taxonomy.replace(';', '_')[:-1] + '.fasta'
        # TIC_LOG.info(old_fasta_file_name)
        new_taxonomy = ';'.join(old_taxonomy.split(';')[:-2]) + ';GOTU' + str(curr_GOTU)
        new_taxonomy += ';' + old_taxonomy.split(';')[-2] + ';\n'
        new_fasta_file_name = '_'.join(old_taxonomy.split(';')[:-2]) + '_GOTU' + str(curr_GOTU)
        new_fasta_file_name += '_' + old_taxonomy.split(';')[-2] + '.fasta'
        old_seqs_dict = fasta2dict(old_fasta_file_name)
        new_centroid_location = '_'.join(old_taxonomy.split(';')[:-2]) + '_GOTU' + str(curr_GOTU)
        new_centroid_dir = create_directory(new_centroid_location, 'species_centroids')
        # TIC_LOG.info(new_centroid_dir)
        new_fasta_file = open(MAIN_DIR + '/' + new_centroid_dir + new_fasta_file_name, 'w+')
        # TIC_LOG.info(MAIN_DIR + '/' + new_centroid_dir + new_fasta_file_name)
        new_fasta_file.write('>' + sequence_identifier + 'tax=' + new_taxonomy)
        new_fasta_file.write(old_seqs_dict[key] + '\n')
        new_fasta_file.close()
        remove(old_fasta_file_name)


def clean_bases():
    unused_bases = glob('*.base.fasta')
    for i in unused_bases:
        remove(i)


def clean_centroids():
    used_centroids = glob('*.centroids.fasta')
    for i in used_centroids:
        remove(i)


def search_unknown_genera_in_known_genera():
    chdir(MAIN_DIR)
    all_genera_queries_fastas = glob('*.genera_queries.fasta')
    # TIC_LOG.info(all_genera_queries_fastas)
    for curr_query_genera in all_genera_queries_fastas:
        vsearch_blast(curr_query_genera, genera_similarity)
        handle_blast_matches_genera(curr_query_genera)
        cluster_not_matched(curr_query_genera, genera_similarity)
        update_genera_stats(curr_query_genera)
    update_genera_centroids()
    clean_bases()
    clean_centroids()


def handle_blast_matches_family(curr_fasta):
    chdir(MAIN_DIR)
    start_of_curr_fasta = curr_fasta.split('.')[0]
    b6_file = start_of_curr_fasta + '.b6'
    if not isfile(b6_file):
        return
    elif not getsize(b6_file):
        pass
    # in order not to repeat the remove command at the end of function
    # and not run anything else but clean the files at the end
    else:
        b6_contents = read_file(b6_file)
        for line in b6_contents:
            # TIC_LOG.info(line)
            query = line.split('\t')[0]
            target = line.split('\t')[1]
            target_taxonomy = target.split('tax=')[1]
            found_family = target_taxonomy.split(';')[4]
            found_genus = target_taxonomy.split(';')[5]
            initial_taxonomy = query.split('tax=')[1]
            initial_dir = '/species_stats/' + '/'.join(initial_taxonomy.split(';')[:-2]) + '/'
            initial_stats_file = MAIN_DIR + initial_dir + initial_taxonomy.replace(';', '_')[:-1] + '.stats'
            # TIC_LOG.info(initial_stats_file)
            new_taxonomy = ';'.join(initial_taxonomy.split(';')[:-2]) + ';' + found_family + ';'
            new_taxonomy += found_genus + ';' + initial_taxonomy.split(';')[-2]
            new_dir = '/species_stats/' + '/'.join(new_taxonomy.split(';')[:-1]) + '/'
            new_stats_file_name = MAIN_DIR + new_dir + new_taxonomy.replace(';', '_') + '.stats'
            # TIC_LOG.info(new_stats_file_name)
            if not isfile(initial_stats_file):
                continue
            old_stats_contents = read_file(initial_stats_file)
            new_stats_file = open(new_stats_file_name, 'w+')
            for line in old_stats_contents:
                line_tokens = line.split('\t')
                sequence_identifier = line_tokens[1].split('tax=')[0]
                new_stats_file.write(line_tokens[0] + '\t' + sequence_identifier + 'tax=' + new_taxonomy + ';\n')
            new_stats_file.close()
            remove(initial_stats_file)
            #
            initial_dir = '/species_centroids/' + '/'.join(initial_taxonomy.split(';')[:-2]) + '/'
            old_centroids_file_name = MAIN_DIR + initial_dir
            old_centroids_file_name += initial_taxonomy.replace(';', '_')[:-1] + '.fasta'
            # TIC_LOG.info(old_centroids_file_name)
            old_centroids_seqs_dict = fasta2dict(old_centroids_file_name)
            new_dir = '/species_centroids/' + '/'.join(new_taxonomy.split(';')[:-1]) + '/'
            new_centroids_file_name = MAIN_DIR + new_dir + new_taxonomy.replace(';', '_') + '.fasta'
            # TIC_LOG.info(new_centroids_file_name)
            new_centroids_file = open(new_centroids_file_name, 'w+')
            query_sequence_identifier = query.split('tax=')[0]
            new_centroids_file.write('>' + query_sequence_identifier + 'tax=' + new_taxonomy + ';\n')
            new_centroids_file.write(old_centroids_seqs_dict[query] + '\n')
            new_centroids_file.close()
            remove(old_centroids_file_name)
    remove(b6_file)


def handle_blast_matches_family_with_class(curr_fasta, level):
    chdir(MAIN_DIR)
    start_of_curr_fasta = curr_fasta.split('.')[0]
    b6_file = start_of_curr_fasta + '.b6'
    if not isfile(b6_file):
        return
    elif not getsize(b6_file):
        pass
    # in order not to repeat the remove command at the end of function
    # and not run anything else but clean the files at the end
    else:
        b6_contents = read_file(b6_file)
        for line in b6_contents:
            # TIC_LOG.info(line)
            query = line.split('\t')[0]
            target = line.split('\t')[1]
            target_taxonomy = target.split('tax=')[1]
            initial_taxonomy = query.split('tax=')[1]
            initial_dir = '/'.join(initial_taxonomy.split(';')[:-2]) + '/'
            initial_stats_file = MAIN_DIR + '/species_stats/' + initial_dir
            initial_stats_file += initial_taxonomy.replace(';', '_')[:-1] + '.stats'
            # TIC_LOG.info(initial_stats_file)
            if level == 'genus':
                new_taxonomy = ';'.join(target_taxonomy.split(';')[:-2]) + ';'
                new_taxonomy += initial_taxonomy.split(';')[-2]
            if level == 'family':
                new_taxonomy = ';'.join(target_taxonomy.split(';')[:-3]) + ';'
                new_taxonomy += ';'.join(initial_taxonomy.split(';')[-3:-1])
            new_dir = '/'.join(new_taxonomy.split(';')[:-1]) + '/'
            if not isdir(MAIN_DIR + '/species_stats/' + new_dir):
                mkdir(MAIN_DIR + '/species_stats/' + new_dir)
            new_stats_file_name = MAIN_DIR + '/species_stats/' + new_dir + new_taxonomy.replace(';', '_') + '.stats'
            # TIC_LOG.info(new_stats_file_name)
            if not isfile(initial_stats_file):
                TIC_LOG.warning('Stats fiel could not be found. Continuing...')
                continue
            old_stats_contents = read_file(initial_stats_file)
            new_stats_file = open(new_stats_file_name, 'w+')
            for line in old_stats_contents:
                line_tokens = line.split('\t')
                sequence_identifier = line_tokens[1].split('tax=')[0]
                new_stats_file.write(line_tokens[0] + '\t' + sequence_identifier + 'tax=' + new_taxonomy + ';\n')
            new_stats_file.close()
            remove(initial_stats_file)
            old_centroids_file_name = MAIN_DIR + '/species_centroids/' + initial_dir
            old_centroids_file_name += initial_taxonomy.replace(';', '_')[:-1] + '.fasta'
            old_centroids_seqs_dict = fasta2dict(old_centroids_file_name)
            if not isdir(MAIN_DIR + '/species_centroids/' + new_dir):
                mkdir(MAIN_DIR + '/species_centroids/' + new_dir)
            new_centroids_file_name = MAIN_DIR + '/species_centroids/' + new_dir
            new_centroids_file_name += new_taxonomy.replace(';', '_') + '.fasta'
            # TIC_LOG.info(old_centroids_file_name)
            # TIC_LOG.info(new_centroids_file_name)
            new_centroids_file = open(new_centroids_file_name, 'w+')
            query_sequence_identifier = query.split('tax=')[0]
            new_centroids_file.write('>' + query_sequence_identifier + 'tax=' + new_taxonomy + ';\n')
            new_centroids_file.write(old_centroids_seqs_dict[query] + '\n')
            new_centroids_file.close()
            remove(old_centroids_file_name)
    remove(b6_file)


def handle_blast_matches_order(curr_fasta):
    chdir(MAIN_DIR)
    start_of_curr_fasta = curr_fasta.split('.')[0]
    b6_file = start_of_curr_fasta + '.b6'
    if not isfile(b6_file):
        return
    elif not getsize(b6_file):
        pass
    # in order not to repeat the remove command at the end of function
    # and not run anything else but clean the files at the end
    else:
        b6_contents = read_file(b6_file)
        for line in b6_contents:
            # TIC_LOG.info(line)
            query = line.split('\t')[0]
            target = line.split('\t')[1]
            target_taxonomy = target.split('tax=')[1]
            found_family = target_taxonomy.split(';')[4]
            initial_taxonomy = query.split('tax=')[1]
            initial_dir = '/species_stats/' + '/'.join(initial_taxonomy.split(';')[:-2]) + '/'
            initial_stats_file = MAIN_DIR + initial_dir + initial_taxonomy.replace(';', '_')[:-1] + '.stats'
            # TIC_LOG.info(initial_stats_file)
            new_taxonomy = ';'.join(initial_taxonomy.split(';')[:-3]) + ';' + found_family + ';'
            new_taxonomy += ';'.join(initial_taxonomy.split(';')[-3:])
            new_dir = '/species_stats/' + '/'.join(new_taxonomy.split(';')[:-2]) + '/'
            create_directory('/'.join(new_taxonomy.split(';')[:-2]), 'species_stats')
            # TIC_LOG.info(new_taxonomy)
            new_stats_file_name = MAIN_DIR + new_dir + new_taxonomy[:-1].replace(';', '_') + '.stats'
            # TIC_LOG.info(new_stats_file_name)
            if not isfile(initial_stats_file):
                continue
            old_stats_contents = read_file(initial_stats_file)
            new_stats_file = open(new_stats_file_name, 'w+')
            for line in old_stats_contents:
                line_tokens = line.split('\t')
                sequence_identifier = line_tokens[1].split('tax=')[0]
                new_stats_file.write(line_tokens[0] + '\t' + sequence_identifier + 'tax=' + new_taxonomy + '\n')
            new_stats_file.close()
            remove(initial_stats_file)
            #
            old_centroids_file_name = MAIN_DIR + '/species_centroids/' + '/'.join(initial_taxonomy.split(';')[:-2]) + '/'
            old_centroids_file_name += initial_taxonomy.replace(';', '_')[:-1] + '.fasta'
            old_centroids_seqs_dict = fasta2dict(old_centroids_file_name)
            new_centroids_file_name = MAIN_DIR + '/species_centroids/' + '/'.join(new_taxonomy.split(';')[:-2]) + '/'
            new_centroids_file_name += new_taxonomy[:-1].replace(';', '_') + '.fasta'
            create_directory('/'.join(new_taxonomy.split(';')[:-2]), 'species_centroids')
            # TIC_LOG.info(old_centroids_file_name)
            # TIC_LOG.info(new_centroids_file_name)
            new_centroids_file = open(new_centroids_file_name, 'w+')
            query_sequence_identifier = query.split('tax=')[0]
            new_centroids_file.write('>' + query_sequence_identifier + 'tax=' + new_taxonomy + '\n')
            new_centroids_file.write(old_centroids_seqs_dict[query] + '\n')
            new_centroids_file.close()
            remove(old_centroids_file_name)
    remove(b6_file)


def update_family_stats(curr_fasta, level='none'):
    chdir(MAIN_DIR)
    global family_number
    global family_dict
    curr_name = curr_fasta.split('.')[0]
    uc_file = curr_name + '.uc'
    if not isfile(uc_file) or not(getsize(uc_file)):
        return
    uc_content = read_file(uc_file)
    for line in uc_content:
        line_tokens = line.split('\t')
        if (line[0] == 'S'):
            family_number += 1
            # TIC_LOG.info(line)
            curr_gotu_species_number = '_GOTU' + line_tokens[8].split('GOTU')[1][:-1].replace(';', '_')
            if level == 'class':
                # TIC_LOG.info(curr_name)
                new_stats_file_name = curr_name + '_UNKCLASS_UNKORDER_FOTU' + str(family_number)
                new_stats_file_name += curr_gotu_species_number + '.stats'
            elif level == 'phylum':
                new_stats_file_name = curr_name + '_UNKPHYLUM_UNKCLASS_UNKORDER_FOTU' + str(family_number)
                new_stats_file_name += curr_gotu_species_number + '.stats'
            elif level == 'order':
                new_stats_file_name = curr_name + '_UNKORDER_FOTU' + str(family_number)
                new_stats_file_name += curr_gotu_species_number + '.stats'
            else:
                new_stats_file_name = curr_name + '_FOTU' + str(family_number) + curr_gotu_species_number + '.stats'
            new_taxonomy_no_species = '_'.join(new_stats_file_name.split('_')[:-1])
            new_dir = create_directory(new_taxonomy_no_species, 'species_stats')
            # TIC_LOG.info(MAIN_DIR + '/' + new_dir + new_stats_file_name)
            new_stats_file = open(MAIN_DIR + '/' + new_dir + new_stats_file_name, 'w+')
            initial_dir = '/species_stats/' + '/'.join(curr_name.split('_')) + '/'
            initial_dir += curr_gotu_species_number.split('_')[1]
            initial_stats_file = MAIN_DIR + initial_dir + '/' + curr_name + curr_gotu_species_number + '.stats'
            initial_stats_contents = read_file(initial_stats_file)
            for stats_line in initial_stats_contents:
                stats_line_tokens = stats_line.split('\t')
                sequence_identifier = stats_line_tokens[1].split('tax=')[0]
                # TIC_LOG.info(stats_line)
                old_taxonomy = stats_line_tokens[1].split('tax=')[1]
                if level == 'order':
                    new_taxonomy = ';'.join(old_taxonomy.split(';')[:-3]) + ';UNKORDER;FOTU'
                    new_taxonomy += str(family_number) + ';' + ';'.join(old_taxonomy.split(';')[-3:]) + '\n'
                elif level == 'class':
                    new_taxonomy = ';'.join(old_taxonomy.split(';')[:-3]) + ';UNKCLASS;UNKORDER;FOTU'
                    new_taxonomy += str(family_number) + ';' + ';'.join(old_taxonomy.split(';')[-3:]) + '\n'
                elif level == 'phylum':
                    new_taxonomy = old_taxonomy.split(';')[0] + ';UNKPHYLUM;UNKCLASS;UNKORDER;FOTU'
                    new_taxonomy += str(family_number) + ';' + ';'.join(old_taxonomy.split(';')[-3:]) + '\n'
                else:
                    new_taxonomy = ';'.join(old_taxonomy.split(';')[:-3]) + ';FOTU' + str(family_number)
                    new_taxonomy += ';' + ';'.join(old_taxonomy.split(';')[-3:]) + '\n'
                new_stats_file.write(stats_line_tokens[0] + '\t' + sequence_identifier + 'tax=' + new_taxonomy)
            new_stats_file.close()
            family_dict[line_tokens[8]] = family_number
            remove(initial_stats_file)
        elif(line[0] == 'H'):
            # TIC_LOG.info(line)
            target_hit_family_num = str(family_dict[line_tokens[9]])
            query_seq_identifier = line_tokens[8].split('tax=')[0]
            curr_gotu_species_number = '_GOTU' + line_tokens[8].split('GOTU')[1][:-1].replace(';', '_')
            old_taxonomy = line_tokens[8].split('tax=')[1]  # It has a ";" at the end
            if level == 'order':  # RECHECK
                new_taxonomy = ';'.join(old_taxonomy.split(';')[:-3]) + ';UNKORDER;FOTU'
                new_taxonomy += str(target_hit_family_num) + ';' + ';'.join(old_taxonomy.split(';')[-3:]) + '\n'
                out_file_name = curr_name + '_UNKORDER_FOTU' + str(target_hit_family_num)
                out_file_name += curr_gotu_species_number + '.stats'
            elif level == 'class':
                new_taxonomy = ';'.join(old_taxonomy.split(';')[:-3]) + ';UNKCLASS;UNKORDER;FOTU'
                new_taxonomy += str(target_hit_family_num) + ';' + ';'.join(old_taxonomy.split(';')[-3:]) + '\n'
                out_file_name = curr_name + '_UNKCLASS_UNKORDER_FOTU' + str(target_hit_family_num)
                out_file_name += curr_gotu_species_number + '.stats'
            elif level == 'phylum':  # RECHECK
                new_taxonomy = old_taxonomy.split(";")[0] + ';UNKPHYLUM;UNKCLASS;UNKORDER;FOTU'
                new_taxonomy += str(target_hit_family_num) + ';' + ';'.join(old_taxonomy.split(';')[-3:]) + '\n'
                out_file_name = curr_name + '_UNKPHYLUM_UNKCLASS_UNKORDER_FOTU' + str(target_hit_family_num)
                out_file_name += curr_gotu_species_number + '.stats'
            else:
                new_taxonomy = ';'.join(old_taxonomy.split(';')[:-3]) + ';FOTU' + str(target_hit_family_num)
                new_taxonomy += ';' + ';'.join(old_taxonomy.split(';')[-3:]) + '\n'
                out_file_name = curr_name + '_FOTU' + str(target_hit_family_num) + curr_gotu_species_number + '.stats'
            new_dir = create_directory('/'.join(out_file_name.split('_')[:-1]), 'species_stats')
            out_file = open(MAIN_DIR + '/' + new_dir + out_file_name, 'a')
            species_hit_with_FOTU = query_seq_identifier + 'tax=' + new_taxonomy
            out_line = 'S' + '\t' + species_hit_with_FOTU + '\n'
            out_file.write(out_line)
            family_dict[line_tokens[8]] = target_hit_family_num
            old_dir = '/'.join(old_taxonomy.split(';')[:-2]) + '/'
            initial_stats_file = MAIN_DIR + '/species_stats/' + old_dir + old_taxonomy.replace(';', '_')[:-1] + '.stats'
            initial_contents = read_file(initial_stats_file)
            for init_line in initial_contents[1:]:
                header = init_line.split('\t')[1].split('tax=')[0]
                nea_grammi = header + 'tax=' + new_taxonomy + '\n'
                out_file.write('H\t' + nea_grammi)
            out_file.close()
            remove(initial_stats_file)
    remove(uc_file)


def update_family_centroids(level='none'):
    # TIC_LOG.info(family_dict)
    for key, value in family_dict.items():
        curr_FOTU = value
        old_taxonomy = key.split('tax=')[1]
        sequence_identifier = key.split('tax=')[0]
        old_dir = MAIN_DIR + '/species_centroids/' + '/'.join(old_taxonomy.split(';')[:-2]) + '/'
        # TIC_LOG.info(old_dir)
        old_fasta_file_name = old_dir + old_taxonomy.replace(';', '_')[:-1] + '.fasta'
        # TIC_LOG.info(old_fasta_file_name)
        old_seqs_dict = fasta2dict(old_fasta_file_name)
        if(level == 'phylum'):
            new_taxonomy = old_taxonomy.split(';')[0] + ';UNKPHYLUM;UNKCLASS;UNKORDER;FOTU' + str(curr_FOTU)
            new_taxonomy += ';' + ';'.join(old_taxonomy.split(';')[-3:][:-1]) + ';\n'
            new_fasta_file_name = new_taxonomy.replace(';', '_')[:-2] + '.fasta'
        elif(level == 'class'):
            new_taxonomy = ';'.join(old_taxonomy.split(';')[:-3]) + ';UNKCLASS;UNKORDER;FOTU' + str(curr_FOTU)
            new_taxonomy += ';' + ';'.join(old_taxonomy.split(';')[-3:][:-1]) + ';\n'
            new_fasta_file_name = '_'.join(old_taxonomy.split(';')[:-3]) + '_UNKCLASS_UNKORDER_FOTU' + str(curr_FOTU)
            new_fasta_file_name += '_' + '_'.join(old_taxonomy.split(';')[-3:][:-1]) + '.fasta'
        elif(level == 'order'):
            new_taxonomy = ';'.join(old_taxonomy.split(';')[:-3]) + ';UNKORDER;FOTU' + str(curr_FOTU)
            new_taxonomy += ';' + ';'.join(old_taxonomy.split(';')[-3:][:-1]) + ';\n'
            new_fasta_file_name = '_'.join(old_taxonomy.split(';')[:-3]) + '_UNKORDER_FOTU' + str(curr_FOTU)
            new_fasta_file_name += '_' + '_'.join(old_taxonomy.split(';')[-3:][:-1]) + '.fasta'
        else:
            new_taxonomy = ';'.join(old_taxonomy.split(';')[:-3]) + ';FOTU' + str(curr_FOTU)
            new_taxonomy += ';' + ';'.join(old_taxonomy.split(';')[-3:][:-1]) + ';\n'
            new_fasta_file_name = '_'.join(old_taxonomy.split(';')[:-3]) + '_FOTU' + str(curr_FOTU)
            new_fasta_file_name += '_' + '_'.join(old_taxonomy.split(';')[-3:][:-1]) + '.fasta'
        # TIC_LOG.info(new_taxonomy)
        # TIC_LOG.info(new_fasta_file_name)
        new_fasta_name_for_dir = '_'.join(new_fasta_file_name.split('_')[:-1])
        # TIC_LOG.info(new_fasta_name_for_dir)
        new_dir = create_directory(new_fasta_name_for_dir, 'species_centroids')
        new_fasta_file = open(MAIN_DIR + '/' + new_dir + new_fasta_file_name, 'w+')
        new_fasta_file.write('>' + sequence_identifier + 'tax=' + new_taxonomy)
        new_fasta_file.write(old_seqs_dict[key] + '\n')
        new_fasta_file.close()
        remove(old_fasta_file_name)


def search_unknown_family_in_known_family():
    chdir(MAIN_DIR)
    unk_family_fastas = get_curr_level_fastas(3)
    for curr_unk_family in unk_family_fastas:
        if 'base' in curr_unk_family:
            continue
        vsearch_blast(curr_unk_family, species_similarity)
        handle_blast_matches_species(curr_unk_family)
        cluster_not_matched(curr_unk_family, species_similarity)
        update_species_stats(curr_unk_family)
        update_species_centroids(curr_unk_family)
    for i in glob('*.base.fasta'):
        remove(i)
    gather_known_genera(3)
    gather_unknown_genera(4, 1)
    chdir(MAIN_DIR)
    genera_dict.clear()
    all_genera_queries_fastas = glob('*.genera_queries.fasta')
    for curr_query_genera in all_genera_queries_fastas:
        vsearch_blast(curr_query_genera, genera_similarity)
        handle_blast_matches_family(curr_query_genera)
        cluster_not_matched(curr_query_genera, genera_similarity)
        update_genera_stats(curr_query_genera)
    update_genera_centroids()
    clean_bases()
    clean_centroids()
    gather_known_genera(3)
    gather_unknown_genera(5, 2)
    chdir(MAIN_DIR)
    all_family_queries_fastas = glob('*.genera_queries.fasta')
    for curr_query_family in all_family_queries_fastas:
        vsearch_blast(curr_query_family, family_similarity)
        handle_blast_matches_order(curr_query_family)
        cluster_not_matched(curr_query_family, family_similarity)
        update_family_stats(curr_query_family)
        curr_dir = curr_query_family.split('.')[0].replace('_', '/') + '/'
        # system('find ' + MAIN_DIR + '/species_stats/' + curr_dir + ' -maxdepth 1 -type d -empty -delete')
        system_sub(
            [
                'find',
                MAIN_DIR + '/species_stats/' + curr_dir,
                '-maxdepth',
                '1',
                '-type',
                'd',
                '-empty',
                '-delete'
            ]
        )
    update_family_centroids()
    for curr_query_family in all_family_queries_fastas:
        curr_dir = curr_query_family.split('.')[0].replace('_', '/') + '/'
        # system('find ' + MAIN_DIR + '/species_centroids/' + curr_dir + ' -maxdepth 1 -type d -empty -delete')
        system_sub(
            [
                'find',
                MAIN_DIR + '/species_centroids/' + curr_dir,
                '-maxdepth',
                '1',
                '-type',
                'd',
                '-empty',
                '-delete'
            ]
        )
    clean_bases()
    clean_centroids()
    chdir(MAIN_DIR)


def search_unknown_order_in_orders(level_fastas_num, known_genera_num, unknown_genera_num,
                                   unknown_genera_cut_num, unknown_families_num, unknown_families_cut_num, curr_level):
    chdir(MAIN_DIR)
    curr_level_fastas = get_curr_level_fastas(level_fastas_num)
    for curr_class_fasta in curr_level_fastas:
        if 'base' in curr_class_fasta:
            continue
        vsearch_blast(curr_class_fasta, species_similarity)
        handle_blast_matches_species(curr_class_fasta)
        cluster_not_matched(curr_class_fasta, species_similarity)
        update_species_stats(curr_class_fasta)
        update_species_centroids(curr_class_fasta)
    for i in glob('*.base.fasta'):
        remove(i)
    gather_known_genera(known_genera_num)
    gather_unknown_genera(unknown_genera_num, unknown_genera_cut_num)
    chdir(MAIN_DIR)
    genera_dict.clear()
    all_genera_queries_fastas = glob('*.genera_queries.fasta')
    for curr_query_genera in all_genera_queries_fastas:
        vsearch_blast(curr_query_genera, genera_similarity)
        handle_blast_matches_family_with_class(curr_query_genera, 'genus')
        cluster_not_matched(curr_query_genera, genera_similarity)
        update_genera_stats(curr_query_genera)
    update_genera_centroids()
    clean_bases()
    clean_centroids()
    gather_known_genera(known_genera_num)
    gather_unknown_genera(unknown_families_num, unknown_families_cut_num)
    chdir(MAIN_DIR)
    family_dict.clear()
    all_family_queries_fastas = glob('*.genera_queries.fasta')
    # TIC_LOG.info(all_family_queries_fastas)
    for curr_query_family in all_family_queries_fastas:
        vsearch_blast(curr_query_family, family_similarity)
        handle_blast_matches_family_with_class(curr_query_family, 'family')
        cluster_not_matched(curr_query_family, family_similarity)
        update_family_stats(curr_query_family, curr_level)
        curr_dir = curr_query_family.split('.')[0].replace('_', '/') + '/'
        # system('find ' + MAIN_DIR + '/species_stats/' + curr_dir + ' -maxdepth 1 -type d -empty -delete')
        system_sub(
            [
                'find',
                MAIN_DIR + '/species_stats/' + curr_dir,
                '-maxdepth',
                '1',
                '-type',
                'd',
                '-empty',
                '-delete'
            ]
        )
    update_family_centroids(curr_level)
    for curr_query_family in all_family_queries_fastas:
        curr_dir = curr_query_family.split('.')[0].replace('_', '/') + '/'
        # system('find ' + MAIN_DIR + '/species_centroids/' + curr_dir + ' -maxdepth 1 -type d -empty -delete')
        system_sub(
            [
                'find',
                MAIN_DIR + '/species_centroids/' + curr_dir,
                '-maxdepth',
                '1',
                '-type',
                'd',
                '-empty',
                '-delete'
            ]
        )
    clean_bases()
    clean_centroids()
    chdir(MAIN_DIR)


def main_complex_TIC(tool: str, data_dir: str,
                     family_sim: float = 0.90,  # float 0.0 -> 1.0
                     genus_sim: float = 0.95,  # float 0.0 -> 1.0
                     species_sim: float = 0.97,  # float 0.0 -> 1.0
                     threads: int = 10):
    '''
    tool: Should be either usearch or vsearch binary paths
    '''
    if species_sim > float(1) or genus_sim > float(1) or family_sim > float(1):
        TIC_LOG.error('Similarity levels > 100')
        exit
    if species_sim <= float(0) or genus_sim <= float(0) or family_sim <= float(0):
        TIC_LOG.error('Similarity levels <= 0')
        exit()
    if species_sim <= genus_sim:
        TIC_LOG.error('Species Similarity <= Genus Similarity, exiting')
        exit()
    if genus_sim <= family_sim:
        TIC_LOG.error('Genus Similarity <= Family Similarity, exiting')
        exit()
    if species_sim <= family_sim:
        TIC_LOG.error('Species Similarity <= Family Similarity, exiting')
        exit()

    global species_similarity, genera_similarity, family_similarity
    global MAIN_DIR, TOOL, THREADS
    global USEARCH_BIN_CLUST, USEARCH_BIN_BLAST
    species_similarity = str(species_sim)
    genera_similarity = str(genus_sim)
    family_similarity = str(family_sim)
    MAIN_DIR = str(data_dir)
    TOOL = str(tool)
    THREADS = str(threads)
    if 'vsearch' in TOOL:
        USEARCH_BIN_CLUST = TOOL + ' --threads ' + THREADS + ' --cluster_fast '
        USEARCH_BIN_BLAST = TOOL + ' --threads ' + THREADS + ' --usearch_global '
    else:
        USEARCH_BIN_CLUST = TOOL + ' --threads ' + THREADS + ' -cluster_fast '
        USEARCH_BIN_BLAST = TOOL + ' --threads ' + THREADS + ' -usearch_global '

    TIC_LOG.info('Known Genera Clustering starting. Similarity Threshold:' + str(species_similarity))
    try:
        genera_level_fastas = get_curr_level_fastas(5)
    except Exception as exc:
        raise ValueError(f"L953: {str(exc)}")
    try:
        create_species_of_known_genera(genera_level_fastas)
    except Exception as exc:
        raise ValueError(f"L957: {str(exc)}")
    TIC_LOG.info('Known Genera Clustering done')
    TIC_LOG.info('Unknown Genera, Known Family Search starting. Similarity Threshold:' + str(genera_similarity))
    try:
        family_level_fastas = get_curr_level_fastas(4)
    except Exception as exc:
        raise ValueError(f"L963: {str(exc)}")
    # navigate to species_centroids directory and create the targets for the blast
    # the 4 says to which level you must have as the last known
    gather_species(5, '.base.fasta')
    try:
        search_unknown_families_in_known_genera(family_level_fastas)
    except Exception as exc:
        raise ValueError(f"L970: {str(exc)}")
    TIC_LOG.info('Unknown Family Search starting. Similarity Threshold:' + str(family_similarity))
    try:
        gather_known_genera(2)
    except Exception as exc:
        raise ValueError(f"L975: {str(exc)}")
    try:
        gather_unknown_genera(5, 1)
    except Exception as exc:
        raise ValueError(f"L979: {str(exc)}")
    try:
        search_unknown_genera_in_known_genera()
    except Exception as exc:
        raise ValueError(f"L983: {str(exc)}")
    try:
        gather_species(4, '.base.fasta')
    except Exception as exc:
        raise ValueError(f"L987: {str(exc)}")
    try:
        search_unknown_family_in_known_family()
    except Exception as exc:
        raise ValueError(f"L991: {str(exc)}")
    try:
        gather_species(3, '.base.fasta')
    except Exception as exc:
        raise ValueError(f"L995: {str(exc)}")
    try:
        search_unknown_order_in_orders(2, 4, 3, 1, 4, 2, 'order')
    except Exception as exc:
        raise ValueError(f"L999: {str(exc)}")
    try:
        gather_species(2, '.base.fasta')
    except Exception as exc:
        raise ValueError(f"L1003: {str(exc)}")
    try:
        search_unknown_order_in_orders(1, 5, 2, 1, 3, 2, 'class')
    except Exception as exc:
        raise ValueError(f"L1007: {str(exc)}")
    try:
        gather_species(1, '.base.fasta')
    except Exception as exc:
        raise ValueError(f"L1011: {str(exc)}")
    try:
        search_unknown_order_in_orders(0, 6, 1, 1, 2, 2, 'phylum')
    except Exception as exc:
        raise ValueError(f"L1015: {str(exc)}")


if __name__ == '__main__':
    # read arguments for the three levels of similarity
    # as well as the absolute path of the MAIN_DIR
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--family_similarity", help="family similarity 0-100", required=True, type=int)
    parser.add_argument("-g", "--genus_similarity", help="genus similarity 0-100", required=True, type=int)
    parser.add_argument("-s", "--species_similarity", help="species similarity 0-100", required=True, type=int)
    parser.add_argument("-t", "--tool", help="usearch or vsearch", required=True, type=str)
    parser.add_argument("-n", "--threads", help="usearch or vsearch", required=True, type=int)
    parser.add_argument("-d", "--data_dir", required=True, help="Data Directory")
    args = parser.parse_args()

# user input 97 transform it to 0.97
# same for genus and family
    species_similarity = float(args.species_similarity) / float(100)
    genera_similarity = float(args.genus_similarity) / float(100)
    family_similarity = float(args.family_similarity) / float(100)
    MAIN_DIR = str(args.data_dir)
    TOOL = str(args.tool)
    THREADS = int(args.threads)
    main_complex_TIC(TOOL, MAIN_DIR,
                     family_similarity,
                     genera_similarity,
                     species_similarity,
                     THREADS)
