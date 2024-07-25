import argparse
from tqdm import tqdm
from os import mkdir
from os import path as ospath
from processing_helper import gimmelogger


global TIC_LOG
TIC_LOG = gimmelogger(
    logger_name="run_imngs2.analysis.TIC",
)


def read_file(filename):
    with open(filename) as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    return content


def check_valid_names(input_taxonomy):
    taxo_tokens = input_taxonomy.split(';')
    valid_names_list = list()
    phylum = ''
    class_s = ''
    order_o = ''
    family = ''
    genus = ''
    if len(taxo_tokens) >= 1:
        try:
            domain = taxo_tokens[0].split()[0]
        except BaseException:
            domain = ''
    if len(taxo_tokens) >= 2:
        try:
            phylum = "-".join(taxo_tokens[1].split())
        except BaseException:
            phylum = ''
    if len(taxo_tokens) >= 3:
        try:
            class_s = "-".join(taxo_tokens[2].split())
        except BaseException:
            class_s = ''
    if len(taxo_tokens) >= 4:
        try:
            order_o = "-".join(taxo_tokens[3].split())
        except BaseException:
            order_o = ''
    if len(taxo_tokens) >= 5:
        try:
            family = "-".join(taxo_tokens[4].split())
        except BaseException:
            family = ''
    if len(taxo_tokens) >= 6:
        try:
            genus = "-".join(taxo_tokens[5].split())
        except BaseException:
            genus = ''
    if "archaea" in domain.lower() or "bacteria" in domain.lower():
        valid_names_list.append(domain)
    else:
        return('error')
    # The list defines the keywords by which we drop the taxa at each level.
    # Some are deliberately left uncomplete to cover more possible combination
    no_goes = [
               "unclassi",
               "uncultur",
               "unknown",
               "incertae",
               "unkdom",
               "unkking",
               "unkphyl",
               "unkclass",
               "unkorder",
               "unkfam",
               "unkgen",
               "unkspec",
               ]
    if phylum and not any([True for ng in no_goes if ng in phylum.lower()]):
        valid_names_list.append(phylum)
    else:
        return('_'.join(valid_names_list) + '.fasta')
    if class_s and not any([True for ng in no_goes if ng in class_s.lower()]):
        valid_names_list.append(class_s)
    else:
        return('_'.join(valid_names_list) + '.fasta')
    if order_o and not any([True for ng in no_goes if ng in order_o.lower()]):
        valid_names_list.append(order_o)
    else:
        return('_'.join(valid_names_list) + '.fasta')
    if family and not any([True for ng in no_goes if ng in family.lower()]):
        valid_names_list.append(family)
    else:
        return('_'.join(valid_names_list) + '.fasta')
    if genus and not any([True for ng in no_goes if ng in genus.lower()]):
        valid_names_list.append(genus)
    return('_'.join(valid_names_list) + '.fasta')


def split_based_on_taxonomy(data_dir: str, input_taxed_zotu: str):
    input_file = ospath.abspath(input_taxed_zotu)
    MAIN_DIR = ospath.abspath(data_dir)
    try:
        mkdir(MAIN_DIR)
    except BaseException:
        TIC_LOG.error(MAIN_DIR + ' already present please select other directory or move the present directory to another location')
        exit()
    # Correcting taxonomy for special characters
    tax_correct_tr_dict = {
        ord("("): '',
        ord(")"): '',
        ord("["): '',
        ord("]"): '',
        ord("_"): ord("-"),
        ord("."): ord("-"),
        ord(" "): '',
    }
    input_contents = read_file(input_file)
    TIC_LOG.info(f"Splitting the input file based on taxonomy")
    for i in tqdm(range(0, len(input_contents), 2)):
        header = input_contents[i]
        sequence = input_contents[i+1]
        curr_taxonomy = str(header.split('tax=')[1]).translate(tax_correct_tr_dict)
        out_file_name = check_valid_names(curr_taxonomy)
        if out_file_name == 'error':
            continue
        tokens = out_file_name.split('_')
        counter = 1
        new_out_file_name_tokens = list()
        for curr_token in tokens:
            if counter < 7:
                new_out_file_name_tokens.append(curr_token)
            counter += 1
        new_out_file_name = '_'.join(new_out_file_name_tokens)
        if '.fasta' not in new_out_file_name:
            new_out_file_name += '.fasta'
        out_file_path = ospath.join(MAIN_DIR, new_out_file_name)
        out_file = open(out_file_path, 'a+')
        new_header = header.split('tax=')[0] + 'tax=' + new_out_file_name.split('.')[0].replace('_', ';') + ';'
        out_file.write(new_header + '\n')
        out_file.write(sequence + '\n')
        out_file.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--data_dir", required=True, help="Data Directory")
    parser.add_argument("-i", "--input", required=True, help="Input File")
    args = parser.parse_args()

    MAIN_DIR = args.data_dir + '/'
    input_file = args.input
    split_based_on_taxonomy(MAIN_DIR, input_file)
