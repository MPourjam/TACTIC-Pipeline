#! /usr/bin/env python3
from zipfile import ZipFile
from glob import glob
from pathlib import Path
import argparse
import re
import os


names_f = "/srv/cfm/inputs/20220309-Asp-missing.txt"


def is_good_fastq(file_path):
    fastq_regex = re.compile(".*fastq(\.gz)+$")
    indices_regex = re.compile(r".*(-|_|\.)+I[1-2]{1}(-|_|\.)+.*")
    f_name = Path(file_path).name
    is_fastq = True if fastq_regex.search(f_name) else False
    is_seq = True if not indices_regex.search(f_name) else False
    if is_seq and is_fastq:
        return True
    return False


def latest_files(files_list):
    '''
    Choosing the latest files from files with the same name(not path)
    '''
    if not isinstance(files_list, list):
        raise TypeError("latest_files function gets list of file paths!")
    files = [Path(f) for f in files_list]
    files_dict = {}
    for f in files:
        sub_list_val = files_dict.get(f.name, [])
        sub_list = sub_list_val
        if not isinstance(sub_list_val, list):
            sub_list = [sub_list_val]
        sub_list.append((f, os.stat(f).st_mtime))
        sub_list = sorted(sub_list, key=lambda x: x[1])
        files_dict[f.name] = sub_list[0]
    files = [str(fil[0]) for fil in list(files_dict.values())]
    return files


def zip_them(file_list, dest_file, flat=False):
    dest_file = Path(dest_file)
    dir_path = dest_file.parent
    if not dir_path.is_dir():
        raise TypeError("Destination directory must be an existing directory!!")
    if not isinstance(file_list, list) or not file_list:
        raise TypeError("file_list must be a list and not empty!")
    print("Writing zip to: {}".format(dest_file))
    # Handling duplicated files. Keeping the most recent one
    if flat:
        msg = "Flat option implies choosing tve latest files amongst "
        msg += "files with the same names while writing to the zip file!"
        print(msg)
        file_list = latest_files(file_list)
    with ZipFile(str(dest_file), "w") as dest_zip:
        for f in file_list:
            f_path = Path(f)
            arc_name = f_path.relative_to(dest_file.parent)
            if flat:
                arc_name = f_path.name
            dest_zip.write(f_path, arcname=arc_name)


def match_score(query, subject):
    '''
    It splits the query by -,_, and . and tries to find number of matched segments
    '''
    if not isinstance(query, str) or not isinstance(subject, str):
        raise("Query and match must be strings.")
    delim = re.compile(r"(\.|-|_)")
    # TODO It needs refinement
    q = delim.split(query)  # if len(delim.split(query)) > 1 else query
    s = delim.split(subject)  # if len(delim.split(subject)) > 1 else subject
    q_m_score = 0
    for qu in q:
        # for su in s:
        if qu in s:
            q_m_score += 1
    s_m_score = 0
    if q_m_score == 0:
        for su in s:
            # for qu in q:
            if su in q:
                s_m_score += 1
    return max(s_m_score, q_m_score)


def main(name_file, search_dir, zip_b=False, use_glob=False, use_walk=True, flat_zip=False, use_match_score=False):
    '''
    name_file must have desired search terms per line.
    '''
    search_dir = Path(search_dir)
    name_file = Path(name_file)
    if not search_dir.is_dir():
        raise TypeError("Search directory must be an existing directory!!")
    if not name_file.is_file():
        raise TypeError("{} is not a file.".format(name_file))

    names_pats = []
    with open(name_file, 'r') as names:
        name = str(names.readline()).strip()
        while name:
            names_pats.append(name)
            name = str(names.readline()).strip()

    files = []
    if use_glob:
        g_pat = str(search_dir) + "/**/*{}*"
        for name in names_pats:
            cur_g_pat = g_pat.format(name)
            found_files = glob(cur_g_pat, recursive=True)
            f_files = [Path(f).resolve() for f in found_files if is_good_fastq(f)]
            files.extend(f_files)
    elif use_walk:
        y_files = []
        for dirpath, dirnames, files_c in os.walk(search_dir, topdown=False):
            for fi in files_c:
                fi_path = Path(dirpath).joinpath(fi)
                add_f = False
                if use_match_score:
                    max_match_list = [match_score(n.lower(), fi.lower()) for n in names_pats]
                    max_sim = max(max_match_list)
                    if max_sim > 0:
                        for i, ma in enumerate(max_match_list):
                            add_f = True if ma == max_sim else False
                else:
                    # Splitting the fi into word list to prevent inclusion of small numbers in bigger numbers.
                    word_reg = re.compile(r"[-_\.\\ ]+")
                    fi_word_list = re.split(word_reg, fi.lower())
                    if any([True for n in names_pats if str(n).lower() in fi_word_list]):
                        add_f = True
                if add_f:
                    y_files.append(fi_path)
        for yf in y_files:
            if is_good_fastq(yf):
                files.append(yf)
    # files = list(set(files))
    if zip_b:
        des_fi = name_file.parent.joinpath(name_file.stem + "_files.zip")
        zip_them(files, des_fi, flat=flat_zip)
    else:
        # TODO Returning a file containing all names for which a file was not found would be grear.
        file_paths = name_file.parent.joinpath(name_file.stem + "_paths.txt")
        with open(file_paths, "w+") as out_file:
            for f in files:
                out_file.write(str(f) + '\n')
    return files


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory",
                        help="The directory to find recursively files",
                        default=".")
    parser.add_argument("-f", "--name-file",
                        help="The file contianing names of files one per line.",
                        )
    parser.add_argument("-z", "--zip",
                        help="Boolean to zip them in the directroy of name_file.",
                        default=False, action='store_true')
    parser.add_argument("-flat", "--flat",
                        help="Boolean to zip files without hierarchy.",
                        default=False, action='store_true')
    args = parser.parse_args()
    lookup_dir = Path(args.directory).resolve()
    name_file_path = Path(args.name_file).resolve()
    zip_it = args.zip
    flat_zip = args.flat
    if flat_zip:
        print("Files with duplicated name but different path would get overwritten.")
    main(name_file_path, lookup_dir, zip_b=zip_it, flat_zip=flat_zip)
