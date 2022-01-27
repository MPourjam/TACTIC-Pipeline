#! /usr/bin/env python3
from zipfile import ZipFile
from glob import glob
from pathlib import Path
import argparse
import re
import os


names_f = "/home/mohsen/CRC1371/Tasks/CFM/Erik_Files/Erik_barIDs.txt"


def is_good_fastq(file_path):
    fastq_regex = re.compile(".*fastq(\.gz)+$")
    indices_regex = re.compile(r".*(-|_|\.)+I[1-2]{1}(-|_|\.)+.*")
    f_name = Path(file_path).name
    is_fastq = True if fastq_regex.search(f_name) else False
    is_seq = True if not indices_regex.search(f_name) else False
    if is_seq and is_fastq:
        return True
    return False


def zip_them(file_list, dest_file, flat=False):
    dest_file = Path(dest_file)
    dir_path = dest_file.parent
    if not dir_path.is_dir():
        raise TypeError("Destination directory must be an existing directory!!")
    if not isinstance(file_list, list) or not file_list:
        raise TypeError("file_list must be a list and not empty!")
    print("Writing zip to: {}".format(dest_file))
    with ZipFile(str(dest_file), "w") as dest_zip:
        for f in file_list:
            f_path = Path(f)
            arc_name = f_path.relative_to(dest_file.parent)
            if flat:
                print("When the flat flag is given, some files with the same name could get overwritten.")
                arc_name = f_path.name
            dest_zip.write(f_path, arcname=arc_name)


def main(name_file, search_dir, zip_b=False, use_glob=False, use_walk=True):
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
        name = names.readline()[:-1]
        while name:
            names_pats.append(name)
            name = names.readline()[:-1]

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
                if any([True for n in names_pats if str(n).lower() in str(fi).lower()]):
                    fi_path = os.path.join(dirpath, fi)
                    y_files.append(fi_path)

        for yf in y_files:
            if is_good_fastq(yf):
                files.append(str(yf))

    if zip_b:
        des_fi = search_dir.parent.joinpath(name_file.stem + "_files.zip")
        zip_them(files, des_fi)
    else:
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
    args = parser.parse_args()
    lookup_dir = Path(args.directory).resolve()
    name_file_path = Path(args.name_file).resolve()
    zip_it = args.zip
    main(name_file_path, lookup_dir, zip_b=zip_it)
