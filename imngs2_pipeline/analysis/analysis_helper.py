import os
from pathlib import Path
import tempfile


def onelinefasta(fastafilepath):
    proper_filepath = Path(fastafilepath).resolve()
    dirpath = proper_filepath.parent
    filename = proper_filepath.name
    os.chdir(dirpath)
    temp_file = tempfile.NamedTemporaryFile(delete=False, suffix=".fasta")
    newfile_temp_name = temp_file.name
    temp_file.close()
    with proper_filepath.open('r') as f, open(newfile_temp_name, "w+") as onelinefa:
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
    try:
        # shutil.copy2(proper_filepath, proper_filepath.with_suffix(proper_filepath.suffix + ".bak"))  # for troubleshooting
        os.replace(newfile_temp_name, filename)
        # print("{} sequence file is now oneline".format(filename))
    finally:
        if Path(newfile_temp_name).exists():
            Path(newfile_temp_name).unlink()
