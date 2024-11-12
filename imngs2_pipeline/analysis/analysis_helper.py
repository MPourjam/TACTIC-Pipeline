import os
from pathlib import Path
import tempfile


def onelinefasta(fastafilepath):
    proper_filepath = Path(fastafilepath).resolve()
    dirpath = proper_filepath.parent
    filename = proper_filepath.name
    temp_file = tempfile.NamedTemporaryFile(
        delete=False,
        suffix=str(proper_filepath.suffix),
        dir=str(dirpath)
    )
    newfile_temp_name = temp_file.name
    temp_file.close()
    with proper_filepath.open('r') as f, open(newfile_temp_name, "w+") as onelinefa:
        line = f.readline()
        sequence = ""
        while line:
            if line[0] == '>':
                if sequence:
                    onelinefa.write(sequence + "\n")
                    sequence = ""
                onelinefa.write(line)
            elif line == '\n' or line == '\r\n':
                pass
            else:
                sequence += line.strip()
            line = f.readline()
        if sequence:
            onelinefa.write(sequence + "\n")
    try:
        os.replace(newfile_temp_name, filename)
    finally:
        if Path(newfile_temp_name).exists():
            Path(newfile_temp_name).unlink()
