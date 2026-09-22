import os
import re
from pathlib import Path, PurePath
import tempfile


def onelinefasta(fastafilepath):
    proper_filepath = Path(fastafilepath).resolve()
    dirpath = proper_filepath.parent
    filename = proper_filepath.name
    with tempfile.NamedTemporaryFile(
            delete=False,
            suffix=str(proper_filepath.suffix),
            dir=str(dirpath)) as temp_file:
        newfile_temp_name = temp_file.name
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
        os.replace(newfile_temp_name, proper_filepath)
    finally:
        if Path(newfile_temp_name).exists():
            Path(newfile_temp_name).unlink()


def keep_seq_id_only(input_fasta: str, output_fasta: str) -> None:
    """
    Removes taxonomy information from fasta headers.
    :param input_fasta: Path to the input fasta file.
    :param output_fasta: Path to the output fasta file without taxonomy info.
    """
    tax_regex = re.compile(r"^>(?P<seq_id>[^\s]+)", re.IGNORECASE)

    with open(input_fasta, 'r') as infile, open(output_fasta, 'w+') as outfile:
        for line in infile:
            if line.startswith('>'):
                cleaned_line = tax_regex.match(line).group(0)
                outfile.write(cleaned_line + '\n')
            else:
                outfile.write(line)


def remove_taxonomy_from_table(input_table: str, output_table: str, sep: str = '\t') -> None:
    """
    Removes taxonomy information from the last column of a tab-separated table.
    :param input_table: Path to the input table file.
    :param output_table: Path to the output table file without taxonomy info.
    """
    with open(input_table, 'r') as infile, open(output_table, 'w+') as outfile:
        header = infile.readline().strip().lower()
        header_l = header.split(sep)
        tax_index = [i for i, col in enumerate(header_l) if col == 'taxonomy']
        if tax_index:
            # Taking everything except the column with taxonomy
            header = sep.join([col for i, col in enumerate(header_l) if i != tax_index[0]])
            header = header.strip()
        outfile.write(header + '\n')
        for line in infile:
            line = line.strip()
            line_l = line.split(sep)
            if tax_index:
                line = sep.join([col for i, col in enumerate(line_l) if i != tax_index[0]])
                line = line.strip()
            outfile.write(line + '\n')


def prepare_picrust2_input(analysis_dir: str) -> None:
    """
    It finds all OTU tables and sequnces fasta files and deletes taxonomy
    from them and only keeps the sequence id in the fastas.
    Any file following the naming pattern '.*[ZS]?OTUs?-Table.*.tab' will be considered as an OTU table
    Any file following the naming pattern '.*[ZS]?OTUs?-Seqs.*.fasta' will be considered as a sequences fasta file
    """
    analysis_dir = Path(PurePath(analysis_dir)).absolute()
    otu_table_pattern = re.compile(r"(?P<prefix>.*[ZS]?OTUs?-Table.*)\.tab$", re.IGNORECASE)
    otu_fasta_pattern = re.compile(r"(?P<prefix>.*[ZS]?OTUs?-Seqs.*)\.fasta$", re.IGNORECASE)
    picrust2_inputs_dir = analysis_dir / "PICRUSt2-Inputs"
    picrust2_inputs_dir.mkdir(exist_ok=True)

    for item in analysis_dir.rglob("*"):
        if item.is_file():
            if picrust2_inputs_dir in item.parents:
                continue
            match_table = otu_table_pattern.match(item.name)
            match_fasta = otu_fasta_pattern.match(item.name)

            if match_table:
                prefix = match_table.group("prefix")
                otu_table_file = item
                rel_parent = item.parent.relative_to(analysis_dir)
                rel_parent_str = str(rel_parent).replace(os.sep, "__")
                rel_parent_str = rel_parent_str.strip(".")
                if rel_parent_str:
                    prefix = f"{rel_parent_str}__{prefix}"
                picrust2_otu_table_file = picrust2_inputs_dir / f"{prefix}-For-PICRUSt2.tab"
                remove_taxonomy_from_table(otu_table_file, picrust2_otu_table_file)

            elif match_fasta:
                prefix = match_fasta.group("prefix")
                otu_fasta_file = item
                rel_parent = item.parent.relative_to(analysis_dir)
                rel_parent_str = str(rel_parent).replace(os.sep, "__")
                rel_parent_str = rel_parent_str.strip(".")
                if rel_parent_str:
                    prefix = f"{rel_parent_str}__{prefix}"
                picrust2_otu_fasta_file = picrust2_inputs_dir / f"{prefix}-For-PICRUSt2.fasta"
                keep_seq_id_only(otu_fasta_file, picrust2_otu_fasta_file)
                onelinefasta(picrust2_otu_fasta_file)
