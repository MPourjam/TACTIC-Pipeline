from ..cutadapt import cutadapt_parameter_sweep
from ..main import main
from pathlib import Path
from pprint import pprint as print


cutadapt_res, stats_res, plot_res = main(
    ref_fastas=[
        "DBs/SILVA/SILVA_138.2_SSURef_NR99_tax_silva_DNA.fasta",
        "DBs/GG2/GG2_backbone_full-length/data/dna-sequences.fasta"
    ],
    primers_tsv=str(Path("pcr_simulation/primers/all_primers.tsv").absolute()),
    n_jobs=22
)
