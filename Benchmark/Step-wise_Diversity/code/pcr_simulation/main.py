#!/usr/bin/env python3
"""
This module gets an input fasta file of sequences and runs cutadapt with desired arguments on it.
It then generate statistics and plots on the trimmed output files.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
from .cutadapt import cutadapt_parameter_sweep
from .fasta_stats import *
from .helpers import *


RUN_PARAMETERS = {
    "ref_fastas": [
        "DBs/SILVA/SILVA_138.2_SSURef_NR99_tax_silva_DNA.fasta",
        "DBs/GG2/GG2_backbone_full-length/data/dna-sequences.fasta"
    ],
    "output_dir": "./results",
    "plot": True,
    "n_jobs": 1,
    "cutadapt_threads": 5,
    "error_rates": [
        0.1,
        0.15,
        # 0.2
    ],
    "min_overlaps": [
        3,
        6,
        9,
        12,
        15
    ],
    "min_amp_lens": [
        # 100,
        150,
        200
    ],
    "max_amp_lens": [
        # 1000,
        1500,
        2000
    ],
    "primers_tsv": str(Path("pcr_simulation/primers/Isi_primers.tsv").absolute())
}

def main(
        ref_fastas: list = [],
        output_dir: str = './results',
        primers_tsv: str = '',
        n_jobs: int = 1):
    """
    It runs parameter sweep on a range of parameters and 
    generates stats and plot for each output amplicon fasta.
    """
    if len(ref_fastas) == 0:
        ref_fastas = RUN_PARAMETERS["ref_fastas"]
    if primers_tsv == '':
        primers_tsv = RUN_PARAMETERS["primers_tsv"]
    print(f"Using primers TSV: {primers_tsv}")
    # exit()
    sweep_results = cutadapt_parameter_sweep(
        ref_fastas=ref_fastas,
        out_dir=output_dir,
        n_jobs=n_jobs,
        cutadapt_threads=RUN_PARAMETERS["cutadapt_threads"],
        error_rates=RUN_PARAMETERS["error_rates"],
        min_overlaps=RUN_PARAMETERS["min_overlaps"],
        min_amp_lens=RUN_PARAMETERS["min_amp_lens"],
        max_amp_lens=RUN_PARAMETERS["max_amp_lens"],
        primers_tsv=primers_tsv
    )


    # write arguments and statistics of trimmed fasta to a file next to it
    stat_results = []
    tsv_writer_results = []
    plot_len_dist_results = []
    with ProcessPoolExecutor(max_workers=n_jobs) as executor:
        futures = {}
        for each_run in sweep_results:
            run_id, params_d, primers_names, results_paths = each_run
            trimmed_fasta, cutadapt_json, untrimmed_fasta = results_paths
            for this_fasta in (trimmed_fasta, untrimmed_fasta):
                future = executor.submit(fasta_stats, this_fasta)
                futures[future] = (run_id, this_fasta, params_d, primers_names)

        for future in as_completed(futures):
            run_id, fasta_path, params_d, primers_names = futures[future]
            try:
                stats = future.result()
                params_d_up = params_d.copy()
                params_d_up.update(
                    {
                        "fwd_name": primers_names[0],
                        "rev_name": primers_names[1]
                    }
                )
                params_d_up.update(stats)
                stat_results.append((run_id, fasta_path, params_d_up))
            except Exception as e:
                print(f"Error processing {fasta_path} in run {run_id}: {e}")
        
        # write the stats to files next to the
        tsv_writer_futures = []
        for res_ in stat_results:
            run_id, fasta_path, params_d_up = res_
            fasta_path = Path(fasta_path).with_suffix('.stats.tsv')
            future = executor.submit(DictWriter.write_tsv, params_d_up, str(fasta_path))
            tsv_writer_futures.append(future)

        for wrtier_future in as_completed(tsv_writer_futures):
            try:
                tsv_writer_results.append(wrtier_future.result())
            except Exception as e:
                print(f"Error writing TSV file: {e}")
                tsv_writer_results.append(None)

        # Write the plots next to each
        plot_len_dist_futures = []
        for res_ in stat_results:
            run_id, fasta_path, params_d_up = res_
            future = executor.submit(plot_fasta_length_histogram, fasta_path)
            plot_len_dist_futures.append(future)

        for plot_future in as_completed(plot_len_dist_futures):
            try:
                plot_len_dist_results.append(plot_future.result())
            except Exception as e:
                print(f"Error writing plot file: {e}")
                plot_len_dist_results.append(None)

    return stat_results, tsv_writer_results, plot_len_dist_results

