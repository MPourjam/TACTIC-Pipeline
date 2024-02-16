#!/usr/bin/env python3
import argparse
import subprocess
import math
import statistics
import pandas as pd
from imngs2_pipeline.processing_helper import gimmelogger
from sys import version_info
if version_info[0] < 3:
    from pathlib2 import Path  # pip2 install pathlib2
else:
    from pathlib import Path


global NORM_LOG
NORM_LOG = gimmelogger(
    logger_name="run_imngs2.analysis.normalization",
    log_file=Path.cwd().joinpath("Analysis_log.txt"),
    only_file=True
)


def normalize_otu_table(otu_table_path: str, spikes_stats_path: str, norm_methods: set = {"Spike", "SampleMinCount"}) -> set:
    valid_norm_methods = {"Spike", "SampleMinCount"}
    # Only normalization methods in valid_norm_methods or numeric values are allowed
    norm_methods = {nm for nm in norm_methods if nm in valid_norm_methods or nm.isnumeric()}
    norm_methods = {str(nm) for nm in norm_methods}
    otu_table_path = Path(otu_table_path).absolute()
    spikes_stats_path = Path(spikes_stats_path).absolute()
    if not otu_table_path.is_file() or not spikes_stats_path.is_file():
        raise FileNotFoundError("One of the given files for normalization does not exist. Please check the paths.")

    samples = {}
    observed_weights = []
    observed_amounts = []
    observed_counts = []
    with open(spikes_stats_path, 'r') as stats_h:
        for line in stats_h:
            if line.startswith("#"):
                assert('#SampleID\tSpikeReads\tspikes_total_weight_in_g\tspike_amount' in line.strip())
            elif line == "\n":
                continue
            else:
                fields = line.strip().split("\t")
                sample_id = fields[0]
                try:
                    spike_reads = int(fields[1])
                except Exception:
                    spike_reads = 0
                try:
                    original_total_weight_in_g = float(fields[2])
                except ValueError:
                    original_total_weight_in_g = float("nan")
                try:
                    amount = float(fields[3])
                except Exception:
                    amount = float("nan")
                if not (math.isclose(amount, 0, abs_tol=1e-5) or math.isnan(amount)):
                    observed_amounts.append(amount)
                if not math.isclose(spike_reads, 0):
                    observed_counts.append(spike_reads)
                if not math.isnan(original_total_weight_in_g):
                    observed_weights.append(original_total_weight_in_g)
                samples[sample_id] = (spike_reads, original_total_weight_in_g, amount)

    sum = 0
    n = 0
    for values in samples.values():
        count, _, amount = values
        if not ((math.isclose(amount, 0, abs_tol=1e-5) or math.isnan(amount)) and math.isclose(count, 0, abs_tol=1e-5)):
            sum += count
            n += 1

    otu_table = pd.read_csv(otu_table_path, delimiter="\t", index_col=0)

    if n == 0 or len(observed_amounts) == 0 or len(observed_counts) == 0:
        # delete "Spike" from norm_methods set if amount and count are not inferralbe from the rest of entries
        NORM_LOG.warning(f"Normalization method 'Spike' is not applicable. No spike reads or amount found in spike stats file.")
        norm_methods.discard("Spike")
        mean = otu_table.iloc[:, :-1].sum().min()
    else:
        mean = (sum / n)
    norm_methods_to_return = set()
    for norm_method in norm_methods:
        otu_table_copy = otu_table.copy()
        try:
            for file_id, values in samples.items():
                count, weight, amount = values
                factor = None
                if norm_method == "Spike":
                    thismean = mean
                    if math.isnan(weight):
                            if len(observed_weights) == 0:
                                weight = 1  # fallback to 1
                            else:
                                weight = statistics.median(observed_weights)  # fallback to median weight
                    # At this point we are sure that either observed_amounts or observed_counts are not empty
                    if math.isclose(amount, 0, abs_tol=1e-5) or math.isnan(amount):
                        amount = statistics.median(observed_amounts)
                    if math.isclose(count, 0, abs_tol=1e-5):
                        count = statistics.mean(observed_counts)
                    molarity_factor = 600 / (amount * 100)
                    factor = thismean / (count * weight * molarity_factor)
                    normalized_otu_path = otu_table_path.parent.joinpath(f"SpikeNormalized-{str(otu_table_path.name)}")
                elif norm_method == "SampleMinCount":
                    thismean = mean
                    factor = thismean / otu_table[file_id].sum()
                elif norm_method.isnumeric():
                    factor = float(norm_method) / otu_table[file_id].sum()
                else:
                    NORM_LOG.warning(f"Normalization method {norm_method} is not implemented")

                normalized_column = otu_table_copy[file_id] * factor
                otu_table_copy[file_id] = normalized_column.round(5)
        except Exception as e:
            NORM_LOG.warning(f"Error while normalizing {file_id} with normalization method '{norm_method} on file {otu_table_path.name}': {e}")
            # continue the outer for loop
            continue
        else:
            normalized_otu_path = otu_table_path.parent.joinpath(f"{norm_method}Normalized-{str(otu_table_path.name)}")
            otu_table_copy.to_csv(normalized_otu_path, sep="\t")
            # add norm_method to return set
            norm_methods_to_return.add(norm_method)
            # correct rights of output folders
            subprocess.call(['chmod', '777', '-R', str(normalized_otu_path)])

    return norm_methods_to_return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("otu_table_path", help="Path of original OTU Table. Output will be written to same folder", type=str)
    parser.add_argument("spikes_stats_path", help="Path to spike stats file.", type=str)
    args = parser.parse_args()

    otu_table_path = args.otu_table_path
    spikes_stats_path = args.spikes_stats_path

    # call the function
    normalized_otu_path = normalize_otu_table(otu_table_path, spikes_stats_path)
