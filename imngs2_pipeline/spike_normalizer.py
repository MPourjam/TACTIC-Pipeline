#!/usr/bin/env python3
import argparse
import subprocess
import math
import statistics
import pandas as pd
from enum import Enum
# from imngs2_pipeline.processing_helper import gimmelogger
from .processing_helper import gimmelogger

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


class SpikeNormValidity(str, Enum):
    NS = "NotSpiked"
    NA = "NotApplicable"
    PC = "PossibleContamination"
    PS = "PossiblySpiked"


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
    sum = 0
    n = 0
    with open(spikes_stats_path, 'r') as stats_h:
        for line in stats_h:
            if line.startswith("#"):
                assert('#SampleID\tSpikeReads\tspikes_total_weight_in_g\tspike_amount' in line.strip())
            elif line == "\n":
                continue
            else:
                fields = line.strip().split("\t")
                sample_id = fields[0]
                # Inferring spike reads, weight and amount from the spike stats file
                try:
                    spike_reads = str(fields[1]).replace(",", ".").split(".")[0]
                    spike_reads = int(spike_reads)
                except Exception:
                    NORM_LOG.warning(f"Could not parse spike reads for sample {sample_id}. Setting to 0.")
                    spike_reads = 0
                try:
                    original_total_weight_in_g = float(fields[2])
                except ValueError:
                    original_total_weight_in_g = float("nan")
                try:
                    amount = float(fields[3])
                except Exception:
                    amount = float("nan")
                # Categorizing entries based on the validity of the spike normalization
                if spike_reads == 0 and amount > 0.0:
                    # This is a possible mistake so that normalization with spike is not applicable
                    spike_nrom_validity = SpikeNormValidity.NA
                elif spike_reads == 0 and (math.isclose(amount, 0, abs_tol=1e-5) or math.isnan(amount)):
                    # The sample is not spiked
                    spike_nrom_validity = SpikeNormValidity.NS
                elif spike_reads > 0 and math.isclose(amount, 0, abs_tol=1e-5):
                    # The sample is not spiked but has spiked reads count
                    # This is a possible contamination so that we force spike_reads to be 0
                    # Then SpikeNormValidity.PC is equal to SpikeNormValidity.NS
                    spike_nrom_validity = SpikeNormValidity.PC
                    spike_reads = 0
                else:
                    # The sample is spiked but might have missed values for weight or amount
                    spike_nrom_validity = SpikeNormValidity.PS
                # Collecting observed valid values for spike normalization
                if spike_reads > 0:
                    sum += spike_reads
                    n += 1
                    observed_counts.append(spike_reads)
                if not (math.isclose(amount, 0, abs_tol=1e-5) or math.isnan(amount)):
                    observed_amounts.append(amount)
                if not math.isnan(original_total_weight_in_g):
                    observed_weights.append(original_total_weight_in_g)
                samples[sample_id] = (spike_reads, original_total_weight_in_g, amount, spike_nrom_validity)

    otu_table = pd.read_csv(otu_table_path, delimiter="\t", index_col=0)

    spike_norm_invalid_samples = [sample_id for sample_id, values in samples.items() if values[3] not in {SpikeNormValidity.PS}]
    all_samples_not_spiked = all([math.isnan(values[2]) or math.isclose(values[2], 0, abs_tol=1e-5) for values in samples.values()])

    if all_samples_not_spiked:
        # delete "Spike" from norm_methods set if all samples are not spiked
        NORM_LOG.info(
            f"Normalization method 'Spike' is not applicable."
            " No spiked sample was found.")
        norm_methods.discard("Spike")

    elif len(observed_amounts) == 0 or len(observed_counts) == 0 or any(spike_norm_invalid_samples):
        # delete "Spike" from norm_methods set if amount and count are not inferralbe from the rest of entries
        NORM_LOG.info(
            f"Normalization method 'Spike' is not applicable."
            " Some of entires in spike stats file are not valid"
            " for spike normalization. Possible mixture of spiked and non-spiked samples!")
        norm_methods.discard("Spike")
    # Raise index error if not all otu_tables columns are in samples
    # Do not include the first and last column in the check
    if not all([col in samples for col in otu_table.columns[1:-1]]):
        raise IndexError(f"Columns in OTU table {otu_table_path.name} are not in spike stats file {spikes_stats_path.name}")
    factor = 1
    norm_methods_to_return = set()
    for norm_method in norm_methods:
        otu_table_copy = otu_table.copy()
        try:
            for file_id, values in samples.items():
                count, weight, amount, spike_norm_validity = values
                factor = None
                if norm_method == "Spike":
                    thismean = sum / n  # mean of spike reads
                    # Interpolation of missing values for weight and amount
                    if math.isnan(amount):
                        amount = statistics.median(observed_amounts)
                    if math.isnan(weight) or math.isclose(weight, 0, abs_tol=1e-5):
                            if len(observed_weights) == 0:
                                weight = 1  # fallback to 1
                            else:
                                weight = statistics.median(observed_weights)  # fallback to median weight
                    molarity_factor = 600 / (amount * 100)
                    spike_count_factor = thismean / count
                    # NOTE:TODO Having molarity factor and taking 600 Microliter of all samples
                    # in the lab make normalization to weight meaningless?
                    # factor = 600 / (amount * 100)
                    # otu_table[file_id] = (otu_table[file_id] * thismean) / (count * weight * factor)
                    factor = spike_count_factor / (weight * molarity_factor)
                elif norm_method == "SampleMinCount":
                    thismean = otu_table.iloc[:, :-1].sum().min()
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
