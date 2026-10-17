"""Preserve the species-label definition used by the older SILVA/LTP summaries."""
from pathlib import Path


def species_name(header):
    parts = header.split(maxsplit=1)
    levels = [part.strip() for part in parts[1].split(";") if part.strip()] if len(parts) == 2 else []
    if len(levels) != 7:
        return ""
    words = [word.strip(".") for word in levels[-1].split() if word.strip()]
    return " ".join(words[:2]).lower() if len(words) >= 2 else ""


def configure(workflow, profile):
    if not hasattr(workflow, "_bundle_original_filter"):
        workflow._bundle_original_filter = workflow._filter_fasta_by_valid_species
        workflow._bundle_original_counter = workflow._count_unique_taxonomies
    workflow._filter_fasta_by_valid_species = workflow._bundle_original_filter
    workflow._count_unique_taxonomies = workflow._bundle_original_counter
    if profile == "strict_full_taxonomy":
        return
    if profile != "legacy_species_name":
        raise ValueError("Unknown taxonomy profile: " + profile)

    def filter_preserving_headers(input_fasta, output_fasta, reuse, discarded_file=None):
        output_fasta = Path(output_fasta)
        output_fasta.parent.mkdir(parents=True, exist_ok=True)
        if reuse and output_fasta.exists():
            return output_fasta
        rejected = []
        kept = 0
        with output_fasta.open("w") as output:
            for header, sequence in workflow._iter_fasta_records(Path(input_fasta)):
                valid, reason = workflow._validate_species_taxonomy(header)
                if valid:
                    output.write(f">{header}\n{sequence}\n")
                    kept += 1
                else:
                    rejected.append((header, reason))
        if discarded_file:
            with Path(discarded_file).open("w") as output:
                output.write("header\treason\n")
                for header, reason in rejected:
                    output.write(header + "\t" + reason + "\n")
        if not kept:
            raise RuntimeError("No valid species headers in " + str(input_fasta))
        return output_fasta

    def count_species_names(path):
        return len({species_name(header) for header, _ in workflow._iter_fasta_records(Path(path))
                    if workflow._validate_species_taxonomy(header)[0]})

    workflow._filter_fasta_by_valid_species = filter_preserving_headers
    workflow._count_unique_taxonomies = count_species_names
