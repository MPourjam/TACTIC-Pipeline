#! /bin/env bash
# set -x  # Uncomment for debugging
set -euo pipefail


SILVA_URL="https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_SSURef_tax_silva.fasta.gz"
OUTPUT_FILE="SILVA_138.2_SSURef_tax_silva.fasta.gz"


/bin/env python3 stepwise_16S_diversity.py \
    download \
    --db silva \
    --url "$SILVA_URL" \
    -o "$OUTPUT_FILE" \
    --decompress


/bin/env python3 stepwise_16S_diversity.py simulate \
    --ref SILVA_138.2_SSURef_tax_silva.fasta \
    --primers primers.tsv \
    --forward 27F \
    --reverse-targets 338R,519R,785R,926R,1492R \
    --rc-reverse \
    --min-length 100 \
    --max-length 2000 \
    --error-rate 0.1 \
    --min-overlap 10 \
    --threads 20 \
    --cutadapt cutadapt \
    --outdir amplicons_silva_27F_stepwise