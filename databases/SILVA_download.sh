#! /bin/bash
curr_wd=$(pwd)
database_dir="/base/databases"
version=138
# Get this files parent path and store it in a variable
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
wget --no-check-certificate https://www.arb-silva.de/fileadmin/arb_web_db/release_138_1/ARB_files/SILVA_${version}.1_SSURef_NR99_12_06_20_opt.arb.gz
if [[ ! -f SILVA_${version}.1_SSURef_NR99_12_06_20_opt.arb.gz ]]; then
    echo "SILVA_${version}.1_SSURef_NR99_12_06_20_opt.arb.gz not found"
    exit 1
else
    # Unzip the file
    gunzip -d SILVA_${version}.1_SSURef_NR99_12_06_20_opt.arb.gz -c > SILVA_${version}.1_SSURef_NR99_12_06_20_opt.arb
    echo "SILVA_${version}.1_SSURef_NR99_12_06_20_opt.arb.gz gunzipped"
fi

if [[ -f silva-arc-16s-id95.fasta.gz ]]; then
    gunzip -d silva-arc-16s-id95.fasta.gz -c > silva-arc-16s-id95.fasta
fi
if [[ -f silva-bac-16s-id90.fasta.gz ]]; then
    gunzip -d silva-bac-16s-id90.fasta.gz -c > silva-bac-16s-id90.fasta
fi

# Creating a link in dataset_dir to these three files
ln -s ${parent_path}/SILVA_${version}.1_SSURef_NR99_12_06_20_opt.arb ${database_dir}/SILVA_${version}.1_SSURef_NR99_12_06_20_opt.arb
ln -s ${parent_path}/silva-arc-16s-id95.fasta ${database_dir}/silva-arc-16s-id95.fasta
ln -s ${parent_path}/silva-bac-16s-id90.fasta ${database_dir}/silva-bac-16s-id90.fasta



