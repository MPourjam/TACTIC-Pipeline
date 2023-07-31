#! /bin/bash
curr_wd=$(pwd)
database_dir="/base/databases"
version=138
cd ${database_dir}
wget https://www.arb-silva.de/fileadmin/arb_web_db/release_138_1/ARB_files/SILVA_${version}.1_SSURef_NR99_12_06_20_opt.arb.gz
gunzip -d SILVA_${version}.1_SSURef_NR99_12_06_20_opt.arb.gz
gunzip -d silva-arc-16s-id95.fasta.gz
gunzip -d silva-bac-16s-id90.fasta.gz
cd ${curr_wd}


