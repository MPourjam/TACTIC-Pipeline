#! /bin/bash
curr_wd=$(pwd)
database_dir="/base/databases"
version=138
SILVA_DB_FILE="SILVA_${version}.1_SSURef_NR99_12_06_20_opt.arb"
SILVA_DB_INDEX_FILE="${SILVA_DB_FILE%.arb}.sidx"
SINA_BIN="/base/binaries/sina/bin/sina"
# Get this files parent path and store it in a variable
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd ${parent_path}
if [[ ! -f ${SILVA_DB_FILE}.gz ]]; then
    echo "SILVA_${version}.1_SSURef_NR99_12_06_20_opt.arb.gz not found"
    # wget --no-check-certificate https://www.arb-silva.de/fileadmin/arb_web_db/release_138_1/ARB_files/${SILVA_DB_FILE}.gz
    wget https://www.arb-silva.de/fileadmin/arb_web_db/release_138_1/ARB_files/${SILVA_DB_FILE}.gz
fi

# Gunzip the file
gunzip -d ${SILVA_DB_FILE}.gz -c > ${SILVA_DB_FILE} && rm ${SILVA_DB_FILE}.gz
echo "${SILVA_DB_FILE}.gz gunzipped"

# Creating database index
touch fake1.fasta
${SINA_BIN} --db ${SILVA_DB_FILE} --in fake1.fasta --out fake2.fasta
rm fake{1..2}.fasta

# Unzipping the SILVA database index
if [[ ! -f ${SILVA_DB_INDEX_FILE} ]]; then
    echo "Could not build index for ${SILVA_DB_FILE}"
    exit 1
fi

if [[ -f silva-arc-16s-id95.fasta.gz ]]; then
    gunzip -d silva-arc-16s-id95.fasta.gz -c > silva-arc-16s-id95.fasta && rm silva-arc-16s-id95.fasta.gz
fi
if [[ -f silva-bac-16s-id90.fasta.gz ]]; then
    gunzip -d silva-bac-16s-id90.fasta.gz -c > silva-bac-16s-id90.fasta && rm silva-bac-16s-id90.fasta.gz
fi

# Creating a link in dataset_dir to these three files
# if the files are the same don't create the link
if [[ -f ${database_dir}/${SILVA_DB_FILE} ]]; then
    if [[ $(diff ${parent_path}/${SILVA_DB_FILE} ${database_dir}/${SILVA_DB_FILE}) ]]; then
        echo "${SILVA_DB_FILE} is different"
        rm ${database_dir}/${SILVA_DB_FILE}
        ln -sf ${parent_path}/${SILVA_DB_FILE} ${database_dir}/${SILVA_DB_FILE}
    fi
else
    ln -sf ${parent_path}/${SILVA_DB_FILE} ${database_dir}/${SILVA_DB_FILE}
fi

if [[ -f ${database_dir}/${SILVA_DB_INDEX_FILE} ]]; then
    if [[ $(diff ${parent_path}/${SILVA_DB_INDEX_FILE} ${database_dir}/${SILVA_DB_INDEX_FILE}) ]]; then
        echo "${SILVA_DB_INDEX_FILE} is different"
        rm ${database_dir}/${SILVA_DB_INDEX_FILE}
        ln -sf ${parent_path}/${SILVA_DB_INDEX_FILE} ${database_dir}/${SILVA_DB_INDEX_FILE}
    fi
else
    ln -sf ${parent_path}/${SILVA_DB_INDEX_FILE} ${database_dir}/${SILVA_DB_INDEX_FILE}
fi

# ln -sf ${parent_path}/SILVA_${version}.1_SSURef_NR99_12_06_20_opt.arb ${database_dir}/SILVA_${version}.1_SSURef_NR99_12_06_20_opt.arb
if [[ -f ${database_dir}/silva-arc-16s-id95.fasta ]]; then
    if [[ $(diff ${parent_path}/silva-arc-16s-id95.fasta ${database_dir}/silva-arc-16s-id95.fasta) ]]; then
        echo "silva-arc-16s-id95.fasta is different"
        rm ${database_dir}/silva-arc-16s-id95.fasta
        ln -sf ${parent_path}/silva-arc-16s-id95.fasta ${database_dir}/silva-arc-16s-id95.fasta
    fi
else
    ln -sf ${parent_path}/silva-arc-16s-id95.fasta ${database_dir}/silva-arc-16s-id95.fasta
fi
# ln -sf ${parent_path}/silva-arc-16s-id95.fasta ${database_dir}/silva-arc-16s-id95.fasta
if [[ -f ${database_dir}/silva-bac-16s-id90.fasta ]]; then
    if [[ $(diff ${parent_path}/silva-bac-16s-id90.fasta ${database_dir}/silva-bac-16s-id90.fasta) ]]; then
        echo "silva-bac-16s-id90.fasta is different"
        rm ${database_dir}/silva-bac-16s-id90.fasta
        ln -sf ${parent_path}/silva-bac-16s-id90.fasta ${database_dir}/silva-bac-16s-id90.fasta
    fi
else
    ln -sf ${parent_path}/silva-bac-16s-id90.fasta ${database_dir}/silva-bac-16s-id90.fasta
fi
