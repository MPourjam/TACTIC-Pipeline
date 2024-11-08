#! /bin/bash
curr_wd=$(pwd)
database_dir="/base/databases"
major_version=138
minor_version=2
version=${major_version}.${minor_version}
version_date=03_07_24
SILVA_DB_FILE="SILVA_${version}_SSURef_NR99_${version_date}_opt.arb.gz"
SILVA_ARB_API_ROOT="https://www.arb-silva.de/fileadmin/silva_databases"
SILVA_DB_URL="${SILVA_ARB_API_ROOT}/release_${major_version}_${minor_version}/ARB_files/${SILVA_DB_FILE}"
SILVA_ARB_FILE_LATEST="SILVA_LATEST.arb"
SILVA_DB_INDEX_FILE="${SILVA_ARB_FILE_LATEST%.arb}.sidx"
SINA_BIN="/base/binaries/sina/bin/sina"
# Get this files parent path and store it in a variable
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd ${parent_path}
if [[ ! -f ${SILVA_DB_FILE} ]]; then
    echo "${SILVA_DB_FILE} not found"
    # wget --no-check-certificate https://www.arb-silva.de/fileadmin/arb_web_db/release_138_1/ARB_files/${SILVA_DB_FILE}.gz
    wget "${SILVA_DB_URL}"
fi

# Gunzip the file
gunzip -d ${SILVA_DB_FILE} -c > ${SILVA_ARB_FILE_LATEST} && rm ${SILVA_DB_FILE}
echo "${SILVA_DB_FILE} gunzipped"

# Creating database index
touch fake1.fasta
${SINA_BIN} --db ${SILVA_ARB_FILE_LATEST} --in fake1.fasta --out fake2.fasta
rm fake{1..2}.fasta

# Unzipping the SILVA database index
if [[ ! -f ${SILVA_DB_INDEX_FILE} ]]; then
    echo "Could not build index for ${SILVA_ARB_FILE_LATEST}"
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
if [[ -f ${database_dir}/${SILVA_ARB_FILE_LATEST} ]]; then
    if [[ $(diff ${parent_path}/${SILVA_ARB_FILE_LATEST} ${database_dir}/${SILVA_ARB_FILE_LATEST}) ]]; then
        echo "${SILVA_ARB_FILE_LATEST} is different"
        rm ${database_dir}/${SILVA_ARB_FILE_LATEST}
        ln -sf ${parent_path}/${SILVA_ARB_FILE_LATEST} ${database_dir}/${SILVA_ARB_FILE_LATEST}
    fi
else
    ln -sf ${parent_path}/${SILVA_ARB_FILE_LATEST} ${database_dir}/${SILVA_ARB_FILE_LATEST}
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

# ln -sf ${parent_path}/${SILVA_ARB_FILE_LATEST} ${database_dir}/${SILVA_ARB_FILE_LATEST}
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
