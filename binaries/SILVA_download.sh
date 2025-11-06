#! /bin/bash

# Enable strict error handling
set -euo pipefail
IFS=$'\n\t'

log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - TACTIC.Databases:$(basename "${BASH_SOURCE[0]}"):${LINENO} - INFO - $*"
}

curr_wd=$(pwd)
database_dir="/base/databases_raw"
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
parent_path="/base/inputs/databases" #  It should be sub-directory of /base/inputs as we mount /base/inputs to /base/databases in the container
# Create the parent path if it does not exist
mkdir -p "${parent_path}"
# cd to the parent path
cd "${parent_path}"
if [[ ! -s "${SILVA_ARB_FILE_LATEST}" || ! -s "${SILVA_DB_INDEX_FILE}" ]]; then
    log "${SILVA_ARB_FILE_LATEST} or ${SILVA_DB_INDEX_FILE} not found"
    # wget --no-check-certificate https://www.arb-silva.de/fileadmin/arb_web_db/release_138_1/ARB_files/${SILVA_DB_FILE}.gz
    log "Downloading SILVA database ${SILVA_DB_URL}"
    wget "${SILVA_DB_URL}" -O "${SILVA_DB_FILE}" > download_SILVA.log 2>&1
    gunzip -d "${SILVA_DB_FILE}" -c > "${SILVA_ARB_FILE_LATEST}" && rm "${SILVA_DB_FILE}"
    # Creating database index
    touch fake1.fasta
    log "Building SILVA database index ${SILVA_DB_INDEX_FILE}"
    "${SINA_BIN}" --db "${SILVA_ARB_FILE_LATEST}" --in fake1.fasta --out fake2.fasta > SILVA_index_build.log 2>&1
    rm fake{1..2}.fasta
else
    log "SILVA database and index found, skipping download and index creation"
fi

log "Checking for silva-arc-16s-id95.fasta and silva-bac-16s-id90.fasta"
if [[ ! -s "${parent_path}/silva-arc-16s-id95.fasta" ]]; then
    # We delete after unzipping to as we force link /base/inputs/databases to /base/databases
    log "Placing silva-arc-16s-id95.fasta.gz in ${parent_path}"
    gunzip -k -d "${database_dir}/silva-arc-16s-id95.fasta.gz" -c > "${parent_path}/silva-arc-16s-id95.fasta"
else
    log "silva-arc-16s-id95.fasta found, skipping download"
fi

if [[ ! -s "${parent_path}/silva-bac-16s-id90.fasta" ]]; then
    # We delete after unzipping to as we force link /base/inputs/databases to /base/databases
    log "Placing silva-bac-16s-id90.fasta.gz in ${parent_path}"
    gunzip -k -d "${database_dir}/silva-bac-16s-id90.fasta.gz" -c > "${parent_path}/silva-bac-16s-id90.fasta"
else
    log "silva-bac-16s-id90.fasta found, skipping download"
fi

# cd back to the current working directory
cd "${curr_wd}"
exit 0
