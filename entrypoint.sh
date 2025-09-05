#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

printf "________________  __________  __________.____________  \n"
printf "\__    ___/  _  \ \_   ___ \ /__    ___/|   \_   ___ \ \n"
printf " |    | /  /_\  \/    \  \/  |    |   |   /    \  \/  \n"
printf " |    |/    |    \     \___  |    |   |   \     \____ \n"
printf " |____|\____|__  /\______    |____|   |___|\______  / \n"
printf "               \/        \/                       \/  \n"



# Call the SILVA database download script to ensure the database is present
/base/binaries/SILVA_download.sh

exec /base/run_imngs2.py "$@"
