#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

printf "________________  __________  __________.____________  \n"
printf "\__    ___/  _  \ \_   ___ \ /__    ___/|   \_   ___ \ \n"
printf " |    | /  /_\  \/    \  \/  |    |   |   /    \  \/  \n"
printf " |    |/    |    \     \___  |    |   |   \     \____ \n"
printf " |____|\____|__  /\______    |____|   |___|\______  / \n"
printf "               \/        \/                       \/  \n"


exec /base/run_tactic.py "$@"
