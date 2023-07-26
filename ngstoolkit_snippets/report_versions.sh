#!/usr/bin/env bash

echo -e "* WSL distro build version"
cat /usr/local/bin/wsl_distro_version.txt
echo -e "\n* USEARCH"
usearch  --version
echo -e "\n* SINA"
sina --version
echo -e "\n* SortMeRNA"
sortmerna --version
echo -e "\n* FastTree"
FastTree
echo -e "\n* bowtie2"
bowtie2 --version
