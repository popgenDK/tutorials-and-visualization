#!/usr/bin/env bash
set -euo pipefail

DATA_URL=https://popgen.dk/albrecht/open/tutorial_data/admixture

mkdir -p data

for EXT in bed bim fam
do
  if [ ! -s "data/example.${EXT}" ]
  then
    wget -nc -P data "${DATA_URL}/example.${EXT}"
  fi
done
