#!/usr/bin/env bash
set -euo pipefail

DATA_URL=https://popgen.dk/albrecht/open/tutorial_data/admixture

mkdir -p data

if [ ! -s data/example.bed ]
then
  wget -nc -P data "${DATA_URL}/example.bed"
  wget -nc -P data "${DATA_URL}/example.bim"
  wget -nc -P data "${DATA_URL}/example.fam"
fi
