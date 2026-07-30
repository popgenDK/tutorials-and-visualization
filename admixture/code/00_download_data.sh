#!/usr/bin/env bash
set -euo pipefail

DATA_DIR=../tutorial_data/admixture
DATA_URL=https://popgen.dk/albrecht/open/tutorial_data/admixture

mkdir -p "${DATA_DIR}"

wget -nc -P "${DATA_DIR}" "${DATA_URL}/example.bed"
wget -nc -P "${DATA_DIR}" "${DATA_URL}/example.bim"
wget -nc -P "${DATA_DIR}" "${DATA_URL}/example.fam"
