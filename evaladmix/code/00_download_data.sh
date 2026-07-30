#!/usr/bin/env bash
set -euo pipefail

DATA_DIR=../tutorial_data/ngsadmix
DATA_URL=https://popgen.dk/albrecht/open/tutorial_data/ngsadmix

mkdir -p "${DATA_DIR}"

wget -nc -P "${DATA_DIR}" "${DATA_URL}/1000G5pops.inputgl.beagle.gz"
wget -nc -P "${DATA_DIR}" "${DATA_URL}/1000G5pops.pop.info"
