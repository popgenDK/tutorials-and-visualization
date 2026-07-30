#!/usr/bin/env bash
set -euo pipefail

DATA_URL=https://popgen.dk/albrecht/open/tutorial_data/ngsadmix

mkdir -p data

wget -nc -P data "${DATA_URL}/1000G5pops.inputgl.beagle.gz"
wget -nc -P data "${DATA_URL}/1000G5pops.pop.info"
