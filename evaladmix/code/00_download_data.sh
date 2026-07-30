#!/usr/bin/env bash
set -euo pipefail

DATA_URL=https://popgen.dk/albrecht/open/tutorial_data/ngsadmix

mkdir -p ../ngsadmix/data

wget -nc -P ../ngsadmix/data "${DATA_URL}/1000G5pops.inputgl.beagle.gz"
wget -nc -P ../ngsadmix/data "${DATA_URL}/1000G5pops.pop.info"
