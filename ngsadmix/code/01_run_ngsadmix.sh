#!/usr/bin/env bash
set -euo pipefail

DATA_DIR=../tutorial_data/ngsadmix
RESULTS_DIR=results
THREADS=${THREADS:-4}

mkdir -p "${RESULTS_DIR}"

BEAGLE=${DATA_DIR}/1000G5pops.inputgl.beagle.gz
K=${K:-3}
SEEDS=${SEEDS:-20}

for SEED in $(seq 1 "${SEEDS}")
do
  NGSadmix \
    -likes "${BEAGLE}" \
    -K "${K}" \
    -P "${THREADS}" \
    -minMaf 0.05 \
    -seed "${SEED}" \
    -o "${RESULTS_DIR}/1000G5pops.ngsadmix.K${K}.seed${SEED}"
done

rm -f "${RESULTS_DIR}/ngsadmix.K${K}.likes"

for SEED in $(seq 1 "${SEEDS}")
do
  grep "best like" "${RESULTS_DIR}/1000G5pops.ngsadmix.K${K}.seed${SEED}.log" \
    | awk -v seed="${SEED}" '{print seed, $NF}' \
    >> "${RESULTS_DIR}/ngsadmix.K${K}.likes"
done

sort -k2,2gr "${RESULTS_DIR}/ngsadmix.K${K}.likes"
