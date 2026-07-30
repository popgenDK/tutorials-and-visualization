#!/usr/bin/env bash
set -euo pipefail

NGSADMIX_RESULTS_DIR=../ngsadmix/results
RESULTS_DIR=results
THREADS=${THREADS:-4}

mkdir -p "${RESULTS_DIR}"

K=${K:-3}
SEED=${SEED:-3}
BEAGLE=../ngsadmix/data/1000G5pops.inputgl.beagle.gz
QOPT=${NGSADMIX_RESULTS_DIR}/1000G5pops.ngsadmix.K${K}.seed${SEED}.qopt
FOPT=${NGSADMIX_RESULTS_DIR}/1000G5pops.ngsadmix.K${K}.seed${SEED}.fopt.gz

evalAdmix \
  -beagle "${BEAGLE}" \
  -fname "${FOPT}" \
  -qname "${QOPT}" \
  -o "${RESULTS_DIR}/1000G5pops.K${K}.seed${SEED}.corres" \
  -P "${THREADS}"
