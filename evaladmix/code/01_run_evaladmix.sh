#!/usr/bin/env bash
set -euo pipefail

TUTORIAL_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
NGSADMIX_DIR=${TUTORIAL_DIR}/../ngsadmix
DATA_ROOT=${TUTORIAL_DIR}/../tutorial_data
DATA_DIR=${DATA_ROOT}/ngsadmix
NGSADMIX_RESULTS_DIR=${NGSADMIX_DIR}/results
RESULTS_DIR=${TUTORIAL_DIR}/results
THREADS=${THREADS:-4}

mkdir -p "${RESULTS_DIR}"

K=${K:-3}
SEED=${SEED:-3}
BEAGLE=${DATA_DIR}/1000G5pops.inputgl.beagle.gz
QOPT=${NGSADMIX_RESULTS_DIR}/1000G5pops.ngsadmix.K${K}.seed${SEED}.qopt
FOPT=${NGSADMIX_RESULTS_DIR}/1000G5pops.ngsadmix.K${K}.seed${SEED}.fopt.gz

evalAdmix \
  -beagle "${BEAGLE}" \
  -fname "${FOPT}" \
  -qname "${QOPT}" \
  -o "${RESULTS_DIR}/1000G5pops.K${K}.seed${SEED}.corres" \
  -P "${THREADS}"
