#!/usr/bin/env bash
set -euo pipefail

TUTORIAL_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
DATA_ROOT=${TUTORIAL_DIR}/../tutorial_data
DATA_DIR=${DATA_ROOT}/admixture
RESULTS_DIR=${TUTORIAL_DIR}/results
THREADS=${THREADS:-4}

mkdir -p "${RESULTS_DIR}"
cd "${TUTORIAL_DIR}"

PLINK_PREFIX=${DATA_DIR}/example
PRUNED_PREFIX=${RESULTS_DIR}/example.pruned

plink \
  --bfile "${PLINK_PREFIX}" \
  --indep-pairwise 50 10 0.1 \
  --out "${RESULTS_DIR}/example.ld"

plink \
  --bfile "${PLINK_PREFIX}" \
  --extract "${RESULTS_DIR}/example.ld.prune.in" \
  --make-bed \
  --out "${PRUNED_PREFIX}"

for K in 2 3 4 5
do
  admixture -j"${THREADS}" --cv "${PRUNED_PREFIX}.bed" "${K}" \
    | tee "${RESULTS_DIR}/admixture.K${K}.log"

  mv "example.pruned.${K}.Q" "${RESULTS_DIR}/example.pruned.K${K}.Q"
  mv "example.pruned.${K}.P" "${RESULTS_DIR}/example.pruned.K${K}.P"
done

grep -h "CV error" "${RESULTS_DIR}"/admixture.K*.log \
  > "${RESULTS_DIR}/admixture.cv_errors.txt"
