#!/usr/bin/env bash
set -euo pipefail

THREADS=${THREADS:-4}

mkdir -p results

PLINK_PREFIX=data/example
PRUNED_PREFIX=results/example.pcaone_pruned

if [ ! -f "${PRUNED_PREFIX}.bed" ]; then
  echo "Missing ${PRUNED_PREFIX}.bed"
  echo "Create this PCAone-pruned PLINK data set before running ADMIXTURE."
  exit 1
fi

for K in 2 3 4 5
do
  admixture -j"${THREADS}" "${PRUNED_PREFIX}.bed" "${K}" \
    | tee "results/admixture.K${K}.log"

  mv "example.pcaone_pruned.${K}.Q" "results/example.pcaone_pruned.K${K}.Q"
  mv "example.pcaone_pruned.${K}.P" "results/example.pcaone_pruned.K${K}.P"
done
