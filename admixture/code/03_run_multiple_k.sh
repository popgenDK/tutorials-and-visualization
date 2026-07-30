#!/usr/bin/env bash
set -euo pipefail

FORCE=${FORCE:-0}

if [ ! -s results/runtime.tsv ]
then
  printf "step\tseconds\n" > results/runtime.tsv
fi

record_time() {
  STEP=$1
  START=$2
  END=$(date +%s)
  printf "%s\t%s\n" "${STEP}" "$((END - START))" >> results/runtime.tsv
}

for K in 2 3 4 5
do
  QOUT="results/example.pcaone_pruned.K${K}.seed1.Q"
  POUT="results/example.pcaone_pruned.K${K}.seed1.P"
  if [ "${FORCE}" != 1 ] && [ -s "${QOUT}" ] && [ -s "${POUT}" ]
  then
    echo "Skipping K=${K} seed=1; output already exists."
  else
    START=$(date +%s)
    software/dist/admixture_linux-1.3.0/admixture -s 1 results/example.pcaone_pruned.bed "${K}" \
      | tee "results/admixture.K${K}.seed1.log"

    mv "example.pcaone_pruned.${K}.Q" "${QOUT}"
    mv "example.pcaone_pruned.${K}.P" "${POUT}"
    record_time "admixture_k${K}_seed1" "${START}"
  fi

  if [ "${FORCE}" = 1 ] || [ ! -s "results/evaladmix.K${K}.seed1.corres" ]
  then
    START=$(date +%s)
    software/evalAdmix/evalAdmix \
      -plink results/example.pcaone_pruned \
      -fname "${POUT}" \
      -qname "${QOUT}" \
      -o "results/evaladmix.K${K}.seed1.corres" \
      -P 4
    record_time "evaladmix_k${K}_seed1" "${START}"
  fi
done
