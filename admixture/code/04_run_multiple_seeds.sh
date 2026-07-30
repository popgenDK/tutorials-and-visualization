#!/usr/bin/env bash
set -euo pipefail

MAX_SEEDS=${MAX_SEEDS:-10}
FORCE=${FORCE:-0}
THREADS=${THREADS:-4}

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
  START=$(date +%s)
  FORCE="${FORCE}" code/multiConv_ADMIXTURE.sh \
    results/example.pcaone_pruned "${MAX_SEEDS}" "${THREADS}" results "${K}" 1
  record_time "convergence_admixture_k${K}" "${START}"

  BEST_SEED=$(awk 'NR == 1 {print $1}' "results/admixture.K${K}.likes")
  START=$(date +%s)
  software/evalAdmix/evalAdmix \
    -plink results/example.pcaone_pruned \
    -fname "results/example.pcaone_pruned.K${K}.seed${BEST_SEED}.P" \
    -qname "results/example.pcaone_pruned.K${K}.seed${BEST_SEED}.Q" \
    -o "results/evaladmix.K${K}.best.corres" \
    -P 4
  record_time "convergence_evaladmix_k${K}_best_seed${BEST_SEED}" "${START}"
done
