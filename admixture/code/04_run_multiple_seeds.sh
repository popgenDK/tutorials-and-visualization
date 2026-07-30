#!/usr/bin/env bash
set -euo pipefail

MAX_SEEDS=${MAX_SEEDS:-10}
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
  LIKE_FILE="results/admixture.K${K}.likes"
  if [ "${FORCE}" = 1 ]
  then
    : > "${LIKE_FILE}"
  fi
  touch "${LIKE_FILE}"

  for SEED in $(seq 1 "${MAX_SEEDS}")
  do
    QOUT="results/example.pcaone_pruned.K${K}.seed${SEED}.Q"
    POUT="results/example.pcaone_pruned.K${K}.seed${SEED}.P"
    LOG="results/admixture.K${K}.seed${SEED}.log"

    if [ "${FORCE}" != 1 ] && [ -s "${QOUT}" ] && [ -s "${POUT}" ] && [ -s "${LOG}" ]
    then
      echo "Skipping K=${K} seed=${SEED}; output already exists."
    else
      START=$(date +%s)
      software/dist/admixture_linux-1.3.0/admixture -s "${SEED}" results/example.pcaone_pruned.bed "${K}" \
        | tee "${LOG}"

      mv "example.pcaone_pruned.${K}.Q" "${QOUT}"
      mv "example.pcaone_pruned.${K}.P" "${POUT}"
      record_time "convergence_admixture_k${K}_seed${SEED}" "${START}"
    fi

    LOG_LIK=$(awk '/^Loglikelihood:/ {ll=$2} END {print ll}' "${LOG}")
    awk -v seed="${SEED}" -v ll="${LOG_LIK}" '$1 != seed {print} END {print seed, ll}' \
      "${LIKE_FILE}" 2>/dev/null | sort -k2,2nr > "${LIKE_FILE}.tmp"
    mv "${LIKE_FILE}.tmp" "${LIKE_FILE}"

    CONV=$(awk 'NR == 1 {best=$2} best - $2 < 5 {n++} END {print n + 0}' "${LIKE_FILE}")
    if [ "${CONV}" -ge 3 ]
    then
      break
    fi
  done

  BEST_SEED=$(awk 'NR == 1 {print $1}' "${LIKE_FILE}")
  START=$(date +%s)
  software/evalAdmix/evalAdmix \
    -plink results/example.pcaone_pruned \
    -fname "results/example.pcaone_pruned.K${K}.seed${BEST_SEED}.P" \
    -qname "results/example.pcaone_pruned.K${K}.seed${BEST_SEED}.Q" \
    -o "results/evaladmix.K${K}.best.corres" \
    -P 4
  record_time "convergence_evaladmix_k${K}_best_seed${BEST_SEED}" "${START}"
done
