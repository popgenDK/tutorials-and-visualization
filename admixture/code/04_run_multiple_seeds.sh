#!/usr/bin/env bash
set -euo pipefail

MAX_SEEDS=${MAX_SEEDS:-10}
FORCE=${FORCE:-0}
CONV_TIMES=${CONV_TIMES:-3}
LL_DIFF=${LL_DIFF:-3}
Q_DIFF=${Q_DIFF:-0.01}

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
  LIKE_TMP="results/admixture.K${K}.likes.tmp"
  QLIST="results/admixture.K${K}.Qlist"
  if [ "${FORCE}" = 1 ]
  then
    : > "${LIKE_FILE}"
    : > "${LIKE_TMP}"
    : > "${QLIST}"
  fi
  touch "${LIKE_FILE}"
  touch "${LIKE_TMP}"
  touch "${QLIST}"

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
        > "${LOG}"

      mv "example.pcaone_pruned.${K}.Q" "${QOUT}"
      mv "example.pcaone_pruned.${K}.P" "${POUT}"
      record_time "convergence_admixture_k${K}_seed${SEED}" "${START}"
    fi

    LOG_LIK=$(awk '/^Loglikelihood:/ {ll=$2} END {print ll}' "${LOG}")
    awk -v seed="${SEED}" -v ll="${LOG_LIK}" '$1 != seed {print} END {print seed, ll}' \
      "${LIKE_TMP}" 2>/dev/null > "${LIKE_TMP}.tmp"
    mv "${LIKE_TMP}.tmp" "${LIKE_TMP}"
    awk -v qout="${QOUT}" '$1 != qout {print} END {print qout}' \
      "${QLIST}" 2>/dev/null > "${QLIST}.tmp"
    mv "${QLIST}.tmp" "${QLIST}"

    awk -v seed="${SEED}" -v ll="${LOG_LIK}" '$1 != seed {print} END {print seed, ll}' \
      "${LIKE_FILE}" 2>/dev/null | sort -k2,2nr > "${LIKE_FILE}.tmp"
    mv "${LIKE_FILE}.tmp" "${LIKE_FILE}"

    awk '{print $2}' "${LIKE_TMP}" > "results/ll.K${K}.${SEED}.txt"
    CONV=$(awk -v diff="${LL_DIFF}" 'NR == 1 {best=$2} best - $2 < diff {n++} END {print n + 0}' "${LIKE_FILE}")
    CONV2=$(Rscript code/testQconv.R "results/ll.K${K}.${SEED}.txt" "${QLIST}" "${Q_DIFF}")
    echo "K=${K} seed=${SEED} likelihood_converged=${CONV} q_converged=${CONV2}"
    if [ "${CONV}" -gt "${CONV_TIMES}" ] || [ "${CONV2}" -gt "${CONV_TIMES}" ]
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
