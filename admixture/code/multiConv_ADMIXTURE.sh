#!/usr/bin/env bash
set -euo pipefail

file=$1
num=$2
P=$3
out=$4
K=$5
star=$6
bfile=$(basename "${file}")

# Based on popgenDK/analysis Admix/multiConv_ADMIXTURE.sh, with repo-local paths.
ADM=${ADMIXTURE:-software/dist/admixture_linux-1.3.0/admixture}
CONV_TIMES=${CONV_TIMES:-3}
LL_DIFF=${LL_DIFF:-3}
Q_DIFF=${Q_DIFF:-0.01}
FORCE=${FORCE:-0}
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

echo "num = ${num}"
echo "p = ${P}"
echo "out = ${out}"
echo "K = ${K}"

mkdir -p "${out}"
like_tmp="${out}/admixture.K${K}.likes.tmp"
like_file="${out}/admixture.K${K}.likes"
q_list="${out}/admixture.K${K}.Qlist"

if [ "${FORCE}" = 1 ]
then
  : > "${like_tmp}"
  : > "${like_file}"
  : > "${q_list}"
fi
touch "${like_tmp}" "${like_file}" "${q_list}"
: > "${like_tmp}"
: > "${like_file}"
: > "${q_list}"

expected_snps=$(awk 'END {print NR}' "${file}.bim")
expected_individuals=$(awk 'END {print NR}' "${file}.fam")

outputs_match_input() {
  qout=$1
  pout=$2
  log=$3
  [ -s "${qout}" ] || return 1
  [ -s "${pout}" ] || return 1
  [ -s "${log}" ] || return 1
  [ "$(awk 'END {print NR}' "${qout}")" = "${expected_individuals}" ] || return 1
  [ "$(awk 'END {print NR}' "${pout}")" = "${expected_snps}" ] || return 1
}

for f in $(seq "${star}" "${num}")
do
  qout="${out}/${bfile}.K${K}.seed${f}.Q"
  pout="${out}/${bfile}.K${K}.seed${f}.P"
  log="${out}/admixture.K${K}.seed${f}.log"

  if [ "${FORCE}" != 1 ] && outputs_match_input "${qout}" "${pout}" "${log}"
  then
    echo "Skipping K=${K} seed=${f}; output already exists."
  else
    "${ADM}" "${file}.bed" "${K}" -s "${f}" -j"${P}" > "${log}"
    mv "${bfile}.${K}.Q" "${qout}"
    mv "${bfile}.${K}.P" "${pout}"
  fi

  log_lik=$(awk '/^Loglikelihood/ {ll=$2} END {print ll}' "${log}")
  awk -v seed="${f}" -v ll="${log_lik}" '$1 != seed {print} END {print seed, ll}' \
    "${like_tmp}" > "${like_tmp}.tmp"
  mv "${like_tmp}.tmp" "${like_tmp}"
  sort -k2,2nr "${like_tmp}" > "${like_file}"

  awk -v qout="${qout}" '$1 != qout {print} END {print qout}' \
    "${q_list}" > "${q_list}.tmp"
  mv "${q_list}.tmp" "${q_list}"

  awk '{print $2}' "${like_tmp}" > "${out}/ll.K${K}.${f}.txt"
  conv=$(awk -v diff="${LL_DIFF}" 'NR == 1 {best=$2} best - $2 < diff {n++} END {print n + 0}' "${like_file}")
  echo "conv ${conv} with chosen ${CONV_TIMES}"

  conv2=$(Rscript "${script_dir}/testQconv.R" "${out}/ll.K${K}.${f}.txt" "${q_list}" "${Q_DIFF}")
  echo "second criteria conv: ${conv2} with chosen ${CONV_TIMES}"

  if [ "${conv}" -gt "${CONV_TIMES}" ] || [ "${conv2}" -gt "${CONV_TIMES}" ]
  then
    echo "k: ${K} first conv criteria ${conv} second conv criteria ${conv2}"
    break
  fi
done
