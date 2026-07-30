#!/usr/bin/env bash
set -euo pipefail

mkdir -p results

FORCE=${FORCE:-0}
if [ "${FORCE}" = 1 ] || [ ! -s results/runtime.tsv ]
then
  printf "step\tseconds\n" > results/runtime.tsv
fi

record_time() {
  STEP=$1
  START=$2
  END=$(date +%s)
  printf "%s\t%s\n" "${STEP}" "$((END - START))" >> results/runtime.tsv
}

for FILE in data/example.bed data/example.bim data/example.fam software/dist/admixture_linux-1.3.0/admixture software/plink/plink software/PCAone software/evalAdmix/evalAdmix
do
  if [ ! -e "${FILE}" ]
  then
    echo "Missing ${FILE}. Run code/00_download_data.sh and code/00_install_software.sh first." >&2
    exit 1
  fi
done

if [ "${FORCE}" = 1 ] || [ ! -s results/example.qc.bed ]
then
  START=$(date +%s)
  software/plink/plink --bfile data/example \
    --maf 0.05 \
    --geno 0.01 \
    --make-bed \
    --out results/example.qc
  record_time "plink_maf_geno_qc" "${START}"
fi

if [ "${FORCE}" = 1 ] || [ ! -s results/example.pcaone_plot.eigvecs2 ]
then
  START=$(date +%s)
  software/PCAone -b results/example.qc -k 10 -V -o results/example.pcaone_plot
  record_time "pcaone_plot_10pcs" "${START}"
fi

if [ "${FORCE}" = 1 ] || [ ! -s results/example.hwe.hwe ]
then
  START=$(date +%s)
  software/PCAone -b results/example.qc \
    -k 5 \
    --inbreed 1 \
    -o results/example.hwe
  record_time "pcaone_hwe_test" "${START}"
fi

if [ "${FORCE}" = 1 ] || [ ! -s results/example.hwe.keep ]
then
  START=$(date +%s)
  awk 'NR > 1 && !($2 <= 1e-6 && ($4 <= -0.05 || $4 >= 0.05)) {print $1}' \
    results/example.hwe.hwe > results/example.hwe.keep
  record_time "hwe_filter_snplist" "${START}"
fi

if [ "${FORCE}" = 1 ] || [ ! -s results/example.hwe_filtered.bed ]
then
  START=$(date +%s)
  software/plink/plink --bfile results/example.qc \
    --extract results/example.hwe.keep \
    --make-bed \
    --out results/example.hwe_filtered
  record_time "plink_make_hwe_filtered" "${START}"
fi

if [ "${FORCE}" = 1 ] || [ ! -s results/example.adjld.residuals ]
then
  START=$(date +%s)
  software/PCAone -b results/example.hwe_filtered -k 5 --ld -o results/example.adjld
  record_time "pcaone_adjusted_ld_matrix" "${START}"
fi

if [ "${FORCE}" = 1 ] || [ ! -s results/example.adjld.ld.prune.in ]
then
  START=$(date +%s)
  software/PCAone -B results/example.adjld.residuals \
    --match-bim results/example.adjld.mbim \
    --ld-r2 0.2 \
    --ld-bp 1000000 \
    -o results/example.adjld
  record_time "pcaone_ld_prune" "${START}"
fi

if [ "${FORCE}" = 1 ] || [ ! -s results/example.adjld.ld.prune.ids ]
then
  START=$(date +%s)
  awk '{print $2}' results/example.adjld.ld.prune.in > results/example.adjld.ld.prune.ids
  record_time "ld_prune_extract_ids" "${START}"
fi

if [ "${FORCE}" = 1 ] || [ ! -s results/example.pcaone_pruned.bed ]
then
  START=$(date +%s)
  software/plink/plink --bfile results/example.hwe_filtered \
    --extract results/example.adjld.ld.prune.ids \
    --make-bed \
    --out results/example.pcaone_pruned
  record_time "plink_make_admixture_input" "${START}"
fi

QOUT=results/example.pcaone_pruned.K3.seed1.Q
POUT=results/example.pcaone_pruned.K3.seed1.P

if [ "${FORCE}" = 1 ] || [ ! -s "${QOUT}" ] || [ ! -s "${POUT}" ]
then
  START=$(date +%s)
  software/dist/admixture_linux-1.3.0/admixture -s 1 results/example.pcaone_pruned.bed 3 \
    | tee results/admixture.K3.seed1.log

  mv example.pcaone_pruned.3.Q "${QOUT}"
  mv example.pcaone_pruned.3.P "${POUT}"
  record_time "admixture_k3_seed1" "${START}"
fi

if [ "${FORCE}" = 1 ] || [ ! -s results/evaladmix.K3.seed1.corres ]
then
  START=$(date +%s)
  software/evalAdmix/evalAdmix \
    -plink results/example.pcaone_pruned \
    -fname results/example.pcaone_pruned.K3.seed1.P \
    -qname results/example.pcaone_pruned.K3.seed1.Q \
    -o results/evaladmix.K3.seed1.corres \
    -P 4
  record_time "evaladmix_k3_seed1" "${START}"
fi
