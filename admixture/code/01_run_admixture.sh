#!/usr/bin/env bash
set -euo pipefail

mkdir -p results

ADMIXTURE=software/dist/admixture_linux-1.3.0/admixture
PLINK=software/plink/plink
PCAONE=software/PCAone

for FILE in data/example.bed data/example.bim data/example.fam "${ADMIXTURE}" "${PLINK}" "${PCAONE}"
do
  if [ ! -e "${FILE}" ]
  then
    echo "Missing ${FILE}. Run code/00_download_data.sh and code/00_install_software.sh first." >&2
    exit 1
  fi
done

if [ ! -s results/example.qc.bed ]
then
  "${PLINK}" --bfile data/example \
    --maf 0.05 \
    --geno 0.01 \
    --make-bed \
    --out results/example.qc
fi

if [ ! -s results/example.pcaone.eigvecs2 ]
then
  "${PCAONE}" -b results/example.qc -k 10 -V -o results/example.pcaone
fi

if [ ! -s results/example.hwe.hwe ]
then
  "${PCAONE}" -b results/example.qc \
    --USV results/example.pcaone \
    -k 10 \
    --inbreed 1 \
    -o results/example.hwe
fi

if [ ! -s results/example.hwe.keep ]
then
  awk 'NR > 1 && $2 >= 1e-6 {print $1}' results/example.hwe.hwe > results/example.hwe.keep
fi

if [ ! -s results/example.adjld.residuals ]
then
  "${PCAONE}" -b results/example.qc -k 10 --ld -o results/example.adjld
fi

if [ ! -s results/example.adjld.ld.prune.in ]
then
  "${PCAONE}" -B results/example.adjld.residuals \
    --match-bim results/example.adjld.mbim \
    --ld-r2 0.2 \
    --ld-bp 1000000 \
    -o results/example.adjld
fi

if [ ! -s results/example.adjld.ld.prune.ids ]
then
  awk '{print $2}' results/example.adjld.ld.prune.in > results/example.adjld.ld.prune.ids
fi

if [ ! -s results/example.keep.snps ]
then
  grep -Fxf results/example.adjld.ld.prune.ids results/example.hwe.keep > results/example.keep.snps
fi

if [ ! -s results/example.pcaone_pruned.bed ]
then
  "${PLINK}" --bfile results/example.qc \
    --extract results/example.keep.snps \
    --make-bed \
    --out results/example.pcaone_pruned
fi

for K in 2 3 4 5
do
  for SEED in $(seq 1 10)
  do
    QOUT="results/example.pcaone_pruned.K${K}.seed${SEED}.Q"
    POUT="results/example.pcaone_pruned.K${K}.seed${SEED}.P"
    if [ -s "${QOUT}" ] && [ -s "${POUT}" ]
    then
      echo "Skipping K=${K} seed=${SEED}; output already exists."
      continue
    fi

    "${ADMIXTURE}" -s "${SEED}" results/example.pcaone_pruned.bed "${K}" \
      | tee "results/admixture.K${K}.seed${SEED}.log"

    mv "example.pcaone_pruned.${K}.Q" "${QOUT}"
    mv "example.pcaone_pruned.${K}.P" "${POUT}"
  done
done
