#!/usr/bin/env bash
set -euo pipefail

mkdir -p results

plink --bfile data/example \
  --maf 0.05 \
  --geno 0.01 \
  --make-bed \
  --out results/example.qc

PCAone -b results/example.qc -k 10 -V -o results/example.pcaone

PCAone -b results/example.qc \
  --USV results/example.pcaone \
  -k 10 \
  --inbreed 1 \
  -o results/example.hwe

awk 'NR > 1 && $2 >= 1e-6 {print $1}' results/example.hwe.hwe > results/example.hwe.keep

PCAone -b results/example.qc -k 10 --ld -o results/example.adjld

PCAone -B results/example.adjld.residuals \
  --match-bim results/example.adjld.mbim \
  --ld-r2 0.2 \
  --ld-bp 1000000 \
  -o results/example.adjld

grep -Fxf results/example.adjld.ld.prune.in results/example.hwe.keep > results/example.keep.snps

plink --bfile results/example.qc \
  --extract results/example.keep.snps \
  --make-bed \
  --out results/example.pcaone_pruned

for K in 2 3 4 5
do
  for SEED in $(seq 1 10)
  do
    admixture -s "${SEED}" results/example.pcaone_pruned.bed "${K}" \
      | tee "results/admixture.K${K}.seed${SEED}.log"

    mv "example.pcaone_pruned.${K}.Q" "results/example.pcaone_pruned.K${K}.seed${SEED}.Q"
    mv "example.pcaone_pruned.${K}.P" "results/example.pcaone_pruned.K${K}.seed${SEED}.P"
  done
done
