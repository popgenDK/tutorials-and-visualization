#!/usr/bin/env bash
set -euo pipefail

mkdir -p software

ADMIXTURE=software/dist/admixture_linux-1.3.0/admixture
PLINK=software/plink/plink
PCAONE=software/PCAone

if [ ! -x "${ADMIXTURE}" ]
then
  curl -L https://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz \
    -o software/admixture_linux-1.3.0.tar.gz
  tar -xzf software/admixture_linux-1.3.0.tar.gz -C software
fi

if [ ! -x "${PLINK}" ]
then
  curl -L https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20241022.zip \
    -o software/plink_linux_x86_64_20241022.zip
  unzip -n software/plink_linux_x86_64_20241022.zip -d software/plink
fi

if [ ! -x "${PCAONE}" ]
then
  curl -L https://github.com/Zilong-Li/PCAone/releases/latest/download/PCAone-Linux.zip \
    -o software/PCAone-Linux.zip
  unzip -n software/PCAone-Linux.zip -d software
  chmod +x "${PCAONE}"
fi

"${ADMIXTURE}" --help >/dev/null || true
"${PLINK}" --help >/dev/null || true
"${PCAONE}" --help >/dev/null || true
