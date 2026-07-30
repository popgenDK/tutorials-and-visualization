#!/usr/bin/env bash
set -euo pipefail

mkdir -p software

if [ ! -x software/dist/admixture_linux-1.3.0/admixture ]
then
  curl -L https://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz \
    -o software/admixture_linux-1.3.0.tar.gz
  tar -xzf software/admixture_linux-1.3.0.tar.gz -C software
fi

if [ ! -x software/plink/plink ]
then
  curl -L https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20241022.zip \
    -o software/plink_linux_x86_64_20241022.zip
  unzip -n software/plink_linux_x86_64_20241022.zip -d software/plink
fi

if [ ! -x software/PCAone ]
then
  curl -L https://github.com/Zilong-Li/PCAone/releases/latest/download/PCAone-Linux.zip \
    -o software/PCAone-Linux.zip
  unzip -n software/PCAone-Linux.zip -d software
  chmod +x software/PCAone
fi

if [ ! -x software/evalAdmix/evalAdmix ]
then
  git clone https://github.com/GenisGE/evalAdmix.git software/evalAdmix
  make -C software/evalAdmix
fi

software/dist/admixture_linux-1.3.0/admixture --help >/dev/null || true
software/plink/plink --help >/dev/null || true
software/PCAone --help >/dev/null || true
software/evalAdmix/evalAdmix 2>&1 | head >/dev/null || true
