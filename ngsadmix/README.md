# NGSadmix tutorial

This tutorial estimates ancestry proportions from genotype likelihoods using NGSadmix. It is adapted from the population-structure exercises in `popgenDK/courses`, especially `advBinf/exercises/admixture.md` and `summer2025/exercises/Day3_Morning_Admixture.ipynb`.

## Learning goals

- Create or inspect Beagle-format genotype likelihood input.
- Run NGSadmix for one or more values of `K`.
- Compare repeated runs using log likelihoods.
- Plot ancestry proportions and save the figure from reproducible R code.

## Input data

The course tutorials use low-depth 1000 Genomes data from ASW, CEU, CHB, YRI, and MXL. For this repository, large input files should be staged locally under `../tutorial_data/ngsadmix/` and mirrored publicly at:

```text
https://popgen.dk/albrecht/open/tutorial_data/ngsadmix/
```

Relevant course sources:

- https://github.com/popgenDK/courses/blob/main/advBinf/exercises/admixture.md
- https://github.com/popgenDK/courses/blob/main/summer2025/exercises/Day3_Morning_Admixture.ipynb

## Set paths

Run all commands from the `ngsadmix/` folder.

```bash
DATA_DIR=../tutorial_data/ngsadmix
DATA_URL=https://popgen.dk/albrecht/open/tutorial_data/ngsadmix
RESULTS_DIR=results
FIGURES_DIR=figures
CODE_DIR=code

THREADS=${THREADS:-4}
mkdir -p "${DATA_DIR}" "${RESULTS_DIR}" "${FIGURES_DIR}"

BEAGLE=${DATA_DIR}/1000G5pops.inputgl.beagle.gz
POPINFO=${DATA_DIR}/1000G5pops.pop.info
```

Download the data if they are not already present:

```bash
mkdir -p "${DATA_DIR}"
wget -nc -P "${DATA_DIR}" "${DATA_URL}/1000G5pops.inputgl.beagle.gz"
wget -nc -P "${DATA_DIR}" "${DATA_URL}/1000G5pops.pop.info"
```

<details>
<summary>Install software</summary>

NGSadmix documentation: https://popgen.dk/software/index.php/NgsAdmix

Install locally under `../software/`:

```bash
SOFTWARE_DIR=../software
mkdir -p "${SOFTWARE_DIR}/ngsadmix"
cd "${SOFTWARE_DIR}/ngsadmix"

wget -nc http://www.popgen.dk/software/download/NGSadmix/NGSadmix32.cpp
g++ -O3 NGSadmix32.cpp -lz -o NGSadmix

cd -
export PATH="$(pwd)/${SOFTWARE_DIR}/ngsadmix:${PATH}"

which NGSadmix
NGSadmix 2>&1 | head
```

If starting from BAM files, ANGSD is also needed to generate Beagle-format genotype likelihoods:

```bash
SOFTWARE_DIR=../software
git clone https://github.com/ANGSD/angsd.git "${SOFTWARE_DIR}/angsd"
cd "${SOFTWARE_DIR}/angsd"
make
cd -
export PATH="$(pwd)/${SOFTWARE_DIR}/angsd:${PATH}"

which angsd
angsd --help 2>&1 | head
```

</details>

## Inspect genotype likelihood input

```bash
zcat "${BEAGLE}" | wc -l
zcat "${BEAGLE}" | head -n 7 | cut -f1-9 | column -t
cut -f 1 -d " " "${POPINFO}" | sort | uniq -c
```

## Run NGSadmix

Run several seeds for a fixed `K`. The tutorial should start with `K=3` and then repeat with `K=4`.

```bash
K=3

for SEED in $(seq 1 20)
do
  NGSadmix \
    -likes "${BEAGLE}" \
    -K "${K}" \
    -P "${THREADS}" \
    -minMaf 0.05 \
    -seed "${SEED}" \
    -o "${RESULTS_DIR}/1000G5pops.ngsadmix.K${K}.seed${SEED}"
done
```

## Pick the best run

```bash
K=3
rm -f "${RESULTS_DIR}/ngsadmix.K${K}.likes"

for SEED in $(seq 1 20)
do
  grep "best like" "${RESULTS_DIR}/1000G5pops.ngsadmix.K${K}.seed${SEED}.log" \
    | awk -v seed="${SEED}" '{print seed, $NF}' \
    >> "${RESULTS_DIR}/ngsadmix.K${K}.likes"
done

sort -k2,2gr "${RESULTS_DIR}/ngsadmix.K${K}.likes"
```

## Make figures

Run the plotting script:

```bash
Rscript code/02_plot_ngsadmix.R
```

![NGSadmix K3 ancestry proportions](figures/ngsadmix_k3.png)

Save this as `code/02_plot_ngsadmix.R` and run it from `ngsadmix/`.

```r
source("https://raw.githubusercontent.com/GenisGE/evalAdmix/master/visFuns.R")

data_dir <- "../tutorial_data/ngsadmix"
results_dir <- "results"
figures_dir <- "figures"

pop <- read.table(file.path(data_dir, "1000G5pops.pop.info"), as.is = TRUE)
q <- read.table(file.path(results_dir, "1000G5pops.ngsadmix.K3.seed3.qopt"))

ord <- orderInds(pop = pop[, 1], q = q)

png(file.path(figures_dir, "ngsadmix_k3.png"), width = 1400, height = 500, res = 160)
par(mar = c(5, 4, 1, 1))
barplot(
  t(q)[, ord],
  col = 2:10,
  space = 0,
  border = NA,
  xlab = "Individuals",
  ylab = "Admixture proportions"
)
text(sort(tapply(seq_len(nrow(pop)), pop[ord, 1], mean)), -0.05, unique(pop[ord, 1]), xpd = TRUE)
abline(v = cumsum(sapply(unique(pop[ord, 1]), function(x) sum(pop[ord, 1] == x))), col = 1, lwd = 1.2)
dev.off()
```

## Interpret the result

The final tutorial should ask users to compare the inferred clusters with the population labels, identify admixed populations, and explain why more sites usually give cleaner ancestry profiles.

## Sources

- https://github.com/popgenDK/courses/blob/main/advBinf/exercises/admixture.md
- https://github.com/popgenDK/courses/blob/main/summer2025/exercises/Day3_Morning_Admixture.ipynb
- https://popgen.dk/software/index.php/NgsAdmix
