# ADMIXTURE tutorial

This tutorial estimates ancestry proportions from called genotypes using ADMIXTURE. It complements the NGSadmix tutorial, which uses genotype likelihoods from low-depth sequencing data.

## Learning goals

- Prepare PLINK input for ADMIXTURE.
- LD prune genotype data before ancestry inference.
- Run ADMIXTURE for several values of `K`.
- Plot ancestry proportions from ADMIXTURE `Q` files.
- Compare genotype-based ADMIXTURE output with NGSadmix-style output.

## Input data

The final tutorial should use a PLINK binary data set staged locally under `../tutorial_data/admixture/` and mirrored publicly at:

```text
https://popgen.dk/albrecht/open/tutorial_data/admixture/
```

Course material to mine for commands and data choices:

- https://github.com/popgenDK/courses/blob/main/summer2025/exercises/Day3_Admixture_structure_bonus.ipynb
- https://github.com/popgenDK/courses/blob/main/bgi23/Admixture.ipynb
- https://github.com/popgenDK/courses/blob/main/chinaCourse2025/Day4_Morning_admixture_genotype.ipynb

## Set paths

Run all commands from the `admixture/` folder.

```bash
TUTORIAL_DIR=$(pwd)
DATA_ROOT=${TUTORIAL_DIR}/../tutorial_data
DATA_DIR=${DATA_ROOT}/admixture
DATA_URL=https://popgen.dk/albrecht/open/tutorial_data/admixture
RESULTS_DIR=${TUTORIAL_DIR}/results
FIGURES_DIR=${TUTORIAL_DIR}/figures
CODE_DIR=${TUTORIAL_DIR}/code

THREADS=${THREADS:-4}
mkdir -p "${DATA_DIR}" "${RESULTS_DIR}" "${FIGURES_DIR}"

PLINK_PREFIX=${DATA_DIR}/example
PRUNED_PREFIX=${RESULTS_DIR}/example.pruned
```

<details>
<summary>Install software</summary>

ADMIXTURE documentation:

- https://dalexander.github.io/admixture/

PLINK documentation:

- https://www.cog-genomics.org/plink/

Example checks:

```bash
which admixture
admixture --help | head

which plink
plink --help | head
```

</details>

## Prepare PLINK input

Start by checking the input data.

```bash
plink --bfile "${PLINK_PREFIX}" --freq --out "${RESULTS_DIR}/example.freq"
plink --bfile "${PLINK_PREFIX}" --missing --out "${RESULTS_DIR}/example.missing"
```

LD prune before running ADMIXTURE.

```bash
plink \
  --bfile "${PLINK_PREFIX}" \
  --indep-pairwise 50 10 0.1 \
  --out "${RESULTS_DIR}/example.ld"

plink \
  --bfile "${PLINK_PREFIX}" \
  --extract "${RESULTS_DIR}/example.ld.prune.in" \
  --make-bed \
  --out "${PRUNED_PREFIX}"
```

## Run ADMIXTURE

```bash
for K in 2 3 4 5
do
  admixture --cv "${PRUNED_PREFIX}.bed" "${K}" \
    | tee "${RESULTS_DIR}/admixture.K${K}.log"

  mv "example.pruned.${K}.Q" "${RESULTS_DIR}/example.pruned.K${K}.Q"
  mv "example.pruned.${K}.P" "${RESULTS_DIR}/example.pruned.K${K}.P"
done
```

Extract cross-validation errors:

```bash
grep -h "CV error" "${RESULTS_DIR}"/admixture.K*.log \
  > "${RESULTS_DIR}/admixture.cv_errors.txt"
```

## Make figures

![ADMIXTURE K3 ancestry proportions](figures/admixture_k3.png)

Save this as `code/02_plot_admixture.R` and run it from `admixture/`.

```r
data_dir <- "../tutorial_data/admixture"
results_dir <- "results"
figures_dir <- "figures"

q <- read.table(file.path(results_dir, "example.pruned.K3.Q"))

fam_file <- file.path(results_dir, "example.pruned.fam")
if (!file.exists(fam_file)) {
  fam_file <- file.path(data_dir, "example.fam")
}
fam <- read.table(fam_file, as.is = TRUE)

pop <- fam[, 1]
ord <- order(pop, q[, 1])

png(file.path(figures_dir, "admixture_k3.png"), width = 1400, height = 500, res = 160)
par(mar = c(5, 4, 1, 1))
barplot(
  t(q)[, ord],
  col = 2:10,
  space = 0,
  border = NA,
  xlab = "Individuals",
  ylab = "Ancestry proportion"
)
text(sort(tapply(seq_len(nrow(fam)), pop[ord], mean)), -0.05, unique(pop[ord]), xpd = TRUE)
abline(v = cumsum(sapply(unique(pop[ord]), function(x) sum(pop[ord] == x))), col = 1, lwd = 1.2)
dev.off()
```

Optional CV plot:

```r
cv <- readLines(file.path("results", "admixture.cv_errors.txt"))
k <- as.integer(sub(".*K=([0-9]+).*", "\\1", cv))
err <- as.numeric(sub(".*: ", "", cv))

png(file.path("figures", "admixture_cv.png"), width = 800, height = 600, res = 150)
plot(k, err, type = "b", pch = 19, xlab = "K", ylab = "Cross-validation error")
dev.off()
```

## Interpret the result

The tutorial should make users compare `K` values using both CV error and ancestry barplots. It should also mention that ADMIXTURE assumes called genotypes, while NGSadmix is designed for genotype likelihoods and is therefore better suited to low-depth sequencing data.

## Sources

- https://github.com/popgenDK/courses/blob/main/summer2025/exercises/Day3_Admixture_structure_bonus.ipynb
- https://github.com/popgenDK/courses/blob/main/bgi23/Admixture.ipynb
- https://github.com/popgenDK/courses/blob/main/chinaCourse2025/Day4_Morning_admixture_genotype.ipynb
- https://dalexander.github.io/admixture/
- https://www.cog-genomics.org/plink/
