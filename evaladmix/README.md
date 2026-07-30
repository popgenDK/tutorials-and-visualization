# evalAdmix tutorial

This tutorial evaluates admixture model fit using evalAdmix. It is designed to follow either the `ngsadmix/` tutorial or an ADMIXTURE tutorial with compatible `Q`, allele-frequency, and genotype/genotype-likelihood input files.

## Learning goals

- Run evalAdmix on an ancestry model.
- Generate a correlation-of-residuals plot.
- Use residual structure to decide whether a chosen `K` is adequate.
- Compare model fit across `K=3` and `K=4`.

## Input data

This starter tutorial assumes the NGSadmix tutorial produced:

- `../ngsadmix/results/1000G5pops.ngsadmix.K3.seed3.qopt`
- `../ngsadmix/results/1000G5pops.ngsadmix.K3.seed3.fopt.gz`
- `../tutorial_data/ngsadmix/1000G5pops.inputgl.beagle.gz`
- `../tutorial_data/ngsadmix/1000G5pops.pop.info`

Relevant course source:

- https://github.com/popgenDK/courses/blob/main/summer2025/exercises/Day3_Morning_Admixture.ipynb

## Set paths

Run all commands from the `evaladmix/` folder. If using outputs from `ngsadmix/`, point `NGSADMIX_DIR` at that folder.

```bash
TUTORIAL_DIR=$(pwd)
NGSADMIX_DIR=${TUTORIAL_DIR}/../ngsadmix

DATA_ROOT=${TUTORIAL_DIR}/../tutorial_data
DATA_DIR=${DATA_ROOT}/ngsadmix
DATA_URL=https://popgen.dk/albrecht/open/tutorial_data/ngsadmix
NGSADMIX_RESULTS_DIR=${NGSADMIX_DIR}/results
RESULTS_DIR=${TUTORIAL_DIR}/results
FIGURES_DIR=${TUTORIAL_DIR}/figures
CODE_DIR=${TUTORIAL_DIR}/code

THREADS=${THREADS:-4}
mkdir -p "${RESULTS_DIR}" "${FIGURES_DIR}"

BEAGLE=${DATA_DIR}/1000G5pops.inputgl.beagle.gz
POPINFO=${DATA_DIR}/1000G5pops.pop.info
QOPT=${NGSADMIX_RESULTS_DIR}/1000G5pops.ngsadmix.K3.seed3.qopt
FOPT=${NGSADMIX_RESULTS_DIR}/1000G5pops.ngsadmix.K3.seed3.fopt.gz
```

<details>
<summary>Install software</summary>

evalAdmix repository:

- https://github.com/GenisGE/evalAdmix

Example checks:

```bash
which evalAdmix
evalAdmix 2>&1 | head
```

The plotting functions used by the course material are available from:

```r
source("https://raw.githubusercontent.com/GenisGE/evalAdmix/master/visFuns.R")
```

</details>

## Run evalAdmix

```bash
evalAdmix \
  -beagle "${BEAGLE}" \
  -fname "${FOPT}" \
  -qname "${QOPT}" \
  -o "${RESULTS_DIR}/1000G5pops.K3.corres" \
  -P "${THREADS}"
```

## Make figures

![evalAdmix residual correlation plot](figures/evaladmix_k3_residuals.png)

Save this as `code/02_plot_evaladmix.R` and run it from `evaladmix/`.

```r
source("https://raw.githubusercontent.com/GenisGE/evalAdmix/master/visFuns.R")

data_dir <- "../tutorial_data/ngsadmix"
ngsadmix_results_dir <- "../ngsadmix/results"
results_dir <- "results"
figures_dir <- "figures"

pop <- read.table(file.path(data_dir, "1000G5pops.pop.info"), as.is = TRUE)
q <- read.table(file.path(ngsadmix_results_dir, "1000G5pops.ngsadmix.K3.seed3.qopt"))
r <- as.matrix(read.table(file.path(results_dir, "1000G5pops.K3.corres")))

ord <- orderInds(pop = pop[, 1], q = q)

png(file.path(figures_dir, "evaladmix_k3_residuals.png"), width = 900, height = 800, res = 150)
plotCorRes(r, pop = pop[, 1], ord = ord, max_z = 0.2)
dev.off()
```

## Interpret the result

The final tutorial should explain that large blocks of positive or negative residual correlation can indicate structure not captured by the chosen admixture model. Users should compare residual plots across `K` values rather than treating the ancestry barplot alone as proof of a good model.

## Sources

- https://github.com/popgenDK/courses/blob/main/summer2025/exercises/Day3_Morning_Admixture.ipynb
- https://github.com/GenisGE/evalAdmix
