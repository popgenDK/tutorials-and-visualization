# Software tutorial suggestions

This repository should act as a curated entry point for tutorials, examples, code snippets, and small reproducible data sets for software associated with the Albrechtsen/popgen ecosystem.

Primary discovery sources:

- Albrechtsen Lab software page: https://website.popgen.dk/software/
- Legacy popgen.dk software wiki: https://popgen.dk/software/index.php/Main_Page
- popgenDK GitHub organization: https://github.com/popgenDK
- Anders Albrechtsen GitHub profile: https://github.com/aalbrechtsen
- External software documentation for tools commonly used with these workflows

## Recommended front page structure

The repository front page (`README.md`) should give a method-oriented overview first, then link to one folder per software package.

Suggested sections:

1. Population structure and admixture
2. Low-depth sequencing and genotype likelihood workflows
3. Relatedness and pedigree inference
4. PCA, LD, pruning, and model evaluation
5. Runs of homozygosity, inbreeding, and genome scan workflows
6. Imputation, structural variants, and specialized workflows
7. External software used in lab workflows
8. Utilities and reusable workflow components

Each software folder should contain:

- `README.md`: tutorial overview, installation notes, expected runtime, and links to upstream documentation.
- `code/`: runnable scripts and notebooks.
- `tutorial_data/software-name/`: large input data staged outside git and mirrored to `https://popgen.dk/albrecht/open/tutorial_data/software-name/`.
- `results/`: ignored or lightweight example outputs, depending on size.

Example layout:

```text
software-name/
  README.md
  code/
  results/

tutorial_data/
  software-name/
```

## Priority 1: core tutorials

These should be included first because they are central, widely useful, and good anchors for connected tutorials.

| Software | Suggested folder | Method area | Tutorial angle | Upstream/source |
| --- | --- | --- | --- | --- |
| ANGSD | `angsd/` | Low-depth NGS, genotype likelihoods | From BAM files to genotype likelihoods, allele frequencies, PCA/admixture-ready output | https://www.popgen.dk/angsd/ |
| NGSadmix | `ngsadmix/` | Admixture inference from genotype likelihoods | Run K values, choose K, plot ancestry proportions, compare with ADMIXTURE-style workflows | https://popgen.dk/software/index.php/NgsAdmix |
| PCAngsd | `pcangsd/` | PCA and structure from low-depth NGS | Build covariance matrices, PCA plots, selection scans, admixture mode | https://github.com/Rosemeis/pcangsd |
| fastNGSadmix | `fastngsadmix/` | Single-sample ancestry and PCA | Analyze one low-depth sample against a reference panel; include BAM-to-Beagle path via ANGSD | https://popgen.dk/software/index.php/FastNGSadmix |
| PCAone | `pcaone/` | Large-scale PCA and LD-aware workflows | Scalable PCA, projection, LD corrected for structure, pruning/clumping examples | https://github.com/Zilong-Li/PCAone |
| winSFS | `winsfs/` | Site frequency spectrum inference | Estimate 1D/2D SFS from low-coverage data; compare runtime/memory choices | https://github.com/malthesr/winsfs |
| evalAdmix | `evaladmix/` | Model checking for admixture | Evaluate NGSadmix/ADMIXTURE fit and visualize residual structure | https://website.popgen.dk/software/ |

## Priority 2: relatedness, pedigree, and admixed samples

These are important as a second block because they make the repository more than a population-structure tutorial collection.

| Software | Suggested folder | Method area | Tutorial angle | Upstream/source |
| --- | --- | --- | --- | --- |
| relate | `relate/` | IBD and relatedness mapping | Identify local tracts of relatedness while accounting for LD | https://popgen.dk/software/index.php/Relate |
| relateAdmix | `relateadmix/` | Relatedness in admixed populations | Estimate pairwise relatedness or inbreeding using ancestry proportions and allele frequencies | https://popgen.dk/software/index.php/RelateAdmix |
| NGSremix | `ngsremix/` | Relatedness for admixed NGS samples | Pairwise relatedness from low-depth admixed individuals, including F1/recent admixture examples | https://github.com/KHanghoj/NGSremix |
| IBSrelate | `ibsrelate/` | Close relationships without allele frequencies | Allele-frequency-free inference of close familial relationships | https://website.popgen.dk/software/ |
| APOH | `apoh/` | Hybrid pedigree inference | Infer admixture pedigrees of recent hybrids without a contiguous reference genome | https://github.com/popgenDK/apoh |

## Priority 3: specialized and newer methods

These can be added once the core tutorials have a stable style and data layout.

| Software | Suggested folder | Method area | Tutorial angle | Upstream/source |
| --- | --- | --- | --- | --- |
| QUILT2 | `quilt2/` | Read-aware genotype imputation | Low-coverage, long-read, ancient DNA, and NIPT imputation examples | https://github.com/rwdavies/QUILT |
| SVUPP | `svupp/` | Structural variant genotyping | Genotyping structural variants with pre-phased long reads | https://website.popgen.dk/software/ |
| HaploNet | `haplonet/` | Neural-network haplotype/population structure | Population structure inference from simulations and whole-genome data | https://website.popgen.dk/software/ |
| EMU | `emu/` | PCA with high missingness | EM-PCA for ultra-low coverage or non-overlapping sample data | https://github.com/Rosemeis/emu |
| evalPCA | `evalpca/` | PCA/admixture interpretation | Evaluate whether PCA distances reflect genetic ancestry | https://website.popgen.dk/software/ |
| ASAmap | `asamap/` | Ancestry-specific association mapping | Association mapping when local ancestry is unavailable | https://github.com/e-jorsboe/asaMap |
| SATC | `satc/` | Sex assignment and sex-linked scaffolds | Determine sex and identify sex-linked scaffolds from coverage | https://github.com/popgenDK/SATC |

## Priority 4: helper repositories and workflow components

These are useful, but they should probably be presented as utilities rather than primary tutorial tracks.

| Repository | Suggested folder | Role | Upstream/source |
| --- | --- | --- | --- |
| sites_filters | `sites-filters/` | Reusable genome site filtering scripts for population genomics analyses | https://github.com/popgenDK/sites_filters |
| popgenDK analysis scripts | `analysis-workflows/` | General pipeline examples and reusable plotting/reporting patterns | https://github.com/popgenDK/analysis |
| courses | `courses-crosslinks/` | Cross-links to teaching material, not duplicated tutorial content | https://github.com/popgenDK/courses |

## External software used in workflows

These were not developed by us, but should be included because they are standard tools in workflows around the lab software. They should be clearly labelled as external dependencies or companion methods on the front page.

| Software | Suggested folder | Method area | Tutorial angle | Upstream/source |
| --- | --- | --- | --- | --- |
| ADMIXTURE | `admixture/` | Genotype-based ancestry inference | Prepare PLINK input, LD prune, run multiple K values, cross-validation, plot Q matrices, compare with NGSadmix | https://dalexander.github.io/admixture/ |
| PLINK | `plink/` | Genotype data management and QC | Convert file formats, filter SNPs/samples, LD pruning, PCA helper steps, and prepare inputs for ADMIXTURE/PCAone | https://www.cog-genomics.org/plink/ |
| PLINK ROH detection | `plink-roh/` | Runs of homozygosity | Detect ROH with `--homozyg`, tune SNP/window/length thresholds, summarize `F_ROH`, and plot ROH burden | https://www.cog-genomics.org/plink/1.9/ibd#homozyg |
| ROHan | `rohan/` | ROH and heterozygosity from BAM files | Infer ROH and genome-wide heterozygosity directly from modern or ancient BAMs, including damage-aware ancient DNA settings | https://grenaud.github.io/ROHan/ |
| bcftools/RoH | `bcftools-roh/` | HMM-based ROH inference | Use genotype likelihoods or called genotypes to infer autozygosity tracts; compare with PLINK ROH and ROHan | https://samtools.github.io/bcftools/bcftools.html |
| ADMIXTOOLS | `admixtools/` | f-statistics and admixture graph workflows | Convert inputs, run qp3Pop/qpDstat/qpAdm-style analyses, and connect results to PCA/admixture interpretation | https://github.com/DReichLab/AdmixTools |
| EIGENSOFT/smartpca | `eigensoft/` | PCA and population-genetic utilities | Run smartpca, compare PCA results with PCAone/PCAngsd, and handle EIGENSTRAT conversion | https://github.com/DReichLab/EIG |

## Cross-software workflow suggestions

Some tutorials should intentionally combine in-house and external tools:

- `angsd/` -> `ngsadmix/` -> `evaladmix/`: low-depth genotype likelihoods, ancestry inference, model checking.
- `plink/` -> `admixture/` -> `evaladmix/`: called-genotype ancestry workflow.
- `plink/` -> `pcaone/` -> `evalpca/`: scalable PCA and PCA interpretation.
- `plink-roh/` -> `rohan/`: compare ROH calls from called genotypes and BAM-aware likelihood-based inference.
- `angsd/` -> `pcangsd/` -> `winsfs/`: genotype likelihood population structure and site frequency spectrum inference.

## Suggested first build order

1. Create the front-page `README.md` with the method overview and software index.
2. Add `angsd/`, `ngsadmix/`, and `pcangsd/` because they share low-depth NGS inputs and can reuse example Beagle/BAM data.
3. Add `fastngsadmix/` as a single-sample extension of the same input ecosystem.
4. Add `pcaone/`, `winsfs/`, and `evaladmix/` to cover PCA, SFS, and model evaluation.
5. Add external basics that many workflows depend on: `plink/`, `admixture/`, and `plink-roh/`.
6. Add `rohan/` so ROH tutorials cover both called-genotype and BAM-aware approaches.
7. Add relatedness tutorials: `relate/`, `relateadmix/`, `ngsremix/`, `ibsrelate/`, and `apoh/`.
8. Add newer/specialized tutorials after the common folder template is stable.

## Data policy suggestion

Keep the repository small while allowing substantial tutorial data outside git. For each tutorial:

- Put large tutorial data under `tutorial_data/software-name/`.
- Mirror the same subfolder to `https://popgen.dk/albrecht/open/tutorial_data/software-name/`.
- Publish staged data with `rsync -a -P tutorial_data/ kelly.popgen.dk:/kellyData/home/albrecht/public/open/tutorial_data`.
- Track only lightweight `README.md` files from `tutorial_data/`.
- Put generated outputs under `results/` and decide case by case whether to track them.
- Add checksums for downloaded data where possible.

## Naming policy suggestion

Use lower-case folder names with no punctuation except hyphens where needed. Use the upstream software capitalization in page titles and prose.

Examples:

- Folder: `ngsadmix/`; title: `NGSadmix`
- Folder: `fastngsadmix/`; title: `fastNGSadmix`
- Folder: `relateadmix/`; title: `relateAdmix`
- Folder: `plink-roh/`; title: `PLINK ROH detection`
