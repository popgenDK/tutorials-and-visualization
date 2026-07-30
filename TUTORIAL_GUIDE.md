# Tutorial guide

This document defines the standard format for tutorials in this repository. Each software package should have its own folder with the tutorial markdown, code used in the tutorial, local data download commands, and figures.

## Folder layout

Use this layout for each tutorial:

```text
software-name/
  README.md
  code/
    00_download_data.sh
    01_run_analysis.sh
    02_plot_results.R
  data/
  figures/
    figure-name.svg
    figure-name.png
  results/

tutorial_data/
  software-name/
    README.md
    large-input-file
```

`README.md` is the tutorial page. All commands shown in the tutorial should also exist in scripts under `code/`, so every analysis step and figure can be regenerated.

Large tutorial data should be downloaded into each tutorial's local `data/` folder. The source mirror is staged from the top-level `tutorial_data/` folder and published at:

```text
https://popgen.dk/albrecht/open/tutorial_data/
```

Publish or update the public mirror with:

```bash
rsync -a -P tutorial_data/ kelly.popgen.dk:/kellyData/home/albrecht/public/open/tutorial_data
```

Use this existing tutorial as a style reference:

- https://github.com/aalbrechtsen/fastNGSadmix/blob/master/PROJECTION_TUTORIAL.md

It is a useful model because it builds the software, downloads data, defines only the paths needed for the next command, explains expected outputs, embeds figures, and states what the figures show.

## Required tutorial structure

Each tutorial `README.md` should follow this order:

1. Title and short purpose
2. Learning goals
3. Input data
4. Set up folders and download data
5. Installation inside a folded `<details>` section
6. Run the analysis
7. Explain the program output
8. Inspect key output files
9. Visualizations
10. Code used to make the visualizations
11. Interpret the visualizations
12. Troubleshooting
13. Sources and related tutorials

## Required setup

Every tutorial should start the runnable part with simple local folders. Do not define extra path variables unless they remove real repetition.

Use this pattern:

```bash
DATA_URL=https://popgen.dk/albrecht/open/tutorial_data/software-name

mkdir -p data results figures
```

Each tutorial should include a download block and a matching `code/00_download_data.sh` script:

```bash
DATA_URL=https://popgen.dk/albrecht/open/tutorial_data/software-name

mkdir -p data
wget -nc -P data "${DATA_URL}/input-file-name"
```

Avoid course-server-only paths in the main workflow. Mention them only in notes when they are useful for tracing where an older course exercise came from.

## Explanatory text

Tutorials should explain each major step before showing the command. The reader should understand what the command does and why it is needed before running it.

Use this pattern:

````markdown
We first create the output folders and download the input files. The `-nc` option tells `wget` not to download a file again if it is already present.

```bash
mkdir -p data results figures
wget -nc -P data "${DATA_URL}/input-file-name"
```
````

For analysis commands, explain:

- What input file is being used.
- What output files will be created.
- Why the important parameters are chosen.
- What assumption or limitation matters for interpretation.

Prefer short, local explanations tied to the command immediately below them. For example, the fastNGSadmix projection tutorial first says that the next command estimates admixture proportions for a projected sample, then shows the command, then lists the expected `.qopt` and `.log` files.

## Program output

Every tutorial must explain the important output files before asking users to interpret results. Do not only list filenames. Explain what each file contains, how rows and columns map to samples or variants, and which files are used for downstream plots.

Use a short table where possible:

| File | Created by | What it contains | Used for |
| --- | --- | --- | --- |
| `results/example.K3.Q` | ADMIXTURE | One row per individual and one column per ancestry component | Ancestry barplot |
| `results/example.K3.P` | ADMIXTURE | Allele-frequency estimates for each ancestry component | Model parameters |
| `results/example.K3.log` | ADMIXTURE | Runtime information and likelihood/convergence messages | Checking whether the run behaved as expected |

For NGSadmix-style tutorials, explain at least:

- `.qopt`: ancestry proportions; one row per individual and one column per ancestry component.
- `.fopt.gz`: allele-frequency estimates; one row per site and one column per ancestry component.
- `.log`: runtime information, likelihood, and convergence messages.

For evalAdmix tutorials, explain at least:

- `.corres`: pairwise correlation of residuals between individuals.
- Positive residual correlation means two individuals are more similar than expected under the fitted model.
- Blocks of residual correlation can indicate structure not captured by the chosen `K` or by the reference populations.

The expected-output style should follow the fastNGSadmix projection tutorial:

```markdown
This writes:

- `results/sample.qopt`: estimated ancestry proportions.
- `results/sample.log`: run settings, overlap counts, and convergence information.
```

## Folded installation section

Installation instructions should be present but folded so the tutorial stays readable for users who already have the software.

Use this exact markdown pattern:

````markdown
<details>
<summary>Install software</summary>

```bash
# installation commands here
```

</details>
````

Inside the folded section, include:

- Upstream URL.
- Tested version or commit when known.
- Actual install commands, such as `git clone`, `make`, binary `wget`, or conda/module commands.
- A short check such as `program --help` or `which program`.

## Visualizations and plotting code

Every tutorial must contain visible visualizations in the markdown. Do not leave broken image links in a tutorial.

Use tracked SVG figures for the first draft if the real tutorial data are not yet available. Once the analysis can be run, replace or supplement them with PNG figures generated by `code/02_plot_results.R`.

Each visualization should have:

- A figure file under `figures/`.
- Code under `code/` that creates or updates the figure.
- A caption explaining what is plotted.
- A short paragraph explaining what pattern the reader should look for.
- A note naming the output file used to make the figure.

In the markdown:

```markdown
![Admixture proportions](figures/admixture_k3.svg)

Figure 1. Ancestry proportions for K=3. Each vertical bar is one individual, and colors show inferred ancestry components. The figure is generated from `results/example.K3.Q` by `code/02_plot_results.R`.
```

In `code/02_plot_results.R`:

```r
dir.create("figures", showWarnings = FALSE)

png("figures/admixture_k3.png", width = 1400, height = 500, res = 160)
# plotting code
dev.off()
```

Rules:

- Do not paste screenshots as the only source of a figure.
- Every image linked in the markdown must exist in the repository or be generated by a command shown immediately above it.
- Use stable output filenames.
- Include enough code to regenerate each plot from `data/` and `results/`.
- Keep plot labels short and readable.
- Prefer tracked SVG for draft figures and generated PNG for final analysis figures. PDF can be optional extra output.

For admixture tutorials, include at least:

- An ancestry barplot for each `K` being discussed.
- A convergence summary plot or table across seeds.
- An evalAdmix residual-correlation heatmap for the `K` values being compared.

For projection tutorials, include at least:

- An admixture barplot for the projected sample.
- A PCA projection plot showing the projected sample together with the reference panel.
- Text explaining whether the projected position agrees with the estimated ancestry proportions.

For PCA tutorials, include at least:

- A PCA scatter plot with population/sample labels or colors.
- A plot or table explaining variance explained or eigenvalues, when relevant.

For ROH tutorials, include at least:

- A genome-wide ROH track or segment plot.
- A per-sample summary of total ROH length or `F_ROH`.

## LD pruning

Do not use PLINK LD pruning for data with population structure. PLINK's LD pruning does not correct for structure, so it can remove or retain variants for the wrong reason when ancestry differences create long-range correlations.

For structured data, use PCAone for LD estimation, pruning, and clumping because PCAone can account for population structure. Tutorials that need LD-pruned data for PCA, ADMIXTURE, or related analyses should explain this explicitly before the command.

Use this wording as the default explanation:

```markdown
This data set has population structure, so we do not use PLINK for LD pruning. Standard LD estimates can be confounded by ancestry differences. Instead, we use PCAone, which estimates LD after accounting for population structure, and then use the retained SNPs in the downstream analysis.
```

If a tutorial uses PLINK for file conversion or basic genotype QC, that is fine. The restriction is specifically about LD pruning/clumping in structured data.

If a tutorial demonstrates PLINK LD pruning only as a comparison or as a warning example, label it clearly as inappropriate for structured data and explain why PCAone is preferred.

## Choosing K in admixture analyses

Do not present cross-validation error as the way to choose the meaningful value of `K` in ADMIXTURE, NGSadmix, or related ancestry tutorials.

Cross-validation can be useful as a narrow prediction diagnostic, but it is a bad guide to the biologically meaningful `K` for these tutorials. At higher `K`, optimization can fail to converge reliably across independent runs, so a cross-validation score may reflect a poor local optimum rather than a better model. Even when optimization is stable, the `K` with the best cross-validation score is not necessarily the most interpretable or meaningful description of population structure. A lower or higher `K` can be more useful depending on the question, sampling scheme, reference populations, and the structure visible in the residuals.

Instead, tutorials should:

- Run multiple independent seeds for each `K`.
- Check convergence by comparing log likelihoods across seeds.
- Plot the best runs for several `K` values.
- Use evalAdmix to evaluate model fit from residual correlations.
- Explain which `K` is useful for the biological question, rather than calling one `K` the true or best value.

Use this wording as the default explanation:

```markdown
We do not choose `K` only from ADMIXTURE cross-validation error. Cross-validation measures predictive fit, but the lowest score is not necessarily the most meaningful population-genetic model. Also, at larger `K`, admixture algorithms can fail to converge for some seeds, so cross-validation scores can be affected by local optima. We therefore compare repeated runs, inspect the ancestry plots, and use evalAdmix residual correlations to assess whether important structure remains unexplained.
```

If cross-validation is shown, label it as an optional diagnostic and explicitly state that it is not the criterion for choosing the meaningful `K`.

## Code style

Use shell scripts for command-line workflow and R scripts for plotting unless a tutorial clearly needs another language.

Shell scripts should start with:

```bash
#!/usr/bin/env bash
set -euo pipefail
```

Use variables for paths and thread counts:

```bash
THREADS=${THREADS:-4}
```

Avoid commands that only work on one course server unless clearly marked as course-server-only.

## Data policy

Large data are allowed for the tutorials, but they should not be committed to GitHub.

- Stage large data for upload under `tutorial_data/software-name/`.
- Mirror that folder to `https://popgen.dk/albrecht/open/tutorial_data/software-name/`.
- Sync the mirror from the repository root with `rsync -a -P tutorial_data/ kelly.popgen.dk:/kellyData/home/albrecht/public/open/tutorial_data`.
- Document expected filenames in `tutorial_data/software-name/README.md`.
- In the tutorial, download files into the local tutorial folder with `wget -nc -P data "${DATA_URL}/filename"`.
- Ignore local `data/` folders in git.
- Add checksums when possible.
- Do not commit large input data or generated outputs.

## Source attribution

Each tutorial should list the sources it was adapted from. For the first admixture tutorials, useful sources include:

- `popgenDK/courses/advBinf/exercises/admixture.md`
- `popgenDK/courses/summer2025/exercises/Day3_Morning_Admixture.ipynb`
- `popgenDK/courses/summer2025/exercises/Day3_Admixture_structure_bonus.ipynb`
- `aalbrechtsen/fastNGSadmix/PROJECTION_TUTORIAL.md`
- NGSadmix documentation: https://popgen.dk/software/index.php/NgsAdmix
- evalAdmix repository/documentation: https://github.com/GenisGE/evalAdmix
- ADMIXTURE documentation: https://dalexander.github.io/admixture/

## Checklist before committing a tutorial

- The tutorial starts with path setup.
- Installation is in a folded `<details>` section.
- Major commands are introduced with text explaining what is being done and why.
- Important output files are explained before interpretation.
- Tutorials include visible visualizations with captions and code to generate them.
- LD pruning with population structure uses PCAone, not PLINK.
- Admixture tutorials use evalAdmix and convergence checks to evaluate `K`, not cross-validation as the main criterion.
- Every displayed final figure is generated by code in `code/`; draft SVG figures are allowed until the full analysis can be run.
- Data source and licensing are documented.
- Commands have been run or clearly marked as not yet tested.
- The tutorial has a short interpretation section, not only commands.
