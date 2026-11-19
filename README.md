# RNA-seq QC & Analysis in R: yeast heat-shock workshop

[![DOI](https://zenodo.org/badge/1460704566.svg)](https://doi.org/10.5281/zenodo.24381827)


A 1-hour hands-on workshop (20-min demo + 30-min exercises) on RNA-seq QC, normalization, differential expression and publication-ready figures in R, using GEO GSE135568.

**Participants:** follow this page before the workshop. 
> Pre-workshop handbook · last updated November 16, 2026


## Before you arrive: checklist

Please finish these steps before the workshop. They take 30–45 minutes, mostly waiting for packages to install. The workshop itself is 90 minutes: a 20-minute live demo, and 30 minutes of hands-on R. The rest of the time will be dedicated to troubleshooting and answering questions.

- [ ] R 4.6.x installed (see [Install R and RStudio](#install-r-and-rstudio))
- [ ] RStudio Desktop installed
- [ ] Workshop repository downloaded and opened as an RStudio project
- [ ] `source("scripts/install_workshop.R")` finished with "All core packages installed. You're ready!"
- [ ] Count matrix present: `data/demo_counts.csv` (it ships with the repo; optionally re-download it from GEO)
- [ ] The app opens: `shiny::runApp()`

Bring a laptop with at least 8 GB RAM and 3 GB free disk space, and its charger. If anything fails, send the full error text during the workshop or Q/A session.

## Install R and RStudio

Install R first, then RStudio. The current R release is 4.6, and the Bioconductor release built for it is 3.23 ([Bioconductor install page](https://bioconductor.org/install/)). If you already have R 4.5 or older, update: Bioconductor packages are built for one R version at a time.

| System | R | RStudio Desktop | Notes |
| --- | --- | --- | --- |
| Windows 10/11 | [CRAN › Download R for Windows › base](https://cloud.r-project.org/bin/windows/base/) | [posit.co/download/rstudio-desktop](https://posit.co/download/rstudio-desktop/) | Rtools is not needed: all workshop packages install as binaries |
| macOS (Apple silicon M1–M4) | [CRAN › macOS](https://cloud.r-project.org/bin/macosx/), the **arm64** .pkg | same link, macOS .dmg | Pick arm64 on Apple silicon, x86_64 only on Intel Macs |
| macOS (Intel) | same page, the **x86_64** .pkg | same link | — |
| Ubuntu / Debian | [CRAN Ubuntu instructions](https://cloud.r-project.org/bin/linux/ubuntu/) | same link, .deb | Use [r2u](https://eddelbuettel.github.io/r2u/) or Posit Package Manager for fast binary packages |

Check the install by opening RStudio and typing in the Console:

```r
R.version.string    # should print R version 4.6.x
```


## Get the workshop repository

The repository holds the scripts, the Shiny app and a copy of the count matrix. Repository link: this one (you are here).

1. **Without git:** on the GitHub page click **Code › Download ZIP**, then unzip it somewhere without spaces or accents in the path (for example `Documents/rnaseq-yeast-workshop`).
2. **With git:** in RStudio choose **File › New Project › Version Control › Git**, paste the repository URL and pick a folder.
3. Open the project by double-clicking **`rnaseq-yeast-workshop.Rproj`**. The RStudio title bar should show the project name. Every script assumes this folder is the working directory.

What is inside:

| Path | What it is |
| --- | --- |
| `scripts/install_workshop.R` | Installs and checks every package |
| `data/demo_counts.csv` | The count matrix (7,126 features × 24 samples) |
| `scripts/00_get_data.R` | Re-downloads the counts from GEO |
| `scripts/01_qc_normalization_subsampling.R` | Depth and replicate subsampling, normalization, irregular-replicate experiments |
| `scripts/02_demo_walkthrough.R` | The 20-minute live demo |
| `scripts/03_exercises.R` | The 30-minute hands-on (fill the blanks) |
| `app.R` | The Shiny dashboard |

## Install the R packages

With the project open, run this once in the RStudio Console. It installs about 25 CRAN and Bioconductor packages plus their dependencies and takes 10–25 minutes on a fresh R:

```r
source("scripts/install_workshop.R")
```

If R asks *"Update all/some/none?"*, type `n`. If it asks *"Do you want to install from sources the packages which need compilation?"*, type `no`.

| Group | Packages | Used for |
| --- | --- | --- |
| Wrangling and plots (CRAN) | tidyverse, readxl, ggrepel, patchwork, pheatmap, ggdist, UpSetR, scales, matrixStats | dplyr pipelines, jitter and raincloud plots, multi-panel figures |
| App (CRAN) | shiny, bslib, plotly, DT | the interactive dashboard |
| Statistics (Bioconductor) | DESeq2, edgeR, limma, apeglm, ashr | normalization, differential expression, LFC shrinkage |
| Data and annotation (Bioconductor) | GEOquery, org.Sc.sgd.db | GEO download; ORF → gene name (YHR015W → MIP6) |
| Optional | TxDb.Scerevisiae.UCSC.sacCer3.sgdGene, clusterProfiler, enrichplot | degradation simulation; GO/KEGG after the workshop |

You are done when the Console ends with **"All N core packages installed … You're ready!"** and a bar chart of library sizes appears in the Plots pane. A line listing *optional* packages is fine.

## Download the dataset (GSE135568)

We use the RNA-seq part of [GEO GSE135568](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135568) from [Nuño-Cabanes et al. 2020, *Scientific Data* 7:69](https://doi.org/10.1038/s41597-020-0412-z): *S. cerevisiae* wild type and *mip6*Δ, at 30 °C and after 20 or 120 minutes at 39 °C, 4 biological replicates each (24 libraries, Illumina 100 bp paired-end, HTSeq counts on sacCer3).

The repository already contains the counts as `data/demo_counts.csv`. Downloading the original once is optional. Pick one route:

1. **In R (optional):** `source("scripts/00_get_data.R")`. It uses `GEOquery::getGEOSuppFiles()` to fetch `GSE135568_CountMatrixRNAseq.xlsx` into `data/raw/`, falls back to a direct HTTPS download, checks the 24 columns, and rewrites `data/GSE135568_counts.csv` plus a sample table.
2. **In the browser:** open the GEO page, scroll to **Supplementary file**, download `GSE135568_CountMatrixRNAseq.xlsx` (click *http* in the Download column), and save it in `data/raw/`. Then run step 1, which will find the file.

What to expect in the file:

| Item | Value |
| --- | --- |
| Rows | 7,126 features: 6,685 nuclear protein-coding ORFs plus mitochondrial genes, tRNAs, snoRNAs, rRNAs and other ncRNAs (SCR1, LSR1 …) |
| Columns | 24 samples named `genotype_temperature_minutes_replicate`, e.g. `mip6_39_20_4` = *mip6*Δ, 39 °C, 20 min, replicate 4 |
| Values | Raw integer read counts (11.5–23.1 million counted reads per sample) |
| Quirk | The first sample name has leading spaces; the scripts trim them |

Raw reads (FASTQ) are in SRA under the same GEO series, e.g. run SRR9929263, a WT 0-min library. You do **not** need them: the workshop starts from the count matrix.

## Final check: run the app

In the RStudio Console, with the project open:

```r
source("scripts/workshop_helpers.R")
counts <- read_counts("data/demo_counts.csv")
coldata <- read_metadata("data/demo_metadata.csv", colnames(counts))
dim(counts) # 7126 24
table(coldata$condition) # 4 samples in each of 6 groups
shiny::runApp() # opens the dashboard; first load takes ~20 s
```

In the app, type **HSP104** in the gene search box. It should be highlighted on the volcano plot (Differential expression tab) and plotted in the Gene explorer tab. Close the app with the red **Stop** button in the Console. Optional: run `source("scripts/02_demo_walkthrough.R")` to produce the demo figures in `figures/02_demo/` (about 2 minutes).

## Troubleshooting

| Symptom | Likely cause | Fix |
| --- | --- | --- |
| `cannot open file 'R/helpers.R'` or `Could not find demo_counts.csv` | Working directory is not the repository | Open `rnaseq-yeast-workshop.Rproj`, or *Session › Set Working Directory › To Project Directory* |
| `Bioconductor version '3.xx' requires R version '4.x'` | Old R | Install R 4.6.x, restart RStudio, re-run `install_packages.R` |
| `package 'X' is not available for this version of R` | Same, or a typo | Check `R.version.string`; then `BiocManager::install("X")` |
| `lazy-load database ... is corrupt` | A package was updated while loaded | *Session › Restart R*, then run again |
| Windows: `cannot remove prior installation` / permission denied | Library folder locked or not writable | Close other R sessions; run RStudio once as administrator, or accept the personal library R offers |
| macOS: `xcrun: error` or compiler errors | R tried to build from source | Answer `no` to compiling from source; check you installed the right arm64 / x86_64 R |
| Linux: compilation of `xml2`, `curl` or `ragg` fails | Missing system libraries | `sudo apt install libxml2-dev libcurl4-openssl-dev libssl-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libpng-dev` |
| `GEOquery failed` / timeout in `00_get_data.R` | Network or firewall | Use the browser route, or skip it: the CSV in `data/` is identical |
| App opens but plots are blank | Browser blocks pop-ups or old RStudio viewer | Click *Open in Browser* at the top of the app window |

Still stuck? Send the output of `sessionInfo()` and the full error message to the organizers.

## Optional: 30 minutes of preparation

None of this is required, but it makes the hands-on smoother.

- **R refresher:** the first chapters of [R for Data Science (2e)](https://r4ds.hadley.nz) on data transformation (dplyr) and data visualization (ggplot2).
- **The workflow we follow:** skim the introduction and the "exploratory analysis and visualization" section of Love, Anders, Kim & Huber, [RNA-Seq workflow: gene-level exploratory analysis and differential expression](https://doi.org/10.12688/f1000research.7035.1), *F1000Research* 4:1070 (2015). The [current Bioconductor version](https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html) is kept up to date.
- **A second workflow, with edgeR and limma:** Law CW, Alhamdoosh M, Su S, Dong X, Tian L, Smyth GK, Ritchie ME (2016). [RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR](https://doi.org/10.12688/f1000research.9005.3). *F1000Research* 5:1408 (version 3, 2018). The code runs as the Bioconductor workflow package [RNAseq123](https://bioconductor.org/packages/RNAseq123): `BiocManager::install("RNAseq123")`.
- **The dataset:** read the abstract and the RNA-seq paragraphs of [Nuño-Cabanes et al. 2020](https://doi.org/10.1038/s41597-020-0412-z).
- **One question to think about:** heat shock makes a handful of genes (HSP12, HSP26, HSP104) hundreds of times more abundant. What could that do to a normalization that simply divides by total reads?
