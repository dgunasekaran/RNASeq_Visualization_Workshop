# =============================================================================
# install_packages.R
# Run ONCE before the workshop (10-25 min on a fresh R install).
#   In RStudio: open this file and click "Source", or run
#   source("install_packages.R")
# =============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_pkgs <- c(
  # data wrangling + plotting
  "tidyverse",    # dplyr, tidyr, ggplot2, readr, purrr, stringr, forcats, tibble
  "readxl",       # read the GEO .xlsx count matrix
  "ggrepel",      # non-overlapping gene labels
  "patchwork",    # multi-panel figures
  "pheatmap",     # static, publication-quality heatmaps
  "ggdist",       # raincloud / half-eye plots
  "RColorBrewer", "scales",
  "matrixStats",  # fast row/column statistics
  "UpSetR",       # set-overlap plots (better than Venn diagrams)
  "ashr",         # LFC shrinkage for arbitrary contrasts
  # interactive app
  "shiny", "bslib", "plotly", "DT",
  # utilities
  "BiocManager", "R.utils"
)

bioc_pkgs <- c(
  "DESeq2",        # negative-binomial GLM for differential expression
  "edgeR",         # TMM normalization, filterByExpr
  "limma",         # voom, voomWithQualityWeights, removeBatchEffect
  "GEOquery",      # download from GEO
  "apeglm",        # LFC shrinkage (Zhu, Ibrahim & Love 2019)
  "org.Sc.sgd.db"  # yeast annotation: systematic ORF -> gene name, GO
)

# Optional: gene lengths (for the RNA-degradation simulation) and post-workshop
# enrichment analysis. The workshop runs fine without them.
optional_bioc <- c("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene", "GenomicFeatures",
                   "clusterProfiler", "enrichplot")

missing_cran <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(missing_cran)) install.packages(missing_cran)

missing_bioc <- setdiff(bioc_pkgs, rownames(installed.packages()))
if (length(missing_bioc)) BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)

missing_opt <- setdiff(optional_bioc, rownames(installed.packages()))
if (length(missing_opt)) {
  try(BiocManager::install(missing_opt, update = FALSE, ask = FALSE))
}

# ---- Verify -----------------------------------------------------------------
core <- c(cran_pkgs, bioc_pkgs)
ok <- vapply(core, requireNamespace, logical(1), quietly = TRUE)
opt_ok <- vapply(optional_bioc, requireNamespace, logical(1), quietly = TRUE)
if (all(ok)) {
  message("\nAll ", length(ok), " core packages installed. R ", getRversion(),
          ", Bioconductor ", as.character(BiocManager::version()), ". You're ready!")
} else {
  message("\nThese core packages FAILED: ", paste(names(ok)[!ok], collapse = ", "),
          "\nSee the troubleshooting section of the pre-workshop handbook.")
}
if (!all(opt_ok)) message("Optional (not needed on the day): ",
                          paste(names(opt_ok)[!opt_ok], collapse = ", "))

# Smoke test: load the data and draw one plot
if (all(ok) && file.exists("app/R/helpers.R")) {
  source("app/R/helpers.R")
  d <- load_workshop_data()
  message("Loaded ", nrow(d$counts), " features x ", ncol(d$counts), " samples. Test plot shown in Plots pane.")
  print(ggplot(tibble(lib = colSums(d$counts) / 1e6, s = colnames(d$counts)), aes(lib, s)) +
          geom_col() + labs(x = "Library size (M)", y = NULL) + theme_pub())
}
