# =============================================================================
# install_workshop.R: everything the app AND the workshop scripts need.
# Run ONCE before the workshop (10-25 min on a fresh R install), from the repo
# folder (the one with app.R):   source("scripts/install_workshop.R")
# =============================================================================

source("install.R")    # the app's packages (shiny, DESeq2, edgeR, limma, ...)

options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_extra <- c(
  "stringr", "readr", "forcats",  # tidyverse pieces used by the scripts
  "ggrepel",      # non-overlapping gene labels
  "patchwork",    # multi-panel figures
  "ggdist",       # raincloud / half-eye plots
  "UpSetR",       # set-overlap plots (better than Venn diagrams)
  "R.utils"       # gunzip for the GEO download
)
bioc_extra <- c(
  "GEOquery",      # download from GEO (script 00)
  "org.Sc.sgd.db"  # yeast annotation: systematic ORF -> gene name, GO
)
# Optional: gene lengths (RNA-degradation simulation) and post-workshop
# enrichment analysis (script 04). The workshop runs fine without them.
optional_bioc <- c("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene", "GenomicFeatures",
                   "clusterProfiler", "enrichplot")

missing_cran <- setdiff(cran_extra, rownames(installed.packages()))
if (length(missing_cran)) install.packages(missing_cran)
missing_bioc <- setdiff(bioc_extra, rownames(installed.packages()))
if (length(missing_bioc)) BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
missing_opt <- setdiff(optional_bioc, rownames(installed.packages()))
if (length(missing_opt)) try(BiocManager::install(missing_opt, update = FALSE, ask = FALSE))

# ---- Verify -----------------------------------------------------------------
core <- c(cran_extra, bioc_extra)
ok <- vapply(core, requireNamespace, logical(1), quietly = TRUE)
opt_ok <- vapply(optional_bioc, requireNamespace, logical(1), quietly = TRUE)
if (all(ok)) {
  message("\nWorkshop packages installed. R ", getRversion(),
          ", Bioconductor ", as.character(BiocManager::version()), ".")
} else {
  message("\nThese packages FAILED: ", paste(names(ok)[!ok], collapse = ", "),
          "\nSee the troubleshooting section of the pre-workshop handbook.")
}
if (!all(opt_ok)) message("Optional (not needed on the day): ",
                          paste(names(opt_ok)[!opt_ok], collapse = ", "))

# ---- Smoke test: load the demo data and draw one plot -------------------------
if (all(ok)) {
  source("scripts/workshop_helpers.R")
  counts  <- read_counts("data/demo_counts.csv")
  coldata <- read_metadata("data/demo_metadata.csv", colnames(counts)) |> as_workshop_coldata()
  counts  <- counts[, coldata$sample]
  message("Loaded ", nrow(counts), " features x ", ncol(counts),
          " samples. Test plot shown in the Plots pane. You're ready!")
  print(ggplot(tibble(lib = colSums(counts) / 1e6, s = colnames(counts)), aes(lib, s)) +
          geom_col() + labs(x = "Library size (M)", y = NULL) + theme_pub())
}
