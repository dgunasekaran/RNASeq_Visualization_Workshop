# =============================================================================
# 00_get_data.R — download the GSE135568 count matrix from GEO
#
# Dataset: Nuño-Cabanes et al. (2020) Sci Data 7:69. "A multi-omics dataset of
# heat-shock response in the yeast RNA binding protein Mip6."  GEO: GSE135568
#   Design: 2 genotypes (WT, mip6Δ) x 3 conditions (30 °C = 0 min; 39 °C for
#   20 min; 39 °C for 120 min) x 4 biological replicates = 24 RNA-seq libraries.
#   100 bp paired-end Illumina; TopHat2 -> HTSeq counts on sacCer3.
#
# The repo already ships the same data as data/demo_counts.csv and
# data/demo_metadata.csv (the app's demo files), so this script is only needed
# if you want to fetch the original file yourself (recommended once, so you
# know where the data comes from!). It writes the counts and metadata in the
# app's upload format and checks they match the demo files.
#
# Run from the repo folder (the one with app.R):  source("scripts/00_get_data.R")
# =============================================================================

dir.create("data/raw", showWarnings = FALSE, recursive = TRUE)

geo_file <- "GSE135568_CountMatrixRNAseq.xlsx"
dest <- file.path("data/raw", geo_file)

if (!file.exists(dest) && !file.exists(paste0(dest, ".gz"))) {
  # Option A: GEOquery (handles the GEO FTP layout for you)
  ok <- tryCatch({
    GEOquery::getGEOSuppFiles("GSE135568", baseDir = "data/raw", makeDirectory = FALSE,
                              filter_regex = "CountMatrixRNAseq")
    TRUE
  }, error = function(e) { message("GEOquery failed: ", conditionMessage(e)); FALSE })

  # Option B: direct HTTPS download from the GEO supplementary-file folder
  if (!ok) {
    url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135568/suppl/", geo_file)
    tryCatch(download.file(url, dest, mode = "wb"),
             error = function(e) message("Direct download failed too. Download '", geo_file,
               "' manually from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135568",
               " (Supplementary file section) and put it in data/raw/."))
  }
}

# GEO sometimes gzips supplementary files
if (!file.exists(dest) && file.exists(paste0(dest, ".gz"))) {
  R.utils_ok <- requireNamespace("R.utils", quietly = TRUE)
  if (R.utils_ok) R.utils::gunzip(paste0(dest, ".gz"), remove = FALSE) else
    stop("Please unzip ", dest, ".gz (double-click it) and re-run.")
}

source("scripts/workshop_helpers.R")

# Counts in the app's format: first column gene IDs, one column per sample.
# (read_counts() also trims the leading space in one GEO column name.)
counts <- read_counts(dest)
cat("Read", nrow(counts), "features x", ncol(counts), "samples\n")
stopifnot(ncol(counts) == 24)

# Metadata in the app's format: sample, condition, genotype, time, replicate, batch
metadata <- metadata_from_names(colnames(counts))
print(dplyr::count(metadata, genotype, time))

out <- data.frame(gene_id = rownames(counts), counts, check.names = FALSE)
readr::write_csv(out, "data/GSE135568_counts.csv")
readr::write_csv(metadata, "data/GSE135568_metadata.csv")
cat("Wrote data/GSE135568_counts.csv and data/GSE135568_metadata.csv\n",
    "(upload both in the app, or use them in place of the demo files)\n")

# Same numbers as the demo files that ship with the app?
demo_counts <- read_counts("data/demo_counts.csv")
demo_meta   <- read_metadata("data/demo_metadata.csv", colnames(demo_counts))
same_counts <- identical(dim(demo_counts), dim(counts)) &&
  setequal(colnames(demo_counts), colnames(counts)) &&
  setequal(rownames(demo_counts), rownames(counts)) &&
  all(demo_counts[rownames(counts), colnames(counts)] == counts)
same_meta <- all(demo_meta$condition == metadata$condition[match(demo_meta$sample, metadata$sample)])
cat("Matches demo_counts.csv:", same_counts, " | matches demo_metadata.csv:", same_meta, "\n")
