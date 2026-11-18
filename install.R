# Install everything the app needs (run once): source("install.R")
options(repos = c(CRAN = "https://cloud.r-project.org"))
cran <- c("shiny", "bslib", "plotly", "DT", "dplyr", "tidyr", "tibble", "purrr", "ggplot2", "readr",
          "readxl", "pheatmap", "RColorBrewer", "viridisLite", "matrixStats", "scales", "ashr", "BiocManager")
bioc <- c("DESeq2", "edgeR", "limma", "apeglm")
optional <- c("IHW", "org.Sc.sgd.db")   # weighted FDR; yeast gene names
install.packages(setdiff(cran, rownames(installed.packages())))
BiocManager::install(setdiff(bioc, rownames(installed.packages())), update = FALSE, ask = FALSE)
try(BiocManager::install(setdiff(optional, rownames(installed.packages())), update = FALSE, ask = FALSE))
ok <- vapply(c(cran, bioc), requireNamespace, logical(1), quietly = TRUE)
if (all(ok)) message("Ready: shiny::runApp()") else message("Failed: ", paste(names(ok)[!ok], collapse = ", "))
