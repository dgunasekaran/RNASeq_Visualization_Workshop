# =============================================================================
# workshop_helpers.R: functions for the workshop scripts (scripts/00-04).
#
# Data are read with the APP'S OWN readers (read_counts(), read_metadata() from
# R/helpers.R), so the scripts and the app use exactly the same format:
#   counts  = integer matrix, genes x samples, gene IDs as row names
#             (data/demo_counts.csv: first column gene_id, one column per sample)
#   coldata = one row per sample, in the order of colnames(counts)
#             (data/demo_metadata.csv: sample, condition, genotype, time,
#              replicate, batch)
#
# Usage, from the repo folder (the one that contains app.R):
#   source("scripts/workshop_helpers.R")
#   counts  <- read_counts("data/demo_counts.csv")
#   coldata <- read_metadata("data/demo_metadata.csv", colnames(counts)) |>
#     as_workshop_coldata()
#   counts  <- counts[, coldata$sample]
#
# Sections
#   1. Reading counts + metadata (app format)
#   2. Feature classes (protein-coding, mito, rRNA, tRNA, snoRNA, ncRNA)
#   3. Plot theme + palettes for publication-ready figures
#   4. Synthetic QC samples generated from real counts (no read mapping)
#   5. Count-level QC metrics + automatic flags (reference = real samples)
#   6. Normalization strategies
#   7. Differential expression + annotation helpers
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(purrr)
  library(readr)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# 1. Reading counts + metadata (app format)
# -----------------------------------------------------------------------------
if (!file.exists("R/helpers.R")) {
  stop("Run the workshop scripts from the repo folder (the one that contains app.R). ",
       "In RStudio: Session > Set Working Directory, or open the project.", call. = FALSE)
}

# Load the app's helpers into their own environment and take only the readers,
# so the app's QC functions don't clash with the workshop versions below.
.app_helpers <- new.env()
sys.source("R/helpers.R", envir = .app_helpers)
read_counts   <- .app_helpers$read_counts     # gene x sample integer matrix
read_metadata <- .app_helpers$read_metadata   # rows matched to colnames(counts)

condition_levels <- paste(rep(c("WT", "mip6"), each = 3),
                          rep(c("0min", "20min", "120min"), times = 2), sep = "_")

#' Turn the uploaded-style metadata into the factors the scripts model with,
#' add the columns used for synthetic samples, and sort samples by design.
as_workshop_coldata <- function(cd) {
  need <- c("sample", "condition", "genotype", "time", "replicate", "batch")
  miss <- setdiff(need, names(cd))
  if (length(miss)) stop("Metadata is missing column(s): ", paste(miss, collapse = ", "))
  cd <- as_tibble(cd) |>
    mutate(condition = factor(condition, levels = condition_levels),
           genotype  = factor(genotype, levels = c("WT", "mip6")),
           time      = factor(time, levels = c("0min", "20min", "120min")),
           replicate = as.integer(replicate),
           batch     = factor(batch))
  bad <- is.na(cd$condition) | is.na(cd$genotype) | is.na(cd$time)
  if (any(bad)) stop("Unexpected condition/genotype/time values for: ",
                     paste(cd$sample[bad], collapse = ", "))
  if (!"synthetic" %in% names(cd)) cd$synthetic <- FALSE
  if (!"qc_issue" %in% names(cd))  cd$qc_issue  <- "none (real sample)"
  if (!"source" %in% names(cd))    cd$source    <- cd$sample
  arrange(cd, genotype, time, replicate)
}

#' Build demo_metadata.csv-style metadata from GSE135568 sample names,
#' e.g. "mip6_39_20_4" = genotype_temperature_minutes_replicate (script 00).
metadata_from_names <- function(samples) {
  p <- stringr::str_match(samples, "^(wt|mip6)_(\\d+)_(\\d+)_(\\d+)$")
  if (anyNA(p[, 1])) stop("Sample names not in genotype_temp_time_rep format: ",
                          paste(samples[is.na(p[, 1])], collapse = ", "))
  genotype <- ifelse(p[, 2] == "wt", "WT", "mip6")
  time <- paste0(p[, 4], "min")
  tibble(sample    = samples,
         condition = paste(genotype, time, sep = "_"),
         genotype  = genotype,
         time      = time,
         replicate = as.integer(p[, 5]),
         # Data-driven batch: residual correlations (after removing condition
         # means) pair replicates {1,3} and {2,4}; see script 01, part A.
         batch     = ifelse(as.integer(p[, 5]) %in% c(1L, 3L), "B1", "B2"))
}

# -----------------------------------------------------------------------------
# 2. Feature classes
# -----------------------------------------------------------------------------
classify_features <- function(ids) {
  dplyr::case_when(
    grepl("^Y[A-P][LR][0-9]{3}[WC](-[A-Z])?$", ids)          ~ "protein_coding",
    grepl("^Q[0-9]{4}$", ids) | ids %in% c("15S_rRNA", "21S_rRNA") ~ "mitochondrial",
    grepl("^RDN", ids)                                        ~ "rRNA",
    grepl("^t[A-Z]\\(", ids)                                  ~ "tRNA",
    grepl("^snR", ids)                                        ~ "snoRNA",
    TRUE                                                      ~ "other_ncRNA" # SCR1, LSR1, NME1, ...
  )
}

MIP6_ORF <- "YHR015W"   # deleted in mip6Δ: should have ~0 reads -> genotype check

# Reference genes reported as stable in S. cerevisiae RT-qPCR (Teste et al. 2009)
REF_GENES <- c(TAF10 = "YDR167W", TFC1 = "YBR123C", UBC6 = "YER100W", ALG9 = "YNL219C")

# -----------------------------------------------------------------------------
# 3. Plot theme + palettes
# -----------------------------------------------------------------------------
pal_genotype <- c(WT = "#0072B2", mip6 = "#D55E00")              # Okabe-Ito
pal_time     <- c("0min" = "#999999", "20min" = "#E69F00", "120min" = "#CC79A7")
pal_condition <- setNames(c("#9ECAE1", "#4292C6", "#08519C", "#FDAE6B", "#F16913", "#A63603"),
                          condition_levels)
pal_de       <- c(Up = "#B2182B", Down = "#2166AC", NS = "grey75")
pal_status   <- c(real = "grey55", synthetic = "#D55E00")

#' Clean, journal-friendly theme: 8-9 pt text, thin black axes, no grid.
theme_pub <- function(base_size = 9, base_family = "") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.text        = element_text(colour = "black"),
      axis.line        = element_line(linewidth = 0.3),
      axis.ticks       = element_line(colour = "black", linewidth = 0.3),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold"),
      legend.key.size  = unit(3, "mm"),
      legend.title     = element_text(face = "bold"),
      plot.title       = element_text(face = "bold", size = base_size + 1),
      plot.tag         = element_text(face = "bold", size = base_size + 3)
    )
}

#' TRUE if cairo_pdf actually works here. capabilities("cairo") can say TRUE
#' while the device still fails to load (macOS without XQuartz:
#' "failed to load cairo DLL"), so test it once and remember the answer.
cairo_pdf_ok <- local({
  ok <- NULL
  function() {
    if (is.null(ok)) {
      ok <<- isTRUE(tryCatch({
        f <- tempfile(fileext = ".pdf")
        grDevices::cairo_pdf(f)
        grDevices::dev.off()
        unlink(f)
        TRUE
      }, error = function(e) FALSE, warning = function(w) FALSE))
      if (!ok) message("cairo_pdf unavailable (on macOS, install XQuartz); ",
                       "saving PDFs with the standard pdf device instead.")
    }
    ok
  }
})

#' Save at journal sizes: 89 mm = single column, 183 mm = double column.
save_pub <- function(plot, file, width_mm = 89, height_mm = 70, dpi = 600) {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  dev <- NULL
  if (grepl("\\.pdf$", file) && cairo_pdf_ok()) dev <- grDevices::cairo_pdf
  ggsave(file, plot, width = width_mm, height = height_mm, units = "mm",
         dpi = dpi, device = dev)
  invisible(file)
}

log_cpm <- function(counts, prior = 1) {
  log2(t(t(counts) / colSums(counts)) * 1e6 + prior)
}

# -----------------------------------------------------------------------------
# 4. Synthetic QC samples — simulated from REAL counts, no read mapping needed
# -----------------------------------------------------------------------------
# Binomial thinning = sampling each read independently with probability p.
# It is exactly what you would get by sequencing the same library less deeply.
thin_counts <- function(x, p) {
  x <- as.integer(round(x))
  as.integer(stats::rbinom(length(x), size = x, prob = pmin(pmax(p, 0), 1)))
}
thin_matrix <- function(m, p) {
  prob <- if (length(p) == 1) p else rep(p, each = nrow(m))
  out <- matrix(stats::rbinom(length(m), size = as.vector(m), prob = prob),
                nrow = nrow(m), dimnames = dimnames(m))
  storage.mode(out) <- "integer"
  out
}

#' Low sequencing depth.
sim_low_depth <- function(x, target_reads = 0.75e6) {
  thin_counts(x, min(1, target_reads / sum(x)))
}

#' Low-complexity library (too little input RNA, over-amplified by PCR):
#' only `unique_frac` of molecules are captured, then each molecule is copied a
#' highly variable number of times so the library SIZE still looks normal.
#' Sum of u negative-binomial copy numbers = NB(mu = u/f, size = u*k).
sim_low_complexity <- function(x, unique_frac = 0.02, k = 0.3) {
  u <- thin_counts(x, unique_frac)
  out <- integer(length(x))
  pos <- u > 0
  out[pos] <- stats::rnbinom(sum(pos), size = u[pos] * k, mu = u[pos] / unique_frac)
  as.integer(out)
}

#' Genomic DNA contamination: reads land on every feature in proportion to its
#' length (uniform if lengths unknown) -> "silent" genes suddenly look expressed.
sim_gdna <- function(x, frac = 0.15, lengths = NULL) {
  w <- if (is.null(lengths)) rep(1, length(x)) else ifelse(is.na(lengths), stats::median(lengths, na.rm = TRUE), lengths)
  lam <- w / sum(w) * sum(x) * frac / (1 - frac)
  as.integer(x + stats::rpois(length(x), lam))
}

#' Poor poly(A) selection / rRNA depletion: extra reads on non-coding RNAs.
sim_ncrna <- function(x, is_noncoding, frac = 0.35) {
  w <- x * is_noncoding
  lam <- w / sum(w) * sum(x) * frac / (1 - frac)
  as.integer(x + stats::rpois(length(x), lam))
}

#' RNA degradation in a poly(A)-selected library: long transcripts lose more
#' reads (3' bias summed to gene level). Needs gene lengths.
sim_degradation <- function(x, lengths, strength = 1) {
  L <- ifelse(is.na(lengths), stats::median(lengths, na.rm = TRUE), lengths)
  keep <- exp(-strength * L / stats::median(L))
  thin_counts(x, keep / max(keep))
}

#' Irregular replicate: gene-wise technical noise (log-normal, mean-preserving).
sim_noise <- function(x, sdlog = 0.6) {
  mult <- exp(stats::rnorm(length(x), mean = -sdlog^2 / 2, sd = sdlog))
  as.integer(stats::rpois(length(x), x * mult))
}

#' Cross-contamination / mixed sample: frac_x of sample x plus (1 - frac_x) of y.
sim_mixture <- function(x, y, frac_x = 0.7) {
  thin_counts(x, frac_x) + thin_counts(y, 1 - frac_x)
}

#' Optional: exonic gene lengths from the sacCer3 SGD TxDb (if installed).
get_gene_lengths <- function(ids) {
  pkg <- "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene"
  if (!requireNamespace(pkg, quietly = TRUE) ||
      !requireNamespace("GenomicFeatures", quietly = TRUE)) {
    message("Gene lengths unavailable (install ", pkg, "); skipping the degradation sample.")
    return(NULL)
  }
  txdb <- getExportedValue(pkg, pkg)
  ex <- GenomicFeatures::exonsBy(txdb, by = "gene")
  len <- sum(GenomicRanges::width(GenomicRanges::reduce(ex)))
  out <- unname(len[ids])
  names(out) <- ids
  out
}

#' Append a panel of clearly named synthetic samples to the real data.
#' Each synthetic sample is labelled with the group a lab would THINK it is.
make_synthetic_panel <- function(counts, coldata, lengths = NULL, seed = 2026,
                                 which = c("low_depth", "low_complexity", "gDNA", "ncRNA",
                                           "degraded", "noisy", "mixture", "swap")) {
  set.seed(seed)
  fclass <- classify_features(rownames(counts))
  is_nc <- !(fclass %in% c("protein_coding"))
  specs <- tribble(
    ~id,              ~name,                      ~source,          ~label_as,        ~issue,
    "low_depth",      "SYN01_lowDepth",           "wt_30_0_1",      "wt_30_0_1",      "Low depth (0.75 M reads)",
    "low_complexity", "SYN02_lowComplexity",      "wt_30_0_2",      "wt_30_0_2",      "Low-complexity library (2% unique molecules, PCR-amplified)",
    "gDNA",           "SYN03_gDNA",               "wt_39_20_1",     "wt_39_20_1",     "Genomic DNA contamination (15% of reads)",
    "ncRNA",          "SYN04_ncRNA",              "wt_39_20_2",     "wt_39_20_2",     "Poor poly(A) selection: +35% ncRNA/rRNA reads",
    "degraded",       "SYN05_degraded",           "mip6_30_0_1",    "mip6_30_0_1",    "RNA degradation (length-dependent loss)",
    "noisy",          "SYN06_noisyReplicate",     "mip6_39_20_1",   "mip6_39_20_1",   "Irregular replicate: gene-wise technical noise",
    "mixture",        "SYN07_mixture",            "mip6_39_120_1",  "mip6_39_120_1",  "Cross-contamination: 70% mip6 + 30% WT",
    "swap",           "SYN08_swap",               "wt_39_20_3",     "mip6_30_0_3",    "Sample swap: WT 20 min labelled as mip6 0 min"
  ) |> filter(id %in% which)
  if (is.null(lengths)) specs <- filter(specs, id != "degraded")

  new <- lapply(seq_len(nrow(specs)), function(i) {
    s <- specs[i, ]
    x <- counts[, s$source]
    switch(s$id,
      low_depth      = sim_low_depth(x),
      low_complexity = sim_low_complexity(x),
      gDNA           = sim_gdna(x, lengths = lengths),
      ncRNA          = sim_ncrna(x, is_nc),
      degraded       = sim_degradation(x, lengths),
      noisy          = sim_noise(x),
      mixture        = sim_mixture(x, counts[, "wt_39_120_1"]),
      swap           = thin_counts(x, 0.9)
    )
  })
  syn_counts <- do.call(cbind, new)
  colnames(syn_counts) <- specs$name
  rownames(syn_counts) <- rownames(counts)

  # Each synthetic sample copies the metadata row of the sample it pretends to be
  syn_cd <- coldata[match(specs$label_as, coldata$sample), ] |>
    mutate(sample = specs$name, synthetic = TRUE, qc_issue = specs$issue,
           source = specs$source)

  list(counts = cbind(counts, syn_counts),
       coldata = bind_rows(coldata, syn_cd))
}

# -----------------------------------------------------------------------------
# 5. Count-level QC metrics (what you can check WITHOUT BAM files)
# -----------------------------------------------------------------------------
#' Gene-detection saturation: genes detected after thinning to 25% vs 100%.
#' A normal library keeps discovering genes with depth; a low-complexity
#' library is "saturated" early (few new genes with more reads).
saturation_gain <- function(x, fracs = c(0.25, 1), min_count = 1) {
  det <- vapply(fracs, function(f) sum(thin_counts(x, f) >= min_count), numeric(1))
  det[length(det)] - det[1]
}

robust_z <- function(v, ref) {
  m <- stats::median(ref, na.rm = TRUE)
  s <- stats::mad(ref, na.rm = TRUE)
  s <- max(s, 0.01 * abs(m), 1e-6)   # floor avoids division by ~0
  (v - m) / s
}

#' Compute per-sample QC metrics. `coldata$synthetic` marks samples that are NOT
#' used as the reference distribution (so real samples define "normal").
qc_metrics <- function(counts, coldata, marker = MIP6_ORF, z_cut = 3) {
  stopifnot(identical(colnames(counts), coldata$sample))
  fclass <- classify_features(rownames(counts))
  pc <- fclass == "protein_coding"
  if (sum(pc) < 100) pc <- rep(TRUE, nrow(counts))   # non-yeast IDs: use all features
  real <- !coldata$synthetic
  lib <- colSums(counts)
  lc <- log_cpm(counts[pc, , drop = FALSE])
  expressed <- rowSums(counts[pc, real, drop = FALSE] >= 10) == sum(real)
  high <- rowMeans(lc[, real, drop = FALSE]) > log2(50 + 1)
  cors <- stats::cor(lc[expressed, ], method = "pearson")

  per_sample <- lapply(seq_len(ncol(counts)), function(j) {
    s <- coldata$sample[j]
    grp_real <- coldata$sample[real & coldata$condition == coldata$condition[j] & coldata$sample != s]
    others_real <- setdiff(coldata$sample[real], s)
    silent <- rowSums(counts[pc, others_real, drop = FALSE] > 0) == 0
    # nearest group by median correlation to real samples of each group
    grp_cor <- vapply(levels(coldata$condition), function(g) {
      ss <- setdiff(coldata$sample[real & coldata$condition == g], s)
      if (!length(ss)) NA_real_ else stats::median(cors[s, ss])
    }, numeric(1))
    ref_med <- matrixStats::rowMedians(lc[high, grp_real, drop = FALSE])
    tibble(
      sample          = s,
      lib_size_M      = lib[j] / 1e6,
      genes_detected  = sum(counts[pc, j] >= 1),
      genes_ge10      = sum(counts[pc, j] >= 10),
      pct_top50       = 100 * sum(sort(counts[, j], decreasing = TRUE)[1:50]) / lib[j],
      pct_ncRNA       = 100 * sum(counts[fclass %in% c("rRNA", "tRNA", "snoRNA", "other_ncRNA"), j]) / lib[j],
      pct_mito        = 100 * sum(counts[fclass == "mitochondrial", j]) / lib[j],
      ghost_genes     = sum(counts[pc, j][silent] > 0),
      cor_to_group    = stats::median(cors[s, grp_real]),
      nearest_group   = names(which.max(grp_cor)),
      noise_mad       = stats::median(abs(lc[high, j] - ref_med)),
      saturation_gain = saturation_gain(counts[pc, j]),
      marker_cpm      = if (marker %in% rownames(counts)) 1e6 * counts[marker, j] / lib[j] else NA_real_
    )
  })
  m <- bind_rows(per_sample) |> left_join(coldata, by = "sample")

  ref <- m[!m$synthetic, ]
  m |>
    mutate(
      z_cor        = robust_z(cor_to_group, ref$cor_to_group),
      z_noise      = robust_z(noise_mad, ref$noise_mad),
      z_ncRNA      = robust_z(pct_ncRNA, ref$pct_ncRNA),
      z_ghost      = robust_z(ghost_genes, ref$ghost_genes),
      z_saturation = robust_z(saturation_gain, ref$saturation_gain),
      flags = pmap_chr(list(lib_size_M, z_cor, z_noise, z_ncRNA, z_ghost, z_saturation,
                            nearest_group, as.character(condition), as.character(genotype), marker_cpm),
        function(lib, zc, zn, znc, zg, zs, near, grp, gt, mk) {
          ok <- function(v) isTRUE(v)          # NA-safe (e.g. a group with no other replicate)
          f <- c(
            if (lib < 5) "low depth (<5 M)",
            # gDNA also "saturates" detection (every gene gets reads), so only
            # call low complexity when there is no ghost-gene signal
            if (ok(zs < -z_cut) && !ok(zg > z_cut)) "low complexity (saturated)",
            if (ok(zc < -z_cut)) "low concordance with replicates",
            if (ok(zn > z_cut)) "excess technical noise",
            if (ok(znc > z_cut)) "ncRNA/rRNA contamination",
            if (ok(zg > z_cut)) "reads on silent genes (gDNA?)",
            # WT vs mip6 at the same time point are nearly identical, so only a
            # TIME-point mismatch is informative here; genotype is checked with
            # the MIP6 marker gene below.
            if (!is.na(near) && !is.na(grp) && sub(".*_", "", near) != sub(".*_", "", grp))
              paste0("clusters with ", near, " (swap?)"),
            if (ok(gt == "mip6") && ok(mk > 0.5)) "MIP6 reads in mip6Δ (contamination/mislabel)",
            if (ok(gt == "WT") && ok(mk < 0.5)) "no MIP6 reads in WT (mislabel?)"
          )
          if (length(f)) paste(f, collapse = "; ") else "PASS"
        })
    )
}

# -----------------------------------------------------------------------------
# 6. Normalization strategies
# -----------------------------------------------------------------------------
#' Size factors (scaled to geometric mean 1) for several strategies.
#'   TC  = total count (CPM)        UQ  = upper quartile (Bullard 2010)
#'   TMM = trimmed mean of M (Robinson & Oshlack 2010)
#'   MoR = median-of-ratios (Anders & Huber 2010; DESeq2)
size_factors <- function(counts, method = c("TC", "UQ", "TMM", "MoR")) {
  method <- match.arg(method)
  counts <- counts[rowSums(counts) > 0, , drop = FALSE]
  lib <- colSums(counts)
  sf <- switch(method,
    TC  = lib,
    UQ  = lib * edgeR::calcNormFactors(counts, method = "upperquartile"),
    TMM = lib * edgeR::calcNormFactors(counts, method = "TMM"),
    MoR = DESeq2::estimateSizeFactorsForMatrix(counts)
  )
  sf / exp(mean(log(sf)))
}

#' log2 normalized expression for a strategy (incl. quantile + VST).
normalize_counts <- function(counts, method = c("none", "TC", "UQ", "TMM", "MoR", "quantile", "VST")) {
  method <- match.arg(method)
  switch(method,
    none     = log2(counts + 1),
    quantile = limma::normalizeQuantiles(log_cpm(counts)),
    # for a plain matrix DESeq2 returns a matrix (median-of-ratios inside)
    VST      = {
      v <- DESeq2::varianceStabilizingTransformation(as.matrix(counts), blind = TRUE)
      dimnames(v) <- dimnames(counts)
      v
    },
    {
      sf <- size_factors(counts, method)
      log2(t(t(counts) / sf) + 1)
    }
  )
}

#' Relative log expression (RLE): each gene minus its median across samples.
#' Well-normalized data -> boxes centered on 0 with similar spread.
rle_long <- function(logmat) {
  (logmat - matrixStats::rowMedians(logmat)) |>
    as.data.frame() |>
    rownames_to_column("gene") |>
    pivot_longer(-gene, names_to = "sample", values_to = "rle")
}

# -----------------------------------------------------------------------------
# 7. Differential expression + annotation helpers
# -----------------------------------------------------------------------------
#' Filter lowly expressed protein-coding genes (edgeR::filterByExpr logic).
filter_genes <- function(counts, group) {
  pc <- classify_features(rownames(counts)) == "protein_coding"
  keep <- edgeR::filterByExpr(counts, group = group) & pc
  counts[keep, , drop = FALSE]
}

#' DESeq2 with design ~ batch + condition; returns a tidy results table for A vs B.
run_deseq_contrast <- function(counts, coldata, A, B, use_batch = TRUE,
                               shrink = c("none", "ashr", "normal"), alpha = 0.05,
                               dds = NULL) {
  shrink <- match.arg(shrink)
  if (is.null(dds)) {
    cd <- as.data.frame(coldata)
    rownames(cd) <- cd$sample
    cd$condition <- droplevels(cd$condition)
    design <- if (use_batch && nlevels(droplevels(cd$batch)) > 1) ~ batch + condition else ~ condition
    dds <- DESeq2::DESeqDataSetFromMatrix(counts, cd, design)
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
  }
  res <- DESeq2::results(dds, contrast = c("condition", A, B), alpha = alpha)
  if (shrink != "none") {
    res <- DESeq2::lfcShrink(dds, contrast = c("condition", A, B), res = res,
                             type = shrink, quiet = TRUE)
  }
  as.data.frame(res) |>
    rownames_to_column("gene_id") |>
    as_tibble() |>
    mutate(contrast = paste(A, "vs", B))
}

#' Map systematic ORF IDs (YHR015W) to common names (MIP6) with org.Sc.sgd.db.
gene_symbols <- function(ids) {
  sym <- rep(NA_character_, length(ids))
  if (requireNamespace("org.Sc.sgd.db", quietly = TRUE)) {
    sym <- tryCatch(unname(suppressMessages(AnnotationDbi::mapIds(
             org.Sc.sgd.db::org.Sc.sgd.db, keys = ids, keytype = "ORF",
             column = "GENENAME", multiVals = "first"))),
           error = function(e) rep(NA_character_, length(ids)))
  }
  ifelse(is.na(sym) | sym == "", ids, sym)
}

label_de <- function(res, padj_cut = 0.05, lfc_cut = 1) {
  res |>
    mutate(status = case_when(
      !is.na(padj) & padj < padj_cut & log2FoldChange >=  lfc_cut ~ "Up",
      !is.na(padj) & padj < padj_cut & log2FoldChange <= -lfc_cut ~ "Down",
      TRUE ~ "NS"
    ), status = factor(status, levels = c("Up", "Down", "NS")))
}

#' Volcano y-axis. padj = 0 means the p-value underflowed (smaller than the
#' smallest number R can store, ~1e-308), so -log10(padj) would be Inf and the
#' gene would vanish from the plot. Draw those genes one log unit above the most
#' significant gene that has a finite padj, and mark them as capped.
volcano_y <- function(padj) {
  pos <- padj[!is.na(padj) & padj > 0]
  cap <- if (length(pos)) min(pos) / 10 else 1e-300
  -log10(ifelse(!is.na(padj) & padj == 0, cap, padj))
}
capped_note <- function(n, width = 50) {
  if (!n) return(NULL)
  txt <- sprintf("Triangles: %d gene%s with padj = 0 (below numerical precision), drawn at the top",
                 n, if (n > 1) "s" else "")
  paste(strwrap(txt, width), collapse = "\n")   # ggplot does not wrap titles: wrap here
}

#' Publication-style volcano plot with ggrepel labels.
plot_volcano <- function(res, padj_cut = 0.05, lfc_cut = 1, n_label = 15, title = NULL) {
  if (!"symbol" %in% names(res)) res$symbol <- res$gene_id
  res <- label_de(res, padj_cut, lfc_cut) |> filter(!is.na(padj)) |>
    mutate(capped = padj == 0, y = volcano_y(padj))
  lab <- res |> filter(status != "NS") |> arrange(padj, desc(abs(log2FoldChange))) |> head(n_label)
  n_by <- table(res$status)
  legend_lab <- c(Up = paste0("Up (", n_by[["Up"]], ")"),
                  Down = paste0("Down (", n_by[["Down"]], ")"), NS = "NS")
  ggplot(res, aes(log2FoldChange, y, colour = status)) +
    geom_point(data = filter(res, !capped), size = 0.6, alpha = 0.7, stroke = 0) +
    geom_point(data = filter(res, capped), shape = 17, size = 1.4, show.legend = FALSE) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 2, linewidth = 0.25) +
    geom_hline(yintercept = -log10(padj_cut), linetype = 2, linewidth = 0.25) +
    ggrepel::geom_text_repel(data = lab, aes(label = symbol), size = 2.3,
                             colour = "black", max.overlaps = 30, segment.size = 0.2,
                             min.segment.length = 0) +
    scale_colour_manual(values = pal_de, labels = legend_lab, drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    labs(x = expression(log[2]~fold~change), y = expression(-log[10]~adjusted~italic(P)),
         colour = NULL, title = title, caption = capped_note(sum(res$capped))) +
    theme_pub()
}
