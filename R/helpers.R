# =============================================================================
# helpers.R: functions for the RNA-seq QC Explorer app
# (Shiny sources every file in R/ automatically before app.R.)
#
#   1. Palettes & themes (colour-blind safe)
#   2. Reading counts, metadata, annotation
#   3. Feature classes, per-sample QC metrics and flags
#   4. Normalization strategies
#   5. Differential expression with reportable choices
#   6. Plain-language interpretation of QC, PCA, normalization and the model
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# 1. Palettes & themes
# -----------------------------------------------------------------------------
OKABE_ITO    <- c("#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7",
                  "#56B4E9", "#000000", "#F0E442")
GROUP_SHAPES <- c(16, 17, 15, 18, 1, 2, 0, 5, 6, 3, 4, 8)
DE_COLOURS   <- c(Up = "#D55E00", Down = "#0072B2", NS = "#BBBBBB")

group_palette <- function(levels, scheme = "okabe") {
  n <- length(levels)
  cols <- switch(scheme,
    okabe   = if (n <= 8) OKABE_ITO[seq_len(n)] else viridisLite::viridis(n, end = 0.9),
    viridis = viridisLite::viridis(n, end = 0.9),
    cividis = viridisLite::cividis(n, end = 0.95),
    grey    = grDevices::grey.colors(n, start = 0.1, end = 0.65))
  stats::setNames(cols, levels)
}
group_shapes <- function(levels) stats::setNames(rep_len(GROUP_SHAPES, length(levels)), levels)

heat_palette <- function(scheme = "PuOr", n = 101) {
  switch(scheme,
    PuOr = rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "PuOr"))(n)),
    BuOr = grDevices::colorRampPalette(c("#2166AC", "#92C5DE", "#F7F7F7", "#FDB863", "#B35806"))(n))
}

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(axis.text = element_text(colour = "black"),
          axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(colour = "black", linewidth = 0.3),
          strip.background = element_blank(), strip.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"), plot.title = element_text(face = "bold"))
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

save_plot <- function(plot, file, width_mm = 89, height_mm = 75, dpi = 600) {
  dev <- if (grepl("\\.pdf$", file) && cairo_pdf_ok()) grDevices::cairo_pdf else NULL
  ggsave(file, plot, width = width_mm, height = height_mm, units = "mm", dpi = dpi, device = dev)
}

# -----------------------------------------------------------------------------
# 2. Input
# -----------------------------------------------------------------------------
read_table_any <- function(path) {
  ext <- tolower(tools::file_ext(sub("\\.gz$", "", path)))
  df <- switch(ext,
    xlsx = , xls = readxl::read_excel(path),
    tsv = , txt = readr::read_tsv(path, show_col_types = FALSE),
    readr::read_csv(path, show_col_types = FALSE))
  as.data.frame(df, check.names = FALSE)
}

read_counts <- function(path) {
  df <- read_table_any(path)
  ids <- as.character(df[[1]])
  m <- as.matrix(df[, -1, drop = FALSE])
  suppressWarnings(storage.mode(m) <- "numeric")
  colnames(m) <- trimws(colnames(m))
  rownames(m) <- ids
  if (anyDuplicated(ids)) stop("Gene IDs in the first column must be unique.")
  if (any(m < 0, na.rm = TRUE)) stop("Counts cannot be negative.")
  non_int <- any(abs(m - round(m)) > 1e-6, na.rm = TRUE)
  m[is.na(m)] <- 0
  m <- round(m); storage.mode(m) <- "integer"
  attr(m, "non_integer") <- non_int
  m
}

read_metadata <- function(path, samples) {
  cd <- read_table_any(path)
  id_col <- if ("sample" %in% names(cd)) "sample" else names(cd)[1]
  cd[[id_col]] <- trimws(as.character(cd[[id_col]]))
  miss <- setdiff(samples, cd[[id_col]])
  if (length(miss)) stop("These count-matrix columns are missing from the metadata: ",
                         paste(head(miss, 5), collapse = ", "))
  cd <- cd[match(samples, cd[[id_col]]), , drop = FALSE]
  names(cd)[names(cd) == id_col] <- "sample"
  cd
}

#' Factor with syntactically valid levels (needed for DESeq2 coefficient names)
safe_factor <- function(x) {
  x <- as.character(x); lv <- unique(x)
  factor(make.names(lv, unique = TRUE)[match(x, lv)], levels = make.names(lv, unique = TRUE))
}

is_yeast_ids <- function(ids) mean(grepl("^Y[A-P][LR][0-9]{3}[WC]", ids)) > 0.5

annotate_genes <- function(ids, annotation = NULL) {
  sym <- rep(NA_character_, length(ids))
  if (!is.null(annotation) && ncol(annotation) >= 2) {
    sym <- as.character(annotation[[2]])[match(ids, as.character(annotation[[1]]))]
  } else if (is_yeast_ids(ids) && requireNamespace("org.Sc.sgd.db", quietly = TRUE)) {
    sym <- tryCatch(unname(suppressMessages(AnnotationDbi::mapIds(
      org.Sc.sgd.db::org.Sc.sgd.db, keys = ids, keytype = "ORF",
      column = "GENENAME", multiVals = "first"))), error = function(e) sym)
  }
  tibble(gene_id = ids, symbol = ifelse(is.na(sym) | sym == "", ids, sym))
}

# -----------------------------------------------------------------------------
# 3. Feature classes & QC metrics
# -----------------------------------------------------------------------------
#' Yeast: SGD naming classes. Other organisms: mitochondrial genes by symbol
#' (MT-/mt-) when an annotation is given; everything else = "gene".
classify_features <- function(ids, symbols = ids) {
  if (is_yeast_ids(ids)) {
    return(case_when(
      grepl("^Y[A-P][LR][0-9]{3}[WC](-[A-Z])?$", ids)                  ~ "protein_coding",
      grepl("^Q[0-9]{4}$", ids) | ids %in% c("15S_rRNA", "21S_rRNA")   ~ "mitochondrial",
      grepl("^RDN", ids)                                               ~ "rRNA",
      grepl("^t[A-Z]\\(", ids)                                         ~ "tRNA",
      grepl("^snR", ids)                                               ~ "snoRNA",
      TRUE                                                             ~ "other_ncRNA"))
  }
  ifelse(grepl("^(MT|mt|Mt)-", symbols), "mitochondrial", "gene")
}

log_cpm <- function(counts, prior = 1) log2(t(t(counts) / colSums(counts)) * 1e6 + prior)

thin_counts <- function(x, p) as.integer(stats::rbinom(length(x), size = as.integer(x), prob = p))

#' Genes gained between 25% and 100% of reads (a count-level complexity proxy)
saturation_gain <- function(x) sum(x >= 1) - sum(thin_counts(x, 0.25) >= 1)

robust_z <- function(v) {
  m <- stats::median(v, na.rm = TRUE)
  s <- max(stats::mad(v, na.rm = TRUE), 0.01 * abs(m), 1e-6)
  (v - m) / s
}

#' Per-sample QC table with flags. All flags are RELATIVE to the other samples.
qc_metrics <- function(counts, group, fclass, marker = NULL, marker_absent_in = NULL,
                       z_cut = 3, swap_margin = 0.015) {
  set.seed(1)
  lib <- colSums(counts)
  gene_rows <- fclass %in% c("protein_coding", "gene")
  if (sum(gene_rows) < 100) gene_rows <- rep(TRUE, nrow(counts))
  g <- counts[gene_rows, , drop = FALSE]
  lc <- log_cpm(g)
  expressed <- rowSums(g >= 10) == ncol(g)
  if (sum(expressed) < 50) expressed <- rowMeans(g) >= 10
  cors <- stats::cor(lc[expressed, , drop = FALSE])
  high <- rowMeans(lc) > log2(51)
  if (sum(high) < 50) high <- rowMeans(lc) >= stats::quantile(rowMeans(lc), 0.75)
  group <- as.character(group)
  noncoding <- fclass %in% c("rRNA", "tRNA", "snoRNA", "other_ncRNA")
  mito <- fclass == "mitochondrial"

  rows <- lapply(seq_len(ncol(g)), function(j) {
    mates <- setdiff(which(group == group[j]), j)
    grp_cor <- vapply(unique(group), function(G) {
      ss <- setdiff(which(group == G), j)
      if (!length(ss)) NA_real_ else stats::median(cors[j, ss])
    }, numeric(1))
    own <- grp_cor[group[j]]
    others <- grp_cor[names(grp_cor) != group[j]]
    nearest <- if (all(is.na(others))) NA_character_ else names(which.max(others))
    silent <- rowSums(g[, -j, drop = FALSE] > 0) == 0
    tibble(
      sample          = colnames(g)[j],
      group           = group[j],
      reads_M         = lib[j] / 1e6,
      genes_detected  = sum(g[, j] >= 1),
      genes_10plus    = sum(g[, j] >= 10),
      pct_top50       = 100 * sum(sort(counts[, j], decreasing = TRUE)[1:min(50, nrow(counts))]) / lib[j],
      pct_noncoding   = if (any(noncoding)) 100 * sum(counts[noncoding, j]) / lib[j] else NA_real_,
      pct_mito        = if (any(mito)) 100 * sum(counts[mito, j]) / lib[j] else NA_real_,
      reads_on_silent_genes = sum(g[silent, j] > 0),
      r_to_replicates = if (length(mates)) unname(own) else NA_real_,
      nearest_other_group = nearest,
      r_margin        = if (length(mates) && !is.na(nearest)) unname(own - others[nearest]) else NA_real_,
      noise_mad       = if (length(mates)) stats::median(abs(lc[high, j] -
                          matrixStats::rowMedians(lc[high, mates, drop = FALSE]))) else NA_real_,
      saturation_gain = saturation_gain(g[, j]),
      marker_cpm      = if (!is.null(marker) && marker %in% rownames(counts)) 1e6 * counts[marker, j] / lib[j] else NA_real_)
  })
  q <- bind_rows(rows)

  z <- list(r = robust_z(q$r_to_replicates), noise = robust_z(q$noise_mad),
            nc = robust_z(q$pct_noncoding), mito = robust_z(q$pct_mito),
            ghost = robust_z(q$reads_on_silent_genes), sat = robust_z(q$saturation_gain))
  ok <- function(v) !is.na(v) & v
  med_lib <- stats::median(q$reads_M)
  q$flags <- vapply(seq_len(nrow(q)), function(i) {
    f <- c(
      if (q$reads_M[i] < 0.25 * med_lib) "very low depth",
      if (ok(z$sat[i] < -z_cut) && !ok(z$ghost[i] > z_cut)) "low complexity",
      if (ok(z$r[i] < -z_cut)) "low agreement with replicates",
      if (ok(z$noise[i] > z_cut)) "excess noise",
      if (ok(z$nc[i] > z_cut)) "high non-coding RNA",
      if (ok(z$mito[i] > z_cut)) "high mitochondrial",
      if (ok(z$ghost[i] > z_cut)) "reads on silent genes",
      if (ok(q$r_margin[i] < -swap_margin)) paste0("more similar to ", q$nearest_other_group[i]),
      if (!is.null(marker_absent_in) && !is.na(q$marker_cpm[i]) && q$group[i] %in% marker_absent_in && q$marker_cpm[i] > 0.5)
        "marker gene expressed where it should be absent",
      if (!is.null(marker_absent_in) && !is.na(q$marker_cpm[i]) && !(q$group[i] %in% marker_absent_in) && q$marker_cpm[i] < 0.5)
        "marker gene missing where expected")
    if (length(f)) paste(f, collapse = "; ") else "PASS"
  }, character(1))
  q
}

# -----------------------------------------------------------------------------
# 4. Normalization
# -----------------------------------------------------------------------------
size_factors <- function(counts, method = c("TC", "UQ", "TMM", "MoR")) {
  method <- match.arg(method)
  counts <- counts[rowSums(counts) > 0, , drop = FALSE]
  lib <- colSums(counts)
  sf <- switch(method,
    TC  = lib,
    UQ  = lib * edgeR::calcNormFactors(counts, method = "upperquartile"),
    TMM = lib * edgeR::calcNormFactors(counts, method = "TMM"),
    MoR = DESeq2::estimateSizeFactorsForMatrix(counts))
  sf / exp(mean(log(sf)))
}

normalize_counts <- function(counts, method) {
  switch(method,
    none     = log2(counts + 1),
    quantile = limma::normalizeQuantiles(log_cpm(counts)),
    VST      = { v <- DESeq2::varianceStabilizingTransformation(as.matrix(counts), blind = TRUE)
                 dimnames(v) <- dimnames(counts); v },
    { sf <- size_factors(counts, method); log2(t(t(counts) / sf) + 1) })
}

rle_long <- function(m) {
  (m - matrixStats::rowMedians(m)) |> as.data.frame() |> rownames_to_column("gene") |>
    pivot_longer(-gene, names_to = "sample", values_to = "rle")
}

norm_diagnostics <- function(norm_list, group, A, B) {
  imap_dfr(norm_list, function(m, nm) tibble(
    method = nm,
    median_log2FC_A_vs_B = stats::median(rowMeans(m[, group == A, drop = FALSE]) -
                                         rowMeans(m[, group == B, drop = FALSE])),
    RLE_median_IQR = stats::median(apply(m - matrixStats::rowMedians(m), 2, stats::IQR)),
    within_group_SD = stats::median(map_dbl(unique(group), \(lv)
                        stats::median(matrixStats::rowSds(m[, group == lv, drop = FALSE])))),
    max_abs_RLE_median = max(abs(apply(m - matrixStats::rowMedians(m), 2, stats::median)))))
}

# -----------------------------------------------------------------------------
# 5. Differential expression
# -----------------------------------------------------------------------------
run_de <- function(dds, A, B, alpha = 0.05, lfc = 1, shrink = "none", padj_method = "BH",
                   threshold_mode = "posthoc", indep_filter = TRUE, cooks = TRUE) {
  notes <- character(0)
  if (shrink == "apeglm" && !requireNamespace("apeglm", quietly = TRUE)) { shrink <- "ashr"; notes <- c(notes, "apeglm not installed: used ashr") }
  if (shrink == "ashr" && !requireNamespace("ashr", quietly = TRUE)) { shrink <- "normal"; notes <- c(notes, "ashr not installed: used normal prior") }
  if (padj_method == "IHW" && !requireNamespace("IHW", quietly = TRUE)) { padj_method <- "BH"; notes <- c(notes, "IHW not installed: used BH") }
  args <- list(alpha = alpha, lfcThreshold = if (threshold_mode == "tested") lfc else 0,
               altHypothesis = "greaterAbs", independentFiltering = indep_filter)
  if (!cooks) args$cooksCutoff <- FALSE
  if (padj_method == "IHW") args$filterFun <- IHW::ihw else args$pAdjustMethod <- padj_method

  if (shrink == "apeglm") {
    dds$group <- stats::relevel(dds$group, ref = B)
    dds <- DESeq2::nbinomWaldTest(dds, quiet = TRUE)
    coef <- paste0("group_", A, "_vs_", B)
    res <- do.call(DESeq2::results, c(list(object = dds, name = coef), args))
    res_s <- DESeq2::lfcShrink(dds, coef = coef, res = res, type = "apeglm", quiet = TRUE)
  } else {
    res <- do.call(DESeq2::results, c(list(object = dds, contrast = c("group", A, B)), args))
    res_s <- if (shrink == "none") res else
      DESeq2::lfcShrink(dds, contrast = c("group", A, B), res = res, type = shrink, quiet = TRUE)
  }
  tab <- tibble(gene_id = rownames(res), baseMean = res$baseMean, log2FC_MLE = res$log2FoldChange,
                lfcSE = res$lfcSE, log2FoldChange = res_s$log2FoldChange,
                stat = res$stat, pvalue = res$pvalue, padj = res$padj)
  sig <- !is.na(tab$padj) & tab$padj < alpha
  dir_lfc <- if (threshold_mode == "tested") tab$log2FC_MLE else tab$log2FoldChange
  cut <- if (threshold_mode == "tested") 0 else lfc
  tab$status <- factor(ifelse(sig & dir_lfc >= cut & dir_lfc > 0, "Up",
                       ifelse(sig & dir_lfc <= -cut & dir_lfc < 0, "Down", "NS")),
                       levels = c("Up", "Down", "NS"))
  meta <- list(n_input = nrow(dds), n_tested = sum(!is.na(res$pvalue)),
               n_zero = sum(res$baseMean == 0), n_cooks = sum(is.na(res$pvalue) & res$baseMean > 0),
               n_indep = sum(!is.na(res$pvalue) & is.na(res$padj)),
               filter_thresh = tryCatch(unname(S4Vectors::metadata(res)$filterThreshold), error = function(e) NA),
               n_up = sum(tab$status == "Up"), n_down = sum(tab$status == "Down"), notes = notes)
  list(res = tab, meta = meta,
       params = list(A = A, B = B, alpha = alpha, lfc = lfc, shrink = shrink, padj_method = padj_method,
                     threshold_mode = threshold_mode, indep_filter = indep_filter, cooks = cooks,
                     n_A = sum(dds$group == A), n_B = sum(dds$group == B)))
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

#' DT column renderer for p-values: 3 significant digits, and "< 1e-300" for 0.
#' Sorting and filtering still use the underlying numbers.
pvalue_render <- function() DT::JS(
  "function(d, type) { if (type !== 'display' || d === null) return d;",
  "  return d === 0 ? '&lt; 1e-300' : Number(d).toPrecision(3); }")

de_grid <- function(res, alphas = c(0.01, 0.05, 0.1), lfcs = c(0, 0.585, 1, 2)) {
  expand_grid(padj = alphas, lfc = lfcs) |>
    mutate(Up   = map2_int(padj, lfc, \(a, l) sum(!is.na(res$padj) & res$padj < a & res$log2FoldChange >=  l)),
           Down = map2_int(padj, lfc, \(a, l) sum(!is.na(res$padj) & res$padj < a & res$log2FoldChange <= -l)),
           FC = round(2^lfc, 2)) |>
    transmute(`padj <` = padj, `|log2FC| ≥` = lfc, `fold change ≥` = FC, Up, Down)
}

SHRINK_TEXT <- c(
  apeglm = "apeglm (heavy-tailed Cauchy prior; Zhu, Ibrahim & Love 2019)",
  ashr   = "ashr (adaptive shrinkage; Stephens 2017)",
  normal = "the DESeq2 normal prior (Love, Huber & Anders 2014)",
  none   = "maximum-likelihood estimates without shrinkage")
PADJ_TEXT <- c(
  BH = "Benjamini–Hochberg false discovery rate", IHW = "independent hypothesis weighting (Ignatiadis et al. 2016)",
  BY = "Benjamini–Yekutieli false discovery rate", bonferroni = "Bonferroni (family-wise error)",
  holm = "Holm (family-wise error)")

methods_text <- function(p, meta, design_txt, n_genes_total, prefilter_txt) {
  thr <- if (p$threshold_mode == "tested")
    sprintf("Wald tests were performed against |log2 fold change| > %s (lfcThreshold, altHypothesis = \"greaterAbs\")", p$lfc)
  else sprintf("Wald tests were performed against a log2 fold change of 0 and genes were additionally required to have |log2 fold change| ≥ %s", p$lfc)
  paste0(
    "Differential expression was analysed with DESeq2 (v", as.character(utils::packageVersion("DESeq2")),
    "; Love, Huber & Anders 2014) in R. ", prefilter_txt, " (", format(meta$n_input, big.mark = ","), " of ",
    format(n_genes_total, big.mark = ","), " features kept). Size factors were estimated with the median-of-ratios method ",
    "and gene-wise dispersions were shrunk towards a fitted mean–dispersion trend. The model was ", design_txt, ". ",
    p$A, " (n = ", p$n_A, ") was compared with ", p$B, " (n = ", p$n_B, "). ", thr, ". ",
    "P-values were adjusted with the ", PADJ_TEXT[[p$padj_method]], " procedure and genes with adjusted p < ", p$alpha,
    " were considered significant. ",
    if (p$indep_filter && p$padj_method != "IHW") "Independent filtering on mean normalized count was applied. " else "",
    if (p$cooks) "Count outliers were flagged by Cook's distance. " else "Cook's distance filtering was disabled. ",
    "Reported log2 fold changes use ", SHRINK_TEXT[[p$shrink]], ". ",
    "This identified ", meta$n_up, " up-regulated and ", meta$n_down, " down-regulated genes among ",
    format(meta$n_tested, big.mark = ","), " tested.")
}

# -----------------------------------------------------------------------------
# 6. Interpretation (returns markdown)
# -----------------------------------------------------------------------------
FLAG_ADVICE <- c(
  "very low depth" = "DESeq2 corrects for depth, but a library with < 25% of the median depth adds noise. Keep it if it agrees with its replicates; otherwise consider excluding it.",
  "low complexity" = "Few new genes appear as reads increase: often low input or PCR over-amplification. Check duplication rates (Picard/preseq) and consider excluding or down-weighting.",
  "low agreement with replicates" = "Look at it on the PCA: next to another group suggests a swap or mislabel; alone suggests a failed library or a biological outlier.",
  "excess noise" = "Gene-level variation against its replicates is unusually high. It often coincides with low depth or degraded RNA (check RIN).",
  "high non-coding RNA" = "Reads are spent on rRNA/tRNA/ncRNA: incomplete poly(A) selection or rRNA depletion. Fewer informative reads but usually still usable.",
  "high mitochondrial" = "A high mitochondrial share can mean cell stress or lysis during preparation. Compare with the biology you expect.",
  "reads on silent genes" = "Reads on genes silent in every other sample suggest genomic DNA contamination. Check DNase treatment and intronic/intergenic read rates.",
  "more similar to" = "This sample correlates more with another group than with its own: verify the sample sheet, labels and marker genes before analysis.",
  "marker gene" = "The marker gene contradicts the label: possible contamination, mislabel or a failed genotype.")

advice_for <- function(flag) {
  key <- names(FLAG_ADVICE)[vapply(names(FLAG_ADVICE), function(k) startsWith(flag, k), logical(1))]
  if (length(key)) FLAG_ADVICE[[key[1]]] else ""
}

fmt <- function(x, d = 1) formatC(x, format = "f", digits = d, big.mark = ",")

qc_interpretation <- function(q) {
  n <- nrow(q); lib <- q$reads_M
  out <- c(
    sprintf("**Depth.** Libraries hold %s–%s million reads (median %s; largest/smallest = %s×). %s",
            fmt(min(lib)), fmt(max(lib)), fmt(stats::median(lib)), fmt(max(lib) / min(lib)),
            if (max(lib) / min(lib) > 3) "A spread above 3× is worth noting: normalization corrects it, but the shallowest samples carry more noise." else "This spread is easily handled by normalization."),
    sprintf("**Detection.** %s–%s features have ≥ 10 reads per sample. Samples well below the others detect fewer genes, usually because of depth or complexity.",
            fmt(min(q$genes_10plus), 0), fmt(max(q$genes_10plus), 0)))
  rr <- q |> filter(!is.na(r_to_replicates)) |> group_by(group) |>
    summarise(med = stats::median(r_to_replicates), low = min(r_to_replicates),
              who = sample[which.min(r_to_replicates)], .groups = "drop")
  if (nrow(rr)) out <- c(out, sprintf(
    "**Replicate agreement.** Median correlation with same-group replicates is %s–%s across groups; the lowest single sample is %s (r = %s). Values above ~0.95 on log-CPM of expressed genes are typical for well-controlled experiments.",
    fmt(min(rr$med), 3), fmt(max(rr$med), 3), rr$who[which.min(rr$low)], fmt(min(rr$low), 3)))
  if (any(!is.na(q$r_margin)) && mean(q$r_margin < 0, na.rm = TRUE) > 0.3)
    out <- c(out, "**Similar groups.** Many samples correlate about as well with another group as with their own. Those groups differ only slightly, so expect few or small differences between them.")
  if (any(!is.na(q$pct_noncoding))) out <- c(out, sprintf("**Non-coding share.** %s–%s%% of reads fall on rRNA/tRNA/snoRNA/ncRNA features.", fmt(min(q$pct_noncoding)), fmt(max(q$pct_noncoding))))
  if (any(!is.na(q$pct_mito))) out <- c(out, sprintf("**Mitochondrial share.** %s–%s%% of reads.", fmt(min(q$pct_mito)), fmt(max(q$pct_mito))))
  fl <- q |> filter(flags != "PASS")
  if (!nrow(fl)) {
    out <- c(out, sprintf("**Flags.** None of the %d samples stands out on any metric (robust z-score beyond ±3 against the other samples).", n))
  } else {
    items <- vapply(seq_len(nrow(fl)), function(i) {
      fs <- strsplit(fl$flags[i], "; ")[[1]]
      sprintf("- **%s** (%s): %s. %s", fl$sample[i], fl$group[i], fl$flags[i],
              paste(unique(vapply(fs, advice_for, character(1))), collapse = " "))
    }, character(1))
    out <- c(out, sprintf("**Flags.** %d of %d samples stand out:", nrow(fl), n), items,
             "Flags are relative to your other samples. Decide on removal before looking at DE results, and report what you removed and why. You can exclude samples in the sidebar.")
  }
  if (n < 6) out <- c(out, "_With fewer than 6 samples the robust z-scores are unstable: treat flags as prompts to look, not verdicts._")
  paste(out, collapse = "\n\n")
}

pca_associations <- function(pc, coldata, factors, k = 4) {
  k <- min(k, ncol(pc$x))
  expand_grid(PC = paste0("PC", seq_len(k)), factor = factors) |>
    mutate(R2 = map2_dbl(PC, factor, \(p, f) {
      y <- pc$x[, p]; x <- factor(coldata[[f]])
      if (nlevels(x) < 2 || nlevels(x) >= length(y)) return(NA_real_)
      summary(stats::lm(y ~ x))$r.squared
    }))
}

pca_interpretation <- function(pc, assoc, coldata) {
  ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  best <- assoc |> filter(!is.na(R2)) |> group_by(PC) |> slice_max(R2, n = 1, with_ties = FALSE) |> ungroup()
  out <- sprintf("**Main axes.** PC1 explains %s%%, PC2 %s%% and PC3 %s%% of the variance among the top variable genes.",
                 fmt(ve[1]), fmt(ve[2]), fmt(ve[3]))
  lines <- vapply(seq_len(nrow(best)), function(i)
    sprintf("- %s (%s%%) is best explained by **%s** (R² = %s).", best$PC[i],
            fmt(ve[as.integer(sub("PC", "", best$PC[i]))]), best$factor[i], fmt(best$R2[i], 2)), character(1))
  out <- c(out, lines)
  g <- assoc |> filter(factor == "group", PC %in% c("PC1", "PC2"))
  if (nrow(g) && max(g$R2, na.rm = TRUE) < 0.3)
    out <- c(out, "**Weak group signal.** Your groups explain < 30% of PC1 and PC2. Other factors or noise dominate, so expect fewer or smaller differences.")
  if ("covariate" %in% assoc$factor) {
    cv <- assoc |> filter(factor == "covariate", PC %in% c("PC1", "PC2", "PC3"))
    if (nrow(cv) && max(cv$R2, na.rm = TRUE) > 0.3)
      out <- c(out, sprintf("**Covariate effect.** The covariate explains up to %s%% of a main axis. Keep it in the design (*~ covariate + group*); use covariate removal only for plots.",
                            fmt(100 * max(cv$R2, na.rm = TRUE), 0)))
  }
  # samples far from their group centroid on PC1-2
  xy <- pc$x[, 1:min(2, ncol(pc$x)), drop = FALSE]
  d <- vapply(seq_len(nrow(xy)), function(i) {
    same <- which(coldata$group == coldata$group[i])
    sqrt(sum((xy[i, ] - colMeans(xy[same, , drop = FALSE]))^2))
  }, numeric(1))
  far <- rownames(xy)[robust_z(d) > 3]
  if (length(far)) out <- c(out, sprintf("**Possible outliers on PC1–PC2:** %s sit far from their group centre. Cross-check with the QC tab before excluding.", paste(far, collapse = ", ")))
  paste(out, collapse = "\n\n")
}

norm_interpretation <- function(nd, A, B) {
  get <- function(m) nd$median_log2FC_A_vs_B[nd$method == m]
  out <- sprintf("**Assumption check (%s vs %s).** If most genes do not change, the median log2 fold change should be close to 0.", A, B)
  if (all(c("TC", "MoR") %in% nd$method)) {
    tc <- get("TC"); mor <- get("MoR")
    out <- c(out, sprintf("- Total count (CPM): median log2FC = %+.2f; median-of-ratios: %+.2f%s.",
                          tc, mor, if ("TMM" %in% nd$method) sprintf("; TMM: %+.2f", get("TMM")) else ""))
    out <- c(out, if (abs(tc - mor) > 0.1)
      "- **Composition bias detected.** A few highly changed genes take up a large share of reads, which makes other genes look shifted under total-count scaling. Use median-of-ratios (DESeq2) or TMM (edgeR), as the DE tab does."
    else "- Total-count and median-of-ratios agree, so composition bias is small here.")
  }
  best <- nd |> filter(method %in% c("TC", "UQ", "TMM", "MoR")) |> arrange(max_abs_RLE_median)
  if (nrow(best)) out <- c(out, sprintf("- Most centred RLE boxes: **%s** (largest sample offset %s log2 units). Boxes should sit on 0 with similar widths; a wide box marks a noisy sample.",
                                       best$method[1], fmt(best$max_abs_RLE_median[1], 2)))
  out <- c(out, "- Remember: every method assumes most genes do not change. If total mRNA per cell changes between groups (e.g. global transcription shutdown), only spike-in controls can reveal it.")
  paste(out, collapse = "\n")
}

model_interpretation <- function(x, dds, design_txt, cov_dropped, crosstab_ok, pca_cov_r2 = NA) {
  p <- x$params; m <- x$meta; r <- x$res
  bcv <- sqrt(stats::median(S4Vectors::mcols(dds)$dispersion, na.rm = TRUE))
  sf <- DESeq2::sizeFactors(dds)
  pv <- r$pvalue[!is.na(r$pvalue) & r$baseMean > 1]
  pi0 <- if (length(pv)) min(1, 2 * mean(pv > 0.5)) else NA
  hi_tail <- mean(pv > 0.9); mid <- mean(pv > 0.5 & pv <= 0.6)
  out <- c(
    sprintf("**Design.** %s. %s (n = %d) vs %s (n = %d).%s", design_txt, p$A, p$n_A, p$B, p$n_B,
            if (min(p$n_A, p$n_B) < 3) " **With fewer than 3 replicates per group, dispersion estimates and p-values are unreliable.**"
            else if (min(p$n_A, p$n_B) < 6) " Power is limited: in a 48-replicate yeast study, 3 replicates found only 20–40% of DE genes and ≥ 6 per group is recommended (Schurch et al. 2016)." else ""),
    if (cov_dropped) "**Covariate dropped.** It is confounded with the groups (every level of one maps onto the other), so its effect cannot be separated from the group effect. Redesign or interpret with caution." else NULL,
    if (!is.na(pca_cov_r2) && pca_cov_r2 > 0.3 && !grepl("+", design_txt, fixed = TRUE)) "**Consider adding the covariate.** It explains a main PCA axis but is not in the model." else NULL,
    sprintf("**Normalization.** Median-of-ratios size factors range %s–%s (1 = average depth).", fmt(min(sf), 2), fmt(max(sf), 2)),
    sprintf("**Biological variability.** The typical biological coefficient of variation is %s. Around 0.1 is common for genetically identical lab organisms; 0.2–0.4 or more is typical for human or field samples (edgeR User's Guide). Higher variation needs more replicates.", fmt(bcv, 2)),
    sprintf("**Filtering.** %s features were modelled; %s had zero counts, %s were flagged as count outliers (Cook's distance) and %s were removed by independent filtering%s.",
            fmt(m$n_input, 0), fmt(m$n_zero, 0), fmt(m$n_cooks, 0), fmt(m$n_indep, 0),
            if (!is.na(m$filter_thresh)) sprintf(" (mean count below %s)", fmt(m$filter_thresh, 1)) else ""),
    sprintf("**Fold changes.** %s %s",
            switch(p$shrink,
              none = "Unshrunken (MLE) fold changes are noisy for low-count genes; prefer shrinkage for ranking and plots.",
              apeglm = "apeglm shrinks noisy low-count fold changes towards 0 but keeps large, well-supported changes. A good default.",
              ashr = "ashr shrinks with a flexible unimodal prior and works for any contrast. A good default.",
              normal = "The normal prior can over-shrink very large true changes; apeglm or ashr are usually preferred."),
            if (p$threshold_mode == "tested") sprintf("Genes are tested against |log2FC| > %s, the correct way to ask for a minimum effect size.", p$lfc)
            else sprintf("Genes must pass adjusted p < %s and then |log2FC| ≥ %s. Filtering after testing is common, but its error rate is not controlled for the fold-change cut-off; choose \"Test against threshold\" for that.", p$alpha, p$lfc)),
    sprintf("**Multiple testing.** %s at %s.", PADJ_TEXT[[p$padj_method]], p$alpha),
    sprintf("**p-value histogram.** About %s%% of genes appear unchanged (π₀ from p > 0.5), so roughly %s%% carry some signal. %s",
            fmt(100 * pi0, 0), fmt(100 * (1 - pi0), 0),
            if (hi_tail < 0.5 * mid) "The right tail is depleted (hump or slope): p-values may be conservative. Check for unmodelled batch or outlier samples."
            else "The shape (flat with a peak near 0) is what a well-fitted model looks like."),
    sprintf("**Result.** %d up and %d down at the chosen cut-offs.%s", m$n_up, m$n_down,
            if (m$n_up + m$n_down == 0) " Nothing passed: check power (replicates), effect sizes (MA plot) and whether the groups separate on the PCA." else ""))
  paste(out[!vapply(out, is.null, logical(1))], collapse = "\n\n")
}
