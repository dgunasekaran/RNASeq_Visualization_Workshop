# =============================================================================
# 01_qc_normalization_subsampling.R
#
# How do sequencing depth, number of replicates, normalization choice and ONE
# irregular replicate change what we conclude?  Everything starts from the
# COUNT MATRIX (no FASTQ/BAM needed): failures are simulated from real counts.
#
#   Part A  Know your matrix: feature composition, genotype marker, hidden batch
#   Part B  Synthetic QC panel: low depth, low complexity, gDNA, ncRNA,
#           degradation, noisy replicate, contamination, sample swap
#   Part C  Depth subsampling (binomial thinning): detection, concordance, DE
#   Part D  Replicate subsampling (n = 2, 3, 4): DE recall and precision
#   Part E  Normalization strategies: TC, UQ, TMM, median-of-ratios, quantile, VST
#   Part F  One irregular replicate: effect on DE, and three remedies
#
# Run from the repo folder (the one with app.R). FAST = TRUE takes ~5-8 min on a laptop;
# FAST = FALSE adds more random draws (~20-30 min).
# Figures -> figures/01_qc/   Tables -> results/
# =============================================================================

source("scripts/workshop_helpers.R")
suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(limma)
  library(patchwork)
})

FAST <- TRUE
set.seed(2026)
fig_dir <- "figures/01_qc"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# Same files (and format) as the app's demo data / upload
counts  <- read_counts("data/demo_counts.csv")          # genes x samples (integer matrix)
coldata <- read_metadata("data/demo_metadata.csv", colnames(counts)) |>
  as_workshop_coldata()                                 # one row per sample, factors set
counts  <- counts[, coldata$sample]                     # columns in design order
fclass  <- classify_features(rownames(counts))
pc_ids  <- rownames(counts)[fclass == "protein_coding"]

# Fit DESeq2 with ~ batch + condition, dropping batch if it is not estimable
# (can happen after subsampling replicates).
fit_dds <- function(cnt, cd) {
  cd <- as.data.frame(cd)
  rownames(cd) <- cd$sample
  cd$condition <- droplevels(cd$condition)
  cd$batch <- droplevels(cd$batch)
  design <- ~ batch + condition
  mm <- model.matrix(~ batch + condition, cd)
  if (nlevels(cd$batch) < 2 || qr(mm)$rank < ncol(mm)) design <- ~ condition
  dds <- DESeqDataSetFromMatrix(cnt, cd, design)
  DESeq(dds, quiet = TRUE)
}
de_calls <- function(dds, A, B, alpha = 0.05) {
  r <- results(dds, contrast = c("condition", A, B), alpha = alpha)
  tibble(gene_id = rownames(r), log2FoldChange = r$log2FoldChange,
         pvalue = r$pvalue, padj = r$padj,
         sig = !is.na(r$padj) & r$padj < alpha)
}
compare_to_truth <- function(test, truth) {
  j <- inner_join(test, truth, by = "gene_id", suffix = c("", "_truth"))
  tp <- sum(j$sig & j$sig_truth)
  tibble(n_de      = sum(test$sig),
         recall    = tp / max(1, sum(truth$sig)),
         precision = if (sum(test$sig) > 0) tp / sum(test$sig) else NA_real_,
         lfc_cor   = cor(j$log2FoldChange, j$log2FoldChange_truth, use = "complete.obs"))
}

# =============================================================================
# Part A — Know your matrix
# =============================================================================
message("Part A: feature composition, genotype marker, batch")

comp <- as.data.frame(counts) |>
  rownames_to_column("gene_id") |>
  mutate(class = fclass) |>
  pivot_longer(-c(gene_id, class), names_to = "sample", values_to = "count") |>
  group_by(sample, class) |>
  summarise(reads = sum(count), .groups = "drop_last") |>
  mutate(pct = 100 * reads / sum(reads)) |>
  ungroup() |>
  mutate(sample = factor(sample, levels = coldata$sample))

p_comp <- ggplot(comp, aes(sample, pct, fill = class)) +
  geom_col(width = 0.85) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = NULL, y = "% of counted reads", fill = "Feature class",
       title = "Where do the reads go?") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))

top_features <- tibble(gene_id = rownames(counts), class = fclass,
                       share = 100 * rowSums(counts) / sum(counts)) |>
  arrange(desc(share)) |> head(10)
print(top_features)   # SCR1 (signal recognition particle RNA) tops the list

# Genotype check: the deleted gene must have ~0 reads in mip6Δ
marker <- tibble(sample = coldata$sample, genotype = coldata$genotype,
                 cpm = 1e6 * counts[MIP6_ORF, ] / colSums(counts))
p_marker <- ggplot(marker, aes(factor(sample, levels = sample), cpm, fill = genotype)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = pal_genotype) +
  labs(x = NULL, y = "MIP6 (YHR015W) CPM", title = "Genotype check: deleted gene") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))

# Hidden batch: remove each condition's mean, then look for structure that is left.
lc_pc <- log_cpm(counts[pc_ids, ])
expressed <- rowSums(counts[pc_ids, ] >= 10) == ncol(counts)
lc_pc <- lc_pc[expressed, ]
grp_means <- vapply(levels(coldata$condition),
                    function(g) rowMeans(lc_pc[, coldata$condition == g, drop = FALSE]),
                    numeric(nrow(lc_pc)))
resid <- lc_pc - grp_means[, as.character(coldata$condition)]
pca_res <- prcomp(t(resid))
var_res <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
p_batch <- as_tibble(pca_res$x[, 1:2]) |>
  mutate(sample = coldata$sample, replicate = factor(coldata$replicate)) |>
  ggplot(aes(PC1, PC2, colour = replicate)) +
  geom_point(size = 2) +
  scale_colour_manual(values = c("1" = "#0072B2", "3" = "#56B4E9", "2" = "#D55E00", "4" = "#E69F00")) +
  labs(x = paste0("Residual PC1 (", var_res[1], "%)"), y = paste0("Residual PC2 (", var_res[2], "%)"),
       title = "Within-condition residuals: replicates 1+3 vs 2+4",
       colour = "Replicate") +
  theme_pub()

save_pub((p_comp / p_marker) | p_batch, file.path(fig_dir, "A_matrix_overview.pdf"),
         width_mm = 183, height_mm = 120)
# Teaching point: the paper reports two culture days. The data say replicates
# {1,3} and {2,4} pair up -> we model batch = B1 (1,3) / B2 (2,4).

# =============================================================================
# Part B — Synthetic QC panel (clearly named SYN## samples)
# =============================================================================
message("Part B: synthetic QC samples")

lens  <- get_gene_lengths(rownames(counts))       # NULL if TxDb not installed
panel <- make_synthetic_panel(counts, coldata, lengths = lens)
qc    <- qc_metrics(panel$counts, panel$coldata)

qc |>
  select(sample, qc_issue, lib_size_M, genes_detected, pct_ncRNA, ghost_genes,
         cor_to_group, noise_mad, saturation_gain, marker_cpm, flags) |>
  write_csv("results/01_qc_metrics.csv")
print(qc |> select(sample, flags), n = Inf)
# Note: one REAL sample, mip6_39_20_4 (lowest depth, 11.5 M), is flagged for
# excess noise -> a naturally irregular replicate. Keep it in mind for Part F.

qc_long <- qc |>
  select(sample, synthetic, `Library size (M)` = lib_size_M,
         `Genes detected` = genes_detected, `% ncRNA/rRNA` = pct_ncRNA,
         `Reads on silent genes` = ghost_genes, `Median r to replicates` = cor_to_group,
         `Noise (MAD log2)` = noise_mad, `Saturation gain` = saturation_gain,
         `MIP6 CPM` = marker_cpm) |>
  pivot_longer(-c(sample, synthetic), names_to = "metric") |>
  mutate(sample = factor(sample, levels = rev(panel$coldata$sample)),
         status = ifelse(synthetic, "synthetic", "real"))

p_qc <- ggplot(qc_long, aes(value, sample, fill = status)) +
  geom_col(width = 0.75) +
  facet_wrap(~ metric, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = pal_status) +
  labs(x = NULL, y = NULL, fill = NULL, title = "Count-level QC dashboard") +
  theme_pub(base_size = 7) +
  theme(legend.position = "bottom")
save_pub(p_qc, file.path(fig_dir, "B_qc_dashboard.pdf"), width_mm = 260, height_mm = 140)

# Saturation (rarefaction) curves: does detection keep rising with depth?
sat_samples <- c("wt_30_0_1", "wt_30_0_2", "SYN01_lowDepth", "SYN02_lowComplexity", "SYN03_gDNA")
sat_samples <- intersect(sat_samples, colnames(panel$counts))
fracs <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1)
sat <- expand_grid(sample = sat_samples, frac = fracs) |>
  mutate(
    reads_M  = map2_dbl(sample, frac, \(s, f) f * sum(panel$counts[pc_ids, s]) / 1e6),
    detected = map2_dbl(sample, frac, \(s, f) sum(thin_counts(panel$counts[pc_ids, s], f) >= 1))
  )
p_sat <- ggplot(sat, aes(reads_M, detected, colour = sample)) +
  geom_line(linewidth = 0.5) + geom_point(size = 1) +
  scale_x_log10() +
  scale_colour_manual(values = c("#0072B2", "#56B4E9", "#009E73", "#D55E00", "#CC79A7")) +
  labs(x = "Reads (millions, log scale)", y = "Protein-coding genes detected (≥1 read)",
       title = "Saturation curves from counts", colour = NULL) +
  theme_pub()

# PCA: fit on REAL samples, then PROJECT the synthetic ones onto the same axes
lc_all <- log_cpm(panel$counts[pc_ids, ])[expressed, ]
real <- !panel$coldata$synthetic
top <- order(matrixStats::rowVars(lc_all[, real]), decreasing = TRUE)[1:500]
pca <- prcomp(t(lc_all[top, real]), center = TRUE)
proj <- predict(pca, t(lc_all[top, ]))
ve <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pca_df <- as_tibble(proj[, 1:2]) |> bind_cols(panel$coldata)
p_pca_syn <- ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(data = filter(pca_df, !synthetic), aes(colour = condition), size = 1.8) +
  geom_point(data = filter(pca_df, synthetic), shape = 4, size = 2.5, stroke = 0.8) +
  ggrepel::geom_text_repel(data = filter(pca_df, synthetic), aes(label = sample),
                           size = 2.2, min.segment.length = 0) +
  scale_colour_manual(values = pal_condition) +
  labs(x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)"),
       title = "Synthetic samples projected onto real-sample PCA") +
  theme_pub()
save_pub(p_sat | p_pca_syn, file.path(fig_dir, "B_saturation_pca.pdf"),
         width_mm = 183, height_mm = 80)

# Sample-to-sample correlation heatmap (real + synthetic)
ann <- as.data.frame(panel$coldata[, c("condition", "synthetic")])
rownames(ann) <- panel$coldata$sample
ann$synthetic <- factor(ifelse(ann$synthetic, "yes", "no"))
pheatmap::pheatmap(cor(lc_all), annotation_col = ann,
                   annotation_colors = list(condition = pal_condition,
                                            synthetic = c(no = "grey85", yes = "#D55E00")),
                   color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdYlBu")))(100),
                   fontsize = 6, border_color = NA,
                   filename = file.path(fig_dir, "B_correlation_heatmap.pdf"),
                   width = 8, height = 7)

# =============================================================================
# Part C — Depth subsampling
# =============================================================================
message("Part C: depth subsampling (DESeq2 fits — takes a few minutes)")

sub_conditions <- c("WT_0min", "WT_20min", "mip6_20min")
cd_sub  <- coldata |> filter(condition %in% sub_conditions) |> mutate(condition = droplevels(condition))
cnt_sub <- counts[pc_ids, cd_sub$sample]
keep    <- filterByExpr(cnt_sub, group = cd_sub$condition)
cnt_sub <- cnt_sub[keep, ]
contrasts <- list(heat = c("WT_20min", "WT_0min"), genotype = c("mip6_20min", "WT_20min"))

dds_full <- fit_dds(cnt_sub, cd_sub)
truth <- map(contrasts, \(ct) de_calls(dds_full, ct[1], ct[2]))
map_int(truth, \(t) sum(t$sig))   # heat: thousands; genotype: far fewer

replicate_cor <- function(cnt, cd) {
  lc <- log_cpm(cnt + 0L)
  map_dbl(levels(cd$condition), \(g) {
    s <- cd$sample[cd$condition == g]
    cm <- cor(lc[, s]); median(cm[upper.tri(cm)])
  }) |> median()
}

depth_fracs <- c(0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1)
n_seeds <- if (FAST) 2 else 5
depth_res <- expand_grid(frac = depth_fracs, seed = seq_len(n_seeds)) |>
  filter(!(frac == 1 & seed > 1)) |>
  mutate(out = map2(frac, seed, \(f, s) {
    set.seed(s)
    cnt <- if (f < 1) thin_matrix(cnt_sub, f) else cnt_sub
    dds <- fit_dds(cnt, cd_sub)
    de <- imap_dfr(contrasts, \(ct, nm) compare_to_truth(de_calls(dds, ct[1], ct[2]), truth[[nm]]) |>
                     mutate(contrast = nm))
    tibble(median_reads_M = median(colSums(cnt)) / 1e6,
           detected_ge1   = median(colSums(cnt >= 1)),
           detected_ge10  = median(colSums(cnt >= 10)),
           replicate_r    = replicate_cor(cnt, cd_sub),
           de = list(de))
  })) |>
  unnest(out)

depth_de <- depth_res |> select(frac, seed, median_reads_M, de) |> unnest(de)
write_csv(select(depth_res, -de), "results/01_depth_detection.csv")
write_csv(depth_de, "results/01_depth_de.csv")

p_d1 <- depth_res |>
  pivot_longer(c(detected_ge1, detected_ge10), names_to = "threshold", values_to = "genes") |>
  mutate(threshold = recode(threshold, detected_ge1 = "≥1 read", detected_ge10 = "≥10 reads")) |>
  ggplot(aes(median_reads_M, genes, colour = threshold)) +
  geom_point(size = 1) + stat_summary(fun = mean, geom = "line") +
  scale_x_log10() + scale_colour_manual(values = c("#0072B2", "#D55E00")) +
  labs(x = "Reads per library (M)", y = "Genes detected", colour = NULL,
       title = "Detection vs depth") + theme_pub()
p_d2 <- ggplot(depth_res, aes(median_reads_M, replicate_r)) +
  geom_point(size = 1) + stat_summary(fun = mean, geom = "line") +
  scale_x_log10() +
  labs(x = "Reads per library (M)", y = "Median Pearson r (replicates)",
       title = "Replicate concordance") + theme_pub()
p_d3 <- ggplot(depth_de, aes(median_reads_M, recall, colour = contrast)) +
  geom_point(size = 1) + stat_summary(fun = mean, geom = "line") +
  scale_x_log10() + scale_y_continuous(limits = c(0, 1)) +
  scale_colour_manual(values = c(heat = "#E69F00", genotype = "#009E73")) +
  labs(x = "Reads per library (M)", y = "Recall of full-depth DE genes",
       title = "DE recovery") + theme_pub()
p_d4 <- ggplot(depth_de, aes(median_reads_M, n_de, colour = contrast)) +
  geom_point(size = 1) + stat_summary(fun = mean, geom = "line") +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c(heat = "#E69F00", genotype = "#009E73")) +
  labs(x = "Reads per library (M)", y = "DE genes (padj < 0.05)", title = "Number of DE calls") +
  theme_pub()
save_pub((p_d1 | p_d2) / (p_d3 | p_d4) + plot_annotation(tag_levels = "a"),
         file.path(fig_dir, "C_depth_subsampling.pdf"), width_mm = 183, height_mm = 130)
# Teaching point: large effects (heat) are found at 1-2 M reads; subtle effects
# (genotype) need depth AND replicates.

# =============================================================================
# Part D — Replicate subsampling
# =============================================================================
message("Part D: replicate subsampling")

n_draws <- if (FAST) 5 else 15
draw_reps <- function(n) {
  cd_sub |> group_by(condition) |> slice_sample(n = n) |> ungroup()
}
rep_res <- expand_grid(n_rep = 2:3, draw = seq_len(n_draws)) |>
  bind_rows(tibble(n_rep = 4L, draw = 1L)) |>
  mutate(de = map2(n_rep, draw, \(n, d) {
    set.seed(100 * n + d)
    cd <- if (n == 4) cd_sub else draw_reps(n)
    dds <- fit_dds(cnt_sub[, cd$sample], cd)
    imap_dfr(contrasts, \(ct, nm) compare_to_truth(de_calls(dds, ct[1], ct[2]), truth[[nm]]) |>
               mutate(contrast = nm))
  })) |>
  unnest(de)
write_csv(rep_res, "results/01_replicate_subsampling.csv")

p_r <- rep_res |>
  pivot_longer(c(n_de, recall, precision), names_to = "metric") |>
  ggplot(aes(factor(n_rep), value, colour = contrast)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(0.75)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
             size = 0.8, alpha = 0.7) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_colour_manual(values = c(heat = "#E69F00", genotype = "#009E73")) +
  labs(x = "Biological replicates per condition", y = NULL, colour = NULL,
       title = "Replicates matter more than depth for subtle effects") +
  theme_pub()
save_pub(p_r, file.path(fig_dir, "D_replicate_subsampling.pdf"), width_mm = 183, height_mm = 70)
# Caveat: "truth" = the 4-replicate result, itself uncertain. Schurch et al.
# (2016, RNA) used 48 replicates per condition to get a better gold standard.

# =============================================================================
# Part E — Normalization strategies
# =============================================================================
message("Part E: normalization strategies")

cnt_n <- counts[pc_ids, ]
cnt_n <- cnt_n[filterByExpr(cnt_n, group = coldata$condition), ]
methods <- c("none", "TC", "UQ", "TMM", "MoR", "quantile", "VST")
norm_list <- set_names(map(methods, \(m) normalize_counts(cnt_n, m)), methods)

# E1. RLE plots (relative log expression)
rle_df <- imap_dfr(norm_list, \(m, nm) rle_long(m) |> mutate(method = nm)) |>
  left_join(coldata, by = "sample") |>
  mutate(method = factor(method, levels = methods),
         sample = factor(sample, levels = coldata$sample))
p_rle <- ggplot(rle_df, aes(sample, rle, fill = time)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.2) +
  geom_hline(yintercept = 0, colour = "red", linewidth = 0.3) +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  facet_wrap(~ method, ncol = 2) +
  scale_fill_manual(values = pal_time) +
  labs(x = NULL, y = "Relative log expression", title = "RLE by normalization") +
  theme_pub(base_size = 7) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# E2. Metrics: global shift between conditions + replicate tightness
norm_metrics <- imap_dfr(norm_list, \(m, nm) {
  g <- coldata$condition
  mean_by <- \(lv) rowMeans(m[, g == lv, drop = FALSE])
  within_sd <- map_dbl(levels(g), \(lv) median(matrixStats::rowSds(m[, g == lv, drop = FALSE])))
  ref <- intersect(REF_GENES, rownames(m))
  tibble(method = nm,
         median_lfc_20v0  = median(mean_by("WT_20min")  - mean_by("WT_0min")),
         median_lfc_120v0 = median(mean_by("WT_120min") - mean_by("WT_0min")),
         rle_iqr          = median(apply(m - matrixStats::rowMedians(m), 2, IQR)),
         within_group_sd  = median(within_sd),
         ref_gene_sd      = mean(matrixStats::rowSds(m[ref, , drop = FALSE])))
})
print(norm_metrics)
write_csv(norm_metrics, "results/01_normalization_metrics.csv")
# Under "most genes do not change", median log2FC should be ~0. Total-count
# scaling gives about +0.3 (0.28-0.29 depending on the filter) at 20 min: a few HSP mRNAs soak up reads (composition
# bias), so every other gene looks UP. TMM / median-of-ratios fix this.
# (If heat shock changed TOTAL mRNA per cell, only spike-ins could tell:
#  Lovén et al. 2012, Cell.)

# E3. Size factors: total count vs median-of-ratios, colored by time
sf_df <- tibble(sample = colnames(cnt_n),
                TC = size_factors(cnt_n, "TC"), TMM = size_factors(cnt_n, "TMM"),
                MoR = size_factors(cnt_n, "MoR")) |>
  left_join(coldata, by = "sample")
p_sf <- ggplot(sf_df, aes(TC, MoR, colour = time, shape = genotype)) +
  geom_abline(linetype = 2, linewidth = 0.3) +
  geom_point(size = 2) +
  scale_colour_manual(values = pal_time) +
  labs(x = "Size factor: total count", y = "Size factor: median-of-ratios",
       title = "Composition bias moves size factors") +
  theme_pub()

# E4. Consequence for DE: same DESeq2 model, different size factors
dds_n <- DESeqDataSetFromMatrix(cnt_n, data.frame(coldata, row.names = coldata$sample),
                                ~ batch + condition)
de_by_norm <- map_dfr(c("TC", "UQ", "TMM", "MoR"), \(m) {
  sizeFactors(dds_n) <- size_factors(cnt_n, m)
  d <- DESeq(dds_n, quiet = TRUE)
  r <- results(d, contrast = c("condition", "WT_20min", "WT_0min"), alpha = 0.05)
  tibble(method = m,
         up   = sum(r$padj < 0.05 & r$log2FoldChange > 0, na.rm = TRUE),
         down = sum(r$padj < 0.05 & r$log2FoldChange < 0, na.rm = TRUE))
})
print(de_by_norm)
p_de_norm <- de_by_norm |>
  pivot_longer(c(up, down), names_to = "direction", values_to = "n") |>
  mutate(n = ifelse(direction == "down", -n, n),
         method = factor(method, levels = c("TC", "UQ", "TMM", "MoR"))) |>
  ggplot(aes(method, n, fill = direction)) +
  geom_col(width = 0.7) + geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c(up = "#B2182B", down = "#2166AC")) +
  scale_y_continuous(labels = abs) +
  labs(x = NULL, y = "DE genes (WT 20 vs 0 min)", fill = NULL,
       title = "Normalization shifts the up/down balance") +
  theme_pub()

save_pub(p_rle | (p_sf / p_de_norm), file.path(fig_dir, "E_normalization.pdf"),
         width_mm = 183, height_mm = 150)

# =============================================================================
# Part F — One irregular replicate
# =============================================================================
message("Part F: irregular replicate scenarios")

target  <- "mip6_39_20_1"     # replicate we corrupt (condition mip6_20min)
cnt_f   <- counts[pc_ids, ]
cnt_f   <- cnt_f[filterByExpr(cnt_f, group = coldata$condition), ]
contr_f <- list(genotype_20 = c("mip6_20min", "WT_20min"), heat_mip6 = c("mip6_20min", "mip6_0min"))

dds_clean <- fit_dds(cnt_f, coldata)
truth_f   <- map(contr_f, \(ct) de_calls(dds_clean, ct[1], ct[2]))

make_irregular <- function(type) {
  x <- counts[rownames(cnt_f), target]
  switch(type,
    noisy          = sim_noise(x, sdlog = 0.6),
    low_complexity = sim_low_complexity(x),
    contamination  = sim_mixture(x, counts[rownames(cnt_f), "wt_39_20_1"], frac_x = 0.6),
    swap           = thin_counts(counts[rownames(cnt_f), "mip6_30_0_1"], 0.95)
  )
}

limma_calls <- function(cnt, cd, quality_weights = FALSE) {
  cd <- as.data.frame(cd)
  cd$condition <- droplevels(cd$condition)
  design <- model.matrix(~ 0 + condition + batch, cd)
  colnames(design) <- make.names(colnames(design))
  y <- calcNormFactors(DGEList(cnt))
  v <- if (quality_weights) voomWithQualityWeights(y, design) else voom(y, design)
  fit <- lmFit(v, design)
  sw <- if (quality_weights) {
    if (!is.null(v$targets$sample.weights)) v$targets$sample.weights else v$sample.weights
  } else NULL
  out <- imap(contr_f, \(ct, nm) {
    cm <- makeContrasts(contrasts = paste0("condition", ct[1], " - condition", ct[2]), levels = design)
    tt <- topTable(eBayes(contrasts.fit(fit, cm)), number = Inf, sort.by = "none")
    tibble(gene_id = rownames(tt), log2FoldChange = tt$logFC, pvalue = tt$P.Value,
           padj = tt$adj.P.Val, sig = tt$adj.P.Val < 0.05)
  })
  list(calls = out, weights = sw)
}

types <- if (FAST) c("noisy", "contamination", "swap") else
                   c("noisy", "low_complexity", "contamination", "swap")
irr <- map_dfr(types, \(type) {
  set.seed(7)
  cnt_bad <- cnt_f
  cnt_bad[, target] <- make_irregular(type)

  # strategy 1: keep everything (naive)
  dds_bad <- fit_dds(cnt_bad, coldata)
  naive <- imap_dfr(contr_f, \(ct, nm) {
    r <- results(dds_bad, contrast = c("condition", ct[1], ct[2]), alpha = 0.05)
    compare_to_truth(de_calls(dds_bad, ct[1], ct[2]), truth_f[[nm]]) |>
      mutate(contrast = nm, cooks_flagged = sum(is.na(r$pvalue) & r$baseMean > 0))
  }) |> mutate(strategy = "1. keep irregular (DESeq2)")

  # strategy 2: QC-driven removal of the flagged sample
  keep_s <- coldata$sample != target
  dds_rm <- fit_dds(cnt_bad[, keep_s], coldata[keep_s, ])
  removed <- imap_dfr(contr_f, \(ct, nm) compare_to_truth(de_calls(dds_rm, ct[1], ct[2]), truth_f[[nm]]) |>
                        mutate(contrast = nm)) |>
    mutate(strategy = "2. remove sample (DESeq2)")

  # strategy 3: down-weight instead of remove (limma-voom sample quality weights)
  lq <- limma_calls(cnt_bad, coldata, quality_weights = TRUE)
  weighted <- imap_dfr(contr_f, \(ct, nm) compare_to_truth(lq$calls[[nm]], truth_f[[nm]]) |>
                         mutate(contrast = nm)) |>
    mutate(strategy = "3. voom quality weights",
           target_weight = if (is.null(lq$weights)) NA_real_ else
                             lq$weights[which(coldata$sample == target)])

  bind_rows(naive, removed, weighted) |> mutate(type = type)
})

# reference: removing a GOOD replicate (cost of n = 3 alone)
keep_s <- coldata$sample != target
dds_n3 <- fit_dds(cnt_f[, keep_s], coldata[keep_s, ])
n3 <- imap_dfr(contr_f, \(ct, nm) compare_to_truth(de_calls(dds_n3, ct[1], ct[2]), truth_f[[nm]]) |>
                 mutate(contrast = nm, strategy = "0. clean, n = 3 (reference)", type = "none"))
irr <- bind_rows(irr, n3)
write_csv(irr, "results/01_irregular_replicate.csv")
print(irr |> select(type, strategy, contrast, n_de, recall, precision, target_weight, cooks_flagged))

p_irr <- irr |>
  filter(type != "none") |>
  pivot_longer(c(recall, precision), names_to = "metric") |>
  ggplot(aes(type, value, fill = strategy)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_hline(data = n3 |> pivot_longer(c(recall, precision), names_to = "metric"),
             aes(yintercept = value), linetype = 2, linewidth = 0.3) +
  facet_grid(metric ~ contrast) +
  scale_fill_manual(values = c("#BBBBBB", "#0072B2", "#009E73")) +
  labs(x = "Type of irregular replicate", y = "vs clean 4-replicate result", fill = NULL,
       title = "One bad replicate: keep, remove, or down-weight?",
       caption = "Dashed line: dropping a GOOD replicate (n = 3)") +
  theme_pub() + theme(legend.position = "bottom")

# PCA showing where the irregular replicate falls (noisy example)
set.seed(7)
cnt_bad <- cnt_f
cnt_bad[, target] <- make_irregular("noisy")
vs <- assay(vst(DESeqDataSetFromMatrix(cnt_bad, data.frame(coldata, row.names = coldata$sample), ~ condition),
                blind = TRUE))
pc_b <- prcomp(t(vs[order(matrixStats::rowVars(vs), decreasing = TRUE)[1:500], ]))
p_irr_pca <- as_tibble(pc_b$x[, 1:2]) |> bind_cols(coldata) |>
  mutate(irregular = sample == target) |>
  ggplot(aes(PC1, PC2, colour = condition)) +
  geom_point(aes(shape = irregular), size = 2) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 8)) +
  scale_colour_manual(values = pal_condition) +
  labs(title = paste("Noisy", target, "on PCA")) + theme_pub()

save_pub(p_irr_pca | p_irr, file.path(fig_dir, "F_irregular_replicate.pdf"),
         width_mm = 240, height_mm = 110)

# The REAL borderline sample: what happens to the genotype contrast without it?
real_check <- {
  ks <- coldata$sample != "mip6_39_20_4"
  d  <- fit_dds(cnt_f[, ks], coldata[ks, ])
  compare_to_truth(de_calls(d, "mip6_20min", "WT_20min"), truth_f$genotype_20)
}
message("Dropping real sample mip6_39_20_4 -> genotype_20 DE genes: ", real_check$n_de,
        " (with it: ", sum(truth_f$genotype_20$sig), ")")

sessionInfo()
