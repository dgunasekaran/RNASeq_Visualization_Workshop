# =============================================================================
# 02_demo_walkthrough.R — the 20-minute live demo
#
# From count matrix to publication-ready figures with dplyr + ggplot2 + DESeq2.
# Run section by section (Ctrl/Cmd + Enter). Timings are for the instructor.
#
#   [0-3 min]   1. Load + tidy (dplyr / tidyr)
#   [3-6 min]   2. Quick QC: library sizes, detection, composition
#   [6-9 min]   3. DESeq2 model: batch + genotype x time
#   [9-12 min]  4. Exploratory: VST, PCA, sample distances
#   [12-16 min] 5. Results: volcano, MA (shrinkage), p-value histogram, heatmap
#   [16-20 min] 6. Beyond the volcano: raincloud gene plots, UpSet overlaps,
#                  response-concordance quadrant, trajectory clusters,
#                  gene-module scores; multi-panel figure; launch the app
# =============================================================================

source("scripts/workshop_helpers.R")
suppressPackageStartupMessages({
  library(DESeq2)
  library(patchwork)
})
set.seed(2026)
fig_dir <- "figures/02_demo"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Load + tidy -----------------------------------------------------------------
# -----------------------------------------------------------------------------
# Same files (and format) as the app's demo data / upload
counts  <- read_counts("data/demo_counts.csv")          # genes x samples (integer matrix)
coldata <- read_metadata("data/demo_metadata.csv", colnames(counts)) |>
  as_workshop_coldata()                                 # one row per sample, factors set
counts  <- counts[, coldata$sample]                     # columns in design order
dim(counts)


# Long ("tidy") format: one row per gene x sample — what ggplot2 wants
counts_long <- as.data.frame(counts) |>
  rownames_to_column("gene_id") |>
  pivot_longer(-gene_id, names_to = "sample", values_to = "count") |>
  left_join(coldata, by = "sample")
counts_long

# dplyr verbs answer quick questions:
counts_long |>
  filter(gene_id %in% c("YFL014W", "YLL026W")) |>          # HSP12, HSP104
  group_by(gene_id, genotype, time) |>
  summarise(mean_count = mean(count), .groups = "drop")

# -----------------------------------------------------------------------------
# 2. Quick QC ----------------------------------------------------------------------
# -----------------------------------------------------------------------------
qc <- counts_long |>
  group_by(sample, genotype, time, replicate) |>
  summarise(lib_M = sum(count) / 1e6, detected = sum(count >= 10), .groups = "drop")

p_lib <- ggplot(qc, aes(time, lib_M, colour = genotype)) +
  geom_jitter(width = 0.15, height = 0, size = 1.8) +       # jitter: show every replicate
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.4, linewidth = 0.3,
               position = position_dodge(0.6)) +
  scale_colour_manual(values = pal_genotype) +
  labs(x = NULL, y = "Library size (M reads)", title = "Sequencing depth") +
  theme_pub()
p_lib
# -> mip6_39_20_4 has ~11.5 M reads, about half of the rest. Remember it.

# -----------------------------------------------------------------------------
# 3. DESeq2 model --------------------------------------------------------------------
# -----------------------------------------------------------------------------
pc <- classify_features(rownames(counts)) == "protein_coding"
cnt <- counts[pc, ]
cnt <- cnt[edgeR::filterByExpr(cnt, group = coldata$condition), ]
nrow(cnt)

cd <- data.frame(coldata, row.names = coldata$sample)

# (a) one factor per condition: easiest for pairwise contrasts
dds <- DESeqDataSetFromMatrix(cnt, cd, design = ~ batch + condition)
dds <- DESeq(dds)

# (b) factorial model: is the heat response DIFFERENT in mip6Δ? (interaction)
dds_int <- DESeqDataSetFromMatrix(cnt, cd, design = ~ batch + genotype * time)
dds_int <- DESeq(dds_int)
resultsNames(dds_int)

annot <- tibble(gene_id = rownames(cnt), symbol = gene_symbols(rownames(cnt)))

get_res <- function(A, B, shrink = TRUE) {
  r0 <- results(dds, contrast = c("condition", A, B), alpha = 0.05)
  r <- r0
  if (shrink) {
    r <- lfcShrink(dds, contrast = c("condition", A, B), res = r0, type = "ashr", quiet = TRUE)
    r$stat <- r0$stat          # keep the Wald statistic for ranking / GSEA (script 04)
  }
  as.data.frame(r) |> rownames_to_column("gene_id") |> as_tibble() |>
    left_join(annot, by = "gene_id")
}
res_heat_wt   <- get_res("WT_20min", "WT_0min")
res_heat_mip6 <- get_res("mip6_20min", "mip6_0min")
res_geno_20   <- get_res("mip6_20min", "WT_20min")
res_int <- results(dds_int, name = "genotypemip6.time20min", alpha = 0.05) |>
  as.data.frame() |> rownames_to_column("gene_id") |> as_tibble() |> left_join(annot, by = "gene_id")
summary(results(dds_int, name = "genotypemip6.time20min", alpha = 0.05))

# -----------------------------------------------------------------------------
# 4. Exploratory: VST, PCA, sample distances ----------------------------------------
# -----------------------------------------------------------------------------
vsd <- vst(dds, blind = FALSE)
vmat <- assay(vsd)

pca_data <- plotPCA(vsd, intgroup = c("genotype", "time", "batch"), ntop = 500, returnData = TRUE)
pv <- round(100 * attr(pca_data, "percentVar"))
p_pca <- ggplot(pca_data, aes(PC1, PC2, colour = time, shape = genotype)) +
  stat_ellipse(aes(group = interaction(time, genotype)), linewidth = 0.25,
               linetype = 2, show.legend = FALSE) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = pal_time) +
  labs(x = paste0("PC1 (", pv[1], "%)"), y = paste0("PC2 (", pv[2], "%)"),
       title = "PCA (top 500 variable genes, VST)") +
  coord_fixed(ratio = sqrt(pv[2] / pv[1])) +     # axes scaled by variance explained
  theme_pub()
p_pca

# Same PCA after removing the batch effect (for VISUALIZATION only; the DE
# model keeps batch as a covariate instead of altering the data)
vmat_nb <- limma::removeBatchEffect(vmat, batch = vsd$batch,
                                    design = model.matrix(~ condition, data = cd))

dists <- as.matrix(dist(t(vmat)))
ann <- data.frame(genotype = cd$genotype, time = cd$time, batch = cd$batch, row.names = cd$sample)
pheatmap::pheatmap(dists, annotation_col = ann, annotation_row = ann,
  annotation_colors = list(genotype = pal_genotype, time = pal_time,
                           batch = c(B1 = "grey30", B2 = "grey80")),
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(100),
  fontsize = 6, border_color = NA, show_colnames = FALSE,
  filename = file.path(fig_dir, "sample_distances.pdf"), width = 6.5, height = 5.5)

# -----------------------------------------------------------------------------
# 5. Results: volcano, MA, p-value histogram, heatmap --------------------------------
# -----------------------------------------------------------------------------
# For volcano we use the UNshrunken p-values with shrunken LFC (ashr keeps padj)
p_volc <- plot_volcano(res_heat_wt, title = "WT: 20 min 39 °C vs 30 °C")
p_volc$layers[[1]]$aes_params$size <- 1.5   # regular points (default 0.6)
p_volc$layers[[2]]$aes_params$size <- 2.5   # triangles for padj = 0 (default 1.4)
p_volc

# MA plots: why shrink fold changes? (Zhu, Ibrahim & Love 2019)
res_raw <- results(dds, contrast = c("condition", "mip6_20min", "WT_20min")) |>
  as.data.frame() |> mutate(type = "MLE (no shrinkage)")
res_shr <- as.data.frame(res_geno_20) |> mutate(type = "ashr shrinkage")
p_ma <- bind_rows(res_raw, res_shr) |>
  filter(!is.na(padj)) |>
  ggplot(aes(baseMean, log2FoldChange, colour = padj < 0.05)) +
  geom_point(size = 0.4, alpha = 0.6, stroke = 0) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_x_log10() + facet_wrap(~ type) +
  scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "#B2182B"), guide = "none") +
  labs(x = "Mean normalized count", y = expression(log[2]~FC),
       title = "mip6Δ vs WT at 20 min: noisy low-count LFCs shrink to 0") +
  theme_pub()
p_ma

# p-value histogram: a diagnostic every analysis should look at
p_phist <- ggplot(filter(res_raw, baseMean > 1), aes(pvalue)) +
  geom_histogram(breaks = seq(0, 1, 0.025), fill = "grey40", colour = "white", linewidth = 0.1) +
  labs(x = "p-value", y = "Genes",
       title = "Uniform + spike at 0 = healthy; hump or U-shape = check model") +
  theme_pub()

# Heatmap: top 40 genes by interaction significance, z-scored VST, k-means rows
top_int <- res_int |> filter(!is.na(padj)) |> arrange(padj) |> head(40)
hm <- vmat_nb[top_int$gene_id, ]
hm <- t(scale(t(hm)))
rownames(hm) <- make.unique(top_int$symbol)
ord <- order(cd$genotype, cd$time)
pheatmap::pheatmap(hm[, ord], cluster_cols = FALSE, annotation_col = ann[ord, c("genotype", "time")],
  annotation_colors = list(genotype = pal_genotype, time = pal_time),
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(101),
  breaks = seq(-2.5, 2.5, length.out = 102), kmeans_k = NA, cutree_rows = 3,
  gaps_col = cumsum(table(interaction(cd$time, cd$genotype)[ord]))[-6],
  fontsize = 6, border_color = NA, show_colnames = FALSE,
  main = "Genes whose 20-min heat response depends on MIP6 (z-score, batch-corrected VST)",
  filename = file.path(fig_dir, "heatmap_interaction.pdf"), width = 6, height = 6.5)

# -----------------------------------------------------------------------------
# 6. Beyond the volcano ------------------------------------------------------------
# -----------------------------------------------------------------------------
norm <- counts(dds, normalized = TRUE)
gene_long <- function(ids) {
  as.data.frame(norm[ids, , drop = FALSE]) |>
    rownames_to_column("gene_id") |>
    pivot_longer(-gene_id, names_to = "sample", values_to = "norm_count") |>
    left_join(coldata, by = "sample") |>
    left_join(annot, by = "gene_id")
}

# 6a. Raincloud plots (Allen et al. 2021): distribution + box + every replicate
show_genes <- c("YFL014W", "YBR072W", "YLL026W", "YER103W")   # HSP12, HSP26, HSP104, SSA4
p_rain <- gene_long(show_genes) |>
  ggplot(aes(time, norm_count + 1, fill = genotype, colour = genotype)) +
  ggdist::stat_halfeye(adjust = 1, width = 0.5, .width = 0, justification = -0.35,
                       alpha = 0.5, point_colour = NA, position = position_dodge(0.8)) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.3, linewidth = 0.25,
               position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.8),
             size = 0.9) +
  scale_y_log10(labels = scales::label_comma()) +
  facet_wrap(~ symbol, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = pal_genotype) + scale_colour_manual(values = pal_genotype) +
  labs(x = NULL, y = "Normalized count (log)", title = "Heat-shock genes, every replicate shown") +
  theme_pub()
p_rain

# 6b. UpSet plot (Lex et al. 2014; Conway et al. 2017) instead of a Venn diagram
sig_set <- function(r) r$gene_id[!is.na(r$padj) & r$padj < 0.05 & abs(r$log2FoldChange) > 1]
sets <- list(`WT heat 20'`    = sig_set(res_heat_wt),
             `mip6 heat 20'`  = sig_set(res_heat_mip6),
             `WT heat 120'`   = sig_set(get_res("WT_120min", "WT_0min")),
             `mip6 vs WT 20'` = sig_set(res_geno_20))
pdf(file.path(fig_dir, "upset_de_sets.pdf"), width = 6, height = 3.5, onefile = FALSE)
print(UpSetR::upset(UpSetR::fromList(sets), order.by = "freq", nsets = 4,
                    main.bar.color = "grey25", sets.bar.color = "#0072B2", text.scale = 1))
dev.off()

# 6c. Response-concordance (quadrant) plot: heat response in WT vs in mip6Δ.
#     Points off the diagonal = genotype-dependent heat response.
quad <- inner_join(select(res_heat_wt, gene_id, symbol, lfc_wt = log2FoldChange),
                   select(res_heat_mip6, gene_id, lfc_mip6 = log2FoldChange), by = "gene_id") |>
  left_join(select(res_int, gene_id, padj_int = padj), by = "gene_id") |>
  mutate(interaction_sig = !is.na(padj_int) & padj_int < 0.05)
p_quad <- ggplot(quad, aes(lfc_wt, lfc_mip6)) +
  geom_hline(yintercept = 0, linewidth = 0.2) + geom_vline(xintercept = 0, linewidth = 0.2) +
  geom_abline(linetype = 2, linewidth = 0.3, colour = "grey40") +
  geom_point(data = filter(quad, !interaction_sig), size = 0.4, colour = "grey75", stroke = 0) +
  geom_point(data = filter(quad, interaction_sig), size = 0.9, colour = "#D55E00") +
  ggrepel::geom_text_repel(data = quad |> filter(interaction_sig) |> arrange(padj_int) |> head(12),
                           aes(label = symbol), size = 2.2, min.segment.length = 0) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.6, size = 2.5,
           label = paste0("r = ", round(cor(quad$lfc_wt, quad$lfc_mip6, use = "complete.obs"), 2))) +
  coord_equal() +
  labs(x = expression(log[2]~FC~heat~(WT)), y = expression(log[2]~FC~heat~(mip6*Delta)),
       title = "Heat response: WT vs mip6Δ",
       subtitle = "orange = significant genotype x time interaction") +
  theme_pub()
p_quad

# -----------------------------------------------------------------------------
# Assemble a publication figure (183 mm = double column) ------------------------
# -----------------------------------------------------------------------------
fig1 <- (p_pca / p_quad / p_volc / p_rain) +
  plot_annotation(tag_levels = "a") +
  plot_layout(heights = c(0.8, 1.1, 1, 1.1))
save_pub(fig1, file.path(fig_dir, "Figure1_heat_shock.pdf"), width_mm = 230, height_mm = 230)
save_pub(fig1, file.path(fig_dir, "Figure1_heat_shock.png"), width_mm = 230, height_mm = 230, dpi = 300)
save_pub(p_ma / p_phist, file.path(fig_dir, "FigureS_shrinkage_pvalues.pdf"), width_mm = 183, height_mm = 120)

# Save results tables for the hands-on and for the app
dir.create("results", showWarnings = FALSE)
write_csv(res_heat_wt, "results/DE_WT_20min_vs_0min.csv")
write_csv(res_geno_20, "results/DE_mip6_vs_WT_20min.csv")
write_csv(res_int,     "results/DE_interaction_genotype_x_20min.csv")

# -----------------------------------------------------------------------------
# Launch the interactive dashboard -----------------------------------------------
# -----------------------------------------------------------------------------
# shiny::runApp()      # from the repo folder: Upload tab -> "Use demo data"
