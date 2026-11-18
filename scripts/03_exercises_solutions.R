# =============================================================================
# 03_exercises_solutions.R — worked answers for 03_exercises.R
# =============================================================================
source("scripts/workshop_helpers.R")
suppressPackageStartupMessages(library(DESeq2))
set.seed(42)
# Same files (and format) as the app's demo data / upload
counts  <- read_counts("data/demo_counts.csv")          # genes x samples (integer matrix)
coldata <- read_metadata("data/demo_metadata.csv", colnames(counts)) |>
  as_workshop_coldata()                                 # one row per sample, factors set
counts  <- counts[, coldata$sample]                     # columns in design order

counts_long <- as.data.frame(counts) |>
  rownames_to_column("gene_id") |>
  pivot_longer(-gene_id, names_to = "sample", values_to = "count") |>
  left_join(coldata, by = "sample")

# ---- Exercise 1: dplyr warm-up (5 min) --------------------------------------
# 1a. Library size per sample, smallest first
lib <- counts_long |>
  group_by(sample) |>
  summarise(lib_M = sum(count) / 1e6) |>
  arrange(lib_M)
lib                                    # mip6_39_20_4 ~ 11.5 M

# 1b. % of reads in the 50 most-abundant features of each sample
counts_long |>
  group_by(sample) |>
  slice_max(count, n = 50, with_ties = FALSE) |>
  summarise(top50 = sum(count)) |>
  left_join(lib, by = "sample") |>
  mutate(pct_top50 = 100 * top50 / (lib_M * 1e6)) |>
  arrange(desc(pct_top50))

# 1c. HSP12 (YFL014W) mean CPM per condition and fold change 20 vs 0 min in WT
hsp12 <- counts_long |>
  group_by(sample) |> mutate(cpm = 1e6 * count / sum(count)) |> ungroup() |>
  filter(gene_id == "YFL014W") |>
  group_by(genotype, time) |>
  summarise(mean_cpm = mean(cpm), .groups = "drop")
hsp12
with(hsp12, mean_cpm[genotype == "WT" & time == "20min"] / mean_cpm[genotype == "WT" & time == "0min"])

# ---- Exercise 2: QC detective (6 min) ----------------------------------------
panel <- make_synthetic_panel(counts, coldata)        # SYN01..SYN08 appended
qc <- qc_metrics(panel$counts, panel$coldata)
qc |> select(sample, lib_size_M, genes_detected, pct_ncRNA, cor_to_group,
             noise_mad, marker_cpm, flags) |> print(n = Inf)

ggplot(qc, aes(cor_to_group, noise_mad, colour = synthetic)) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(data = filter(qc, flags != "PASS"), aes(label = sample),
                           size = 2.5, min.segment.length = 0) +
  scale_colour_manual(values = c(`FALSE` = "grey50", `TRUE` = "#D55E00")) +
  labs(x = "Median correlation with replicates", y = "Noise (MAD of log2 ratio)") +
  theme_pub()
# Q: Which failure is invisible to correlation/noise but caught by the genotype
#    marker?  A: SYN07_mixture (30% WT in a mip6Δ sample) -> MIP6 CPM > 0.
# Q: Which REAL sample gets flagged?  A: mip6_39_20_4 (excess noise, low depth).

# ---- Exercise 3: sequencing depth (5 min) -------------------------------------
s1 <- "wt_30_0_1"; s2 <- "wt_30_0_3"
pc <- classify_features(rownames(counts)) == "protein_coding"
full <- counts[pc, c(s1, s2)]
low  <- thin_matrix(full, 0.05)                 # 5% of reads ~ 1 M per library
colSums(full >= 10); colSums(low >= 10)

rep_df <- bind_rows(
  tibble(depth = "100% (~20 M)", a = log_cpm(full)[, 1], b = log_cpm(full)[, 2]),
  tibble(depth = "5% (~1 M)",    a = log_cpm(low)[, 1],  b = log_cpm(low)[, 2])
)
ggplot(rep_df, aes(a, b)) +
  geom_jitter(size = 0.3, alpha = 0.3, width = 0.05, height = 0.05) +   # jitter un-stacks ties at low counts
  geom_abline(colour = "red", linewidth = 0.3) +
  facet_wrap(~ depth) +
  labs(x = paste(s1, "log2 CPM"), y = paste(s2, "log2 CPM"),
       title = "Low depth inflates noise for weakly expressed genes") +
  theme_pub()
rep_df |> group_by(depth) |> summarise(r = cor(a, b))

# ---- Exercise 4: normalization (5 min) ----------------------------------------
cnt <- counts[pc, ]
cnt <- cnt[rowSums(cnt >= 10) >= 4, ]
sf <- tibble(sample = colnames(cnt),
             TC  = size_factors(cnt, "TC"),
             MoR = size_factors(cnt, "MoR")) |>
  left_join(coldata, by = "sample")
ggplot(sf, aes(TC, MoR, colour = time)) +
  geom_abline(linetype = 2) + geom_point(size = 2) +
  scale_colour_manual(values = pal_time) + theme_pub()

ma_df <- map_dfr(c("TC", "MoR"), \(m) {
  x <- normalize_counts(cnt, m)
  w0  <- rowMeans(x[, coldata$condition == "WT_0min"])
  w20 <- rowMeans(x[, coldata$condition == "WT_20min"])
  tibble(method = m, A = (w0 + w20) / 2, M = w20 - w0)
})
ma_df |> group_by(method) |> summarise(median_M = median(M))   # TC about +0.29, MoR about 0
ggplot(ma_df, aes(A, M)) +
  geom_point(size = 0.3, alpha = 0.3) +
  geom_hline(yintercept = 0, colour = "red") +
  geom_hline(data = ma_df |> group_by(method) |> summarise(m = median(M)),
             aes(yintercept = m), colour = "blue", linetype = 2) +
  facet_wrap(~ method) + theme_pub()

# ---- Exercise 5: DE + volcano, then add ONE bad replicate (7 min) -------------
keep_s <- coldata$condition %in% c("WT_0min", "WT_20min", "mip6_20min")
cd5 <- droplevels(as.data.frame(coldata[keep_s, ]))
rownames(cd5) <- cd5$sample
cnt5 <- counts[pc, cd5$sample]
cnt5 <- cnt5[edgeR::filterByExpr(cnt5, group = cd5$condition), ]

run5 <- function(m) {
  dds <- DESeq(DESeqDataSetFromMatrix(m, cd5, ~ batch + condition), quiet = TRUE)
  results(dds, contrast = c("condition", "mip6_20min", "WT_20min"), alpha = 0.05) |>
    as.data.frame() |> rownames_to_column("gene_id") |> as_tibble() |>
    mutate(symbol = gene_symbols(gene_id))
}
res_clean <- run5(cnt5)
plot_volcano(res_clean, lfc_cut = 0.5, title = "mip6Δ vs WT, 20 min")

cnt_bad <- cnt5
cnt_bad[, "mip6_39_20_1"] <- sim_noise(cnt5[, "mip6_39_20_1"], sdlog = 0.6)
res_bad <- run5(cnt_bad)
c(clean = sum(res_clean$padj < 0.05, na.rm = TRUE),
  with_bad_replicate = sum(res_bad$padj < 0.05, na.rm = TRUE))
# Extra: how many genes lost their p-value to Cook's distance outlier filtering?
sum(is.na(res_bad$pvalue) & res_bad$baseMean > 0)

# ---- Exercise 6: a publication-ready gene plot (2 min) ------------------------
g <- "TPS1"                                          # trehalose-6-P synthase
annot <- tibble(gene_id = rownames(counts), symbol = gene_symbols(rownames(counts)))
gid <- annot$gene_id[annot$symbol == g][1]
nc <- t(t(counts) / size_factors(counts[pc, ], "MoR"))
p_gene <- tibble(sample = colnames(nc), value = nc[gid, ]) |>
  left_join(coldata, by = "sample") |>
  ggplot(aes(time, value, colour = genotype)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6), size = 1.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2,
               position = position_dodge(0.6), linewidth = 0.4) +
  stat_summary(fun = mean, geom = "point", shape = 95, size = 6,
               position = position_dodge(0.6)) +
  scale_colour_manual(values = pal_genotype) +
  labs(x = NULL, y = "Normalized count", title = g) +
  theme_pub()
p_gene
save_pub(p_gene, "figures/my_gene_TPS1.pdf", width_mm = 60, height_mm = 50)
