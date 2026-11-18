# =============================================================================
# 03_exercises.R — 30-minute hands-on
#
# Fill in every ___ and run each block. Stuck? Answers: 03_exercises_solutions.R
# Work in pairs; one types, one reads the hints. Swap at every exercise.
#
#   Ex1  dplyr warm-up ................................ 5 min
#   Ex2  QC detective: find the 8 synthetic failures .. 6 min
#   Ex3  Sequencing depth and replicate scatter ........ 5 min
#   Ex4  Normalization: total count vs median-of-ratios  5 min
#   Ex5  DE + volcano, then ONE bad replicate .......... 7 min
#   Ex6  Publication-ready gene plot ................... 2 min
#   Then: shiny::runApp() and explore your own gene
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
# 1a. Library size (millions) per sample, smallest first.
lib <- counts_long |>
  group_by(___) |>
  summarise(lib_M = ___ / 1e6) |>
  arrange(___)
lib
# Q: Which sample is the shallowest? By how much?

# 1b. % of reads in the 50 most-abundant features of each sample
#     Hint: slice_max(count, n = 50, with_ties = FALSE)
counts_long |>
  group_by(sample) |>
  ___ |>
  summarise(top50 = sum(count)) |>
  left_join(lib, by = "sample") |>
  mutate(pct_top50 = ___) |>
  arrange(desc(pct_top50))

# 1c. HSP12 (YFL014W): mean CPM per genotype x time.
hsp12 <- counts_long |>
  group_by(sample) |> mutate(cpm = 1e6 * count / sum(count)) |> ungroup() |>
  filter(gene_id == ___) |>
  group_by(___, ___) |>
  summarise(mean_cpm = mean(cpm), .groups = "drop")
hsp12
# Q: Roughly what is the fold change at 20 min in WT?

# ---- Exercise 2: QC detective (6 min) ----------------------------------------
# make_synthetic_panel() appends 8 simulated failures (SYN01..SYN08) built from
# real counts: low depth, low complexity, gDNA, ncRNA, degradation*, noise,
# contamination, sample swap.  (*needs the TxDb package)
panel <- make_synthetic_panel(counts, coldata)
qc <- qc_metrics(panel$counts, panel$coldata)
qc |> select(sample, lib_size_M, genes_detected, pct_ncRNA, cor_to_group,
             noise_mad, marker_cpm, flags) |> print(n = Inf)

# Plot correlation-to-replicates vs noise, colour by synthetic, label flagged
ggplot(qc, aes(x = ___, y = ___, colour = synthetic)) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(data = filter(qc, flags != "PASS"), aes(label = sample),
                           size = 2.5, min.segment.length = 0) +
  theme_pub()
# Q1: Which failure is invisible on this plot but caught by the MIP6 marker?
# Q2: One REAL sample is flagged too. Which one, and why might that be?
# Q3: SYN02 has a NORMAL library size. What gives it away?

# ---- Exercise 3: sequencing depth (5 min) -------------------------------------
s1 <- "wt_30_0_1"; s2 <- "wt_30_0_3"
pc <- classify_features(rownames(counts)) == "protein_coding"
full <- counts[pc, c(s1, s2)]
low  <- thin_matrix(full, ___)                 # keep 5% of reads
colSums(full >= 10); colSums(low >= 10)

rep_df <- bind_rows(
  tibble(depth = "100%", a = log_cpm(full)[, 1], b = log_cpm(full)[, 2]),
  tibble(depth = "5%",   a = log_cpm(low)[, 1],  b = log_cpm(low)[, 2])
)
ggplot(rep_df, aes(a, b)) +
  ___(size = 0.3, alpha = 0.3, width = 0.05, height = 0.05) +   # un-stack tied low counts
  geom_abline(colour = "red") +
  facet_wrap(~ ___) +
  theme_pub()
# Q: Where on the plot does the noise appear: low or high expression? Why?

# ---- Exercise 4: normalization (5 min) ----------------------------------------
cnt <- counts[pc, ]
cnt <- cnt[rowSums(cnt >= 10) >= 4, ]
sf <- tibble(sample = colnames(cnt),
             TC  = size_factors(cnt, "TC"),
             MoR = size_factors(cnt, ___)) |>
  left_join(coldata, by = "sample")
ggplot(sf, aes(TC, MoR, colour = time)) +
  geom_abline(linetype = 2) + geom_point(size = 2) +
  scale_colour_manual(values = pal_time) + theme_pub()

ma_df <- map_dfr(c("TC", "MoR"), \(m) {
  x <- normalize_counts(cnt, m)
  w0  <- rowMeans(x[, coldata$condition == "WT_0min"])
  w20 <- rowMeans(x[, coldata$condition == ___])
  tibble(method = m, A = (w0 + w20) / 2, M = w20 - w0)
})
ma_df |> group_by(method) |> summarise(median_M = ___)
# Q: Under "most genes don't change", median M should be ~0. Which method
#    passes? Why do a handful of heat-shock genes bias total-count scaling?

# ---- Exercise 5: DE + volcano, then ONE bad replicate (7 min) -------------------
keep_s <- coldata$condition %in% c("WT_0min", "WT_20min", "mip6_20min")
cd5 <- droplevels(as.data.frame(coldata[keep_s, ]))
rownames(cd5) <- cd5$sample
cnt5 <- counts[pc, cd5$sample]
cnt5 <- cnt5[edgeR::filterByExpr(cnt5, group = cd5$condition), ]

run5 <- function(m) {
  dds <- DESeq(DESeqDataSetFromMatrix(m, cd5, ~ batch + condition), quiet = TRUE)
  results(dds, contrast = c("condition", ___, ___), alpha = 0.05) |>   # mip6 vs WT at 20 min
    as.data.frame() |> rownames_to_column("gene_id") |> as_tibble() |>
    mutate(symbol = gene_symbols(gene_id))
}
res_clean <- run5(cnt5)
plot_volcano(res_clean, lfc_cut = 0.5)

# Now corrupt ONE replicate with gene-wise noise and re-run
cnt_bad <- cnt5
cnt_bad[, "mip6_39_20_1"] <- sim_noise(cnt5[, "mip6_39_20_1"], sdlog = ___)
res_bad <- run5(cnt_bad)
c(clean = sum(res_clean$padj < 0.05, na.rm = TRUE),
  with_bad_replicate = sum(res_bad$padj < 0.05, na.rm = TRUE))
# Q: How many DE genes did one bad sample cost you? Try sdlog = 0.3 and 1.0.

# ---- Exercise 6: publication-ready gene plot (2 min) ----------------------------
g <- "___"                       # try TPS1, HSP104, RPL3, CTT1 or your favourite
annot <- tibble(gene_id = rownames(counts), symbol = gene_symbols(rownames(counts)))
gid <- annot$gene_id[annot$symbol == g][1]
nc <- t(t(counts) / size_factors(counts[pc, ], "MoR"))
p_gene <- tibble(sample = colnames(nc), value = nc[gid, ]) |>
  left_join(coldata, by = "sample") |>
  ggplot(aes(time, value, colour = genotype)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6), size = 1.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2,
               position = position_dodge(0.6)) +
  scale_colour_manual(values = pal_genotype) +
  labs(x = NULL, y = "Normalized count", title = g) +
  theme_pub()
p_gene
save_pub(p_gene, paste0("figures/my_gene_", g, ".pdf"), width_mm = ___, height_mm = 50)

# Done early?  shiny::runApp()  -> Gene explorer tab -> type your gene
