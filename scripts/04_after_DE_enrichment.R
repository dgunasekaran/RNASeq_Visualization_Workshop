# =============================================================================
# 04_after_DE_enrichment.R — what to do AFTER differential expression
# (post-workshop; needs the optional packages in scripts/install_workshop.R and an
#  internet connection for KEGG)
#
#   1. Over-representation (ORA): GO Biological Process, with the RIGHT universe
#   2. KEGG pathways (organism code "sce"; yeast KEGG IDs are ORF names)
#   3. Rank-based GSEA (no arbitrary cut-off) on the DESeq2 Wald statistic
#   4. compareCluster: enrichment of each temporal trajectory cluster
#   5. Redundancy reduction + publication plots (dotplot, enrichment map)
# Run scripts/02_demo_walkthrough.R first (it writes results/*.csv).
# =============================================================================

source("scripts/workshop_helpers.R")
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(org.Sc.sgd.db)
})
dir.create("figures/04_enrichment", recursive = TRUE, showWarnings = FALSE)

res <- read_csv("results/DE_WT_20min_vs_0min.csv", show_col_types = FALSE)
universe <- res$gene_id[!is.na(res$padj)]            # genes actually TESTED
up   <- res$gene_id[which(res$padj < 0.05 & res$log2FoldChange >  1)]
down <- res$gene_id[which(res$padj < 0.05 & res$log2FoldChange < -1)]

# ---- 1. GO over-representation ------------------------------------------------
# Always pass the tested genes as the universe: using the whole genome inflates
# significance (expressed genes are already biased towards some functions).
ego_up <- enrichGO(up, universe = universe, OrgDb = org.Sc.sgd.db, keyType = "ORF",
                   ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
ego_up <- simplify(ego_up, cutoff = 0.6)              # merge near-duplicate GO terms
head(as.data.frame(ego_up), 10)

# ---- 2. KEGG ------------------------------------------------------------------------
ekegg_down <- enrichKEGG(down, organism = "sce", universe = universe, keyType = "kegg")
head(as.data.frame(ekegg_down), 10)                   # expect Ribosome among top terms

# ---- 3. GSEA on the ranked list ----------------------------------------------------
ranks <- res |> filter(!is.na(stat)) |> arrange(desc(stat))
geneList <- setNames(ranks$stat, ranks$gene_id)
gsea_go <- gseGO(geneList, OrgDb = org.Sc.sgd.db, keyType = "ORF", ont = "BP",
                 minGSSize = 15, maxGSSize = 500, eps = 0, seed = TRUE)
p_gsea <- dotplot(gsea_go, showCategory = 10, split = ".sign") + facet_grid(. ~ .sign)

# ---- 4. Enrichment per trajectory cluster (compareCluster) -------------------------
# Save km$cluster from script 02 section 6d, e.g.:
#   write_csv(tibble(gene_id = names(km$cluster), cluster = km$cluster), "results/traj_clusters.csv")
if (file.exists("results/traj_clusters.csv")) {
  cl <- read_csv("results/traj_clusters.csv", show_col_types = FALSE)
  cc <- compareCluster(gene_id ~ cluster, data = cl, fun = "enrichGO",
                       OrgDb = org.Sc.sgd.db, keyType = "ORF", ont = "BP", universe = universe)
  save_pub(dotplot(cc, showCategory = 4) + theme_pub(),
           "figures/04_enrichment/compareCluster_trajectories.pdf", 183, 150)
}

# ---- 5. Figures ---------------------------------------------------------------------
save_pub(dotplot(ego_up, showCategory = 12) + theme_pub() + ggtitle("GO BP: up at 20 min (WT)"),
         "figures/04_enrichment/GO_up_dotplot.pdf", 140, 110)
save_pub(emapplot(pairwise_termsim(ego_up), showCategory = 25),
         "figures/04_enrichment/GO_up_emap.pdf", 183, 160)
save_pub(p_gsea, "figures/04_enrichment/GSEA_GO_dotplot.pdf", 183, 120)
save_pub(gseaplot2(gsea_go, geneSetID = 1:3, pvalue_table = TRUE),
         "figures/04_enrichment/GSEA_running_score.pdf", 183, 120)
