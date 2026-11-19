# =============================================================================
# RNA-seq QC Explorer
# Upload a count matrix + sample metadata; get QC, PCA, normalization checks,
# DESeq2 differential expression, heatmaps and gene plots, each with a
# plain-language interpretation of what the data and the model choices mean.
#
#   shiny::runApp()                               (from this folder)
#   shiny::runGitHub("<user>/rnaseq-qc-explorer")
# =============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(plotly)
  library(DT)
  library(DESeq2)
})
options(shiny.maxRequestSize = 300 * 1024^2)

DEMO_COUNTS <- "data/demo_counts.csv"
DEMO_META   <- "data/demo_metadata.csv"

interp_card <- function(id, title = "Interpretation") {
  card(class = "border-info", card_header(class = "bg-info-subtle", icon_text(title)), uiOutput(id))
}
icon_text <- function(x) HTML(paste0("&#9432;&nbsp;", x))

# ---- UI ------------------------------------------------------------------------
ui <- page_navbar(
  title = "RNA-seq QC Explorer",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  fillable = FALSE,
  sidebar = sidebar(
    width = 330,
    selectizeInput("genes", "Find genes (highlighted in every tab)", choices = NULL, multiple = TRUE,
                   options = list(placeholder = "Gene name or ID…", maxItems = 8)),
    accordion(
      open = c("Data", "Design"),
      accordion_panel("Data",
        radioButtons("source", NULL, c("Upload my data" = "upload",
                                       "Demo: yeast heat shock (GSE135568)" = "demo"), selected = "demo"),
        conditionalPanel("input.source == 'upload'",
          fileInput("f_counts", "Count matrix (csv / tsv / xlsx)", accept = c(".csv", ".tsv", ".txt", ".xlsx")),
          fileInput("f_meta", "Sample metadata (csv / tsv / xlsx)", accept = c(".csv", ".tsv", ".txt", ".xlsx")),
          fileInput("f_annot", "Optional: gene annotation (gene_id, symbol)", accept = c(".csv", ".tsv", ".txt")))),
      accordion_panel("Design",
        selectInput("group_col", "Groups to compare (metadata column)", choices = NULL),
        selectInput("cov_col", "Covariate to adjust for (e.g. batch)", choices = "none"),
        selectizeInput("drop", "Exclude samples", choices = NULL, multiple = TRUE,
                       options = list(placeholder = "None excluded")),
        radioButtons("prefilter", "Gene pre-filter",
                     c("edgeR::filterByExpr (recommended)" = "fbe",
                       "≥ 10 counts in ≥ smallest-group samples" = "min10"))),
      accordion_panel("Differential expression",
        selectInput("A", "Numerator group (A)", NULL),
        selectInput("B", "Reference group (B)", NULL),
        checkboxInput("use_cov", "Include covariate in the model", TRUE),
        numericInput("alpha", "Adjusted p-value <", 0.05, 0.001, 0.25, 0.005),
        numericInput("lfc", "|log2 fold change| ≥", 1, 0, 5, 0.25),
        selectInput("thr_mode", "Fold-change cut-off", c("Filter after testing" = "posthoc",
                                                         "Test against threshold" = "tested")),
        selectInput("shrink", "Fold-change estimation", c("ashr shrinkage" = "ashr", "apeglm shrinkage" = "apeglm",
                                                          "normal-prior shrinkage" = "normal", "None (MLE)" = "none")),
        selectInput("padj", "P-value adjustment", c("Benjamini–Hochberg (FDR)" = "BH", "IHW (weighted FDR)" = "IHW",
                                                    "Benjamini–Yekutieli" = "BY", "Bonferroni" = "bonferroni", "Holm" = "holm")),
        checkboxInput("indep", "Independent filtering", TRUE),
        checkboxInput("cooks", "Cook's distance outlier flagging", TRUE),
        actionButton("run", "Run / update model", class = "btn-primary w-100")),
      accordion_panel("Marker-gene check (optional)",
        helpText("E.g. a knocked-out gene, a sex-specific gene or a transgene. Flags samples where it contradicts the label."),
        selectizeInput("marker", "Marker gene", choices = NULL, options = list(placeholder = "none")),
        selectizeInput("marker_absent", "Groups where it should be ABSENT", choices = NULL, multiple = TRUE)),
      accordion_panel("Colours",
        selectInput("pal", "Groups", c("Okabe–Ito (colour-blind safe)" = "okabe", "viridis" = "viridis",
                                       "cividis" = "cividis", "Greyscale (print)" = "grey")),
        selectInput("hpal", "Heatmaps", c("Purple–Orange" = "PuOr", "Blue–Orange" = "BuOr")),
        helpText("Groups are always shown by shape as well as colour.")))
  ),

  # ------------------------------------------------------------ Overview
  nav_panel("Overview",
    layout_columns(
      value_box("Samples", textOutput("vb_n")), value_box("Features", textOutput("vb_g")),
      value_box("Median depth", textOutput("vb_d")), value_box("Samples flagged in QC", textOutput("vb_f"), theme = "warning")),
    uiOutput("data_alerts"),
    layout_columns(col_widths = c(7, 5),
      card(card_header("Sample metadata"), DTOutput("t_meta")),
      card(card_header("Design check: samples per group × covariate"), tableOutput("t_cross"),
           uiOutput("design_note"))),
    card(card_header("Input format"), markdown("
**Count matrix**: first column = unique gene IDs, then one column per sample with **raw integer counts** (from featureCounts, HTSeq, or Salmon via tximport). Not CPM/TPM/FPKM.
**Metadata**: one row per sample; a `sample` column (or the first column) matching the count-matrix column names, plus any columns describing samples (condition, genotype, time, batch…).
**Annotation** (optional): `gene_id, symbol` to show readable names. Yeast ORF IDs are translated automatically."))
  ),

  # ------------------------------------------------------------ QC
  nav_panel("QC",
    interp_card("qc_interp", "What the QC says"),
    card(card_header("Per-sample metrics and flags", downloadButton("dl_qc", "CSV", class = "btn-sm float-end")),
         DTOutput("t_qc")),
    card(full_screen = TRUE, card_header("QC dashboard"), plotOutput("p_qc", height = "520px")),
    layout_columns(col_widths = c(6, 6),
      card(full_screen = TRUE, card_header("Read composition by feature class"), plotlyOutput("p_comp", height = "420px")),
      card(full_screen = TRUE, card_header("Sample correlation (log CPM, expressed genes)"), plotlyOutput("p_cor", height = "420px"))),
    card(card_header("How to read each metric"), markdown("
| Metric | What it measures | Watch for |
|---|---|---|
| reads_M | Counted reads (millions) | < 25% of the median |
| genes_detected / genes_10plus | Features with ≥ 1 / ≥ 10 reads | Clearly fewer than replicates |
| pct_top50 | Share of reads in the 50 most abundant features | Very high = low diversity or a dominant transcript |
| pct_noncoding / pct_mito | Reads on rRNA/tRNA/ncRNA / mitochondrial genes | Far above the other samples |
| reads_on_silent_genes | Genes with reads here but zero in all other samples | Many = genomic DNA contamination |
| r_to_replicates | Median Pearson r with same-group samples | Low = outlier, failed library or swap |
| r_margin | Own-group r minus best other-group r | Strongly negative = looks like another group |
| noise_mad | Median absolute log2 deviation from replicates (expressed genes) | High = noisy library |
| saturation_gain | Genes gained going from 25% to 100% of reads | Near 0 with normal depth = low complexity |
| marker_cpm | Expression of the marker gene you chose | Contradicts the label |
"))
  ),

  # ------------------------------------------------------------ PCA
  nav_panel("PCA",
    layout_columns(col_widths = c(3, 9),
      card(
        sliderInput("ntop", "Most variable genes", 100, 5000, 500, 100),
        selectInput("pcx", "x axis", paste0("PC", 1:6), "PC1"),
        selectInput("pcy", "y axis", paste0("PC", 1:6), "PC2"),
        selectInput("pca_col", "Colour / shape by", choices = "group"),
        checkboxInput("pca_rmcov", "Remove covariate before PCA (display only)", FALSE),
        downloadButton("dl_pca", "PDF")),
      card(full_screen = TRUE, plotlyOutput("p_pca", height = "480px"))),
    interp_card("pca_interp", "What the PCA says"),
    layout_columns(col_widths = c(4, 4, 4),
      card(card_header("Variance explained"), plotOutput("p_scree", height = "280px")),
      card(card_header("Which metadata explains each PC (R²)"), tableOutput("t_assoc")),
      card(card_header("Top loading genes"), DTOutput("t_load")))
  ),

  # ------------------------------------------------------------ Normalization
  nav_panel("Normalization",
    layout_columns(col_widths = c(3, 9),
      card(checkboxGroupInput("norm_methods", "Methods",
             c("none", "TC", "UQ", "TMM", "MoR", "quantile", "VST"), selected = c("none", "TC", "TMM", "MoR")),
           helpText("TC = total count (CPM) · UQ = upper quartile · TMM = trimmed mean of M (edgeR) ·",
                    "MoR = median-of-ratios (DESeq2, used for DE) · VST = variance-stabilizing transform")),
      card(full_screen = TRUE, card_header("Relative log expression (RLE)"), plotOutput("p_rle", height = "520px"))),
    interp_card("norm_interp", "What the normalization check says"),
    card(card_header("Diagnostics for groups A vs B (sidebar)"), tableOutput("t_norm"))
  ),

  # ------------------------------------------------------------ DE
  nav_panel("Differential expression",
    uiOutput("de_headline"),
    interp_card("model_interp", "What the model and your choices mean"),
    checkboxInput("show_stats", strong("Show detailed statistics & methods paragraph"), FALSE),
    conditionalPanel("input.show_stats",
      layout_columns(col_widths = c(4, 4, 4),
        card(card_header("What happened to each gene"), tableOutput("t_flow")),
        card(card_header("Up / down at other cut-offs"), tableOutput("t_grid")),
        card(card_header("Effect sizes of called genes"), tableOutput("t_eff"))),
      card(card_header("Methods paragraph (edit before publishing)",
                       downloadButton("dl_methods", "TXT", class = "btn-sm float-end")),
           verbatimTextOutput("methods", placeholder = TRUE))),
    layout_columns(col_widths = c(7, 5),
      card(full_screen = TRUE, card_header("Volcano plot", downloadButton("dl_volcano", "PDF", class = "btn-sm float-end")),
           plotlyOutput("p_volcano", height = "520px")),
      card(full_screen = TRUE, card_header("MA plot"), plotlyOutput("p_ma", height = "250px"),
           card_header("p-value histogram"), plotOutput("p_phist", height = "210px"))),
    card(card_header("Results (searchable)", downloadButton("dl_de", "CSV", class = "btn-sm float-end")), DTOutput("t_de"))
  ),

  # ------------------------------------------------------------ Heatmap
  nav_panel("Heatmap",
    layout_columns(col_widths = c(3, 9),
      card(
        radioButtons("hm_src", "Genes", c("Top DE (current contrast)" = "de", "Most variable" = "var",
                                         "Searched genes + list" = "list")),
        sliderInput("hm_n", "Number of genes", 10, 200, 50, 5),
        textAreaInput("hm_list", "Gene list (names or IDs)", "", rows = 3),
        checkboxInput("hm_z", "Row z-score", TRUE),
        checkboxInput("hm_cov", "Remove covariate for display", TRUE),
        downloadButton("dl_hm", "PDF")),
      card(full_screen = TRUE, plotlyOutput("p_hm", height = "640px")))
  ),

  # ------------------------------------------------------------ Gene explorer
  nav_panel("Gene explorer",
    layout_columns(col_widths = c(3, 9),
      card(
        selectInput("gx_x", "x axis", choices = "group"),
        selectInput("gx_col", "Colour / shape by", choices = "group"),
        radioButtons("gp_type", "Plot", c("Points + mean ± SE" = "pts", "Box + points" = "box",
                                          "Violin + points" = "vio")),
        checkboxInput("gp_log", "log10 y axis", FALSE),
        numericInput("gp_w", "Export width (mm)", 89, 40, 183),
        numericInput("gp_h", "Export height (mm)", 70, 30, 200),
        downloadButton("dl_gp_pdf", "PDF"), downloadButton("dl_gp_png", "PNG")),
      card(full_screen = TRUE, plotOutput("p_gene", height = "520px"))),
    card(card_header("Statistics for these genes (current contrast)"), DTOutput("t_gene"))
  ),

  # ------------------------------------------------------------ Response map
  nav_panel("Response map",
    layout_columns(col_widths = c(3, 9),
      card(
        selectInput("rm_A1", "Contrast 1: A", NULL), selectInput("rm_B1", "Contrast 1: B", NULL),
        selectInput("rm_A2", "Contrast 2: A", NULL), selectInput("rm_B2", "Contrast 2: B", NULL),
        helpText("Each point is a gene. On the diagonal = same response in both contrasts;",
                 "off-diagonal = context-dependent response (e.g. a treatment effect that differs by genotype).")),
      card(full_screen = TRUE, plotlyOutput("p_rm", height = "560px"))),
    interp_card("rm_interp", "What the response map says")
  ),

  # ------------------------------------------------------------ About
  nav_panel("About", card(markdown("
### Methods and references
| Step | Method | Reference |
|---|---|---|
| Pre-filter | edgeR::filterByExpr | Chen, Lun & Smyth 2016, *F1000Research* 5:1438 |
| Normalization (DE) | Median-of-ratios | Anders & Huber 2010, *Genome Biol* 11:R106 |
| Other normalizations | TMM; upper quartile; quantile | Robinson & Oshlack 2010; Bullard et al. 2010; Bolstad et al. 2003 |
| Model | Negative-binomial GLM, Wald test, dispersion shrinkage, VST | Love, Huber & Anders 2014, *Genome Biol* 15:550 |
| Fold-change shrinkage | apeglm / ashr | Zhu, Ibrahim & Love 2019; Stephens 2017 |
| Multiple testing | Benjamini–Hochberg; IHW | Benjamini & Hochberg 1995; Ignatiadis et al. 2016 |
| Threshold test | H0: abs(log2FC) ≤ t | McCarthy & Smyth 2009, *Bioinformatics* 25:765 |
| Replicates | Power vs replicate number | Schurch et al. 2016, *RNA* 22:839 |
| Workflow | Exploratory analysis & DE | Love et al. 2015, *F1000Research* 4:1070 |

### Accessibility
Groups use the Okabe–Ito palette (or viridis, cividis, greys) plus a unique point shape each.
Up = vermillion, down = blue; heatmaps use purple–orange / blue–orange and cividis, all readable with red–green colour-vision deficiency.

### Interpretation text
Interpretations are generated from your data with simple, stated rules (robust z-scores beyond ±3, R² from linear models of PC scores, the π₀ estimate from p > 0.5).
They are prompts for judgement, not automatic decisions.

Demo data: Nuño-Cabanes et al. 2020, *Sci Data* 7:69 (GEO GSE135568).
")))
)

# ---- Server --------------------------------------------------------------------------
server <- function(input, output, session) {

  # ---------- input ----------
  raw <- reactive({
    if (input$source == "demo") {
      cnt <- read_counts(DEMO_COUNTS); cd <- read_metadata(DEMO_META, colnames(cnt)); ann <- NULL
    } else {
      validate(need(input$f_counts, "Upload a count matrix and a metadata file (sidebar ▸ Data)."),
               need(input$f_meta, "Upload the sample metadata file."))
      # Build the message only when reading failed: need() evaluates its message
      # eagerly, and conditionMessage() of a successfully read matrix is an error.
      cnt <- tryCatch(read_counts(input$f_counts$datapath), error = function(e) e)
      if (inherits(cnt, "error")) validate(need(FALSE, paste("Count matrix:", conditionMessage(cnt))))
      cd <- tryCatch(read_metadata(input$f_meta$datapath, colnames(cnt)), error = function(e) e)
      if (inherits(cd, "error")) validate(need(FALSE, paste("Metadata:", conditionMessage(cd))))
      ann <- if (!is.null(input$f_annot)) read_table_any(input$f_annot$datapath) else NULL
    }
    annot <- annotate_genes(rownames(cnt), ann)
    list(counts = cnt, meta = cd, annot = annot, fclass = classify_features(rownames(cnt), annot$symbol))
  })

  observe({
    r <- raw(); cd <- r$meta
    cols <- setdiff(names(cd), "sample")
    cand <- cols[vapply(cd[cols], function(x) { k <- length(unique(x)); k >= 2 && k < nrow(cd) }, logical(1))]
    updateSelectInput(session, "group_col", choices = cand,
                      selected = if ("condition" %in% cand) "condition" else cand[1])
    updateSelectInput(session, "cov_col", choices = c("none", cand),
                      selected = if ("batch" %in% cand) "batch" else "none")
    updateSelectizeInput(session, "drop", choices = cd$sample, selected = character(0))
    a <- r$annot
    ch <- stats::setNames(a$gene_id, ifelse(a$symbol == a$gene_id, a$gene_id, paste0(a$symbol, " (", a$gene_id, ")")))
    updateSelectizeInput(session, "genes", choices = ch, server = TRUE)
    updateSelectizeInput(session, "marker", choices = c(none = "", ch), server = TRUE,
                         selected = if (input$source == "demo") "YHR015W" else "")
    meta_cols <- c("group", cols)
    updateSelectInput(session, "pca_col", choices = meta_cols, selected = "group")
    updateSelectInput(session, "gx_x", choices = meta_cols, selected = if ("time" %in% cols) "time" else "group")
    updateSelectInput(session, "gx_col", choices = meta_cols, selected = if ("genotype" %in% cols) "genotype" else "group")
  })

  dat <- reactive({
    r <- raw(); req(input$group_col %in% names(r$meta))
    keep <- !(r$meta$sample %in% input$drop)
    validate(need(sum(keep) >= 3, "Keep at least 3 samples."))
    cd <- r$meta[keep, , drop = FALSE]
    cd$group <- droplevels(safe_factor(cd[[input$group_col]]))
    cd$covariate <- if (input$cov_col != "none") droplevels(factor(cd[[input$cov_col]])) else factor(rep("all", nrow(cd)))
    rownames(cd) <- cd$sample
    list(counts = r$counts[, keep, drop = FALSE], meta = cd, annot = r$annot, fclass = r$fclass)
  })

  observe({
    g <- dat()$meta$group; lv <- levels(g)
    sel <- function(pref, i) if (pref %in% lv) pref else lv[min(i, length(lv))]
    updateSelectInput(session, "A", choices = lv, selected = sel("WT_20min", 2))
    updateSelectInput(session, "B", choices = lv, selected = sel("WT_0min", 1))
    updateSelectInput(session, "rm_A1", choices = lv, selected = sel("WT_20min", 2))
    updateSelectInput(session, "rm_B1", choices = lv, selected = sel("WT_0min", 1))
    updateSelectInput(session, "rm_A2", choices = lv, selected = sel("mip6_20min", min(4, length(lv))))
    updateSelectInput(session, "rm_B2", choices = lv, selected = sel("mip6_0min", min(3, length(lv))))
    updateSelectizeInput(session, "marker_absent", choices = lv,
                         selected = if (input$source == "demo") intersect(lv, c("mip6_0min", "mip6_20min", "mip6_120min")) else character(0))
  })

  pal <- reactive({ g <- dat()$meta$group; group_palette(levels(g), input$pal) })
  shp <- reactive({ g <- dat()$meta$group; group_shapes(levels(g)) })
  pal_for <- function(x) { lv <- levels(factor(x)); list(col = group_palette(lv, input$pal), shp = group_shapes(lv)) }

  filtered <- reactive({
    d <- dat()
    keep <- if (input$prefilter == "fbe") edgeR::filterByExpr(d$counts, group = d$meta$group)
            else rowSums(d$counts >= 10) >= min(table(d$meta$group))
    d$counts[keep, , drop = FALSE]
  })

  vst_m <- reactive({
    cnt <- filtered(); cd <- dat()$meta; d0 <- DESeqDataSetFromMatrix(cnt, cd, ~ 1)
    SummarizedExperiment::assay(if (nrow(cnt) >= 1000) vst(d0, blind = TRUE) else varianceStabilizingTransformation(d0, blind = TRUE))
  })
  remove_cov <- function(v) {
    cd <- dat()$meta
    if (nlevels(cd$covariate) < 2) return(v)
    limma::removeBatchEffect(v, batch = cd$covariate, design = model.matrix(~ group, cd))
  }

  # ---------- Overview ----------
  # Reactives are always read into a local variable FIRST. DESeq2 makes colSums,
  # median, nrow, paste, order, ... Bioconductor (S4) generics; if a reactive
  # that is still waiting for input (req/validate) is evaluated inside one of
  # them, Shiny's "not ready yet" signal becomes a real error
  # ("error in evaluating the argument 'x' in selecting a method for ...").
  output$vb_n <- renderText({ cnt <- dat()$counts; ncol(cnt) })
  output$vb_g <- renderText({ cnt <- dat()$counts; format(nrow(cnt), big.mark = ",") })
  output$vb_d <- renderText({ cnt <- dat()$counts; paste0(round(median(colSums(cnt)) / 1e6, 1), " M reads") })
  output$vb_f <- renderText({ q <- qc(); sum(q$flags != "PASS") })
  output$t_meta <- renderDT({ m <- raw()$meta; datatable(m, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE)) })
  output$t_cross <- renderTable({
    cd <- dat()$meta
    as.data.frame.matrix(table(group = cd$group, covariate = cd$covariate)) |> rownames_to_column("group")
  }, striped = TRUE)
  output$design_note <- renderUI({
    cd <- dat()$meta; n <- table(cd$group)
    msg <- c(
      if (any(n < 2)) "<b>At least one group has a single sample:</b> differential expression needs replicates.",
      if (any(n == 2)) "Some groups have only 2 replicates: estimates will be unstable.",
      if (nlevels(cd$covariate) > 1 && any(rowSums(table(cd$group, cd$covariate) > 0) == 1))
        "Some groups appear in only one covariate level: the covariate is partly confounded with groups.")
    if (length(msg)) div(class = "alert alert-warning", HTML(paste(msg, collapse = "<br>")))
    else p(class = "text-success", "Every group has replicates", if (nlevels(cd$covariate) > 1) "and spans several covariate levels." else ".")
  })
  output$data_alerts <- renderUI({
    cnt <- raw()$counts
    if (isTRUE(attr(cnt, "non_integer")))
      div(class = "alert alert-danger", "The count matrix contains non-integer values and was rounded. DESeq2 needs RAW counts, not CPM/TPM/FPKM or normalized data.")
  })

  # ---------- QC ----------
  qc <- reactive({
    d <- dat()
    mk <- if (nzchar(input$marker %||% "")) input$marker else NULL
    qc_metrics(d$counts, d$meta$group, d$fclass, marker = mk,
               marker_absent_in = if (!is.null(mk) && length(input$marker_absent)) input$marker_absent else NULL)
  })
  output$qc_interp <- renderUI({ q <- qc(); markdown(qc_interpretation(q)) })
  output$t_qc <- renderDT({
    q <- qc() |> select(where(~ !all(is.na(.x))))
    num <- names(q)[vapply(q, is.numeric, logical(1))]
    datatable(q, rownames = FALSE, filter = "top", options = list(pageLength = 12, scrollX = TRUE)) |>
      formatRound(num, 3) |>
      formatStyle("flags", backgroundColor = styleEqual("PASS", "#E8F2FA", default = "#FCE9DD"))
  })
  output$dl_qc <- downloadHandler("qc_metrics.csv", function(f) { q <- qc(); readr::write_csv(q, f) })

  output$p_qc <- renderPlot({
    q <- qc()
    lab <- c(reads_M = "Reads (M)", genes_10plus = "Genes ≥ 10 reads", pct_top50 = "% reads in top 50",
             pct_noncoding = "% non-coding", pct_mito = "% mitochondrial", reads_on_silent_genes = "Reads on silent genes",
             r_to_replicates = "r to replicates", noise_mad = "Noise (MAD)", saturation_gain = "Saturation gain",
             marker_cpm = "Marker CPM")
    keep <- names(lab)[vapply(names(lab), function(k) k %in% names(q) && !all(is.na(q[[k]])), logical(1))]
    q |> select(sample, group, flags, all_of(keep)) |>
      pivot_longer(all_of(keep), names_to = "metric") |>
      mutate(metric = factor(lab[metric], levels = lab[keep]),
             sample = factor(sample, levels = rev(q$sample)),
             status = ifelse(flags == "PASS", "Pass", "Flagged")) |>
      ggplot(aes(value, sample, fill = status)) +
      geom_col(width = 0.75) +
      facet_wrap(~ metric, scales = "free_x", nrow = 1) +
      scale_fill_manual(values = c(Pass = "#9AA5B1", Flagged = "#D55E00")) +
      labs(x = NULL, y = NULL, fill = NULL) + theme_pub(10) + theme(legend.position = "bottom")
  })
  output$p_comp <- renderPlotly({
    cnt <- dat()$counts; fc <- dat()$fclass
    validate(need(length(unique(fc)) > 1, "Feature classes are available for yeast IDs, or mitochondrial genes when an annotation with MT- symbols is uploaded."))
    comp <- rowsum(cnt, fc); comp <- sweep(comp, 2, colSums(comp), "/") * 100
    df <- as.data.frame(comp) |> rownames_to_column("class") |>
      pivot_longer(-class, names_to = "sample", values_to = "pct") |>
      mutate(sample = factor(sample, levels = colnames(cnt)))
    p <- ggplot(df, aes(sample, pct, fill = class, text = paste0(sample, "<br>", class, ": ", round(pct, 2), "%"))) +
      geom_col() + scale_fill_manual(values = unname(OKABE_ITO)) + labs(x = NULL, y = "% reads") + theme_pub(10) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7))
    ggplotly(p, tooltip = "text")
  })
  output$p_cor <- renderPlotly({
    cnt <- filtered(); cm <- cor(log_cpm(cnt)); o <- hclust(as.dist(1 - cm))$order
    plot_ly(x = colnames(cm)[o], y = colnames(cm)[o], z = cm[o, o], type = "heatmap",
            colors = viridisLite::cividis(100), hovertemplate = "%{x}<br>%{y}<br>r = %{z:.4f}<extra></extra>") |>
      layout(xaxis = list(tickfont = list(size = 8)), yaxis = list(tickfont = list(size = 8)))
  })

  # ---------- PCA ----------
  pca <- reactive({
    v <- vst_m(); if (isTRUE(input$pca_rmcov)) v <- remove_cov(v)
    top <- order(matrixStats::rowVars(v), decreasing = TRUE)[seq_len(min(input$ntop, nrow(v)))]
    prcomp(t(v[top, ]))
  })
  assoc <- reactive({
    cd <- dat()$meta
    facs <- unique(c("group", if (nlevels(cd$covariate) > 1) "covariate",
                     setdiff(names(cd), c("sample", "group", "covariate", input$group_col, input$cov_col))))
    p <- pca(); pca_associations(p, cd, facs)
  })
  pca_gg <- reactive({
    p <- pca(); ve <- round(100 * p$sdev^2 / sum(p$sdev^2), 1)
    k <- min(6, ncol(p$x)); xi <- as.integer(sub("PC", "", input$pcx)); yi <- as.integer(sub("PC", "", input$pcy))
    validate(need(xi <= k && yi <= k, "Not enough samples for that PC."))
    cd <- dat()$meta; by <- if (input$pca_col %in% names(cd)) input$pca_col else "group"
    df <- as_tibble(p$x[, seq_len(k), drop = FALSE]) |> mutate(sample = rownames(p$x), by = factor(cd[[by]]))
    pp <- pal_for(df$by)
    ggplot(df, aes(.data[[input$pcx]], .data[[input$pcy]], colour = by, shape = by,
                   text = paste0(sample, "<br>", by))) +
      geom_point(size = 3.2) +
      scale_colour_manual(values = pp$col) + scale_shape_manual(values = pp$shp) +
      labs(x = paste0(input$pcx, " (", ve[xi], "%)"), y = paste0(input$pcy, " (", ve[yi], "%)"),
           colour = by, shape = by) + theme_pub()
  })
  output$p_pca <- renderPlotly({ g <- pca_gg(); ggplotly(g, tooltip = "text") })
  output$dl_pca <- downloadHandler("PCA.pdf", function(f) { g <- pca_gg(); save_plot(g, f, 120, 90) })
  output$p_scree <- renderPlot({
    p <- pca(); ve <- 100 * p$sdev^2 / sum(p$sdev^2)
    tibble(PC = factor(paste0("PC", seq_along(ve)), levels = paste0("PC", seq_along(ve))), ve = ve) |> head(10) |>
      ggplot(aes(PC, ve)) + geom_col(fill = "#0072B2") +
      geom_text(aes(label = sprintf("%.0f%%", ve)), vjust = -0.3, size = 3) + labs(x = NULL, y = "% variance") + theme_pub()
  })
  output$t_assoc <- renderTable({
    assoc() |> mutate(R2 = round(R2, 2)) |> pivot_wider(names_from = PC, values_from = R2)
  }, striped = TRUE, na = "–")
  output$pca_interp <- renderUI({ p <- pca(); a <- assoc(); cd <- dat()$meta; markdown(pca_interpretation(p, a, cd)) })
  output$t_load <- renderDT({
    p <- pca(); an <- dat()$annot
    tibble(gene_id = rownames(p$rotation), loading = p$rotation[, input$pcx]) |>
      left_join(an, by = "gene_id") |> arrange(desc(abs(loading))) |> head(50) |>
      datatable(rownames = FALSE, options = list(pageLength = 8, dom = "ftip")) |> formatRound("loading", 4)
  })

  # ---------- Normalization ----------
  norm_r <- reactive({
    req(input$norm_methods)
    cnt <- filtered()
    set_names(lapply(input$norm_methods, \(m) normalize_counts(cnt, m)), input$norm_methods)
  })
  norm_diag <- reactive({
    req(input$A, input$B); nr <- norm_r(); g <- dat()$meta$group
    norm_diagnostics(nr, g, input$A, input$B)
  })
  output$p_rle <- renderPlot({
    cd <- dat()$meta; nr <- norm_r(); pl <- pal()
    imap_dfr(nr, \(m, nm) rle_long(m) |> mutate(method = nm)) |>
      left_join(select(cd, sample, group), by = "sample") |>
      mutate(sample = factor(sample, levels = cd$sample), method = factor(method, levels = input$norm_methods)) |>
      ggplot(aes(sample, rle, fill = group)) +
      geom_boxplot(outlier.shape = NA, linewidth = 0.2) +
      geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
      coord_cartesian(ylim = c(-1.5, 1.5)) + facet_wrap(~ method, ncol = 2) +
      scale_fill_manual(values = pl) + labs(x = NULL, y = "Relative log expression", fill = "Group") +
      theme_pub(10) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))
  })
  output$t_norm <- renderTable(norm_diag(), digits = 3, striped = TRUE)
  output$norm_interp <- renderUI({ nd <- norm_diag(); markdown(norm_interpretation(nd, input$A, input$B)) })

  # ---------- Differential expression ----------
  dds <- eventReactive(list(input$run, filtered()), {
    cd <- dat()$meta; cnt <- filtered()
    validate(need(nlevels(cd$group) >= 2, "Need at least two groups."))
    use_cov <- isTRUE(input$use_cov) && nlevels(cd$covariate) > 1
    dropped <- FALSE
    if (use_cov) {
      mm <- model.matrix(~ covariate + group, cd)
      if (qr(mm)$rank < ncol(mm)) { use_cov <- FALSE; dropped <- TRUE }
    }
    withProgress(message = "Fitting DESeq2…", value = 0.4, {
      d <- DESeqDataSetFromMatrix(cnt, cd, if (use_cov) ~ covariate + group else ~ group)
      d <- DESeq(d, quiet = TRUE)
    })
    S4Vectors::metadata(d)$use_cov <- use_cov
    S4Vectors::metadata(d)$cov_dropped <- dropped
    d
  }, ignoreNULL = FALSE)

  de <- eventReactive(list(input$run, dds(), input$A, input$B), {
    req(input$A, input$B)
    validate(need(input$A != input$B, "Choose two different groups for A and B."))
    d <- dds()
    validate(need(all(c(input$A, input$B) %in% levels(d$group)), "Press 'Run / update model'."))
    x <- withProgress(message = "Testing and estimating fold changes…", value = 0.6,
      run_de(d, input$A, input$B, alpha = input$alpha, lfc = input$lfc, shrink = input$shrink,
             padj_method = input$padj, threshold_mode = input$thr_mode,
             indep_filter = input$indep, cooks = input$cooks))
    if (length(x$meta$notes)) showNotification(paste(x$meta$notes, collapse = "; "), type = "warning")
    an <- dat()$annot
    x$res <- left_join(x$res, an, by = "gene_id")
    x
  }, ignoreNULL = FALSE)

  design_txt <- reactive({
    d <- dds()
    if (isTRUE(S4Vectors::metadata(d)$use_cov)) sprintf("~ %s + %s (covariate + groups)", input$cov_col, input$group_col)
    else sprintf("~ %s", input$group_col)
  })

  output$de_headline <- renderUI({
    x <- de(); p <- x$params
    HTML(sprintf("<p class='lead'><b>%s vs %s</b>: <span style='color:#B34700'><b>%d up</b></span> and <span style='color:#005A8C'><b>%d down</b></span> of %s genes tested (adjusted p &lt; %s, |log2 FC| ≥ %s%s).</p>",
                 p$A, p$B, x$meta$n_up, x$meta$n_down, format(x$meta$n_tested, big.mark = ","), p$alpha, p$lfc,
                 if (p$threshold_mode == "tested") ", tested" else ""))
  })
  output$model_interp <- renderUI({
    a <- assoc(); cv <- a$R2[a$factor == "covariate" & a$PC %in% c("PC1", "PC2")]
    x <- de(); d <- dds(); dt <- design_txt()
    markdown(model_interpretation(x, d, dt, isTRUE(S4Vectors::metadata(d)$cov_dropped),
                                  TRUE, if (length(cv)) max(cv, na.rm = TRUE) else NA))
  })
  output$t_flow <- renderTable({
    m <- de()$meta; cnt <- dat()$counts
    data.frame(Step = c("Features in matrix", "Kept by pre-filter", "Zero counts", "Count outliers (Cook's)",
                        "Removed by independent filtering", "With adjusted p-value", "Up", "Down"),
               Genes = format(c(nrow(cnt), m$n_input, m$n_zero, m$n_cooks, m$n_indep,
                                m$n_tested - m$n_indep, m$n_up, m$n_down), big.mark = ","))
  }, striped = TRUE)
  output$t_grid <- renderTable({ r <- de()$res; de_grid(r) }, striped = TRUE, digits = 3)
  output$t_eff <- renderTable({
    r <- de()$res
    r |> filter(status != "NS") |> group_by(Direction = as.character(status)) |>
      summarise(Genes = n(), `Median log2FC` = median(log2FoldChange), `Max |log2FC|` = max(abs(log2FoldChange)),
                `Median mean count` = median(baseMean), .groups = "drop")
  }, striped = TRUE, digits = 2)
  methods_txt <- reactive({
    pf <- if (input$prefilter == "fbe") "Lowly expressed features were removed with edgeR::filterByExpr"
          else "Features with ≥ 10 counts in at least as many samples as the smallest group were kept"
    excl <- if (length(input$drop)) sprintf(" Samples excluded after QC: %s.", paste(input$drop, collapse = ", ")) else ""
    x <- de(); dt <- design_txt(); cnt <- dat()$counts
    paste0(methods_text(x$params, x$meta, dt, nrow(cnt), pf), excl)
  })
  output$methods <- renderText(methods_txt())
  output$dl_methods <- downloadHandler("methods.txt", function(f) { txt <- methods_txt(); writeLines(txt, f) })

  volcano_gg <- reactive({
    x <- de(); p <- x$params
    r <- filter(x$res, !is.na(padj)) |>
      mutate(capped = padj == 0, y = volcano_y(padj),
             padj_txt = ifelse(capped, "0 (below numerical precision; drawn at the top)",
                               as.character(signif(padj, 3))))
    n_cap <- sum(r$capped)
    hl <- filter(r, gene_id %in% input$genes)
    lab <- r |> filter(status != "NS") |> arrange(padj, desc(abs(log2FoldChange))) |> head(10)
    ggplot(r, aes(log2FoldChange, y, colour = status,
                  text = paste0(symbol, "<br>log2FC = ", round(log2FoldChange, 2), "<br>padj = ", padj_txt))) +
      geom_point(data = filter(r, !capped), size = 1, alpha = 0.75) +
      geom_point(data = filter(r, capped), shape = 17, size = 2.4, show.legend = FALSE) +
      geom_hline(yintercept = -log10(p$alpha), linetype = 2, linewidth = 0.3) +
      geom_vline(xintercept = c(-p$lfc, p$lfc), linetype = 2, linewidth = 0.3) +
      geom_point(data = hl, colour = "black", fill = "#F0E442", shape = 21, size = 3.2) +
      geom_text(data = bind_rows(lab, hl), aes(label = symbol), colour = "black", size = 3, vjust = -0.9, check_overlap = TRUE) +
      scale_colour_manual(values = DE_COLOURS, drop = FALSE) +
      # room above the top points and at the sides so gene labels are not clipped
      scale_x_continuous(expand = expansion(mult = 0.08)) +
      scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
      labs(x = "log2 fold change", y = "-log10 adjusted p", colour = NULL,
           title = paste(p$A, "vs", p$B), subtitle = capped_note(n_cap)) +
      theme_pub() + theme(legend.position = "bottom")
  })
  output$p_volcano <- renderPlotly({
    g <- volcano_gg()
    ttl <- g$labels$title                       # ggplotly drops subtitles: add it to the title
    if (!is.null(g$labels$subtitle)) ttl <- paste0(ttl, "<br><sup>", gsub("\n", " ", g$labels$subtitle), "</sup>")
    ggplotly(g, tooltip = "text") |> layout(title = list(text = ttl)) |> toWebGL()
  })
  output$dl_volcano <- downloadHandler(function() paste0("volcano_", input$A, "_vs_", input$B, ".pdf"),
                                       function(f) { g <- volcano_gg(); save_plot(g, f, 130, 120) })
  output$p_ma <- renderPlotly({
    res <- de()$res; r <- filter(res, !is.na(padj))
    p <- ggplot(r, aes(baseMean, log2FoldChange, colour = status,
                       text = paste0(symbol, "<br>mean = ", round(baseMean), "<br>log2FC = ", round(log2FoldChange, 2)))) +
      geom_point(size = 0.8, alpha = 0.7) + geom_hline(yintercept = 0, linewidth = 0.3) +
      scale_x_log10() + scale_colour_manual(values = DE_COLOURS, drop = FALSE) +
      labs(x = "Mean normalized count", y = "log2 FC", colour = NULL) + theme_pub(10)
    ggplotly(p, tooltip = "text") |> toWebGL()
  })
  output$p_phist <- renderPlot({
    res <- de()$res
    ggplot(filter(res, !is.na(pvalue), baseMean > 1), aes(pvalue)) +
      geom_histogram(breaks = seq(0, 1, 0.025), fill = "#56B4E9", colour = "white", linewidth = 0.1) +
      labs(x = "Raw p-value", y = "Genes") + theme_pub(10)
  })
  output$t_de <- renderDT({
    res <- de()$res
    r <- res |> select(symbol, gene_id, baseMean, log2FoldChange, log2FC_MLE, lfcSE, pvalue, padj, status)
    pcols <- which(names(r) %in% c("pvalue", "padj")) - 1          # 0-based for DT
    datatable(r, rownames = FALSE, filter = "top",
              options = list(pageLength = 15, scrollX = TRUE, order = list(list(7, "asc")),
                             columnDefs = list(list(targets = pcols, render = pvalue_render())))) |>
      formatRound(c("baseMean", "log2FoldChange", "log2FC_MLE", "lfcSE"), 2)
  })
  output$dl_de <- downloadHandler(function() paste0("DE_", input$A, "_vs_", input$B, ".csv"),
                                  function(f) { res <- de()$res; readr::write_csv(res, f) })

  # ---------- Heatmap ----------
  hm <- reactive({
    v <- vst_m(); if (isTRUE(input$hm_cov)) v <- remove_cov(v)
    a <- dat()$annot
    genes <- switch(input$hm_src,
      de   = { res <- de()$res; res |> filter(status != "NS") |> arrange(padj) |> head(input$hm_n) |> pull(gene_id) },
      var  = rownames(v)[order(matrixStats::rowVars(v), decreasing = TRUE)][seq_len(min(input$hm_n, nrow(v)))],
      list = { tk <- strsplit(trimws(input$hm_list), "[,;[:space:]]+")[[1]]
               unique(c(input$genes, a$gene_id[a$gene_id %in% tk | toupper(a$symbol) %in% toupper(tk)])) })
    genes <- intersect(genes, rownames(v))
    validate(need(length(genes) >= 2, "Fewer than 2 genes to show (no DE genes, or names not found)."))
    g <- dat()$meta$group; ord <- order(g)
    m <- v[genes, ord, drop = FALSE]
    if (isTRUE(input$hm_z)) { m <- t(scale(t(m))); m[is.na(m)] <- 0 }
    rownames(m) <- make.unique(a$symbol[match(rownames(m), a$gene_id)])
    m[hclust(dist(m))$order, , drop = FALSE]
  })
  output$p_hm <- renderPlotly({
    m <- hm(); z <- isTRUE(input$hm_z)
    lim <- if (z) c(-2.5, 2.5) else range(m)
    plot_ly(x = colnames(m), y = rownames(m), z = pmin(pmax(m, lim[1]), lim[2]), type = "heatmap",
            zmin = lim[1], zmax = lim[2], colors = if (z) heat_palette(input$hpal) else viridisLite::cividis(100),
            colorbar = list(title = if (z) "z-score" else "VST"),
            hovertemplate = "%{y}<br>%{x}<br>%{z:.2f}<extra></extra>") |>
      layout(xaxis = list(tickfont = list(size = 8), categoryorder = "array", categoryarray = colnames(m)),
             yaxis = list(tickfont = list(size = 8), autorange = "reversed"))
  })
  output$dl_hm <- downloadHandler("heatmap.pdf", function(f) {
    m <- hm(); cd <- dat()$meta; pl <- pal(); z <- isTRUE(input$hm_z)
    pheatmap::pheatmap(m, cluster_rows = FALSE, cluster_cols = FALSE,
      annotation_col = data.frame(group = cd$group, row.names = cd$sample), annotation_colors = list(group = pl),
      color = if (z) heat_palette(input$hpal) else viridisLite::cividis(101),
      breaks = if (z) seq(-2.5, 2.5, length.out = 102) else NA,
      border_color = NA, fontsize = 7, show_colnames = FALSE, filename = f, width = 6.5,
      height = max(3, 0.12 * nrow(m) + 1.5))
  })

  # ---------- Gene explorer ----------
  norm_counts <- reactive({
    cnt_f <- filtered(); cnt <- dat()$counts
    sf <- DESeq2::estimateSizeFactorsForMatrix(cnt_f)
    t(t(cnt) / sf[colnames(cnt)])
  })
  gene_gg <- reactive({
    validate(need(length(input$genes) > 0, "Search for genes in the sidebar."))
    a <- dat()$annot; cd <- dat()$meta
    xv <- if (input$gx_x %in% names(cd)) input$gx_x else "group"
    cv <- if (input$gx_col %in% names(cd)) input$gx_col else "group"
    nc <- norm_counts()
    df <- as.data.frame(nc[input$genes, , drop = FALSE]) |> rownames_to_column("gene_id") |>
      pivot_longer(-gene_id, names_to = "sample", values_to = "value") |>
      left_join(tibble(sample = cd$sample, x = factor(cd[[xv]], levels = unique(cd[[xv]])),
                       col = factor(cd[[cv]], levels = unique(cd[[cv]]))), by = "sample") |>
      left_join(a, by = "gene_id") |>
      mutate(symbol = factor(symbol, levels = unique(a$symbol[match(input$genes, a$gene_id)])))
    pp <- pal_for(df$col)
    dg <- position_dodge(0.7)
    jd <- position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7, seed = 1)
    p <- ggplot(df, aes(x, value, colour = col, shape = col))
    p <- switch(input$gp_type,
      pts = p + geom_point(position = jd, size = 2) +
        stat_summary(aes(group = col), fun.data = mean_se, geom = "errorbar", width = 0.25, position = dg, linewidth = 0.45) +
        stat_summary(aes(group = col), fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.45, position = dg, linewidth = 0.3),
      box = p + geom_boxplot(aes(fill = col), alpha = 0.15, outlier.shape = NA, position = dg, width = 0.6) +
        geom_point(position = jd, size = 1.8),
      vio = p + geom_violin(aes(fill = col), alpha = 0.15, position = dg, scale = "width") +
        geom_point(position = jd, size = 1.8))
    if (isTRUE(input$gp_log)) p <- p + scale_y_log10(labels = scales::label_comma())
    p + facet_wrap(~ symbol, scales = "free_y") +
      scale_colour_manual(values = pp$col) + scale_fill_manual(values = pp$col) + scale_shape_manual(values = pp$shp) +
      labs(x = xv, y = "Normalized count (median-of-ratios)", colour = cv, shape = cv, fill = cv) +
      theme_pub() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  })
  output$p_gene <- renderPlot(gene_gg())
  output$dl_gp_pdf <- downloadHandler("gene_plot.pdf", function(f) { g <- gene_gg(); save_plot(g, f, input$gp_w, input$gp_h) })
  output$dl_gp_png <- downloadHandler("gene_plot.png", function(f) { g <- gene_gg(); save_plot(g, f, input$gp_w, input$gp_h) })
  output$t_gene <- renderDT({
    req(length(input$genes) > 0)
    x <- de()
    x$res |> filter(gene_id %in% input$genes) |> select(symbol, gene_id, baseMean, log2FoldChange, padj, status) |>
      mutate(contrast = paste(x$params$A, "vs", x$params$B)) |>
      datatable(rownames = FALSE, options = list(dom = "t",
                columnDefs = list(list(targets = 4, render = pvalue_render())))) |>   # padj = 5th column
      formatRound(c("baseMean", "log2FoldChange"), 2)
  })

  # ---------- Response map ----------
  rm_df <- reactive({
    req(input$rm_A1, input$rm_B1, input$rm_A2, input$rm_B2)
    validate(need(input$rm_A1 != input$rm_B1 && input$rm_A2 != input$rm_B2, "Each contrast needs two different groups."))
    d <- dds(); an <- dat()$annot
    validate(need(all(c(input$rm_A1, input$rm_B1, input$rm_A2, input$rm_B2) %in% levels(d$group)), "Press 'Run / update model'."))
    get <- function(A, B) {
      r <- results(d, contrast = c("group", A, B), alpha = input$alpha)
      tibble(gene_id = rownames(r), lfc = r$log2FoldChange,
             sig = !is.na(r$padj) & r$padj < input$alpha & abs(r$log2FoldChange) >= input$lfc)
    }
    inner_join(get(input$rm_A1, input$rm_B1), get(input$rm_A2, input$rm_B2), by = "gene_id", suffix = c("1", "2")) |>
      filter(!is.na(lfc1), !is.na(lfc2)) |>
      left_join(an, by = "gene_id") |>
      mutate(class = case_when(sig1 & sig2 ~ "Both", sig1 ~ "Contrast 1 only", sig2 ~ "Contrast 2 only", TRUE ~ "Neither"))
  })
  output$p_rm <- renderPlotly({
    df <- rm_df(); hl <- filter(df, gene_id %in% input$genes)
    l1 <- paste(input$rm_A1, "vs", input$rm_B1); l2 <- paste(input$rm_A2, "vs", input$rm_B2)
    p <- ggplot(df, aes(lfc1, lfc2, colour = class, shape = class,
                        text = paste0(symbol, "<br>", l1, ": ", round(lfc1, 2), "<br>", l2, ": ", round(lfc2, 2)))) +
      geom_hline(yintercept = 0, linewidth = 0.2) + geom_vline(xintercept = 0, linewidth = 0.2) +
      geom_abline(linetype = 2, colour = "grey40") + geom_point(size = 1, alpha = 0.75) +
      geom_point(data = hl, colour = "black", fill = "#F0E442", shape = 21, size = 3) +
      geom_text(data = hl, aes(label = symbol), colour = "black", vjust = -1) +
      scale_colour_manual(values = c(Both = "#CC79A7", `Contrast 1 only` = "#E69F00", `Contrast 2 only` = "#0072B2", Neither = "#C8C8C8")) +
      scale_shape_manual(values = c(Both = 17, `Contrast 1 only` = 16, `Contrast 2 only` = 15, Neither = 16)) +
      labs(x = paste("log2 FC:", l1), y = paste("log2 FC:", l2), colour = NULL, shape = NULL) + theme_pub()
    ggplotly(p, tooltip = "text") |> toWebGL()
  })
  output$rm_interp <- renderUI({
    df <- rm_df()
    r <- cor(df$lfc1, df$lfc2)
    resp <- df[df$sig1 | df$sig2, ]
    sl <- if (nrow(resp) > 10) unname(coef(stats::lm(lfc2 ~ lfc1, resp))[2]) else NA
    n <- table(df$class); getn <- function(k) if (k %in% names(n)) as.integer(n[[k]]) else 0L
    markdown(paste0(
      sprintf("**Agreement.** Fold changes correlate with r = %.2f across %s genes. ", r, format(nrow(df), big.mark = ",")),
      sprintf("%d genes respond in both contrasts, %d only in contrast 1 and %d only in contrast 2. ",
              getn("Both"), getn("Contrast 1 only"), getn("Contrast 2 only")),
      if (!is.na(sl)) sprintf("Among responding genes the slope is %.2f: %s", sl,
        if (sl < 0.8) "responses are weaker in contrast 2." else if (sl > 1.25) "responses are stronger in contrast 2." else "responses have similar magnitude.") else "",
      "\n\nGenes that respond in one contrast only are candidates for context-dependent regulation. ",
      "To test that formally, fit an interaction model (e.g. *~ genotype × treatment*) rather than comparing lists."))
  })
}

shinyApp(ui, server)
