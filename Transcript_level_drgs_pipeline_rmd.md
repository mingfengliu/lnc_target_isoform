---
title: "Temporal Transcriptomics of lncRNA trans-regulation in hiPSCs — Analysis Pipeline"
author: "Dr. Izabela Mamede, iza.mamede@gmail.com"
date: "07-19-2026"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)
```

This document reproduces the computational analyses from raw Salmon quantifications
through the published figures. Steps run in order; intermediate objects are passed
between chunks. Adjust `edger_dir` and `salmon_dir` to your paths.

```{r global-config}
suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(tximport)
  library(ComplexHeatmap)
  library(circlize)
  library(ggrepel)
  library(cowplot)
  library(patchwork)
  library(ggnewscale)
})

salmon_dir <- "salmon_quant"
edger_dir  <- "edgeR_by_lncrna"
fig_dir    <- "figs"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Source-of-truth: transgene -> ENST
tx_map <- tribble(
  ~transgene_name,   ~ensg,            ~enst,
  "BANCR",           "ENSG00000278910","ENST00000624238.3",
  "CRNDE",           "ENSG00000245694","ENST00000502066.7",
  "DANCR",           "ENSG00000226950","ENST00000441504.2",
  "FENDRR",          "ENSG00000268388","ENST00000595886.1",
  "HAGLR",           "ENSG00000224189","ENST00000416928.8",
  "HEIH",            "ENSG00000278970","ENST00000623091.2",
  "HULC",            "ENSG00000285219","ENST00000503668.2",
  "LINC00667",       "ENSG00000263753","ENST00000668807.1",
  "LINC00847",       "ENSG00000245060","ENST00000653584.1",
  "LINC01547",       "ENSG00000183250","ENST00000667826.1",
  "LINC-PINT",       "ENSG00000231721","ENST00000451786.6",
  "LINC-ROR",        "ENSG00000258609","ENST00000553704.3",
  "LNCPRESS1",       "ENSG00000232301","ENST00000429254.3",
  "PNKY",            "ENSG00000283010","ENST00000635423.1",
  "RP11-1055B8.4",   "ENSG00000262877","ENST00000571724.3",
  "TUG1",            "ENSG00000253352","ENST00000643071.1",
  "TUG1",            "ENSG00000253352","ENST00000644773.3"
)

# manuscript cell-line display order (used across figures)
manuscript_order <- c(
  "WT","GFPc4_B1","GFPc5_B1","GFPc5_B2","GFPc5_B3",
  "BANCR","CRNDE","DANCR203","FENDRR","HAGLR","HEIH","HULC",
  "LINC00667c1_B1","LINC00667c2_B1","LINC00667c4_B1","LINC00667c4_B2","LINC00667c5_B1",
  "LINC00847","LINC01547","LINCPINT","LINCROR","LNCPRESS1","PNKY","RP11",
  "TUG1210c3","TUG1210c9","TUG1217"
)

# map a file prefix (clone-level) to its canonical transgene key
extract_transgene_key <- function(prefix) {
  p <- toupper(prefix)
  if (str_detect(p, "^TUG1"))       return("TUG1")
  if (str_detect(p, "^RP11"))       return("RP11-1055B8.4")
  if (str_detect(p, "^LINC\\d+"))   return(str_extract(p, "^LINC\\d+"))
  if (str_detect(p, "^LINC-"))      return(str_extract(p, "^LINC-[A-Z]+"))
  if (str_detect(p, "^LINCPINT"))   return("LINC-PINT")
  if (str_detect(p, "^LINCROR"))    return("LINC-ROR")
  if (str_detect(p, "^LNCPRESS"))   return("LNCPRESS1")
  str_remove(p, "\\d+[A-Z]?\\d*$")
}
```

# 01 — Differential expression across the time course (~time)

Per lncRNA line, import Salmon counts and test for time-dependent expression with
the edgeR quasi-likelihood F-test on the `time` coefficient. One CSV per line.

```{r de-time, eval=FALSE}
parse_sample_id <- function(sample_id) {
  stem  <- sub("_quant$", "", sample_id)
  parts <- strsplit(stem, "_", fixed = TRUE)[[1]]
  rep    <- as.integer(tail(parts, 1))
  time   <- as.double(tail(parts, 2)[1])
  lncrna <- paste(head(parts, -2), collapse = "_")
  list(lncrna_name = lncrna, timepoint = time, replicate = rep)
}

fit_edger_time <- function(counts_mat, meta_df) {
  y      <- DGEList(counts = counts_mat)
  design <- model.matrix(~ time, data = meta_df)
  keep   <- filterByExpr(y, design = design)
  y      <- y[keep, , keep.lib.sizes = FALSE]
  y      <- calcNormFactors(y)
  y      <- estimateDisp(y, design)
  fit    <- glmQLFit(y, design)
  qlf    <- glmQLFTest(fit, coef = "time")
  tab    <- topTags(qlf, n = Inf)$table
  tab$transcript_id <- rownames(tab)
  tab
}

quant_sf   <- file.path(list.dirs(salmon_dir, recursive = FALSE), "quant.sf")
quant_sf   <- quant_sf[file.exists(quant_sf)]
sample_ids <- basename(dirname(quant_sf))

meta_all <- map_dfr(seq_along(sample_ids), function(i) {
  p <- parse_sample_id(sample_ids[i])
  tibble(sample_id = sample_ids[i], files = quant_sf[i],
         lncrna_name = p$lncrna_name, time = p$timepoint, replicate = p$replicate)
})

for (target in sort(unique(meta_all$lncrna_name))) {
  meta <- meta_all %>% filter(lncrna_name == target) %>% arrange(time, replicate)
  files_named <- setNames(meta$files, meta$sample_id)
  txi    <- tximport(files_named, type = "salmon", txOut = TRUE, countsFromAbundance = "no")
  counts <- round(txi$counts)
  meta2  <- meta[match(colnames(counts), meta$sample_id), ]
  tab    <- fit_edger_time(counts, meta2)
  write.csv(tab, file.path(edger_dir, paste0(target, "__time_vs_time0.csv")), row.names = FALSE)
}
```

# 02 — Transcript-level DE with inferential replicates (edgeR v4)

Transcript-level modelling using Salmon Gibbs inferential replicates. The keep-set
is fixed once from observed counts, then reused across all bootstraps.

```{r de-infreps, eval=FALSE}
nboot_use <- 100

rank_reduce_design <- function(X) {
  qr_x <- qr(X); X[, qr_x$pivot[seq_len(qr_x$rank)], drop = FALSE]
}

fit_edger_nofilter <- function(counts_mat, meta_df) {
  y   <- DGEList(counts = counts_mat)
  X   <- rank_reduce_design(model.matrix(~ time + time:treatment, data = meta_df))
  y   <- calcNormFactors(y)
  y   <- estimateDisp(y, X)
  fit <- glmQLFit(y, X)
  tt  <- glmQLFTest(fit, coef = ncol(X))
  tab <- topTags(tt, n = Inf)$table
  tab$transcript_id <- rownames(tab)
  tab
}

for (target in sort(unique(meta_all$lncrna_name))) {
  meta <- meta_all %>% filter(lncrna_name == target) %>% arrange(time, replicate) %>%
    mutate(treatment = factor(ifelse(time == 0, "baseline", "treated"),
                              levels = c("baseline", "treated")))
  files_named <- setNames(meta$files, meta$sample_id)
  txi <- tximport(files_named, type = "salmon", txOut = TRUE,
                  countsFromAbundance = "no", dropInfReps = FALSE)
  counts_obs <- round(txi$counts)
  meta2 <- meta[match(colnames(counts_obs), meta$sample_id), ]

  X    <- rank_reduce_design(model.matrix(~ time + time:treatment, data = meta2))
  keep <- filterByExpr(DGEList(counts = counts_obs), design = X)
  keep_ids <- rownames(counts_obs)[keep]
  if (length(keep_ids) == 0) next
  counts_keep <- counts_obs[keep_ids, , drop = FALSE]

  obs_tab <- fit_edger_nofilter(counts_keep, meta2)
  write.csv(obs_tab, file.path(edger_dir, paste0(target, "__observed_edgeR.csv")), row.names = FALSE)

  inf <- txi$infReps
  if (is.null(inf) || length(inf) == 0) next
  B <- min(nboot_use, ncol(inf[[1]]))

  logFC_mat <- matrix(NA_real_, length(keep_ids), B, dimnames = list(keep_ids, NULL))
  p_mat     <- matrix(NA_real_, length(keep_ids), B, dimnames = list(keep_ids, NULL))
  for (b in seq_len(B)) {
    cb_all <- sapply(colnames(counts_obs), function(s) round(inf[[s]][, b]))
    rownames(cb_all) <- rownames(inf[[1]])
    tab_b  <- fit_edger_nofilter(cb_all[keep_ids, , drop = FALSE], meta2)
    common <- intersect(tab_b$transcript_id, keep_ids)
    logFC_mat[common, b] <- tab_b$logFC[match(common, tab_b$transcript_id)]
    p_mat[common, b]     <- tab_b$PValue[match(common, tab_b$transcript_id)]
  }

  p_med <- apply(p_mat, 1, median, na.rm = TRUE)
  boot_summary <- tibble(
    transcript_id = keep_ids,
    logFC_med = apply(logFC_mat, 1, median, na.rm = TRUE),
    p_med = p_med,
    fdr_med = p.adjust(p_med, method = "BH"),
    frac_fdr_lt_0.05 = apply(p_mat, 1, function(p) mean(p.adjust(p, "BH") < 0.05, na.rm = TRUE))
  )
  write.csv(boot_summary, file.path(edger_dir, paste0(target, "__infreps_summary_B", B, ".csv")), row.names = FALSE)
}
```

# 03 — Import results and annotate

Import the time + timeXtreatment tables, join a transcript dictionary, and keep
the fields used downstream.

```{r import-results}
library(isoformic)

files_timextreat <- list.files(edger_dir, pattern = "__time_plus_timeXtreatment\\.csv$", full.names = TRUE)

timextreat_tbl <- files_timextreat %>%
  set_names() %>%
  map_dfr(~ read_csv(.x, show_col_types = FALSE) %>%
            mutate(lncRNA = sub("__.*$", "", basename(.x)),
                   source_file = basename(.x)),
          .id = "file_path")

fasta_path <- download_reference(version = "43", file_type = "fasta")
tx_to_gene <- make_tx_to_gene(file_path = fasta_path, file_type = "fasta")
tx_to_gene_sel <- tx_to_gene %>% dplyr::select(transcript_id, transcript_name, transcript_type, gene_name)

timextreat_tbl_sel_dic <- timextreat_tbl %>%
  dplyr::select(transcript_id, gene_id, logFC, logCPM, PValue, FDR, lncRNA) %>%
  left_join(tx_to_gene_sel, by = "transcript_id")

saveRDS(timextreat_tbl_sel_dic, "time_timextreat_all.rds")
```

# 04 — DEG consistency (Figure 1D / 2B)

Runs the `~time` LRT for any missing lines (WT, TUG1), collapses to gene level,
and quantifies cross-line sharing to build the consistency bar chart.

```{r deg-consistency, eval=FALSE}
library(tximeta); library(fishpond)

PVAL_CUTOFF <- 0.05; MIN_LOGFC <- 1
TRANSCRIPT_TYPE_FILTER <- "protein_coding"
edger_dir_time <- "edgeR_by_lncrna_time"
dir.create(edger_dir_time, showWarnings = FALSE, recursive = TRUE)

sanitize <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
safe_write_csv <- function(df, path) {
  df <- as.data.frame(df, check.names = FALSE)
  for (nm in names(df)) if (is.list(df[[nm]]))
    df[[nm]] <- vapply(df[[nm]], function(z) paste(z, collapse = ";"), character(1))
  write.csv(df, path, row.names = FALSE, quote = TRUE)
}

coldata_tx <- read.csv("samples.csv", stringsAsFactors = FALSE)
existing   <- str_remove(list.files(edger_dir_time, pattern = "__time_model_LRT\\.csv$"),
                         "__time_model_LRT\\.csv$")
missing    <- setdiff(unique(coldata_tx$lncrna_name), existing)

for (target in missing) {
  cd <- subset(coldata_tx, lncrna_name == target)
  if (nrow(cd) == 0 || !all(file.exists(cd$files)) || length(unique(cd$timepoint)) < 2) next
  se <- tximeta(cd); y <- makeDGEList(se)
  y  <- y[rowSums(y$counts) >= 100, , keep.lib.sizes = FALSE]
  meta <- as.data.frame(colData(se)); meta$time <- as.numeric(meta$timepoint)
  meta <- meta[colnames(y), , drop = FALSE]
  design <- model.matrix(~ time, data = meta)
  y   <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  res <- topTags(glmQLFTest(fit, coef = 2), n = Inf)$table
  res$transcript_id <- rownames(res)
  safe_write_csv(res, file.path(edger_dir_time, paste0(sanitize(target), "__time_model_LRT.csv")))
}

tx_to_gene <- as_tibble(readRDS("tx_to_gene.rds")) %>%
  filter(transcript_type == TRANSCRIPT_TYPE_FILTER)
tx2gene <- tx_to_gene %>% dplyr::select(transcript_id, gene_name) %>% distinct()

all_dets <- list.files(edger_dir_time, pattern = "__time_model_LRT\\.csv$", full.names = TRUE) %>%
  map_dfr(function(f) {
    read.csv(f, stringsAsFactors = FALSE) %>%
      mutate(cell_line = str_remove(basename(f), "__time_model_LRT\\.csv$"),
             significant = PValue < PVAL_CUTOFF & abs(logFC) > MIN_LOGFC)
  }) %>%
  left_join(tx2gene, by = "transcript_id") %>%
  filter(!is.na(gene_name))

cell_lines   <- unique(all_dets$cell_line)
n_cell_lines <- length(cell_lines)

sig_genes_per_line <- all_dets %>%
  filter(significant) %>%
  group_by(gene_name, cell_line) %>%
  slice_max(abs(logFC), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(gene_name, cell_line, logFC)

gene_sharing <- sig_genes_per_line %>%
  count(gene_name, name = "n_lines_sig") %>%
  mutate(share_fraction = n_lines_sig / n_cell_lines,
         share_category = case_when(
           share_fraction >= 0.40 ~ "High Shared (>=75%)",
           share_fraction >= 0.15 ~ "Medium Shared (50-75%)",
           share_fraction >= 0.02 ~ "Low Shared (25-50%)",
           TRUE ~ "Unique (<25%)"))

fig2b_data <- sig_genes_per_line %>%
  left_join(gene_sharing %>% dplyr::select(gene_name, share_category), by = "gene_name") %>%
  count(cell_line, share_category, name = "n_genes") %>%
  group_by(cell_line) %>% mutate(percentage = 100 * n_genes / sum(n_genes)) %>% ungroup() %>%
  mutate(share_category = factor(share_category,
           levels = c("High Shared (>=75%)","Medium Shared (50-75%)","Low Shared (25-50%)","Unique (<25%)")),
         cell_line = factor(cell_line, levels = manuscript_order[manuscript_order %in% cell_lines]))

fig2b <- ggplot(fig2b_data, aes(cell_line, percentage, fill = share_category)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = c("High Shared (>=75%)"="#1a1a1a","Medium Shared (50-75%)"="#666666",
                               "Low Shared (25-50%)"="#b3b3b3","Unique (<25%)"="#4a90d9"), name = "Consistency") +
  labs(x = NULL, y = "Percentage of DEs (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7), legend.position = "bottom")

ggsave2(file.path(fig_dir, "Fig2B_DEG_consistency.pdf"), fig2b, width = 7, height = 3.5)
```

# 05 — DRGS heatmap (Figure 2A)

Loads all per-timepoint contrasts, filters to the 209 DRGS transcripts, and draws
the synchronized log2FC heatmap across all lines and timepoints.

```{r drgs-heatmap, eval=FALSE}
KEEP_TIMEPOINTS <- c(0, 2, 4, 8, 16, 24, 48, 96)

drgs_genes <- unique(read_csv("dox_209.csv", show_col_types = FALSE)$gene_name)
drgs_transcripts <- as_tibble(readRDS("tx_to_gene.rds")) %>%
  filter(gene_name %in% drgs_genes, transcript_type == "protein_coding") %>%
  dplyr::select(transcript_id, gene_name) %>% distinct()

all_results <- list.files(edger_dir, pattern = "__ALL_timepoints_vs_time0\\.csv$", full.names = TRUE) %>%
  map_dfr(function(f) read.csv(f, stringsAsFactors = FALSE) %>%
            mutate(cell_line = str_remove(basename(f), "__ALL_timepoints_vs_time0\\.csv$"))) %>%
  mutate(timepoint_vs_0 = as.numeric(timepoint_vs_0)) %>%
  filter(timepoint_vs_0 %in% setdiff(KEEP_TIMEPOINTS, 0)) %>%
  inner_join(drgs_transcripts, by = "transcript_id")

heatmap_data <- all_results %>%
  mutate(sample_col = paste0(cell_line, "_", timepoint_vs_0)) %>%
  dplyr::select(transcript_id, sample_col, logFC) %>%
  pivot_wider(names_from = sample_col, values_from = logFC)

mat <- as.matrix(heatmap_data %>% dplyr::select(-transcript_id))
rownames(mat) <- heatmap_data$transcript_id

col_info <- tibble(col_name = colnames(mat)) %>%
  mutate(cell_line = str_remove(col_name, "_[0-9.]+$"),
         timepoint = as.numeric(str_extract(col_name, "[0-9.]+$")),
         cell_line_order = match(cell_line, manuscript_order)) %>%
  arrange(cell_line_order, timepoint)

mat <- mat[, col_info$col_name, drop = FALSE]
mat[is.na(mat)] <- 0
mat <- mat[hclust(dist(mat))$order, ]

ht <- Heatmap(mat, name = "log2FC",
  col = colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#b2182b")),
  cluster_rows = FALSE, cluster_columns = FALSE,
  column_split = factor(col_info$cell_line, levels = unique(col_info$cell_line)),
  column_gap = unit(0.5, "mm"), show_row_names = FALSE, show_column_names = FALSE,
  column_title_rot = 90, column_title_gp = gpar(fontsize = 7),
  row_title = paste0("DRGS (n=", nrow(mat), ")"),
  width = unit(16, "cm"), use_raster = TRUE, raster_quality = 2)

pdf(file.path(fig_dir, "Fig2A_DRGS_heatmap.pdf"), width = 16, height = 6)
draw(ht, heatmap_legend_side = "right"); dev.off()
```

# 06 — Transgene–transcript correlations (Figure 3, filter step 1)

For each lncRNA counts file, correlate the transgene ENST against all transcripts
(>100 counts) using Pearson and Spearman; keep pairs significant by either.

```{r correlations, eval=FALSE}
read_counts_matrix <- function(csv_file) {
  df <- read_csv(csv_file, show_col_types = FALSE)
  id_col <- if ("transcript_id" %in% names(df)) "transcript_id" else names(df)[1]
  ids <- as.character(df[[id_col]])
  num <- df %>% dplyr::select(where(is.numeric))
  m <- as.matrix(num); rownames(m) <- ids; m
}

correlate_vs_target <- function(expr_mat, target_tx) {
  if (!target_tx %in% rownames(expr_mat)) return(NULL)
  M <- log2(expr_mat + 1)
  x <- as.numeric(M[target_tx, ])
  if (sd(x, na.rm = TRUE) == 0) return(NULL)
  A <- t(M); n <- nrow(A)
  pr <- as.numeric(cor(A, x, use = "pairwise.complete.obs", method = "pearson"))
  sr <- as.numeric(cor(A, x, use = "pairwise.complete.obs", method = "spearman"))
  tstat <- function(r) { r <- pmin(pmax(r, -0.999999999), 0.999999999); r*sqrt((n-2)/(1-r^2)) }
  tibble(transcript_id = colnames(A),
         pearson_r = pr,  pearson_pvalue  = 2*pt(-abs(tstat(pr)), df = n-2),
         spearman_rho = sr, spearman_pvalue = 2*pt(-abs(tstat(sr)), df = n-2),
         n_samples = n) %>%
    filter(transcript_id != target_tx) %>%
    mutate(pearson_padj_BH  = p.adjust(pearson_pvalue,  "BH"),
           spearman_padj_BH = p.adjust(spearman_pvalue, "BH"))
}

count_files <- list.files(edger_dir, pattern = "__counts_FILTERED_minTotal100\\.csv$", full.names = TRUE)
transgene_to_ensts <- tx_map %>% mutate(key = toupper(transgene_name)) %>%
  group_by(key) %>% summarise(ensts = list(unique(enst)), .groups = "drop")

all_results <- map_dfr(count_files, function(f) {
  fn <- basename(f)
  if (str_detect(toupper(fn), "^GFP")) return(tibble())
  key <- extract_transgene_key(str_split_fixed(fn, "__counts", 2)[, 1])
  ensts <- transgene_to_ensts %>% filter(key == !!key) %>% pull(ensts) %>% unlist()
  if (length(ensts) == 0) return(tibble())
  expr_mat <- read_counts_matrix(f)
  present  <- intersect(rownames(expr_mat), ensts)
  if (length(present) == 0) return(tibble())
  map_dfr(present, function(enst) {
    ct <- correlate_vs_target(expr_mat, enst)
    if (is.null(ct)) return(tibble())
    ct %>% mutate(target_transgene_enst = enst,
                  target_transgene_name = tx_map$transgene_name[match(enst, tx_map$enst)],
                  counts_file = fn)
  }) %>% filter(pearson_padj_BH < 0.05 | spearman_padj_BH < 0.05)
})

dir.create("rds", showWarnings = FALSE)
saveRDS(all_results, "rds/all_corr_res.rds")
```

# 07 — GFP subtraction and dictionary (Figure 3, filter steps 2–3)

Remove correlated pairs also significant in GFP controls, then keep high-confidence
pairs (BH < 0.01 and |Spearman| > 0.7).

```{r gfp-filter, eval=FALSE}
all_results <- readRDS("rds/all_corr_res.rds")

standardize_de <- function(df) {
  nms <- names(df)
  id  <- intersect(c("transcript_id","tx","feature","Geneid"), nms)[1] %||% nms[1]
  lfc <- intersect(c("log2FoldChange","logFC","log2FC"), nms)[1]
  pv  <- intersect(c("pvalue","PValue","p.value"), nms)[1]
  transmute(df, transcript_id = as.character(.data[[id]]),
            log2FC = as.numeric(.data[[lfc]]), pvalue = as.numeric(.data[[pv]]))
}

gfp_files <- list.files(edger_dir, pattern = "GFP.*_time\\d+\\w*_vs_time0\\.csv$", full.names = TRUE)
gfp_exclude_ids <- map_dfr(gfp_files, ~ read_csv(.x, show_col_types = FALSE) %>%
    standardize_de() %>%
    filter(!is.na(transcript_id), pvalue < 0.05, abs(log2FC) > 1, log2FC > 0)) %>%
  distinct(transcript_id) %>% pull(transcript_id)

all_results_filtered <- all_results %>% filter(!transcript_id %in% gfp_exclude_ids)

tx_dic <- tx_to_gene_sel %>% dplyr::select(transcript_id, transcript_name)
all_results_filtered_dic <- all_results_filtered %>%
  left_join(tx_dic %>% rename(transcript_name_trgt_RNA = transcript_name), by = "transcript_id") %>%
  left_join(tx_dic %>% rename(transcript_name_transgene = transcript_name),
            by = c("target_transgene_enst" = "transcript_id")) %>%
  filter(pearson_padj_BH < 0.01, abs(spearman_rho) > 0.7)

saveRDS(all_results_filtered_dic, "rds/all_results_filtered_dic.rds")
```



Filtered pairs plotted as curve-deviation vs GFP per transgene library; selected
targets labelled.

```{r fig3a-violin, eval=FALSE}
all_slope_results_plot  <- readRDS("rds/all_pairs_with_slopes.rds")
all_slope_filtered_plot <- all_slope_results_plot %>% filter(curve_deviation_vs_gfp > 5)

add_display <- function(df) df %>% mutate(
  clone_info = str_extract(counts_file, "[A-Z0-9-]+c[0-9]+(_B[0-9]+)?"),
  target_transgene_name_display = coalesce(clone_info, as.character(target_transgene_name)))

all_slope_results_plot  <- add_display(all_slope_results_plot)
all_slope_filtered_plot <- add_display(all_slope_filtered_plot)

order_disp <- intersect(
  c("BANCR","CRNDE","DANCR","FENDRR","HAGLR","HEIH","HULC",
    "LINC00667c1_B1","LINC00667c2_B1","LINC00667c4_B1","LINC00667c4_B2","LINC00667c5_B1",
    "LINC00847","LINC01547","LINC-PINT","LINC-ROR","LNCPRESS1","PNKY","RP11-1055B8.4",
    "TUG1210c3","TUG1210c9","TUG1217"),
  unique(all_slope_results_plot$target_transgene_name_display))

all_slope_results_plot  <- all_slope_results_plot  %>% mutate(target_transgene_name_display = factor(target_transgene_name_display, levels = order_disp))
all_slope_filtered_plot <- all_slope_filtered_plot %>% mutate(target_transgene_name_display = factor(target_transgene_name_display, levels = order_disp))

TO_MARK <- c("HHLA-202","MACROH2A1-205","MACROH2A1-218","PCDHB5-201","EIF4A3-207",
             "MT-ND1-201","MT-ND2-201","ESRG-201","NEFM-201","HNRNPD-202","PTMS-201",
             "PPL-201","SART1-201","CLIP3-201","GCN-201","ZC3H18-201","MYH9-201","SAFB2-201")
to_label <- all_slope_filtered_plot %>% filter(transcript_name_trgt_RNA %in% TO_MARK)

p_slope <- ggplot() +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 5, ymax = Inf, alpha = 0.3, fill = "#C6BBDD") +
  geom_violin(data = all_slope_results_plot,
              aes(target_transgene_name_display, curve_deviation_vs_gfp),
              fill = "gray95", color = "gray70", alpha = 0.5, scale = "width", width = 0.8) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_hline(yintercept = 0, color = "gray20", linewidth = 0.3) +
  geom_point(data = all_slope_filtered_plot,
             aes(target_transgene_name_display, curve_deviation_vs_gfp, size = neg_log_pval_slope),
             color = "#C6BBDD", alpha = 0.7, shape = 16,
             position = position_jitter(width = 0.3, height = 0, seed = 42)) +
  geom_text_repel(data = to_label,
             aes(target_transgene_name_display, curve_deviation_vs_gfp, label = transcript_name_trgt_RNA),
             size = 2.5, fontface = "bold", color = "black", segment.color = "gray40", max.overlaps = Inf) +
  scale_size_continuous(name = "-log10(p-value)", range = c(2, 8)) +
  labs(x = "Target Transgene", y = "Curve Deviation vs. GFP") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5))

ggsave2(file.path(fig_dir, "Fig3A_curve_deviation.pdf"), p_slope, height = 5, width = 8)
```

# 09 — Genome-wide target distribution (Figure 3C)

Max log2FC per transcript across the genome, targets coloured by inducing lncRNA,
H11 locus shaded. Requires `tx_coord` (transcript_id, transcript_name, chr, start, end).

```{r fig3c-manhattan, eval=FALSE}
H11_START <- 31399560; H11_END <- 31528740
chr_order <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

target_lists <- list(
  CRNDE     = c("CRNDE-220","HHLA1-202","MT-ND6-201","ENST00000674403"),
  DANCR     = c("DDX5-208","MACROH2A1-205","MACROH2A1-218","XACT-203","LRRK1-205"),
  FENDRR    = c("PCDHB5-201"),
  LINC00667 = c("EIF4A3-207"),
  LINC00847 = c("FXYD6-204","ANP32B-201","PPP1R1A-201","PRKACA-201","EMILIN1-201","AEBP1-201",
                "MBD6-201","KCTD17-206","PREX1-201","PRR14-202","HNRNPD-202","PTMS-201",
                "CABLES2-201","EHD2-201","NACC1-202","CCDC9-207","PPL-201","HNRNPUL2-BSCL2-201",
                "SLC25A23-202","PABPC1-207","PPP1R9B-202","TAF15-210","CLIP3-201","SART1-201",
                "NUMA1-233","GCN-201","ZC3H18-201","ARID3A-201","ZC3H4-201","INCENP-201",
                "RNASEK-210","GNAS-211","PRKD2-201","PRKCSH-209","TAF15-214","HNRNPU-228",
                "SCAF1-201","MYH9-201","HNRNPD-203","SUPT6H-202","RPL36A-206","CEP131-203",
                "SAFB2-201","RPSA-204","SLC2A1-215","GSK3A-201","SAFB-211","HNRNPH3-202",
                "TRAF7-208","RUVBL2-212","HNRNPD-214","TTYH3-201","RPL3-202","PGK1-203"),
  LNCPRESS1 = c("MT-ND2-201","MT-ND1-201","MTND1P23-201","MT-CYB-201"),
  PNKY      = c("LINC01013-262","ENST00000691011","ENST00000691049","ESRG-201","ENST00000688847",
                "HHLA1-202","ENST00000588548","ENST00000618966","LRRK1-205","SAFB2-201",
                "LINC00665-228","OLMALINC-224","CCDC88C-203","NEFM-201"))

target_lookup <- imap_dfr(target_lists, ~ tibble(transcript_name = .x, lncrna_group = .y)) %>% distinct()

all_raw <- list.files(edger_dir, pattern = "__ALL_timepoints_vs_time0\\.csv$", full.names = TRUE) %>%
  map_dfr(~ read.csv(.x, stringsAsFactors = FALSE) %>%
            mutate(cell_line = str_remove(basename(.x), "__ALL_timepoints_vs_time0\\.csv$")))

coord_lookup <- tx_coord %>% dplyr::select(transcript_id, transcript_name, chr, start, end) %>%
  distinct(transcript_id, .keep_all = TRUE)

plot_data <- all_raw %>%
  left_join(coord_lookup, by = "transcript_id") %>%
  filter(!is.na(chr), !is.na(transcript_name)) %>%
  mutate(pos = (start + end) / 2) %>%
  group_by(transcript_id, transcript_name, chr, pos) %>%
  slice_max(abs(logFC), n = 1, with_ties = FALSE) %>% ungroup() %>%
  mutate(chr_std = ifelse(str_detect(chr, "^chr"), chr, paste0("chr", chr))) %>%
  filter(chr_std %in% chr_order) %>%
  mutate(chr_std = factor(chr_std, levels = chr_order))

chr_lengths <- plot_data %>% group_by(chr_std) %>%
  summarise(chr_max = max(pos, na.rm = TRUE), .groups = "drop") %>% arrange(chr_std) %>%
  mutate(cumulative_offset = lag(cumsum(chr_max), default = 0),
         chr_centre = cumulative_offset + chr_max / 2)

plot_data <- plot_data %>%
  left_join(chr_lengths %>% dplyr::select(chr_std, cumulative_offset), by = "chr_std") %>%
  mutate(pos_cumulative = pos + cumulative_offset) %>%
  left_join(target_lookup, by = "transcript_name") %>%
  mutate(is_selected = !is.na(lncrna_group),
         lncrna_group = ifelse(is_selected, lncrna_group, "background"))

chr22_offset <- chr_lengths$cumulative_offset[chr_lengths$chr_std == "chr22"]
ymax <- max(abs(plot_data$logFC), na.rm = TRUE)

lncrna_colours <- c(CRNDE="#E41A1C", DANCR="#FF7F00", FENDRR="#4DAF4A", LINC00667="#984EA3",
                    LINC00847="#377EB8", LNCPRESS1="#A65628", PNKY="#F781BF")

p <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.25) +
  geom_point(data = filter(plot_data, !is_selected),
             aes(pos_cumulative, logFC, colour = chr_std), size = 0.4, alpha = 0.35) +
  scale_colour_manual(values = setNames(rep(c("#444444","#AAAAAA"), length.out = length(chr_order)), chr_order), guide = "none") +
  new_scale_colour() +
  geom_point(data = filter(plot_data, is_selected),
             aes(pos_cumulative, logFC, colour = lncrna_group), size = 3) +
  scale_colour_manual(values = lncrna_colours, name = "lncRNA") +
  annotate("rect", xmin = chr22_offset + H11_START, xmax = chr22_offset + H11_END,
           ymin = -Inf, ymax = Inf, fill = "#FFD700", alpha = 0.5) +
  geom_vline(xintercept = chr22_offset + (H11_START + H11_END)/2,
             colour = "firebrick", linewidth = 0.5, linetype = "dashed") +
  scale_x_continuous(breaks = chr_lengths$chr_centre,
                     labels = str_remove(chr_order, "chr"), expand = c(0.01, 0)) +
  labs(x = "Chromosome", y = "Max log2FC vs 0 h") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(size = 8, colour = "grey30"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3))

ggsave2(file.path(fig_dir, "Fig3C_manhattan.pdf"), p, width = 15, height = 5)
```

# 10 — Cis-window check at H11 (Figure 3C inset / S8)

Per cell line, plot log2FC of transcripts in a ±1 Mb (and ±2 Mb) window around the
H11 insertion; confirms absence of local cis effects.

```{r cis-window, eval=FALSE}
H11_CHR <- "chr22"; H11_MID <- 31464150
all_sel <- unique(unlist(target_lists))

all_raw <- all_raw %>% mutate(transcript_id_clean = sub("\\..*", "", transcript_id))

make_cis_plots <- function(window_bp) {
  tx_window <- tx_coord %>%
    mutate(transcript_id_clean = sub("\\..*", "", transcript_id),
           chr_std = ifelse(str_detect(chr, "^chr"), chr, paste0("chr", chr)),
           pos = (start + end) / 2) %>%
    filter(chr_std == H11_CHR, pos >= H11_MID - window_bp, pos <= H11_MID + window_bp) %>%
    dplyr::select(transcript_id_clean, transcript_name, gene_name, pos) %>%
    distinct(transcript_id_clean, .keep_all = TRUE)

  chr22_expr <- all_raw %>%
    inner_join(tx_window, by = "transcript_id_clean") %>%
    group_by(gene_name, transcript_name, pos, cell_line) %>%
    slice_max(abs(logFC), n = 1, with_ties = FALSE) %>% ungroup() %>%
    mutate(dist_kb = (pos - H11_MID) / 1e3, is_selected = transcript_name %in% all_sel)

  plots <- map(sort(unique(chr22_expr$cell_line)), function(cl) {
    d <- chr22_expr %>% filter(cell_line == cl); if (nrow(d) == 0) return(NULL)
    ymax <- ceiling(max(abs(d$logFC), na.rm = TRUE) * 1.1); if (!is.finite(ymax) || ymax == 0) ymax <- 2
    ggplot(d, aes(dist_kb, logFC)) +
      geom_hline(yintercept = 0, linewidth = 0.3) +
      geom_vline(xintercept = 0, colour = "#4D7FBF", linewidth = 0.8) +
      geom_point(data = filter(d, !is_selected), colour = "black", size = 1.8, alpha = 0.7) +
      geom_point(data = filter(d,  is_selected), colour = "#CC0000", size = 2.5) +
      { if (any(d$is_selected)) geom_text_repel(data = filter(d, is_selected),
          aes(label = transcript_name), colour = "#CC0000", fontface = "italic",
          size = 2.8, max.overlaps = 30) } +
      scale_y_continuous(limits = c(-ymax, ymax)) +
      labs(title = cl, x = "distance from H11 (kb)", y = "log2FC") +
      theme_classic(base_size = 10) +
      theme(plot.title = element_text(face = "italic", size = 10, hjust = 0.5),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4))
  }) %>% compact()

  plot_grid(plotlist = plots, ncol = 5)
}

ggsave2(file.path(fig_dir, "Fig3C_cis_H11_1Mb.pdf"), make_cis_plots(1e6), width = 20, height = 3.5 * ceiling(length(unique(all_raw$cell_line))/5))
ggsave2(file.path(fig_dir, "Fig3C_cis_H11_2Mb.pdf"), make_cis_plots(2e6), width = 20, height = 3.5 * ceiling(length(unique(all_raw$cell_line))/5))
```

# 11 — LINC00847 targets: heatmap + trajectories (Figure 4)

Per-transgene log2FC heatmap of high-confidence targets, and count-trajectory plots
(transgene vs GFP vs WT) for representative targets.

```{r fig4-linc00847, eval=FALSE}
all_slope_results_filt <- readRDS("rds/transgene_specific_pairs_filtered.rds")

timepoint_order <- c(0, 2, 4, 8, 16, 24, 48, 96)

read_counts_long <- function(file_path, transcripts_needed, label) {
  read_csv(file_path, show_col_types = FALSE) %>%
    filter(transcript_id %in% transcripts_needed) %>%
    pivot_longer(-transcript_id, names_to = "sample", values_to = "count") %>%
    mutate(timepoint = as.numeric(str_extract(sample, "(?<=_)\\d+\\.?\\d*(?=_)")), source = label) %>%
    filter(!is.na(timepoint))
}

gfp_files <- list.files(edger_dir, pattern = "^GFP.*counts_FILTERED", full.names = TRUE)
wt_files  <- list.files(edger_dir, pattern = "^WT.*counts_FILTERED",  full.names = TRUE)

counts_index <- tibble(path = list.files(edger_dir, pattern = "counts_FILTERED", full.names = TRUE)) %>%
  mutate(file = basename(path),
         prefix = str_split_fixed(file, "__counts", 2)[, 1],
         transgene_key = sapply(prefix, extract_transgene_key)) %>%
  filter(!str_detect(toupper(file), "^GFP|^WT"))

# --- Figure 4A: LINC00847 target heatmap (z-scored trajectories) ---
transgene <- "LINC00847"
pairs_tg  <- all_slope_results_filt %>% filter(target_transgene_name == transgene)
counts_file <- counts_index %>% filter(toupper(transgene_key) == toupper(transgene)) %>% pull(path) %>% .[1]

log2fc_file <- list.files(edger_dir, pattern = paste0("^", transgene, ".*__ALL_timepoints_vs_time0\\.csv$"), full.names = TRUE)[1]
log2fc <- read_csv(log2fc_file, show_col_types = FALSE) %>%
  filter(transcript_id %in% c(unique(pairs_tg$target_transgene_enst), unique(pairs_tg$transcript_id))) %>%
  left_join(pairs_tg %>% dplyr::select(transcript_id, transcript_name_trgt_RNA) %>% distinct(), by = "transcript_id") %>%
  mutate(display_name = coalesce(transcript_name_trgt_RNA, transcript_id))

mat47 <- log2fc %>%
  group_by(display_name, timepoint_vs_0) %>% summarise(logFC = mean(logFC, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = timepoint_vs_0, values_from = logFC) %>%
  column_to_rownames("display_name") %>% as.matrix()
mat47[is.na(mat47)] <- 0
mat47 <- mat47[, order(as.numeric(colnames(mat47)))]

pdf(file.path(fig_dir, "Fig4A_LINC00847_heatmap.pdf"), width = 8, height = max(6, 2 + nrow(mat47) * 0.2))
draw(Heatmap(t(scale(t(mat47))), name = "z-score",
             col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
             cluster_columns = FALSE, column_title = paste0(transgene, " targets")))
dev.off()

# --- Figure 4F / 3D-E: per-target count trajectories vs GFP + WT ---
plot_target <- function(target_id, target_name, counts_file, transgene) {
  plot_data <- bind_rows(
    read_counts_long(counts_file, target_id, paste0(transgene, " Induced")),
    map_dfr(gfp_files, ~ read_counts_long(.x, target_id, "GFP Control")),
    map_dfr(wt_files,  ~ read_counts_long(.x, target_id, "WT Control"))
  ) %>%
    filter(timepoint %in% timepoint_order) %>%
    mutate(timepoint_factor = factor(timepoint, levels = timepoint_order),
           source = factor(source, levels = c(paste0(transgene, " Induced"), "GFP Control", "WT Control"))) %>%
    filter(!is.na(timepoint_factor))
  if (nrow(plot_data) == 0) return(NULL)

  ggplot(plot_data, aes(timepoint_factor, count, color = source, fill = source, group = source)) +
    geom_point(size = 2, alpha = 0.6, position = position_jitter(width = 0.1, height = 0)) +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.2, span = 1) +
    scale_color_manual(values = c(setNames("#D62828", paste0(transgene, " Induced")),
                                  "GFP Control" = "gray50", "WT Control" = "#457B9D")) +
    scale_fill_manual(values = c(setNames("#D62828", paste0(transgene, " Induced")),
                                 "GFP Control" = "gray50", "WT Control" = "#457B9D")) +
    labs(x = "Time (h)", y = "Counts", title = target_name, color = NULL, fill = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
}

target_map <- tx_to_gene_sel %>% dplyr::select(transcript_id, transcript_name)
targets_47 <- pairs_tg %>% left_join(target_map, by = "transcript_id")

plots_47 <- pmap(list(targets_47$transcript_id, targets_47$transcript_name), 
                 ~ plot_target(..1, ..2, counts_file, transgene)) %>% compact()

ggsave2(file.path(fig_dir, "Fig4F_LINC00847_trajectories.pdf"),
        wrap_plots(plots_47, ncol = 3, guides = "collect") & theme(legend.position = "bottom"),
        width = 12, height = 3 * ceiling(length(plots_47) / 3), limitsize = FALSE)
```

```{r session-info}
sessionInfo()
```
