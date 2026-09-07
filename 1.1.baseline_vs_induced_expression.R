## ============================================================
## Supplementary Figure: baseline (WT/GFP 0h vs own transgene 0h)
## vs. induced (own transgene, max timepoint) expression
## for the 17 lncRNA transgenes
##
## Run this locally where DESeq2 + tidyverse are installed.
## Needs: col_data.csv, combined_counts_plus.csv (same folder as this script,
## or edit the paths below)
## ============================================================

library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)
library(ggbeeswarm)   # install.packages("ggbeeswarm") if missing, for geom_quasirandom

## ---- 1. Load data -------------------------------------------------------
col_data <- read_csv("col_data.csv", show_col_types = FALSE) %>%
  mutate(time = as.numeric(time)) %>%
  column_to_rownames("sample_id")

counts_raw <- read_csv("combined_counts_plus.csv", show_col_types = FALSE)
gene_name  <- setNames(counts_raw$gene_name, counts_raw$gene_id)
counts_mat <- counts_raw %>%
  select(-gene_name) %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

# make sure column order matches col_data row order
counts_mat <- counts_mat[, rownames(col_data)]

## ---- 2. DESeq2 size factors / normalized counts ------------------------
dds <- DESeqDataSetFromMatrix(countData = round(counts_mat),
                               colData   = col_data,
                               design    = ~1)     # only used for normalization here
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

## ---- 3. Gene + cell-line mapping (edit if your naming differs) --------
gene_map <- c(
  "BANCR"="ENSG00000278910", "CRNDE"="ENSG00000245694", "DANCR"="ENSG00000226950",
  "FENDRR"="ENSG00000268388", "HAGLR"="ENSG00000224189", "HEIH"="ENSG00000278970",
  "HULC"="ENSG00000285219", "LINC00667"="ENSG00000263753", "LINC00847"="ENSG00000245060",
  "LINC01547"="ENSG00000183250", "LINC-PINT"="ENSG00000231721", "LINC-ROR"="ENSG00000258609",
  "LNCPRESS1"="ENSG00000232301", "PNKY"="ENSG00000283010", "RP11-1055B8.4"="ENSG00000262877",
  "TUG1"="ENSG00000253352"
)

line_map <- list(
  BANCR          = list(BANCR = c("BANCR")),
  CRNDE          = list(CRNDE = c("CRNDE")),
  DANCR          = list(DANCR = c("DANCR203")),
  FENDRR         = list(FENDRR = c("FENDRR")),
  HAGLR          = list(HAGLR = c("HAGLR")),
  HEIH           = list(HEIH = c("HEIH")),
  HULC           = list(HULC = c("HULC")),
  LINC00667      = list(LINC00667 = c("LINC00667c1b1","LINC00667c2b1","LINC00667c4b1",
                                       "LINC00667c4b2","LINC00667c5b1")),
  LINC00847      = list(LINC00847 = c("LINC00847")),
  LINC01547      = list(LINC01547 = c("LINC01547")),
  `LINC-PINT`    = list(`LINC-PINT` = c("LINCPINT")),
  `LINC-ROR`     = list(`LINC-ROR` = c("LINCROR")),
  LNCPRESS1      = list(LNCPRESS1 = c("LNCPRESS1")),
  PNKY           = list(PNKY = c("PNKY")),
  `RP11-1055B8.4`= list(`RP11-1055B8.4` = c("RP11")),
  TUG1           = list(`TUG1-210` = c("TUG1210c3","TUG1210c9"),
                         `TUG1-217` = c("TUG1217"))
)

gfp_conds <- c("GFPc4b1","GFPc5b1","GFPc5b2","GFPc5b3","GFPc18")

## Pooled "no transgene induced" reference: WT + GFP at 0h.
## Rationale: WT alone has only n=3 at 0h in this dataset -- too small to be
## a stable baseline on its own. GFP and WT both lack the test lncRNA's
## Dox-inducible cassette, so pooling them (n=20) gives a much more robust
## baseline for "what does this locus look like with no transgene present."
## We deliberately do NOT fold in the other 15 lncRNA lines' own 0h samples,
## because those cell lines DO carry a (different, uninduced) Tet-On cassette
## and keeping them out preserves a clean, minimal-confound reference group.
wtgfp_0h_samples <- rownames(col_data)[col_data$time == 0 &
                                          col_data$condition %in% c("WT", gfp_conds)]

## ---- 4. Build long-format table for plotting ---------------------------
results_list <- list()
for (sym in names(gene_map)) {
  ensg <- gene_map[[sym]]
  gid  <- grep(paste0("^", ensg, "(\\.|$)"), rownames(norm_counts), value = TRUE)[1]
  if (is.na(gid)) { warning("Gene not found: ", sym); next }

  for (line_name in names(line_map[[sym]])) {
    conds <- line_map[[sym]][[line_name]]
    own_samples <- rownames(col_data)[col_data$condition %in% conds]
    own_0h      <- own_samples[col_data[own_samples, "time"] == 0]

    # empirically pick the timepoint (t>0) with the highest mean normalized
    # count for this gene in its own line -> "max induction" timepoint
    own_pos_t <- own_samples[col_data[own_samples, "time"] > 0]
    if (length(own_pos_t) > 0) {
      by_t <- tapply(norm_counts[gid, own_pos_t], col_data[own_pos_t, "time"], mean)
      max_t <- as.numeric(names(by_t)[which.max(by_t)])
    } else {
      max_t <- NA
    }
    own_max <- own_samples[col_data[own_samples, "time"] == max_t]

    df <- bind_rows(
      tibble(sample_id = wtgfp_0h_samples, group = "WT+GFP 0h"),
      tibble(sample_id = own_0h,           group = "own 0h"),
      tibble(sample_id = own_max,          group = paste0("own max (t=", max_t, "h)"))
    ) %>%
      mutate(gene_symbol = sym, lncRNA_line = line_name, gene_id = gid,
             normalized_count = norm_counts[gid, sample_id],
             time_h = col_data[sample_id, "time"],
             condition = col_data[sample_id, "condition"])
    results_list[[paste(sym, line_name)]] <- df
  }
}

plot_df <- bind_rows(results_list) %>%
  mutate(group = factor(group, levels = c("WT+GFP 0h", "own 0h",
                                           unique(group[str_starts(group, "own max")]))))

write_csv(plot_df, "baseline_vs_induced_expression_data_R.csv")

## ---- 5. Plot: all points shown (quasirandom/jitter) + median crossbar --
## one facet per lncRNA line, log10 y-axis, no error bars -- every replicate
## is plotted individually as requested.
plot_df <- plot_df %>%
  mutate(group_simple = case_when(
    group == "WT+GFP 0h" ~ "WT+GFP 0h",
    group == "own 0h" ~ "own 0h",
    TRUE ~ "own max"
  ) %>% factor(levels = c("WT+GFP 0h", "own 0h", "own max")))

p <- ggplot(plot_df, aes(x = group_simple, y = normalized_count + 1, color = group_simple)) +
  geom_quasirandom(width = 0.15, size = 1.6, alpha = 0.85) +
  stat_summary(fun = median, geom = "crossbar", width = 0.4, color = "black", linewidth = 0.3) +
  scale_y_log10() +
  scale_color_manual(values = c("WT+GFP 0h" = "grey50", "own 0h" = "#4C72B0", "own max" = "#C44E52")) +
  facet_wrap(~ paste0(gene_symbol, " (", lncRNA_line, ")"), scales = "free_y", ncol = 4) +
  labs(x = NULL, y = "DESeq2 normalized count + 1", color = NULL,
       title = "Baseline (WT+GFP 0h vs. own 0h) vs. induced (own max) expression") +
  theme_bw(base_size = 10) +
  theme(legend.position = "top", strip.background = element_blank(),
        axis.text.x = element_text(angle = 0))

ggsave("Supplementary_Figure_baseline_vs_induced.pdf", p, width = 14, height = 14)
ggsave("Supplementary_Figure_baseline_vs_induced.png", p, width = 14, height = 14, dpi = 300)

## ---- 6. Quick leakiness check: own 0h vs pooled WT+GFP 0h --------------
leak_check <- plot_df %>%
  filter(group_simple %in% c("WT+GFP 0h", "own 0h")) %>%
  group_by(gene_symbol, lncRNA_line, group_simple) %>%
  summarise(mean_count = mean(normalized_count), n = n(), .groups = "drop") %>%
  pivot_wider(names_from = group_simple, values_from = c(mean_count, n)) %>%
  mutate(fold_own0h_vs_baseline = (`mean_count_own 0h` + 1) / (`mean_count_WT+GFP 0h` + 1)) %>%
  arrange(desc(fold_own0h_vs_baseline))

print(leak_check, n = 20)
write_csv(leak_check, "leakiness_check_own0h_vs_pooled_baseline.csv")
