# ==============================================================================
# MAIN ANALYSIS SCRIPT: HPC / POSIT RSTUDIO VERSION (OPTIMIZED)
# ==============================================================================

# 1. Load Libraries & Utils
library(foreach)
library(doParallel)
library(doRNG)
library(DESeq2)
library(dplyr)
library(tidyr)
library(bnlearn) 
library(limma)

# 清理旧的并行环境
try({ stopCluster(cl) }, silent = TRUE) 
registerDoSEQ() 

# Load your utility functions 
source("/scratch/Shares/rinn/ML/RNAseq_FULL_MODEL/combined time/RF_results/dbn_utils.R") 

# ------------------------------------------------------------------------------
# 2. Load Data & Prepare Groups
# ------------------------------------------------------------------------------
setwd("/scratch/Shares/rinn/ML/RNAseq_FULL_MODEL/combined time/RF_results")

# 加载数据
if(file.exists("/scratch/Shares/rinn/ML/RNAseq_FULL_MODEL/combined time/results/combined_counts_with_gene_name.RData")) {
  load("/scratch/Shares/rinn/ML/RNAseq_FULL_MODEL/combined time/results/combined_counts_with_gene_name.RData")
} else { stop("Data file not found: combined_counts_with_gene_name.RData") }

if(file.exists("/scratch/Shares/rinn/ML/RNAseq_FULL_MODEL/combined time/figures/fig4/sig_lncRNA_18_to_24.RData")) {
  load("/scratch/Shares/rinn/ML/RNAseq_FULL_MODEL/combined time/figures/fig4/sig_lncRNA_18_to_24.RData")
} else { stop("Data file not found: sig_lncRNA_18_to_24.RData") }

if(file.exists("/scratch/Shares/rinn/ML/RNAseq_FULL_MODEL/combined time/RF_results/importance_df.RData")) {
  load("/scratch/Shares/rinn/ML/RNAseq_FULL_MODEL/combined time/RF_results/importance_df.RData")
} else { stop("Data file not found: importance_df.RData") }

# 数据预处理
counts_matrix <- combined_counts_plus[, -665]
rownames(counts_matrix) <- counts_matrix$gene_id
counts_matrix <- counts_matrix[, -1]
counts_matrix <- as.matrix(counts_matrix)
lncRNAs <- unique(sig_lncRNA_18_to_24$gene_id) 

dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData   = data.frame(row.names = colnames(counts_matrix)),
  design    = ~ 1
)

vst_mat <- assay(vst(dds, blind=TRUE)) 
sel_samples <- grep("_(0|2|4|8|16|24|48|96)_", colnames(vst_mat), value = TRUE)
expr <- vst_mat[, sel_samples]

gene_var <- apply(expr, 1, var)
top_genes <- union(lncRNAs, names(sort(gene_var, decreasing=TRUE))[1:5000]) 

expr <- expr[top_genes, ]
expr_t <- t(expr)

# 分组逻辑
set.seed(123)
dbn_input_list_all_groups <- list()

cat("Generating groups based on RF importance...\n")
for (lnc in unique(importance_df$lncRNA)) {
  df <- importance_df %>%
    filter(lncRNA == lnc) %>%
    arrange(desc(IncMSE)) %>% 
    head(2000)
  
  all_predictors <- df$feature
  if (length(all_predictors) == 0) next
  
  total_genes <- length(all_predictors)
  
  # === [关键参数] ===
  # 你设置了 200。如果跑得慢，可以改回 100。
  group_size <- 50
  min_coverage <- 20
  
  required_groups <- ceiling((total_genes * min_coverage) / group_size * 1.2)
  actual_groups <- required_groups
  
  coverage_count <- rep(0, total_genes)
  names(coverage_count) <- all_predictors
  groups <- list()
  
  for (i in 1:actual_groups) {
    undercovered_genes <- all_predictors[coverage_count < min_coverage]
    if (length(undercovered_genes) >= group_size) {
      selected <- sample(undercovered_genes, group_size)
    } else {
      n_undercovered <- length(undercovered_genes)
      n_additional <- group_size - n_undercovered
      covered_genes <- all_predictors[coverage_count >= min_coverage]
      max_allowed_coverage <- max(min_coverage + 2, 10) 
      eligible_covered <- covered_genes[coverage_count[covered_genes] < max_allowed_coverage]
      if (length(eligible_covered) >= n_additional) {
        additional <- sample(eligible_covered, n_additional)
      } else { additional <- eligible_covered }
      selected <- unique(c(undercovered_genes, additional))
      if (length(selected) < group_size) {
        needed <- group_size - length(selected)
        supplemental <- sample(all_predictors, needed)
        selected <- unique(c(selected, supplemental))
      }
    }
    coverage_count[selected] <- coverage_count[selected] + 1
    groups[[i]] <- selected
  }
  
  for (i in 1:length(groups)) {
    group_name <- paste0(lnc, "_group", i)
    genes <- c(lnc, groups[[i]])
    genes <- genes[genes %in% colnames(expr_t)]
    if (length(genes) > 1) {
      dbn_input_list_all_groups[[group_name]] <- expr_t[, genes, drop = FALSE]
      attr(dbn_input_list_all_groups[[group_name]], "lnc_name") <- lnc
      attr(dbn_input_list_all_groups[[group_name]], "group_number") <- i
      attr(dbn_input_list_all_groups[[group_name]], "predictors") <- groups[[i]]
    }
  }
}
cat("Groups generated. Total groups:", length(dbn_input_list_all_groups), "\n")


# ------------------------------------------------------------------------------
# 3. SETUP PARALLEL PROCESSING (FORK Mode)
# ------------------------------------------------------------------------------

# 优先读取 Slurm 变量
slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK")

if (slurm_cpus != "") {
  num_cores <- as.integer(slurm_cpus)
  cat("Detected SLURM allocation. Using", num_cores, "cores.\n")
} else {
  # 如果没有 Slurm，手动设置
  # 请确保这个数字不要超过你申请的核数！
  num_cores <- 120
  cat("No SLURM allocation detected. Using manual setting:", num_cores, "cores.\n")
}

# === [关键修复] ===
# 必须使用 type = "FORK" 来避免 'all 128 connections are in use' 报错
# FORK 模式在 Linux HPC 上更高效，且没有连接数限制
cl <- makeCluster(num_cores, type = "FORK", outfile = "")
registerDoParallel(cl)

# ------------------------------------------------------------------------------
# 4. RUN PARALLEL DBN PIPELINE
# ------------------------------------------------------------------------------

log_dir <- "logs_monitor"
if(!dir.exists(log_dir)) dir.create(log_dir)

# 参数设置
# 注意：group_size=200 配合 100次 bootstrap 计算量巨大
# 如果发现速度慢，可以考虑将 N_BOOTSTRAP 降为 50
N_BOOTSTRAP <- 25
TIME_POINTS <- c(0, 2, 4, 8, 16, 24, 48, 96) 
group_names <- names(dbn_input_list_all_groups)
total_groups <- length(group_names)

if(!exists("dbn_input_list_all_groups")) stop("Missing: dbn_input_list_all_groups")
if(!exists("g2s")) stop("Missing: g2s mapping object")

cat("Starting Parallel Analysis on HPC (FORK mode)...\n")
cat("Total Groups:", total_groups, "\n")
cat("Bootstrap Iterations:", N_BOOTSTRAP, "\n")

res_dir <- "results_parts"
if(!dir.exists(res_dir)) dir.create(res_dir)

# 执行并行循环
results_list <- foreach(
  expr_mat = dbn_input_list_all_groups,
  grp_name = group_names,
  .packages = c("bnlearn", "dplyr", "tidyr", "limma", "stringr"),
  # FORK 模式下 .export 其实是可选的，因为会复制父环境，但为了保险保留
  .export = c("run_temporal_DBN_with_bootstrap_v2", 
              "run_temporal_DBN_v2_with_groups",
              "build_robust_dynamic_bayesian_network_final",
              "extract_robust_inter_time_edges_completely_fixed",
              "perform_time_differential_analysis_v2",
              "analyze_lag_correlations_v2",
              "integrate_temporal_results_fixed_v3",
              "match_samples_by_cellline_v2",
              "convert_gene_ids_to_names_local",
              "align_time_series",
              "extract_dbn_evidence_for_lnc",
              "infer_relationship_type",
              "aggregate_bootstrap_results_v2",
              "g2s", "TIME_POINTS", "N_BOOTSTRAP"),
  .options.RNG = 123 
) %dorng% {
  
  # --- [断点续传逻辑] ---
  # 如果这个组之前已经跑完了，直接跳过
  finished_flag <- file.path(log_dir, paste0(grp_name, ".finished"))
  if(file.exists(finished_flag)) {
    # 返回 NULL 或者之前的占位符，最后我们会过滤掉 NULL
    return(NULL)
  }
  
  # --- [LOGGING START] ---
  my_log_file <- file.path(log_dir, paste0(grp_name, ".log"))
  start_time <- Sys.time()
  
  # 状态标记
  running_flag <- file.path(log_dir, paste0(grp_name, ".running"))
  file.create(running_flag)
  
  cat(paste0(">>> [START] ", grp_name, " at ", start_time, "\n"), 
      file = my_log_file, append = FALSE)
  
  # --- [CORE CALCULATION] ---
  lnc_name <- attr(expr_mat, "lnc_name")
  if(length(dim(expr_mat)) > 2) expr_mat <- expr_mat[,,1]
  expr_mat <- as.data.frame(expr_mat)
  
  res <- tryCatch({
    run_temporal_DBN_with_bootstrap_v2(
      expr_mat = expr_mat,
      lnc_name = lnc_name,
      g2s_mapping = g2s, 
      time_points = TIME_POINTS,
      bootstrap_iterations = N_BOOTSTRAP,
      group_name = grp_name
    )
  }, error = function(e) {
    cat(paste0("[ERROR] ", e$message, "\n"), file = my_log_file, append = TRUE)
    return(NULL)
  })
  
  if(!is.null(res)) {
    save(res, file = file.path(res_dir, paste0(grp_name, ".RData")))
  }
  
  # --- [LOGGING END] ---
  end_time <- Sys.time()
  duration <- round(difftime(end_time, start_time, units = "mins"), 2)
  
  cat(paste0("<<< [DONE] ", grp_name, " at ", end_time, 
             " | Duration: ", duration, " mins\n"), 
      file = my_log_file, append = TRUE)
  
  # 切换标记
  file.remove(running_flag)
  file.create(finished_flag)
  
  return(res)
}

# 关闭集群
stopCluster(cl)
cat("All tasks completed.\n")

# 筛选有效结果 (去掉因为断点续传跳过的 NULL)
valid_results <- results_list[sapply(results_list, is.list)]
save(valid_results, file = "DBN_Optimized_Results_HPC.RData")


# ------------------------------------------------------------------------------
# 5. ORGANIZE & AGGREGATE RESULTS (Load from Disk)
# ------------------------------------------------------------------------------
cat("Starting aggregation from disk files...\n")

# 1. 找到所有保存好的小文件
result_files <- list.files("results_parts", pattern = "\\.RData$", full.names = TRUE)
cat("Found", length(result_files), "completed group files.\n")

if(length(result_files) == 0) stop("No result files found in 'results_parts'!")

# 2. 循环读取并合并
# 使用 lapply 快速读取所有文件
all_edges_list <- lapply(result_files, function(f) {
  # 加载 RData，它会把变量 'res' 加载到局部环境
  env <- new.env()
  load(f, envir = env)
  res <- env$res
  
  if(!is.null(res$bootstrap_stats) && nrow(res$bootstrap_stats) > 0){
    stats <- res$bootstrap_stats
    stats$group_id <- res$group_name
    stats$lncRNA <- res$lnc_display_name 
    return(stats)
  } else {
    return(NULL)
  }
})

# 3. 合并大表格
raw_network_table <- do.call(rbind, all_edges_list)

cat("Raw table rows:", nrow(raw_network_table), "\n")

# 4. Ensemble Voting (和之前一样)
cat("Performing Ensemble Voting...\n")

final_network_table <- raw_network_table %>%
  group_by(lncRNA, gene) %>%
  summarise(
    n_groups_tested = n_distinct(group_id),
    n_groups_found = n(),
    ensemble_score = mean(bootstrap_frequency),
    max_group_freq = max(bootstrap_frequency),
    final_relationship = names(sort(table(relationship_type), decreasing=TRUE))[1],
    .groups = 'drop'
  ) %>%
  filter(ensemble_score >= 0.7 | (n_groups_found >= 2 & max_group_freq >= 0.8)) %>%
  arrange(desc(ensemble_score))

# 5. Save Final
save(final_network_table, file = "Final_Robust_LncRNA_Network_Ensemble_HPC.RData")
cat("Done.\n")