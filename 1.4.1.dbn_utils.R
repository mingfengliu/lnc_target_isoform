# ==============================================================================
# UTILITY FUNCTIONS FOR DYNAMIC BAYESIAN NETWORK ANALYSIS
# ==============================================================================

library(bnlearn)
library(visNetwork)
library(dplyr)
library(tidyr)
library(limma)
library(pROC)

# ------------------------------------------------------------------------------
# 1. Gene ID Conversion Helper
# ------------------------------------------------------------------------------
convert_gene_ids_to_names_local <- function(gene_ids, g2s_mapping){
  base_ids <- gsub("\\..*","", gene_ids)
  g2s_base_ids <- gsub("\\..*","", g2s_mapping$gene_id)
  mapping <- setNames(g2s_mapping$gene_name, g2s_base_ids)
  
  converted <- sapply(base_ids, function(id){
    if(id %in% names(mapping) && !is.na(mapping[[id]]) && mapping[[id]]!="") {
      return(mapping[[id]]) 
    } else {
      return(id)
    }
  })
  names(converted) <- gene_ids
  # cat(sprintf("Gene ID conversion: %.1f%% successful\n", sum(converted!=base_ids)/length(converted)*100))
  return(converted)
}

# ------------------------------------------------------------------------------
# 2. Sample Matching (CRITICAL UPDATE: Preserving Replicates)
# ------------------------------------------------------------------------------
# Updated to NOT average replicates. It pairs replicates to maximize N.
match_samples_by_cellline_v2 <- function(data_prev, data_curr, t_prev, t_curr){
  samples_prev <- rownames(data_prev)
  samples_curr <- rownames(data_curr)
  
  # Extract cell line and replicate info (assuming format: Cell_Time_Rep)
  # Adjust the split logic if your naming convention differs
  info_prev <- data.frame(
    sample_id = samples_prev,
    cell = sapply(strsplit(samples_prev, "_"), `[`, 1),
    rep = sapply(strsplit(samples_prev, "_"), `[`, 3), # Assuming Rep is 3rd part
    stringsAsFactors = FALSE
  )
  
  info_curr <- data.frame(
    sample_id = samples_curr,
    cell = sapply(strsplit(samples_curr, "_"), `[`, 1),
    rep = sapply(strsplit(samples_curr, "_"), `[`, 3),
    stringsAsFactors = FALSE
  )
  
  # Merge based on Cell line AND Replicate to pair them correctly
  # If exact replicate matching isn't possible, matching by Cell line alone is okay,
  # but strictly pairing Rep1->Rep1 is statistically cleaner for time-series.
  merged_map <- inner_join(info_prev, info_curr, by = c("cell", "rep"), suffix = c("_prev", "_curr"))
  
  if(nrow(merged_map) == 0) return(NULL)
  
  # Construct the combined matrix
  mat_prev <- data_prev[merged_map$sample_id_prev, , drop=FALSE]
  mat_curr <- data_curr[merged_map$sample_id_curr, , drop=FALSE]
  
  colnames(mat_prev) <- paste0(colnames(mat_prev), "_t", t_prev)
  colnames(mat_curr) <- paste0(colnames(mat_curr), "_t", t_curr)
  
  combined <- cbind(mat_prev, mat_curr)
  return(as.data.frame(combined))
}

# ------------------------------------------------------------------------------
# 3. Time-Course Differential Expression Analysis
# ------------------------------------------------------------------------------
perform_time_differential_analysis_v2 <- function(expr_mat, sample_info, lnc_name){
  results <- list()
  time_points <- sort(unique(sample_info$time))
  
  for(i in 2:length(time_points)){
    t_prev <- time_points[i-1]
    t_curr <- time_points[i]
    
    samples_prev <- sample_info$sample[sample_info$time == t_prev]
    samples_curr <- sample_info$sample[sample_info$time == t_curr]
    
    if(length(samples_prev) < 2 | length(samples_curr) < 2) next
    
    tryCatch({
      # Calculate means
      prev_expr_means <- colMeans(expr_mat[samples_prev, , drop = FALSE], na.rm = TRUE)
      curr_expr_means <- colMeans(expr_mat[samples_curr, , drop = FALSE], na.rm = TRUE)
      fold_changes <- curr_expr_means - prev_expr_means
      
      # Limma analysis
      expr_subset <- expr_mat[c(samples_prev, samples_curr), ]
      gene_vars <- apply(expr_subset, 2, var, na.rm = TRUE)
      valid_genes <- which(gene_vars > 0)
      
      if(length(valid_genes) < 5) {
        # Fallback to simple FC if not enough variance
        fc_threshold <- quantile(abs(fold_changes), 0.75, na.rm = TRUE)
        sig_genes <- names(which(abs(fold_changes) > fc_threshold))
        results[[paste(t_prev, "to", t_curr)]] <- list(fold_changes = fold_changes, significant_genes = sig_genes)
        next
      }
      
      expr_filtered <- expr_subset[, valid_genes, drop = FALSE]
      design <- model.matrix(~ factor(c(rep("prev", length(samples_prev)), rep("curr", length(samples_curr)))))
      
      fit <- lmFit(t(expr_filtered), design)
      fit <- eBayes(fit)
      pvals <- fit$p.value[, 2]
      fdr <- p.adjust(pvals, method = "fdr")
      sig_genes <- names(fdr)[fdr < 0.05 & !is.na(fdr)]
      
      results[[paste(t_prev, "to", t_curr)]] <- list(
        fold_changes = fold_changes,
        significant_genes = sig_genes,
        fdr = fdr
      )
      
    }, error = function(e) {
      # cat("    DE failed for", t_prev, "to", t_curr, "\n")
      NULL
    })
  }
  return(results)
}

# ------------------------------------------------------------------------------
# 4. Extract Edges from BN Object (Helper)
# ------------------------------------------------------------------------------
extract_robust_inter_time_edges_completely_fixed <- function(bn, lnc_name, t_prev, t_curr) {
  edges <- data.frame(from = character(), to = character(), from_time = character(), to_time = character(), stringsAsFactors = FALSE)
  nodes <- bnlearn::nodes(bn)
  
  prev_pattern <- paste0("_t", t_prev)
  curr_pattern <- paste0("_t", t_curr)
  
  for (node in nodes) {
    parents <- bn$nodes[[node]]$parents
    if (length(parents) > 0) {
      for (parent in parents) {
        # Check if edge is inter-time (prev -> curr)
        if (grepl(prev_pattern, parent) && grepl(curr_pattern, node)) {
          
          clean_parent <- gsub(prev_pattern, "", parent)
          clean_node <- gsub(curr_pattern, "", node)
          
          new_edge <- data.frame(
            from = clean_parent,
            to = clean_node,
            from_time = t_prev,
            to_time = t_curr,
            stringsAsFactors = FALSE
          )
          edges <- rbind(edges, new_edge)
        }
      }
    }
  }
  return(edges)
}

# ==============================================================================
# 5. Build Dynamic Bayesian Network (OPTIMIZED with Blacklist)
# ==============================================================================
build_robust_dynamic_bayesian_network_final <- function(time_slice_data, lnc_name, time_points){
  results <- list()
  
  for(i in 2:length(time_points)){
    t_prev <- as.character(time_points[i-1])
    t_curr <- as.character(time_points[i])
    
    if(!t_prev %in% names(time_slice_data) | !t_curr %in% names(time_slice_data)) next
    
    data_prev <- time_slice_data[[t_prev]]
    data_curr <- time_slice_data[[t_curr]]
    
    # 1. Identify common genes
    common_genes <- intersect(colnames(data_prev), colnames(data_curr))
    if(!lnc_name %in% common_genes) next
    
    # 2. Match samples by cell line (preserving replicates)
    combined_data <- match_samples_by_cellline_v2(
      data_prev[, common_genes, drop = FALSE], 
      data_curr[, common_genes, drop = FALSE], 
      t_prev, t_curr
    )
    
    if(is.null(combined_data) || nrow(combined_data) < 10) next 
    
    # Critical Optimization: Create Blacklist to constrain search space
    # Logic: Only allow edges from Prev_Time -> Curr_Time.
    # Ban: Curr -> Prev (Time travel), Prev -> Prev (Intra-slice), Curr -> Curr (Intra-slice)
    
    all_nodes <- colnames(combined_data)
    prev_cols <- grep(paste0("_t", t_prev, "$"), all_nodes, value = TRUE)
    curr_cols <- grep(paste0("_t", t_curr, "$"), all_nodes, value = TRUE)
    
    # Ban all edges originating from Current (No time travel, no intra-slice loops)
    from_curr <- expand.grid(from = curr_cols, to = all_nodes, stringsAsFactors = FALSE)
    
    # Ban all edges pointing to Previous (Previous can only be parents, not children)
    to_prev <- expand.grid(from = all_nodes, to = prev_cols, stringsAsFactors = FALSE)
    
    bl <- unique(rbind(from_curr, to_prev))
    
    tryCatch({
      # 4. Structure Learning (with Blacklist)
      # maxp=5 restricts parent set size to prevent overfitting and improve speed
      bn <- hc(combined_data, score = "bge", blacklist = bl, maxp = 3)
      
      # 5. Extract robust inter-time edges
      edges <- extract_robust_inter_time_edges_completely_fixed(bn, lnc_name, t_prev, t_curr)
      
      results[[paste(t_prev, "->", t_curr)]] <- list(
        network = bn, 
        inter_time_edges = edges, 
        n_edges = nrow(edges)
      )
      
    }, error = function(e) {
      # Errors are caught in the main loop logging mechanism
    })
  }
  return(results)
}

# ------------------------------------------------------------------------------
# 6. Lag Correlation Analysis
# ------------------------------------------------------------------------------
align_time_series <- function(lnc_expr, gene_expr, sample_info, time_points, lag) {
  aligned <- data.frame(lnc = numeric(), gene = numeric(), time_pair = character(), stringsAsFactors = FALSE)
  
  for (i in (1 + lag):length(time_points)) {
    lnc_time <- time_points[i - lag]
    gene_time <- time_points[i]
    
    lnc_samples <- sample_info$sample[sample_info$time == lnc_time]
    gene_samples <- sample_info$sample[sample_info$time == gene_time]
    
    # Simple matching for correlation (can be improved, but sufficient for lag check)
    lnc_cells <- sapply(strsplit(lnc_samples, "_"), `[`, 1)
    gene_cells <- sapply(strsplit(gene_samples, "_"), `[`, 1)
    common_cells <- intersect(lnc_cells, gene_cells)
    
    for (cell in common_cells) {
      # Here we average just for correlation calculation purpose
      lnc_val <- mean(lnc_expr[lnc_samples[lnc_cells == cell], , drop = TRUE], na.rm = TRUE)
      gene_val <- mean(gene_expr[gene_samples[gene_cells == cell], , drop = TRUE], na.rm = TRUE)
      if (!is.na(lnc_val) && !is.na(gene_val)) {
        aligned <- rbind(aligned, data.frame(lnc = lnc_val, gene = gene_val, time_pair = paste(lnc_time, "->", gene_time)))
      }
    }
  }
  return(aligned)
}

analyze_lag_correlations_v2 <- function(expr_mat, sample_info, lnc_name, time_points){
  results <- list()
  if(!lnc_name %in% colnames(expr_mat)) return(results)
  lnc_expr <- expr_mat[, lnc_name, drop = FALSE]
  if (var(lnc_expr, na.rm = TRUE) == 0) return(results)
  
  for(gene in setdiff(colnames(expr_mat), lnc_name)){
    gene_expr <- expr_mat[, gene, drop = FALSE]
    if (var(gene_expr, na.rm = TRUE) == 0) next
    
    best_corr <- 0
    best_lag <- 0
    best_p <- 1
    
    for(lag in 1:(length(time_points)-1)){
      aligned <- align_time_series(lnc_expr, gene_expr, sample_info, time_points, lag)
      if(nrow(aligned) > 5){
        ct <- tryCatch(cor.test(aligned$lnc, aligned$gene), error = function(e) NULL)
        if(!is.null(ct) && !is.na(ct$p.value)){
          if(ct$p.value < 0.05 && abs(ct$estimate) > abs(best_corr)){
            best_corr <- ct$estimate; best_lag <- lag; best_p <- ct$p.value
          }
        }
      }
    }
    if (best_lag > 0) results[[gene]] <- list(best_lag = best_lag, max_correlation = best_corr, pvalue = best_p)
  }
  return(results)
}

# ------------------------------------------------------------------------------
# 7. Integration & Scoring
# ------------------------------------------------------------------------------
extract_dbn_evidence_for_lnc <- function(dbn_results, lnc_name) {
  temporal_evidence <- list()
  if (is.null(dbn_results)) return(temporal_evidence)
  
  for (trans_name in names(dbn_results)) {
    edges <- dbn_results[[trans_name]]$inter_time_edges
    if (!is.null(edges) && nrow(edges) > 0) {
      # Lnc Regulates Gene
      lnc_regs <- edges[edges$from == lnc_name & edges$to != lnc_name, ]
      for (i in seq_len(nrow(lnc_regs))) {
        g <- lnc_regs$to[i]
        if(is.na(g)) next
        if(is.null(temporal_evidence[[g]])) temporal_evidence[[g]] <- list()
        temporal_evidence[[g]]$lnc_regulates <- c(temporal_evidence[[g]]$lnc_regulates, trans_name)
      }
      # Gene Regulates Lnc
      reg_lnc <- edges[edges$to == lnc_name & edges$from != lnc_name, ]
      for (i in seq_len(nrow(reg_lnc))) {
        g <- reg_lnc$from[i]
        if(is.na(g)) next
        if(is.null(temporal_evidence[[g]])) temporal_evidence[[g]] <- list()
        temporal_evidence[[g]]$regulates_lnc <- c(temporal_evidence[[g]]$regulates_lnc, trans_name)
      }
    }
  }
  return(temporal_evidence)
}

infer_relationship_type <- function(evidence, n_lnc_regulates, n_regulates_lnc) {
  if (n_lnc_regulates > 0 && n_regulates_lnc > 0) {
    if (n_lnc_regulates >= n_regulates_lnc) return("lnc_regulates_gene") else return("gene_regulates_lnc")
  } else if (n_lnc_regulates > 0) return("lnc_regulates_gene")
  else if (n_regulates_lnc > 0) return("gene_regulates_lnc")
  else if (!is.null(evidence$lag)) {
    if (evidence$lag$max_correlation > 0.5) return("associated_pos")
    if (evidence$lag$max_correlation < -0.5) return("associated_neg")
  }
  return("associated")
}

integrate_temporal_results_fixed_v3 <- function(diff_results, dbn_results, lag_results, lnc_name, time_points) {
  temporal_evidence <- extract_dbn_evidence_for_lnc(dbn_results, lnc_name)
  
  # Merge Lag info
  for(g in names(lag_results)){
    if(is.null(temporal_evidence[[g]])) temporal_evidence[[g]] <- list()
    temporal_evidence[[g]]$lag <- lag_results[[g]]
  }
  
  # Score
  scores_df <- data.frame()
  
  for (g in names(temporal_evidence)) {
    ev <- temporal_evidence[[g]]
    score <- 0
    n_lnc_reg <- length(ev$lnc_regulates)
    n_reg_lnc <- length(ev$regulates_lnc)
    
    # Weighted Scoring System
    if(n_lnc_reg > 0) score <- score + (n_lnc_reg * 8)
    if(n_reg_lnc > 0) score <- score + (n_reg_lnc * 8)
    
    lag_corr <- 0; lag_p <- 1
    if(!is.null(ev$lag)) {
      lag_corr <- ev$lag$max_correlation
      lag_p <- ev$lag$pvalue
      if(abs(lag_corr) > 0.6 && lag_p < 0.05) score <- score + 4
    }
    
    # Diversity Bonus
    if((n_lnc_reg + n_reg_lnc > 0) && (!is.null(ev$lag))) score <- score + 3
    
    rel_type <- infer_relationship_type(ev, n_lnc_reg, n_reg_lnc)
    
    scores_df <- rbind(scores_df, data.frame(
      gene = g,
      temporal_score = score,
      relationship_type = rel_type,
      n_lnc_regulates = n_lnc_reg,
      n_regulates_lnc = n_reg_lnc,
      lag_correlation = lag_corr,
      stringsAsFactors = FALSE
    ))
  }
  
  if(nrow(scores_df) > 0) scores_df <- scores_df[order(-scores_df$temporal_score), ]
  
  return(list(temporal_evidence = temporal_evidence, integrated_scores = scores_df))
}

# ------------------------------------------------------------------------------
# 8. Main Wrapper Function (Per Group)
# ------------------------------------------------------------------------------
run_temporal_DBN_v2_with_groups <- function(expr_mat, lnc_name, g2s_mapping, time_points=c(0,24,48,96), group_name = NULL, group_number = NULL) {
  
  # Clean IDs
  all_gene_ids <- colnames(expr_mat)
  gene_name_mapping <- convert_gene_ids_to_names_local(all_gene_ids, g2s_mapping)
  colnames(expr_mat) <- gene_name_mapping
  lnc_display_name <- convert_gene_ids_to_names_local(lnc_name, g2s_mapping)[1]
  
  # Sample Info
  sample_info <- data.frame(sample=rownames(expr_mat)) %>%
    separate(sample, into=c("cell","time","rep"), sep="_", remove=FALSE)
  sample_info$time <- as.numeric(sample_info$time)
  detected_time_points <- sort(unique(sample_info$time))
  
  # Time slices
  time_slice_data <- list()
  for(t in detected_time_points){
    tsamps <- sample_info$sample[sample_info$time==t]
    if(length(tsamps)>1) time_slice_data[[as.character(t)]] <- expr_mat[tsamps,,drop=FALSE]
  }
  
  # Run Modules
  diff_res <- tryCatch(perform_time_differential_analysis_v2(expr_mat, sample_info, lnc_display_name), error=function(e) list())
  dbn_res <- tryCatch(build_robust_dynamic_bayesian_network_final(time_slice_data, lnc_display_name, detected_time_points), error=function(e) list())
  lag_res <- tryCatch(analyze_lag_correlations_v2(expr_mat, sample_info, lnc_display_name, detected_time_points), error=function(e) list())
  
  # Integrate
  final_res <- tryCatch(integrate_temporal_results_fixed_v3(diff_res, dbn_res, lag_res, lnc_display_name, detected_time_points), error=function(e) NULL)
  
  if(is.null(final_res)) return(NULL)
  
  return(list(
    lnc_display_name = lnc_display_name,
    original_lnc_id = lnc_name,
    integrated_results = final_res,
    expr_data = expr_mat,
    group_name = group_name,
    group_number = group_number
  ))
}

# ------------------------------------------------------------------------------
# 9. Bootstrap Aggregation Helper (Added Missing Function)
# ------------------------------------------------------------------------------
aggregate_bootstrap_results_v2 <- function(bootstrap_results_list, lnc_name) {
  if (length(bootstrap_results_list) == 0) return(NULL)
  
  # Combine all results into one dataframe
  all_edges <- do.call(rbind, lapply(1:length(bootstrap_results_list), function(i) {
    df <- bootstrap_results_list[[i]]
    if(is.null(df) || nrow(df) == 0) return(NULL)
    # Filter for valid edges (Score > 0 implies some evidence found)
    df <- df[df$temporal_score > 0, c("gene", "relationship_type", "temporal_score")]
    if(nrow(df) > 0) df$iter <- i
    return(df)
  }))
  
  if(is.null(all_edges)) return(NULL)
  
  # Calculate statistics
  edge_stats <- all_edges %>%
    group_by(gene) %>%
    summarise(
      bootstrap_frequency = n_distinct(iter) / length(bootstrap_results_list),
      avg_score = mean(temporal_score),
      consensus_relationship = names(which.max(table(relationship_type))),
      consistency = max(table(relationship_type)) / n(), # How stable is the direction?
      .groups = 'drop'
    ) %>%
    arrange(desc(bootstrap_frequency))
  
  return(edge_stats)
}

# ------------------------------------------------------------------------------
# 10. Bootstrap Wrapper (Can be run inside parallel loop)
# ------------------------------------------------------------------------------
run_temporal_DBN_with_bootstrap_v2 <- function(expr_mat, lnc_name, g2s_mapping, 
                                               time_points, bootstrap_iterations = 50,
                                               group_name = NULL) {
  
  # 1. Run on Original Data
  original_result <- run_temporal_DBN_v2_with_groups(expr_mat, lnc_name, g2s_mapping, time_points, group_name)
  
  if(bootstrap_iterations < 1) return(original_result)
  
  # 2. Bootstrap Loop
  sample_info <- data.frame(sample=rownames(expr_mat)) %>%
    separate(sample, into=c("cell","time","rep"), sep="_", remove=FALSE)
  unique_cells <- unique(sample_info$cell)
  
  boot_res_list <- list()
  
  for(i in 1:bootstrap_iterations){
    # Resample Cell Lines (Preserving time structure)
    boot_cells <- sample(unique_cells, length(unique_cells), replace = TRUE)
    
    # Construct bootstrap matrix
    boot_samples <- c()
    for(cl in boot_cells){
      # We must handle duplicated cell names by making them unique in rownames if needed,
      # but simple row extraction allows duplicates in R matrices.
      # However, for downstream mapping, unique rownames are safer.
      # Simplified strategy: extract indices
      indices <- which(sample_info$cell == cl)
      boot_samples <- c(boot_samples, indices)
    }
    
    boot_expr <- expr_mat[boot_samples, , drop=FALSE]
    # Fix rownames to be unique for limma/bnlearn
    rownames(boot_expr) <- make.unique(rownames(boot_expr))
    
    res <- tryCatch({
      run_temporal_DBN_v2_with_groups(boot_expr, lnc_name, g2s_mapping, time_points, paste0(group_name, "_b", i))
    }, error = function(e) NULL)
    
    if(!is.null(res) && !is.null(res$integrated_results)) {
      boot_res_list[[i]] <- res$integrated_results$integrated_scores
    }
  }
  
  # 3. Aggregate
  original_result$bootstrap_stats <- aggregate_bootstrap_results_v2(boot_res_list, lnc_name)
  
  return(original_result)
}

#nohup Rscript dbn_analysis.R > main_output.log 2>&1 &
#pkill -u mili7948 R
#ls -l logs_monitor/*.finished | wc -l
