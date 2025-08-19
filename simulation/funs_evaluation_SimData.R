library(clue)

#' @param data data.frame with annotations

evaluateRareCloneAccuracy <- function(data, truth_col, pred_col, 
  rare_labels = NULL,rare_prop=0.1, dominance_threshold = 0.8) {

  # if (!requireNamespace("caret", quietly = TRUE)) install.packages("caret")
  # if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  #if (!requireNamespace("mclust", quietly = TRUE)) install.packages("mclust")
  library(mclust)
  library(dplyr)
  library(caret)


  if (!(truth_col %in% names(data)) || !(pred_col %in% names(data))) {
    stop("Error: Specified columns not found in the input data.")
  }
  
  truth <- as.character(data[[truth_col]])
  pred <- as.character(data[[pred_col]])

  # table(truth,pred)
  num_truth_clusters <- length(unique(truth))
  num_pred_clusters <- length(unique(pred))
  if (num_truth_clusters != num_pred_clusters) {
    warning(sprintf("Cluster count mismatch: %d (truth) vs %d (pred). Evaluation skipped.",
                    num_truth_clusters, num_pred_clusters))
    #return(NULL)
  }

  # Step 1: Automatically select rare labels (if not specified)
  if (is.null(rare_labels)) {
    freq_table <- table(truth)
    rare_candidates <- names(freq_table[freq_table / sum(freq_table) <= rare_prop])
    if (length(rare_candidates) == 0) {
      warning("No rare labels found under the specified rare_prop.")
      return(NULL)
    }
    if (length(rare_candidates) > 1) {
      min_freq <- min(freq_table[rare_candidates])
      rare_labels <- names(freq_table[rare_candidates][freq_table[rare_candidates] == min_freq])[1]
    } else {
      rare_labels <- rare_candidates
    }
  }
  
  # Step 2: Construct binary_truth from the entire data and evaluate which clusters in the predicted labels represent rare clones
  all_pred <- pred
  all_truth <- truth
  pred_labels <- unique(all_pred)
  pred_labels <- na.omit(pred_labels)
  rare_clone_cells <- which(all_truth %in% rare_labels)


  # Step 3: The proportion of rare clones in each predicted label
  pred_to_rare_ratio <- sapply(pred_labels, function(label) {
    pred_cells <- which(all_pred %in% label)
    overlap <- length(intersect(pred_cells, rare_clone_cells))
    total <- length(pred_cells)
    if (total == 0) return(0)
    return(overlap / total)
  })

  # Step 4: Select predicted labels representing rare clones
  rare_pred_labels <- names(pred_to_rare_ratio[pred_to_rare_ratio >= dominance_threshold])

  # Step 5: Construct a binary classification of predictions and true labels
  binary_truth <- all_truth %in% rare_labels
  binary_pred <- all_pred %in% rare_pred_labels
  # Step 6: Evaluating performance
  cm_bin <- confusionMatrix(
    factor(binary_pred, levels = c(TRUE, FALSE)),
    factor(binary_truth, levels = c(TRUE, FALSE)),
    positive = "TRUE"
  )
  
  ari_score <- adjustedRandIndex(truth, pred)

  results <- list(
    rare_label = rare_labels[1],
    TP = sum(binary_pred & binary_truth),
    FP = sum(binary_pred & !binary_truth),
    FN = sum(!binary_pred & binary_truth),
    TN = sum(!binary_pred & !binary_truth),
    accuracy = cm_bin$overall["Accuracy"],
    precision = cm_bin$byClass["Precision"],
    recall = cm_bin$byClass["Recall"],
    f1 = cm_bin$byClass["F1"],
    balanced_accuracy = cm_bin$byClass["Balanced Accuracy"],
    ARI = ari_score,
    confusion_matrix = cm_bin$table,
    rare_predicted_labels = rare_pred_labels,
    dominance_ratios = pred_to_rare_ratio
  )

  return(results)
}

evaluate_Binary <- function(truth,pred){
   cm_bin <- caret::confusionMatrix(
        factor(truth, levels = c(1, 0)),
        factor(pred, levels = c(1, 0)),
        positive = "1"
      )

      results <- list(
        accuracy = cm_bin$overall["Accuracy"],
        precision = cm_bin$byClass["Precision"],
        recall = cm_bin$byClass["Recall"],
        f1 = cm_bin$byClass["F1"],
        balanced_accuracy = cm_bin$byClass["Balanced Accuracy"],
        confusion_matrix = cm_bin$table
      )
      return(results)
}

##truth and pred label clusters do not match
evaluate_MultiClassification <- function(truth,pred){
  truth <- factor(truth)
  pred <- factor(pred)
  # Ensure levels are aligned
  common_levels <- union(levels(truth), levels(pred))
  truth <- factor(truth, levels = common_levels)
  pred  <- factor(pred,  levels = common_levels)

  # Check how many unique classes are in truth
  unique_truth <- unique(truth)
  unique_pred <- unique(pred)

  if(length(unique_pred) <= length(unique_truth)){
    conf_mat <- table(truth = truth, pred = pred)
    max_dim <- max(nrow(conf_mat), ncol(conf_mat))
    conf_square <- matrix(0, max_dim, max_dim)
    conf_square[1:nrow(conf_mat), 1:ncol(conf_mat)] <- conf_mat
    mapping <- solve_LSAP(conf_square, maximum = TRUE) ## Optimal label matching
    # Create a new predicted label: map the predicted number to the actual label number
    pred_labels <- factor(pred)
    aligned_pred <- mapping[pred_labels] 
    aligned_pred <- factor(aligned_pred, levels = levels(truth))
    pred <- aligned_pred
  }else{

    cluster_counts <- sort(table(pred), decreasing = TRUE)
    new_ids <- setNames(seq_along(cluster_counts), names(cluster_counts))
    pred_labels_reindexed <- new_ids[as.character(pred)]
    pred_labels_reindexed <- factor(pred_labels_reindexed,levels = common_levels)
    pred <- pred_labels_reindexed

  }
  


  results <- list()

  if (length(unique_truth) == 1) {
    # Single-class case
    class_label <- unique_truth[1]
    correct <- sum(pred == truth)
    total <- length(truth)
    TP <- sum(truth == class_label & pred == class_label)
    FP <- sum(truth != class_label & pred == class_label)
    FN <- sum(truth == class_label & pred != class_label)
    TN <- sum(truth != class_label & pred != class_label)
    accuracy <- (TP + TN) / (TP + TN + FP + FN)
    precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA
    recall <- if ((TP + FN) > 0) TP / (TP + FN) else NA
    f1 <- if (!is.na(precision) & !is.na(recall) & (precision + recall) > 0) {
      2 * precision * recall / (precision + recall)
    } else {
      NA
    }
    specificity <- if ((TN + FP) > 0) TN / (TN + FP) else NA
    balanced_accuracy <- if (!is.na(recall) & !is.na(specificity)) {
      (recall + specificity) / 2
    } else {
      recall
    }
    cm_table <- table(truth, pred)
    
    results <- list(
      note            = "Only one class in truth. Binary/multiclass metrics not applicable.",
      class_label     = class_label,
      total           = total,
      correct         = correct,
      accuracy        = accuracy,
      balanced_accuracy = balanced_accuracy,
      precision       = precision,
      recall          = recall,
      f1              = f1,
      confusion       = cm_table

    )
    
  }else{


    cm_bin <- caret::confusionMatrix(pred, truth)
    support <- as.numeric(table(truth))
    if(!is.null(nrow(cm_bin$byClass))){
      class_names <- rownames(cm_bin$byClass)
      metrics <- cm_bin$byClass[, c("Precision", "Recall", "F1","Balanced Accuracy")]
      rownames(metrics) <- class_names
       macro_avg <- colMeans(metrics, na.rm = TRUE)
        # 6. Weighted averages
      weighted_avg <- colSums(metrics * support, na.rm = TRUE) / sum(support)

    }else{
      weighted_avg <- cm_bin$byClass
    }
      
    results <- list(
        accuracy = cm_bin$overall["Accuracy"],
        precision = weighted_avg["Precision"],
        recall    = weighted_avg["Recall"],
        f1        = weighted_avg["F1"],
        balanced_accuracy = weighted_avg["Balanced Accuracy"],
        confusion_matrix = cm_bin$table
      )
  }
  return(results)
}



#' @param data data.frame with annotations
evaluateCNVPerformance <- function(data,
                                    truth_col = "GroundTruth",
                                    pred_col = "prediction",
                                    normal_value = 2,
                                    evaluate_CNVregion = FALSE,
                                    evaluate_CNVstates = TRUE) {

  library(caret)
  library(mclust)   # for adjustedRandIndex
  library(dplyr)

  if (!(truth_col %in% names(data)) || !(pred_col %in% names(data))) {
    stop("Specified column names do not exist in the input data.")
  }
  df <- data %>% dplyr::filter(!is.na(.data[[truth_col]]))
  truth <- df[[truth_col]]
  pred  <- df[[pred_col]]
  # Check how many unique classes are in truth
  unique_truth <- unique(truth)
  unique_pred <- unique(pred)

  results <- list()

  if (length(unique_truth) == 1) {
    # Single-class case
    class_label <- unique_truth[1]
    correct <- sum(pred == truth)
    total <- length(truth)
    TP <- sum(truth == class_label & pred == class_label)
    FP <- sum(truth != class_label & pred == class_label)
    FN <- sum(truth == class_label & pred != class_label)
    TN <- sum(truth != class_label & pred != class_label)
    accuracy <- (TP + TN) / (TP + TN + FP + FN)
    precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA
    recall <- if ((TP + FN) > 0) TP / (TP + FN) else NA
    f1 <- if (!is.na(precision) & !is.na(recall) & (precision + recall) > 0) {
      2 * precision * recall / (precision + recall)
    } else {
      NA
    }
    cm_table <- table(truth, pred)
    
    results <- list(
      note            = "Only one class in truth. Binary/multiclass metrics not applicable.",
      class_label     = class_label,
      total           = total,
      correct         = correct,
      accuracy        = accuracy,
      precision       = precision,
      recall          = recall,
      f1              = f1,
      confusion       = cm_table

    )
    
  }else{

    # ---- 1. Binary classification: Is the region a CNV (vs. normal) ----
    if (evaluate_CNVregion) {
      bin_truth <- ifelse(truth != normal_value, 1, 0)
      bin_pred  <- ifelse(pred != normal_value & !is.na(pred), 1, 0)

      cm_bin <- caret::confusionMatrix(
        factor(bin_pred, levels = c(1, 0)),
        factor(bin_truth, levels = c(1, 0)),
        positive = "1"
      )

      results$CNV_region_detection <- list(
        accuracy = cm_bin$overall["Accuracy"],
        precision = cm_bin$byClass["Precision"],
        recall = cm_bin$byClass["Recall"],
        f1 = cm_bin$byClass["F1"],
        balanced_accuracy = cm_bin$byClass["Balanced Accuracy"],
        confusion_matrix = cm_bin$table
      )
    }
    # ---- 2. Multiclass classification: Are CNV states correctly predicted ----
    if (evaluate_CNVstates) {
        df <- data %>% dplyr::filter(!is.na(.data[[truth_col]]),!is.na(.data[[pred_col]]))
        truth <- df[[truth_col]]
        pred  <- df[[pred_col]]
        #平均绝对误差（MAE, Mean Absolute Error）
        mae <- mean(abs(as.integer(truth) - as.integer(pred)))
        rmse <- sqrt(mean((as.integer(truth) - as.integer(pred))^2))

        truth <- factor(df[[truth_col]])
        pred  <- factor(df[[pred_col]])
          # Ensure levels are aligned
        common_levels <- union(levels(truth), levels(pred))
        truth <- factor(truth, levels = common_levels)
        pred  <- factor(pred,  levels = common_levels)
        cm <- caret::confusionMatrix(pred, truth)
        ari <- mclust::adjustedRandIndex(as.integer(truth), as.integer(pred))
         # 3. Support per class
        support <- as.numeric(table(truth))
        class_names <- rownames(cm$byClass)
        metrics <- cm$byClass[, c("Precision", "Recall", "F1","Balanced Accuracy")]
        rownames(metrics) <- class_names
         macro_avg <- colMeans(metrics, na.rm = TRUE)
          # 6. Weighted averages
        weighted_avg <- colSums(metrics * support, na.rm = TRUE) / sum(support)
        


        results$CNV_states_detection  <- list(
          accuracy = cm$overall["Accuracy"],
          precision = weighted_avg["Precision"],
          recall    = weighted_avg["Recall"],
          f1        = weighted_avg["F1"],
          balanced_accuracy = weighted_avg["Balanced Accuracy"],
          kappa = cm$overall["Kappa"],  # Cohen’s Kappa
          ARI = ari,                  # Adjusted Rand Index
          MAE = mae,
          RMSE = rmse,
          confusion_matrix = cm$table,
          per_class = metrics
        )
    }
  }

 


  return(results)

}


extract_runtime_info <- function(file_path) {
  log_lines <- readLines(file_path)

  start_time <- NA
  end_time <- NA
  total_runtime <- NA
  memory_used <- NA

  for (line in log_lines) {
    if (grepl("^Start time: ", line)) {
      start_time <- gsub("^Start time: (.*)$", "\\1", line)
    }
    if (grepl("^End time: ", line)) {
      end_time <- gsub("^End time: (.*)$", "\\1", line)
    }
    if (grepl("^Total runtime \\(minutes\\): ", line)) {
      total_runtime <- as.numeric(gsub("^Total runtime \\(minutes\\): (.*)$", "\\1", line))
    }
    if (grepl("^Memory used \\(MB\\): ", line)) {
      memory_used <- as.numeric(gsub("^Memory used \\(MB\\): (.*)$", "\\1", line))
    }
  }
  result <- data.frame(
    StartTime = start_time,
    EndTime = end_time,
    Runtime_Minutes = total_runtime,
    Memory_MB = memory_used
  )
  
  return(result)
}









expand_matrix_by_row_sampling <- function(mat, M) {
  N <- nrow(mat)
  
  if (M <= N) {
    sampled_indices <- sample(1:N, M, replace = FALSE)
  } else {
    sampled_indices <- sample(1:N, M, replace = TRUE)
  }
  
  expanded_mat <- mat[sampled_indices, , drop = FALSE]
  
  return(expanded_mat)
}

expand_matrix_by_col_sampling <- function(mat, M) {
  N <- ncol(mat)
  
  if (M <= N) {
    sampled_indices <- sample(1:N, M, replace = FALSE)
  } else {
    sampled_indices <- sample(1:N, M, replace = TRUE)
  }
  
  expanded_mat <- mat[,sampled_indices, drop = FALSE]
  
  return(expanded_mat)
}


# Function to randomly generate CNV segments
generate_cnv_column <- function(bin_table, total_cnv_bins, cnv_values, min_cnv_len = 100, min_gap_bins = 500) {
  total_bins <- nrow(bin_table)
  bin_table$cn <- 2
  
  current_cnv_bins <- 0
  cursor <- 1

  last_end_idx <- -Inf  # The last bin index of the previous CNV

  while (cursor + min_cnv_len - 1 <= total_bins && current_cnv_bins < total_cnv_bins) {

    # Whether the bin interval limit is met
    if (cursor < last_end_idx + min_gap_bins) {
      cursor <- cursor + 1
      next
    }
    max_len <- min(500, total_bins - cursor + 1, total_cnv_bins - current_cnv_bins)
    if (max_len < min_cnv_len) break
    # The length of the current segment
    cnv_len <- max_len

    start_idx <- cursor
    end_idx <- cursor + cnv_len - 1
    # Actual genomic endpoint
    current_end_bp <- bin_table$end[end_idx]


    cn <- sample(cnv_values, 1)
    bin_table$cn[start_idx:end_idx] <- cn

    current_cnv_bins <- current_cnv_bins + cnv_len
    last_end_idx <- end_idx
    cursor <- end_idx + 1  
  }

  

  return(bin_table$cn)
}

simulate_arm_level_cnv <- function(bin_df, cnv_chrs = NULL, cnv_values = c(1, 3, 4)) {
  # bin_df: Data frame of genome bins with columns: seqnames, start, end, bin
  # cnv_chrs: Optional. If NULL, use all chromosomes in bin_df
  # cnv_values: Vector of CNV values to assign to arms. Values may repeat.

  # Get valid chromosomes
  available_chrs <- unique(bin_df$seqnames)
  if (is.null(cnv_chrs)) {
    cnv_chrs <- available_chrs
  } else {
    cnv_chrs <- intersect(cnv_chrs, available_chrs)
  }

  # Create all possible chromosome-arm combinations (chr_p and chr_q)
  chrom_arm_combinations <- unlist(lapply(cnv_chrs, function(chr) {
    c(paste0(chr, "_p"), paste0(chr, "_q"))
  }))

  if (length(chrom_arm_combinations) == 0) {
    stop("No valid chromosome arms available.")
  }

  selected_arms <- chrom_arm_combinations
  
  # Sample CNV values (with replacement if needed)
  assigned_cnv_values <- sample(cnv_values, length(selected_arms), replace = TRUE)

  assigned_cnv <- data.frame(
    chr = sub("_(p|q)$", "", selected_arms),
    arm = sub("^.*_(p|q)$", "\\1", selected_arms),
    cnv_value = assigned_cnv_values
  )

  # Initialize CNV values
  bin_df$cnv_value <- 2  # default diploid

  for (i in seq_len(nrow(assigned_cnv))) {
    chr <- assigned_cnv$chr[i]
    arm <- assigned_cnv$arm[i]
    cnv_val <- assigned_cnv$cnv_value[i]

    chr_bins <- bin_df[bin_df$seqnames == chr, ]
    mid_pos <- median(chr_bins$start)

    if (arm == "p") {
      target_bins <- chr_bins[chr_bins$end <= mid_pos, ]
    } else {
      target_bins <- chr_bins[chr_bins$start > mid_pos, ]
    }

    bin_df$cnv_value[bin_df$bin %in% target_bins$bin] <- cnv_val
    message(sprintf("Assigned CNV value %d to %s %s-arm", cnv_val, chr, arm))
  }

  return(bin_df)
}









