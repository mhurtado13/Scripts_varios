
library(caret)
library(dplyr)
library(mlbench)
library(h2o)
library(ggplot2)
library(patchwork)   
library(stringr)
library(tidyr)

compute_features.training.ML = function(features_train, task_type = c("classification", "survival"), target_var = NULL, trait.positive = NULL,
                                        time_var = NULL, event_var = NULL, metric = "Accuracy", stack, k_folds = 10, n_rep = 5, LODO = FALSE,
                                        batch_id = NULL, file_name = NULL, return = FALSE,
                                        fold_construction_fun = NULL, fold_construction_args_fixed = NULL, fold_construction_args_tunable = NULL){
  
  # ---------------------------------------------------------------------------
  # Validate task_type and required arguments
  # ---------------------------------------------------------------------------
  
  if (task_type == "classification") {
    if (is.null(target_var) || is.null(trait.positive)) {
      stop("For classification, both `target_var` and `trait.positive` must be provided.")
    }
  } else if (task_type == "survival") {
    if (is.null(time_var) || is.null(event_var)) {
      stop("For survival, both `time_var` and `event_var` must be provided.")
    }
  } else {
    stop("Invalid task_type. Must be either 'classification' or 'survival'.")
  }
  
  # ---------------------------------------------------------------------------
  # === CASE 1: CLASSIFICATION TASK ==========================================
  # ---------------------------------------------------------------------------
  
  if (task_type == "classification") {
    
    #Set training set for classification
    train_data = features_train %>%
      data.frame() %>%
      dplyr::mutate(Trait = target_var,
                    target = as.factor(ifelse(Trait == trait.positive, 'yes', 'no'))) %>%
      dplyr::select(-Trait)
    
    train_data$target <- factor(train_data$target, levels = c("no", "yes"))  # Order, just in case to ensure positive class is not well defined
    
    if(LODO == T){
      train_data = train_data %>%
        dplyr::mutate(dataset = traitData_train[,batch_id])
    }
    
    #Cross-validation training
    training = compute_k_fold_CV(train_data, k_folds = k_folds, n_rep = n_rep, metric = metric, stacking = stack,
                                 file_name = file_name, LODO = LODO, return= return,
                                 fold_construction_fun = fold_construction_fun, fold_construction_args_fixed = fold_construction_args_fixed,
                                 fold_construction_args_tunable = fold_construction_args_tunable)
    
  }
  
  # ---------------------------------------------------------------------------
  # === CASE 2: SURVIVAL ANALYSIS TASK =======================================
  # ---------------------------------------------------------------------------
  else if (task_type == "survival"){
    
    #Prepare training data
    train_data <- features_train %>%
      data.frame() %>%
      dplyr::mutate(
        time  = time_var,
        event = as.numeric(event_var == trait.positive)
      )
    
    if (LODO == TRUE) {
      train_data = train_data %>%
        dplyr::mutate(dataset = traitData_train[,batch_id])
    }
    
    #Split into features + outcomes
    df_outcome  <- train_data %>% dplyr::select(time, event)
    df_features <- train_data %>% dplyr::select(-time, -event)
    
    #Run survival cross-validation and tuning
    training = compute_k_fold_CV_survival(
      df_features = df_features,
      df_outcome  = df_outcome,
      outcome_col = "time",
      event_col   = "event",
      k_folds = k_folds,
      n_rep = n_rep,
      ncores = ncores,
      file_name   = file_name,
      fold_construction_fun = fold_construction_fun,
      fold_construction_args_fixed = fold_construction_args_fixed,
      fold_construction_args_tunable = fold_construction_args_tunable
    )
    
  }
  
  ####################################################Predicting
  if(length(training)!=0){
    return(training)
  }else{  #No features are selected as predictive
    message("No features selected as predictive after Boruta runs. No model returned.")
    return(NULL)
  }
  
}

compute_k_fold_CV = function(model, k_folds, n_rep, stacking = FALSE, metric = "Accuracy", file_name = NULL, LODO = FALSE,
                             return = FALSE, fold_construction_fun = NULL,
                             fold_construction_args_fixed = NULL,
                             fold_construction_args_tunable = NULL){
  
  if(!(metric %in% c("AUROC", "AUPRC","Accuracy"))){
    stop("The metric assigned is not supported. Choose either accuracy or AUC.")
  }
  
  if(is.null(fold_construction_fun)){ ### Preprocessing (remove collinear variables and low variance)
    if(LODO == TRUE){
      train_data = preprocess_features(model %>% dplyr::select(-dataset), cor_thresh = 0.9, target_col = "target") %>%
        dplyr::mutate(dataset = model$dataset)
    }else{
      train_data = preprocess_features(model, cor_thresh = 0.9, target_col = "target")
    }
  }else{
    train_data = model
  }
  
  rm(model) #Clean memory
  gc()
  
  ######### Machine Learning models
  
  ######### Stratify K fold cross-validation
  if(LODO == T){
    multifolds = construct_stratified_cohort_folds(train_data, 'dataset', 'target', k_folds = k_folds, n_rep = n_rep)
    train_data = train_data %>% dplyr::select(-dataset) #After creating multifolds we remove this variable to be able to train
  }else{
    multifolds = caret::createMultiFolds(train_data[,'target'], k = k_folds, times = n_rep) #repeated folds
  }
  
  if(!is.null(fold_construction_fun)){ #Custom function provided, using normal CV

    # Custom fold construction (is running in parallel)
    do.call(fold_construction_fun, c(list(data = train_data, folds = multifolds), fold_construction_args_fixed, fold_construction_args_tunable))
    
    ### Extract the file names of the folds
    result_files <- list.files("Results", pattern = "^fold_.*\\.rds$", full.names = TRUE)
    fold_data = vector("list", length(result_files))
    
    # Initialize master list to store everything in memory
    models_all_folds <- vector("list", length(result_files))
    
    for (fold_i in seq_along(result_files)) {
      
      result = readRDS(result_files[[fold_i]])
      
      models_all_params <- vector("list", length(result))
      
      for (parameter_i in seq_along(result)) {
        
        train_df <- result[[parameter_i]][["train_data"]]
        test_df  <- result[[parameter_i]][["test_data"]]
        
        test_df = cbind(test_df, target = result[[parameter_i]][["obs_test"]])
        
        res <- run_h2o_fold(train_df, test_df, "target", 
                            nfolds_inner = 5, max_runtime_secs = 120)
        
        # Store results
        models_all_params[[parameter_i]] <- res
        
      }
      
      models_all_folds[[fold_i]] <- models_all_params

    }
    
    nested_result <- select_best_model_nested(models_all_folds)
    
    grDevices::pdf(paste0("Results/AUROC_CV_", file_name, ".pdf"), width = 10)
    plot_metric_summary_nested(nested_result, metric = "auc", statistic = "median")
    dev.off()
    
    grDevices::pdf(paste0("Results/AUPRC_CV_", file_name, ".pdf"), width = 10)
    plot_metric_summary_nested(nested_result, metric = "auprc", statistic = "median")
    dev.off()
    
    #fam_summary <- aggregate_nested_families(models_all_folds)
    #plot_metric_summary(fam_summary, metric = "auc", statistic = "median")
    #best_info <- pick_best_model_outer(fam_summary)
    custom_output <- do.call(fold_construction_fun,
                             c(list(data = train_data, bestune = result[[nested_result$best_parameter]][["params"]]), fold_construction_args_fixed))
    
    output = list(nested_result, custom_output)
  }
  
  return(output)
  
}


run_h2o_fold <- function(train_df, test_df, outcome_col, 
                         max_runtime_secs = 60, seed = 1234,
                         nfolds_inner = 5, balance_classes = TRUE) {
  
  # Ensure H2O is running
  tryCatch({ h2o::h2o.init(nthreads = -1) }, 
           error = function(e) h2o::h2o.init(nthreads = -1))
  
  # Convert to H2O frames
  train_h2o <- as.h2o(train_df)
  test_h2o  <- as.h2o(test_df)
  
  x <- setdiff(names(train_df), outcome_col)
  y <- outcome_col
  
  train_h2o[, y] <- as.factor(train_h2o[, y])
  test_h2o[,  y] <- as.factor(test_h2o[, y])
  
  # ======================================================
  # 1) Run AutoML (inner CV)
  # ======================================================
  aml <- h2o.automl(
    x = x,
    y = y,
    training_frame = train_h2o,
    nfolds = nfolds_inner,
    max_runtime_secs = max_runtime_secs,
    balance_classes = balance_classes,
    seed = seed,
    keep_cross_validation_predictions = TRUE,
    keep_cross_validation_models = TRUE
  )
  
  # ======================================================
  # 2) Extract leaderboard (all models)
  # ======================================================
  lb <- h2o.get_leaderboard(aml, extra_columns = "ALL")
  model_ids <- as.data.frame(lb$model_id)[,1]
  model_families <- sub("_.*", "", model_ids)
  
  # ======================================================
  # 3) Collect INNER CV metrics per model
  # ======================================================
  inner_cv <- lapply(seq_along(model_ids), function(i) {
    
    mid <- model_ids[i]
    fam <- model_families[i]
    
    m <- h2o.getModel(mid)
    perf <- m@model$cross_validation_metrics_summary
    
    # Convert perf to a clean dataframe
    df <- data.frame(
      metric = rownames(perf),
      mean   = perf[, "mean"],
      sd     = perf[, "sd"],
      row.names = NULL
    )
    
    # Remove GLM-only metrics
    df <- df %>% 
      dplyr::filter(!metric %in% c("null_deviance", "residual_deviance")) %>%
      dplyr::mutate(sd_name = paste0(metric, "_sd"))
    
    # Reshape to wide format: mean columns + sd columns
    df_mean <- df %>%
      select(metric, mean) %>%
      tidyr::pivot_wider(names_from = metric, values_from = mean)
    
    df_sd <- df %>%
      select(sd_name, sd) %>%
      tidyr::pivot_wider(names_from = sd_name, values_from = sd)
    
    df_wide <- dplyr::bind_cols(df_mean, df_sd)
    
    # Combine mean + sd metrics
    df_wide <- dplyr::bind_cols(
      data.frame(
        model_id = mid,
        family   = fam,
        stringsAsFactors = FALSE
      ),
      df_wide
    )
    
    df_wide
  })
  
  inner_cv <- dplyr::bind_rows(inner_cv) 
  
  # ======================================================
  # 4) For each FAMILY → pick best model_id by INNER CV AUC
  # ======================================================
  best_per_family <- inner_cv %>%
    group_by(family) %>%
    slice_max(auc, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # ======================================================
  # 5) Evaluate *each family’s best model* on the OUTER test set
  # ======================================================
  outer_eval <- lapply(seq_len(nrow(best_per_family)), function(i) {
    mid <- best_per_family$model_id[i]
    fam <- best_per_family$family[i]
    
    model <- h2o.getModel(mid)
    perf  <- h2o.performance(model, newdata = test_h2o)
    
    auc_val      <- as.numeric(h2o.auc(perf))
    auprc_val    <- as.numeric(h2o.aucpr(perf))
    logloss_val  <- as.numeric(h2o.logloss(perf))
    mse_val      <- as.numeric(h2o.mse(perf))
    rmse_val     <- as.numeric(h2o.rmse(perf))
    
    data.frame(
      model_id = mid,
      family   = fam,
      auc     = auc_val,
      auprc   = auprc_val,
      logloss = logloss_val,
      mse     = mse_val,
      rmse    = rmse_val
    )
    
  })
  
  outer_eval <- dplyr::bind_rows(outer_eval)
  
  # Return ALL families performance for this fold
  list(
    inner_cv_metrics   = inner_cv,
    best_models_family = best_per_family,
    outer_family_perf  = outer_eval
  )
}

aggregate_families <- function(all_results) {
  
  # ============================================================
  # 1) Aggregate OUTER-FOLD performance (per family)
  # ============================================================
  outer_df <- bind_rows(
    lapply(seq_along(all_results), function(i) {
      x <- all_results[[i]]$outer_family_perf
      x$fold <- i
      x
    })
  )
  
  outer_summary <- outer_df %>%
    group_by(family) %>%
    summarise(
      auc_mean      = mean(auc, na.rm = TRUE),
      auc_median    = median(auc, na.rm = TRUE),
      auc_sd        = sd(auc, na.rm = TRUE),
      
      auprc_mean    = mean(auprc, na.rm = TRUE),
      auprc_median  = median(auprc, na.rm = TRUE),
      auprc_sd      = sd(auprc, na.rm = TRUE),
      
      logloss_mean  = mean(logloss, na.rm = TRUE),
      logloss_median = median(logloss, na.rm = TRUE),
      logloss_sd    = sd(logloss, na.rm = TRUE),
      
      mse_mean      = mean(mse, na.rm = TRUE),
      mse_median    = median(mse, na.rm = TRUE),
      mse_sd        = sd(mse, na.rm = TRUE),
      
      rmse_mean     = mean(rmse, na.rm = TRUE),
      rmse_median   = median(rmse, na.rm = TRUE),
      rmse_sd       = sd(rmse, na.rm = TRUE),
      
      n_folds = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(auc_median))
  
  
  
  # ============================================================
  # 2) Aggregate INNER CV metrics (all metrics + SD)
  # ============================================================
  inner_df <- bind_rows(
    lapply(seq_along(all_results), function(i) {
      x <- all_results[[i]]$inner_cv_metrics
      x$fold <- i
      x
    })
  )
  
  # ---------------------
  # Identify only NON-SD metric columns
  # ---------------------
  metric_cols <- names(inner_df)
  metric_cols <- metric_cols[
    !(metric_cols %in% c("model_id", "family", "fold")) &
      !stringr::str_detect(metric_cols, "_sd$")
  ]
  
  inner_summary <- inner_df %>%
    group_by(family) %>%
    summarise(
      across(
        all_of(metric_cols),
        list(
          mean   = ~mean(.x, na.rm = TRUE),
          median = ~median(.x, na.rm = TRUE),
          sd     = ~sd(.x, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      n_models = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(auc_median))
  
  
  # ============================================================
  # 3) Return everything cleanly
  # ============================================================
  list(
    outer_per_fold = outer_df,
    outer_summary  = outer_summary,
    inner_per_fold = inner_df,
    inner_summary  = inner_summary
  )
}

plot_metric_summary_nested <- function(nested_obj, metric = "auc", statistic = "median") {
  
  # ---------- Determine expected names ----------
  metric_col  <- paste0(metric, "_", statistic)
  sd_col      <- paste0(metric, "_sd")
  
  # ---------- INNER ----------
  has_inner <- metric_col %in% colnames(nested_obj$inner_family_summary)
  
  df_inner <- NULL
  if (has_inner) {
    df_inner <- nested_obj$inner_family_summary %>%
      select(family, value = all_of(metric_col), sd = all_of(sd_col)) %>%
      mutate(type = "CV training")
  }
  
  # ---------- OUTER ----------
  has_outer <- metric_col %in% colnames(nested_obj$outer_family_summary)
  
  df_outer <- NULL
  if (has_outer) {
    df_outer <- nested_obj$outer_family_summary %>%
      select(family, value = all_of(metric_col), sd = all_of(sd_col)) %>%
      mutate(type = "CV test")
  }
  
  if (!has_inner & !has_outer) {
    stop("Metric ", metric, " not found in nested summaries.")
  }
  
  # ---------- Combine ----------
  df <- dplyr::bind_rows(df_inner, df_outer)
  
  # ---------- Plot ----------
  print(
    ggplot(df, aes(x = family, y = value, fill = type)) +
      geom_col(position = position_dodge(width = 0.8)) +
      geom_errorbar(
        aes(ymin = value - sd, ymax = value + sd),
        width = 0.15,
        position = position_dodge(width = 0.8)
      ) +
      scale_y_continuous(
        breaks = seq(0, 1, by = 0.1),
        limits = c(0, 1)
      ) +
      scale_fill_manual(values = c(
        "CV training" = "#1f77b4",
        "CV test"     = "#ff7f0e"
      )) +
      ggtitle(paste("Nested CV performance –", toupper(metric), statistic)) +
      xlab("Family") +
      ylab(metric) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank()
      )
  )
}


pick_best_model_outer <- function(fam_summary) {
  
  outer_per_fold  <- fam_summary$outer_per_fold
  outer_summary   <- fam_summary$outer_summary
  
  # 1) Best family by *outer* performance
  best_family <- outer_summary %>%
    arrange(desc(auc_median)) %>%
    slice(1) %>%
    pull(family)
  
  cat("Winning family based on OUTER AUC:", best_family, "\n")
  
  # 2) Subset models from that family
  fam_models <- outer_per_fold %>%
    filter(family == best_family)
  
  # 3) Pick best model (outer AUC)
  best_model_row <- fam_models %>%
    arrange(desc(auc)) %>%
    slice(1)
  
  best_model_id <- best_model_row$model_id
  
  list(
    winning_family = best_family,
    best_model_id  = best_model_id,
    best_row       = best_model_row
  )
}

aggregate_nested_families <- function(all_results) {

  
  # ============================================================
  # 1) Extract ALL inner + outer metrics across folds×parameters
  # ============================================================
  
  inner_list <- list()
  outer_list <- list()
  idx <- 1
  
  for (fold_i in seq_along(all_results)) {
    
    fold_results <- all_results[[fold_i]]
    
    for (param_i in seq_along(fold_results)) {
      
      res <- fold_results[[param_i]]
      
      # INNER
      inner_df <- res$inner_cv_metrics
      inner_df$fold <- fold_i
      inner_df$parameter <- param_i
      
      # OUTER
      outer_df <- res$outer_family_perf
      outer_df$fold <- fold_i
      outer_df$parameter <- param_i
      
      inner_list[[idx]] <- inner_df
      outer_list[[idx]] <- outer_df
      
      idx <- idx + 1
    }
  }
  
  inner_all <- bind_rows(inner_list)
  outer_all <- bind_rows(outer_list)
  
  
  
  # ============================================================
  # 2) Aggregate OUTER metrics (group by family)
  # ============================================================
  
  outer_metric_cols <- c("auc", "auprc", "logloss", "mse", "rmse")
  
  outer_summary <- outer_all %>%
    group_by(family) %>%
    summarise(
      across(
        all_of(outer_metric_cols),
        list(
          mean   = ~mean(.x, na.rm = TRUE),
          median = ~median(.x, na.rm = TRUE),
          sd     = ~sd(.x, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      n = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(auc_median))
  
  
  
  # ============================================================
  # 3) Aggregate INNER metrics (group by family)
  #    EXCLUDING SD columns (they represent CV fold SD per model)
  # ============================================================
  
  inner_metric_cols <- names(inner_all)
  inner_metric_cols <- inner_metric_cols[
    !(inner_metric_cols %in% c("model_id", "family", "fold", "parameter")) &
      !str_detect(inner_metric_cols, "_sd$")
  ]
  
  inner_summary <- inner_all %>%
    group_by(family) %>%
    summarise(
      across(
        all_of(inner_metric_cols),
        list(
          mean   = ~mean(.x, na.rm = TRUE),
          median = ~median(.x, na.rm = TRUE),
          sd     = ~sd(.x, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      n = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(auc_median))
  
  
  
  # ============================================================
  # 4) Return everything
  # ============================================================
  
  list(
    inner_per_fold     = inner_all,
    inner_summary = inner_summary,
    outer_per_fold     = outer_all,
    outer_summary = outer_summary
  )
}

aggregate_by_parameter <- function(all_results) {
  
  outer_list <- list()
  inner_list <- list()
  
  # ============================================================
  # 1) Collect ALL inner + outer metrics across folds × params
  # ============================================================
  
  for (fold in seq_along(all_results)) {
    for (param in seq_along(all_results[[fold]])) {
      
      res <- all_results[[fold]][[param]]
      
      ## ---- OUTER ----
      df_outer <- res$outer_family_perf
      df_outer$fold      <- fold
      df_outer$parameter <- param
      outer_list[[length(outer_list) + 1]] <- df_outer
      
      ## ---- INNER ----
      df_inner <- res$inner_cv_metrics
      df_inner$fold      <- fold
      df_inner$parameter <- param
      inner_list[[length(inner_list) + 1]] <- df_inner
    }
  }
  
  outer_df <- dplyr::bind_rows(outer_list)
  inner_df <- dplyr::bind_rows(inner_list)
  
  
  # ============================================================
  # 2) Summary per PARAMETER — OUTER METRICS
  # ============================================================
  
  outer_summary <- outer_df %>%
    group_by(parameter) %>%
    summarise(
      auc_mean   = mean(auc, na.rm = TRUE),
      auc_median = median(auc, na.rm = TRUE),
      auc_sd     = sd(auc, na.rm = TRUE),
      
      auprc_mean   = mean(auprc, na.rm = TRUE),
      auprc_median = median(auprc, na.rm = TRUE),
      auprc_sd     = sd(auprc, na.rm = TRUE),
      
      logloss_mean   = mean(logloss, na.rm = TRUE),
      logloss_median = median(logloss, na.rm = TRUE),
      logloss_sd     = sd(logloss, na.rm = TRUE),
      
      mse_mean   = mean(mse, na.rm = TRUE),
      mse_median = median(mse, na.rm = TRUE),
      mse_sd     = sd(mse, na.rm = TRUE),
      
      rmse_mean   = mean(rmse, na.rm = TRUE),
      rmse_median = median(rmse, na.rm = TRUE),
      rmse_sd     = sd(rmse, na.rm = TRUE),
      
      n = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(auc_median))
  
  
  # ============================================================
  # 3) Summary per PARAMETER — INNER CV METRICS
  #    (excluding *_sd columns)
  # ============================================================
  
  # Identify metrics to aggregate (ignore sd, ids, fold, parameter)
  metric_cols_inner <- names(inner_df)
  metric_cols_inner <- metric_cols_inner[
    !(metric_cols_inner %in% c("model_id", "family", "fold", "parameter")) &
      !stringr::str_detect(metric_cols_inner, "_sd$")
  ]
  
  inner_summary <- inner_df %>%
    group_by(parameter) %>%
    summarise(
      across(
        all_of(metric_cols_inner),
        list(
          mean   = ~mean(.x, na.rm = TRUE),
          median = ~median(.x, na.rm = TRUE),
          sd     = ~sd(.x, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      n = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(auc_median))
  
  
  # ============================================================
  # 4) Return all aggregated components
  # ============================================================
  
  return(list(
    outer_per_fold = outer_df,
    outer_summary  = outer_summary,
    inner_per_fold = inner_df,
    inner_summary  = inner_summary
  ))
}

select_best_model_nested <- function(all_results) {
  
  message("==== Selecting Best Parameter → Best Family → Best Model ID ====\n")
  
  # ================================================================
  # (1) Aggregate by parameter → choose BEST PARAMETER
  # ================================================================
  param_aggr <- aggregate_by_parameter(all_results)
  
  param_summary <- param_aggr$outer_summary
  
  best_param <- param_summary %>%
    arrange(desc(auc_median)) %>%
    slice(1) %>%
    pull(parameter)
  
  message("👉 Best PARAMETER: ", best_param, "\n")
  
  
  # ================================================================
  # (2) Subset OUTER & INNER (per-fold) for the best parameter ONLY
  # ================================================================
  outer_df_best_param <- param_aggr$outer_per_fold %>%
    filter(parameter == best_param)
  
  inner_df_best_param <- param_aggr$inner_per_fold %>%
    filter(parameter == best_param)
  
  
  # ================================================================
  # (3) Aggregate OUTER metrics by FAMILY (best parameter only)
  # ================================================================
  outer_family_summary <- outer_df_best_param %>%
    group_by(family) %>%
    summarise(
      auc_mean   = mean(auc, na.rm = TRUE),
      auc_median = median(auc, na.rm = TRUE),
      auc_sd     = sd(auc, na.rm = TRUE),
      
      auprc_mean   = mean(auprc, na.rm = TRUE),
      auprc_median = median(auprc, na.rm = TRUE),
      auprc_sd     = sd(auprc, na.rm = TRUE),
      
      logloss_mean   = mean(logloss, na.rm = TRUE),
      logloss_median = median(logloss, na.rm = TRUE),
      logloss_sd     = sd(logloss, na.rm = TRUE),
      
      mse_mean   = mean(mse, na.rm = TRUE),
      mse_median = median(mse, na.rm = TRUE),
      mse_sd     = sd(mse, na.rm = TRUE),
      
      rmse_mean   = mean(rmse, na.rm = TRUE),
      rmse_median = median(rmse, na.rm = TRUE),
      rmse_sd     = sd(rmse, na.rm = TRUE),
      
      n_folds = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(auc_median))
  
  best_family <- outer_family_summary$family[1]
  message("👉 Best FAMILY for this parameter: ", best_family, "\n")
  
  
  # ================================================================
  # (4) Aggregate INNER CV metrics by FAMILY (best parameter only)
  # ================================================================
  # Identify NON-SD inner metric columns
  inner_metric_cols <- names(inner_df_best_param)
  inner_metric_cols <- inner_metric_cols[
    !(inner_metric_cols %in% c("model_id", "family", "fold", "parameter")) &
      !stringr::str_detect(inner_metric_cols, "_sd$")
  ]
  
  inner_family_summary <- inner_df_best_param %>%
    group_by(family) %>%
    summarise(
      across(
        all_of(inner_metric_cols),
        list(
          mean   = ~mean(.x, na.rm = TRUE),
          median = ~median(.x, na.rm = TRUE),
          sd     = ~sd(.x, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      n_models = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(auc_median))
  
  
  # ================================================================
  # (5) Select Best MODEL ID inside best family & parameter
  # ================================================================
  best_model_row <- outer_df_best_param %>%
    filter(family == best_family) %>%
    arrange(desc(auc)) %>%
    slice(1)
  
  best_model_id <- best_model_row$model_id
  
  message("👉 Best MODEL ID: ", best_model_id, "\n")
  message("============================================================\n")
  
  
  # ================================================================
  # (6) Return everything
  # ================================================================
  list(
    best_parameter = best_param,
    best_family = best_family,
    best_model_id = best_model_id,
    
    parameter_summary = param_summary,
    
    outer_family_summary = outer_family_summary,
    inner_family_summary = inner_family_summary,
    
    outer_df_best_param = outer_df_best_param,
    inner_df_best_param = inner_df_best_param,
    
    best_model_row = best_model_row
  )
}

