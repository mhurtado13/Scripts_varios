
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
    test_predictions_all_folds = list()
    
    for (fold_i in seq_along(result_files)) {
      
      result = readRDS(result_files[[fold_i]])
      
      models_all_params <- vector("list", length(result))
      test_predictions_all_params = list()
      for (parameter_i in seq_along(result)) {
        
        train_df <- result[[parameter_i]][["train_data"]]
        test_df  <- result[[parameter_i]][["test_data"]]
        
        test_df = cbind(test_df, target = result[[parameter_i]][["obs_test"]])
        
        res <- run_h2o_fold(train_df, test_df, "target", 
                            nfolds_inner = 0, max_runtime_secs = 60)
        
        # Store results
        models_all_params[[parameter_i]] <- res
        test_predictions_all_params[[parameter_i]] <- res$test_predictions
      }
      
      models_all_folds[[fold_i]] <- models_all_params
      test_predictions_all_folds[[fold_i]] <- test_predictions_all_params
    }
    
    nested_result <- select_best_parameter_family_model(models_all_folds, result_files, test_predictions_all_folds)

    grDevices::pdf(paste0("Results/AUROC_CV_", file_name, ".pdf"), width = 10)
    plot_metric_summary(nested_result, metric = "auc", statistic = "median")
    dev.off()
    
    grDevices::pdf(paste0("Results/AUPRC_CV_", file_name, ".pdf"), width = 10)
    plot_metric_summary(nested_result, metric = "auprc", statistic = "median")
    dev.off()

    h2o::h2o.init(nthreads = -1, bind_to_localhost = TRUE)
    custom_output <- do.call(fold_construction_fun,
                             c(list(data = train_data, bestune = result[[nested_result$best_parameter]][["params"]]), fold_construction_args_fixed))
    h2o::h2o.shutdown(prompt = FALSE)
    
    output = list(nested_result, custom_output)
  }
  
  return(output)
  
}


run_h2o_fold <- function(train_df, test_df, outcome_col, 
                         max_runtime_secs = 60, seed = 1234,
                         nfolds_inner = 0, balance_classes = TRUE) {
  
  train_preds <- list()
  test_preds  <- list()
  
  h2o::h2o.init(
    nthreads = -1,
    bind_to_localhost = TRUE
  )
  
  on.exit({
    h2o::h2o.shutdown(prompt = FALSE)
    Sys.sleep(3)
  }, add = TRUE)
  
  train_h2o <- as.h2o(train_df)
  test_h2o  <- as.h2o(test_df)
  
  x <- setdiff(names(train_df), outcome_col)
  y <- outcome_col
  
  train_h2o[, y] <- as.factor(train_h2o[, y])
  test_h2o[,  y] <- as.factor(test_h2o[, y])
  
  aml <- h2o.automl(
    x = x,
    y = y,
    training_frame = train_h2o,
    nfolds = nfolds_inner,
    max_runtime_secs = max_runtime_secs,
    balance_classes = balance_classes,
    seed = seed,
    keep_cross_validation_predictions = (nfolds_inner > 1),
    keep_cross_validation_models = (nfolds_inner > 1)
  )
  
  lb <- h2o.get_leaderboard(aml, extra_columns = "ALL")
  lb <- as.data.frame(lb)
  
  best_per_family <- lb %>%
    mutate(family = sub("_.*", "", model_id)) %>%
    group_by(family) %>%
    slice(1) %>%
    ungroup() %>%
    select(model_id, family)
  
  # -------- Save each best model per family --------
  dir.create("Results/ML_models", recursive = TRUE, showWarnings = FALSE)
  for (i in seq_len(nrow(best_per_family))) {
    model <- h2o.getModel(best_per_family$model_id[i])
    h2o.saveModel(model, path = "Results/ML_models", force = TRUE)
  }
  
  # ----------- TEST predictions (meta-testing features) ----------
  for (i in seq_len(nrow(best_per_family))) {
    mid  <- best_per_family$model_id[i]
    fam  <- best_per_family$family[i]
    model <- h2o.getModel(mid)
    
    pred_test <- as.data.frame(h2o.predict(model, test_h2o))[ , 3]
    test_preds[[fam]] <- pred_test
  }
  
  # ----------- OUTER METRICS ------------
  outer_eval <- lapply(seq_len(nrow(best_per_family)), function(i) {
    mid <- best_per_family$model_id[i]
    fam <- best_per_family$family[i]
    
    model <- h2o.getModel(mid)
    perf  <- h2o.performance(model, newdata = test_h2o)
    
    data.frame(
      model_id = mid,
      family   = fam,
      auc      = as.numeric(h2o.auc(perf)),
      auprc    = as.numeric(h2o.aucpr(perf)),
      logloss  = as.numeric(h2o.logloss(perf)),
      mse      = as.numeric(h2o.mse(perf)),
      rmse     = as.numeric(h2o.rmse(perf))
    )
  })
  
  outer_eval <- dplyr::bind_rows(outer_eval)
  
  return(list(
    outer_family_perf  = outer_eval,
    test_predictions   = test_preds     # <-- used to evaluate meta-learner
  ))
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

plot_metric_summary <- function(nested_obj, metric = "auc", statistic = "median") {
  
  # Expected column names such as "auc_median", "auc_sd"
  metric_col <- paste0(metric, "_", statistic)
  sd_col     <- paste0(metric, "_sd")
  
  # -------------------------------------------------------------
  #   Only OUTER summary exists now → nested_obj$family_summary
  # -------------------------------------------------------------
  df_outer <- nested_obj$family_summary
  
  # Validate required columns
  if (!metric_col %in% colnames(df_outer) || !sd_col %in% colnames(df_outer)) {
    stop("Metric ", metric, " not found in nested_obj$family_summary.")
  }
  
  # Prepare DF for plotting
  df <- df_outer %>%
    dplyr::select(
      family,
      value = dplyr::all_of(metric_col),
      sd    = dplyr::all_of(sd_col)
    ) %>%
    dplyr::mutate(type = "Performance")
  
  # -------------------------------------------------------------
  #   Plot
  # -------------------------------------------------------------
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
      scale_fill_manual(values = c("Performance" = "#1f77b4")) +
      ggtitle(paste0("Performance by ML models family - ", metric_col)) +
      xlab("ML model family") +
      ylab(metric) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        legend.position = "none"
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
  
  # ============================================================
  # 1) Collect OUTER metrics across folds × params
  # ============================================================
  for (fold in seq_along(all_results)) {
    for (param in seq_along(all_results[[fold]])) {
      
      res <- all_results[[fold]][[param]]
      
      ## ---- OUTER ----
      df_outer <- res$outer_family_perf
      df_outer$fold      <- fold
      df_outer$parameter <- param
      outer_list[[length(outer_list) + 1]] <- df_outer
      
      ## ---- NO INNER CV ----
      ## nfolds_inner = 0 → inner_cv_metrics = NULL
    }
  }
  
  outer_df <- dplyr::bind_rows(outer_list)
  
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

      .groups = "drop"
    ) %>%
    arrange(desc(auc_median))
  
  # ============================================================
  # 3) Return only outer summaries
  # ============================================================
  return(list(
    outer_per_fold = outer_df,
    outer_summary  = outer_summary
  ))
}


select_best_parameter_family_model <- function(models, result_files, test_predictions_all_folds) {
  
  ############################################## Aggregate by parameter
  param_aggr <- aggregate_by_parameter(models)
  
  ############################################## Ensemble stacking models
  
  test_preds <- extract_predictions_all_params(test_predictions_all_folds) ## Extract predictions from base models
   
  stack_res <- build_stacking_matrices_all_params(test_preds, result_files) ## Build stacking matrices across params
  
  stacking_outer <- train_and_evaluate_meta_learner_all_params( ## Train and test metalearners
    stacking_df_list = stack_res,
    families = names(test_preds[[1]][[1]])
  )
  
  param_aggr$outer_summary = rbind(param_aggr$outer_summary, stacking_outer$outer_summary)  ## Join CV results with base models
  param_aggr$outer_per_fold = rbind(param_aggr$outer_per_fold, stacking_outer$outer_per_fold) ## Join CV results with base models
  
  ############################################## Choose best parameter
  
  param_summary <- param_aggr$outer_summary
  
  best_param <- param_summary %>%
    arrange(desc(auc_median)) %>%
    slice(1) %>%
    pull(parameter)
  
  ############################################## Aggregate OUTER metrics by FAMILY

  outer_df_best_param <- param_aggr$outer_per_fold %>%
    filter(parameter == best_param)
  
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
  
  ############################################## Select Best MODEL ID inside best family

  best_model_row <- outer_df_best_param %>%
    filter(family == best_family) %>%
    arrange(desc(auc)) %>%
    slice(1)
  
  best_model_id <- best_model_row$model_id
  
  message("============================================================\n")
  
  message("==== Selecting Best Parameter → Best Family → Best Model ID ====\n")
  message("Best PARAMETER: ", best_param, "\n")
  message("Best FAMILY for this parameter: ", best_family, "\n")
  message("Best MODEL ID: ", best_model_id, "\n")
  
  return(list(
    best_parameter = best_param,
    best_family    = best_family,
    best_model_id  = best_model_id,
    
    parameter_summary  = param_summary,
    family_summary = outer_family_summary,
    
    per_fold_results_best_param = outer_df_best_param
  ))
}


get_free_h2o_port <- function(min = 15000, max = 60000) {
  repeat {
    p <- sample(min:max, 1)
    
    ok <- try({
      # H2O needs p and p+1 FREE
      con1 <- socketConnection("127.0.0.1", port = p,     open = "r+", server = TRUE)
      con2 <- socketConnection("127.0.0.1", port = p + 1, open = "r+", server = TRUE)
      close(con1); close(con2)
      TRUE
    }, silent = TRUE)
    
    if (identical(ok, TRUE)) return(p)
  }
}

extract_predictions_all_params <- function(predictions_all_folds) {
  
  # Helper: extract model family from a model_id
  extract_family <- function(name) sub("_.*", "", name)
  
  out <- vector("list", length(predictions_all_folds))
  
  for (fold_i in seq_along(predictions_all_folds)) {
    
    fold_pred_list <- predictions_all_folds[[fold_i]]
    n_params <- length(fold_pred_list)
    
    # Store predictions per parameter
    param_list <- vector("list", n_params)
    
    for (param_i in seq_len(n_params)) {
      
      preds_param <- fold_pred_list[[param_i]]
      
      # Convert model IDs to families
      new_names <- sapply(names(preds_param), extract_family)
      names(preds_param) <- new_names
      
      # Store
      param_list[[param_i]] <- preds_param
    }
    
    out[[fold_i]] <- param_list
  }
  
  return(out)
}



build_stacking_matrices_all_params <- function(preds_all_params, result_files) {
  
  K <- length(preds_all_params)                      # number of folds
  P <- length(preds_all_params[[1]])                 # number of parameters
  
  out <- vector("list", P)
  
  for (param_i in seq_len(P)) {
    
    stacking_rows <- list()
    
    # Get consistent family names (from first fold)
    families <- names(preds_all_params[[1]][[param_i]])
    
    # Loop over folds
    for (fold_i in seq_len(K)) {
      
      preds_fold_param <- preds_all_params[[fold_i]][[param_i]]
      
      # Ensure consistent ordering
      pred_mat <- sapply(families, function(f) preds_fold_param[[f]])
      colnames(pred_mat) <- families
      
      # Extract correct target for this fold + param
      result <- readRDS(result_files[[fold_i]])
      target <- result[[param_i]][["obs_test"]]
      
      # Build DF
      df_fold <- data.frame(
        fold = fold_i,
        pred_mat,
        target = target,
        row.names = NULL
      )
      
      stacking_rows[[fold_i]] <- df_fold
    }
    
    # Bind all folds
    out[[param_i]] <- dplyr::bind_rows(stacking_rows)
  }
  
  return(out)  # list of stacking_df per parameter
}

train_meta_learner <- function(stacking_df, families) {
  
  # Start a fresh H2O cluster
  h2o.init(
    nthreads = -1,
    bind_to_localhost = TRUE
  )
  
  # Ensure clean shutdown when finished
  on.exit({
    h2o::h2o.shutdown(prompt = FALSE)
    Sys.sleep(2)
  }, add = TRUE)
  
  # Convert to H2O
  df_h2o <- as.h2o(stacking_df)
  
  x <- families
  y <- "target"
  
  # Ensure target is categorical
  df_h2o[, y] <- as.factor(df_h2o[, y])
  
  # Train meta-learner (GLM0 = simple logistic regression)
  meta_model <- h2o.glm(
    x = x,
    y = y,
    training_frame = df_h2o,
    family = "binomial",
    lambda = 0,
    alpha = 0
  )
  
  return(meta_model)
}

evaluate_meta_learner <- function(meta_model, preds_best, result_files, best_param) {
  
  h2o::h2o.init(
    nthreads = -1,
    bind_to_localhost = TRUE
  )
  
  on.exit({
    h2o::h2o.shutdown(prompt = FALSE)
    Sys.sleep(3)
  }, add = TRUE)
  
  # ------ Collect ALL OOF predictions & labels ------
  all_pred_rows <- list()
  all_true_rows <- list()
  
  for (fold_i in seq_along(preds_best)) {
    
    preds_fold <- preds_best[[fold_i]]
    families   <- names(preds_fold)
    
    # prediction matrix for this fold
    pred_mat <- sapply(families, function(f) preds_fold[[f]])
    all_pred_rows[[fold_i]] <- as.data.frame(pred_mat)
    
    # true labels for this fold
    result <- readRDS(result_files[[fold_i]])
    true_y <- result[[best_param]]$obs_test
    
    all_true_rows[[fold_i]] <- data.frame(y = true_y)
  }
  
  # ---- Bind ALL folds together (OOF dataset) ----
  pred_df <- do.call(rbind, all_pred_rows)
  true_df <- do.call(rbind, all_true_rows)
  colnames(true_df) <- 'target'
  
  # Convert to H2O
  pred_h2o <- as.h2o(pred_df)
  true_h2o <- as.h2o(true_df)
  
  # ------------------------------------------------
  # **Single evaluation on ALL OOF data**
  # ------------------------------------------------
  perf <- h2o.performance(
    meta_model,
    newdata = h2o.cbind(pred_h2o, true_h2o)
  )
  
  # ---- Return AutoML-style summary row ----
  data.frame(
    model_id = "STACK_META_LEARNER",
    family   = "Stacking",
    auc      = as.numeric(h2o.auc(perf)),
    auprc    = as.numeric(h2o.aucpr(perf)),
    logloss  = as.numeric(h2o.logloss(perf)),
    mse      = as.numeric(h2o.mse(perf)),
    rmse     = as.numeric(h2o.rmse(perf)),
    parameter = best_param
  )
}

train_and_evaluate_meta_learner_all_params <- function(stacking_df_list, families){
  
  h2o.init(nthreads = -1, bind_to_localhost = TRUE)
  on.exit({ h2o.shutdown(prompt = FALSE); Sys.sleep(2) }, add = TRUE)
  
  P <- length(stacking_df_list)  # number of parameters
  
  outer_per_fold  <- list()
  outer_summary   <- list()
  
  # Create folder for saving metalearners
  save_dir <- "Results/ML_models"
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  #====================================================================
  #   LOOP OVER PARAMETERS
  #====================================================================
  for (param_i in seq_len(P)) {
    
    stacking_df <- stacking_df_list[[param_i]]
    K <- length(unique(stacking_df$fold))
    
    fold_results <- list()
    
    #---------------------------------------------------------
    #    1) TRAIN/TEST META-LEARNER FOR EACH FOLD
    #---------------------------------------------------------
    for (fold_i in seq_len(K)) {
      
      # ---- TRAIN (fold != i) ----
      train_df  <- stacking_df[stacking_df$fold != fold_i, ]
      train_h2o <- as.h2o(train_df)
      train_h2o[, "target"] <- as.factor(train_h2o[, "target"])
      
      meta_model <- h2o.glm(
        x = families,
        y = "target",
        training_frame = train_h2o,
        family = "binomial",
        lambda = 0,
        alpha = 0
      )
      
      # ---- SAVE METALERNER MODEL ----
      save_name <- paste0("StackedEnsemble_param_", param_i, "_fold", fold_i)
      h2o.saveModel(meta_model, path = save_dir, force = TRUE, filename = save_name)
      
      # ---- TEST (fold == i) ----
      test_df  <- stacking_df[stacking_df$fold == fold_i, ]
      test_h2o <- as.h2o(test_df)
      
      perf <- h2o.performance(meta_model, newdata = test_h2o)
      
      fold_results[[fold_i]] <- data.frame(
        model_id  = save_name,
        family    = "StackedEnsemble",
        auc       = as.numeric(h2o.auc(perf)),
        auprc     = as.numeric(h2o.aucpr(perf)),
        logloss   = as.numeric(h2o.logloss(perf)),
        mse       = as.numeric(h2o.mse(perf)),
        rmse      = as.numeric(h2o.rmse(perf)),
        fold      = fold_i,
        parameter = param_i
      )
    }
    
    # Combine into a single DF for this parameter
    df_param <- dplyr::bind_rows(fold_results)
    outer_per_fold[[param_i]] <- df_param
    
    #---------------------------------------------------------
    #    2) SUMMARY (matches nested_result$outer_family_summary)
    #---------------------------------------------------------
    outer_summary[[param_i]] <- df_param %>%
      summarise(
        parameter     = param_i,
        
        auc_mean      = mean(auc, na.rm = TRUE),
        auc_median    = median(auc, na.rm = TRUE),
        auc_sd        = sd(auc, na.rm = TRUE),
        
        auprc_mean    = mean(auprc, na.rm = TRUE),
        auprc_median  = median(auprc, na.rm = TRUE),
        auprc_sd      = sd(auprc, na.rm = TRUE),
        
        logloss_mean   = mean(logloss, na.rm = TRUE),
        logloss_median = median(logloss, na.rm = TRUE),
        logloss_sd     = sd(logloss, na.rm = TRUE),
        
        mse_mean       = mean(mse, na.rm = TRUE),
        mse_median     = median(mse, na.rm = TRUE),
        mse_sd         = sd(mse, na.rm = TRUE),
        
        rmse_mean      = mean(rmse, na.rm = TRUE),
        rmse_median    = median(rmse, na.rm = TRUE),
        rmse_sd        = sd(rmse, na.rm = TRUE)
        
      )
  }
  
  return(list(
    outer_per_fold = dplyr::bind_rows(outer_per_fold),
    outer_summary  = dplyr::bind_rows(outer_summary)
  ))
}
