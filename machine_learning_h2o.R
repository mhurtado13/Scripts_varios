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
    
    fam_summary <- aggregate_nested_families(models_all_folds)
    plot_metric_summary(fam_summary, metric = "auc", statistic = "median")
    best_info <- pick_best_model_outer(fam_summary)
    
    output = list(fam_summary, best_info)
  }
  
  return(output)
  
}
