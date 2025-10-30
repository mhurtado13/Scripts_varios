
suppressPackageStartupMessages({
  library(censored)
  library(tidyverse)
  library(tidymodels)
  library(survival)
  library(cowplot)
  library(future)
  library(recipes)
  library(doParallel)
})

#' Get default hyperparameter grids for supported survival models
#'
#' This function generates default hyperparameter grids for various survival models
#' compatible with the tidymodels framework. It uses the default parameter ranges
#' from the {dials} package and produces a sequence of evenly spaced values
#' within those ranges for each tunable hyperparameter.
#'
#' The function supports models such as Cox proportional hazards (regular and penalized),
#' parametric survival regression, decision trees, bagging, random forests, and gradient boosting.
#' For models without tunable parameters (e.g., classic Cox or AFT models),
#' the function returns \code{NULL}.
#'
#' @param model_name Character string specifying the model name.
#'   Supported options include:
#'   \itemize{
#'     \item \code{"cox_ph_survival"} – Classic Cox proportional hazards model
#'     \item \code{"proportional_hazards_glmnet"} – Penalized Cox (LASSO / Elastic Net)
#'     \item \code{"survreg_flexsurv"} – Parametric accelerated failure time (AFT)
#'     \item \code{"decision_tree_partykit"} – Single survival tree
#'     \item \code{"bag_tree_rpart"} – Bagged CART survival trees
#'     \item \code{"rand_forest_partykit"} – Random survival forest (ctree-based)
#'     \item \code{"rand_forest_aorsf"} – Oblique random survival forest
#'     \item \code{"boost_tree_mboost"} – Gradient boosting for survival
#'   }
#' @param train_x Optional data frame or matrix of training predictors.
#'   This is required for parameters that depend on the number of features,
#'   such as \code{mtry}, which determines the number of variables randomly
#'   sampled at each split in tree-based models.
#' @param levels Integer specifying how many values to generate per hyperparameter.
#'   Defaults to \code{5}. Must be at least 2.
#'
#' @return A named list of hyperparameter grids.
#'   Each element is a numeric vector of sampled values for that parameter.
#'   Returns \code{NULL} for models without tunable hyperparameters.
#'
#' @details
#' The helper function \code{vs()} internally calls \code{dials::value_seq()}
#' to generate evenly spaced sequences of parameter values across their default ranges.
#' For data-dependent parameters (like \code{mtry}), \code{dials::finalize()}
#' is used to compute appropriate limits based on \code{train_x}.
#'
#' @examples
#' # Example 1: Random forest with feature-dependent hyperparameter
#' get_default_hyperparams("rand_forest_partykit", train_x = iris[, 1:4])
#'
#' # Example 2: Penalized Cox model with default parameter ranges
#' get_default_hyperparams("proportional_hazards_glmnet")
#'
get_default_hyperparams <- function(model_name, train_x = NULL, levels = 5) {
  stopifnot(levels >= 2)
  
  # --- helpers -----------------------------------------------------------
  vs <- function(param) dials::value_seq(param, n = levels)  # evenly spaced across default range
  p  <- if (!is.null(train_x)) ncol(train_x) else NA_integer_
  
  # --- model-specific parameter grids -----------------------------------
  if (model_name == "cox_ph_survival") {
    # Classic Cox model (no tunable hyperparameters)
    return(NULL)
    
  } else if (model_name == "proportional_hazards_glmnet") {
    # Penalized Cox regression (LASSO/Elastic Net)
    return(list(
      penalty = vs(dials::penalty()),  # regularization strength (log10 scale)
      mixture = vs(dials::mixture())   # elastic net mixing parameter (0 = ridge, 1 = lasso)
    ))
    
  } else if (model_name == "survreg_flexsurv") {
    # Parametric accelerated failure time (AFT) model
    # No tunable hyperparameters in parsnip interface
    return(NULL)
    
  } else if (model_name == "decision_tree_partykit") {
    # Decision tree model hyperparameters
    return(list(
      cost_complexity = vs(dials::cost_complexity()),  # complexity/pruning penalty
      tree_depth      = vs(dials::tree_depth()),       # maximum depth of the tree
      min_n           = vs(dials::min_n())             # minimum number of samples per node
    ))
    
  } else if (model_name == "bag_tree_rpart") {
    # Bagging of decision trees (rpart engine)
    return(list(
      trees = vs(dials::trees())  # number of bootstrap trees
    ))
    
  } else if (model_name == "rand_forest_partykit") {
    # Random survival forest using the partykit engine
    vals <- list(
      trees = vs(dials::trees()),  # number of trees in the forest
      min_n = vs(dials::min_n())   # minimum samples per node
    )
    # mtry depends on number of predictors, so finalize using train_x
    if (!is.na(p)) vals$mtry <- vs(dials::finalize(dials::mtry(), train_x))
    return(vals)
    
  } else if (model_name == "rand_forest_aorsf") {
    # Oblique random survival forest (aorsf engine)
    vals <- list(
      trees = vs(dials::trees()),  # number of trees
      min_n = vs(dials::min_n())   # minimum samples per node
    )
    # finalize mtry if training data provided
    if (!is.na(p)) vals$mtry <- vs(dials::finalize(dials::mtry(), train_x))
    return(vals)
    
  } else if (model_name == "boost_tree_mboost") {
    # Gradient boosting model for survival analysis
    return(list(
      trees          = vs(dials::trees()),           # number of boosting iterations
      min_n          = vs(dials::min_n()),           # minimum number of samples per node
      tree_depth     = vs(dials::tree_depth()),      # depth of individual trees
      learn_rate     = vs(dials::learn_rate()),      # learning rate (step size)
      loss_reduction = vs(dials::loss_reduction()),  # minimum loss reduction per split
      sample_size    = vs(dials::sample_size(range = c(0, 1))),  # fraction of samples per iteration
      stop_iter      = vs(dials::stop_iter())        # early stopping iterations
    ))
    
  } else {
    stop("Unsupported model name: ", model_name)
  }
}


unregister_dopar <- function() {
  if (!is.null(foreach::getDoParRegistered())) {
    # switch back to sequential backend
    foreach::registerDoSEQ()
    gc()
  }
}

# =============================================================================
# Function: machine_learning_custom()
# Purpose : Train and evaluate a single survival model (Cox, AFT, Random Forest,
#           Boosting, etc.) using the tidymodels framework.
#           This function supports multiple model engines, user-defined
#           hyperparameters, and computes the Concordance Index (C-index)
#           as the main performance metric.
# =============================================================================

machine_learning_custom <- function(df_train, df_test,
                                    outcome_col, event_col = NULL,
                                    model, models_hyperparameters){
  # ---------------------------------------------------------------------------
  # Step 3a: Define the model specification based on the model name
  # ---------------------------------------------------------------------------
  # Each model is created using the {parsnip} interface, which provides a unified
  # syntax for specifying models across various engines.
  # All models are defined with mode = "censored regression", which is required
  # for survival analysis tasks.
  #
  # The engine determines the computational backend used for fitting:
  #   - "survival" → Classical Cox or AFT model from {survival}
  #   - "glmnet"   → Penalized regression (LASSO / Elastic Net)
  #   - "partykit" → Tree-based survival models (ctree framework)
  #   - "aorsf"    → Oblique random survival forest (fast, optimized)
  #   - "rpart"    → CART-style bagged trees
  #   - "mboost"   → Gradient boosting for survival analysis
  #
  if (model == "cox_ph_survival") {
    # Standard Cox proportional hazards model (no tunable hyperparameters)
    model_spec <- parsnip::proportional_hazards(
      mode = "censored regression",
      engine = "survival"
    )
    
  } else if (model == "proportional_hazards_glmnet") {
    # Penalized Cox model — introduces two tunable parameters:
    #   penalty (λ): regularization strength
    #   mixture (α): elastic net mixing (0 = ridge, 1 = lasso)
    model_spec <- parsnip::proportional_hazards(
      penalty = tune(),
      mixture = tune(),
      mode = "censored regression",
      engine = "glmnet"
    )
    
  } else if (model == "survreg_flexsurv") {
    # Parametric accelerated failure time (AFT) model via flexsurv engine
    model_spec <- parsnip::survival_reg(
      mode = "censored regression",
      engine = "flexsurv"
    )
    
  } else if (model == "rand_forest_partykit") {
    # Random survival forest using conditional inference trees
    model_spec <- parsnip::rand_forest(
      trees = tune(),
      mode = "censored regression",
      engine = "partykit"
    )
    
  } else if (model == "rand_forest_aorsf") {
    # Random survival forest using oblique splits (aorsf package)
    model_spec <- parsnip::rand_forest(
      trees = tune(),
      mode = "censored regression",
      engine = "aorsf"
    )
    
  } else if (model == "decision_tree_partykit") {
    # Single conditional inference survival tree
    model_spec <- parsnip::decision_tree(
      mode = "censored regression",
      engine = "partykit"
    )
    
  } else if (model == "bag_tree_rpart") {
    # Bagged ensemble of CART survival trees
    model_spec <- parsnip::bag_tree(
      mode = "censored regression",
      engine = "rpart"
    )
    
  } else if (model == "boost_tree_mboost") {
    # Gradient boosting model for censored survival data
    model_spec <- parsnip::boost_tree(
      trees = tune(),
      mode = "censored regression",
      engine = "mboost"
    )
    
  } else {
    # Catch any unsupported model names
    stop("Unsupported model: ", model)
  }
  
  
  # ---------------------------------------------------------------------------
  # Step 3b: Apply user-provided hyperparameters
  # ---------------------------------------------------------------------------
  # Hyperparameters are provided as a list of parameter sets, typically created
  # from a tuning grid. Example:
  #   models_hyperparameters = list(list(trees = 500, min_n = 10))
  #
  # The parsnip::set_args() function updates the model specification with
  # these custom values. Wrapping the call in do.call() allows passing an
  # arbitrary number of named arguments dynamically.
  #
  if (!is.null(models_hyperparameters)) {
    model_spec <- do.call(
      parsnip::set_args,
      c(list(object = model_spec), models_hyperparameters[[1]])
    )
  }
  
  # ---------------------------------------------------------------------------
  # Step 3c: Fit the model using a tidymodels workflow
  # ---------------------------------------------------------------------------
  # A workflow combines a model specification with a formula or recipe.
  # For survival models, the response must be defined as:
  #   Surv(time, event) ~ predictors
  #
  # This structure is recognized by survival engines and allows them to
  # compute censored regression models correctly.
  #
  formula_model <- as.formula(
    paste0("Surv(", outcome_col, ", ", event_col, ") ~ .")
  )
  
  wf <- workflows::workflow() %>%
    workflows::add_model(model_spec) %>%
    workflows::add_formula(formula_model)
  
  # Fit the model to the training data
  fitted <- parsnip::fit(wf, data = df_train)
  
  # ---------------------------------------------------------------------------
  # Step 3d: Generate predictions on test data
  # ---------------------------------------------------------------------------
  # The code tries different prediction "types" depending on what the fitted model supports.
  #
  # 1. "linear_pred"  → Risk score or log hazard (e.g., Cox model)
  #      - Higher values = higher risk (shorter survival)
  #      - Direction: higher = worse outcome
  #
  # 2. "time"         → Expected survival time (e.g., parametric AFT models)
  #      - Higher values = longer expected lifetime
  #      - Direction: higher = better outcome (will be reversed later)
  #
  # 3. "survival"     → Survival probability at a given evaluation time
  #      - Higher values = more likely to survive (less risk)
  #      - Direction: higher = better outcome (will be reversed later)
  #
  # The try–catch structure ensures the code works for any survival model:
  # - Try risk-based prediction first
  # - If not supported, try time-based prediction
  # - If that fails, fall back to survival probability at a fixed eval_time
  # ---------------------------------------------------------------------------
  
  pred_type <- NULL  # track which type succeeded
  
  preds <- tryCatch({
    pred_type <<- "linear_pred"
    stats::predict(fitted, df_test, type = "linear_pred")
  }, error = function(e1) {
    tryCatch({
      pred_type <<- "time"
      stats::predict(fitted, df_test, type = "time")
    }, error = function(e2) {
      # Fallback: use survival probabilities at the median observed time
      eval_time <- stats::median(df_train[[outcome_col]], na.rm = TRUE)
      pred_type <<- "survival"
      stats::predict(fitted, df_test, type = "survival", eval_time = eval_time)
    })
  })
  
  # ---------------------------------------------------------------------------
  # Step 3e: Standardize prediction output into a numeric tibble
  # ---------------------------------------------------------------------------
  # Prediction outputs vary by engine:
  #   - Some return numeric vectors (risk scores)
  #   - Others return matrices (e.g., survival curves across time)
  #   - Some may return lists of probabilities
  #
  # This section standardizes all prediction outputs into a tibble with a
  # single numeric column `.pred`, computing the median when multiple
  # predictions per sample exist.
  #
  
  if ("matrix" %in% class(preds[[1]])) {
    preds <- tibble::tibble(.pred = apply(as.data.frame(preds[[1]]), 1, median, na.rm = TRUE))
  } else if (!is.numeric(preds[[1]])) {
    preds <- tibble::tibble(.pred = apply(as.matrix(preds), 1, median, na.rm = TRUE))
  } else {
    pred_col <- names(preds)[1]
    preds <- preds %>% dplyr::rename(.pred = dplyr::all_of(pred_col))
  }
  
  # ---------------------------------------------------------------------------
  # Step 3f: Ensure direction consistency (higher = higher risk)
  # ---------------------------------------------------------------------------
  
  # Reverse predictions if type implies "more = better survival"
  if (is.null(pred_type)) {
    warning("⚠️ No successful prediction type for model: ", model)
  } else if (pred_type %in% c("time", "survival")) {
    preds$.pred <- -preds$.pred
  }
  
  # ---------------------------------------------------------------------------
  # Step 3g: Evaluate model performance using the Concordance Index (C-index)
  # ---------------------------------------------------------------------------
  # The C-index measures how well the model ranks survival times relative to
  # true outcomes (i.e., discrimination ability). It ranges from 0.5 (random)
  # to 1 (perfect discrimination).
  #
  # The function `censored::concordance_survival_vec()` computes the metric.
  # If the computation fails, the value is safely returned as NA.
  #
  metric_val <- tryCatch({
    yardstick::concordance_survival_vec(
      truth = survival::Surv(df_test[[outcome_col]], df_test[[event_col]]),
      estimate = preds$.pred
    )
  }, error = function(e) {
    message("Concordance calculation failed: ", e$message)
    NA_real_
  })
  
  # ---------------------------------------------------------------------------
  # Step 3h: Return the model name and computed C-index
  # ---------------------------------------------------------------------------
  # The function returns a tibble containing the model name and the computed
  # C-index, making it easy to aggregate and compare across models and folds.
  #
  return(tibble::tibble(model = model, c_index = metric_val))
}

# =============================================================================
# Function: cross_validation_custom()
# Purpose : Perform *nested cross-validation* with hyperparameter tuning for
#           multiple survival models using tidymodels.
#           For each model and hyperparameter configuration, the function:
#             1. Trains across multiple folds
#             2. Computes the C-index (concordance)
#             3. Aggregates median and MAD across folds
#           Finally, it identifies the best-performing model and hyperparameters.
# =============================================================================

cross_validation_custom <- function(df_features, df_outcome, 
                                    outcome_col, event_col, 
                                    ml_options = list(nb_folds = 5,
                                                      nb_repeats = 1,
                                                      ncores = parallel::detectCores() - 1,
                                                      LODO = FALSE,
                                                      batch_id = NULL),
                                    file_name = NULL){
  
  ### Implement parallelization over folds
  # Set the number of CPU cores for parallel processing, as specified in ml_options.
  ncores <- ml_options$ncores
  
  # ---------------------------------------------------------------------------
  # Step 2: Combine predictors and outcomes
  # ---------------------------------------------------------------------------
  # Merge the feature matrix (X) and the outcome data (Y) into a single dataset.
  # This ensures that both survival time and event indicators are correctly aligned
  # for each individual prior to creating cross-validation splits.
  #
  df_all <- df_features %>% 
    dplyr::bind_cols(df_outcome)
  
  # ---------------------------------------------------------------------------
  # Step 3: Create v-fold cross-validation splits  (with optional LODO stratification)
  # ---------------------------------------------------------------------------
  # Use the {rsample} package to generate stratified K-fold CV partitions.
  # Stratification ensures that the proportion of events (1) vs. censored (0)
  # is roughly preserved across folds, improving stability of C-index estimates.
  # Each fold contains:
  #   - an analysis set (training data)
  #   - an assessment set (validation data)
  # If `nb_repeats` > 1, the v-fold partitioning is repeated multiple times
  # for more stable performance estimates.
  # If LODO = FALSE → standard stratified by event rate
  # If LODO = TRUE  → stratified by both cohort (batch_id) and event indicator
  
  if (isTRUE(ml_options$LODO)) {
    if (is.null(ml_options$batch_id) || !ml_options$batch_id %in% names(df_all)) {
      stop("When LODO = TRUE, you must provide 'batch_id' as a column name in df_all.")
    }
    
    batch_col <- ml_options$batch_id

    # Create composite stratification variable: cohort × event
    df_all <- df_all %>%
      dplyr::mutate(strata = interaction(.data[[batch_col]], .data[[event_col]], drop = TRUE))
    
    folds <- rsample::vfold_cv(
      df_all,
      v = ml_options$nb_folds,
      repeats = ml_options$nb_repeats,
      strata = "strata"
    )
    
  } else {
    cat("Creating stratified v-fold CV (stratified by event rate only)\n")
    
    folds <- rsample::vfold_cv(
      df_all,
      v = ml_options$nb_folds,
      repeats = ml_options$nb_repeats,
      strata = df_all[[event_col]]  # stratify by event indicator only
    )
  }
  
  # ---------------------------------------------------------------------------
  # Step 4: Define or load hyperparameter grids
  # ---------------------------------------------------------------------------
  # If no explicit hyperparameter grids are provided, automatically generate
  # default grids for supported survival models using the helper function
  # `get_default_hyperparams()`, which leverages the {dials} package.
  #
  
  # List of survival models to evaluate
  model_list <- c(
    "cox_ph_survival",              # classical Cox proportional hazards model
    "proportional_hazards_glmnet",  # penalized Cox (LASSO / elastic net)
    "survreg_flexsurv",             # parametric AFT regression
    "decision_tree_partykit",       # single survival tree
    "bag_tree_rpart",               # bagged decision trees
    #"rand_forest_partykit",        # random survival forest (ctree engine)
    "rand_forest_aorsf"             # oblique random survival forest
    #"boost_tree_mboost"            # gradient boosting for survival (optional)
  )
  
  # Create model-specific parameter grids.
  # Each element in model_grids is a list of hyperparameter value vectors.
  model_grids <- purrr::map(model_list, ~ get_default_hyperparams(.x, train_x = df_features))
  names(model_grids) <- model_list
  
  # ---------------------------------------------------------------------------
  # Step 5: Nested loop over models and hyperparameter configurations
  # ---------------------------------------------------------------------------
  # Outer loop: iterate over each model.
  # Inner loop: iterate over all hyperparameter combinations for that model.
  # Innermost loop: perform v-fold cross-validation for each configuration.
  #
  # This structure implements nested CV: tuning parameters inside cross-validation.
  #
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  
  # Export necessary global objects to all worker processes
  parallel::clusterExport(
    cl,
    varlist = c("folds", "outcome_col", "event_col", "ml_options"),
    envir = environment()
  )
  
  # ---------------------------------------------------------------------------
  # Step 6: Model training and evaluation loop
  # ---------------------------------------------------------------------------
  # For each model in the list, train and evaluate all configurations in parallel.
  # Results are aggregated into a single tibble using map_dfr().
  #
  all_results <- purrr::map_dfr(model_list, function(current_model) {
    
    cat("Training with model", current_model, "\n")
    
    # Retrieve hyperparameter grid for the current model
    hyperparams <- model_grids[[current_model]]
    
    # Fallback: if the model has no tunable parameters, create a single configuration
    if (is.null(hyperparams)) {
      param_grid <- tibble::tibble(.config_id = 1)
    } else {
      param_grid <- tidyr::expand_grid(!!!hyperparams) %>%
        dplyr::mutate(.config_id = row_number())
    }
    
    # -----------------------------------------------------------------------
    # Inner hyperparameter tuning loop
    # -----------------------------------------------------------------------
    # For each hyperparameter configuration, perform full v-fold CV evaluation.
    #
    model_results <- purrr::map_dfr(1:nrow(param_grid), function(g) {
      current_params <- param_grid[g, , drop = FALSE]
      
      # Export model-specific variables to all cluster workers
      parallel::clusterExport(
        cl,
        varlist = c("current_model", "hyperparams", "current_params", "g"),
        envir = environment()
      )
      
      # ---------------------------------------------------------------------
      # Inner-most loop: parallelized cross-validation folds
      # ---------------------------------------------------------------------
      # Train the model and compute C-index for each fold in parallel.
      #
      results <- foreach::foreach(i = seq_along(folds$splits), 
                                  .combine = dplyr::bind_rows,
                                  .packages = c("dplyr")) %dopar% {
                                    # Source the model training function in worker environments
                                    source("machine_learning_survival.R")
                                    
                                    # Split into training and validation sets
                                    split <- folds$splits[[i]]
                                    train_df <- rsample::analysis(split)
                                    test_df  <- rsample::assessment(split)
                                    
                                    # Remove auxiliary columns only in LODO mode
                                    if (isTRUE(ml_options$LODO)) {
                                      drop_cols <- c("strata", ml_options$batch_id)
                                      train_df <- train_df %>% dplyr::select(-dplyr::any_of(drop_cols))
                                      test_df  <- test_df %>% dplyr::select(-dplyr::any_of(drop_cols))
                                    }
                                    
                                    # Train model and evaluate C-index on the held-out fold
                                    trained <- machine_learning_custom(
                                      df_train = train_df,
                                      df_test = test_df,
                                      outcome_col = outcome_col,
                                      event_col = event_col,
                                      model = current_model,
                                      models_hyperparameters = if (is.null(hyperparams)) NULL else list(
                                        current_params %>% dplyr::select(-.config_id)
                                      )
                                    )
                                    
                                    # Return fold-specific results
                                    tibble::tibble(
                                      model = current_model,
                                      .config_id = g,
                                      fold = folds$id[i],
                                      c_index = trained$c_index
                                    )
                                  }
      
      # Aggregate fold results: compute median and MAD (robust dispersion)
      results %>%
        dplyr::group_by(model, .config_id) %>%
        dplyr::summarise(
          median = stats::median(c_index, na.rm = TRUE),
          mad = stats::mad(c_index, na.rm = TRUE),
          .groups = "drop"
        )
      
    })
    
    # Return summarized cross-validation results for the current model
    model_results
    
  })
  
  # Stop the parallel cluster and clean up
  parallel::stopCluster(cl)
  unregister_dopar() # Ensure no background parallel backend remains
  
  # ---------------------------------------------------------------------------
  # Step 7: Identify the best hyperparameter configuration
  # ---------------------------------------------------------------------------
  # After all configurations have been evaluated:
  #   - Identify the best-performing configuration per model
  #   - Select the top model overall based on the median C-index
  #
  best_configs_per_model <- all_results %>%
    dplyr::group_by(model) %>%
    dplyr::slice_max(median, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  # Select the best model overall based on highest median C-index
  best_model <- best_configs_per_model %>%
    dplyr::slice_max(median, n = 1, with_ties = FALSE)
  
  # Extract the best model name and configuration ID
  best_model_name <- best_model$model
  best_config_id <- best_model$.config_id
  
  # Retrieve corresponding hyperparameter values for the winning configuration
  best_hyperparams <- tidyr::expand_grid(!!!get_default_hyperparams(best_model_name)) %>%
    dplyr::slice(best_config_id)  
  
  # Visualize and summarize performance
  plot_summary <- compute_cv_CINDEX(best_configs_per_model, file_name = file_name)
  
  print(plot_summary$CINDEX_summary)
  cat("Top model:", plot_summary$Top_model, "\n")
  
  # ---------------------------------------------------------------------------
  # Step 8: Return structured results
  # ---------------------------------------------------------------------------
  # The output is a comprehensive list containing:
  #   - all_results: C-index performance for each model/configuration/fold
  #   - best_configs_per_model: top configuration per model
  #   - best_model: name of the best-performing model
  #   - best_score: median C-index of the top model
  #   - best_sd: MAD (robust SD) of that configuration
  #
  list(
    all_results = all_results,
    best_configs_per_model = best_configs_per_model,
    best_model = best_model$model,
    best_score = best_model$median,
    best_sd = best_model$mad
  )
}

# =============================================================================
# Function: compute_cv_CINDEX()
# Purpose : Summarize and visualize C-index results across survival models
#           obtained from the output of cross_validation_custom().
#           Produces a summary table and generates a bar plot showing
#           median C-index ± MAD (robust SD) for each model.
# =============================================================================

compute_cv_CINDEX = function(surv_results, file_name = NULL, return = TRUE) {
  
  # ---------------------------------------------------------------------------
  # Step 1: Prepare summarized C-index results
  # ---------------------------------------------------------------------------
  # The input `surv_results` should contain one row per model with:
  #   - model: model identifier
  #   - median: median C-index from cross-validation
  #   - mad: median absolute deviation of C-index (robust SD estimate)
  #
  # This section renames and reorders columns for readability and plotting.
  #
  res_cindex <- surv_results %>%
    dplyr::select(model, median, mad) %>%
    dplyr::rename(Median_CINDEX = median, MAD_CINDEX = mad) %>%
    dplyr::arrange(desc(Median_CINDEX))
  
  # ---------------------------------------------------------------------------
  # Step 2: Plot median C-index ± MAD per model
  # ---------------------------------------------------------------------------
  # Generates a bar plot summarizing the cross-validation performance
  # of each survival model. Bars represent the median C-index, and error bars
  # correspond to ± MAD_CINDEX (robust variability).
  #
  # The figure is saved as a PDF to the "Results/" folder, named according
  # to the provided `file_name` argument.
  #
  # The layout uses a consistent ggplot2 style (matching compute_cv_AUC()).
  #
  if (return) {
    grDevices::pdf(paste0("Results/CINDEX_CV_methods_", file_name, ".pdf"), width = 10)
    plot(
      ggplot2::ggplot(res_cindex, ggplot2::aes(x = model, y = Median_CINDEX, fill = model)) +
        ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(), width = 0.6) +
        ggplot2::geom_errorbar(
          ggplot2::aes(ymin = Median_CINDEX - MAD_CINDEX, ymax = Median_CINDEX + MAD_CINDEX),
          width = 0.2,
          position = ggplot2::position_dodge(0.6)
        ) +
        ggplot2::labs(
          title = "Performance of Survival Models",
          x = "Model",
          y = "Median C-index"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          legend.position = "none",
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        ) +
        ggplot2::scale_y_continuous(breaks = seq(0, 1, by = 0.05))
    )
    grDevices::dev.off()
  }
  
  # ---------------------------------------------------------------------------
  # Step 3: Identify top-performing model
  # ---------------------------------------------------------------------------
  # Determines the model with the highest median C-index.
  # Ties are resolved by taking the first occurrence (with_ties = FALSE).
  #
  top_model <- res_cindex %>%
    dplyr::slice_max(Median_CINDEX, n = 1, with_ties = FALSE) %>%
    dplyr::pull(model)
  
  # ---------------------------------------------------------------------------
  # Step 4: Return structured summary
  # ---------------------------------------------------------------------------
  # Returns a list containing:
  #   - CINDEX_summary: table of models and corresponding C-index statistics
  #   - Top_model: the name of the model with the best median C-index
  #
  list(
    CINDEX_summary = res_cindex,
    Top_model = top_model
  )
}


#' Plot and save survival performance of a model on test data
#'
#' This function groups individuals into risk strata (e.g., Low/Medium/High)
#' based on their predicted risk scores from a fitted survival model.
#' It then plots Kaplan–Meier survival curves for each risk group,
#' including a log-rank test for group separation and displays the
#' concordance index (C-index) of the model on the test data.
#'
#' @param df_test A data frame containing at least the following columns:
#'   \describe{
#'     \item{time}{Observed survival or follow-up time (numeric).}
#'     \item{event}{Event indicator (1 = event occurred, 0 = censored).}
#'     \item{.pred}{Predicted risk score or linear predictor from the model
#'       (higher values indicate higher risk).}
#'   }
#' @param c_index Numeric. The concordance index (C-index) computed on the test data.
#' @param n_groups Integer. Number of risk groups to stratify by (default = 3).
#'   Typically 3 groups correspond to "Low", "Medium", and "High" risk strata.
#' @param file_name Character (optional). If provided, the Kaplan–Meier plot will be
#'   saved to \code{"Results/Survival_KM_<file_name>.pdf"}.
#'
#' @details
#' Risk groups are defined by quantile-based cut points on the predicted risk scores.
#' The log-rank test is used to assess whether survival curves differ significantly
#' between risk strata. The C-index is displayed for interpretability.
#'
#' @return Invisibly returns the \code{ggsurvplot} object for further customization,
#'   and saves a PDF of the plot in the "Results/" directory if \code{file_name} is provided.
#'
#' @examples
#' \dontrun{
#' plot_survival_performance(
#'   df_test = test_data,
#'   c_index = 0.74,
#'   n_groups = 3,
#'   file_name = "cox_model_test"
#' )
#' }
#'
#' @export
plot_survival_performance <- function(df_test, c_index = NULL, n_groups = 3, file_name = NULL) {
  suppressPackageStartupMessages({
    library(survival)
    library(ggplot2)
    library(survminer)
    library(dplyr)
  })
  
  # Check required columns
  required_cols <- c("time", "event", ".pred")
  if (!all(required_cols %in% names(df_test))) {
    stop("df_test must include columns: time, event, and .pred")
  }
  
  # Group patients into risk categories based on predicted risk (.pred)
  df_test <- df_test %>%
    dplyr::mutate(
          risk_group = cut(
            .pred,
            breaks = stats::quantile(.pred, probs = seq(0, 1, length.out = n_groups + 1), na.rm = TRUE),
            include.lowest = TRUE,
            labels = paste(c("Low", "Medium", "High")[1:n_groups], "risk")
          )
    )
  
  # Fit Kaplan–Meier survival curves directly from formula
  fit_km <- survival::survfit(survival::Surv(time, event) ~ risk_group, data = df_test)
  
  # Log-rank test
  logrank <- survival::survdiff(survival::Surv(time, event) ~ risk_group, data = df_test)
  p_val <- 1 - stats::pchisq(logrank$chisq, df = length(logrank$n) - 1)
  
  subtitle_text <- paste0("C-index: ", round(c_index, 3),
                          " | Log-rank p = ", format.pval(p_val, digits = 3, eps = .001))
  
  # Plot survival curves
  plt <- survminer::ggsurvplot(
    fit_km,
    data = df_test,
    risk.table = TRUE,
    pval = FALSE,
    ggtheme = ggplot2::theme_minimal(),
    palette = c("#1B9E77", "#7570B3", "#D95F02")[1:n_groups],
    legend.title = "Risk Group",
    legend.labs = levels(df_test$risk_group),
    title = paste0("Test-set Survival Performance_", file_name),
    subtitle = subtitle_text
  )
  
  grDevices::pdf(paste0("Results/Survival_KM_", file_name, ".pdf"), width = 8, height = 6)
  print(plt) 
  grDevices::dev.off()
  
}
