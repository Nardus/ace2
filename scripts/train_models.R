## Train models to predict susceptibility to sarbecovirus infection

# ---- Input args ---------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(argparse)
})

parser <- ArgumentParser(description = "Train ACE2 models using varying subsets of data and features")

parser$add_argument("response_var", type = "character", 
                    choices = c("infection", "shedding"),
                    help = "response variable to train on.")

parser$add_argument("output_path", type = "character",
                    help = "location for output files. Missing folders will be created.")


features_group <- parser$add_argument_group("Feature sets")
features_group$add_argument("--aa_categorical", action = "store_const", const = TRUE, default = FALSE,
                            help = "include features representing variable amino acids directly as categorical ('N/T/S/...').")

features_group$add_argument("--aa_distance", action = "store_const", const = TRUE, default = FALSE,
                            help = "include features representing variable amino acids as a distance to the closest positive.")

features_group$add_argument("--aa_properties", action = "store_const", const = TRUE, default = FALSE,
                            help = "include features representing variable amino acids by their physico-chemical properties.")

features_group$add_argument("--distance_to_humans", action = "store_const", const = TRUE, default = FALSE,
                            help = "include a feature measuring overall amino acid distance to human ACE2")

features_group$add_argument("--distance_to_positive", action = "store_const", const = TRUE, default = FALSE,
                            help = paste("include features measuring overall amino acid distance to",
                                         "the closest positive (at each evidence level)"))

features_group$add_argument("--binding_affinity", action = "store_const", const = TRUE, default = FALSE,
                            help = "include features measuring binding affinity to the SARS-CoV-2 spike protein")


data_group <- parser$add_argument_group("Dataset options")
data_group$add_argument("--evidence_min", type = "integer", choices = 1L:4L, default = 1L,
                        help = paste("lowest evidence level to include. Higher numbers include less robust evidence:",
                                     "1 = natural infection, 2 = experimental infection, 3 = cell culture",
                                     "4 = cells modified to express ACE2 (default = 1)."))

data_group$add_argument("--evidence_max", type = "integer", choices = 1L:4L, default = 4L,
                        help = "maximum evidence level to include (default = 4).")


other_opts_group <- parser$add_argument_group("Other options")
other_opts_group$add_argument("--random_seed", type = "integer", default = trunc(runif(1, max = 1e5)),
                              help = "random seed to use (default: a random integer between 0 and 1e5)")

other_opts_group$add_argument("--n_threads", type = "integer", default = 16,
                              help = "number of parallel threads allowed (default: 16)")


## Check input
INPUT <- parser$parse_args()

if (!any(INPUT$aa_categorical, INPUT$aa_distance, INPUT$aa_properties, INPUT$distance_to_humans, 
         INPUT$distance_to_positive, INPUT$binding_affinity))
  stop("No features selected. Run train_models.R --help for available feature sets")


# Response variable:
metadata_path <- sprintf("data/calculated/cleaned_%s_data.rds", INPUT$response_var)
response <- switch(INPUT$response_var,
                   "infection" = "infected",
                   "shedding" = "shedding")

# Features
feature_prefixes <- c()

if (INPUT$aa_categorical)
  feature_prefixes <- c(feature_prefixes, "variable_site_")

if (INPUT$aa_distance)
  feature_prefixes <- c(feature_prefixes, "dist_variable_site_")

if (INPUT$aa_properties)
  feature_prefixes <- c(feature_prefixes, "property_")

if (INPUT$distance_to_humans)
  feature_prefixes <- c(feature_prefixes, "distance_to_humans", "distance_to_rhinolophid")

if (INPUT$distance_to_positive)
  feature_prefixes <- c(feature_prefixes, "closest_positive_")

if (INPUT$binding_affinity)
  feature_prefixes <- c(feature_prefixes, "haddock_score")


# ---- Setup --------------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  
  library(xgboost)
  library(recipes)
  library(themis)
  library(rsample)
  library(workflows)
  library(dials)
  library(tune)
  library(yardstick)
  library(parsnip)
  
  library(purrr)
  library(future)
  library(doFuture)
  library(furrr)
  
  source("scripts/utils/aa_distance_utils.R")
  source("scripts/utils/feature_calc_utils.R")
  source("scripts/utils/feature_recipe_utils.R")
})


# Constants:
N_HYPER_PARAMS <- 100      # Number of hyper-parameter combinations to try

set.seed(INPUT$random_seed)
plan(strategy = multisession(workers = INPUT$n_threads))
registerDoFuture()


# ---- Data ---------------------------------------------------------------------------------------
# Metadata
metadata <- read_rds(metadata_path) %>% 
  rename(label = all_of(response))

# Filter evidence levels
# - If filtering, also need to adjust the response variable to match the best (lowest) evidence
#   level available
# - Note: the which.min() step below:
#     - accurately handles NA's (caused by empty evidence string)
#     - returns 1 (True) when there is a tie - as in "prepare_data.R", conflicting evidence is 
#       ignored because we are simply asking "has it _ever_ been reported to be infected/shedding/etc?"
if (!(INPUT$evidence_min == 1 & INPUT$evidence_max == 4)) {
  removed_levels <- c(seq(0, (INPUT$evidence_min - 1)),
                      seq((INPUT$evidence_max + 1), 5))
  removed_levels <- paste(removed_levels, collapse = "|")
  
  metadata <- metadata %>% 
    mutate(all_evidence_true = str_remove_all(.data$all_evidence_true, removed_levels),
           all_evidence_false = str_remove_all(.data$all_evidence_false, removed_levels),
           all_evidence_true = str_remove(.data$all_evidence_true, "^,*"),
           all_evidence_false = str_remove(.data$all_evidence_false, "^,*")) %>% 
    filter(!(.data$all_evidence_true == "" & .data$all_evidence_false == "")) %>% 
    
    group_by(.data$species) %>% 
    mutate(evidence_true = substring(.data$all_evidence_true, 1, 1),  # Levels sorted, so this gets the lowest remaining evidence level
           evidence_false = substring(.data$all_evidence_false, 1, 1),
           new_response = which.min(c(.data$evidence_true, .data$evidence_false)),
           label = c("True", "False")[.data$new_response],
           evidence_level = c(.data$evidence_true, .data$evidence_false)[.data$new_response],
           evidence_level = as.integer(.data$evidence_level)) %>% 
    ungroup() %>% 
    select(-.data$evidence_true, -.data$evidence_false, -.data$new_response)
}


# Feature data
pairwise_dist_data <- read_rds("data/calculated/features_pairwise_dists.rds")
dist_to_humans <- read_rds("data/calculated/features_dist_to_humans.rds")
variable_sites <- read_rds("data/calculated/features_variable_sites.rds")
site_properties <- read_rds("data/calculated/features_site_properties.rds")
haddock_scores <- read_rds("data/calculated/features_haddock_scores.rds")

# Combine
final_data <- metadata %>% 
  left_join(dist_to_humans, by = "ace2_accession") %>% 
  left_join(variable_sites, by = "ace2_accession") %>% 
  left_join(site_properties, by = "ace2_accession") %>% 
  left_join(haddock_scores, by = "species")


stopifnot(nrow(final_data) == n_distinct(metadata$species))

# Final processing
final_data <- final_data %>% 
  mutate(label = factor(.data$label, levels = c("True", "False"))) %>% 
  mutate_if(is.character, as.factor)


# ---- Pre-processing -----------------------------------------------------------------------------
# Specify base recipe
id_columns <- c("species", "all_evidence_true", "all_evidence_false", 
                "evidence_level", "ace2_accession")

preprocessing_recipe <- 
  recipe(label ~ ., data = final_data) %>% 
  update_role(all_of(id_columns), new_role = "ID") %>% 
  step_distance_features(everything(),
                         options = list(all_pairwise_dists = pairwise_dist_data, 
                                        all_variable_sites = variable_sites, 
                                        metadata = metadata),
                         role = "predictor")


# Subset features to match input options
preprocessing_recipe <- 
  preprocessing_recipe %>% 
  step_rm(-has_role("ID"), -all_outcomes(), -starts_with(feature_prefixes))


# Convert character values and remove invariant columns:
preprocessing_recipe <-
  preprocessing_recipe %>% 
  step_zv(all_predictors(), -has_role("ID"), -all_outcomes()) %>% 
  step_unknown(all_nominal(), -has_role("ID"), -all_outcomes()) %>% 
  step_other(all_nominal(), -has_role("ID"), -all_outcomes(), threshold = 0.1) %>%  # Rare amino acids (frequency < 10%) collapsed 
  step_dummy(all_nominal(), -has_role("ID"), -all_outcomes(), one_hot = TRUE)


# Objects/settings needed to evaluate this recipe in parallel:
recipe_opts <- furrr_options(globals = c(recipe_globals, "feature_prefixes"),
                             packages = c("dplyr", "tidyr", "themis", "tune", "yardstick", "rsample"),
                             seed = TRUE)


# ---- Model setup --------------------------------------------------------------------------------
xgboost_model <- 
  boost_tree(mode = "classification",
             learn_rate = tune(),
             tree_depth = tune(),
             sample_size = tune(),
             mtry = tune(),
             trees = tune(),
             min_n = tune(),
             loss_reduction = tune()) %>%
  set_engine("xgboost", eval_metric = "logloss")


tuning_parameters <- parameters(learn_rate(range = c(-5, -0.7)),          # eta
                                tree_depth(),                             # max_depth
                                sample_prop(range = c(0.6, 1.0)),         # subsample (sample_size(), but as proportion)
                                mtry(),                                   # colsample_bynode, limits set based on number of features below
                                trees(range = c(1L, 250L)),               # nrounds, keeping this low to prevent over-fitting (more boosting rounds also increases number of features used)
                                min_n(range = c(5L, 10L)),                # min_child_weight, somewhat high - we don't want features used to predict just one or two cases (so at least 5)
                                loss_reduction(range = c(-10, 10)))       # gamma, allow very high values, which would make model extremely conservative


# ---- Training -----------------------------------------------------------------------------------
# Set up cross-validation
cv_folds <- nested_cv(data = final_data,
                      outside = loo_cv(), 
                      inside = validation_split(strata = label, prop = 0.75))

# Record training data by preparing the pre-processing_recipe for each split 
#  - Only re-calculating features for outer folds, since this takes a while
cv_folds$recipes <- future_map(cv_folds$splits, prepper, recipe = preprocessing_recipe,
                               .options = recipe_opts)

# Apply to tuning to the inner splits of each (outer) fold:
tune_inner <- function(inner_splits, fold_recipe, model, param_grid) {
  test_cutoff <- function(split_result, cutoff) {
    # split_result should represent tuning results from a single param combination and data split
    split_result %>% 
      mutate(new_prediction = if_else(.data$.pred_True > cutoff, "True", "False"),
             new_prediction = factor(new_prediction, levels = c("True", "False"))) %>% 
      group_by(.data$label, .add = TRUE) %>% 
      summarise(class_acc = accuracy_vec(truth = .data$label, estimate = .data$new_prediction),
                .groups = "drop_last") %>% 
      summarise(.estimate = sum(.data$class_acc)/2,
                .groups = "keep") %>% 
      mutate(.metric = "balanced_accuracy",
             .estimator = "binary",
             cutoff = cutoff)
  }
  
  find_best_cutoff <- function(split_result, cutoffs = seq(0.15, 0.95, by = 0.01)) {
    sapply(cutoffs, test_cutoff, split_result = split_result, simplify = FALSE) %>% 
      bind_rows() %>% 
      slice_max(.data$.estimate, with_ties = TRUE) %>% 
      slice_sample(n = 1) # Choose randomly in case of ties
  }
  
  tuning_results <- tune_grid(object = model,
                              preprocessor = fold_recipe,
                              resamples = inner_splits,
                              grid = param_grid,
                              metrics = metric_set(roc_auc), # Not used, but triggers return of quantitative predictions for cutoff optimisation below
                              control = control_grid(save_pred = TRUE,
                                                     allow_par = FALSE)) # Will parallelise the outer loop instead (so we can specify recipe_opts)
  
  # Test cutoffs on the same training/validation splits:
  for (i in 1:nrow(tuning_results)) {  # Loop over splits
    cutoff_metrics <- tuning_results[i, ] %>% 
      unnest(.data$.predictions) %>% 
      group_by(.data$.config) %>% 
      find_best_cutoff()
    
    tuning_results$.metrics[[i]] <- tuning_results$.metrics[[i]] %>% 
      select(-starts_with("."), .config) %>% 
      full_join(cutoff_metrics, by = ".config")
  } 
  
  tuning_results %>% 
    select_best()
}

parameter_combos <- tuning_parameters %>% 
  finalize(final_data) %>% 
  grid_max_entropy(size = N_HYPER_PARAMS)

best_params <- future_map2(cv_folds$inner_resamples, cv_folds$recipes, tune_inner,
                           model = xgboost_model, 
                           param_grid = parameter_combos,
                           .options = recipe_opts)


# Fit final modelling workflow on all training data from each outer fold
fit_final <- function(outer_split, fold_recipe, hyperparams, model) {
  tuned_model <- finalize_model(model, hyperparams)
  fold_training_data <- analysis(outer_split)
  
  workflow() %>% 
    add_recipe(fold_recipe) %>% 
    add_model(tuned_model) %>% 
    fit(fold_training_data)
}

final_models <- future_pmap(.l = list(outer_split = cv_folds$splits, 
                                      fold_recipe = cv_folds$recipes,
                                      hyperparams = best_params), 
                            .f = fit_final,
                            model = xgboost_model,
                            .options = recipe_opts)

# Get predictions for the holdout in each outer fold
predict_final <- function(outer_split, fold_recipe, fold_id, trained_workflow, params) {
  cutoff <- params$cutoff
  fold_test_data <- assessment(outer_split)
  
  predict(trained_workflow, fold_test_data, type = "prob") %>% 
    bind_cols(fold_test_data) %>% 
    mutate(cv_fold = fold_id,
           cutoff = cutoff,
           prediction = if_else(.data$.pred_True > cutoff, "True", "False"),
           prediction = factor(.data$prediction, levels = c("True", "False"))) %>% 
    select(species, label, cv_fold, prediction, cutoff,
           p_true = .pred_True,
           p_false = .pred_False)
}

predictions <- future_pmap_dfr(.l = list(outer_split = cv_folds$splits, 
                                         fold_recipe = cv_folds$recipes,
                                         fold_id = cv_folds$id,
                                         trained_workflow = final_models,
                                         params = best_params), 
                               .f = predict_final,
                               .options = recipe_opts)

# Fit a final model on all data:
# - Steps above gave an estimate of how well the *entire fitting procedure* works
# - Use best params found across all outer replicates (and repeats) to fit final model
outer_performance <- predictions %>% 
  group_by(.data$cv_fold) %>% 
  summarise(accuracy = accuracy_vec(truth = .data$label, estimate = .data$prediction),
            .groups = "drop")

best_id <- outer_performance %>% 
  slice_max(accuracy) %>% 
  slice_sample(n = 1) # If there are ties, choose randomly (many ties expected with LOOCV)

best_id <- which(cv_folds$id == best_id$cv_fold)

final_params <- best_params[[best_id]]

final_workflow <- workflow() %>% 
  add_recipe(preprocessing_recipe) %>% 
  add_model(xgboost_model) %>%
  finalize_workflow(final_params) %>% 
  fit(final_data)


# ---- Feature importance --------------------------------------------------------------------------
shap_values <- predict(final_workflow, final_data, type = "raw", opts = list(predcontrib = TRUE))

feature_importance <- shap_values %>% 
  data.frame() %>% 
  pivot_longer(everything(), names_to = "feature", values_to = "shap") %>% 
  group_by(feature) %>% 
  summarise(importance = mean(abs(shap))) %>% 
  filter(.data$feature != "BIAS") %>% 
  arrange(-importance)


# ---- Output -------------------------------------------------------------------------------------
dir.create(INPUT$output_path, recursive = TRUE)

# Predictions
write_rds(predictions, file.path(INPUT$output_path, "predictions.rds"))

# Models associated with these predictions
names(final_models) <- cv_folds$id
write_rds(final_models, file.path(INPUT$output_path, "cv_models.rds"),
          compress = "gz", compression = 9)

# Final model workflow
write_rds(final_workflow, file.path(INPUT$output_path, "trained_model_workflow.rds"), 
          compress = "gz", compression = 9)

# Feature usage in final workflow
write_rds(shap_values, file.path(INPUT$output_path, "shap_values.rds"),
          compress = "gz", compression = 9)
write_rds(feature_importance, file.path(INPUT$output_path, "feature_importance.rds"))
