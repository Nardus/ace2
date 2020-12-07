## Train models to predict susceptibility to sarbecovirus infection

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(xgboost)
library(caret)
library(parallel)
library(doParallel)

source("scripts/utils/aa_distance_utils.R")
source("scripts/utils/feature_calc_utils.R")
source("scripts/utils/training_utils.R")


# Training iterations (repeated test/train splits)
#  - Doing 10 rounds of training with all features, then keeping only the most important features,
#    which are used for a final 100 iterations. Only this final set gets saved.
N_SELECTION_ROUNDS <- 10
N_ITERATIONS <- 100
MAX_FEATURES <- 50


# Hyper-parameter search (performed for each iteration):
CV_K <- 5                  # Number of folds for k-fold cross-validation
CV_REPS <- 5               # Number of cross-validation rounds (repeated cross-validation)
N_HYPER_PARAMS <- 100      # Number of hyper-parameter combinations to try in each iteration

# Tuning parameter values to try:
# - N = N_HYPER_PARAMS random combinations of these will be tested 
TUNING_PARAMETERS <- list(eta = c(0.001, 0.005, seq(0.01, 0.2, by = 0.02)),
                          max_depth = seq(6, 15, by = 1),
                          subsample = seq(0.6, 1.0, by = 0.1),
                          colsample_bytree = seq(0.5, 1.0, by = 0.1),
                          nrounds = seq(50, 250, by = 10), # Keeping this somewhat low to prevent over-fitting
                          min_child_weight = c(5, 6, 8, 10),
                          gamma = seq(0, 7, by = 0.5))


set.seed(41210412)
registerDoParallel(8)


# ---- Data ---------------------------------------------------------------------------------------
# Metadata
infection_data <- read_rds("data/calculated/cleaned_infection_data.rds")

# Feature data
pairwise_dist_data <- read_rds("data/calculated/features_pairwise_dists.rds")
dist_to_humans <- read_rds("data/calculated/features_dist_to_humans.rds")
variable_sites <- read_rds("data/calculated/features_variable_sites.rds")

# Combine
final_data <- infection_data %>% 
  left_join(dist_to_humans, by = "ace2_accession") %>% 
  left_join(variable_sites, by = "ace2_accession")


stopifnot(nrow(final_data) == n_distinct(infection_data$species))


# Output directories:
dir.create("output/infection/xgboost_models/", recursive = TRUE)


# ---- Training -----------------------------------------------------------------------------------

# Tuning using 3 replicates of 5-fold CV in each iteration
train_setup <- trainControl(method = "repeatedcv",
                            number = CV_K,
                            repeats = CV_REPS,
                            classProbs = TRUE,
                            summaryFunction = twoClassSummary,
                            search = "random")


feature_usage <- tibble()
predictions <- data.frame()
training_results <- list()

for (iteration in 1:(N_SELECTION_ROUNDS + N_ITERATIONS)) {
  
  if (iteration <= N_SELECTION_ROUNDS) {
    message(sprintf("Feature selection: %i of %i\r", iteration, N_SELECTION_ROUNDS))
  } else {
    message(sprintf("Training: %i of %i\r", iteration - N_SELECTION_ROUNDS, N_ITERATIONS))
  }
  
  # Train/test split
  # - Some species point to the same sequence: all references to a given accession number should be
  #   in the same partition to avoid a data leak
  train_index <- createDataPartition(final_data$infected, p = 0.7, times = 1)[[1]]
  train_accessions <- unique(final_data$ace2_accession[train_index])
  
  
  # Calculate additional (training set-specific) features
  closest_infectable <- pairwise_dist_data %>% 
    filter(.data$other_seq %in% train_accessions) %>% 
    get_dist_to_closest_infectable(infection_data)
  
  consensus_dists <- variable_sites %>% 
    filter(.data$ace2_accession %in% train_accessions) %>% 
    get_consensus_dist(variable_sites, infection_data)
  
  iteration_data <- final_data %>% 
    left_join(closest_infectable, by = "ace2_accession") %>% 
    left_join(consensus_dists, by = "ace2_accession") %>% 
    mutate(across(where(is.character), as.factor))
  
  
  # Feature selection
  if (iteration == (N_SELECTION_ROUNDS + 1)) {
    # First real training round, so define which features to keep
    message("\n")
    
    final_features <- feature_usage %>% 
      group_by(.data$feature) %>% 
      summarise(mean_importance = mean(.data$Overall), .groups = "drop") %>% 
      filter(.data$mean_importance > 0) %>% 
      top_n(n = MAX_FEATURES, wt = .data$mean_importance) %>% 
      mutate(feature = if_else(str_starts(.data$feature, "variable_site_"), # These columns dummy-coded in model,
                               str_remove(.data$feature, "[A-Z]$"),        # so remove amino acid at end
                               .data$feature)) %>% 
      pull(.data$feature) %>% 
      unique()
    
  } 
  
  if (iteration > N_SELECTION_ROUNDS) {
    # A training round, so reduce number of features
    iteration_data <- iteration_data %>% 
      select(.data$species, .data$ace2_accession, .data$evidence_level, .data$infected,
             all_of(final_features))
  }
  
  
  # Split and prepare data
  training_data <- iteration_data %>% 
    filter(.data$ace2_accession %in% train_accessions) %>% 
    as.data.frame() %>% 
    mutate(infected = if_else(.data$infected, "True", "False"),
           infected = factor(.data$infected, levels = c("True", "False"))) # Caret's twoClassSummary treats level 1 as "TRUE"
  
  test_data <- iteration_data %>% 
    filter(!.data$ace2_accession %in% train_accessions) %>%
    as.data.frame() %>% 
    mutate(infected = if_else(.data$infected, "True", "False"),
           infected = factor(.data$infected, levels = c("True", "False")))
  
  
  # Train
  parameter_combos <- lapply(TUNING_PARAMETERS, sample, size = N_HYPER_PARAMS, replace = TRUE) %>% 
    bind_rows() %>% 
    as.data.frame()
  
  train_data_fo <- training_data %>% 
    select(-.data$species, -.data$ace2_accession, -.data$evidence_level) # Remove non-feature columns
  
  trained_model <- train(infected ~ ., 
                         data = train_data_fo,
                         method = "xgbTree",
                         metric = "ROC",
                         trControl = train_setup,
                         tuneGrid = parameter_combos,
                         na.action = na.pass,
                         nthread = 1)
  
  
  # Record results
  if (iteration <= N_SELECTION_ROUNDS) {
    # Pre-selection: record feature usage only
    feature_usage <- rbind(feature_usage,
                           as_tibble(varImp(trained_model)$imp, rownames = "feature"))
    
  } else {
    # Training rounds: record predictions and models
    iteration <- iteration - N_SELECTION_ROUNDS
    
    # Predict
    train_preds <- record_predictions(trained_model, training_data, 
                                      data_name = "train", 
                                      iteration = iteration,
                                      data_cols = c("species", "ace2_accession", "infected", "evidence_level"))
    
    test_preds <- record_predictions(trained_model, test_data, 
                                     data_name = "test", 
                                     iteration = iteration,
                                     data_cols = c("species", "ace2_accession", "infected", "evidence_level"))
    
    predictions <- rbind(predictions, train_preds, test_preds)
    
    
    # Record this iteration's data and model
    training_results[[as.character(iteration)]] <- list(iteration = iteration, 
                                                        finalModel = trained_model$finalModel,
                                                        trainingData = training_data,
                                                        test_data = test_data)
    
    # Save binary version of the current model, compatible with future xgboost releases:
    xgb.save(trained_model$finalModel, sprintf("output/infection/xgboost_models/%i.model", iteration))
  }
}


message("\nDone\n")


# ---- Output -------------------------------------------------------------------------------------
write_rds(predictions, "output/infection/predictions.rds")

write_rds(training_results, "output/infection/training_results.rds", 
          compress = "gz", compression = 9)
