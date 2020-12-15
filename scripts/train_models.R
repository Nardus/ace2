## Train models to predict susceptibility to sarbecovirus infection

# ---- Input args ---------------------------------------------------------------------------------
# Takes one argument, specifying the dataset
dataset_name = commandArgs(trailingOnly = TRUE)

if (length(dataset_name) != 1 | !(dataset_name %in% c("infection", "shedding", "transmission")))
  stop(paste("A single argument specifying the dataset to predict is required", 
             "(either 'infection', 'shedding', or 'transmission')"), 
       call. = FALSE)

metadata_path <- sprintf("data/calculated/cleaned_%s_data.rds", dataset_name)
response <- switch(dataset_name,
                   "infection" = "infected",
                   "shedding" = "shedding",
                   "transmission" = "transmission")


# ---- Setup --------------------------------------------------------------------------------------
library(tibble)
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



# Hyper-parameter search:
N_HYPER_PARAMS <- 100      # Number of hyper-parameter combinations to try

# Tuning parameter values to try:
# - N = N_HYPER_PARAMS random combinations of these will be tested 
TUNING_PARAMETERS <- list(eta = c(0.001, 0.005, seq(0.01, 0.2, by = 0.02)),
                          max_depth = seq(6, 15, by = 1),
                          subsample = seq(0.6, 1.0, by = 0.1),
                          colsample_bytree = seq(0.5, 1.0, by = 0.1),
                          nrounds = seq(50, 250, by = 10), # Keeping this somewhat low to prevent over-fitting
                          min_child_weight = seq(0, 10, by = 2),
                          gamma = seq(0, 7, by = 0.5))


set.seed(41210412)
registerDoParallel(8)


# ---- Data ---------------------------------------------------------------------------------------
# Metadata
metadata <- read_rds(metadata_path) %>% 
  rename(label = all_of(response))

# Feature data
pairwise_dist_data <- read_rds("data/calculated/features_pairwise_dists.rds")
dist_to_humans <- read_rds("data/calculated/features_dist_to_humans.rds")
variable_sites <- read_rds("data/calculated/features_variable_sites.rds")

# Combine
final_data <- metadata %>% 
  left_join(dist_to_humans, by = "ace2_accession") %>% 
  left_join(variable_sites, by = "ace2_accession")


stopifnot(nrow(final_data) == n_distinct(metadata$species))

# Remove features which do not vary much in the current dataset:
# - happens among "variable site" features in particular
# - these columns will often be zero-variance once data are split, causing an error in caret
remove_cols <- nearZeroVar(final_data)
final_data <- final_data[, -remove_cols]


# ---- Training -----------------------------------------------------------------------------------

# Tuning using 3 replicates of 5-fold CV in each iteration
train_setup <- trainControl(method = "LOOCV",
                            classProbs = TRUE,
                            search = "random",
                            savePredictions = "final")
  
# Calculate additional (training set-specific) features
#  - These depend on the particular test set, but are correct by default for leave-one-out CV, 
#    since the current virus is not included when summarising its neighbours
closest_positive <- get_dist_to_closest_positive(pairwise_dist_data, metadata)
consensus_dists <- get_consensus_dist(variable_sites, variable_sites, metadata)
  
final_data <- final_data %>% 
  left_join(closest_positive, by = "ace2_accession") %>% 
  left_join(consensus_dists, by = "ace2_accession") %>% 
  mutate(across(where(is.character), as.factor))
  
# Prepare data for caret
final_data <- final_data %>% 
  mutate(label = factor(.data$label, levels = c("True", "False"))) %>%  # Caret's twoClassSummary treats level 1 as "TRUE"
  as.data.frame()

train_data <- final_data %>% 
  select(-.data$species, -.data$ace2_accession, -.data$evidence_level) # Remove non-feature columns
  
# Train
parameter_combos <- lapply(TUNING_PARAMETERS, sample, size = N_HYPER_PARAMS, replace = TRUE) %>% 
  bind_rows() %>% 
  as.data.frame()
  
trained_model <- train(label ~ .,
                       data = train_data,
                       method = "xgbTree",
                       metric = "Accuracy",
                       trControl = train_setup,
                       tuneGrid = parameter_combos,
                       na.action = na.pass,
                       nthread = 1)

# Correct predictions
# - Caret uses a hard-coded cutoff of 0.5, but that's far from optimal
cutoff <- sum(final_data$label == "True")/nrow(final_data)

predictions <- trained_model$pred %>% 
  mutate(pred = if_else(.data$True > cutoff, "True", "False"))


# ---- Output -------------------------------------------------------------------------------------
dir.create(sprintf("output/%s/", dataset_name), 
           recursive = TRUE)

# Predictions
predictions <- final_data %>% 
  rowid_to_column("rowIndex") %>% 
  full_join(predictions, by = "rowIndex") %>% 
  select(.data$species, .data$ace2_accession, .data$label, .data$evidence_level,
         prediction = .data$pred,
         prob = .data$True)

stopifnot(nrow(predictions) == nrow(final_data))

write_rds(predictions, sprintf("output/%s/predictions.rds", dataset_name))


# Model
training_results <- list(finalModel = trained_model,
                         trainingData = final_data)

write_rds(training_results, sprintf("output/%s/training_results.rds", dataset_name), 
          compress = "gz", compression = 9)

# Binary version of model, compatible with future xgboost releases:
xgb.save(trained_model$finalModel, sprintf("output/%s/xgboost_model.model", dataset_name))