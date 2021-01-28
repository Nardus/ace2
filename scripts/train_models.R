## Train models to predict susceptibility to sarbecovirus infection

# ---- Input args ---------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(argparse)
})

parser <- ArgumentParser(description = "Train ACE2 models using varying subsets of data and features")

parser$add_argument("response_var", type = "character", 
                    choices = c("infection", "shedding", "transmission"),
                    help = "response variable to train on.")

parser$add_argument("output_path", type = "character",
                    help = "location for output files. Missing folders will be created.")


features_group <- parser$add_argument_group("Feature sets")
features_group$add_argument("--aa_categorical", action = "store_const", const = TRUE, default = FALSE,
                            help = "include features representing variable amino acids directly as categorical ('N/T/S/...').")

features_group$add_argument("--aa_distance", action = "store_const", const = TRUE, default = FALSE,
                            help = "include features representing variable amino acids as a distance to the closest positive.")

features_group$add_argument("--distance_to_humans", action = "store_const", const = TRUE, default = FALSE,
                            help = "include a feature measuring overall amino acid distance to human ACE2")

features_group$add_argument("--binding_affinity", action = "store_const", const = TRUE, default = FALSE,
                            help = "include features measuring binding affinity to the SARS-CoV-2 spike protein")


data_group <- parser$add_argument_group("Dataset options")
data_group$add_argument("--max_evidence_level", type = "integer", choices = 1L:4L, default = 4L,
                        help = paste("data to include. Higher numbers include less robust evidence:",
                                     "1 = natural infection, 2 = experimental infection, 3 = cell culture",
                                     "4 = cells modified to express ACE2 (default = 4)."))


other_opts_group <- parser$add_argument_group("Other options")
other_opts_group$add_argument("--random_seed", type = "integer", default = trunc(runif(1, max = 1e5)),
                              help = "random seed to use (default: a random integer between 0 and 1e5)")

other_opts_group$add_argument("--n_threads", type = "integer", default = 8,
                              help = "number of parallel threads allowed (default: 8)")


## Check input
INPUT <- parser$parse_args()

if (!any(INPUT$aa_categorical, INPUT$aa_distance, INPUT$distance_to_humans, INPUT$binding_affinity))
  stop("No features selected. Run train_models.R --help for available feature sets")


# Response variable:
metadata_path <- sprintf("data/calculated/cleaned_%s_data.rds", INPUT$response_var)
response <- switch(INPUT$response_var,
                   "infection" = "infected",
                   "shedding" = "shedding",
                   "transmission" = "transmission")

# Features
feature_prefixes <- c()

if (INPUT$aa_categorical)
  feature_prefixes <- c(feature_prefixes, "variable_site_")

if (INPUT$aa_distance)
  feature_prefixes <- c(feature_prefixes, "dist_variable_site_")

if (INPUT$distance_to_humans)
  feature_prefixes <- c(feature_prefixes, "distance_to_humans")

if (INPUT$binding_affinity)
  stop("Not yet implemented") # TODO


# ---- Setup --------------------------------------------------------------------------------------
suppressPackageStartupMessages({
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
})


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


set.seed(INPUT$random_seed)
registerDoParallel(INPUT$n_threads)


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
  left_join(consensus_dists, by = "ace2_accession")

# Final dataset: subset to match input options
final_data <- final_data %>% 
  mutate(across(where(is.character), as.factor)) %>% 
  select(.data$species, .data$ace2_accession, .data$evidence_level, .data$label,
         starts_with(feature_prefixes)) %>% 
  filter(.data$evidence_level <= INPUT$max_evidence_level)
  
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
# TODO: may need to be optimized
cutoff <- 0.5 

predictions <- trained_model$pred %>% 
  mutate(pred = if_else(.data$True > cutoff, "True", "False"))


# ---- Output -------------------------------------------------------------------------------------
dir.create(INPUT$output_path, recursive = TRUE)

# Predictions
predictions <- final_data %>% 
  rowid_to_column("rowIndex") %>% 
  full_join(predictions, by = "rowIndex") %>% 
  select(.data$species, .data$ace2_accession, .data$label, .data$evidence_level,
         prediction = .data$pred,
         prob = .data$True)

stopifnot(nrow(predictions) == nrow(final_data))

write_rds(predictions, file.path(INPUT$output_path, "predictions.rds"))


# Model
training_results <- list(finalModel = trained_model,
                         trainingData = final_data)

write_rds(training_results, file.path(INPUT$output_path, "training_results.rds"), 
          compress = "gz", compression = 9)

# Binary version of model, compatible with future xgboost releases:
xgb.save(trained_model$finalModel, file.path(INPUT$output_path, "xgboost_model.model"))