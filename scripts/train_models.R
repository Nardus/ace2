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

features_group$add_argument("--distance_to_humans", action = "store_const", const = TRUE, default = FALSE,
                            help = "include a feature measuring overall amino acid distance to human ACE2")

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
other_opts_group$add_argument("--select_features", type = "integer",
                              help = paste("number of features to retain in model. If not specified, all",
                                           "features will be used. If specified, --feature_importance",
                                           "is also required."))

other_opts_group$add_argument("--feature_importance", type = "character",
                              help = paste("location of a variable importance table to be used for",
                                           "feature selection prior to training. Ignored if", 
                                           "--select_features is not specified"))

other_opts_group$add_argument("--random_seed", type = "integer", default = trunc(runif(1, max = 1e5)),
                              help = "random seed to use (default: a random integer between 0 and 1e5)")

other_opts_group$add_argument("--n_threads", type = "integer", default = 16,
                              help = "number of parallel threads allowed (default: 16)")


## Check input
INPUT <- parser$parse_args()

if (!any(INPUT$aa_categorical, INPUT$aa_distance, INPUT$distance_to_humans, INPUT$binding_affinity))
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

if (INPUT$distance_to_humans)
  feature_prefixes <- c(feature_prefixes, "distance_to_humans")

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
  library(caret)
  library(parallel)
  library(doParallel)
  
  source("scripts/utils/aa_distance_utils.R")
  source("scripts/utils/feature_calc_utils.R")
  source("scripts/utils/training_utils.R")
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
haddock_scores <- read_rds("data/calculated/features_haddock_scores.rds")

# Combine
final_data <- metadata %>% 
  left_join(dist_to_humans, by = "ace2_accession") %>% 
  left_join(variable_sites, by = "ace2_accession") %>% 
  left_join(haddock_scores, by = "species")


stopifnot(nrow(final_data) == n_distinct(metadata$species))

# Remove features which do not vary much in the current dataset:
# - happens among "variable site" features in particular
# - these columns will often be zero-variance once data are split, causing an error in caret
remove_cols <- nearZeroVar(final_data, names = TRUE)
remove_cols <- remove_cols[!remove_cols %in% c("evidence_level", "all_evidence_true", "all_evidence_false")]
final_data <- final_data %>% 
  select(-all_of(remove_cols))


# ---- Feature selection data ---------------------------------------------------------------------
if (!is.null(INPUT$select_features)) {
  feature_importance <- readRDS(INPUT$feature_importance)
  
  if (nrow(feature_importance) <= INPUT$select_features)
    stop("Feature importance list contains less than the requested number of features")
  
  final_features <- feature_importance %>% 
    slice_max(n = INPUT$select_features, order_by = .data$importance, with_ties = FALSE) %>% 
    filter(!.data$feature %in% remove_cols) %>% 
    pull(.data$feature)
}


# ---- Training -----------------------------------------------------------------------------------

# Tuning using 3 replicates of 5-fold CV in each iteration
train_setup <- trainControl(method = "LOOCV",
                            classProbs = TRUE,
                            search = "random",
                            savePredictions = "final")#,
                            #summaryFunction = GMSummary)
  
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
         starts_with(feature_prefixes))

if (!is.null(INPUT$select_features)) {  # Need to do feature selection
  if (!all(final_features %in% colnames(final_data)))
    stop("Not all selected features found in data. Does feature selection list match input options?")
  
  final_data <- final_data %>% 
    select(.data$species, .data$ace2_accession, .data$evidence_level, .data$label,
           all_of(final_features))
}

  
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
                       metric = "Accuracy",#"GM",  # Geometric mean of sensitivity and specificity
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