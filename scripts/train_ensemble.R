# Evaluate an ensemble model combining the predictions from two models
# - Since we have trained models for each CV fold (and these match since we used LOO-CV), 
#   there's no need to retrain the base models
# - However, we do need to test our cutoff-finding method on the averaged predictions,
#   using the same folds

# NOTE: this script uses a lot of memory - reduce N_CORES if needed to control this

# ---- Input args ----------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(argparse)
})

parser <- ArgumentParser(description = "Train an ensemble averaging predictions accross two models")

parser$add_argument("--m1", type = "character", 
                    help = "Path to the output of an existing model")
                    
parser$add_argument("--m2", type = "character", 
                    help = "Path to the output of another model")

parser$add_argument("--output_path", type = "character",
                    help = "location for output files. Missing folders will be created.")
                    
parser$add_argument("--random_seed", type = "integer", 
                    help = "Random seed")

INPUT <- parser$parse_args()


# ---- Other setup ---------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(workflows)
  library(yardstick)
  library(parallel)
  
  source("scripts/utils/aa_distance_utils.R")
  source("scripts/utils/feature_calc_utils.R")
  source("scripts/utils/feature_recipe_utils.R")
  source("scripts/utils/cutoff_utils.R")
})

N_CORES = 12  # Number of threads run in parallel
set.seed(INPUT$random_seed)


# ---- Data & models -------------------------------------------------------------------------------
# Data
metadata <- read_rds("data/calculated/cleaned_infection_data.rds") %>% 
  mutate(label = factor(.data$infected, levels = c("True", "False")))

pairwise_dist_data <- read_rds("data/calculated/features_pairwise_dists.rds")
dist_to_humans <- read_rds("data/calculated/features_dist_to_humans.rds")
variable_sites <- read_rds("data/calculated/features_variable_sites.rds")
site_properties <- read_rds("data/calculated/features_site_properties.rds")
binding_affinity <- read_rds("data/calculated/features_binding_affinity.rds")
phylogeny_features <- read_rds("data/calculated/features_phylogeny_eigenvectors.rds")

# Combine
final_data <- metadata %>% 
  left_join(dist_to_humans, by = "ace2_accession") %>% 
  left_join(variable_sites, by = "ace2_accession") %>% 
  left_join(site_properties, by = "ace2_accession") %>% 
  left_join(binding_affinity, by = "species") %>% 
  left_join(phylogeny_features, by = "species")

stopifnot(nrow(final_data) == n_distinct(metadata$species))


# Read base models
base_model_m1 <- read_rds(file.path(INPUT$m1, "cv_models.rds"))
base_model_m2 <- read_rds(file.path(INPUT$m2, "cv_models.rds"))

# Base predictions
base_preds_m1 <- read_rds(file.path(INPUT$m1, "predictions.rds"))
base_preds_m2 <- read_rds(file.path(INPUT$m2, "predictions.rds"))


# ---- Pre-processing ------------------------------------------------------------------------------
# Combine predictions
base_preds_combined <- base_preds_m1 %>% 
  select(.data$species, .data$label,
         m1_pred = .data$p_true) %>% 
  left_join(base_preds_m2, by = c("species", "label")) %>% 
  select(.data$species, .data$label, .data$m1_pred,
         m2_pred = .data$p_true) 

# Find holdout virus of each fold
# - Order may differ between models
inds_m1 <- base_preds_m1$cv_fold
names(inds_m1) <- base_preds_m1$species

inds_m2 <- base_preds_m2$cv_fold
names(inds_m2) <- base_preds_m2$species

# Match models to holdout viruses
stopifnot(all(names(inds_m1) %in% names(inds_m2)))

base_models <- lapply(names(inds_m1), function(sp) list(species = sp,
                                                        model_1 = base_model_m1[[inds_m1[sp]]],
                                                        model_2 = base_model_m2[[inds_m2[sp]]]))


# ---- Cross-validate ensemble model ---------------------------------------------------------------
# Training (find optimal cutoff in each fold)
train_ensemble <- function(models, all_data = final_data, base_preds = base_preds_combined) {
  test_species <- models$species
  m1_model <- models$model_1
  m2_model <- models$model_2
  
  # Get best cutoff
  training_data <- all_data %>% 
    filter(.data$species != test_species)
  
  combined_predictions <- training_data %>% 
    select(.data$species, .data$label) %>% 
    mutate(m1_pred = predict(m1_model, new_data = training_data, type = "raw"),
           m2_pred = predict(m2_model, new_data = training_data, type = "raw")) %>% 
    group_by(.data$species, .data$label) %>% 
    summarise(probability = mean(c(.data$m1_pred, .data$m2_pred)), 
              .groups = "drop")
  
  cutoff <- find_best_cuttof(combined_predictions)
  
  # Return test prediction
  base_preds %>% 
    filter(.data$species == test_species) %>% 
    mutate(probability = mean(c(.data$m1_pred, .data$m2_pred)),
           cutoff = cutoff,
           prediction = if_else(.data$probability > cutoff, "True", "False"))
}

# LOO-CV, so can re-use splits from base models
test_predictions <- mclapply(base_models, train_ensemble, 
                             mc.cores = N_CORES)

test_predictions <- test_predictions %>% 
  bind_rows()


# ---- Output --------------------------------------------------------------------------------------
# Make a version matching the predictions output by train_models.R:
cleaned_predictions <- test_predictions %>% 
  mutate(cv_fold = .data$species, 
         p_true = .data$probability,
         p_false = 1 - .data$p_true) %>% 
  select(.data$species, .data$label, .data$cv_fold, .data$prediction, 
         .data$cutoff, .data$p_true, .data$p_false)

dir.create(INPUT$output_path, recursive = TRUE)
saveRDS(test_predictions, file.path(INPUT$output_path, "raw_predictions.rds"))
saveRDS(cleaned_predictions, file.path(INPUT$output_path, "predictions.rds"))
