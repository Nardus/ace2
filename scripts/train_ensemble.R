# Evaluate an ensemble model combining the ACE2-based and phylogeny-based predictions
# - Since we have trained models for each CV fold (and these match since we used LOO-CV), 
#   there's no need to retrain the base models
# - However, we do need to test our cutoff-finding method on the averaged predictions,
#   using the same folds

# NOTE: this script uses a lot of memory - reduce N_CORES if needed to control this

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
set.seed(10012022)

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
base_model_ace2 <- read_rds("output/all_data/infection/all_features/cv_models.rds")
base_model_phylo <- read_rds("output/all_data/infection/phylogeny/cv_models.rds")

# Base predictions
base_preds_ace2 <- read_rds("output/all_data/infection/all_features/predictions.rds")
base_preds_phylo <- read_rds("output/all_data/infection/phylogeny/predictions.rds")


# ---- Pre-processing ------------------------------------------------------------------------------
# Combine predictions
base_preds_combined <- base_preds_ace2 %>% 
  select(.data$species, .data$label,
         ace2_pred = .data$p_true) %>% 
  left_join(base_preds_phylo, by = c("species", "label")) %>% 
  select(.data$species, .data$label, .data$ace2_pred,
         phylo_pred = .data$p_true) 

# Find holdout virus of each fold
# - Order may differ between models
inds_ace2 <- base_preds_ace2$cv_fold
names(inds_ace2) <- base_preds_ace2$species

inds_phylo <- base_preds_phylo$cv_fold
names(inds_phylo) <- base_preds_phylo$species

# Match models to holdout viruses
stopifnot(all(names(inds_ace2) %in% names(inds_phylo)))

base_models <- lapply(names(inds_ace2), function(sp) list(species = sp,
                                                          ace2 = base_model_ace2[[inds_ace2[sp]]],
                                                          phylo = base_model_phylo[[inds_phylo[sp]]]))


# ---- Cross-validate ensemble model ---------------------------------------------------------------
# Training (find optimal cutoff in each fold)
train_ensemble <- function(models, all_data = final_data, base_preds = base_preds_combined) {
  test_species <- models$species
  ace2_model <- models$ace2
  phylo_model <- models$phylo
  
  # Get best cutoff
  training_data <- all_data %>% 
    filter(.data$species != test_species)
  
  combined_predictions <- training_data %>% 
    select(.data$species, .data$label) %>% 
    mutate(ace2_pred = predict(ace2_model, new_data = training_data, type = "raw"),
           phylo_pred = predict(phylo_model, new_data = training_data, type = "raw")) %>% 
    group_by(.data$species, .data$label) %>% 
    summarise(probability = mean(c(.data$ace2_pred, .data$phylo_pred)), 
              .groups = "drop")
  
  cutoff <- find_best_cuttof(combined_predictions)
  
  # Return test prediction
  base_preds %>% 
    filter(.data$species == test_species) %>% 
    mutate(probability = mean(c(.data$ace2_pred, .data$phylo_pred)),
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

dir.create("output/all_data/infection/ensemble/", recursive = TRUE)
saveRDS(test_predictions, "output/all_data/infection/ensemble/raw_predictions.rds")
saveRDS(cleaned_predictions, "output/all_data/infection/ensemble/predictions.rds")
