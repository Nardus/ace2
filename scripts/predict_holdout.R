## Predict additional species using trained models

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(workflows)
  library(yardstick)
  
  source("scripts/utils/aa_distance_utils.R")
  source("scripts/utils/feature_calc_utils.R")
  source("scripts/utils/feature_recipe_utils.R")
  source("scripts/utils/cutoff_utils.R")
})

set.seed(2312532)


# ---- Model ---------------------------------------------------------------------------------------
ace2_workflow <- read_rds("output/all_data/infection/all_features/trained_model_workflow.rds")
phylo_workflow <- read_rds("output/all_data/infection/phylogeny/trained_model_workflow.rds")


# ---- Metadata ------------------------------------------------------------------------------------
# Used for training
training_metadata <- read_rds("data/calculated/cleaned_infection_data.rds") %>% 
  rename(label = .data$infected)

# Other species
# - Removing both matched species names AND accessions already used to represent related species
additional_metadata <- read_csv("data/internal/NCBI_ACE2_orthologs.csv",
                                col_types = cols(.default = "c")) %>% 
  select(species = .data$`Scientific name`, 
         ace2_accession = .data$`RefSeq Protein accessions`) %>% 
  filter(!.data$species %in% training_metadata$species & 
           !.data$ace2_accession %in% training_metadata$ace2_accession)

final_metadata <- training_metadata %>% 
  bind_rows(additional_metadata) %>% 
  mutate(label = factor(.data$label, levels = c("True", "False")))


# ---- Feature data --------------------------------------------------------------------------------
pairwise_dist_data <- read_rds("data/calculated/features_pairwise_dists.rds")
dist_to_humans <- read_rds("data/calculated/features_dist_to_humans.rds")
variable_sites <- read_rds("data/calculated/features_variable_sites.rds")
site_properties <- read_rds("data/calculated/features_site_properties.rds")
haddock_scores <- read_rds("data/calculated/features_haddock_scores.rds")
phylogeny_features <- read_rds("data/calculated/features_phylogeny_eigenvectors.rds")

# Combine 
final_data <- final_metadata %>% 
  left_join(dist_to_humans, by = "ace2_accession") %>% 
  left_join(variable_sites, by = "ace2_accession") %>% 
  left_join(site_properties, by = "ace2_accession") %>% 
  left_join(haddock_scores, by = "species") %>% 
  left_join(phylogeny_features, by = "species")

stopifnot(nrow(final_data) == n_distinct(final_metadata$species))

# Expanded dataset for phylogeny model
# - Here we can generate predictions for all species, not just those with ACE2 sequences available
expanded_phylo_data <- final_metadata %>% 
  left_join(dist_to_humans, by = "ace2_accession") %>% 
  left_join(variable_sites, by = "ace2_accession") %>% 
  left_join(site_properties, by = "ace2_accession") %>% 
  left_join(haddock_scores, by = "species") %>% 
  full_join(phylogeny_features, by = "species")  # Different from "final_data" above
  
# Final processing
final_data <- final_data %>% 
  mutate_if(is.character, as.factor)

expanded_phylo_data <- expanded_phylo_data %>% 
  mutate_if(is.character, as.factor)


# ---- Continuous predictions ----------------------------------------------------------------------
ace2_predictions <- final_metadata
ace2_predictions$probability <- predict(ace2_workflow, final_data, type = "raw")

phylo_predictions <- expanded_phylo_data %>% 
  select(.data$species) %>% 
  full_join(final_metadata, by = "species")
  
phylo_predictions$probability <- predict(phylo_workflow, expanded_phylo_data, type = "raw")


# ---- Find best cutoff ----------------------------------------------------------------------------
# Use training data to find optimal cut-off for novel species
ace2_cutoff <- find_best_cuttof(ace2_predictions)
phylo_cutoff <- find_best_cuttof(phylo_predictions)


# ---- Discrete predictions ------------------------------------------------------------------------
add_discrete_predictions <- function(predictions, best_cutoff) {
  predictions %>% 
    mutate(prediction_type = if_else(!is.na(.data$label), "Fitted value", "Prediction"),
           predicted_label = if_else(.data$probability > best_cutoff, "True", "False"),
           cutoff = best_cutoff)
}

ace2_predictions <- ace2_predictions %>% 
  add_discrete_predictions(ace2_cutoff) %>% 
  select(.data$species, .data$ace2_accession, .data$label, 
         .data$prediction_type, .data$predicted_label, .data$probability, .data$cutoff)

phylo_predictions <- phylo_predictions %>% 
  add_discrete_predictions(phylo_cutoff) %>% 
  select(.data$species, .data$label, 
         .data$prediction_type, .data$predicted_label, .data$probability, .data$cutoff)


# ---- Output ---------------------------------------------------------------------------
saveRDS(ace2_predictions, "output/all_data/infection/all_features/holdout_predictions.rds")
saveRDS(phylo_predictions, "output/all_data/infection/phylogeny/holdout_predictions.rds")
