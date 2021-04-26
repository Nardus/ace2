## Predict additional species using trained models

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

set.seed(2312532)
registerDoParallel(cores = 16)

# ---- Models --------------------------------------------------------------------------------------
infection_model_path <- "output/all_data/infection/feature_selection_2/"
shedding_model_path <- "output/all_data/shedding/feature_selection_2/"

model_infection <- readRDS(file.path(infection_model_path, "training_results.rds"))
model_shedding <- readRDS(file.path(shedding_model_path, "training_results.rds"))


# ---- Data ----------------------------------------------------------------------------------------
# Species with known responses
infection_data <- read_rds("data/calculated/cleaned_infection_data.rds")
shedding_data <- read_rds("data/calculated/cleaned_shedding_data.rds")

training_metadata <- infection_data %>% 
  full_join(shedding_data, by = c("species", "ace2_accession")) %>% 
  select(.data$species, .data$ace2_accession, .data$infected, .data$shedding)

# Other species
# - Removing both matched species names AND accessions already used to represent related species
additional_metadata <- read_csv("data/internal/NCBI_ACE2_orthologs.csv",
                                col_types = cols(.default = "c")) %>% 
  select(species = .data$`Scientific name`, 
         ace2_accession = .data$`RefSeq Protein accessions`) %>% 
  filter(!.data$species %in% training_metadata$species & 
           !.data$ace2_accession %in% training_metadata$ace2_accession)

final_metadata <- training_metadata %>% 
  bind_rows(additional_metadata)


# ---- Features ------------------------------------------------------------------------------------
# General features
pairwise_dist_data <- read_rds("data/calculated/features_pairwise_dists.rds")
dist_to_humans <- read_rds("data/calculated/features_dist_to_humans.rds")
variable_sites <- read_rds("data/calculated/features_variable_sites.rds")
haddock_scores <- read_rds("data/calculated/features_haddock_scores.rds")


# ---- Predict infection ---------------------------------------------------------------------------
predict_additional <- function(model, label, all_data) {
  training_data <- model$trainingData
  
  label_unknown <- is.na(all_data[[label]])
  
  new_data <- all_data[label_unknown, ] %>% 
    select(.data$species, .data$ace2_accession)
  
  # Training set-specific features
  stopifnot(all(new_data$ace2_accession %in% pairwise_dist_data$ace2_accession))
  stopifnot(all(new_data$ace2_accession %in% variable_sites$ace2_accession))
  
  closest_positive <- get_dist_to_closest_positive(pairwise_dist_data, training_data)
  consensus_dists <- get_consensus_dist(variable_sites, variable_sites, training_data)
  
  # Combine
  new_features <- new_data %>% 
    left_join(dist_to_humans, by = "ace2_accession") %>% 
    left_join(variable_sites, by = "ace2_accession") %>% 
    left_join(haddock_scores, by = "species") %>% 
    left_join(closest_positive, by = "ace2_accession") %>% 
    left_join(consensus_dists, by = "ace2_accession")
  
  # Predict
  new_data$predicted_prob <- predict(model$finalModel, newdata = new_features, type = "raw")
  
  new_data %>% 
    mutate(predicted_label = if_else(.data$predicted_prob > cutoff, "True", "False"))
}


preds_infection <- predict_additional(model_infection,
                                      label = "infected",
                                      all_data = final_metadata)

preds_shedding <- predict_additional(model_shedding,
                                     label = "shedding",
                                     all_data = final_metadata)

# ---- Output ---------------------------------------------------------------------------
saveRDS(preds_infection, file.path(infection_model_path, "additional_preds_infection.rds"))
saveRDS(preds_shedding, file.path(shedding_model_path, "additional_preds_shedding.rds"))
