# Create a dendrogram clustering our predictions alongside published ones

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(readr)
  library(scales)
  library(cluster)
  library(ape)
})

existing_predictions <- read_csv("data/internal/existing_predictions.csv",
                                 col_types = cols(raw_score = "d",
                                                  reverse_score = "l",
                                                  .default = "c"))

taxonomy <- readRDS("data/calculated/taxonomy.rds")

# Using holdout predictions and "fitted values", to match earlier studies (which had no cross-validation)
ace2_preds <- readRDS("output/all_data/infection/all_features/holdout_predictions.rds")
phylo_preds <- readRDS("output/all_data/infection/phylogeny/holdout_predictions.rds")
ensemble_preds <- readRDS("output/all_data/infection/ensemble_aa_distance_self/holdout_predictions.rds")


# ---- Fix taxonomy --------------------------------------------------------------------------------
known_invalid <- c("Nyctophilus timoriensis", "Pseudopotto martini") # Invalid/dubious species in Fischoff et al. data (see ITIS records)
stopifnot(all(existing_predictions$species %in% c(taxonomy$internal_name, known_invalid)))

taxonomy <- taxonomy %>% 
  select(.data$internal_name, .data$species)

existing_predictions <- existing_predictions %>%
  rename(internal_name = .data$species) %>% 
  left_join(taxonomy, by = "internal_name") %>% 
  filter(!is.na(.data$species))


# Summarise over all subspecies in the same study
# - For categorical predictions, the most common value (generally the same across subspecies anyway)
# - For numeric scores, the mean 
categorical_preds <- existing_predictions %>% 
  group_by(.data$citation_key, .data$species, .data$prediction) %>% 
  summarise(n = n(), .groups = "drop_last") %>% 
  filter(.data$n == max(.data$n)) %>% 
  sample_n(size = 1) %>%  # in case of ties
  ungroup() %>% 
  select(-.data$n)
  
existing_predictions <- existing_predictions %>% 
  group_by(.data$citation_key, .data$prediction_type, .data$predictor, .data$reverse_score,
           .data$species) %>% 
  summarise(raw_score = mean(.data$raw_score),  
            .groups = "drop") %>% 
  left_join(categorical_preds, by = c("citation_key", "species"))


# ---- Add our values ------------------------------------------------------------------------------
# Unify column names
ace2_preds <- ace2_preds %>% 
  select(.data$species, 
         prediction = .data$predicted_label,
         raw_score = .data$probability) %>% 
  mutate(citation_key = "This study (All ACE2 representations)",
         predictor = "All ACE2 representations",
         prediction_type = "sequence-based")


phylo_preds <- phylo_preds %>% 
  select(.data$species, 
         prediction = .data$predicted_label,
         raw_score = .data$probability) %>% 
  mutate(citation_key = "This study (host phylogeny)",
         predictor = "host phylogeny",
         prediction_type = "phylogeny")


ensemble_preds <- ensemble_preds %>% 
  select(.data$species, 
         prediction = .data$predicted_label,
         raw_score = .data$probability) %>% 
  mutate(citation_key = "This study (ACE2 ensemble)",
         predictor = "ACE2 ensemble",
         prediction_type = "ensemble")


# Merge
internal_predictions <- ace2_preds %>% 
  bind_rows(phylo_preds) %>% 
  bind_rows(ensemble_preds) %>% 
  mutate(reverse_score = FALSE)

all_predictions <- existing_predictions %>% 
  bind_rows(internal_predictions)


# ---- Make scores comparable ----------------------------------------------------------------------
# - All scores should be ascending, meaning a higher value = more likely to be positive
# - Rescale values from each study to lie in [0, 1], so we can use a single colour scale in plots
all_predictions <- all_predictions %>% 
  group_by(.data$citation_key) %>% 
  mutate(scaled_score = if_else(.data$reverse_score, .data$raw_score * -1, .data$raw_score),
         scaled_score = rescale(.data$scaled_score, to = c(0, 1))) %>% 
  ungroup()


# ---- Correlation between studies -----------------------------------------------------------------
score_mat <- all_predictions %>% 
  select(.data$species, .data$citation_key, .data$scaled_score) %>% 
  pivot_wider(id_cols = "species", names_from = "citation_key", values_from = "scaled_score") %>% 
  column_to_rownames("species")

cor_mat <- cor(score_mat, use = "pairwise.complete", method = "spearman")

# Assume studies with no overlap have no correlation:
cor_mat <- replace(cor_mat, is.na(cor_mat), 0)

# Hierarchical clustering by correlation
prediction_clusters <- agnes(1 - cor_mat, diss = TRUE)
prediction_dendrogram <- as.phylo(as.hclust(prediction_clusters))


# ---- Output --------------------------------------------------------------------------------------
saveRDS(prediction_dendrogram, "output/plots/intermediates/prediction_dendrogram.rds")
saveRDS(all_predictions, "output/plots/intermediates/unified_predictions.rds")
saveRDS(cor_mat, "output/plots/intermediates/study_correlation_matrix.rds")
