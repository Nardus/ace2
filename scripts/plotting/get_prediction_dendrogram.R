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

ace2_cv_preds <- readRDS("output/all_data/infection/all_features/predictions.rds")
phylo_cv_preds <- readRDS("output/all_data/infection/phylogeny/predictions.rds")
ensemble_cv_preds <- readRDS("output/all_data/infection/ensemble/predictions.rds")

ace2_holdout_preds <- readRDS("output/all_data/infection/all_features/holdout_predictions.rds") %>% 
  filter(.data$prediction_type != "Fitted value")

phylo_holdout_preds <- readRDS("output/all_data/infection/phylogeny/holdout_predictions.rds") %>% 
  filter(.data$prediction_type != "Fitted value")

ensemble_holdout_preds <- readRDS("output/all_data/infection/ensemble/holdout_predictions.rds") %>% 
  filter(.data$prediction_type != "Fitted value")


# ---- Fix taxonomy --------------------------------------------------------------------------------
stopifnot(all(existing_predictions$species %in% taxonomy$internal_name))

taxonomy <- taxonomy %>% 
  select(.data$internal_name, .data$species)

existing_predictions <- existing_predictions %>%
  rename(internal_name = .data$species) %>% 
  left_join(taxonomy, by = "internal_name")


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
ace2_cv_preds <- ace2_cv_preds %>% 
  select(.data$species, .data$prediction,
         raw_score = .data$p_true)

phylo_cv_preds <- phylo_cv_preds %>% 
  select(.data$species, .data$prediction,
         raw_score = .data$p_true)

ensemble_cv_preds <- ensemble_cv_preds %>% 
  select(.data$species, .data$prediction,
         raw_score = .data$p_true)


ace2_holdout_preds <- ace2_holdout_preds %>% 
  select(.data$species, 
         prediction = .data$predicted_label,
         raw_score = .data$probability)

phylo_holdout_preds <- phylo_holdout_preds %>% 
  select(.data$species, 
         prediction = .data$predicted_label,
         raw_score = .data$probability)

ensemble_holdout_preds <- ensemble_holdout_preds %>% 
  select(.data$species, 
         prediction = .data$predicted_label,
         raw_score = .data$probability)


# Merge
ace2_predictions <- ace2_cv_preds %>% 
  bind_rows(ace2_holdout_preds) %>% 
  mutate(citation_key = "This study (ACE2-based)",
         predictor = "ACE2-based",
         prediction_type = "sequence-based")

phylo_predictions <- phylo_cv_preds %>% 
  bind_rows(phylo_holdout_preds) %>% 
  mutate(citation_key = "This study (host phylogeny)",
         predictor = "host phylogeny",
         prediction_type = "phylogeny")

ensemble_predictions <- ensemble_cv_preds %>% 
  bind_rows(ensemble_holdout_preds) %>% 
  mutate(citation_key = "This study (ensemble)",
         predictor = "ensemble",
         prediction_type = "ensemble")

internal_predictions <- ace2_predictions %>% 
  bind_rows(phylo_predictions) %>% 
  bind_rows(ensemble_predictions) %>% 
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
