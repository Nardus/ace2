## Plot detailed feature usage in the infection model

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(cowplot)
})

source("scripts/utils/plot_utils.R")


# ---- Data ----------------------------------------------------------------------------------------
full_model <- readRDS("output/all_data/infection/all_features/training_results.rds")
final_model <- readRDS("output/all_data/infection/feature_selection_6/training_results.rds")

internal_data <- final_model$trained_models[[1]]$trainingData %>% 
  mutate(.outcome = factor(.data$.outcome))

dummy_coded_data <- internal_data %>% 
  dummyVars(formula = .outcome ~ ., data = ., sep = "") %>% 
  predict(newdata = internal_data) %>% 
  data.frame() %>% 
  select(-any_of(colnames(final_model$trainingData)))

all_training_data <- final_model$trainingData %>% 
  bind_cols(dummy_coded_data)


# Raw features
#metadata <- readRDS("data/calculated/cleaned_infection_data.rds") %>% 
#  rename(label = .data$infected)

#aa_categorical <- readRDS("data/calculated/features_variable_sites.rds")

#aa_categorical <- aa_categorical %>% 
#  left_join(metadata, by = "ace2_accession")



# Readable feature names
feature_label_data <- data.frame(feature = colnames(full_model$trainingData)) %>% 
  filter(!.data$feature %in% c("species", "ace2_accession", "evidence_level", "label",
                               "variable_site_391V")) %>% 
  add_readable_feature_names()

feature_labels <- feature_label_data$feature_label
names(feature_labels) <- feature_label_data$feature


# Constants:
LABEL <- "Infected"

MAIN_VAR_1 <- "dist_variable_site_412"
MAIN_LABEL_1 <- feature_labels[MAIN_VAR_1]

MAIN_VAR_2 <- "variable_site_391V"
MAIN_LABEL_2 <- sprintf("%s = Valine", feature_labels["variable_site_391"])

EXPANDED_VAR_1 <- "variable_site_412"  # Show observed amino acids
EXPANDED_VAR_2 <- "variable_site_391"
EXPANDED_LABEL_1 <- feature_labels[EXPANDED_VAR_1]
EXPANDED_LABEL_2 <- feature_labels[EXPANDED_VAR_2]


# ---- Feature importance --------------------------------------------------------------------------
shap_values <- get_shap_values(final_model)

feature_importance <- shap_values %>% 
  mutate(full_feature_name = .data$feature,
         feature = str_remove(.data$feature, "[[:alpha:]]$")) %>% 
  group_by(.data$feature, .data$full_feature_name) %>% 
  summarise(mean_shap = mean(abs(.data$shap_value)),
            .groups = "drop_last") %>% 
  summarise(total_effect_magnitude = sum(.data$mean_shap)) %>% 
  add_readable_feature_names() %>% 
  arrange(.data$total_effect_magnitude) %>% 
  mutate(feature_label = factor(.data$feature_label, levels = .data$feature_label))


ggplot(feature_importance, aes(x = feature_label, y = total_effect_magnitude)) +
  geom_col() +
  labs(x = NULL, y = "Total effect magnitude") +
  coord_flip() +
  PLOT_THEME


# ---- SHAP values ---------------------------------------------------------------------------------
shap_values <- shap_values %>% 
  mutate(full_feature_name = .data$feature,
         feature = str_remove(.data$full_feature_name, "[[:alpha:]]$"),
         observed_aa = str_extract(.data$full_feature_name, "[[:alpha:]]$")) %>% 
  add_readable_feature_names() %>% 
  mutate(feature_label = if_else(!is.na(.data$observed_aa),
                                 sprintf("%s = %s", .data$feature_label, .data$observed_aa),
                                 .data$feature_label))


ggplot(shap_values, aes(x = feature_value, y = shap_value, colour = label, shape = label)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_jitter(width = 0.1, height = 0, size = 1.8) + 
  scale_colour_brewer(palette = "Set1", na.value = "grey80") +
  scale_shape_manual(values = c(19, 1)) +
  facet_wrap(vars(feature_label), scales = "free") +
  labs(x = "Observed value", y = "Effect on log-odds",
       title = "SHAP values") +
  PLOT_THEME


# ---- SHAP interaction ----------------------------------------------------------------------------
shap_interaction_values <- get_shap_values(final_model, interaction = TRUE)


p_interaction_1 <- plot_shap_interaction(shap_interaction_values, MAIN_VAR_1, MAIN_VAR_1, 
                                         jitter_width = 5) +
  labs(x = MAIN_LABEL_1,
       y = MAIN_LABEL_1,
       shape = LABEL,
       colour = "SHAP interaction value") +
  guides(shape = FALSE,
         colour = guide_colourbar(title.position = "top")) +
  theme(legend.position = "top")


p_interaction_2 <- plot_shap_interaction(shap_interaction_values, MAIN_VAR_2, MAIN_VAR_1,
                                         jitter_height = 5) +
  labs(x = MAIN_LABEL_2,
       y = MAIN_LABEL_1,
       shape = LABEL,
       colour = "SHAP interaction value") +
  guides(shape = FALSE,
         colour = guide_colourbar(title.position = "top")) +
  theme(legend.position = "top")


# Are interactions important?
interaction_importance <- shap_interaction_values %>% 
  group_by(.data$feature_1, .data$feature_2) %>% 
  summarise(effect_magnitude = mean(abs(.data$shap_value)),
            .groups = "drop") %>%
  filter(effect_magnitude > 0) %>% 
  
  mutate(full_feature_name = .data$feature_1,
         feature = str_remove(.data$full_feature_name, "[[:alpha:]]$"),
         observed_aa = str_extract(.data$full_feature_name, "[[:alpha:]]$")) %>% 
  add_readable_feature_names() %>% 
  mutate(feature_label_1 = if_else(!is.na(.data$observed_aa),
                                   sprintf("%s = %s", .data$feature_label, .data$observed_aa),
                                   .data$feature_label)) %>% 
  
  mutate(full_feature_name = .data$feature_2,
         feature = str_remove(.data$full_feature_name, "[[:alpha:]]$"),
         observed_aa = str_extract(.data$full_feature_name, "[[:alpha:]]$")) %>% 
  add_readable_feature_names() %>% 
  mutate(feature_label_2 = if_else(!is.na(.data$observed_aa),
                                   sprintf("%s = %s", .data$feature_label, .data$observed_aa),
                                   .data$feature_label)) %>% 
  
  mutate(is_interaction = .data$feature_1 != .data$feature_2,
         label = if_else(.data$is_interaction, 
                         sprintf("%s : %s", feature_label_1, feature_label_2),
                         .data$feature_label_1)) %>% 
  arrange(.data$effect_magnitude) %>% 
  mutate(label = factor(.data$label, levels = .data$label))
  


ggplot(interaction_importance, aes(x = label, y = effect_magnitude, fill = is_interaction, colour = is_interaction)) +
  geom_col() +
  coord_flip() +
  labs(x = NULL, y = "Effect magnitude", fill = "Interaction", colour = "Interaction") +
  PLOT_THEME

# Tree structure: how deep are trees?
tree_structure <- xgboost::xgb.model.dt.tree(feature_names = final_model$trained_models[[1]]$finalModel$feature_names, 
                                             model = final_model$trained_models[[1]]$finalModel)

feats_per_tree <- data.frame(tree_structure) %>% 
  filter(.data$Feature != "Leaf") %>% 
  group_by(.data$Tree) %>% 
  summarise(n_feats = n_distinct(.data$Feature))

ggplot(feats_per_tree, aes(x = n_feats)) +
  geom_histogram() +
  xlim(0, 4) +
  labs(x = "Number of features", y = "Number of trees") +
  PLOT_THEME


node_gain <- data.frame(tree_structure) %>% 
  filter(.data$Feature != "Leaf") %>% 
  group_by(.data$Tree) %>% 
  mutate(node_index = 1:n(),
         feature_number = as.numeric(factor(.data$Feature, levels = unique(.data$Feature))))

ggplot(node_gain, aes(x = node_index, y = Quality, colour = factor(feature_number), group = Tree)) +
  geom_line(colour = "grey80") +
  geom_point() +
  labs(x = "Node index", y = "Gain", colour = "Feature number") +
  PLOT_THEME



## Try a gam with these features
library(mgcv)

gam_data <- internal_data %>% 
  mutate(.outcome = .data$.outcome == "True")# %>% 
  #mutate(across(where(is.character), as.factor))

# Remove duplicates (feat1:feat2 is the same as feat2:feat1)
terms <- interaction_importance %>% 
  filter(.data$effect_magnitude > 1e-4) %>% 
  rowwise() %>% 
  mutate(combined = paste(sort(c(.data$feature_1, .data$feature_2)), collapse = "/")) %>% 
  distinct(.data$combined) %>% 
  separate(.data$combined, into = c("feature_1", "feature_2"), sep = "/")


gam_formula <- terms %>% 
  mutate(term = case_when(.data$feature_1 == .data$feature_2 & startsWith(.data$feature_1, "dist") ~ sprintf("s(%s, bs = 'tp', k = 6)", .data$feature_1),
                          .data$feature_1 == .data$feature_2 ~ .data$feature_1,
                          startsWith(.data$feature_1, "dist") & startsWith(.data$feature_2, "dist") ~ sprintf("s(%s, %s, bs = 'tp', k = 6)", .data$feature_1, .data$feature_2),
                          startsWith(.data$feature_1, "dist") ~ sprintf("s(%s, by = %s, bs = 'tp', k = 6)", .data$feature_1, .data$feature_2),
                          startsWith(.data$feature_2, "dist") ~ sprintf("s(%s, by = %s, bs = 'tp', k = 6)", .data$feature_2, .data$feature_1),
                          TRUE ~ sprintf("%s:%s", .data$feature_1, .data$feature_2))) %>% 
  pull(.data$term) %>% 
  paste(collapse = " + ") %>% 
  paste("label ~", .) %>% 
  as.formula()

gam_fit <- gam(gam_formula, data = all_training_data)








# ---- Observed values------------------------------------------------------------------------------
human <- all_training_data %>% 
  filter(.data$species == "Homo sapiens")

p_obs_values <- ggplot(all_training_data, aes_string(x = MAIN_VAR_1, y = MAIN_VAR_2, 
                                                     colour = "label", shape = "factor(evidence_level)")) +
  geom_vline(aes_string(xintercept = MAIN_VAR_1), linetype = 2, data = human) +
  geom_hline(aes_string(yintercept = MAIN_VAR_2), linetype = 2, data = human) +
  geom_jitter(width = 5, height = 0.05, size = 1.8) +
  scale_colour_brewer(palette = "Set1", guide = FALSE) + 
  labs(x = feature_labels[MAIN_VAR_1],
       y = feature_labels[MAIN_VAR_2],
       colour = LABEL,
       shape = "Evidence level") +
  guides(shape = guide_legend(title.position = "top")) +
  PLOT_THEME +
  theme(legend.position = "top")


# Actual amino acids - include holdout viruses
human <- aa_categorical %>% 
  filter(.data$species == "Homo sapiens")

training_data <- aa_categorical %>% 
  filter(!is.na(.data$label))

holdout_data <- aa_categorical %>% 
  filter(is.na(.data$label))

p_obs_aminoacids <- ggplot(training_data, aes_string(x = EXPANDED_VAR_1, y = EXPANDED_VAR_2, 
                                                     colour = "label", shape = "factor(evidence_level)")) +
  geom_vline(xintercept = human[[EXPANDED_VAR_1]], linetype = 2) +
  geom_hline(yintercept = human[[EXPANDED_VAR_2]], linetype = 2) +
  
  geom_jitter(width = 0.1, height = 0.1, size = 1.8, data = holdout_data) +
  geom_jitter(width = 0.1, height = 0.1, size = 1.8) +
  
  scale_colour_brewer(palette = "Set1", guide = FALSE, na.value = "grey70") + 
  scale_shape_discrete(na.value = 0) + 
  
  facet_grid(cols = vars(label)) +
  
  labs(x = EXPANDED_LABEL_1,
       y = EXPANDED_LABEL_2,
       colour = LABEL,
       shape = "Evidence level") +
  
  guides(shape = guide_legend(title.position = "top")) +
  PLOT_THEME +
  theme(legend.position = "top")


# ---- Combine-------------------------------------------------------------------------------------
p_detail <- plot_grid(p_shap_1, p_shap_2,
                      p_interaction_1, p_interaction_2,
                      p_obs_values, p_obs_aminoacids,
                      ncol = 2,
                      labels = LETTERS[1:5],
                      align = "hv", axis = "lrt")



ggsave2("output/plots/varimp_detail_infection.png", p_detail,
        width = 6, height = 6)
