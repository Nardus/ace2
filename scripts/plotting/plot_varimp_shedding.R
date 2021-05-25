## Plot feature usage in the infection model

library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(ggplot2)
library(cowplot)

source("scripts/utils/plot_utils.R")


# ---- Data ----------------------------------------------------------------------------------------
full_model <- readRDS("output/all_data/shedding/all_features/training_results.rds")
final_model <- readRDS("output/all_data/shedding/feature_selection_2/training_results.rds")

all_training_data <- full_model$trainingData %>% 
  mutate(variable_site_886G = .data$variable_site_886 == "G",
         variable_site_120E = .data$variable_site_120 == "E",)

feature_importance <- readRDS("output/all_data/shedding/all_features/feature_usage.rds")


# Readable feature names
feature_labels <- feature_importance %>% 
  add_readable_feature_names() %>% 
  pull(.data$feature_label)

names(feature_labels) <- feature_importance$feature


# Constants:
LABEL <- "Shedding"

MAIN_VAR_1 <- "variable_site_886G"
MAIN_LABEL_1 <- sprintf("%s = Glycine", feature_labels["variable_site_886"])

MAIN_VAR_2 <- "variable_site_120E"
MAIN_LABEL_2 <- sprintf("%s = Glutamic acid", feature_labels["variable_site_120"])

EXPANDED_VAR_1 <- "variable_site_886"
EXPANDED_VAR_2 <- "variable_site_120"
EXPANDED_LABEL_1 <- feature_labels[EXPANDED_VAR_1]
EXPANDED_LABEL_2 <- feature_labels[EXPANDED_VAR_2]


# ---- Feature importance --------------------------------------------------------------------------
feature_importance_nz <- feature_importance %>% 
  filter(.data$importance > 0) %>% 
  add_readable_feature_names() %>% 
  arrange(.data$importance) %>% 
  mutate(feature_label = factor(.data$feature_label, levels = .data$feature_label))

p_importance <- ggplot(feature_importance_nz, aes(x = feature_label, y = importance)) +
  geom_col() +
  coord_flip() +
  labs(y = "Importance", x = NULL) +
  PLOT_THEME


# ---- SHAP values ---------------------------------------------------------------------------------
shap_values <- get_shap_values(final_model)

# Plot
p_shap_1 <- plot_shap(shap_values, MAIN_VAR_1, x_is_factor = TRUE) +
  labs(x = MAIN_LABEL_1,
       y = "SHAP value",
       colour = LABEL,
       shape = LABEL) +
  guides(colour = guide_legend(title.position = "top")) +
  theme(legend.position = "top")


p_shap_2 <- plot_shap(shap_values, MAIN_VAR_2, x_is_factor = TRUE) +
  labs(x = MAIN_LABEL_2,
       y = "SHAP value",
       colour = LABEL,
       shape = LABEL) +
  guides(colour = guide_legend(title.position = "top")) +
  theme(legend.position = "top")


# ---- SHAP interaction ----------------------------------------------------------------------------
p_interaction_1 <- plot_shap_interaction(shap_values, MAIN_VAR_1, MAIN_VAR_2, 
                                         x_is_factor = TRUE, y_is_factor = TRUE,
                                         jitter_width = 0.1, jitter_height = 0.1) +
  labs(x = MAIN_LABEL_1,
       y = MAIN_LABEL_2,
       shape = LABEL,
       colour = sprintf("SHAP value: %s", MAIN_LABEL_1)) +
  guides(shape = FALSE,
         colour = guide_colourbar(title.position = "top")) +
  theme(legend.position = "top")


p_interaction_2 <- plot_shap_interaction(shap_values, MAIN_VAR_2, MAIN_VAR_1,
                                         x_is_factor = TRUE, y_is_factor = TRUE,
                                         jitter_width = 0.1, jitter_height = 0.1) +
  labs(x = MAIN_LABEL_2,
       y = MAIN_LABEL_1,
       shape = LABEL,
       colour = sprintf("SHAP value: %s", MAIN_LABEL_2)) +
  guides(shape = FALSE,
         colour = guide_colourbar(title.position = "top")) +
  theme(legend.position = "top")


# ---- Key amino acids------------------------------------------------------------------------------
human <- all_training_data %>% 
  filter(.data$species == "Homo sapiens")

p_obs_values <- ggplot(all_training_data, aes(x = dist_variable_site_886, y = dist_variable_site_120, 
                                              colour = label, shape = factor(evidence_level))) +
  geom_vline(aes(xintercept = dist_variable_site_886), linetype = 2, data = human) +
  geom_hline(aes(yintercept = dist_variable_site_120), linetype = 2, data = human) +
  geom_jitter(width = 5, height = 5, size = 1.8) +
  scale_colour_brewer(palette = "Set1", guide = FALSE) + 
  labs(x = feature_labels["dist_variable_site_886"],
       y = feature_labels["dist_variable_site_120"],
       colour = LABEL,
       shape = "Evidence level") +
  guides(shape = guide_legend(title.position = "top")) +
  PLOT_THEME +
  theme(legend.position = "top")


p_obs_aminoacids <- ggplot(all_training_data, aes_string(x = EXPANDED_VAR_1, y = EXPANDED_VAR_2, 
                                                         colour = "label", shape = "factor(evidence_level)")) +
  geom_vline(aes_string(xintercept = EXPANDED_VAR_1), linetype = 2, data = human) +
  geom_hline(aes_string(yintercept = EXPANDED_VAR_2), linetype = 2, data = human) +
  geom_jitter(width = 0.1, height = 0.1, size = 1.8) +
  scale_colour_brewer(palette = "Set1", guide = FALSE) + 
  labs(x = EXPANDED_LABEL_1,
       y = EXPANDED_LABEL_1,
       colour = LABEL,
       shape = "Evidence level") +
  guides(shape = guide_legend(title.position = "top")) +
  PLOT_THEME +
  theme(legend.position = "top")


# ---- Combine-------------------------------------------------------------------------------------
bottom_rows <- plot_grid(p_shap_1, p_shap_2,
                         p_interaction_1, p_interaction_2,
                         p_obs_values, p_obs_aminoacids,
                         ncol = 2,
                         labels = LETTERS[2:7],
                         align = "hv", axis = "lrt")

p_final <- plot_grid(p_importance,
                     bottom_rows,
                     ncol = 1, rel_heights = c(0.5, 3), 
                     labels = c("A", ""))

ggsave2("output/plots/varimp_shedding.png", p_final,
        width = 6, height = 8)
