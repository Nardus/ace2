## Plot overview of feature usage

suppressPackageStartupMessages({
  library(argparse)
})

parser <- ArgumentParser(description = "Plot an overview of feature usage")

parser$add_argument("varimp_file", type = "character",
                    help = "location of (iteration level) variable importance data")

parser$add_argument("output_name", type = "character",
                    help = "location/name for output file")

INPUT <- parser$parse_args()


# ---- Data ----------------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(cowplot)
  
  source("scripts/utils/plot_utils.R")
})

feature_importance <- readRDS(INPUT$varimp_file)


feature_clusters <- readRDS("output/plots/feature_clusters_infection.rds") # TODO: should come from outside


# ---- Mark sites known to interact with S ---------------------------------------------------------
# Binding sites taken from human ACE2 genbank entry (NP_001358344.1)
s_binding_sites <- tribble(
  ~start_pos, ~stop_pos, ~name,
  30,         41,        "ECO:0000269|PubMed:15791205 1",
  82,         84,        "ECO:0000269|PubMed:15791205 2",
  353,        357,       "ECO:0000269|PubMed:15791205 3"
)

all_s_binding_inds <- mapply(seq,
                             from = s_binding_sites$start_pos, 
                             to = s_binding_sites$stop_pos) %>% 
  unlist()

cluster_binding <- feature_clusters %>% 
  group_by(.data$cluster) %>% 
  summarise(includes_binding_site = any(.data$feature_position_corrected %in% all_s_binding_inds))

feature_locations <- feature_clusters %>% 
  left_join(cluster_binding, by = "cluster") %>% 
  mutate(s_binding = case_when(.data$feature_position_corrected %in% s_binding_sites ~ "S-binding",
                               .data$includes_binding_site ~ "Correlated with\nS-binding site",
                               TRUE ~ "Other"))


# ---- Feature importance --------------------------------------------------------------------------
n_models <- feature_importance %>% 
  group_by(.data$feature) %>% 
  summarise(n_models = sum(.data$importance > 0),
            .groups = "drop") %>% 
  mutate(n_model_rank = rank(-.data$n_models)) %>% 
  filter(.data$n_model_rank < 20) %>% 
  add_readable_feature_names() %>% 
  arrange(-.data$n_model_rank) %>% 
  mutate(feature_label = factor(.data$feature_label, levels = .data$feature_label)) %>% 
  left_join(feature_locations, by = "feature_position_corrected")

top_importance <- feature_importance %>% 
  filter(.data$feature %in% n_models$feature) %>% 
  add_readable_feature_names() %>% 
  mutate(feature_label = factor(.data$feature_label, levels = levels(n_models$feature_label))) %>% 
  left_join(feature_locations, by = "feature_position_corrected")


p_n_models <- ggplot(n_models, aes(x = feature_label, y = n_models, fill = s_binding)) +
  geom_col(colour = "grey20") +
  coord_flip() +
  scale_fill_brewer(palette = "Set1", na.value = "grey60", guide = FALSE) +
  labs(y = "Number of replicates", x = NULL, fill = "S-interaction") +
  PLOT_THEME +
  theme(plot.margin = margin(t = 5.5, r = 2.5, b = 5.5, l = 5.5))

p_overall_importance <- ggplot(top_importance, aes(x = feature_label, y = importance, fill = s_binding)) +
  geom_boxplot(colour = "grey20") +
  coord_flip() +
  scale_fill_brewer(palette = "Set1", na.value = "grey60", guide = FALSE) +
  labs(y = "Effect magnitude\n(all replicates)", x = NULL, fill = "S-interaction") +
  PLOT_THEME +
  theme(axis.text.y = element_blank(),
        plot.margin = margin(t = 5.5, r = 2.5, b = 5.5, l = 0))

p_present_importance <- ggplot(top_importance[top_importance$importance > 0, ], 
                               aes(x = feature_label, y = importance, colour = s_binding)) +
  geom_boxplot() +
  coord_flip() +
  scale_colour_brewer(palette = "Set1", na.value = "grey60") +
  labs(y = "Effect magnitude\n(when present)", x = NULL, colour = "S-interaction") +
  PLOT_THEME +
  theme(axis.text.y = element_blank(),
        plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 0))


p_importance <- plot_grid(p_n_models, p_overall_importance, p_present_importance, 
                          ncol = 3, rel_widths = c(1.4, 1, 1.4),
                          align = "hv", axis = "tb")


# ---- Importance by site --------------------------------------------------------------------------
site_importance <- feature_importance %>% 
  add_readable_feature_names() %>% 
  filter(!is.na(.data$feature_position_corrected)) %>% 
  group_by(.data$iteration, .data$feature_position_corrected) %>% 
  summarise(total_importance = sum(.data$importance),
            .groups = "drop") %>% 
  group_by(.data$feature_position_corrected) %>% 
  summarise(n_models = sum(.data$total_importance > 0),
            mean_importance = mean(.data$total_importance),
            .groups = "drop") %>% 
  left_join(feature_locations, by = "feature_position_corrected")


p_site_count <- ggplot(site_importance) +
  geom_rect(aes(xmin = start_pos, xmax = stop_pos, ymin = -Inf, ymax = Inf), 
            fill = "grey70", data = s_binding_sites) +
  geom_point(aes(x = feature_position_corrected, y = n_models, colour = s_binding)) +
  scale_fill_brewer(palette = "Set1", na.value = "grey60", guide = FALSE) +
  labs(x = NULL, y = "Number of replicates", colour = "S-interaction") +
  PLOT_THEME +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 5.5, r = 5.5, b = 0, l = 5.5))

p_site_importance <- ggplot(site_importance) +
  geom_rect(aes(xmin = start_pos, xmax = stop_pos, ymin = -Inf, ymax = Inf), 
            fill = "grey70", data = s_binding_sites) +
  geom_point(aes(x = feature_position_corrected, y = mean_importance, colour = s_binding)) +
  labs(x = "Sequence position", y = "Mean effect magnitude", colour = "S-interaction") +
  PLOT_THEME +
  theme(plot.margin = margin(t = 0, r = 5.5, b = 5.5, l = 5.5))


p_sites <- plot_grid(p_site_count, p_site_importance, 
                     ncol = 1, align = "v", axis = "lr")


# ---- Combine-------------------------------------------------------------------------------------
p_overview <- plot_grid(p_importance, 
                        p_sites,
                        ncol = 1, rel_heights = c(1, 1), 
                        labels = c("A", "B"))

ggsave2(INPUT$output_name, p_overview, width = 6, height = 6)