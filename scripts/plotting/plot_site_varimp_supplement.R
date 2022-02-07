## Plot variable importance for the best site-specific model
#   - Supplement to "plot_varimp_overview.R"

# ---- Data ----------------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(cowplot)
  library(seqinr)
  library(vegan)
  library(ggbeeswarm)
  library(ggsignif)
  
  source("scripts/utils/plot_utils.R")
  source("scripts/plotting/plotting_constants.R")
})

feature_importance <- readRDS("output/all_data/infection/aa_distance/feature_importance.rds")

feature_clusters <- readRDS("output/plots/intermediates/feature_clusters.rds")

ace2_alignment <- read.alignment("data/calculated/ace2_protein_alignment.fasta",
                                 format = "fasta")

infection_data <- readRDS("data/calculated/cleaned_infection_data.rds")


# ---- Mark sites known to interact with S ---------------------------------------------------------
# Sites defined in plotting_constants.R
cluster_binding <- feature_clusters %>% 
  group_by(.data$cluster) %>% 
  summarise(includes_binding_site = any(.data$feature_position_corrected %in% ALL_S_BINDING_INDS))

feature_locations <- feature_clusters %>% 
  left_join(cluster_binding, by = "cluster") %>% 
  mutate(s_binding = case_when(.data$feature_position_corrected %in% ALL_S_BINDING_INDS ~ "S-binding",
                               .data$includes_binding_site ~ "Correlated with\nS-binding site",
                               TRUE ~ "Other"),
         s_binding = factor(.data$s_binding, levels = c("S-binding", 
                                                        "Correlated with\nS-binding site",
                                                        "Other"))) %>% 
  select(-.data$feature_position_corrected)


# ---- Feature importance --------------------------------------------------------------------------
top_importance <- feature_importance %>% 
  filter(.data$importance > 0) %>% 
  add_readable_feature_names() %>% 
  arrange(.data$importance) %>% 
  mutate(feature_label = factor(.data$feature_position_corrected, 
                                levels = .data$feature_position_corrected),
         feature_type = factor(.data$feature_type)) %>% 
  left_join(feature_locations, by = "feature_position")

# Plot
p_overall_importance <- ggplot(top_importance, aes(x = feature_label, y = importance, 
                                                   fill = s_binding)) +
  geom_col(colour = "grey20", size = 0.2) +
  scale_fill_brewer(palette = "Set1", na.value = "grey60", drop = FALSE) +
  coord_flip() +
  labs(y = "Effect magnitude", x = "Sequence position\n(human ACE2)", fill = "S-interaction")


# ---- Phylogenetic information in these sites -----------------------------------------------------
ace2_alignment <- as.matrix(ace2_alignment)
training_alignment <- ace2_alignment[unique(infection_data$ace2_accession), ]

entropy_training <- get_entropy(training_alignment)

entropy_training <- entropy_training %>% 
  left_join(feature_locations, by = c("position" = "feature_position")) %>% 
  mutate(selected = if_else(.data$position %in% top_importance$feature_position, "Yes", "No"),
         selected = factor(.data$selected, levels = c("Yes", "No"))) 


# Remove sites never considered in the model (sites removed by the near zero variance recipe filter)
included_sites <- feature_importance %>% 
  add_readable_feature_names() %>% 
  pull(.data$feature_position) %>% 
  unique()

entropy_training <- entropy_training %>% 
  filter(.data$position %in% included_sites)

# Plot
p_entropy <- ggplot(entropy_training, aes(x = selected, y = entropy)) +
  geom_quasirandom(aes(colour = s_binding), size = 0.4) +
  geom_boxplot(fill = NA, outlier.colour = NA, colour = "grey20") +
  
  geom_signif(test = "wilcox.test", comparisons = list(c("Yes", "No")),
              size = 0.3, tip_length = 0.014, textsize = 1.8, vjust = -0.2) +
  
  scale_y_continuous(expand = expansion(add = c(0.1, 0.15))) +
  scale_colour_brewer(palette = "Set1", na.value = "grey60", guide = "none") +
  labs(x = "Position retained", 
       y = "Phylogenetic informativeness\n(Shannon entropy)", 
       colour = "S-interaction")


# ---- Combine--------------------------------------------------------------------------------------
shared_legend <- get_legend(p_overall_importance)

p_overall_importance <- p_overall_importance +
  guides(fill = "none")

p_combined <- plot_grid(p_overall_importance, p_entropy, shared_legend,
                        ncol = 3, rel_widths = c(1, 1, 0.3),
                        labels = c("A", "B", ""))

ggsave2("output/plots/site_varimp_supplement.pdf", p_combined, width = 7, height = 5)
