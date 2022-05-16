## Plot entropy and S-binding status of alignment sites used by ACE2 models

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(cowplot)
  library(seqinr)
  library(vegan)
  library(vipor)
  library(ggsignif)
  
  source("scripts/utils/plot_utils.R")
  source("scripts/plotting/plotting_constants.R")
})

set.seed(2742022)


# ---- Data ----------------------------------------------------------------------------------------
feature_clusters <- readRDS("output/plots/intermediates/feature_clusters.rds")

ace2_alignment <- read.alignment("data/calculated/ace2_protein_alignment.fasta",
                                 format = "fasta")

infection_data <- readRDS("data/calculated/cleaned_infection_data.rds")


varimp_full_model <- readRDS("output/all_data/infection/all_features/feature_importance.rds")
varimp_dist_model <- readRDS("output/all_data/infection/aa_distance/feature_importance.rds")


# ---- Phylogenetic information in each site -------------------------------------------------------
ace2_alignment <- as.matrix(ace2_alignment)
training_alignment <- ace2_alignment[unique(infection_data$ace2_accession), ]

entropy_training <- get_entropy(training_alignment)


# ---- Mark sites known to interact with S ---------------------------------------------------------
# Sites defined in plotting_constants.R
cluster_binding <- feature_clusters %>% 
  group_by(.data$cluster) %>% 
  summarise(includes_binding_site = any(.data$feature_position_corrected %in% ALL_S_BINDING_INDS))

position_data <- feature_clusters %>%
  left_join(cluster_binding, by = "cluster") %>%
  mutate(s_binding = case_when(.data$feature_position_corrected %in% ALL_S_BINDING_INDS ~ "S-binding",
                               .data$includes_binding_site ~ "Correlated with\nS-binding site",
                               TRUE ~ "Other"))


# ---- Plot ----------------------------------------------------------------------------------------
prepare_data <- function(varimps, model_label, 
                         entropy = entropy_training, 
                         pos_data = position_data) {
  varimps <- add_readable_feature_names(varimps)
  
  # Remove sites never considered in the model (sites removed by the near zero variance recipe filter)
  included_sites <- varimps %>%
    pull(.data$feature_position) %>%
    unique()
  
  # Mark selected sites
  varimps <- varimps %>%
    filter(.data$importance > 0)

  entropy <- entropy %>%
    mutate(selected = if_else(.data$position %in% varimps$feature_position, "Yes", "No"),
           selected = factor(.data$selected, levels = c("Yes", "No"))) %>% 
    filter(.data$position %in% included_sites) %>%
    left_join(pos_data, by = c("position" = "feature_position")) %>%
    mutate(s_binding = if_else(is.na(.data$s_binding), "Other", .data$s_binding),
           s_binding = factor(.data$s_binding, levels = c("S-binding", 
                                                          "Correlated with\nS-binding site",
                                                          "Other")),
           position_label = if_else(.data$selected == "Yes", 
                                    .data$feature_position_corrected, 
                                    NA_integer_),
           model = model_label) %>% 
    mutate(x_offset = offsetX(.data$entropy, .data$selected), 
           adjusted_x = as.numeric(.data$selected) + .data$x_offset)
           
  entropy
}

entropy_full <- prepare_data(varimp_full_model, "All ACE2 representations combined")
entropy_dist <- prepare_data(varimp_dist_model, "AA consensus distance")

entropy_combined <- bind_rows(entropy_full, entropy_dist) %>%
  mutate(model = factor(.data$model,
                        levels = c("All ACE2 representations combined", "AA consensus distance")))


p <- ggplot(entropy_combined, aes(x = selected, y = entropy)) +
  geom_blank() +
  geom_point(aes(x = adjusted_x, colour = s_binding), size = 1) +
  geom_text(aes(x = adjusted_x, label = position_label), size = 1.8, colour = "grey40") +
  geom_boxplot(fill = NA, outlier.colour = NA, colour = "grey20") +
  geom_signif(test = "wilcox.test", comparisons = list(c("Yes", "No")),
              size = 0.3, tip_length = 0.014, textsize = 1.8, vjust = -0.2) +
  scale_y_continuous(expand = expansion(add = c(0.1, 0.15))) +
  scale_colour_brewer(palette = "Set1", na.value = "grey60") +
  labs(x = "Position retained",
       y = "Phylogenetic informativeness (Shannon entropy)",
       colour = "S-interaction") +
  facet_grid(cols = vars(model))


# ---- Output--------------------------------------------------------------------------------------
ggsave2("output/plots/varimp_overview.pdf", p, width = 7, height = 3)
