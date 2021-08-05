## Plot a comparison of predictions from previous studies alongside current data

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(readr)
  library(scales)
  library(cluster)
  library(ape)
  library(ggplot2)
  library(ggtree)
  library(vipor)
  library(cowplot)
})

source("scripts/plotting/plotting_constants.R")


timetree <- read.tree("data/internal/timetree_amniota.nwk")
taxonomy_table <- readRDS("data/calculated/taxonomy.rds")

infection_data <- readRDS("data/calculated/cleaned_infection_data.rds") %>% 
  mutate(infected = factor(.data$infected, levels = c("True", "False")))

shedding_data <- readRDS("data/calculated/cleaned_shedding_data.rds") %>% 
  mutate(shedding = factor(.data$shedding, levels = c("True", "False")),)

virus_data_infection <- readRDS("output/plots/intermediates/virus_data_infection.rds")

existing_predictions <- read_csv("data/internal/existing_predictions.csv",
                                 col_types = cols(raw_score = "d",
                                                  reverse_score = "l",
                                                  .default = "c"))


# ---- Filter to mammals and birds -----------------------------------------------------------------
# - Restrict to species for which ACE2 sequences are available (fischhoff2021a in particular 
#   contains many species which cannot be predicted in any other study)
# - Of these, use mammals and birds only  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< TODO
existing_predictions <- existing_predictions %>%
  filter(.data$species %in% c(taxonomy_table$species, taxonomy_table$subspecies)) %>% 
  #mutate(species = case_when(.data$species == "Canis lupus familiaris" ~ "Canis familiaris",
  #                           TRUE ~ .data$species),
  #       species = str_extract(.data$species, "^[[:alpha:]]+ [[:alpha:]]+")) %>%  # Remove subspecies info
  filter(.data$species %in% infection_data$species) %>% 
  group_by(.data$species, .data$citation_key, .data$prediction_type, .data$predictor, .data$reverse_score) %>% 
  summarise(raw_score = mean(.data$raw_score), .groups = "drop") # Average across subspecies


# ---- Make scores comparable ----------------------------------------------------------------------
# - All scores should be ascending, meaning a higher value = more likely to be positive
# - Rescale values from each study to lie in [0, 1], so we can use a single colour scale in the plot
existing_predictions <- existing_predictions %>% 
  group_by(.data$citation_key) %>% 
  mutate(scaled_score = if_else(.data$reverse_score, .data$raw_score * -1, .data$raw_score),
         scaled_score = rescale(.data$scaled_score, to = c(0, 1)))


# ---- Correlation between studies -----------------------------------------------------------------
score_mat <- existing_predictions %>% 
  select(.data$species, .data$citation_key, .data$scaled_score) %>% 
  pivot_wider(id_cols = "species", names_from = "citation_key", values_from = "scaled_score") %>% 
  column_to_rownames("species")

cor_mat <- cor(score_mat, use = "pairwise.complete", method = "spearman")

# Assume studies with no overlap have no correlation:
cor_mat <- replace(cor_mat, is.na(cor_mat), 0)

# Hierarchical clustering by correlation
prediction_clusters <- agnes(1 - cor_mat, diss = TRUE)
prediction_dendrogram <- as.phylo(as.hclust(prediction_clusters))


# ---- Prepare phylogeny ---------------------------------------------------------------------------
# Correct names
timetree$tip.label <- str_replace(timetree$tip.label, "_", " ")

name_replacements <- c("Canis lupus" = "Canis familiaris",
                       "Neophocaena phocaenoides" = "Neophocaena asiaeorientalis",
                       "Monachus schauinslandi" = "Neomonachus schauinslandi")

stopifnot(all(names(name_replacements) %in% timetree$tip.label))

timetree$tip.label <- with(timetree,
                           if_else(tip.label %in% names(name_replacements),
                                   as.character(name_replacements[tip.label]),
                                   tip.label))

stopifnot(all(existing_predictions$species %in% timetree$tip.label))

timetree <- keep.tip(timetree, unique(existing_predictions$species))


# ---- Plot main overview --------------------------------------------------------------------------
get_tree_order <- function(tree_plot) {
  tree_plot$data %>% 
    filter(.data$isTip) %>% 
    arrange(.data$y) %>% 
    pull(.data$label) %>% 
    na.omit() %>% 
    as.vector()
}

# Phylogeny - also show observed status here
annotation_data <- infection_data %>% 
  select(.data$species, .data$infected)

annotated_tree <- ggtree(timetree, aes(colour = infected)) %<+% annotation_data

phylo_panel <- annotated_tree +
  geom_rootedge(rootedge = 10) +
  scale_colour_manual(values = c(True = "#EE7733", False = "#009988"), 
                      na.value = "grey30",
                      guide = "none") +
  scale_x_continuous(expand = expansion(add = 0)) +
  scale_y_discrete(expand = expansion(add = 0.5)) +
  coord_cartesian(clip = "off") +
  theme_tree2() +
  theme(line = PLOT_THEME$line,
        axis.ticks.length = PLOT_THEME$axis.ticks.length,
        axis.text.x = element_text(size = 5.6),
        plot.margin = margin(t = 0, r = 0, b = 69.3, l = 5.5))

phylo_panel <- revts(phylo_panel)


# Study correlation dendrogram
dendro_panel <- ggtree(prediction_dendrogram) +
  scale_x_reverse(expand = expansion(add = c(0.002, 0.008))) +
  scale_y_discrete(expand = expansion(add = 0.2)) +
  coord_flip() +
  theme(line = PLOT_THEME$line,
        plot.margin = margin(t = 5.5, r = 92, b = 0, l = 0))


# Data panel
observed_status <- infection_data %>% 
  select(.data$species, .data$infected, .data$evidence_level)

existing_predictions <- existing_predictions %>% 
  left_join(observed_status, by = "species") %>% 
  mutate(species = factor(.data$species, levels = get_tree_order(phylo_panel)),
         citation_key = factor(.data$citation_key, levels = get_tree_order(dendro_panel)))

study_labels <- existing_predictions %>% 
  distinct(.data$citation_key, .data$prediction_type, .data$predictor) %>% 
  mutate(predictor = if_else(.data$predictor == "binding affinity change relative to humans",
                             "binding affinity change\nrelative to humans", .data$predictor),
         study_label = str_to_sentence(.data$citation_key),
         study_label = str_replace(.data$study_label, "([[:digit:]]{4})[[:alpha:]]*$", " et al., \\1")) %>% 
  str_glue_data("{study_label}\n({predictor})")

names(study_labels) <- unique(as.character(existing_predictions$citation_key))


data_panel <- ggplot(existing_predictions) +
  geom_hline(aes(yintercept = species, linetype = infected), colour = "grey92", size = 0.4) +
  geom_tile(aes(x = citation_key, y = species, fill = scaled_score)) +
  
  scale_linetype_manual(values = c("True" = 1, "False" = 2)) +
  scale_fill_viridis_c(direction = -1) +
  scale_x_discrete(labels = study_labels, expand = expansion(add = 0)) +
  scale_y_discrete(expand = expansion(add = 0.5), position="right") +
  
  labs(x = "Prediction", y = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(face = "italic"),
        legend.position="none",
        plot.margin = margin(t = 0, r = 5.5, b = 5.5, l = 0))



# Combine
shared_legend <- get_legend(
  ggplot(existing_predictions, aes(x = citation_key, y = scaled_score, fill = scaled_score,
                                   colour = infected, linetype = infected)) +
    geom_line() +
    geom_point(shape = 21, colour = "grey20") +
    scale_fill_viridis_c(direction = -1, guide = guide_colorbar(order = 1)) +
    scale_linetype_manual(values = c("True" = 1, "False" = 2), 
                          guide = guide_legend(order = 2, keywidth = unit(1.25, "lines"))) +
    scale_colour_manual(values = c(True = "#EE7733", False = "#009988"), 
                        na.value = "grey30", guide = guide_legend(order = 2)) +
    labs(fill = "Predicted score\n(scaled)", colour = "Infected", linetype = "Infected")
)

phylo_panel2 <- phylo_panel +
  annotation_custom(shared_legend, xmin = -500, ymin = 20)

overview_plot <- plot_grid(NULL, dendro_panel,
                           phylo_panel2, data_panel,
                           nrow = 2, rel_heights = c(1, 20), rel_widths = c(1, 5))
  

# ---- Plot summary panels -------------------------------------------------------------------------
evidence_labels <- c("1" = "Observed infection",
                     "2" = "Experimental infection",
                     "3" = "Cell culture",
                     "4" = "Cell culture (het-ACE2)")

jittered_predictions <- existing_predictions %>% 
  group_by(.data$citation_key) %>% 
  mutate(x_offset = offsetX(.data$scaled_score, .data$infected)) %>% 
  ungroup() %>% 
  mutate(adjusted_x = as.numeric(.data$infected) + .data$x_offset) %>% 
  left_join(virus_data_infection, by = "species")

summary_plot <- ggplot(existing_predictions, aes(x = infected, y = scaled_score)) +
  geom_boxplot(outlier.colour = NA) +
  geom_point(aes(x = adjusted_x, shape = virus, fill = factor(evidence_level)), colour = "grey40",
             data = jittered_predictions) +
  facet_wrap(vars(citation_key), labeller = as_labeller(study_labels),
             ncol = 3) +
  scale_fill_brewer(palette = "YlGnBu", direction = - 1, labels = evidence_labels, 
                    guide = guide_legend(order = 1, override.aes = list(shape = 22), ncol = 2)) +
  scale_shape_manual(values = VIRUS_SHAPES, guide = guide_legend(order = 2, ncol = 2)) +
  labs(x = "Infected", y = "Predicted score (scaled)", fill = "Best evidence", shape = "Virus") +
  theme(legend.position = "top",
        legend.box = "vertical",
        strip.text = element_text(size = 5, margin = margin(2.5, 2.5, 2.5, 2.5)))


# ---- Output --------------------------------------------------------------------------------------
final_figure <- plot_grid(overview_plot, summary_plot, 
                          nrow = 1, rel_widths = c(2, 1.15),
                          labels = c("A", "B"))

ggsave2("output/plots/existing_predictions.pdf", 
        final_figure, 
        width = 7, height = 7)
