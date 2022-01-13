## Plot an overview of predictions for non-training species

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(ape)
  library(ggplot2)
  library(ggtree)
  library(cowplot)
  library(grid)
  
  source("scripts/utils/plot_utils.R")
  source("scripts/plotting/plotting_constants.R")
})


# ---- Data ----------------------------------------------------------------------------------------
all_predictions <- readRDS("output/all_data/infection/all_features/holdout_predictions.rds")
taxonomy_table <- readRDS("output/plots/intermediates/taxonomy_table.rds")

time_tree <- read.tree("data/internal/timetree_amniota.nwk")


# ---- Extract holdout -----------------------------------------------------------------------------
holdout_preds <- all_predictions %>% 
  rename(lookup_species = .data$species) %>% 
  left_join(taxonomy_table, by = c("lookup_species")) %>% 
  filter(.data$prediction_type != "Fitted value")


# ---- Prepare tree --------------------------------------------------------------------------------
# Change taxonomy to match tree
holdout_taxonomy <- taxonomy_table %>% 
  filter(.data$species %in% holdout_preds$species) %>% 
  filter(.data$species != "Bos indicus x Bos taurus") %>% 
  filter(!.data$species %in% c("Tinamus guttatus", 
                               "Chiroxiphia lanceolata", 
                               "Apteryx mantelli"))  # TODO: these should be in tree somewhere

holdout_taxonomy <- holdout_taxonomy %>% 
  mutate(species = case_when(species == "Antrostomus carolinensis" ~ "Caprimulgus carolinensis",
                             species == "Corvus cornix" ~ "Corvus corone",
                             species == "Strigops habroptila" ~ "Strigops habroptilus",
                             species == "Cebus imitator" ~ "Cebus capucinus",
                             species == "Dryobates pubescens" ~ "Picoides pubescens",
                             TRUE ~ species))

# Reduce to relevant species
time_tree$tip.label <- str_replace(time_tree$tip.label, "_", " ")

stopifnot(all(holdout_taxonomy$species %in% time_tree$tip.label))
time_tree <- keep.tip(time_tree, holdout_taxonomy$species)


# Reduce phylogeny to order level:
for (ord in unique(holdout_preds$order)) {
  species <- holdout_taxonomy %>% 
    filter(.data$order == ord) %>% 
    pull(.data$species)
  
  time_tree <- drop.tip(time_tree, species, subtree = TRUE)
  label_position <- grepl("\\[", time_tree$tip.label)  # Find ape's indicator of how many tips were dropped
  time_tree$tip.label[label_position] <- ord
}


# ---- Plot tree -----------------------------------------------------------------------------------
tree_panel <- ggtree(time_tree, right = TRUE) +
  geom_rootedge(rootedge = 10) +
  geom_tiplab(linetype = 1, align = TRUE, offset = ) +  # Compensate for overtrimming above
  scale_x_continuous(expand = expansion(add = 0)) +
  scale_y_discrete(expand = expansion(add = c(0.5, 0.5))) +
  coord_cartesian(clip = "off") +
  theme_tree() +
  theme(plot.margin = margin(t = 5.5, r = 0, b = 5.5, l = 5.5))


# Get order of labels for other panels
label_order <- tree_panel$data %>% 
  filter(.data$isTip) %>% 
  arrange(.data$y) %>% 
  pull(.data$label) %>% 
  na.omit() %>% 
  as.vector()


# ---- Predictions by taxonomic order --------------------------------------------------------------
order_importance <- holdout_preds %>% 
  group_by(.data$order) %>% 
  summarise(n_susceptible = sum(.data$predicted_label == "True"), 
            n_species = n(),
            .groups = "drop")

order_preds <- holdout_preds %>% 
  inner_join(order_importance, by = "order") %>% 
  arrange(.data$n_susceptible, .data$n_species) %>% 
  mutate(order = factor(.data$order, levels = label_order))

data_panel <- ggplot(order_preds, aes(x = order, fill = predicted_label)) +
  geom_bar(position = "stack") +
  scale_x_discrete(position = "top") +
  scale_y_continuous(expand = expansion(add = c(0, 1))) +
  scale_fill_manual(values = INFECTION_STATUS_COLOURS) +
  labs(x = NULL, y = "Number of species", fill = "Predicted\nsusceptible") +
  coord_flip() +
  PLOT_THEME +
  theme(axis.title.y = element_blank(),
        plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 0))


# ---- Output --------------------------------------------------------------------------------------
final_plot <- plot_grid(tree_panel, data_panel,
                        nrow = 1, rel_widths = c(1, 3),
                        align = "h", axis = "tb")

ggsave2("output/plots/holdout_predictions.png", final_plot,
        width = 4, height = 6)

# ---- Values mentioned in text --------------------------------------------------------------------
cat("Holdout predictions:\n")

holdout_preds %>% 
  group_by(.data$class) %>% 
  summarise(susceptible = sum(.data$predicted_label == "True"),
            n_species = n(),
            .groups = "drop") %>% 
  mutate(percent = .data$susceptible/.data$n_species * 100) %>% 
  print()