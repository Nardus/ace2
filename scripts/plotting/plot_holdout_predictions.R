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
  library(ggtext)
  
  source("scripts/utils/plot_utils.R")
  source("scripts/plotting/plotting_constants.R")
})


# ---- Data ----------------------------------------------------------------------------------------
ensemble_predictions <- readRDS("output/all_data/infection/ensemble_all_features_phylogeny/holdout_predictions.rds")
phylogeny_predictions <- readRDS("output/all_data/infection/phylogeny/holdout_predictions.rds")

infection_data <- readRDS("data/calculated/cleaned_infection_data.rds")
taxonomy_table <- readRDS("data/calculated/taxonomy.rds")

time_tree <- read.tree("data/internal/timetree_amniota.nwk")


# ---- Prepare tree --------------------------------------------------------------------------------
# Change taxonomy to match tree
tree_taxonomy <- taxonomy_table %>% 
  filter(.data$internal_name %in% phylogeny_predictions$species) %>% 
  filter(.data$species != "Bos indicus x Bos taurus") %>% 
  filter(!.data$species %in% c("Tinamus guttatus", 
                               "Chiroxiphia lanceolata")) %>% # These don't appear in TimeTree
  filter(!is.na(.data$order))

tree_taxonomy <- tree_taxonomy %>% 
  mutate(internal_name = case_when(internal_name == "Antrostomus carolinensis" ~ "Caprimulgus carolinensis",
                                   internal_name == "Corvus cornix" ~ "Corvus corone",
                                   internal_name == "Strigops habroptila" ~ "Strigops habroptilus",
                                   internal_name == "Cebus imitator" ~ "Cebus capucinus",
                                   internal_name == "Dryobates pubescens" ~ "Picoides pubescens",
                                   internal_name == "Canis familiaris" ~ "Canis lupus",
                                   internal_name == "Neomonachus schauinslandi" ~ "Monachus schauinslandi",
                                   internal_name == "Neophocaena asiaeorientalis" ~ "Neophocaena phocaenoides",
                                   TRUE ~ internal_name))

# Reduce to relevant species
time_tree$tip.label <- str_replace_all(time_tree$tip.label, "_", " ")

stopifnot(all(tree_taxonomy$internal_name %in% time_tree$tip.label))
time_tree <- keep.tip(time_tree, tree_taxonomy$internal_name)


# Reduce phylogeny to order level:
for (ord in unique(tree_taxonomy$order)) {
  species <- tree_taxonomy %>% 
    filter(.data$order == ord) %>% 
    pull(.data$internal_name)
  
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


# ---- Merge data ----------------------------------------------------------------------------------
ensemble_predictions <- ensemble_predictions %>% 
  mutate(run_label = "ACE2 / phylogeny ensemble")

phylogeny_predictions <- phylogeny_predictions %>% 
  mutate(run_label = "Phylogeny-only")

infection_data <- infection_data %>% 
  mutate(run_label = "Observed",
         predicted_label = .data$infected)


all_predictions <- bind_rows(ensemble_predictions, phylogeny_predictions, infection_data) %>% 
  mutate(run_label = factor(.data$run_label, 
                            levels = c("Observed",
                                       "ACE2 / phylogeny ensemble", 
                                       "Phylogeny-only"))) %>% 
  left_join(taxonomy_table, by = c("species" = "internal_name")) %>% 
  filter(!is.na(.data$class))


# ---- Predictions by taxonomic order --------------------------------------------------------------
order_summary <- all_predictions %>%
  filter(!is.na(.data$order)) %>% 
  group_by(.data$run_label, .data$class, .data$order) %>%
  summarise(n_species = n(),
            n_susceptible = sum(.data$predicted_label == "True"),
            .groups = "keep") %>% 
  mutate(prop_susceptible = .data$n_susceptible / .data$n_species,
         lower = binom.test(.data$n_susceptible, .data$n_species)$conf.int[1],
         upper = binom.test(.data$n_susceptible, .data$n_species)$conf.int[2]) %>% 
  ungroup() %>% 
  mutate(order = factor(.data$order, levels = label_order))


data_panel <- ggplot(order_summary, aes(x = order, y = prop_susceptible, fill  = class)) +
  geom_col() +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.4,
                color = "grey20") +
  
  scale_x_discrete(position = "top") +
  scale_y_continuous(expand = expansion(add = c(0, 0.05))) +
  scale_fill_brewer(palette = "Dark2") +
  
  facet_grid(cols = vars(run_label), 
             scales = "free", space = "free", switch = "y") +
  
  labs(x = "Order", y = "Proportion susceptible", fill = "Class") +
  coord_flip() +
  theme(axis.title.y = element_blank(),
        plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 0),
        panel.grid.major.y = element_line(colour = "grey90"))


# ---- Counts vs expected -------------------------------------------------------------------------
count_summary <- order_summary %>% 
  filter(.data$run_label != "ACE2 / phylogeny ensemble")

count_labels <- count_summary %>% 
  arrange(.data$order) %>% 
  group_by(.data$run_label, .data$n_species, .data$n_susceptible) %>% 
  filter(n() < 4) %>% 
  summarise(order = paste(.data$order, collapse = ", "),
            .groups = "drop") %>% 
  mutate(angle = if_else(.data$n_susceptible == 0 | .data$n_species < 9, 45, 0))


count_plot <- ggplot(count_summary, aes(x = log10(n_species), y = log10(n_susceptible), colour = class)) +
  geom_abline(linetype = 2, colour = "grey80") +
  geom_point() +
  geom_text(aes(label = order, angle = angle), 
            data = count_labels,
            colour = "grey10",
            size = 2,
            hjust = 0,
            nudge_x = log10(1.1),
            nudge_y = log10(1.05)) +
  
  scale_colour_brewer(palette = "Dark2", guide = "none") +
  
  facet_grid(cols = vars(run_label)) +
  labs(x = "log<sub>10</sub>(Number of species)", 
       y = "log<sub>10</sub>(Number of susceptible species)") +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_markdown())


# ---- Output --------------------------------------------------------------------------------------
tree_plot <- plot_grid(tree_panel, data_panel,
                       nrow = 1, rel_widths = c(1, 3, 3),
                       align = "h", axis = "tb")

final_plot <- plot_grid(tree_plot, count_plot,
                        nrow = 2, rel_heights = c(1.2, 1),
                        labels = c("A", "B"))

ggsave2("output/plots/predictions_by_order_supplement.pdf", final_plot,
        width = 8, height = 9)
