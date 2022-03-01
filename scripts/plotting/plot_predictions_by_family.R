## Plot predictions by family (zoomed to Boreoeutheria only)
#   - Supplement / follow up to "plot_predictions_by_order.R


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


# Orders in the Boreoeutheria clade:
plot_orders <- c("Rodentia", "Lagomorpha", "Scandentia", "Dermoptera", "Primates",
                 "Eulipotyphla", "Chiroptera", "Cetartiodactyla", "Perissodactyla",
                 "Pholidota", "Carnivora")




# ---- Prepare taxonomy ----------------------------------------------------------------------------
# Fix polyphyletic clades, matching tree in plot_predictions_by_order.R
taxonomy_table <- taxonomy_table %>% 
  mutate(order = case_when(family %in% c("Erinaceidae", "Soricidae", "Talpidae", "Solenodontidae") ~ "Eulipotyphla",
                           order %in% c("Cetacea", "Artiodactyla") ~ "Cetartiodactyla",
                           family %in% c("Scopidae", "Balaenicipitidae") ~ "Ciconiiformes",  # According to NCBI taxonomy
                           order %in% c("Suliformes", "Pelecaniformes") ~ "Pelecaniformes/Suliformes",
                           family == "Cathartidae" ~ "Ciconiiformes",
                           family == "Aegothelidae" ~ "Aegotheliformes", 
                           order %in% c("Leptosomiformes", "Coraciiformes") ~ "Coraciiformes/Leptosomiformes",  
                           TRUE ~ order))

stopifnot(all(plot_orders %in% taxonomy_table$order))  # Check spelling


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
  filter(!is.na(.data$class)) %>% 
  filter(!is.na(.data$order)) %>% 
  filter(.data$order %in% plot_orders)


# ---- Predictions by family--------------------------------------------------------------
family_summary <- all_predictions %>%
  group_by(.data$run_label, .data$order, .data$family) %>%
  summarise(n_species = n(),
            n_susceptible = sum(.data$predicted_label == "True"),
            .groups = "keep") %>% 
  mutate(prop_susceptible = .data$n_susceptible / .data$n_species,
         lower = binom.test(.data$n_susceptible, .data$n_species)$conf.int[1],
         upper = binom.test(.data$n_susceptible, .data$n_species)$conf.int[2]) %>% 
  ungroup()


p <- ggplot(family_summary, aes(x = family, y = prop_susceptible)) +
  geom_col(fill = "grey70") +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.4,
                color = "grey40") +
  
  scale_x_discrete(position = "top") +
  scale_y_continuous(expand = expansion(add = c(0, 0.05))) +
  
  facet_grid(rows = vars(order), cols = vars(run_label), 
             scales = "free", space = "free", switch = "y") +
  
  labs(x = "Family", y = "Proportion susceptible") +
  coord_flip() +
  theme(axis.title.y = element_blank(),
        plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 0),
        panel.grid.major.y = element_line(colour = "grey90"),
        strip.text.y.left = element_text(angle = 0))


# ---- Output --------------------------------------------------------------------------------------
ggsave2("output/plots/predictions_by_family_supplement.pdf", p,
        width = 8, height = 9)

