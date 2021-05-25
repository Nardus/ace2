## Plot observed data on a phylogeny

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ape)
  library(ggplot2)
  library(ggtree)
  library(cowplot)
})

timetree <- read.tree("data/internal/training_spp_timetree.nwk")

infection_data <- readRDS("data/calculated/cleaned_infection_data.rds")
shedding_data <- readRDS("data/calculated/cleaned_shedding_data.rds")



# ---- Fix species names ---------------------------------------------------------------------------
timetree$tip.label <- str_replace(timetree$tip.label, "_", " ")

name_replacements <- c("Canis lupus" = "Canis familiaris",
                       "Neophocaena phocaenoides" = "Neophocaena asiaeorientalis",
                       "Monachus schauinslandi" = "Neomonachus schauinslandi")

stopifnot(all(names(name_replacements) %in% timetree$tip.label))

timetree$tip.label <- with(timetree,
                           if_else(tip.label %in% names(name_replacements),
                                   as.character(name_replacements[tip.label]),
                                   tip.label))

stopifnot(all(timetree$tip.label %in% infection_data$species))
stopifnot(all(infection_data$species %in% timetree$tip.label))


# ---- Plot ----------------------------------------------------------------------------------------
# Tree
max_date <- max(cophenetic.phylo(timetree)) / 2
tree_breaks <- seq(max_date %% 100, max_date, by = 100)
tree_labels <- max_date - tree_breaks # reverse axis, showing years since present

tp <- ggtree(timetree) +
  geom_tiplab(size = 2.5, fontface = "italic") +
  scale_x_ggtree(breaks = tree_breaks, labels = tree_labels) +
  theme_tree2()


# Data
evidence_labels <- c("1" = "Observed infection",
                     "2" = "Experimental infection",
                     "3" = "Cell culture",
                     "4" = "Cell culture (het-ACE2)")
 
infection_matrix <- infection_data %>% 
  select(.data$species, .data$infected, .data$evidence_level) %>% 
  mutate(infected = if_else(.data$infected == "True", "Infected", "Not infected"),
         evidence_level = factor(.data$evidence_level)) %>% 
  pivot_wider(names_from = .data$infected, values_from = .data$evidence_level)

shedding_matrix <- shedding_data %>% 
  select(.data$species, .data$shedding, .data$evidence_level) %>% 
  mutate(shedding = if_else(.data$shedding == "True", "Shedding", "Not shedding"),
         evidence_level = factor(.data$evidence_level)) %>% 
  pivot_wider(names_from = .data$shedding, values_from = .data$evidence_level)


data_matrix <- left_join(infection_matrix, shedding_matrix, by = "species") %>% 
  as.data.frame()

rownames(data_matrix) <- data_matrix$species
data_matrix <- data_matrix[, -1]


# Plot
p_final <- gheatmap(tp, data_matrix, offset = 240, width = 0.3, 
                    font.size = 3,
                    colnames_position = "top", colnames_angle = 90, hjust = 0) +
  geom_vline(xintercept = 610) + 
  scale_fill_brewer(palette = "YlGnBu", direction = - 1,
                    labels = c(evidence_labels)) +
  scale_y_continuous(expand = expand_scale(add = c(1, 10)))


# ---- Output --------------------------------------------------------------------------------------
dir.create("output/plots", recursive = TRUE)
ggsave2("output/plots/raw_data_phylogeny.png", 
        p_final,
        width = 6, height = 7)
