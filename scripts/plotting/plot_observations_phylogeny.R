## Plot observed data on a phylogeny

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(ape)
  library(phytools)
  library(ggplot2)
  library(ggtree)
  library(ggbeeswarm)
  library(ggsignif)
  library(grid)
  library(cowplot)
})

source("scripts/plotting/plotting_constants.R")
source("scripts/utils/feature_calc_utils.R")


timetree <- read.tree("data/internal/timetree_amniota.nwk")

infection_data <- readRDS("data/calculated/cleaned_infection_data.rds")
shedding_data <- readRDS("data/calculated/cleaned_shedding_data.rds")


# Ace2 distance data:
ace2_dists <- readRDS("data/calculated/features_pairwise_dists.rds")

ncbi_metadata <- read.csv("data/internal/NCBI_ACE2_orthologs.csv") %>% 
  select(species = .data$Scientific.name, ace2_accession = .data$RefSeq.Protein.accessions)

internal_metadata <- read.csv("data/internal/ace2_accessions.csv") %>% 
  select(.data$species, .data$ace2_accession) %>% 
  filter(!is.na(.data$ace2_accession))

ace2_metadata <- bind_rows(ncbi_metadata, internal_metadata) %>% 
  distinct()


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

stopifnot(all(infection_data$species %in% timetree$tip.label))

# Remove species with no data
timetree <- keep.tip(timetree, infection_data$species)
stopifnot(all(timetree$tip.label %in% infection_data$species))


# ---- Plot data on phylogeny ----------------------------------------------------------------------
evidence_labels <- c("1" = "Observed infection",
                     "2" = "Experimental infection",
                     "3" = "Cell culture",
                     "4" = "Cell culture (het-ACE2)")

tree_midpoint <- max(cophenetic(timetree))/4

# Tree
tree_panel <- ggtree(timetree) +
  geom_rootedge(rootedge = 10) +
  scale_x_continuous(expand = expansion(add = 0)) +
  scale_y_discrete(expand = expansion(add = c(0.5, 0.5))) +
  coord_cartesian(clip = "off") +
  theme_tree2() +
  theme(line = PLOT_THEME$line,
        axis.ticks.length = PLOT_THEME$axis.ticks.length,
        axis.text.x = element_text(size = 5.6),
        panel.background = element_blank(),
        plot.margin = margin(t = 5.5, r = 0, b = -10, l = 5.5))

tree_panel <- revts(tree_panel)

label_order <- tree_panel$data %>% 
  filter(.data$isTip) %>% 
  arrange(.data$y) %>% 
  pull(.data$label) %>% 
  na.omit() %>% 
  as.vector()

# Data panels
infection_panel <- infection_data %>% 
  mutate(species = factor(.data$species, levels = label_order),
         infected = factor(.data$infected, levels = c("True", "False")),
         evidence_level = factor(.data$evidence_level,
                                 levels = c("1", "2", "3", "4"),
                                 labels = evidence_labels)) %>% 
  
  ggplot(aes(x = infected, y = species, fill = factor(evidence_level))) +
  geom_tile() +
  scale_fill_brewer(palette = "YlGnBu", direction = - 1) +
  scale_x_discrete(expand = expansion(add = 0)) +
  scale_y_discrete(expand = expansion(add = 0.5)) +
  labs(x = "Infected", fill = "Best evidence") +
  theme(legend.position="none",
        panel.grid.major.y = element_line(colour = "grey92"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(t = 5.5, r = 0, b = 5.5, l = 0))


shedding_panel <- shedding_data %>% 
  mutate(species = factor(.data$species, levels = label_order),
         shedding = factor(.data$shedding, levels = c("True", "False")),
         evidence_level = factor(.data$evidence_level,
                                 levels = c("1", "2", "3", "4"),
                                 labels = evidence_labels)) %>% 
  
  ggplot(aes(x = shedding, y = species, fill = evidence_level)) +
  geom_tile() +
  scale_fill_brewer(palette = "YlGnBu", direction = - 1, drop = FALSE) +
  scale_x_discrete(expand = expansion(add = 0)) +
  scale_y_discrete(expand = expansion(add = 0.5), drop = FALSE, position = "right") +
  labs(x = "Shedding", fill = "Best evidence") +
  theme(legend.position="none",
        panel.grid.major.y = element_line(colour = "grey92"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic"),
        plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 2))

# Combine
shared_legend <- get_legend(
  infection_panel +
    theme(legend.position = "left")
)

tree_panel2 <- tree_panel +
  annotation_custom(textGrob("Million years (My)", gp = gpar(fontsize = 7)), 
                    xmin = -tree_midpoint, xmax = -tree_midpoint,   # setting xmin = xmax centers label at xmin
                    ymin = -1.95, ymax = -1.95) +
  annotation_custom(shared_legend, xmin = -500, ymin = 20)

overview_plot <- plot_grid(tree_panel2, infection_panel, shedding_panel,
                           nrow = 1, rel_widths = c(6, 1, 3.41),
                           align = "h", axis = "tb")


# ---- Plot clustering vs distance -----------------------------------------------------------------
# Get ACE2 distances
human_accession <- ace2_metadata$ace2_accession[ace2_metadata$species == "Homo sapiens"]
rhinolophid_accessions <- ace2_metadata$ace2_accession[startsWith(ace2_metadata$species, "Rhinolophus")]

dist_to_human_ace2 <- ace2_dists %>% 
  filter(.data$other_seq == human_accession) %>% 
  left_join(ace2_metadata, by = "ace2_accession") %>% 
  select(.data$species, ace2_dist_human = .data$distance)

dist_to_rhinolophid_ace2 <- ace2_dists %>% 
  filter(.data$other_seq %in% rhinolophid_accessions) %>% 
  group_by(.data$ace2_accession) %>% 
  summarise(ace2_dist_rhinolophid = mean(.data$distance)) %>% 
  left_join(ace2_metadata, by = "ace2_accession") %>% 
  select(.data$species, .data$ace2_dist_rhinolophid)


# Get phylogenetic distances on timetree
# - Distances to closest positive/negative, as an indication of clustering
phylo_dists <- cophenetic(timetree) %>% 
  data.frame(check.names = FALSE) %>% 
  rownames_to_column("from_species") %>% 
  pivot_longer(-.data$from_species, names_to = "to_species", values_to = "distance")

phylo_dists <- phylo_dists %>% 
  rename(ace2_accession = .data$from_species,
         other_seq = .data$to_species)

dist_closest_positive <- infection_data %>% 
  select(-.data$ace2_accession) %>% 
  rename(label = .data$infected, ace2_accession = .data$species) %>% 
  get_dist_to_closest_positive(pairwise_dist_data = phylo_dists, 
                               metadata = .) %>% 
  select(species = .data$ace2_accession, .data$closest_positive_overall)

dist_closest_negative <- infection_data %>% 
  select(-.data$ace2_accession) %>% 
  rename(label = .data$infected, ace2_accession = .data$species) %>% 
  mutate(label = if_else(.data$label == "True", "False", "True")) %>%  # Reverse label so "closest_positive" returns closest negative
  get_dist_to_closest_positive(pairwise_dist_data = phylo_dists, 
                               metadata = .) %>% 
  select(species = .data$ace2_accession, 
         closest_negative_overall = .data$closest_positive_overall)

# Merge all distance/clustering measures
all_dists <- dist_to_human_ace2 %>% 
  full_join(dist_to_rhinolophid_ace2, by = "species") %>% 
  left_join(dist_closest_positive, by = "species") %>% 
  left_join(dist_closest_negative, by = "species")

dist_data <- infection_data %>% 
  left_join(all_dists, by = "species") %>% 
  mutate(infected = factor(.data$infected, levels = c("True", "False")))


# Plots
# - ACE distances
dist_h_plot <- ggplot(dist_data, aes(x = infected, y = ace2_dist_human)) +
  geom_boxplot(outlier.colour = NA) + 
  geom_quasirandom(aes(fill = factor(evidence_level)), shape = 21, size = 1.2, colour = "grey30") +
  
  geom_signif(comparisons = list(c("True", "False")), 
              test = "ks.test", test.args = list(exact = FALSE),
              size = 0.3, tip_length = 0.014, textsize = 1.8, vjust = -0.2, y_position = 20800) +
  
  scale_fill_brewer(palette = "YlGnBu", direction = - 1,
                    labels = c(evidence_labels)) +
  scale_y_continuous(limits = c(0, 23800), expand = expansion(add = 200)) +
  labs(x = "Infected", y = "ACE2 amino acid\ndistance to humans", fill = "Best evidence") + 
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(t = 5.5, r = 5.5, b = 0, l = 5.5))

dist_r_plot <- ggplot(dist_data, aes(x = infected, y = ace2_dist_rhinolophid)) +
  geom_boxplot(outlier.colour = NA) + 
  geom_quasirandom(aes(fill = factor(evidence_level)), shape = 21, size = 1.2, colour = "grey30") +
  
  geom_signif(comparisons = list(c("True", "False")), 
              test = "ks.test", test.args = list(exact = FALSE),
              size = 0.3, tip_length = 0.014, textsize = 1.8, vjust = -0.2, y_position = 20800) +
  
  scale_fill_brewer(palette = "YlGnBu", direction = - 1,
                    labels = c(evidence_labels)) +
  scale_y_continuous(limits = c(0, 23800), expand = expansion(add = 200)) +
  labs(x = "Infected", y = "Mean amino acid distance\nto Rhinolophid species", fill = "Best evidence") + 
  theme(legend.position = "none")


# - Phylogenetic clustering
clust_pos_plot <- ggplot(dist_data, aes(x = infected, y = closest_positive_overall)) +
  geom_boxplot(outlier.colour = NA) + 
  geom_quasirandom(aes(fill = factor(evidence_level)), shape = 21, size = 1.2, colour = "grey30") +
  
  geom_signif(comparisons = list(c("True", "False")), 
              test = "ks.test", test.args = list(exact = FALSE),
              size = 0.3, tip_length = 0.014, textsize = 1.8, vjust = -0.2, y_position = 630) +
  
  scale_fill_brewer(palette = "YlGnBu", direction = - 1,
                    labels = c(evidence_labels)) +
  scale_y_continuous(limits = c(0, 700), expand = expansion(add = 20)) +
  labs(x = "Infected", y = "Phylogenetic distance\nto closest infected\nneighbour (My)", 
       fill = "Best evidence") + 
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(t = 5.5, r = 5.5, b = 0, l = 5.5))

clust_neg_plot <- ggplot(dist_data, aes(x = infected, y = closest_negative_overall)) +
  geom_boxplot(outlier.colour = NA) + 
  geom_quasirandom(aes(fill = factor(evidence_level)), shape = 21, size = 1.2, colour = "grey30") +
  
  geom_signif(comparisons = list(c("True", "False")), 
              test = "ks.test", test.args = list(exact = FALSE),
              size = 0.3, tip_length = 0.014, textsize = 1.8, vjust = -0.2, y_position = 630) +
  
  scale_fill_brewer(palette = "YlGnBu", direction = - 1,
                    labels = c(evidence_labels)) +
  scale_y_continuous(limits = c(0, 700), expand = expansion(add = 20)) +
  labs(x = "Infected", y = "Phylogenetic distance\nto closest non-infected\nneighbour (My)", 
       fill = "Best evidence") + 
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(t = 5.5, r = 5.5, b = 0, l = 5.5))


# ---- Plot Pagel's lambda -------------------------------------------------------------------------
# Calculate clustering
get_lambda <- function(feature_name, data, tree) {
  feature <- setNames(data[, feature_name], data$species)
  pagel <- phylosig(tree, feature, method = "lambda", test = TRUE)
  
  data.frame(feature = feature_name,
             lambda = pagel$lambda,
             likelihood_ratio = -2 * (pagel$logL0 - pagel$logL),
             p_value = pagel$P)
}

# - For infection
dist_data_infection <- dist_data %>% 
  mutate(infected_l1 = if_else(.data$evidence_level > 1, FALSE, .data$infected == "True"),
         infected_l2 = if_else(.data$evidence_level > 2, FALSE, .data$infected == "True"),
         infected_l3 = if_else(.data$evidence_level > 3, FALSE, .data$infected == "True"),
         infected_l4 = .data$infected == "True") %>% 
  data.frame()

dist_tree_infection <- di2multi(timetree)  # Collapse 0-length branches

lambdas_infection <- lapply(paste0("infected_l", 1:4), get_lambda,
                            data = dist_data_infection,
                            tree = dist_tree_infection) %>% 
  bind_rows()


# - For shedding
stopifnot(all(shedding_data$evidence_level) %in% c(1, 2)) # Cell culture does not make sense here

dist_data_shedding <- shedding_data %>% 
  left_join(all_dists, by = "species") %>% 
  mutate(shedding_l1 = if_else(.data$evidence_level > 1, FALSE, .data$shedding == "True"),
         shedding_l2 = .data$shedding == "True") %>% 
  data.frame()

dist_tree_shedding <- keep.tip(timetree, shedding_data$species) %>% 
  di2multi()

lambdas_shedding <- lapply(paste0("shedding_l", 1:2), get_lambda,
                           data = dist_data_shedding,
                           tree = dist_tree_shedding) %>% 
  bind_rows()

# Plot
evidence_labels <- c("1" = "Observed\ninfection only",
                     "2" = "Experimental\ninfection added",
                     "3" = "Cell culture\nadded",
                     "4" = "Cell culture\n(het-ACE2)\nadded")

lambdas <- bind_rows(lambdas_infection, lambdas_shedding) %>% 
  mutate(label = str_extract(.data$feature, "^[[:alpha:]]+"),
         label = str_to_sentence(.data$label),
         max_evidence = str_extract(.data$feature, "[[:digit:]]$"),
         max_evidence = as.numeric(.data$max_evidence))

sample_sizes <- list(Infected = data.frame(table(max_evidence = infection_data$evidence_level)),
                     Shedding = data.frame(table(max_evidence = shedding_data$evidence_level))) %>% 
  bind_rows(.id = "label") %>% 
  group_by(.data$label) %>% 
  mutate(max_evidence = as.numeric(as.character(.data$max_evidence)),
         running_total = cumsum(.data$Freq),
         sample_size = sprintf("N = %d", .data$running_total)) %>% 
  left_join(lambdas, by = c("label", "max_evidence")) %>% 
  mutate(lambda = if_else(.data$lambda > 0.75, .data$lambda + 0.1, .data$lambda - 0.1)) # Nudge label


lambda_plot <- ggplot(lambdas, aes(x = max_evidence, y = lambda, group = label)) +
  geom_line(linetype = 3, colour = "grey70") +
  geom_point(aes(colour = label)) +
  geom_point(aes(shape = p_value < 0.01), colour = "grey20") + 
  
  geom_text(aes(label = sample_size), size = 1.5, colour = "grey20", data = sample_sizes) + 
  
  scale_x_continuous(breaks = 1:4, labels = evidence_labels, expand = expansion(add = 0.5)) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(add = 0.01)) +
  scale_shape_manual(values = c("TRUE" = 1, "FALSE" = 4),
                     labels = c("TRUE" = "True", "FALSE" = "False"),
                     guide = "none") +
  labs(x = "Evidence used", y = expression("Pagel's"~lambda),
       colour = NULL, shape = "p-value < 0.01") +
  theme(legend.position = c(0.75, 0.3),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


# ---- Output --------------------------------------------------------------------------------------
side_panels1 <- plot_grid(clust_pos_plot, clust_neg_plot, dist_h_plot, dist_r_plot, 
                          ncol = 1, rel_heights = c(1, 1, 1, 1.2),
                          align = "v", axis = "lr",
                          labels = c("C", "D", "E", "F"), label_size = 9, hjust = 0)

side_panels2 <- plot_grid(lambda_plot, side_panels1, 
                          ncol = 1, rel_heights = c(1, 3),
                          labels = c("B", ""), label_size = 9, hjust = 0)

final_figure <- plot_grid(overview_plot, side_panels2, 
                          ncol = 2, rel_widths = c(3, 1.2),
                          labels = c("A", ""), label_size = 9, hjust = -0.5)


dir.create("output/plots", recursive = TRUE)
ggsave2("output/plots/raw_data_overview.pdf", 
        final_figure,
        width = 7, height = 7)


# ---- Values mentioned in text --------------------------------------------------------------------
bird_spp <- c("Anser anser", "Anser cygnoides", "Anas platyrhynchos", "Gallus gallus", 
              "Meleagris gallopavo", "Coturnix japonica")

stopifnot(all(bird_spp %in% dist_data$species))

dists_no_birds <- dist_data %>% 
  filter(!.data$species %in% bird_spp)

ks_no_birds <- with(dists_no_birds, 
                    ks.test(closest_negative_overall[infected == "True"], 
                            closest_negative_overall[infected == "False"],
                            exact = FALSE))

cat("\n\nAfter excluding birds, negative cases are still significantly closer to other",
    "negative cases: D =",  ks_no_birds$statistic, "; p-value =", ks_no_birds$p.value)
