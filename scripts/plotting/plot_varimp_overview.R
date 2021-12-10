## Plot overview of feature usage

suppressPackageStartupMessages({
  library(argparse)
})

parser <- ArgumentParser(description = "Plot an overview of feature usage")

parser$add_argument("varimp_file", type = "character",
                    help = "location of variable importance measurements")

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
  library(seqinr)
  library(vegan)
  library(ggbeeswarm)
  library(ggsignif)
  
  source("scripts/utils/plot_utils.R")
  source("scripts/plotting/plotting_constants.R")
})

N_CORES <- 20

feature_importance <- readRDS(INPUT$varimp_file)

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
  mutate(feature_label = factor(.data$feature_label, levels = .data$feature_label)) %>% 
  left_join(feature_locations, by = "feature_position")

# Add cluster info
# - Renaming clusters by importance
cluster_labels <- top_importance %>% 
  group_by(.data$cluster) %>%
  summarise(cluster_importance = sum(.data$importance), .groups = "drop") %>% 
  arrange(.data$cluster_importance) %>% 
  mutate(cluster_label = rank(-.data$cluster_importance, ties.method = "random"),
         cluster_label = factor(.data$cluster_label, levels = .data$cluster_label)) %>% 
  select(-.data$cluster_importance)

top_importance <- top_importance %>% 
  left_join(cluster_labels, by = "cluster")

# Plot
p_overall_importance <- ggplot(top_importance, aes(x = feature_label, y = importance, 
                                                   fill = feature_type)) +
  geom_col(colour = "grey20", size = 0.2) +
  coord_flip() +
  scale_fill_brewer(palette = "Set2", na.value = "grey60", guide = "none") +
  labs(y = "Effect magnitude", x = "Selected feature", fill = "Measure")


# ---- Importance by site --------------------------------------------------------------------------
site_labels <- top_importance %>% 
  group_by(.data$feature_position, .data$feature_position_corrected) %>% 
  summarise(total_importance = sum(.data$importance),
            .groups = "drop") %>% 
  arrange(.data$total_importance) %>% 
  mutate(site_label = factor(.data$feature_position_corrected,
                             levels = .data$feature_position_corrected),
         site_label_continuous = as.numeric(.data$site_label))

site_importance <- top_importance %>% 
  left_join(site_labels, by = c("feature_position", "feature_position_corrected"))

p_site_importance <- ggplot(site_importance, aes(x = site_label_continuous, y = importance, 
                                                 fill = feature_type)) +
  geom_col(colour = "grey20", size = 0.2, position = "stack") +
  coord_flip() +
  scale_x_continuous(breaks = site_importance$site_label_continuous,
                     labels = site_importance$site_label,
                     sec.axis = sec_axis(~ ., name = "Correlated cluster",
                                         breaks = site_importance$site_label_continuous,
                                         labels = site_importance$cluster_label)) +
  scale_fill_brewer(palette = "Set2", na.value = "grey60") +
  labs(y = "Total effect magnitude", x = "Sequence position\n(human ACE2)", fill = "Measure")



# ---- Cluster overview ----------------------------------------------------------------------------
# Need to plot these in terms of alignment coordinates, since some positions do not occur in the
# human sequence (gaps)
cluster_overview <- feature_locations %>% 
  right_join(cluster_labels, by = "cluster")

human_coords <- as_human_coord_v(1:max(feature_clusters$feature_position))
location_adjustments <- 1:max(feature_clusters$feature_position)
names(location_adjustments) <- human_coords

adjusted_s_binding_sites <- S_BINDING_SITES
adjusted_s_binding_sites$start_pos <- location_adjustments[as.character(S_BINDING_SITES$start_pos)]
adjusted_s_binding_sites$stop_pos <- location_adjustments[as.character(S_BINDING_SITES$stop_pos)]

p_cluster_overview <- ggplot(cluster_overview) +
  geom_rect(aes(xmin = start_pos, xmax = stop_pos, ymin = -Inf, ymax = Inf), 
            fill = "grey70", data = adjusted_s_binding_sites) +
  geom_tile(aes(x = feature_position, y = cluster_label,  fill = s_binding)) + 
  
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Sequence position (human ACE2)",
                                         breaks = site_labels$feature_position,
                                         labels = site_labels$feature_position_corrected),
                     expand = expansion(mult = 0.002)) +
  scale_fill_brewer(palette = "Set1", na.value = "grey60", guide = "none") +
  labs(x = "Alignment position", y = "Cluster")


# Add combined importance of each cluster:
cluster_importance <- top_importance %>% 
  group_by(.data$cluster_label) %>% 
  summarise(total_shap = sum(.data$importance), .groups = "drop")

p_cluster_importance <- ggplot(cluster_importance, aes(x = cluster_label, y = total_shap)) +
  geom_col(fill = "grey30") +
  scale_y_continuous(expand = expansion(add = c(0, 0.01))) +
  labs(y = "Total effect\nmagnitude") +
  coord_flip() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = -4))

p_cluster_importance <- plot_grid(p_cluster_overview, p_cluster_importance,
                                  ncol = 2, rel_widths = c(12, 1),
                                  align = "h", axis = "tb")


# ---- Phylogenetic information in these sites -----------------------------------------------------
get_entropy <- function(alignment_matrix) {
  # Count amino acids at each site
  counts <- apply(alignment_matrix, MARGIN = 2, FUN = table, simplify = FALSE) %>% 
    sapply(data.frame, simplify = FALSE) %>% 
    bind_rows(.id = "position") %>% 
    pivot_wider(id_cols = "position", names_from = "Var1", values_from = "Freq", values_fill = 0)
  
  counts_mat <- counts %>% 
    select(-.data$position) %>% 
    data.frame()
  
  rownames(counts_mat) <- counts$position
  
  # Get Shannon entropy at each site
  entropy <- diversity(counts_mat, index = "shannon")
  
  tibble(position = as.numeric(names(entropy)),
         entropy = entropy)
}

ace2_alignment <- as.matrix(ace2_alignment)
training_alignment <- ace2_alignment[unique(infection_data$ace2_accession), ]

entropy_training <- get_entropy(training_alignment)

entropy_training <- entropy_training %>% 
  left_join(feature_locations, by = c("position" = "feature_position")) %>% 
  mutate(selected = if_else(.data$position %in% site_labels$feature_position, "Yes", "No"),
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
  scale_colour_brewer(palette = "Set1", na.value = "grey60") +
  labs(x = "Position retained", y = "Entropy", colour = "S-interaction")


# ---- Combine--------------------------------------------------------------------------------------
p_top <- plot_grid(p_overall_importance, p_site_importance, p_entropy,
                   ncol = 3, rel_widths = c(1, 1.35, 1.5),
                   labels = c("A", "B", "C"))

p_combined <- plot_grid(p_top, p_cluster_importance,
                        ncol = 1, rel_heights = c(1, 0.8),
                        labels = c("", "D"))

ggsave2(INPUT$output_name, p_combined, width = 7, height = 4)


# ---- Values mentioned in text --------------------------------------------------------------------
cat("\n\nSequence positions given with reference to accesion:", human_acc, "\n\n")

# S-interacting sites
s_clusters <- feature_clusters %>% 
  left_join(cluster_binding, by = "cluster") %>% 
  filter(!.data$feature_position_corrected %in% ALL_S_BINDING_INDS) %>% 
  mutate(selected_site = .data$feature_position %in% top_importance$feature_position)

stopifnot(n_distinct(s_clusters$feature_position) == nrow(s_clusters))

# Expected ratio
exp_ratio <- sum(s_clusters$includes_binding_site) / nrow(s_clusters)
sprintf("%3.1f%% of positions considered are correlated with an S-interacting site\n",
        exp_ratio * 100) %>% 
  cat()

# Observed ratio
observed <- sum(s_clusters$selected_site & s_clusters$includes_binding_site)
available <- sum(s_clusters$selected_site)
test_result <- binom.test(x = observed,
                          n = available,
                          p = exp_ratio)

sprintf("%i of %i selected sites clustered with S-binding sites (%3.1f%%, p-value compared to expected ratio = %3.3f)",
        observed,
        available,
        test_result$estimate * 100,
        test_result$p.value) %>% 
  cat()
