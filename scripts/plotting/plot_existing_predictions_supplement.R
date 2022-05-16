# Detailed overview of existing predictions (supplement to existing_predictions.pdf)

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(ape)
  library(ggplot2)
  library(ggtree)
  library(vipor)
  library(ggtext)
  library(cowplot)
  
  source("scripts/plotting/plotting_constants.R")
  source("scripts/utils/timetree_constants.R")
})

infection_data <- readRDS("data/calculated/cleaned_infection_data.rds") %>% 
  mutate(infected = factor(.data$infected, levels = c("True", "False")))

virus_data_infection <- readRDS("output/plots/intermediates/virus_data_infection.rds")

all_predictions <- readRDS("output/plots/intermediates/unified_predictions.rds")
prediction_dendrogram <- readRDS("output/plots/intermediates/prediction_dendrogram.rds")

timetree <- read.tree("data/internal/timetree_amniota.nwk")


# ---- Prepare phylogeny ---------------------------------------------------------------------------
# Correct names
timetree$tip.label <- str_replace(timetree$tip.label, "_", " ")

stopifnot(all(names(TIMETREE_TAXONOMY_CORRECTIONS) %in% timetree$tip.label))

timetree$tip.label <- with(timetree,
                           if_else(tip.label %in% names(TIMETREE_TAXONOMY_CORRECTIONS),
                                   as.character(TIMETREE_TAXONOMY_CORRECTIONS[tip.label]),
                                   tip.label))

# ---- Filter to species with known infection data -------------------------------------------------
all_predictions <- all_predictions %>%
  filter(.data$species %in% infection_data$species)

# Remove studies with <10 species remaining (to reduce clutter)
all_predictions <- all_predictions %>% 
  group_by(.data$citation_key) %>% 
  filter(n() >= 10) %>% 
  ungroup()

prediction_dendrogram <- keep.tip(prediction_dendrogram, unique(all_predictions$citation_key))

timetree <- keep.tip(timetree, unique(all_predictions$species))


# ---- Heatmap of all predictions ------------------------------------------------------------------
get_tree_order <- function(tree_plot) {
  tree_plot$data %>% 
    filter(.data$isTip) %>% 
    arrange(.data$y) %>% 
    pull(.data$label) %>% 
    na.omit() %>% 
    as.vector()
}

# Dendrogram / phylogeny annotations
dendro_panel <- ggtree(prediction_dendrogram, right = TRUE) +
  scale_x_reverse(expand = expansion(add = c(0.002, 0.008))) +
  scale_y_discrete(expand = expansion(add = c(0.35, 1.15))) +
  coord_flip() +
  theme(line = PLOT_THEME$line,
        plot.margin = margin(t = 5.5, r = 78, b = 0, l = 0))

# Phylogeny - also show observed status here
annotation_data <- infection_data %>% 
  select(.data$species, .data$infected)

annotated_tree <- ggtree(timetree, aes(colour = infected)) %<+% annotation_data

phylo_panel <- annotated_tree +
  geom_rootedge(rootedge = 10) +
  scale_colour_manual(values = INFECTION_STATUS_COLOURS, 
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

# Get data in correct order
timetree_preds <- all_predictions %>% 
  mutate(species = factor(.data$species, levels = get_tree_order(phylo_panel)),
         citation_key = factor(.data$citation_key, levels = get_tree_order(dendro_panel)))

study_label_df <- all_predictions %>% 
  distinct(.data$citation_key, .data$prediction_type, .data$predictor) %>% 
  mutate(predictor = if_else(.data$predictor == "binding affinity change relative to humans",
                             "binding affinity change<br/>relative to humans", .data$predictor),
         study_label = as.character(.data$citation_key),
         study_label = if_else(startsWith(.data$study_label, "This study"), "**This study**", 
                               .data$study_label),
         study_label = str_to_sentence(.data$study_label),
         study_label = str_replace(.data$study_label, "([[:digit:]]{4})[[:alpha:]]*$", " *et al.*, \\1"),
         study_label_long = paste0(.data$study_label, "<br/>(", .data$predictor, ")"),
         study_label_long = if_else(.data$study_label_long == "**This study**<br/>(All ACE2 representations)", 
                                     "**This study** (All ACE2<br/>representations)", .data$study_label_long))

study_labels <- study_label_df$study_label_long
names(study_labels) <- as.character(study_label_df$citation_key)

# Plot data
data_panel <- ggplot(timetree_preds) +
  geom_tile(aes(x = citation_key, y = species, fill = scaled_score)) +
  
  scale_fill_viridis_c(direction = -1) +
  scale_x_discrete(labels = study_labels, expand = expansion(add = 0)) +
  scale_y_discrete(expand = expansion(add = 0.5), position="right") +
  
  labs(x = "Prediction source", y = NULL) +
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(face = "italic"),
        legend.position="none",
        plot.margin = margin(t = 0, r = 5.5, b = 5.5, l = 0))


# Combine
timetree_preds <- timetree_preds %>% 
    left_join(annotation_data, by = c("species"))

shared_legend <- get_legend(
  ggplot(timetree_preds, aes(x = citation_key, y = scaled_score, 
                             fill = scaled_score, colour = infected)) +
    geom_line() +
    geom_point(shape = 21, colour = "grey20") +
    scale_fill_viridis_c(direction = -1, guide = guide_colorbar(order = 1)) +
    scale_linetype_manual(values = c("True" = 1, "False" = 2), 
                          guide = guide_legend(order = 2, keywidth = unit(1.25, "lines"))) +
    scale_colour_manual(values = c(True = "#EE7733", False = "#009988"), 
                        na.value = "grey30", guide = guide_legend(order = 2)) +
    labs(fill = "Scaled score", colour = "Infected") +
    theme(legend.background = element_blank())
)

phylo_panel2 <- phylo_panel +
  annotation_custom(shared_legend, xmin = -380, ymin = 80)

final_heatmap <- plot_grid(NULL, dendro_panel,
                           phylo_panel2, data_panel,
                           nrow = 2, rel_heights = c(1, 30), rel_widths = c(1, 5))


# ---- Plot summary panels -------------------------------------------------------------------------
# Filter to species with known infection data 
data_predictions <- infection_data %>% 
  left_join(all_predictions, by = "species")

# Mark model sample size in facet headings
fit_labels <- data_predictions %>% 
  group_by(.data$citation_key) %>% 
  summarise(n_species = n(),
            .groups = "drop") %>% 
  mutate(study_label = if_else(!startsWith(.data$citation_key, "This study"),
                               str_to_sentence(.data$citation_key),
                               .data$citation_key),
         study_label = str_replace(.data$study_label, "([[:digit:]]{4})[[:alpha:]]*$", " *et al.*, \\1"),
         study_label = str_replace(.data$study_label, "This study \\(([[:alnum:] ]*)\\)", "**\\1**"),
         study_label = str_replace(.data$study_label, "host phylogeny", "Host phylogeny"),
         study_label = str_replace(.data$study_label, "All ACE2 representations", "All ACE2<br>representations"),
         final_label = str_glue("{study_label}<br/>(N={n_species})"))

fit_label_vec = as.character(fit_labels$final_label)
names(fit_label_vec) <- as.character(fit_labels$citation_key)

# Jitter points to approximate local density
jittered_predictions <- data_predictions %>% 
  group_by(.data$citation_key) %>% 
  mutate(x_offset = offsetX(.data$scaled_score, .data$infected)) %>% 
  ungroup() %>% 
  mutate(adjusted_x = as.numeric(.data$infected) + .data$x_offset) %>% 
  left_join(virus_data_infection, by = "species")  # This duplicates points (so we can plot all viruses at a given species prediction)

# Plot
summary_plot <- ggplot(data_predictions, aes(x = infected, y = scaled_score)) +
  geom_boxplot(outlier.colour = NA, data = data_predictions) +  # Boxplot should not include the duplicates added for display 
  geom_point(aes(x = adjusted_x, shape = virus, fill = factor(evidence_level)), colour = "grey40",
             data = jittered_predictions) +
  facet_wrap(vars(citation_key), labeller = as_labeller(fit_label_vec),
             ncol = 3) +
  scale_fill_brewer(palette = "YlGnBu", direction = - 1, labels = EVIDENCE_LABELS, 
                    guide = guide_legend(order = 1, override.aes = list(shape = 22), ncol = 2)) +
  scale_shape_manual(values = VIRUS_SHAPES, guide = guide_legend(order = 2, ncol = 2)) +
  labs(x = "Infected", y = "Scaled score", fill = "Best\nevidence", shape = "Virus") +
  theme(legend.position = "top",
        legend.box = "vertical",
        legend.title = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
        strip.text = element_markdown(size = 5, margin = margin(2.5, 2.5, 2.5, 2.5)))


# ---- Combine -------------------------------------------------------------------------------------
final_plot <- plot_grid(final_heatmap, summary_plot,
                        rel_widths = c(2, 1.2),
                        labels = c("A", "B"))

ggsave2("output/plots/existing_predictions_supplement.pdf", 
        final_plot, 
        width = 8, height = 9)
