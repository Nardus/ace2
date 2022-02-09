## Plot a comparison of predictions from previous studies alongside current data

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(ape)
  library(ggplot2)
  library(ggtree)
  library(ggtext)
  library(cowplot)
  
  source("scripts/utils/timetree_constants.R")
  source("scripts/plotting/plotting_constants.R")
})


infection_data <- readRDS("data/calculated/cleaned_infection_data.rds") %>% 
  mutate(infected = factor(.data$infected, levels = c("True", "False")))

virus_data_infection <- readRDS("output/plots/intermediates/virus_data_infection.rds")

timetree <- read.tree("data/internal/timetree_amniota.nwk")

prediction_dendrogram <- readRDS("output/plots/intermediates/prediction_dendrogram.rds")
all_predictions <- readRDS("output/plots/intermediates/unified_predictions.rds")


# ---- Filter to species with known infection data -------------------------------------------------
all_predictions <- all_predictions %>%
  filter(.data$species %in% infection_data$species)

# Remove studies with <10 species remaining (to reduce clutter)
all_predictions <- all_predictions %>% 
  group_by(.data$citation_key) %>% 
  filter(n() >= 10) %>% 
  ungroup()

prediction_dendrogram <- keep.tip(prediction_dendrogram, unique(all_predictions$citation_key))


# ---- Prepare phylogeny ---------------------------------------------------------------------------
# Correct names
timetree$tip.label <- str_replace(timetree$tip.label, "_", " ")

stopifnot(all(names(TIMETREE_TAXONOMY_CORRECTIONS) %in% timetree$tip.label))

timetree$tip.label <- with(timetree,
                           if_else(tip.label %in% names(TIMETREE_TAXONOMY_CORRECTIONS),
                                   as.character(TIMETREE_TAXONOMY_CORRECTIONS[tip.label]),
                                   tip.label))

stopifnot(all(all_predictions$species %in% timetree$tip.label))

timetree <- keep.tip(timetree, unique(all_predictions$species))


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


# Study correlation dendrogram
dendro_panel <- ggtree(prediction_dendrogram, ladderize = FALSE) +
  scale_x_reverse(expand = expansion(add = c(0.002, 0.008))) +
  scale_y_discrete(expand = expansion(add = 0.2)) +
  coord_flip() +
  theme(line = PLOT_THEME$line,
        plot.margin = margin(t = 5.5, r = 92, b = 0, l = 0))


# Data panel
observed_status <- infection_data %>% 
  select(.data$species, .data$infected, .data$evidence_level)

all_predictions <- all_predictions %>% 
  left_join(observed_status, by = "species") %>% 
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
         study_label_long = paste0(.data$study_label, "<br/>(", .data$predictor, ")"))

study_labels <- study_label_df$study_label_long
names(study_labels) <- as.character(study_label_df$citation_key)


data_panel <- ggplot(all_predictions) +
  geom_hline(aes(yintercept = species, linetype = infected), colour = "grey92", size = 0.4) +
  geom_tile(aes(x = citation_key, y = species, fill = scaled_score)) +
  
  scale_linetype_manual(values = c("True" = 1, "False" = 2)) +
  scale_fill_viridis_c(direction = -1) +
  scale_x_discrete(labels = study_labels, expand = expansion(add = 0)) +
  scale_y_discrete(expand = expansion(add = 0.5), position="right") +
  
  labs(x = "Prediction source", y = NULL) +
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(face = "italic"),
        plot.background = element_blank(),
        legend.position="none",
        plot.margin = margin(t = 1, r = 5.5, b = 5.5, l = 1))



# Combine
shared_legend <- get_legend(
  ggplot(all_predictions, aes(x = citation_key, y = scaled_score, fill = scaled_score,
                                   colour = infected, linetype = infected)) +
    geom_line() +
    geom_point(shape = 21, colour = "grey20") +
    scale_fill_viridis_c(direction = -1, guide = guide_colorbar(order = 1)) +
    scale_linetype_manual(values = c("True" = 1, "False" = 2), 
                          guide = guide_legend(order = 2, keywidth = unit(1.25, "lines"))) +
    scale_colour_manual(values = c(True = "#EE7733", False = "#009988"), 
                        na.value = "grey30", guide = guide_legend(order = 2)) +
    labs(fill = "Scaled score", colour = "Infected", linetype = "Infected")
)

phylo_panel2 <- phylo_panel +
  annotation_custom(shared_legend, xmin = -420, ymin = 20)

overview_plot <- plot_grid(NULL, dendro_panel,
                           phylo_panel2, data_panel,
                           nrow = 2, rel_heights = c(1, 20), rel_widths = c(1.2, 5))


# ---- Binary predictions --------------------------------------------------------------------------
# Not all studies give binary predictions - we can't assign a cutoff after the fact
# - For Damas et al., consider "Very high" to "medium" as a positive prediction
binary_preds <- all_predictions %>% 
  filter(!is.na(.data$prediction)) %>% 
  filter(.data$citation_key != "kumar2021") %>% 
  mutate(observed = .data$infected == "True",
         binary_pred = case_when(startsWith(as.character(.data$citation_key), "This study") ~ .data$prediction == "True",
                                 .data$citation_key == "damas2020" ~ .data$prediction %in% c("Very High", "High", "Medium"),
                                 .data$citation_key == "melin2020" ~ .data$prediction == "High",
                                 TRUE ~ .data$prediction == "TRUE"))

# Check that all cases were caught:
needs_correction <- binary_preds %>% 
  mutate(citation_key = as.character(.data$citation_key),
         nn = startsWith(.data$citation_key, "This study") | .data$citation_key %in% c("damas2020", "melin2020")) %>% 
  pull(.data$nn)

stopifnot(all(binary_preds$prediction[!needs_correction] %in% c("TRUE", "FALSE")))

# Sample size by study
sample_size <- binary_preds %>% 
  group_by(.data$citation_key) %>% 
  summarise(n = n(),
            .groups = "drop")

# Simpler study labels
simple_label_df <- study_label_df %>% 
  mutate(citation_key = as.character(.data$citation_key),
         study_label = if_else(startsWith(.data$citation_key, "This study"), 
                               sprintf("**%s**", str_to_sentence(.data$predictor)),
                               .data$study_label))

simple_labels <- simple_label_df$study_label
names(simple_labels) <- simple_label_df$citation_key


# ---- Overall accuracy ----------------------------------------------------------------------------
total_accuracies <- binary_preds %>% 
  group_by(.data$citation_key) %>% 
  summarise(n_accurate = sum(.data$observed == .data$binary_pred),
            n_total = n(),
            .groups = "keep") %>% 
  mutate(accuracy = .data$n_accurate/.data$n_total,
         lower = binom.test(.data$n_accurate, .data$n_total)$conf.int[1],
         upper = binom.test(.data$n_accurate, .data$n_total)$conf.int[2]) %>% 
  ungroup()

# Accuracy of a null model which assigns labels in proportion to observed proportions  (based on full dataset)
prop_true <- sum(infection_data$infected == "True") / nrow(infection_data)
null_accuracy <- prop_true^2 + (1-prop_true)^2 


p_acc <- ggplot(total_accuracies, aes(x = citation_key, y = accuracy, fill = n_total)) +
  geom_col(colour = "grey50") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.4, colour = "grey40") + 
  geom_hline(yintercept = null_accuracy, linetype = 2, colour = "grey20") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(0.02)) +
  scale_fill_distiller(palette = "PuBu", direction = 1, guide = "none") +
  labs(y = "Accuracy") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())


# ---- Sensitivity/specificity ---------------------------------------------------------------------
class_accuracies <- binary_preds %>% 
  group_by(.data$citation_key, .data$infected) %>% 
  summarise(n_accurate = sum(.data$observed == .data$binary_pred),
            n_total = n(),
            .groups = "keep") %>% 
  mutate(accuracy = .data$n_accurate/.data$n_total,
         lower = binom.test(.data$n_accurate, .data$n_total)$conf.int[1],
         upper = binom.test(.data$n_accurate, .data$n_total)$conf.int[2]) %>% 
  ungroup() %>% 
  left_join(sample_size, by = "citation_key")


# Sensitivity
sens <- class_accuracies %>% 
  filter(.data$infected == "True")

p_sens <- ggplot(sens, aes(x = citation_key, y = accuracy, fill = n)) +
  geom_col(colour = "grey50") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.4, colour = "grey40") + 
  geom_hline(yintercept = prop_true, linetype = 2, colour = "grey20") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(0.02)) +
  scale_fill_distiller(palette = "PuBu", direction = 1, guide = "none") +
  labs(y = "Sensitivity") +
  theme(strip.text = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

# Specificity
spec <- class_accuracies %>% 
  filter(.data$infected == "False")

p_spec <- ggplot(spec, aes(x = citation_key, y = accuracy, fill = n)) +
  geom_col(colour = "grey50") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.4, colour = "grey40") + 
  geom_hline(yintercept = 1 -prop_true, linetype = 2, colour = "grey20") +
  scale_x_discrete(labels = simple_labels) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(0.02)) +
  scale_fill_distiller(palette = "PuBu", direction = 1, guide = "none") +
  labs(x = "Prediction source", y = "Specificity") +
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5))


# Combined
cl_spread <- sens %>% 
  rename(sens = .data$accuracy,
         sens_lower = .data$lower,
         sens_upper = .data$upper,
         n_pos = .data$n_total) %>% 
  left_join(spec, by = "citation_key") %>% 
  rename(spec = .data$accuracy,
         spec_lower = .data$lower,
         spec_upper = .data$upper,
         n_neg = .data$n_total) %>% 
  mutate(sample_size = .data$n_pos + .data$n_neg)


cl_labels <- cl_spread %>% 
  filter(.data$sample_size >= 60) %>% 
  left_join(simple_label_df, by = "citation_key")


cl_plot <- ggplot(cl_spread, aes(x = sens, y = spec, colour = sample_size)) +
  
  geom_abline(linetype = 2, colour = "grey20") +
  
  geom_point() +
  geom_errorbar(aes(xmin = sens_lower, xmax = sens_upper), width = 0.02) +
  geom_errorbar(aes(ymin = spec_lower, ymax = spec_upper), width = 0.02) +
  
  geom_richtext(aes(label = study_label), data = cl_labels,
                size = 1.5, label.color = NA, 
                nudge_x = -0.02, nudge_y = -0.02, angle = 90, hjust = 1,
                label.padding = unit(0, "lines")) +
  
  labs(x = "Sensitivity", y = "Specificity", colour = "Number of\nspecies") +
  lims(x = c(0, 1), y = c(0, 1)) +
  coord_equal(expand = FALSE) +
  scale_colour_distiller(palette = "PuBu", direction = 1) +
  theme(legend.position = "top",
        plot.margin = margin(t = 5.5, r = 8, b = 5.5, l = 5.5))


# ---- Output --------------------------------------------------------------------------------------
p_top <- plot_grid(cl_plot, p_acc, p_sens, p_spec,
                   ncol = 1, rel_heights = c(3, 1, 1.2, 2),
                   align = "v", axis = "lr",
                   labels = c("B", "C", "", ""),
                   hjust = 0.5, vjust = c(1.5, 0.5))

final_figure <- plot_grid(overview_plot, p_top, 
                          nrow = 1, rel_widths = c(1.9, 1),
                          labels = c("A", ""))

ggsave2("output/plots/existing_predictions.pdf", 
        final_figure, 
        width = 8, height = 8)


# ---- Stats/citations for text --------------------------------------------------------------------
all_predictions <- readRDS("output/plots/intermediates/unified_predictions.rds") # Unfiltered version

# Correlation
score_mat <- all_predictions %>% 
  select(.data$species, .data$citation_key, .data$scaled_score) %>% 
  pivot_wider(id_cols = "species", names_from = "citation_key", values_from = "scaled_score") %>% 
  column_to_rownames("species")

cor_mat <- cor(score_mat, use = "pairwise.complete", method = "spearman")


# - Summarize
dfcor <- cor_mat %>% 
  as_tibble(rownames = "from") %>% 
  pivot_longer(-.data$from, names_to = "to", values_to = "cor") %>% 
  filter(.data$from != .data$to) %>% 
  group_by(.data$from) %>% 
  mutate(min_cor = min(.data$cor, na.rm = TRUE),
         max_cor = max(.data$cor, na.rm = TRUE)) %>% 
  ungroup()

# - Report
cat("\n\nMinimum correlations:\n")
dfcor %>% 
  filter(.data$cor == .data$min_cor) %>% 
  select(-.data$min_cor, -.data$max_cor) %>% 
  arrange(.data$cor) %>% 
  print()

cat("\n\nMaximum correlations:\n")
dfcor %>% 
  filter(.data$cor == .data$max_cor) %>% 
  select(-.data$min_cor, -.data$max_cor) %>% 
  arrange(.data$cor) %>% 
  print()


# - Overlap with our models:
ace2_spp <- all_predictions$species[all_predictions$citation_key == "This study (ACE2-based)"]
phylo_spp <- all_predictions$species[all_predictions$citation_key == "This study (host phylogeny)"]
huang_spp <- all_predictions$species[all_predictions$citation_key == "huang2020"]
alexander_spp <- all_predictions$species[all_predictions$citation_key == "alexander2020"]

cat("\nACE2-based model shares", sum(ace2_spp %in% huang_spp), "species with Huang et al. 2020\n")
cat("Ensemble model shares", sum(ace2_spp %in% alexander_spp), "species with Alexander et al. 2020\n")
cat("Phylogeny-based model shares", sum(phylo_spp %in% alexander_spp), "species with Alexander et al. 2020\n")


# - Accuracy of directly compared to Huang et al.
huang_accuracies <- binary_preds %>% 
  filter(.data$species %in% huang_spp) %>% 
  filter(.data$citation_key == "huang2020" | 
           startsWith(as.character(.data$citation_key), "This study")) %>% 
  group_by(.data$citation_key, .data$infected) %>% 
  summarise(n_accurate = sum(.data$observed == .data$binary_pred),
            n_total = n(),
            .groups = "keep") %>% 
  mutate(accuracy = .data$n_accurate/.data$n_total,
         lower = binom.test(.data$n_accurate, .data$n_total)$conf.int[1],
         upper = binom.test(.data$n_accurate, .data$n_total)$conf.int[2]) %>% 
  ungroup()

cat("\n\nAccuracy when predicting the exact species used to determine performance of", 
    "the Huang et al. model:\n")
print(huang_accuracies)
