## Plot an overview of quantitative predictions from the phylogeny-only model

suppressPackageStartupMessages({
    library(dplyr)
    library(stringr)
    library(ape)
    library(ggplot2)
    library(cowplot)
    
    source("scripts/utils/timetree_constants.R")
    source("scripts/plotting/plotting_constants.R")
})

# ---- Data ----------------------------------------------------------------------------------------
all_predictions <- readRDS("output/all_data/infection/phylogeny/holdout_predictions.rds")

mammal_tree <- read.tree("data/internal/timetree_mammalia.nwk")


# ---- Clean up species names ----------------------------------------------------------------------
mammal_tree$tip.label <- str_replace(mammal_tree$tip.label, "_", " ")

mammal_tree$tip.label <- with(mammal_tree,
                              if_else(tip.label %in% names(TIMETREE_TAXONOMY_CORRECTIONS),
                                      as.character(TIMETREE_TAXONOMY_CORRECTIONS[tip.label]),
                                      tip.label))


# ---- Plot mammalian scores -----------------------------------------------------------------------
mammal_predictions <- all_predictions %>% 
    filter(.data$species %in% mammal_tree$tip.label) %>% 
    arrange(.data$probability) %>%
    mutate(species = factor(.data$species, levels = .data$species))

# Plot separately by label so we can control drawing order
predictions_unknown <- mammal_predictions %>% 
    filter(is.na(.data$label))
    
predictions_true <- mammal_predictions %>% 
    filter(.data$label == "True")
    
predictions_false <- mammal_predictions %>%
    filter(.data$label == "False")


# Add a rectangle around the top X% of predictions
get_top_species <- function(p, preds = mammal_predictions, all_preds = all_predictions) {
    cutoff <- quantile(preds$probability, p = p)
    top_spp_names <- preds$species[preds$probability >= cutoff]
    
    # Proportion of susceptible predictions included
    prop_susceptible <- preds %>% 
        ungroup() %>% 
        mutate(n_predicted = sum(.data$predicted_label == "True")) %>% 
        filter(.data$species %in% top_spp_names) %>% 
        summarise(prop = sum(.data$predicted_label == "True") / unique(.data$n_predicted)) %>% 
        pull(.data$prop)
    
    # Proportion of known susceptibles included
    captured <- all_preds %>% 
        filter(.data$label == "True") %>% 
        summarise(n_captured = sum(.data$species %in% top_spp_names),
                  n_total = n(),
                  prop = .data$n_captured / .data$n_total) %>% 
        pull(.data$prop)
    
    data.frame(lowest_sp = top_spp_names[1],  # Species already in descending order
               highest_sp = top_spp_names[length(top_spp_names)],
               cutoff = cutoff,
               max_score = max(preds$probability),
               proportion_pred_susceptible = prop_susceptible,
               proportion_known_susceptible = captured,
               n_species = length(top_spp_names)) %>% 
        mutate(cutoff_label = sprintf("Top %0.0f%%", (1 - p) * 100),
               size_label = sprintf("N=%0i (%0.0f%% of susceptibles)", .data$n_species, captured * 100))
}

top_species <- lapply(c(0.5, 0.6, 0.7, 0.8), get_top_species) %>% 
    bind_rows()

# Plot
pred_cutoff <- unique(all_predictions$cutoff)

p <- ggplot(mammal_predictions, aes(x = species, y = probability, colour = label)) +
    geom_blank() +
    geom_hline(yintercept = pred_cutoff, colour = "grey20", linetype = 2) +
    
    geom_rect(aes(xmin = lowest_sp, xmax = highest_sp, 
                  ymin = cutoff, ymax = max_score,
                  x = lowest_sp, y = cutoff), 
              data = top_species,
              colour = "firebrick", fill = NA, size = 0.5) +
    geom_text(aes(x = lowest_sp, y = max_score, label = cutoff_label), 
              data = top_species,
              size = 2, colour = "black", 
              angle = 90, hjust = 1.1, vjust = -0.7) +
    geom_text(aes(x = highest_sp, y = cutoff, label = size_label), 
              data = top_species,
              size = 2, colour = "black", 
              hjust = 1.02, vjust = 1.7) +
    
    geom_point(data = predictions_unknown) +
    geom_point(data = predictions_false) +
    geom_point(data = predictions_true) +
    
    scale_y_continuous(expand = c(0, 0)) +
    scale_colour_manual(values = INFECTION_STATUS_COLOURS, na.value = MISSING_DATA_COLOUR) +
    
    labs(x = "Species", y = "Predicted score", colour = "Infected") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())

ggsave2("output/plots/phylogeny_predictions_supplement.pdf", p, width = 7, height = 5)


# ---- Values in text ------------------------------------------------------------------------------
cat("\n\nProportion of known susceptibles captured by each cutoff:\n")

top_species %>% 
    select(.data$cutoff_label, 
           .data$proportion_pred_susceptible, 
           .data$proportion_known_susceptible) %>% 
    print()
