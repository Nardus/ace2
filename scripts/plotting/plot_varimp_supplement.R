# Plot detailed feature usage for the combined ACE2 model

suppressPackageStartupMessages({
    library(dplyr)
    library(tibble)
    library(tidyr)
    library(stringr)
    library(ggplot2)
    library(cowplot)

    source("scripts/utils/plot_utils.R")
    source("scripts/plotting/plotting_constants.R")
})


# ---- Data ----------------------------------------------------------------------------------------
metadata <- readRDS("data/calculated/cleaned_infection_data.rds")

feature_importance <- readRDS("output/all_data/infection/all_features/feature_importance.rds")
shap_values <- readRDS("output/all_data/infection/all_features/shap_values.rds")


# ---- Feature importance --------------------------------------------------------------------------
top_importance <- feature_importance %>%
    filter(.data$importance > 0) %>%
    add_readable_feature_names() %>%
    arrange(.data$importance) %>%
    mutate(feature_label = factor(.data$feature_label, levels = .data$feature_label),
           feature_type = factor(.data$feature_type))

p_overall_importance <- ggplot(top_importance, aes(x = feature_label, y = importance,
                               fill = feature_type)) +
    geom_col(colour = "grey20", size = 0.2) +
    coord_flip() +
    scale_fill_brewer(palette = "Set2", na.value = "grey60", guide = "none") +
    labs(y = "Effect magnitude", x = "Selected feature", fill = "Measure")


# ---- Inset - compare reservoir and other species -------------------------------------------------
spp_groups <- metadata %>% 
  mutate(species_group = if_else(str_starts(.data$species, "Rhinolophus"),
                                 "Rhinolophid bats", 
                                 "Other species")) %>% 
  select(.data$species_group)

group_importance <- shap_values %>% 
  data.frame() %>% 
  bind_cols(spp_groups) %>% 
  pivot_longer(!species_group, names_to = "feature", values_to = "shap") %>% 
  group_by(.data$species_group, .data$feature) %>% 
  summarise(importance = mean(abs(.data$shap)), .groups = "drop") %>% 
  filter(.data$feature != "BIAS") %>% 
  filter(.data$importance > 0) %>%
  add_readable_feature_names() %>% 
  pivot_wider(id_cols=c("feature", "feature_type"),
              names_from = "species_group", values_from = "importance")

group_max <- max(group_importance$`Rhinolophid bats`, group_importance$`Other species`)

inset <- ggplot(group_importance, aes(x = `Rhinolophid bats`, y = `Other species`, colour = feature_type)) +
  geom_abline(linetype = 2, colour = "grey80") +
  geom_point() +
  scale_colour_brewer(palette = "Set2", na.value = "grey60", guide = "none") +
  scale_x_log10(limits = c(0.001, group_max)) +
  scale_y_log10(limits = c(0.001, group_max)) +
  coord_equal()

p_overall_importance <- ggdraw(p_overall_importance) +
  draw_plot(inset, 
            x = 0.45, y = 0.05, 
            width = 0.5, height = 0.5)


# ---- Importance by site --------------------------------------------------------------------------
site_labels <- top_importance %>%
    filter(!is.na(.data$feature_position)) %>%
    group_by(.data$feature_position, .data$feature_position_corrected) %>%
    summarise(total_importance = sum(.data$importance),
              .groups = "drop") %>%
    arrange(.data$total_importance) %>%
    mutate(site_label = factor(.data$feature_position_corrected,
                               levels = .data$feature_position_corrected),
           site_label_continuous = as.numeric(.data$site_label))

site_importance <- top_importance %>%
    left_join(site_labels, by = c("feature_position", "feature_position_corrected"))

site_labels <- site_importance %>% 
  distinct(.data$site_label_continuous, .data$site_label)

p_site_importance <- ggplot(site_importance, aes(x = site_label_continuous, y = importance,
                                                 fill = feature_type)) +
    geom_col(colour = "grey20", size = 0.2, position = "stack") +
    coord_flip() +
    scale_x_continuous(breaks = site_labels$site_label_continuous,
                       labels = site_labels$site_label) +
    scale_fill_brewer(palette = "Set2", na.value = "grey60", drop = FALSE) +
    labs(y = "Total effect magnitude", x = "Sequence position (human ACE2)", fill = "Measure") +
    theme(legend.position = c(0.74, 0.16))


# ---- Combine--------------------------------------------------------------------------------------
p_combined <- plot_grid(p_overall_importance, p_site_importance,
                        ncol = 2, rel_widths = c(1.4, 1),
                        labels = c("A", "B"))

ggsave2("output/plots/varimp_supplement.pdf", p_combined, width = 7, height = 4.7)


# ---- Values mentioned in text --------------------------------------------------------------------
cat("\n\nSequence positions given with reference to accesion:", human_acc, "\n\n")

cat("Model uses", nrow(top_importance), "features, representing", nrow(site_labels), "sites:\n")
print(table(top_importance$feature_type))

cat("Site-specific features represent", n_distinct(site_labels$site_label), "amino acid positions\n")

# S-interacting sites
cat(
    sum(top_importance$feature_position_corrected %in% ALL_S_BINDING_INDS),
    "included sites are known to interact with S\n"
)
