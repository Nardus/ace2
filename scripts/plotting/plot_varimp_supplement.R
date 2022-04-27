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
feature_importance <- readRDS("output/all_data/infection/all_features/feature_importance.rds")

feature_clusters <- readRDS("output/plots/intermediates/feature_clusters.rds")


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
    mutate(feature_label = factor(.data$feature_label, levels = .data$feature_label),
           feature_type = factor(.data$feature_type)) %>%
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

p_site_importance <- ggplot(site_importance, aes(x = site_label_continuous, y = importance,
                                                 fill = feature_type)) +
    geom_col(colour = "grey20", size = 0.2, position = "stack") +
    coord_flip() +
    scale_x_continuous(breaks = site_importance$site_label_continuous,
                       labels = site_importance$site_label) +
    scale_fill_brewer(palette = "Set2", na.value = "grey60", drop = FALSE) +
    labs(y = "Total effect magnitude", x = "Sequence position (human ACE2)", fill = "Measure") +
    theme(legend.position = c(0.66, 0.16))


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

# Expected ratio
s_clusters <- feature_clusters %>%
    left_join(cluster_binding, by = "cluster") %>%
    filter(!.data$feature_position_corrected %in% ALL_S_BINDING_INDS) %>%
    mutate(selected_site = .data$feature_position %in% top_importance$feature_position)

stopifnot(n_distinct(s_clusters$feature_position) == nrow(s_clusters))

exp_ratio <- sum(s_clusters$includes_binding_site) / nrow(s_clusters)
sprintf(
    "%3.1f%% of available positions are correlated with an S-interacting site\n",
    exp_ratio * 100
) %>%
    cat()

# Observed ratio
observed <- sum(s_clusters$selected_site & s_clusters$includes_binding_site)
available <- sum(s_clusters$selected_site)
test_result <- binom.test(
    x = observed,
    n = available,
    p = exp_ratio
)

sprintf(
    "%i of %i selected sites clustered with S-binding sites (%3.1f%%, p-value compared to expected ratio = %3.3f)\n\n",
    observed,
    available,
    test_result$estimate * 100,
    test_result$p.value
) %>%
    cat()
