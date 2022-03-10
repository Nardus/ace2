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

all_predictions <- readRDS("output/plots/intermediates/unified_predictions.rds")

ace2_cv_preds <- readRDS("output/all_data/infection/all_features/predictions.rds")
phylo_cv_preds <- readRDS("output/all_data/infection/phylogeny/predictions.rds")
ensemble_cv_preds <- readRDS("output/all_data/infection/ensemble_all_features_phylogeny/predictions.rds")


# ---- Filter to species with known infection data -------------------------------------------------
all_predictions <- all_predictions %>%
  filter(.data$species %in% infection_data$species)

# Remove studies with <10 species remaining (to reduce clutter)
all_predictions <- all_predictions %>% 
  group_by(.data$citation_key) %>% 
  filter(n() >= 10) %>% 
  ungroup()


# ---- Add our holdout predictions -----------------------------------------------------------------
all_predictions <- all_predictions %>% 
  filter(!startsWith(.data$citation_key, "This study")) # These are fitted values

# Unify column names
ace2_cv_preds <- ace2_cv_preds %>% 
  select(.data$species, .data$prediction,
         raw_score = .data$p_true) %>% 
  mutate(citation_key = "This study (ACE2-based)",
         predictor = "ACE2-based",
         prediction_type = "sequence-based")


phylo_cv_preds <- phylo_cv_preds %>% 
  select(.data$species, .data$prediction,
         raw_score = .data$p_true) %>% 
  mutate(citation_key = "This study (host phylogeny)",
         predictor = "host phylogeny",
         prediction_type = "phylogeny")


ensemble_cv_preds <- ensemble_cv_preds %>% 
  select(.data$species, .data$prediction,
         raw_score = .data$p_true) %>% 
  mutate(citation_key = "This study (ensemble)",
         predictor = "ensemble",
         prediction_type = "ensemble")

# Merge
internal_predictions <- ace2_cv_preds %>% 
  bind_rows(phylo_cv_preds) %>% 
  bind_rows(ensemble_cv_preds) %>% 
  mutate(reverse_score = FALSE)

all_predictions <- all_predictions %>% 
  bind_rows(internal_predictions)


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


# ---- Produce binary predictions ------------------------------------------------------------------
# Not all studies give binary predictions - we can't assign a cutoff after the fact
# - For Damas et al., consider "Very high" to "medium" as a positive prediction
observed_status <- infection_data %>% 
  select(.data$species, .data$infected, .data$evidence_level)

binary_preds <- all_predictions %>% 
  left_join(observed_status, by = "species") %>% 
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


# Study labels
study_label_df <- all_predictions %>% 
  distinct(.data$citation_key, .data$prediction_type, .data$predictor) %>% 
  mutate(predictor = if_else(.data$predictor == "binding affinity change relative to humans",
                             "binding affinity", 
                             .data$predictor),
         predictor_label = case_when(startsWith(.data$citation_key, "This study") ~ "**This study**<br>(holdout predictions)",
                                    .data$predictor == "host features" ~ "Host<br>features",
                                    TRUE ~ str_to_sentence(.data$predictor)),
         
         study_label = if_else(startsWith(.data$citation_key, "This study"), 
                               as.character(.data$predictor), 
                               as.character(.data$citation_key)),
         study_label = if_else(.data$study_label == "ACE2-based", 
                               "ACE2-based",
                               str_to_sentence(.data$study_label)),
         study_label = str_replace(.data$study_label, "([[:digit:]]{4})[[:alpha:]]*$", " *et al.*, \\1")) %>% 
  mutate(predictor_label = factor(.data$predictor_label, 
                                  levels = c("Distance to humans", "Binding affinity", 
                                             "Host<br>features", 
                                             "**This study**<br>(holdout predictions)"))) %>% 
  select(.data$citation_key, .data$study_label, .data$predictor_label)


# ---- Overall accuracy ----------------------------------------------------------------------------
total_accuracies <- binary_preds %>% 
  group_by(.data$citation_key) %>% 
  summarise(n_accurate = sum(.data$observed == .data$binary_pred),
            n_total = n(),
            .groups = "keep") %>% 
  mutate(accuracy = .data$n_accurate/.data$n_total,
         lower = binom.test(.data$n_accurate, .data$n_total)$conf.int[1],
         upper = binom.test(.data$n_accurate, .data$n_total)$conf.int[2]) %>% 
  ungroup() %>% 
  left_join(study_label_df, by = "citation_key") %>% 
  arrange(.data$n_total, .data$accuracy, decreasing = TRUE) %>%
  mutate(study_label = factor(.data$study_label, levels = .data$study_label))

# Accuracy of a null model which assigns labels in proportion to observed proportions  (based on full dataset)
prop_true <- sum(infection_data$infected == "True") / nrow(infection_data)
null_accuracy <- prop_true^2 + (1-prop_true)^2 


p_acc <- ggplot(total_accuracies, aes(x = study_label, y = accuracy, fill = n_total)) +
  geom_col(colour = "grey50") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.4, colour = "grey40") + 
  geom_hline(yintercept = null_accuracy, linetype = 2, colour = "grey20") +
  
  scale_y_continuous(limits = c(0, 1), expand = expansion(0.02)) +
  scale_fill_distiller(palette = "PuBu", direction = 1) +
  
  facet_grid(cols = vars(predictor_label), scale = "free", space = "free") +
  
  labs(x = "Prediction source", y = "Accuracy", fill = "Number of\nspecies") +
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5),
        strip.text.x = element_markdown(),
        legend.position = "top")


# ---- Sensitivity/specificity ---------------------------------------------------------------------
class_accuracies <- binary_preds %>% 
  group_by(.data$citation_key, .data$infected) %>% 
  summarise(n_accurate = sum(.data$observed == .data$binary_pred),
            n_total = n(),
            .groups = "keep") %>% 
  mutate(accuracy = .data$n_accurate/.data$n_total,
         lower = binom.test(.data$n_accurate, .data$n_total)$conf.int[1],
         upper = binom.test(.data$n_accurate, .data$n_total)$conf.int[2]) %>% 
  ungroup()

sens <- class_accuracies %>% 
  filter(.data$infected == "True")

spec <- class_accuracies %>% 
  filter(.data$infected == "False")

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
  left_join(study_label_df, by = "citation_key") %>% 
  mutate(study_label = if_else(startsWith(.data$citation_key, "This study"),
                               paste0("**", .data$study_label, "**"),
                               .data$study_label))


cl_plot <- ggplot(cl_spread, aes(x = sens, y = spec, colour = sample_size)) +
  
  geom_abline(linetype = 2, colour = "grey20") +
  
  geom_point() +
  geom_errorbar(aes(xmin = sens_lower, xmax = sens_upper), width = 0.02) +
  geom_errorbar(aes(ymin = spec_lower, ymax = spec_upper), width = 0.02) +
  
  geom_richtext(aes(label = study_label), data = cl_labels,
                size = 1.9, label.color = NA, 
                nudge_x = -0.02, nudge_y = -0.02, angle = 90, hjust = 1,
                label.padding = unit(0, "lines")) +
  
  scale_x_continuous(expand = expansion(add = c(0, 0.03))) +
  scale_y_continuous(expand = expansion(add = c(0, 0.03))) +
  
  labs(x = "Sensitivity", y = "Specificity") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_colour_distiller(palette = "PuBu", direction = 1, guide = "none") +
  theme(plot.margin = margin(t = 5.5, r = 8, b = 5.5, l = 5.5))


# ---- Output --------------------------------------------------------------------------------------
final_figure <- plot_grid(p_acc, cl_plot, 
                          ncol = 2, rel_widths = c(1.3, 1),
                          align = "v", axis = "tb",
                          labels = c("A", "B"),
                          hjust = -0.2)

ggsave2("output/plots/existing_predictions.pdf", 
        final_figure, 
        width = 7.5, height = 3.5)


# ---- Stats/citations for text --------------------------------------------------------------------
# Correlation (based on fitted values, see get_prediction_dendrogram.R)
cor_mat <- readRDS("output/plots/intermediates/study_correlation_matrix.rds")

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
cat("\n\nMax correlation between our models and previous studies:\n")
dfcor %>% 
  filter(startsWith(.data$from, "This study")) %>%
  filter(!startsWith(.data$to, "This study")) %>%
  group_by(.data$from) %>% 
  mutate(min_cor = min(.data$cor, na.rm = TRUE),
         max_cor = max(.data$cor, na.rm = TRUE)) %>% 
  ungroup() %>% 
  filter(.data$cor == .data$max_cor) %>% 
  select(-.data$min_cor, -.data$max_cor) %>% 
  arrange(.data$cor) %>% 
  print()

ace2_spp <- all_predictions$species[all_predictions$citation_key == "This study (ACE2-based)"]
phylo_spp <- all_predictions$species[all_predictions$citation_key == "This study (host phylogeny)"]
huang_spp <- all_predictions$species[all_predictions$citation_key == "huang2020"]
melin_spp <- all_predictions$species[all_predictions$citation_key == "melin2020"]

cat("\nACE2-based model shares", sum(ace2_spp %in% melin_spp), "species with Melin et al. 2020\n")
cat("Ensemble model shares", sum(ace2_spp %in% melin_spp), "species with Melin et al. 2020\n")
cat("Phylogeny-based model shares", sum(phylo_spp %in% huang_spp), "species with Huang et al. 2020\n")


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
