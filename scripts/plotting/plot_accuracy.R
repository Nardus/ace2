## Plot accuracy for infection models using different feature sets (all_data models only)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(cowplot)
})

source("scripts/utils/output_utils.R")
source("scripts/plotting/plotting_constants.R")


# ---- Read all available predictions --------------------------------------------------------------
test_preds <- load_all_runs("all_data") %>% 
  filter(.data$response_var == "Infection")

infection_data <- readRDS("data/calculated/cleaned_infection_data.rds")

# ---- Plotting order / model names -----------------------------------------------------------------
RUN_LABELS <- c("aa_categorical" = "AA categorical",
                "aa_properties" = "AA properties", 
                "aa_distance" = "AA consensus distance",
                "distance_to_humans" = "Distance to humans",
                "distance_to_positive" = "Distance to susceptible",
                "binding_affinity" = "Binding affinity",
                "phylogeny" = "Phylogenetic eigenvectors",
                "all_features" = "All ACE2 representations\ncombined",
                "ensemble" = "ACE2 / phylogeny\nensemble")

test_preds <- test_preds %>% 
  mutate(run_label = factor(.data$run_id, labels = RUN_LABELS, levels = names(RUN_LABELS)),
         run_group = case_when(.data$run_id == "phylogeny" ~ "",
                               .data$run_id == "ensemble" ~ " ",
                               TRUE ~ "ACE2 sequence representations"),
         run_group = factor(.data$run_group, levels = c("ACE2 sequence representations", "", " ")))


# ---- Overall accuracy ----------------------------------------------------------------------------
# Showing exact (binomial) confidence intervals
overall_accuracies <- test_preds %>% 
  group_by(.data$run_label, .data$run_group) %>% 
  summarise(n_accurate = sum(.data$label == .data$prediction),
            n_total = n(),
            .groups = "keep") %>% 
  mutate(accuracy = .data$n_accurate/.data$n_total,
         lower = binom.test(.data$n_accurate, .data$n_total)$conf.int[1],
         upper = binom.test(.data$n_accurate, .data$n_total)$conf.int[2]) %>% 
  ungroup() %>% 
  arrange(.data$accuracy) %>% 
  mutate(run_label = factor(.data$run_label, levels = .data$run_label))

# Accuracy of a null model which assigns labels in proportion to observed proportions
# - Probability of sampling a given label and then independently assigning the correct label =
#   (p_true * p_labeled_true) + (p_false * p_labeled_false)
# - Calculation confirmed by simulation
prop_true <- sum(infection_data$infected == "True") / nrow(infection_data)
null_accuracy <- prop_true^2 + (1-prop_true)^2 

# Plot
p_overall <- ggplot(overall_accuracies, aes(x = run_label, y = accuracy)) +
  geom_col(colour = "grey20", fill = "grey50") +
  geom_errorbar(aes(ymin = lower, ymax = upper), colour = "grey20", width = 0.4) + 
  geom_hline(yintercept = null_accuracy, linetype = 2, colour = "grey10") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(0.02)) +
  labs(y = "Accuracy") +
  coord_flip() +
  facet_grid(rows = vars(run_group), scales = "free_y", space = "free_y") +
  theme(strip.background = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())


# ---- Sensitivity/Specificity ---------------------------------------------------------------------
class_accuracies <- test_preds %>% 
  group_by(.data$run_label, .data$run_group, .data$label) %>% 
  summarise(n_accurate = sum(.data$label == .data$prediction),
            n_total = n(),
            .groups = "keep") %>% 
  mutate(accuracy = .data$n_accurate/.data$n_total,
         lower = binom.test(.data$n_accurate, .data$n_total)$conf.int[1],
         upper = binom.test(.data$n_accurate, .data$n_total)$conf.int[2]) %>% 
  ungroup() %>% 
  mutate(run_label = factor(.data$run_label, levels = overall_accuracies$run_label))

# Sensitivity
sens <- class_accuracies %>% 
  filter(.data$label == "True")

p_sens <- ggplot(sens, aes(x = run_label, y = accuracy, fill = label)) +
  geom_col(colour = "grey20") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.4, colour = "grey20") + 
  geom_hline(yintercept = prop_true, linetype = 2, colour = "grey10") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(0.02)) +
  scale_fill_manual(values = INFECTION_STATUS_COLOURS, guide = "none") +
  labs(x = "Predictors used", y = "Sensitivity") +
  coord_flip() +
  facet_grid(rows = vars(run_group), scales = "free_y", space = "free_y") +
  theme(strip.text = element_blank())

# Specificity
spec <- class_accuracies %>% 
  filter(.data$label == "False")

p_spec <- ggplot(spec, aes(x = run_label, y = accuracy, fill = label)) +
  geom_col(colour = "grey20") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.4, colour = "grey20") + 
  geom_hline(yintercept = 1-prop_true, linetype = 2, colour = "grey10") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(0.02)) +
  scale_fill_manual(values = INFECTION_STATUS_COLOURS, guide = "none") +
  labs(y = "Specificity") +
  coord_flip() +
  facet_grid(rows = vars(run_group), scales = "free_y", space = "free_y") +
  theme(strip.text = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.y = element_blank())


# ---- Output --------------------------------------------------------------------------------------
p_combined <- plot_grid(p_sens, p_spec, p_overall,
                        nrow = 1, rel_widths = c(1.6, 1, 1.05),
                        align = "h", axis = "tb")

ggsave2("output/plots/accuracy.pdf", p_combined, 
        width = 7, height = 2.2)


# ---- Values quoted in text -----------------------------------------------------------------------
focal_models <- RUN_LABELS[c("all_features", "phylogeny", "ensemble")]

cat("Overall accuracy:\n")
overall_accuracies %>% 
  filter(.data$run_label %in% focal_models) %>% 
  print()

cat("\nClass accuracy:\n")
class_accuracies %>% 
  filter(.data$run_label %in% focal_models) %>% 
  print()
