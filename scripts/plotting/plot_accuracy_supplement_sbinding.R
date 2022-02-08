## Plot accuracy comparison for an S-binding specific model
#   - Supplement to "plot_accuracy.R"

## Plot accuracy for infection models using different feature sets (all_data models only)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(cowplot)

  source("scripts/utils/output_utils.R")
  source("scripts/plotting/plotting_constants.R")
})

# ---- Data ----------------------------------------------------------------------------------------
infection_data <- readRDS("data/calculated/cleaned_infection_data.rds")
prop_true <- sum(infection_data$infected == "True") / nrow(infection_data)

preds_all_sites <- load_single_run("all_data", "infection", "aa_distance")
preds_s_binding <- load_single_run("all_data", "infection", "_supplementary_runs/aa_distance_s_binding")

RUN_LABELS <- c("_supplementary_runs/aa_distance_s_binding" = "S-binding sites only",
                "aa_distance" = "All ACE2 sites")

test_preds <- bind_rows(preds_all_sites, preds_s_binding) %>% 
  mutate(run_label = factor(.data$run_id, labels = RUN_LABELS, levels = names(RUN_LABELS)))
  

# ---- Overall accuracy ----------------------------------------------------------------------------
# Showing exact (binomial) confidence intervals
overall_accuracies <- test_preds %>% 
  group_by(.data$run_label) %>% 
  summarise(n_accurate = sum(.data$label == .data$prediction),
            n_total = n(),
            .groups = "keep") %>% 
  mutate(accuracy = .data$n_accurate/.data$n_total,
         lower = binom.test(.data$n_accurate, .data$n_total)$conf.int[1],
         upper = binom.test(.data$n_accurate, .data$n_total)$conf.int[2]) %>% 
  ungroup()

# Accuracy of a null model which assigns labels in proportion to observed proportions
# - Probability of sampling a given label and then independently assigning the correct label =
#   (p_true * p_labeled_true) + (p_false * p_labeled_false)
# - Calculation confirmed by simulation
null_accuracy <- prop_true^2 + (1-prop_true)^2 

# Plot
p_overall <- ggplot(overall_accuracies, aes(x = run_label, y = accuracy)) +
  geom_col(colour = "grey20", fill = "grey50") +
  geom_errorbar(aes(ymin = lower, ymax = upper), colour = "grey20", width = 0.4) + 
  geom_hline(yintercept = null_accuracy, linetype = 2, colour = "grey10") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(0.02)) +
  labs(y = "Accuracy") +
  coord_flip() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank())


# ---- Sensitivity/Specificity ---------------------------------------------------------------------
class_accuracies <- test_preds %>% 
  group_by(.data$run_label, .data$label) %>% 
  summarise(n_accurate = sum(.data$label == .data$prediction),
            n_total = n(),
            .groups = "keep") %>% 
  mutate(accuracy = .data$n_accurate/.data$n_total,
         lower = binom.test(.data$n_accurate, .data$n_total)$conf.int[1],
         upper = binom.test(.data$n_accurate, .data$n_total)$conf.int[2]) %>% 
  ungroup()

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
  coord_flip()

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
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank())


# ---- Output --------------------------------------------------------------------------------------
p_combined <- plot_grid(p_sens, p_spec, p_overall,
                        nrow = 1, rel_widths = c(1.6, 1, 1.05),
                        align = "h", axis = "tb")

ggsave2("output/plots/accuracy_supplement_sbinding.pdf", p_combined, 
        width = 7, height = 1)
