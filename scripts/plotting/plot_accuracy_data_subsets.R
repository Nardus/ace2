## Plot infection model accuracy on data subsets (i.e., different data qualities)

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


# ---- Read predictions ----------------------------------------------------------------------------
datasets <- c("all_data", "l2_data", "l1+2_data")

test_preds <- lapply(datasets, load_all_runs) %>% 
  bind_rows() %>% 
  filter(.data$response_var == "Infection")


infection_data <- readRDS("data/calculated/cleaned_infection_data.rds")


# ---- Final model data available ------------------------------------------------------------------
run_labels <- c("all_features" = "All ACE2 representations", 
                "phylogeny" = "Phylogenetic eigenvectors", 
                "all_features_phylogeny" = "All ACE2 representations +\nphylogenetic eigenvectors")

final_preds <- test_preds %>% 
  filter(.data$run_id %in% names(run_labels)) %>% 
  filter(.data$dataset == "all_data") %>% 
  left_join(infection_data, by = "species") %>% 
  mutate(run_label = factor(.data$run_id, 
                            levels = names(run_labels),
                            labels = run_labels),
         evidence_label = factor(.data$evidence_level,
                                 levels = rev(c("1", "2", "3", "4")),
                                 labels = rev(EVIDENCE_LABELS)))

data_available <- final_preds %>% 
  filter(.data$run_id == "all_features_phylogeny") %>%  # Same for all runs
  group_by(.data$evidence_label, .data$label) %>% 
  summarise(n = n(), .groups = "drop")

p_data <- ggplot(data_available, aes(x = evidence_label, y = n, fill = label)) +
  geom_col(colour = "grey20") +
  scale_y_continuous(expand = expansion(0.02)) +
  scale_fill_manual(values = INFECTION_STATUS_COLOURS) +
  labs(x = "Best evidence", y = "Number of species", fill = "Infected") +
  coord_flip() +
  theme(legend.position = "top")


# ---- Final model accuracy ------------------------------------------------------------------------
final_accuracies <- final_preds %>% 
  group_by(.data$evidence_label, .data$run_label) %>% 
  summarise(n_accurate = sum(.data$label == .data$prediction),
            n_total = n(),
            .groups = "keep") %>% 
  mutate(accuracy = .data$n_accurate/.data$n_total,
         lower = binom.test(.data$n_accurate, .data$n_total)$conf.int[1],
         upper = binom.test(.data$n_accurate, .data$n_total)$conf.int[2]) %>% 
  ungroup()


p_accuracy <- ggplot(final_accuracies, aes(x = evidence_label, y = accuracy, fill = run_label)) +
  geom_col(colour = "grey20", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.4, colour = "grey20",
                position = position_dodge(width = 0.9)) + 
  
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(0.02)) +
  scale_fill_brewer(palette = "Set2") +
  labs(y = "Accuracy", fill = "") +
  coord_flip() +
  theme(strip.text = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "top")


# ---- Final model sens/spec -----------------------------------------------------------------------
final_class_accuracies <- final_preds %>% 
  group_by(.data$label, .data$evidence_label, .data$run_label) %>% 
  summarise(n_accurate = sum(.data$label == .data$prediction),
            n_total = n(),
            .groups = "keep") %>% 
  mutate(accuracy = .data$n_accurate/.data$n_total,
         lower = binom.test(.data$n_accurate, .data$n_total)$conf.int[1],
         upper = binom.test(.data$n_accurate, .data$n_total)$conf.int[2]) %>% 
  ungroup()

# Sensitivity:
sens <- final_class_accuracies %>% 
  filter(.data$label == "True")

sens_missing <- data.frame(evidence_label = EVIDENCE_LABELS[!EVIDENCE_LABELS %in% sens$evidence_label],
                           accuracy = 0)

p_sens <- ggplot(sens, aes(x = evidence_label, y = accuracy)) +
  geom_col(aes(fill = run_label), colour = "grey20",
           position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = lower, ymax = upper, group = run_label), width = 0.4, colour = "grey20",
                position = position_dodge(width = 0.9)) + 
  
  geom_text(label = "N/A", hjust = 0, size = 2, data = sens_missing) +
  
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(0.02)) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(x = "Best evidence", y = "Sensitivity") +
  coord_flip() +
  theme(strip.text = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.y = element_blank())


# Specificity:
spec <- final_class_accuracies %>% 
  filter(.data$label == "False")

sens_missing <- data.frame(evidence_label = EVIDENCE_LABELS[!EVIDENCE_LABELS %in% spec$evidence_label],
                           accuracy = 0)

p_spec <- ggplot(spec, aes(x = evidence_label, y = accuracy)) +
  geom_col(aes(fill = run_label), colour = "grey20", 
           position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = lower, ymax = upper, group = run_label), width = 0.4, colour = "grey20", 
                position = position_dodge(width = 0.9)) + 
  
  geom_text(label = "N/A", hjust = 0, size = 2, data = sens_missing) +
  
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(0.02)) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(y = "Specificity") +
  coord_flip() +
  theme(strip.text = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.y = element_blank())


# ---- Data subset models (data available) ---------------------------------------------------------
dataset_labels <- c("l1+2_data" = "Observed and\nexperimental\ninfections",
                    "l2_data" = "Experimental\ninfection\nonly")

subset_preds <- test_preds %>% 
  filter(.data$run_id == "all_features_phylogeny") %>% 
  filter(.data$dataset != "all_data") %>% 
  left_join(infection_data, by = "species") %>% 
  mutate(dataset_label = factor(.data$dataset,
                                levels = rev(names(dataset_labels)),
                                labels = rev(dataset_labels)))

subset_data_available <- subset_preds %>% 
  group_by(.data$dataset_label, .data$label) %>% 
  summarise(n = n(), .groups = "drop")

p_data_subset <- ggplot(subset_data_available, aes(x = dataset_label, y = n, fill = label)) +
  geom_col(colour = "grey20") +
  scale_y_continuous(expand = expansion(0.02)) +
  scale_fill_manual(values = INFECTION_STATUS_COLOURS, guide = "none") +
  labs(x = "Dataset", y = "Number of species") +
  coord_flip()


# ---- Data subset models (accuracy) ---------------------------------------------------------------
subset_accuracies <- subset_preds %>% 
  group_by(.data$dataset_label) %>% 
  summarise(n_accurate = sum(.data$label == .data$prediction),
            n_total = n(),
            .groups = "keep") %>% 
  mutate(accuracy = .data$n_accurate/.data$n_total,
         lower = binom.test(.data$n_accurate, .data$n_total)$conf.int[1],
         upper = binom.test(.data$n_accurate, .data$n_total)$conf.int[2]) %>% 
  ungroup()


p_accuracy_subset <- ggplot(subset_accuracies, aes(x = dataset_label, y = accuracy)) +
  geom_col(colour = "grey20", fill = "grey50") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.4, colour = "grey20") + 
  
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(0.02)) +
  labs(y = "Accuracy") +
  coord_flip() +
  theme(strip.text = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.y = element_blank())


# ---- Data subset models (sens/spec) --------------------------------------------------------------
subset_class_accuracies <- subset_preds %>% 
  group_by(.data$label, .data$dataset_label) %>% 
  summarise(n_accurate = sum(.data$label == .data$prediction),
            n_total = n(),
            .groups = "keep") %>% 
  mutate(accuracy = .data$n_accurate/.data$n_total,
         lower = binom.test(.data$n_accurate, .data$n_total)$conf.int[1],
         upper = binom.test(.data$n_accurate, .data$n_total)$conf.int[2]) %>% 
  ungroup()


# Sensitivity:
sens_subset <- subset_class_accuracies %>% 
  filter(.data$label == "True")

p_sens_subset <- ggplot(sens_subset, aes(x = dataset_label, y = accuracy, fill = label)) +
  geom_col(colour = "grey20") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.4, colour = "grey20") + 
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(0.02)) +
  scale_fill_manual(values = INFECTION_STATUS_COLOURS, guide = "none") +
  labs(x = "Best evidence", y = "Sensitivity") +
  coord_flip() +
  theme(strip.text = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.y = element_blank())


# Specificity:
spec_subset <- subset_class_accuracies %>% 
  filter(.data$label == "False")

p_spec_subset <- ggplot(spec_subset, aes(x = dataset_label, y = accuracy, fill = label)) +
  geom_col(colour = "grey20") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.4, colour = "grey20") + 
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(0.02)) +
  scale_fill_manual(values = INFECTION_STATUS_COLOURS, guide = "none") +
  labs(y = "Specificity") +
  coord_flip() +
  theme(strip.text = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.y = element_blank())


# ---- Output --------------------------------------------------------------------------------------
# Legends (along top)
acc_legend <- get_legend(p_accuracy)
p_accuracy <- p_accuracy + guides(fill = "none")

dat_legend <- get_legend(p_data)
p_data <- p_data + guides(fill = "none")


# Plots
column_1 <- plot_grid(dat_legend, p_data, p_data_subset, 
                      ncol = 1, rel_heights = c(1, 4.5, 4.5),
                      align = "hv", axis = "tblr",
                      labels = c("A", "", "B"))

other_columns <- plot_grid(p_sens, p_spec, p_accuracy,
                           p_sens_subset, p_spec_subset, p_accuracy_subset,
                           nrow = 2, align = "h", axis = "tblr")

other_columns <- plot_grid(acc_legend, other_columns,
                           ncol = 1, rel_heights = c(1, 9))

p_final <- plot_grid(column_1, other_columns,
                     nrow = 1, rel_widths = c(1.2, 2))

ggsave2("output/plots/accuracy_data_subsets.pdf", p_final,
        width = 7, height = 5)


# ---- Values quoted in text -----------------------------------------------------------------------
cat("\nFinal model, sensitivity/specificity by evidence level:\n")
print(final_class_accuracies)


cat("\n\nSubset models:\n")
print(subset_accuracies)
