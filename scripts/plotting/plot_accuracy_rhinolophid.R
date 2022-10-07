## Plot infection model accuracy on data subsets (including/excluding rhinolophid bats)

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
datasets <- c("all_data", "no_rhinolophids")

test_preds <- lapply(datasets, load_all_runs) %>% 
  bind_rows() %>% 
  filter(.data$response_var == "Infection")


infection_data <- readRDS("data/calculated/cleaned_infection_data.rds")


# ---- Human-readable labels -----------------------------------------------------------------------
DATA_LABELS <- c( "no_rhinolophids" = "Rhinolophid bats\nexcluded",
                  "all_data" = "All data")

RUN_LABELS <- c("all_features" = "All ACE2 representations", 
                "phylogeny" = "Phylogenetic eigenvectors", 
                "all_features_phylogeny" = "All ACE2 representations +\nphylogenetic eigenvectors")

test_preds <- test_preds %>% 
  filter(.data$run_id %in% names(RUN_LABELS)) %>% 
  left_join(infection_data, by = "species") %>% 
  mutate(data_label = factor(.data$dataset,
                             levels = names(DATA_LABELS),
                             labels = DATA_LABELS),
         run_label = factor(.data$run_id, 
                            levels = names(RUN_LABELS),
                            labels = RUN_LABELS))


# ---- Data available ------------------------------------------------------------------------------
subset_preds <- test_preds %>% 
  filter(.data$dataset == "no_rhinolophids")

data_available <- test_preds %>% 
  filter(.data$run_id == "all_features_phylogeny") %>%  # Same data for all runs
  group_by(.data$data_label, .data$run_label, .data$label) %>% 
  summarise(n = n(), .groups = "drop")

p_data <- ggplot(data_available, aes(x = data_label, y = n, fill = label)) +
  geom_col(colour = "grey20") +
  scale_y_continuous(expand = expansion(0.02)) +
  scale_fill_manual(values = INFECTION_STATUS_COLOURS) +
  labs(x = "Dataset", y = "Number of species", fill = "Infected") +
  coord_flip() +
  theme(legend.position = "top")


# ---- Accuracy ------------------------------------------------------------------------------------
final_accuracies <- test_preds %>% 
  group_by(.data$data_label, .data$run_label) %>% 
  summarise(n_accurate = sum(.data$label == .data$prediction),
            n_total = n(),
            .groups = "keep") %>% 
  mutate(accuracy = .data$n_accurate/.data$n_total,
         lower = binom.test(.data$n_accurate, .data$n_total)$conf.int[1],
         upper = binom.test(.data$n_accurate, .data$n_total)$conf.int[2]) %>% 
  ungroup()


p_accuracy <- ggplot(final_accuracies, aes(x = data_label, y = accuracy, fill=run_label)) +
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


# ---- Sens/spec -----------------------------------------------------------------------------------
class_accuracies <- test_preds %>% 
  group_by(.data$data_label, .data$label, .data$run_label) %>% 
  summarise(n_accurate = sum(.data$label == .data$prediction),
            n_total = n(),
            .groups = "keep") %>% 
  mutate(accuracy = .data$n_accurate/.data$n_total,
         lower = binom.test(.data$n_accurate, .data$n_total)$conf.int[1],
         upper = binom.test(.data$n_accurate, .data$n_total)$conf.int[2]) %>% 
  ungroup()


# Sensitivity:
sens <- class_accuracies %>% 
  filter(.data$label == "True")

p_sens <- ggplot(sens, aes(x = data_label, y = accuracy, fill = run_label)) +
  geom_col(colour = "grey20", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.4, colour = "grey20", 
                position = position_dodge(width = 0.9)) + 
  
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(0.02)) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(y = "Sensitivity") +
  coord_flip() +
  theme(strip.text = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.y = element_blank())


# Specificity:
spec <- class_accuracies %>% 
  filter(.data$label == "False")

p_spec <- ggplot(spec, aes(x = data_label, y = accuracy, fill = run_label)) +
  geom_col(colour = "grey20", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.4, colour = "grey20", 
                position = position_dodge(width = 0.9)) + 
  
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(0.02)) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
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
column_1 <- plot_grid(dat_legend, p_data, 
                      ncol = 1, rel_heights = c(1, 9),
                      align = "hv", axis = "tblr")

other_columns <- plot_grid(p_sens, p_spec, p_accuracy,
                           nrow = 1,
                           align = "h", axis = "tblr")

other_columns <- plot_grid(acc_legend, other_columns,
                           ncol = 1, rel_heights = c(1, 9))

p_final <- plot_grid(column_1, other_columns,
                     nrow = 1, rel_widths = c(1.2, 2))

ggsave2("output/plots/accuracy_rhinolophid.pdf", p_final,
        width = 7, height = 2)


# ---- Values quoted in text -----------------------------------------------------------------------
cat("\nFinal model, sensitivity/specificity by evidence level:\n")
print(final_class_accuracies)


cat("\n\nSubset models:\n")
print(subset_accuracies)
print(subset_class_accuracies)
