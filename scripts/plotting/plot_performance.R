## Plot prediction accuracy / performance measures

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(cowplot)
  library(ModelMetrics)
})

source("scripts/utils/output_utils.R")

# ---- Read all available predictions --------------------------------------------------------------

datasets <- c("all_data", "l2_data", "l3_data", "l1+2_data")

test_preds <- lapply(datasets, load_all_runs) %>% 
  bind_rows()


# ---- Plotting order -----------------------------------------------------------------------------
# TODO: clean labels for each run_id...

RUN_ORDER <- c("combined", "combined+feature_selection_10")

test_preds <- test_preds %>% 
  mutate(run_id = factor(.data$run_id, levels = RUN_ORDER))


# ---- Overall accuracy ---------------------------------------------------------------------------
overall_accuracies <- test_preds %>% 
  group_by(.data$dataset, .data$response_var, .data$run_id) %>% 
  summarise(n = n(),
            accuracy = sum(.data$label == .data$prediction) / n(),
            .groups = "drop") %>% 
  mutate(n = sprintf("N = %d", .data$n))
  
p_overall <- ggplot(overall_accuracies, aes(x = run_id, y = accuracy, fill = run_id)) +
  geom_col(colour = "grey20") +
  geom_text(aes(label = n), nudge_y = 0.04, size = 1) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey10") +
  labs(x = "feature set") +
  ylim(c(0, 1)) +
  facet_grid(cols = vars(dataset), rows = vars(response_var), scales = "free_x", space = "free_x") +
  guides(fill = FALSE) +
  theme_bw(base_size = 5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggtitle("Overall accuracy")


# ---- Accuracy by class only ---------------------------------------------------------------------
# Only makes sense for "all_data" runs
class_accuracies <- test_preds %>% 
  group_by(.data$dataset, .data$response_var, .data$run_id, .data$label) %>% 
  summarise(n = n(),
            accuracy = sum(.data$label == .data$prediction) / n(),
            .groups = "drop") %>% 
  mutate(n = sprintf("N = %d", .data$n))


p_class <- ggplot(class_accuracies, aes(x = factor(run_id), y = accuracy, fill = label)) +
  geom_col(colour = "grey20", position = "dodge") +
  geom_text(aes(label = n, y = accuracy + 0.05), size = 1, position = position_dodge(width = 1)) +
  labs(x = "feature set") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  facet_grid(rows = vars(response_var), cols = vars(dataset), scales = "free_x", space = "free_x") +
  theme_bw(base_size = 5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.key.size = unit(0.5, "lines")) +
  ggtitle("Accuracy by class (sensitivity/specificity)")


# ---- Accuracy by class and evidence level ---------------------------------------------------------
# Only makes sense for "all_data" and "l1+2" runs
level_accuracies <- test_preds %>% 
  group_by(.data$dataset, .data$response_var, .data$run_id, .data$label, .data$evidence_level) %>% 
  summarise(n = n(),
            accuracy = sum(.data$label == .data$prediction) / n(),
            .groups = "drop") %>% 
  mutate(n = sprintf("N = %d", .data$n))

label_data_level <- level_accuracies %>% 
  group_by(.data$dataset, .data$response_var, .data$label, .data$evidence_level, .data$n) %>% 
  summarise(accuracy = max(.data$accuracy), .groups = "drop")

make_level_plot <- function(data_name, acc_data = level_accuracies, lab_data = label_data_level) {
  acc_data <- acc_data %>% 
    filter(.data$dataset == data_name)
  
  lab_data <- lab_data %>% 
    filter(.data$dataset == data_name)
  
  ggplot(acc_data, aes(x = factor(evidence_level), y = accuracy)) +
    geom_col(aes(fill = run_id), colour = "grey20", position = "dodge") +
    geom_text(aes(label = n), nudge_y = 0.15, size = 1, data = lab_data) +
    labs(x = "evidence level") +
    scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
    facet_grid(rows = vars(response_var), cols = vars(label)) +
    guides(fill = FALSE) +
    theme_bw(base_size = 5) +
    ggtitle(sprintf("Accuracy by evidence level\n(%s model)", data_name))
}

p_level_all <- make_level_plot("all_data")
p_level_l1_2 <- make_level_plot("l1+2_data")


# ---- Combine / save -----------------------------------------------------------------------------
p_final <- plot_grid(p_overall, p_class, p_level_l1_2, p_level_all,
                     rel_widths = c(1, 1.5),
                     nrow = 2, align = "h", axis = "tb")

dir.create("output/plots", recursive = TRUE)
ggsave2("output/plots/performance.png", 
        p_final,
        width = 6, height = 8)




# ---- Stats to quote in text ---------------------------------------------------------------------
# Does including cell culture data make the model less reliable?
test_preds %>% 
  filter(.data$dataset %in% c("all_data", "l1+2_data")) %>% 
  filter(.data$response_var == "Infection") %>% 
  filter(.data$run_id == "combined") %>% 
  group_by(.data$dataset, .data$response_var, .data$run_id, .data$label) %>% 
  summarise(n = n(),
            accuracy = sum(.data$label == .data$prediction) / n(),
            .groups = "drop")


test_preds %>% 
  filter(.data$dataset %in% c("all_data", "l1+2_data")) %>% 
  filter(.data$response_var == "Infection") %>% 
  filter(.data$run_id == "combined") %>% 
  group_by(.data$dataset, .data$response_var, .data$run_id, .data$label, .data$evidence_level) %>% 
  summarise(n = n(),
            accuracy = sum(.data$label == .data$prediction) / n(),
            .groups = "drop")
