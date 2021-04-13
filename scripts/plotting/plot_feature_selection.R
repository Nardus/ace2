## Plot accuracy of feature selection runs

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(cowplot)
})

source("scripts/utils/output_utils.R")

# ---- Read all available predictions --------------------------------------------------------------

test_preds <- load_all_runs("all_data") %>% 
  bind_rows() %>% 
  filter(.data$run_id == "all_features" | startsWith(.data$run_id, "feature_selection"))


# ---- Clean labels and set order ------------------------------------------------------------------
plot_order <- c("All features", 
                sprintf("Feature selection\n(top %d)", 1:11))

test_preds <- test_preds %>% 
  mutate(feature_count = str_extract(.data$run_id, "[[:digit:]]+$"),
         run_name = if_else(.data$run_id == "all_features", "All features", 
                             sprintf("Feature selection\n(top %s)", .data$feature_count))) %>% 
  filter(.data$run_name %in% plot_order) %>% 
  mutate(run_name = factor(.data$run_name, levels = plot_order))


# ---- Accuracy by class ---------------------------------------------------------------------------
# Only makes sense for "all_data" runs
class_accuracies <- test_preds %>% 
  group_by(.data$dataset, .data$response_var, .data$run_name, .data$label) %>% 
  summarise(n = n(),
            accuracy = sum(.data$label == .data$prediction) / n(),
            .groups = "drop") %>% 
  mutate(n = sprintf("N = %d", .data$n))


p_class <- ggplot(class_accuracies, aes(x = run_name, y = accuracy)) +
  geom_col(colour = "grey20", position = "dodge") +
  geom_text(aes(label = n, y = accuracy + 0.05), size = 2, position = position_dodge(width = 1)) +
  labs(x = "feature set") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  facet_grid(rows = vars(response_var), cols = vars(label), scales = "free_x", space = "free_x") +
  theme_bw(base_size = 7) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggtitle("Accuracy by class (sensitivity/specificity)")



# ---- Output --------------------------------------------------------------------------------------
dir.create("output/plots", recursive = TRUE)
ggsave2("output/plots/feature_selection.png", 
        p_class,
        width = 6, height = 8)
