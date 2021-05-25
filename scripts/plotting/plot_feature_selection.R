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
all_feature_location <- 15 + 5
plot_order <- c(sprintf("%d", 1:15),
                "All features")

test_preds <- test_preds %>% 
  mutate(feature_count = str_extract(.data$run_id, "[[:digit:]]+$"),
         feature_count = if_else(is.na(.data$feature_count), all_feature_location, as.numeric(.data$feature_count)),
         run_name = if_else(.data$run_id == "all_features", "All features", 
                             sprintf("%s", .data$feature_count))) %>% 
  filter(.data$run_name %in% plot_order) %>% 
  mutate(run_name = factor(.data$run_name, levels = plot_order))


# ---- Accuracy by class ---------------------------------------------------------------------------
# Only makes sense for "all_data" runs
class_accuracies <- test_preds %>% 
  group_by(.data$dataset, .data$response_var, .data$run_name, .data$iteration, 
           .data$label, .data$feature_count) %>% 
  summarise(n = n(),
            accuracy = sum(.data$label == .data$prediction) / n(),
            .groups = "drop")


p_class <- ggplot(class_accuracies[class_accuracies$feature_count != all_feature_location, ], 
                  aes(x = feature_count, y = accuracy, colour = label)) +
  geom_blank() +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey20") +
  geom_point() +
  geom_line(aes(group = label)) +
  geom_jitter(data = class_accuracies[class_accuracies$feature_count == all_feature_location, ]) +
  
  scale_x_continuous(breaks = unique(class_accuracies$feature_count),
                     labels = plot_order) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0, 1)) +
  scale_colour_brewer(name = NULL, palette = "Set1", 
                      labels = c("True" = "Sensitivity", "False" = "Specificity")) +
  facet_grid(rows = vars(response_var), scales = "free_x", space = "free_x") +
  theme_bw(base_size = 7) +
  labs(title = "Feature selection", x = "Number of features", y = "Proportion accurate")



# ---- Output --------------------------------------------------------------------------------------
dir.create("output/plots", recursive = TRUE)
ggsave2("output/plots/feature_selection.png", 
        p_class,
        width = 6, height = 8)
