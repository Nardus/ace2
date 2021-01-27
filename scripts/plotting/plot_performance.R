## Plot prediction accuracy / performance measures

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(cowplot)
  library(plotROC)
  library(ModelMetrics)
})


# ---- Read all available predictions --------------------------------------------------------------
load_single_run <- function(response_var, run_id) {
  file.path("output", response_var, run_id, "predictions.rds") %>% 
    read_rds() %>% 
    mutate(response_var = str_to_sentence(response_var),
           run_id = run_id)
}

load_all_runs <- function(response_var) {
  top_dir <- file.path("output", response_var)
  run_ids <- list.dirs(top_dir, full.names = FALSE, recursive = FALSE)
  
  lapply(run_ids, load_single_run, response_var = response_var) %>% 
    bind_rows()
}

responses <- c("infection", "shedding", "transmission")

test_preds <- lapply(responses, load_all_runs) %>% 
  bind_rows()

# TODO: clean labels for each run_id...


# ---- AUC ----------------------------------------------------------------------------------------
aucs <- test_preds %>% 
  group_by(.data$response_var, .data$run_id) %>% 
  summarise(n = n(),
            auc = auc(actual = .data$label == "True", predicted = .data$prob),
            .groups = "drop") %>% 
  mutate(n = sprintf("N = %d", .data$n))
  
p_auc <- ggplot(aucs, aes(x = run_id, y = auc, fill = run_id)) +
  geom_col(colour = "grey20") +
  geom_text(aes(label = n), nudge_y = 0.04, size = 2) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey10") +
  labs(x = "response") +
  ylim(c(0, 1)) +
  facet_grid(rows = vars(response_var)) +
  guides(fill = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))



# ---- Accuracy by class / evidence level ---------------------------------------------------------
accuracies <- test_preds %>% 
  group_by(.data$response_var, .data$run_id, .data$label, .data$evidence_level) %>% 
  summarise(n = n(),
            accuracy = sum(.data$label == .data$prediction) / n(),
            .groups = "drop") %>% 
  mutate(n = sprintf("N = %d", .data$n))

label_data <- accuracies %>% 
  group_by(.data$response_var, .data$label, .data$evidence_level, .data$n) %>% 
  summarise(accuracy = max(.data$accuracy), .groups = "drop")
  

p_accuracy <- ggplot(accuracies, aes(x = factor(evidence_level), y = accuracy)) +
  geom_col(aes(fill = run_id), colour = "grey20", position = "dodge") +
  geom_text(aes(label = n), nudge_y = 0.15, size = 2, data = label_data) +
  labs(x = "evidence level") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  facet_grid(rows = vars(response_var), cols = vars(label)) +
  guides(fill = FALSE) +
  theme_bw()


# ---- Combine / save -----------------------------------------------------------------------------
p_final <- plot_grid(p_auc, p_accuracy,
                     ncol = 2, rel_widths = c(1, 2),
                     align = "h", axis = "t")

dir.create("output/plots", recursive = TRUE)
ggsave2("output/plots/performance.png", 
        p_final,
        width = 7, height = 3.5)
