## Plot prediction accuracy / performance measures

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)
library(plotROC)
library(ModelMetrics)

predictions_inf <- read_rds("output/infection/predictions.rds")
predictions_shed <- read_rds("output/shedding/predictions.rds")
predictions_trans <- read_rds("output/transmission/predictions.rds")


# Combine test preds
test_preds <- list(Infection = predictions_inf,
                   Shedding = predictions_shed,
                   Transmission = predictions_trans) %>% 
  bind_rows(.id = "dataset_name")


# ---- AUC ----------------------------------------------------------------------------------------
aucs <- test_preds %>% 
  group_by(.data$dataset_name) %>% 
  summarise(n = n(),
            auc = auc(actual = .data$label == "True", predicted = .data$prob),
            .groups = "drop") %>% 
  mutate(n = sprintf("N = %d", .data$n))
  
p_auc <- ggplot(aucs, aes(x = dataset_name, y = auc)) +
  geom_col(colour = "grey20", fill = "grey60") +
  geom_text(aes(label = n), nudge_y = 0.04, size = 2) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey10") +
  labs(x = "response") +
  ylim(c(0, 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


# ---- Confusion matrix --------------------------------------------------------------------------
confusion <- test_preds %>% 
  group_by(.data$dataset_name, .data$label) %>% 
  mutate(lab_count = n()) %>% 
  group_by(.data$dataset_name, .data$label, .data$prediction) %>% 
  summarise(lab_count = unique(.data$lab_count),
            lab_pred_count = n(),
            .groups = "drop") %>% 
  mutate(prop = .data$lab_pred_count / .data$lab_count,
         plot_label = sprintf("%d/%d\n(%3.1f%%)", .data$lab_pred_count, .data$lab_count, .data$prop*100))

p_confusion <- ggplot(confusion, aes(x = label, y = prediction, fill = prop)) +
  geom_tile(colour = "grey20", size = 0.5) +
  geom_text(aes(label = plot_label), size = 2) +
  scale_fill_distiller(palette = "Blues", direction = 1, guide = FALSE) +
  labs(x = "observed", y = "predicted") +
  facet_grid(rows = vars(dataset_name)) +
  theme_bw()


# ---- Accuracy by class / evidence level ---------------------------------------------------------
accuracies <- test_preds %>% 
  group_by(.data$dataset_name, .data$label, .data$evidence_level) %>% 
  summarise(n = n(),
            accuracy = sum(.data$label == .data$prediction) / n(),
            .groups = "drop") %>% 
  mutate(n = sprintf("N = %d", .data$n))

p_accuracy <- ggplot(accuracies, aes(x = factor(evidence_level), y = accuracy)) +
  geom_col(colour = "grey20", fill = "grey60") +
  geom_text(aes(label = n), nudge_y = 0.15, size = 2) +
  labs(x = "evidence level") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  facet_grid(rows = vars(dataset_name), cols = vars(label)) +
  theme_bw()


# ---- Combine / save -----------------------------------------------------------------------------
p_final <- plot_grid(p_auc, p_confusion, p_accuracy,
                     ncol = 3, rel_widths = c(1, 1.2, 2),
                     align = "h", axis = "t")

dir.create("output/plots", recursive = TRUE)
ggsave2("output/plots/performance.png", 
        p_final,
        width = 7, height = 3.5)
