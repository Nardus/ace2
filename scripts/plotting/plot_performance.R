## Plot prediction accuracy / performance measures

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)
library(plotROC)
library(ModelMetrics)

predictions_inf <- read_rds("output/infection/predictions.rds") %>% 
  filter(.data$dataset == "test")

predictions_shed <- read_rds("output/shedding/predictions.rds") %>% 
  filter(.data$dataset == "test")

predictions_trans <- read_rds("output/transmission/predictions.rds") %>% 
  filter(.data$dataset == "test")


# Cutoffs (temp):
cutoffs_inf <- predictions_inf %>% 
  filter(.data$dataset == "train") %>% 
  group_by(.data$iteration) %>% 
  summarise(cutoff = sum(.data$label == "True") / n(),
            .groups = "drop")

# Combine test preds
test_preds <- list(Infection = predictions_inf,
                   Shedding = predictions_shed,
                   Transmission = predictions_trans) %>% 
  bind_rows(.id = "dataset_name")


# ---- AUC ----------------------------------------------------------------------------------------
aucs <- test_preds %>% 
  group_by(.data$dataset_name, .data$iteration) %>% 
  summarise(auc = auc(actual = .data$label == "True", predicted = .data$prob),
            .groups = "drop")
  
p_auc <- ggplot(aucs, aes(x = dataset_name, y = auc)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey60") +
  geom_violin(fill = "grey60", colour = NA) +
  geom_boxplot(width = 0.5, colour = "grey20") +
  labs(x = "response") +
  ylim(c(0, 1)) +
  theme_bw()


# ---- Accuracy by class / evidence level ---------------------------------------------------------
accuracies <- test_preds %>% 
  left_join(cutoffs_inf, by = "iteration") %>%   # TEMP
  mutate(prediction = if_else(.data$prob > .data$cutoff, "True", "False")) %>% 
  
  group_by(.data$dataset_name, .data$iteration, .data$label, .data$evidence_level) %>% 
  summarise(accuracy = sum(.data$label == .data$prediction) / n(),
            .groups = "drop")

p_accuracy <- ggplot(accuracies, aes(x = factor(evidence_level), y = accuracy)) +
  geom_jitter(width = 0.2, height = 0.01, colour = "grey20") +
  #geom_violin(fill = "grey60", colour = NA) +
  #geom_boxplot(width = 0.05, colour = "grey20") +
  labs(x = "evidence level") +
  facet_grid(rows = vars(dataset_name), cols = vars(label)) +
  theme_bw()


# ---- ROC curves ---------------------------------------------------------------------------------
# ggplot(test_preds, aes(d = label, m = prob)) +
#   geom_roc(aes(fill = factor(iteration)), n.cuts = 0, colour = "grey40", size = 0.2) +
#   geom_abline(colour = "grey20", linetype = 2) +
#   
#   guides(fill = FALSE) +
#   coord_equal() +
#   scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
#   scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
#   labs(x = 'proportion false positives', 
#        y = 'proportion true positives') +
#   facet_wrap(vars(dataset_name)) +
#   theme_bw() +
#   theme(plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 10.5))


# ---- Combine / save -----------------------------------------------------------------------------
p_final <- plot_grid(p_auc, p_accuracy,
                     ncol = 2, rel_widths = c(1, 2),
                     align = "h", axis = "tb")

dir.create("output/plots", recursive = TRUE)
ggsave2("output/plots/performance.png", 
        p_final,
        width = 5, height = 3)
