## Utility functions for finding cutoffs

# Test a single cutoff
test_cutoff <- function(preds, cutoff) {
  preds %>% 
    mutate(new_prediction = if_else(.data$probability > cutoff, "True", "False"),
           new_prediction = factor(new_prediction, levels = c("True", "False"))) %>% 
    group_by(.data$label, .add = TRUE) %>% 
    summarise(class_acc = accuracy_vec(truth = .data$label, estimate = .data$new_prediction),
              .groups = "drop_last") %>% 
    summarise(balanced_accuracy = sum(.data$class_acc)/2,
              .groups = "drop") %>% 
    mutate(cutoff = cutoff)
}

# Evaluate a range of cutoffs
find_best_cuttof <- function(predictions, cutoffs = seq(0.15, 0.95, by = 0.01)) {
  train_preds <- predictions %>% 
    filter(!is.na(.data$label))
  
  best_cutoff <- sapply(cutoffs, test_cutoff, preds = train_preds, simplify = FALSE) %>% 
    bind_rows() %>% 
    slice_max(.data$balanced_accuracy, with_ties = TRUE) %>% 
    slice_sample(n = 1) %>%  # Choose randomly in case of ties
    pull(.data$cutoff)
  
  best_cutoff
}
