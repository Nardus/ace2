## Compare model accuracy across different sarbecoviruses
#   - Supplement to "plot_accuracy.R"

suppressPackageStartupMessages({
    library(argparse)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(cowplot)
    
    source("scripts/plotting/plotting_constants.R")
})

parser <- ArgumentParser(description = "Compare performance on non-ACE2 sarbecoviruses")

parser$add_argument("predictions", type = "character", 
                    help = "path to a single model's predictions")

parser$add_argument("output_name", type = "character",
                    help = "location/name for output plot")

INPUT <- parser$parse_args()


# ---- Data ----------------------------------------------------------------------------------------
infection_data <- readRDS("data/calculated/cleaned_infection_data.rds")

predictions <- readRDS(INPUT$predictions)

# Extract virus info
# - Using the viruses matches the final label used in training
viruses <- c("SARS-like CoV", "ACE2-utilizing SARS-like CoV", "SARS-CoV-2", "SARS-CoV-1")

virus_data <- infection_data %>% 
    mutate(viruses = if_else(.data$infected == "True",
                             .data$viruses_true,
                             .data$viruses_false)) %>%
    select(.data$species, .data$viruses)


non_ace2_data <- virus_data %>% 
    mutate(viruses = if_else(.data$viruses == "SARS-like CoV", 
                             "Linked only to\nsarbecoviruses\nsuspected to\nnot use ACE2", 
                             "Linked to at least\none sarbecovirus\nknown to use ACE2"))

viruses <- c("Linked only to\nsarbecoviruses\nsuspected to\nnot use ACE2", 
             "Linked to at least\none sarbecovirus\nknown to use ACE2")



# ---- Accuracy for species linked ONLY to SARS-like CoVs not known to use ACE2 --------------------
get_accuracy <- function(virus, labels = c("True", "False"), preds = predictions, vir_data = non_ace2_data) {
    preds <- preds %>% 
        left_join(vir_data, by = "species") %>%
    
        filter(.data$label %in% labels) %>% 
        filter(grepl(virus, .data$viruses, fixed = TRUE)) %>%
        ungroup()
        
    if (nrow(preds) == 0) {
        empty <- data.frame(n_total = 0,
                            accuracy = NA,
                            lower = NA,
                            upper = NA,
                            virus = virus)
        return(empty)   
    }
    
    preds %>% 
        summarise(n_accurate = sum(.data$label == .data$prediction),
                  n_total = n()) %>% 
        mutate(accuracy = .data$n_accurate/.data$n_total,
                lower = binom.test(.data$n_accurate, .data$n_total)$conf.int[1],
                upper = binom.test(.data$n_accurate, .data$n_total)$conf.int[2],
                virus = virus)
}

    
accuracy <- lapply(viruses, get_accuracy, vir_data = non_ace2_data) %>% 
    bind_rows() %>% 
    mutate(virus = factor(.data$virus, levels = viruses))
    
sens <- lapply(viruses, get_accuracy, labels = "True", vir_data = non_ace2_data) %>% 
    bind_rows() %>% 
    mutate(virus = factor(.data$virus, levels = viruses))
    
spec <- lapply(viruses, get_accuracy, labels = "False", vir_data = non_ace2_data) %>% 
    bind_rows() %>% 
    mutate(virus = factor(.data$virus, levels = viruses))


# ---- Plot performance ----------------------------------------------------------------------------
plot_performance <- function(acc_df, label, colour) {
    na_labels <- acc_df %>% 
        filter(is.na(.data$accuracy)) %>% 
        mutate(label = "N/A")
    
    p <- ggplot(acc_df, aes(x = virus, y = accuracy)) +
        geom_col(colour = "grey20", fill = colour) +
        geom_errorbar(aes(ymin = lower, ymax = upper), colour = "grey20", width = 0.4) + 
        
        scale_y_continuous(limits = c(0, 1), expand = expansion(0.02)) +
        labs(y = label) +
        coord_flip() +
        theme(axis.text.y = element_blank(),
              axis.title.y = element_blank())
              
    if (!any(is.na(acc_df$accuracy))) {
        return(p)
    }
    
    p + geom_text(aes(label = label, y = 0), hjust = 0, size = 2, data = na_labels)
}

p_acc <- plot_performance(accuracy, "Accuracy", colour = "grey50")
p_sens <- plot_performance(sens, "Sensitivity", colour = INFECTION_STATUS_COLOURS["True"])
p_spec <- plot_performance(spec, "Specificity", colour = INFECTION_STATUS_COLOURS["False"])


# ---- Plot sample size ----------------------------------------------------------------------------
group_sizes <- infection_data %>% 
    left_join(non_ace2_data, by = "species") %>%
    rename(virus = .data$viruses) %>% 
    group_by(.data$virus, .data$infected) %>% 
    summarise(n_total = n(), .groups = "drop") %>% 
    mutate(virus = factor(.data$virus, levels = viruses))

p_data <- ggplot(group_sizes, aes(x = virus, y = n_total, fill = infected)) +
    geom_col(colour = "grey20") +
    scale_y_continuous(expand = expansion(0.02)) +
    scale_fill_manual(values = INFECTION_STATUS_COLOURS, guide = "none") +
    labs(x = "Species group", y = "Number of species") +
    coord_flip()


# ---- Combine -------------------------------------------------------------------------------------
combined_plot <- plot_grid(p_data, p_acc, p_sens, p_spec,
                           nrow = 1, rel_widths = c(2, 1, 1, 1))


ggsave2(INPUT$output_name, combined_plot, width = 7, height = 2)
