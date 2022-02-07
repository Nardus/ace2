## Utility functions for plotting

library(stringr)

# ---- Human-readable feature names ----------------------------------------------------------------
# Expects a data.frame with a column named "feature"
add_readable_feature_names <- function(x) {
  x %>% 
    mutate(feature_type = case_when(startsWith(.data$feature, "dist_") ~ "Consensus distance",
                                    startsWith(.data$feature, "variable_site") ~ "Amino acid identity",
                                    startsWith(.data$feature, "property_polarity") ~ "Polarity",
                                    startsWith(.data$feature, "property_hydrophobicity") ~ "Hydrophobicity",
                                    startsWith(.data$feature, "property_volume") ~ "Volume",
                                    startsWith(.data$feature, "closest_positive") ~ "Minimum distance",
                                    .data$feature %in% c("haddock_score", "huang_score") ~ "Binding affinity",
                                    TRUE ~ "Other"),
           extra_info = case_when(.data$feature == "closest_positive_overall" ~ "all data",
                                  .data$feature == "closest_positive_l1" ~ "observed infections",
                                  .data$feature == "closest_positive_l2" ~ "experimental infections",
                                  .data$feature == "closest_positive_l3" ~ "cell culture",
                                  .data$feature == "closest_positive_l4" ~ "het-ACE2",
                                  .data$feature == "haddock_score" ~ "Fischhoff et al.",
                                  .data$feature == "huang_score" ~ "Huang et al.",
                                  TRUE ~ NA_character_),
           feature_position = if_else(.data$feature_type %in% c("Consensus distance", "Amino acid identity",
                                                                "Polarity", "Hydrophobicity", "Volume"),
                                      str_extract(.data$feature, "[[:digit:]]+$"), 
                                      NA_character_),
           feature_position = as.integer(.data$feature_position),
           feature_position_corrected = as_human_coord_v(.data$feature_position),
           feature_label = case_when(.data$feature_type == "Other" ~ .data$feature,
                                     !is.na(.data$extra_info) ~ sprintf("%s\n(%s)",
                                                                        .data$feature_type,
                                                                        .data$extra_info),
                                     !is.na(.data$feature_position) ~ sprintf("%s (%d)", 
                                                                              .data$feature_type, 
                                                                              .data$feature_position_corrected), 
                                     TRUE ~ .data$feature_type))
}


# ---- Jittering -----------------------------------------------------------------------------------
# Adjust jitter to some percentage of the size of the smallest gap in x
# Based on the default calculation in jitter(), but allows changing jitter_percent
get_jitter_amount <- function(x, jitter_percent = 0.20) {
  z <- diff(r <- range(x[is.finite(x)]))
  
  if (z == 0) 
    z <- abs(r[1L])
  if (z == 0) 
    z <- 1
  
  d <- diff(xx <- unique(sort.int(round(x, 3 - floor(log10(z))))))
  d <- if (length(d)) 
    min(d)
  else if (xx != 0) 
    xx/10
  else z/10
  
  jitter_percent * abs(d)
}


# ---- SHAP values ---------------------------------------------------------------------------------

get_shap_values <- function(training_result, iteration = 1, interaction = FALSE) {
  training_data <- training_result$trainingData
  trained_model <- training_result$trained_models[[iteration]]$finalModel
  
  metadata <- training_data %>% 
    select(.data$species, .data$ace2_accession, .data$evidence_level, .data$label)
  
  features <- training_data %>% 
    select(-.data$species, -.data$ace2_accession, -.data$evidence_level, -.data$label)
  
  # Calculate SHAP values
  dummy_encoding <- dummyVars(~ ., data = features)
  feature_mat <- predict(dummy_encoding, features)
  
  colnames(feature_mat) <- str_remove(colnames(feature_mat), "\\.")
  feature_mat <- feature_mat[, trained_model$feature_names]
  
  shap_values <- predict(trained_model, 
                         newdata = feature_mat,
                         predcontrib = !interaction,  # If SHAP interaction vals not requested, get regular SHAP vals
                         predinteraction = interaction)
  
  if (!interaction) {
    # shap_values is a matrix
    shap_values <- shap_values %>% 
      as.data.frame() %>% 
      rownames_to_column("rowid") %>% 
      select(-.data$BIAS) %>% 
      pivot_longer(-.data$rowid, names_to = "feature", values_to = "shap_value")
    
  } else {
    # shap_values is an array with format x[data_row][feature_1][feature_2]
    shap_values <- shap_values %>% 
      as.data.frame() %>% 
      rownames_to_column("rowid") %>% 
      pivot_longer(-.data$rowid, names_to = "feature", values_to = "shap_value") %>% 
      separate(.data$feature, into = c("feature_1", "feature_2"), sep = "\\.") %>% 
      filter(.data$feature_1 != "BIAS", .data$feature_2 != "BIAS")
  }
  
  
  # Add original data
  metadata <- metadata %>% 
    rownames_to_column("rowid")
  
  feature_vals <- feature_mat %>% 
    as.data.frame() %>% 
    rownames_to_column("rowid") %>% 
    pivot_longer(-.data$rowid, names_to = "feature", values_to = "feature_value")
  
  if (!interaction) {
    shap_values <- shap_values %>% 
      left_join(feature_vals, by = c("rowid", "feature"))
  } else {
    feature_vals_1 <- rename(feature_vals, feature_1_value = .data$feature_value)
    feature_vals_2 <- rename(feature_vals, feature_2_value = .data$feature_value)
    
    shap_values <- shap_values %>% 
      left_join(feature_vals_1, by = c("rowid", "feature_1" = "feature")) %>% 
      left_join(feature_vals_2, by = c("rowid", "feature_2" = "feature"))
  }
  
  shap_values %>% 
    left_join(metadata, by = "rowid")
}


# ---- Plot annotations -------------------------------------------------------------------------
# Return spearman correlation and R-squared in a format suitable for geom_text annotation
trendline_annotation <- function(x, y, x_name, y_name, data, ...) {
  data$x <- data[[x_name]]
  data$y <- data[[y_name]]
  
  spearman_cor <- cor(data$x, data$y, method = "spearman")
  fit <- lm(y ~ x, data = data)
  fit_summary <- summary(fit)
  
  p_val <- fit_summary$coefficients[2, "Pr(>|t|)"]
  p_val <- if_else(p_val < 0.001, "<0.001", sprintf("==%3.3f", p_val))
  
  df <- data.frame(label = c(sprintf("textstyle(Spearman)~rho==%3.3f", spearman_cor),
                             sprintf("list(R^2==%3.3f, p%s)", fit_summary$r.squared, p_val))) 
  
  geom_text(aes(label = label), x = x, y = y, parse = TRUE, data = df, ...)
}


# ---- Specific plot types -------------------------------------------------------------------------
# Plot a distance measure, to illustrate separation of classes
plot_dists <- function(y_var, test_position, distances, viruses) {
  evidence_labels = c("1" = "Observed infection",
                      "2" = "Experimental infection",
                      "3" = "Cell culture",
                      "4" = "Cell culture (het-ACE2)")
  
  # Jitter to approximate local density
  distances$adjusted_x = as.numeric(distances$infected) + offsetX(distances[[y_var]], distances$infected)
  
  # Virus symbols
  #  - This duplicates points, but with the same jitter position 
  #    (so we can plot all viruses at a given species prediction)
  viruses <- left_join(viruses, distances, by = "species")  
  
  # Plots
  ggplot(distances, aes_string(x = "infected", y = y_var)) + # Boxplots/p-value using raw data (before duplication)
    geom_boxplot(outlier.colour = NA) + 
    
    geom_point(aes(x = adjusted_x, fill = factor(evidence_level), shape = virus), 
               colour = "grey40", size = 1, stroke = 0.3, data = viruses) +  # Use repeated, overlapping symbols
    
    geom_signif(comparisons = list(c("True", "False")), 
                test = "wilcox.test", test.args = list(exact = FALSE),
                size = 0.3, tip_length = 0.014, textsize = 1.8, vjust = -0.2, 
                y_position = test_position) +
    
    scale_fill_brewer(palette = "YlGnBu", direction = - 1, labels = evidence_labels) +
    scale_shape_manual(values = VIRUS_SHAPES) +
    theme(legend.position = "none")
}


# Raw SHAP values vs input feature value
plot_shap <- function(shapvals, feature_name, x_is_factor = FALSE, jitter_width = 0.1, jitter_height = 0.05) {
  shapvals <- shapvals %>% 
    filter(.data$feature == feature_name)
  
  if (x_is_factor)
    shapvals$feature_value = as.factor(shapvals$feature_value)
  
  ggplot(shapvals, aes(x = feature_value, y = shap_value, colour = label, shape = label)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_jitter(width = jitter_width, height = jitter_height, size = 1.8) + 
    scale_colour_brewer(palette = "Set1", na.value = "grey80") +
    scale_shape_manual(values = c(19, 1)) +
    PLOT_THEME
}


# Raw shap values, coloured by input values of a potential interacting feature
plot_shap_interaction <- function(shap_values, main_feature, interacting_feature, 
                                  jitter_width = 0.1, jitter_height = 0.05) {
  shap_interaction <- shap_values %>% 
    filter(.data$feature_1 == main_feature & .data$feature_2 == interacting_feature)
  
  ggplot(shap_interaction, aes(x = feature_1_value, y = feature_2_value, colour = shap_value, shape = label)) +
    geom_jitter(width = jitter_width, height = jitter_height) +
    scale_colour_gradient2(low = "#e41a1c", mid = "#edf8b1", high = "#377eb8") +
    scale_shape_manual(values = c(19, 1)) +
    labs(x = main_feature, y = interacting_feature, colour = "SHAP value") +
    PLOT_THEME
}
