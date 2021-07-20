## Utility functions for plotting

library(seqinr)
library(readr)
library(stringr)


# ---- Convert feature coordinates -----------------------------------------------------------------
# Find the offset needed to make coordinates match human ace2 (count gaps added in alignment)
alignment <- read.fasta("data/calculated/ace2_protein_alignment.fasta")
seq_metadata <- read_csv("data/internal/ace2_accessions.csv", col_types = cols(.default = "c"))

human_acc <- seq_metadata$ace2_accession[seq_metadata$species == "Homo sapiens"]
human_seq <- alignment[[which(names(alignment) == human_acc)]]

as_human_coord <- function(alignment_coord, human = human_seq) {
  if (is.na(alignment_coord))
    return(NA_integer_)
  
  gaps <- human == "-"
  offset_by <- sum(gaps[1:alignment_coord])
  
  alignment_coord - offset_by
}

# vectorised version:
as_human_coord_v <- function(alignment_coords) {
  sapply(alignment_coords, as_human_coord, simplify = TRUE)
}


# ---- Human-readable feature names ----------------------------------------------------------------
# Expects a data.frame with a column named "feature"
add_readable_feature_names <- function(x) {
  x %>% 
    mutate(feature_type = case_when(startsWith(.data$feature, "dist_") ~ "Consensus distance",
                                    startsWith(.data$feature, "variable_site") ~ "Amino acid identity",
                                    .data$feature == "haddock_score" ~ "HADDOCK binding score",
                                    TRUE ~ "Other"),
           feature_position = if_else(.data$feature_type %in% c("Consensus distance", "Amino acid identity"),
                                      str_extract(.data$feature, "[[:digit:]]+$"), 
                                      NA_character_),
           feature_position_corrected = as_human_coord_v(as.integer(.data$feature_position)),
           feature_label = case_when(.data$feature_type == "Other" ~ .data$feature,
                                     .data$feature == "haddock_score" ~ .data$feature_type, 
                                     TRUE ~ sprintf("%s (%d)", .data$feature_type, .data$feature_position_corrected)))
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



# ---- Specific plot types -------------------------------------------------------------------------

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
