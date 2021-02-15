## Generate a ranked list of feature importance for a trained model 

# ---- Input args ---------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(xgboost)
  library(caret)
})

parser <- ArgumentParser(description = "Train ACE2 models using varying subsets of data and features")

parser$add_argument("model_path", type = "character",
                    help = paste("Path to the folder containing the model to use. Output will be",
                                 "named 'feature_usage.rds' and placed in this folder."))

INPUT <- parser$parse_args()


# ---- Load model ---------------------------------------------------------------------------------
full_model_path <- file.path(INPUT$model_path, "training_results.rds")
training_results <- readRDS(full_model_path)



# ---- Feature importance -------------------------------------------------------------------------
features <- training_results$trainingData %>% 
  select(-.data$species, -.data$ace2_accession, -.data$evidence_level, -.data$label)

dummy_encoding <- dummyVars(~ ., data = features)
feature_mat <- predict(dummy_encoding, features)

colnames(feature_mat) <- str_remove(colnames(feature_mat), "\\.")
feature_mat <- feature_mat[, training_results$finalModel$finalModel$feature_names]

shap_values <- predict(training_results$finalModel$finalModel, 
                       newdata = feature_mat,
                       predcontrib = TRUE)

# Summarise to feature level:
feature_importance <- shap_values %>% 
  as.data.frame() %>% 
  select(-.data$BIAS) %>% 
  pivot_longer(everything(), names_to = "feature", values_to = "shap_value") %>% 
  group_by(.data$feature) %>% 
  summarise(importance = mean(abs(.data$shap_value)),
            .groups = "drop")

# Handle dummy coded features:
# - Defining importance as the total effect across all columns making up the feature
feature_importance <- feature_importance %>% 
  mutate(feature = if_else(startsWith(.data$feature, "variable_site_"), 
                                      substring(.data$feature, 1, nchar(.data$feature) - 1),
                                      .data$feature)) %>% 
  group_by(.data$feature) %>% 
  summarise(importance = sum(.data$importance),
            .groups = "drop")



# ---- Output ---------------------------------------------------------------------------------
out_path <- file.path(INPUT$model_path, "feature_usage.rds")
saveRDS(feature_importance, out_path)
