# Generate a supplemental table of predictions and report accuracy on holdout species 
# (i.e. those reported after data collection was finalized)

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(readr)
  library(writexl)
})

# ---- Data ----------------------------------------------------------------------------------------
train_preds <- read_rds("output/all_data/infection/phylogeny/predictions.rds")
test_preds <- read_rds("output/all_data/infection/phylogeny/holdout_predictions.rds")

taxonomy <- read_rds("data/calculated/taxonomy.rds")

holdout_spp <- read_csv("data/internal/holdout_species.csv", col_types = cols(.default = "c"))


# ---- Training species ----------------------------------------------------------------------------
train_column_order <- c("class", "order", "suborder", "family", "subfamily", "species", 
                        "probability", "cutoff", "prediction", "observed")

train_preds <- train_preds %>% 
  rename(internal_name = .data$species) %>% 
  left_join(taxonomy, by = "internal_name") %>%
  select(-.data$species) %>% # Using internal names to avoid confusion
  rename(species = .data$internal_name,
         observed = .data$label,
         probability = .data$p_true) %>% 
  arrange(desc(.data$probability), .data$species) %>% 
  select(all_of(train_column_order))


# ---- Holdout species -----------------------------------------------------------------------------
test_preds <- test_preds %>% 
  filter(.data$prediction_type != "Fitted value") %>% 
  rename(internal_name = .data$species) %>% 
  left_join(taxonomy, by = "internal_name") %>%
  select(-.data$species) %>%  # Using internal names to avoid confusion
  rename(species = .data$internal_name,
         prediction = .data$predicted_label) %>% 
  arrange(desc(.data$probability), .data$species)

holdout_preds <- test_preds %>% 
  filter(.data$species %in% holdout_spp$species) %>% 
  mutate(observed = "True") %>% 
  select(all_of(train_column_order))


# ---- All other species ---------------------------------------------------------------------------
test_column_order <- head(train_column_order, -1)  # No observed label here

other_preds <- test_preds %>% 
  filter(!.data$species %in% holdout_spp$species) %>% 
  select(all_of(test_column_order))

# Remove invalid species
# (TimeTree contains a few spurious entries like "Myotis sp._7_W_Ecuador")
other_preds <- other_preds %>% 
  filter(!is.na(.data$class))

# Split by mammal/bird for clarity
mammal_preds <- other_preds %>% 
  filter(.data$class == "Mammalia")

bird_preds <- other_preds %>% 
  filter(.data$class == "Aves")


# ---- Output --------------------------------------------------------------------------------------
predictions <- list("(a) Training species" = train_preds,
                    "(b) Recently infected" = holdout_preds,
                    "(c) Other (mammals)" = mammal_preds,
                    "(d) Other (birds)" = bird_preds)

write_xlsx(predictions, "output/si_tables/supplement_predictions.xlsx")


# ---- Values mentioned in text --------------------------------------------------------------------
stopifnot(length(unique(holdout_preds$species)) == nrow(holdout_preds))
stopifnot(nrow(holdout_preds) == length(unique(holdout_spp$species)))

cat(
  sprintf("\nSusceptible species recognized after data collection was complete: %i of %i correctly predicted",
          sum(holdout_preds$prediction == "True"),
          nrow(holdout_preds))
)

cat("\nIncorrect predictions:", holdout_preds$species[holdout_preds$prediction != "True"])
