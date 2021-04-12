## Collect Haddock scores from Fischoff et al. 2021, and match species to our data
## - As in that paper, scores from the 10 best models are averaged

library(dplyr)
library(stringr)
library(readr)
library(readxl)

DATA_PATH <- "data/external/haddock_scores/ace2-orthologs-dataset/refined_models/runs"

ace2_accessions <- read_excel("data/internal/ace2_accessions.xlsx")


# ---- Load scores --------------------------------------------------------------------------------
get_average_score <- function(folder_name) {
  full_path <- file.path(DATA_PATH, folder_name, "ranks.models")
  message(full_path)
  
  species <- folder_name %>% 
    str_remove("^autosubmit-") %>% 
    str_remove("_params$") %>% 
    str_replace_all("_", " ") %>% 
    str_to_sentence()
  
  scores <- readLines(full_path) %>% 
    str_extract("\\{ [-.[:digit:]]+ \\}") %>% 
    str_extract("[-.[:digit:]]+") %>% 
    as.numeric()
  
  stopifnot(length(scores) == 10)
  
  tibble(species = species, haddock_score = mean(scores))
}


species_folders <- list.dirs(DATA_PATH, full.names = FALSE, recursive = FALSE)

haddock_scores <- lapply(species_folders, get_average_score) %>% 
  bind_rows()


# ---- Match species ------------------------------------------------------------------------------
# Not all species have known ACE2 sequences
stopifnot(length(unique(ace2_accessions$species)) == nrow(ace2_accessions))

replacements <- ace2_accessions %>% 
  mutate(haddock_species = if_else(.data$haddock_species == "this_species",
                                   .data$species, .data$haddock_species)) %>% 
  select(.data$species, .data$haddock_species)

# Check for missing species:
if (any(!replacements$haddock_species %in% haddock_scores$species)) {
  missing_spp <- replacements$species[!replacements$haddock_species %in% haddock_scores$species] %>% 
    sort() %>% 
    paste(collapse = ", ")
  warning("Some species do not have HADDOCK scores: ", missing_spp)
}

data_scores <- replacements %>% 
  left_join(haddock_scores, by = c("haddock_species" = "species")) %>% 
  select(-.data$haddock_species)


# ---- Output -------------------------------------------------------------------------------------
write_csv(haddock_scores, "data/calculated/all_haddock_scores.csv")
write_rds(data_scores, "data/calculated/features_haddock_scores.rds")