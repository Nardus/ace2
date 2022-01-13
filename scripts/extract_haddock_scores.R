## Collect Haddock scores from Fischoff et al. 2021, and match species to our data
## - As in that paper, scores from the 10 best models are averaged

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
})

DATA_PATH <- "data/external/haddock_scores/ace2-orthologs-dataset/refined_models/runs"


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


# ---- Output -------------------------------------------------------------------------------------
write_rds(haddock_scores, "data/calculated/all_haddock_scores.rds")
