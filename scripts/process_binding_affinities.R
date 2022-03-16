## Process binding affinity data to match internal taxonomy

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(readxl)
})

ace2_accessions <- read_excel("data/internal/ace2_accessions.xlsx")

huang_scores <- read_csv("data/internal/existing_predictions.csv", 
                         col_types = cols(raw_score = "d", .default = "c")) %>% 
  filter(.data$citation_key == "huang2020")

haddock_scores <- read_rds("data/calculated/all_haddock_scores.rds")

taxonomy <- read_rds("data/calculated/taxonomy.rds")


# ---- Match taxonomy ------------------------------------------------------------------------------
stopifnot(all(huang_scores$species %in% taxonomy$internal_name))
stopifnot(all(haddock_scores$species %in% c(taxonomy$internal_name, "Bindicus x btaurus")))

remove_names <- c("Bos indicus", "Camelus ferus")  # Cause duplications at species level

huang_scores <- huang_scores %>% 
  filter(!.data$species %in% remove_names) %>%
  rename(internal_name = .data$species) %>% 
  left_join(taxonomy, by = "internal_name") %>% 
  select(.data$species, 
         huang_score = .data$raw_score)

haddock_scores <- haddock_scores %>% 
  filter(!.data$species %in% remove_names) %>%
  rename(internal_name = .data$species) %>% 
  left_join(taxonomy, by = "internal_name") %>% 
  select(.data$species, .data$haddock_score)

combined_scores <- huang_scores %>% 
  full_join(haddock_scores, by = "species")


# ---- Match ACE2 usage ----------------------------------------------------------------------------
# Not all species have known ACE2 sequences, so match the "closest neighbour" strategy used in metadata
stopifnot(length(unique(ace2_accessions$species)) == nrow(ace2_accessions))

replacements <- ace2_accessions %>% 
  mutate(ace2_species = if_else(.data$ace2_species == "this_species",
                                .data$species, .data$ace2_species)) %>% 
  select(.data$species, .data$ace2_species)

# Check for missing species:
if (any(!replacements$ace2_species %in% combined_scores$species)) {
  missing_spp <- replacements$species[!replacements$ace2_species %in% combined_scores$species] %>% 
    sort() %>% 
    paste(collapse = ", ")
  warning("Some species do not have binding affinity scores in either Fischhoff or Huang et al.: ", missing_spp)
}

combined_scores <- replacements %>% 
  left_join(combined_scores, by = c("ace2_species" = "species")) %>% 
  select(-.data$ace2_species)


# ---- Output ----------------------------------------------------------------------------
write_rds(combined_scores, "data/calculated/features_binding_affinity.rds")
