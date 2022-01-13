## Build a taxonomy lookup-table for all species with ACE2 data

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(taxize)
})

# Merge species from all datasets
ncbi_metadata <- read.csv("data/internal/NCBI_ACE2_orthologs.csv") %>% 
  pull(.data$Scientific.name)

internal_metadata <- read.csv("data/internal/ace2_accessions.csv") %>% 
  filter(!is.na(.data$ace2_accession)) %>% 
  pull(.data$species)

existing_predictions <- read.csv("data/internal/existing_predictions.csv") %>% 
  select(.data$species, .data$citation_key)

haddock_scores <- readRDS("data/calculated/all_haddock_scores.rds") %>% 
  pull(.data$species)


all_spp <- c(ncbi_metadata, internal_metadata, existing_predictions$species, haddock_scores) %>% 
  unique()

# Get taxonomy
taxonomy_table <- classification(all_spp, db = "ncbi")

# Check missing
missing_spp <- names(taxonomy_table)[is.na(taxonomy_table)]

allowed_missing <- existing_predictions %>% 
  group_by(.data$species) %>% 
  mutate(non_ace2_only = n() == 1 & .data$citation_key == "fischhoff2021a") %>% # species unique to Fischhoff et al. (these might not be on GenBank)
  filter(.data$non_ace2_only) %>% 
  pull(.data$species)

allowed_missing <- c(allowed_missing, "Bindicus x btaurus") # From haddock scores, not needed

stopifnot(all(missing_spp %in% allowed_missing))

# Convert to data frame
final_taxonomy <- taxonomy_table %>% 
  rbind() %>% 
  filter(! .data$rank %in% c("no rank", "clade")) %>% 
  select(-.data$id) %>% 
  pivot_wider(id_cols = "query", names_from = "rank", values_from = "name") %>% 
  rename(internal_name = .data$query) %>% 
  mutate(species = if_else(is.na(.data$species), .data$internal_name, .data$species))

# Ensure at least the species column is available for all records
missing_spp <- data.frame(internal_name = missing_spp,
                          species = missing_spp)

final_taxonomy <- bind_rows(final_taxonomy, missing_spp)

# Fix known issues / synonyms:
# - Changes required to match timetree.org phylogeny and/or IUCN range maps
# - Based on homotypic synonyms listed in NCBI taxonomy
final_taxonomy <- final_taxonomy %>% 
  mutate(species = case_when(internal_name == "Canis lupus familiaris" ~ "Canis familiaris", # Avoid duplicates from wolves
                             species == "Tupaia chinensis" ~ "Tupaia belangeri",
                             species == "Nannospalax galili" ~ "Nannospalax ehrenbergi",
                             species == "Grammomys surdaster" ~ "Grammomys dolichurus",
                             TRUE ~ species),
         order = case_when(species == "Echinops telfairi" ~ "Afrosoricida",
                           species == "Chrysochloris asiatica" ~ "Afrosoricida",
                           TRUE ~ order))

# Output
saveRDS(final_taxonomy, "data/calculated/taxonomy.rds")
