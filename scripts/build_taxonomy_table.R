## Build a taxonomy lookup-table for all species with ACE2 data

suppressPackageStartupMessages({
  library(dplyr)
  library(dbplyr)
  library(tidyr)
  library(stringr)
  library(ape)
  library(taxize)
  library(taxizedb)
})

# ---- Merge species from all datasets -------------------------------------------------------------
ncbi_metadata <- read.csv("data/internal/NCBI_ACE2_orthologs.csv") %>% 
  pull(.data$Scientific.name)

internal_metadata <- read.csv("data/internal/ace2_accessions.csv") %>% 
  filter(!is.na(.data$ace2_accession)) %>% 
  pull(.data$species)

existing_predictions <- read.csv("data/internal/existing_predictions.csv") %>% 
  select(.data$species, .data$citation_key)

haddock_scores <- readRDS("data/calculated/all_haddock_scores.rds") %>% 
  pull(.data$species)

# - TimeTree species require additional cleaning to remove e.g. "Myotis cf._nattereri_IC-2011"
timetree_mammals <- read.tree("data/internal/timetree_mammalia.nwk")$tip.label %>% 
  str_replace_all("_", " ")

timetree_birds <- read.tree("data/internal/timetree_aves.nwk")$tip.label %>% 
  str_replace_all("_", " ")

timetree_species <- c(timetree_mammals, timetree_birds)

removals <- (str_detect(timetree_species, " cf.") | 
               str_detect(timetree_species, " sp.") | 
               str_detect(timetree_species, "\\w x \\w"))

timetree_species <- timetree_species[!removals]

# - Merge
all_spp <- c(ncbi_metadata, internal_metadata, existing_predictions$species, 
             haddock_scores, timetree_species) %>% 
  unique()


# ---- Manual corrections --------------------------------------------------------------------------
corrections <- c("Abrornis humei" = "Phylloscopus humei",
                 "Abrornis inornata" = "Phylloscopus inornatus",
                 "Anser canagica" = "Anser canagicus",
                 "Brachypteryx albiventris" = "Sholicola albiventris",
                 "Catherpes sumichrasti" = "Hylorchilus sumichrasti",
                 "Chamaeza mollisima" = "Chamaeza mollissima",
                 "Chrysomus cyanopus" = "Agelasticus cyanopus",
                 "Chrysomus thilius" = "Agelasticus thilius",
                 "Criniferoides leucogaster" = "Corythaixoides leucogaster",
                 "Fringillaria impetuani" = "Emberiza impetuani",
                 "Fringillaria tahapisi" = "Emberiza tahapisi",
                 "Galeopterus variegates" = "Galeopterus variegatus",
                 "Gulosus aristotelis" = "Phalacrocorax aristotelis",
                 "Hylobates alibarbis" = "Hylobates albibarbis",
                 "Hystrix brachyurus" = "Hystrix brachyura",
                 "Leimacomys buttneri" = "Leimacomys buettneri",
                 "Luscinia hyperythra" = "Tarsiger hyperythrus",
                 "Luscinia indica" = "Tarsiger indicus",
                 "Luscinia johnstoniae" = "Tarsiger johnstoniae",
                 "Minopterus macrocneme" = "Miniopterus macrocneme",
                 "Murexchinus melanurus" = "Murexia melanurus",
                 "Nannopterum auritus" = "Phalacrocorax auritus",
                 "Nannopterum brasilianus" = "Phalacrocorax brasilianus",
                 "Nesolagus timinsi" = "Nesolagus timminsi",
                 "Phoenicoparrus minor" = "Phoeniconaias minor",
                 "Poecilotriccus albifascies" = "Poecilotriccus albifacies",
                 "Proechimys trinitatus" = "Proechimys trinitatis",
                 "Schoeniclus aureolus" = "Emberiza aureola",
                 "Urile penicillatus" = "Phalacrocorax penicillatus",
                 "Urile urile" = "Phalacrocorax urile")

all_species_corrected <- if_else(all_spp %in% names(corrections),
                                 corrections[all_spp],
                                 all_spp)


# ---- Get taxonomy --------------------------------------------------------------------------------
db_download_itis()
db_download_ncbi()

# Attempt 1 - direct matches on ITIS
itis_ids <- name2taxid(all_species_corrected, out_type = "summary", db = "itis")

itis_taxonomy <- classification(unique(itis_ids$id), db = "itis")
itis_invalid <- sapply(itis_taxonomy, nrow) == 0
itis_taxonomy <- itis_taxonomy[!itis_invalid]

# Attempt 2 - get valid names for species known to ITIS (id present, but no taxonomy)
valid_names <- itis_invalid[itis_invalid] %>% 
  names() %>% 
  itis_acceptname()

itis_taxonomy_2 <- classification(unique(valid_names$acceptedtsn), db = "itis")

# Attempt 3 - try finding non-ITIS species in NCBI taxonomy
missing_spp <- all_species_corrected[!all_species_corrected %in% itis_ids$name]

ncbi_ids <- name2taxid(missing_spp, out_type = "summary", db = "ncbi")
ncbi_taxonomy <- classification(unique(ncbi_ids$id), db = "ncbi")


# ---- Convert to data frames ----------------------------------------------------------------------
# Direct ITIS
itis_ids <- itis_ids %>% 
  mutate(id = as.character(.data$id))

itis_taxonomy <- itis_taxonomy %>% 
  bind_rows(.id = "query_id") %>% 
  select(-.data$id) %>% 
  pivot_wider(id_cols = "query_id", names_from = "rank", values_from = "name") %>% 
  left_join(itis_ids, by = c("query_id" = "id")) %>% 
  rename(internal_name = .data$name)

# Corrected ITIS
uncorrected_names <- valid_names %>% 
  left_join(itis_ids, by = c("submittedtsn" = "id")) %>% 
  select(id = .data$acceptedtsn, .data$name)  # Old name mapped to corrected ID

class(itis_taxonomy_2) <- "list"

itis_taxonomy_2 <- itis_taxonomy_2 %>% 
  bind_rows(.id = "query_id") %>% 
  select(-.data$id) %>% 
  pivot_wider(id_cols = "query_id", names_from = "rank", values_from = "name") %>% 
  left_join(uncorrected_names, by = c("query_id" = "id")) %>% 
  rename(internal_name = .data$name)

# NCBI
class(ncbi_taxonomy) <- "list"

ncbi_taxonomy <- ncbi_taxonomy %>% 
  bind_rows(.id = "query_id") %>% 
  select(-.data$id) %>% 
  filter(! .data$rank %in% c("no rank", "clade")) %>% 
  pivot_wider(id_cols = "query_id", names_from = "rank", values_from = "name") %>% 
  left_join(ncbi_ids, by = c("query_id" = "id")) %>% 
  rename(internal_name = .data$name)


# Merge
taxonomy_table <- bind_rows(itis_taxonomy, itis_taxonomy_2, ncbi_taxonomy)

# Reverse name corrections from above
#  - Need to keep both versions since e.g. typos may be unique to one dataset
reverse_corrections <- names(corrections)
names(reverse_corrections) <- unname(corrections)

extra_taxonomy <- taxonomy_table %>% 
  filter(.data$internal_name %in% names(reverse_corrections)) %>% 
  mutate(internal_name = reverse_corrections[.data$internal_name])

taxonomy_table <- bind_rows(taxonomy_table, extra_taxonomy)


# ---- Check missing -------------------------------------------------------------------------------
missing_spp <- all_spp[!all_spp %in% taxonomy_table$internal_name]

warning("Taxonomy could not be retrieved for some species: ", 
        paste(missing_spp, collapse = ", "))

# Ensure at least the species column is available for all records
missing_spp <- data.frame(internal_name = missing_spp,
                          species = missing_spp)

final_taxonomy <- bind_rows(taxonomy_table, missing_spp)


# ---- Fix known issues / synonyms -----------------------------------------------------------------
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


# ---- Output --------------------------------------------------------------------------------------
saveRDS(final_taxonomy, "data/calculated/taxonomy.rds")
