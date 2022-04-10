## Calculate features from ACE2 protein sequences
## - Only features constant across all training sets calculated here, for the rest, see 
##   'feature_calc_utils.R'

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(seqinr)
})

source("scripts/utils/aa_distance_utils.R")
source("scripts/utils/aa_property_utils.R")

N_CORES <- 16  # Number of parallel cores/threads allowed


# ---- Load sequences + metadata ------------------------------------------------------------------
ace2_alignment <- read.fasta("data/calculated/ace2_protein_alignment.fasta", seqtype = "AA")

aa_properties <- read_aaindex_multi(index_names = c("hydrophobicity", "polarity", 
                                                    "charge", "volume"),
                                    base_path = "data/internal/aa_index_properties")

# Data for known species
infection_data <- read_rds("data/calculated/cleaned_infection_data.rds")
shedding_data <- read_rds("data/calculated/cleaned_shedding_data.rds")

sequence_metadata <- infection_data %>% 
  full_join(shedding_data, by = c("species", "ace2_accession")) %>% 
  select(.data$species, .data$ace2_accession)

# Data for holdout species
# - Removing both matched species names AND accessions already used to represent related species
taxonomy <- read_rds("data/calculated/taxonomy.rds")

additional_metadata <- read_csv("data/internal/NCBI_ACE2_orthologs.csv",
                                col_types = cols(.default = "c")) %>% 
  select(internal_name = .data$`Scientific name`, 
         ace2_accession = .data$`RefSeq Protein accessions`) %>% 
  left_join(taxonomy, by = "internal_name") %>% 
  select(-.data$internal_name) %>% 
  filter(!.data$species %in% sequence_metadata$species & 
           !.data$ace2_accession %in% sequence_metadata$ace2_accession) %>% 
  distinct(.data$species, .data$ace2_accession) %>% 
  filter(!.data$ace2_accession %in% c("XP_006194263.1"))  # Remove duplicate for wild camels

sequence_metadata <- sequence_metadata %>% 
  bind_rows(additional_metadata)

stopifnot(all(sequence_metadata$ace2_accession %in% names(ace2_alignment)))


# Ensure alignment matches this selection
ace2_alignment <- ace2_alignment[unique(sequence_metadata$ace2_accession)]


# ---- Pairwise distances -------------------------------------------------------------------------
grantham_dists <- get_all_distances(ace2_alignment, type = "grantham", cores = N_CORES)

dist_data <- tibble(ace2_accession = rownames(grantham_dists),
                    data.frame(grantham_dists)) %>% 
  pivot_longer(-.data$ace2_accession, names_to = "other_seq", values_to = "distance")



# ---- Distance to humans -------------------------------------------------------------------------
dist_to_humans <- dist_data %>% 
  left_join(sequence_metadata, by = c("other_seq" = "ace2_accession")) %>% 
  filter(.data$species == "Homo sapiens") %>% 
  select(.data$ace2_accession,
         distance_to_humans = .data$distance)

stopifnot(nrow(dist_to_humans) == length(ace2_alignment))


# ---- Variable sites (categorical) ---------------------------------------------------------------
# To find variable sites, compare all characters except X
# Gaps are considered another character class IF is is internal (i.e. a deletion). At the edges of
# a sequence, "-" is actually missing data
replace_external_gap_chars <- function(x) {
  seq_starts <- min(which(x != "-"))
  seq_ends <- max(which(x != "-"))
  
  new_seq <- rep(NA_character_, length(x))
  new_seq[seq_starts:seq_ends] <- x[seq_starts:seq_ends]
  
  new_seq
}

alignment_mat <- ace2_alignment %>% 
  bind_cols() %>%                   # Each sequence is a column, rows represent alignment positions
  mutate(across(everything(), replace_external_gap_chars),
         across(everything(), ~ if_else(.x == "X", NA_character_, .x))) %>% 
  
  t()                               # Rotate: each row now a sequence

variant_counts <- apply(alignment_mat, MARGIN = 2, FUN = n_distinct, na.rm = TRUE)
colnames(alignment_mat) <- as.character(seq(1, ncol(alignment_mat)))

variable_sites <- alignment_mat[, variant_counts > 1] %>% 
  as_tibble(rownames = "ace2_accession", .name_repair = ~ paste0("variable_site_", .x))

stopifnot(nrow(variable_sites) == length(ace2_alignment))


# ---- Variable sites (physico-chemical properties) -----------------------------------------------
site_properties <- variable_sites %>% 
  pivot_longer(!ace2_accession, names_to = "site", names_prefix = "variable_site_", 
               values_to = "AA") %>% 
  left_join(aa_properties, by = "AA") %>% 
  select(-.data$AA) %>% 
  pivot_wider(names_from = .data$site, values_from = c(hydrophobicity, polarity, charge, volume),
              names_glue = "property_{.value}_{site}") # Prefix all columns with "property", e.g. property_polarity_2

stopifnot(nrow(site_properties) == length(ace2_alignment))


# ---- Output -------------------------------------------------------------------------------------
write_rds(dist_data, "data/calculated/features_pairwise_dists.rds")
write_rds(dist_to_humans, "data/calculated/features_dist_to_humans.rds")
write_rds(variable_sites, "data/calculated/features_variable_sites.rds")
write_rds(site_properties, "data/calculated/features_site_properties.rds")
