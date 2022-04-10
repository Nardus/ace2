## Report correlation between ACE2-based distance features used in models and 
##  general phylogenetic distance

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(ape)
  library(phytools)
  
  source("scripts/utils/feature_calc_utils.R")
  source("scripts/utils/timetree_constants.R")
})

# Load data
infection_data <- readRDS("data/calculated/cleaned_infection_data.rds")
timetree <- read.tree("data/internal/timetree_amniota.nwk")

# Ace2 distance data
ace2_dists <- readRDS("data/calculated/features_pairwise_dists.rds")
human_dists <- readRDS("data/calculated/features_dist_to_humans.rds")

ncbi_metadata <- read.csv("data/internal/NCBI_ACE2_orthologs.csv") %>% 
  select(species = .data$Scientific.name, ace2_accession = .data$RefSeq.Protein.accessions)

internal_metadata <- read.csv("data/internal/ace2_accessions.csv") %>% 
  select(.data$species, .data$ace2_accession) %>% 
  filter(!is.na(.data$ace2_accession))

ace2_metadata <- bind_rows(ncbi_metadata, internal_metadata) %>% 
  distinct()


# ---- Prepare phylogeny ---------------------------------------------------------------------------
# Correct names
timetree$tip.label <- str_replace(timetree$tip.label, "_", " ")

stopifnot(all(names(TIMETREE_TAXONOMY_CORRECTIONS) %in% timetree$tip.label))

timetree$tip.label <- with(timetree,
                           if_else(tip.label %in% names(TIMETREE_TAXONOMY_CORRECTIONS),
                                   as.character(TIMETREE_TAXONOMY_CORRECTIONS[tip.label]),
                                   tip.label))

stopifnot(all(infection_data$species %in% timetree$tip.label))

# Remove species with no data
timetree <- keep.tip(timetree, infection_data$species)
stopifnot(all(timetree$tip.label %in% infection_data$species))


# ---- Calculate distances -------------------------------------------------------------------------
# Phylogeny
phylo_dists <- cophenetic(timetree) %>% 
  data.frame(check.names = FALSE) %>% 
  rownames_to_column("from_species") %>% 
  pivot_longer(-.data$from_species, names_to = "to_species", values_to = "distance")

phylo_dists <- phylo_dists %>% 
  rename(ace2_accession = .data$from_species,
         other_seq = .data$to_species)

phylo_positive <- infection_data %>% 
  select(-.data$ace2_accession) %>% 
  rename(label = .data$infected, ace2_accession = .data$species) %>% 
  get_dist_to_closest_positive(pairwise_dist_data = phylo_dists, 
                               metadata = .) %>% 
  rename(species = .data$ace2_accession)


# ACE2
ace2_positive <- infection_data %>% 
  rename(label = .data$infected) %>% 
  get_dist_to_closest_positive(pairwise_dist_data = ace2_dists, 
                               metadata = .)

ace2_positive <- infection_data %>% 
  left_join(ace2_positive, by = "ace2_accession")

stopifnot(all(phylo_positive$species == ace2_positive$species))  # Code below assumes order is the same


# ---- Report correlations -------------------------------------------------------------------------
cat("\n\nSpearman correlation between phylogenetic and ACE2-based distances to the",
    "closest positive (l1, l2, l4, overall):\n",
    sprintf("%3.3f", cor(phylo_positive$closest_positive_l1, ace2_positive$closest_positive_l1)),
    sprintf("%3.3f", cor(phylo_positive$closest_positive_l2, ace2_positive$closest_positive_l2)),
    sprintf("%3.3f", cor(phylo_positive$closest_positive_l4, ace2_positive$closest_positive_l4)),
    sprintf("%3.3f", cor(phylo_positive$closest_positive_overall, ace2_positive$closest_positive_overall)),
    "\n")
