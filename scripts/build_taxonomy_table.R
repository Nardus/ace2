## Build a taxonomy lookup-table for all species with ACE2 data

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(taxize)
})

ncbi_metadata <- read.csv("data/internal/NCBI_ACE2_orthologs.csv") %>% 
  select(species = .data$Scientific.name, ace2_accession = .data$RefSeq.Protein.accessions)

internal_metadata <- read.csv("data/internal/ace2_accessions.csv") %>% 
  select(.data$species, .data$ace2_accession) %>% 
  filter(!is.na(.data$ace2_accession))

ace2_metadata <- bind_rows(ncbi_metadata, internal_metadata) %>% 
  distinct()


all_spp <- unique(ace2_metadata$species)
taxonomy_table <- classification(all_spp, db = "ncbi")

final_taxonomy <- taxonomy_table %>% 
  rbind() %>% 
  filter(! .data$rank %in% c("no rank", "clade")) %>% 
  select(-.data$id) %>% 
  pivot_wider(id_cols = "query", names_from = "rank", values_from = "name")

saveRDS(final_taxonomy, "data/calculated/taxonomy.rds")
