# Calculate alignment coordinates for S-binding sites

suppressPackageStartupMessages({
  library(seqinr)
  library(readr)
    
  # S-binding sites in human ACE2:
  source("scripts/utils/site_info.R")
})


# Data
alignment <- read.fasta("data/calculated/ace2_protein_alignment.fasta")
seq_metadata <- read_csv("data/internal/ace2_accessions.csv", col_types = cols(.default = "c"))

# Find the offset needed to make coordinates match human ace2 (count gaps added in alignment)
human_acc <- seq_metadata$ace2_accession[seq_metadata$species == "Homo sapiens"]
human_seq <- alignment[[which(names(alignment) == human_acc)]]

as_human_coord <- function(alignment_coord, human = human_seq) {
  if (is.na(alignment_coord))
    return(NA_integer_)
  
  if (human[alignment_coord] == "-")
    return(NA_integer_) # position does not exist in human sequence
  
  gaps <- human == "-"
  offset_by <- sum(gaps[1:alignment_coord])
  
  alignment_coord - offset_by
}

# vectorised version:
as_human_coord_v <- function(alignment_coords) {
  sapply(alignment_coords, as_human_coord, simplify = TRUE)
}

# Adjust s-binding positions
human_coords <- as_human_coord_v(1:length(human_seq))  # length(human_seq) also gives size of alignment
location_adjustments <- 1:length(human_seq)
names(location_adjustments) <- human_coords

s_binding_positions <- location_adjustments[as.character(ALL_S_BINDING_INDS)]
s_binding_positions <- unname(s_binding_positions)

# Save
write_rds(s_binding_positions, "data/calculated/s_binding_alignment_positions.rds")