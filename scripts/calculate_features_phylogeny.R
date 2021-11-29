## Summarise the timetree phylogeny using phylogenetic eigenvectors

# The MPSEM package cannot be installed by conda, so install it manually if needed:
if (! "MPSEM" %in% installed.packages()[, 1]) {
    install.packages("MPSEM", repo="http://cran.us.r-project.org")
}

suppressPackageStartupMessages({
    library(ape)
    library(MPSEM)
    library(dplyr)
    library(tibble)
    library(stringr)
})

# Data
source("scripts/utils/timetree_constants.R")

timetree_all <- read.tree("data/internal/timetree_amniota.nwk")
timetree_mammals <- read.tree("data/internal/timetree_mammalia.nwk")
timetree_birds <- read.tree("data/internal/timetree_aves.nwk")

known_data <- readRDS("data/calculated/cleaned_infection_data.rds")


# ---- Reduce tree to mammals and birds ------------------------------------------------------------
keep <- c(timetree_mammals$tip.label, timetree_birds$tip.label)
timetree <- keep.tip(timetree_all, keep)


# ---- Clean up species names ----------------------------------------------------------------------
timetree$tip.label <- str_replace(timetree$tip.label, "_", " ")

# Correct taxonomy:
stopifnot(all(names(TIMETREE_TAXONOMY_CORRECTIONS) %in% timetree$tip.label))

timetree$tip.label <- with(timetree,
                           if_else(tip.label %in% names(TIMETREE_TAXONOMY_CORRECTIONS),
                                   as.character(TIMETREE_TAXONOMY_CORRECTIONS[tip.label]),
                                   tip.label))

# Check that all training data will have matching features (should be the case if taxonomy 
# corrections capture all required changes):
stopifnot(all(known_data$species %in% timetree$tip.label))


# ---- Calculate phylogenetic eigenvectors ---------------------------------------------------------
treegraph <- Phylo2DirectedGraph(timetree)

eigenvectors <- PEM.build(treegraph, a=0) # This assumes Brownian motion (i.e. no selection)


# ---- Output --------------------------------------------------------------------------------------
eigen_df <- data.frame(eigenvectors) %>%
    rownames_to_column("species") %>%
    rename_with(~str_replace(., "^V", "phylogeny"), -.data$species)
    
saveRDS(eigen_df, "data/calculated/features_phylogeny_eigenvectors.rds")
