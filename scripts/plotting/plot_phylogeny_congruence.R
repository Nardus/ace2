# Plot and measure congruence between the ACE2 phylogeny and timetree

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ape)
  library(phytools)
  library(phylogram)
  library(dendextend)
  library(parallel)
  library(ggplot2)
  library(grid)
})

source("scripts/plotting/plotting_constants.R")
source("scripts/utils/timetree_constants.R")

N_CORES <- 20

# Phylogenies
time_tree <- read.tree("data/internal/timetree_amniota.nwk")
ace2_tree <- read.tree("data/calculated/gene_tree/ace_genetree.treefile")

as.dendrogram(ace2_tree)

# Metadata:
ncbi_metadata <- read.csv("data/internal/NCBI_ACE2_orthologs.csv") %>% 
  select(species = .data$Scientific.name, ace2_accession = .data$RefSeq.Protein.accessions)

internal_metadata <- read.csv("data/internal/ace2_accessions.csv") %>% 
  select(.data$species, .data$ace2_accession) %>% 
  filter(!is.na(.data$ace2_accession))

ace2_metadata <- bind_rows(ncbi_metadata, internal_metadata) %>% 
  distinct()

# Infection data
infection_data <- readRDS("data/calculated/cleaned_infection_data.rds")


# ---- Clean tip labels ----------------------------------------------------------------------------
# Timetree
# - Using names as in NCBI taxonomy
time_tree$tip.label <- str_replace_all(time_tree$tip.label, "_", " ")

stopifnot(all(names(TIMETREE_TAXONOMY_CORRECTIONS) %in% time_tree$tip.label))
time_tree$tip.label <- if_else(time_tree$tip.label %in% names(TIMETREE_TAXONOMY_CORRECTIONS), 
                              as.character(TIMETREE_TAXONOMY_CORRECTIONS[time_tree$tip.label]), 
                              time_tree$tip.label)

# ACE2
ace2_tree <- drop.tip(ace2_tree, "XP_027389727.1")  # Hybrid (Bos indicus x Bos taurus) won't be matched in timetree

spp_map <- ace2_metadata$species %>% 
  if_else(. == "Canis lupus familiaris", "Canis familiaris", .) %>% 
  str_extract("^[[:alpha:]]+ [[:alpha:]]+") # Drop subspecies info

names(spp_map) <- ace2_metadata$ace2_accession

ace2_tree$tip.label <- as.character(spp_map[ace2_tree$tip.label])

# Reduce timetree to species with ACE2 sequences
stopifnot(all(ace2_tree$tip.label %in% c(time_tree$tip.label, TIMETREE_KNOWN_MISSING)))

time_tree <- keep.tip(time_tree, ace2_tree$tip.label[!ace2_tree$tip.label %in% TIMETREE_KNOWN_MISSING])
ace2_tree <- drop.tip(ace2_tree, TIMETREE_KNOWN_MISSING)


# ---- Calculate correlation -----------------------------------------------------------------------
get_baker_null <- function(i, dend1, dend2) {
  dend1_random <- sample.dendrogram(dend1, replace = FALSE)
  dend2_random <- sample.dendrogram(dend2, replace = FALSE)
  cor_bakers_gamma(dend1_random, dend2_random)
}

# Convert to dendrograms
ace2_tree <- midpoint.root(ace2_tree)

ace2_dendro <- as.dendrogram(ace2_tree)
time_dendro <- as.dendrogram(time_tree)

# Calculate correlation
obs_gamma <- cor_bakers_gamma(ace2_dendro, time_dendro)

# Permutation test:
n_sims <- 1000
gamma_null <- mclapply(1:n_sims, 
                       FUN = get_baker_null, 
                       dend1 = ace2_dendro,
                       dend2 = time_dendro,
                       mc.cores = N_CORES)

gamma_null <- unlist(gamma_null)
p_val <- sum(gamma_null >= obs_gamma)/n_sims
p_val <- if_else(p_val == 0, 
                 sprintf("< %3.3f", 1/n_sims), # Smallest p-value possible with current n_sims
                 sprintf("%3.3f", p_val))

cat(sprintf("\n\nBaker's gamma: %3.3f, exact p-value from permutation test: %s\n", obs_gamma, p_val))


# ---- Correlation in pairwise phylogenetic distances ----------------------------------------------
ace2_dists <- cophenetic(ace2_tree) %>% 
  data.frame() %>% 
  rownames_to_column("from_species") %>% 
  pivot_longer(-.data$from_species, names_to = "to_species", values_to = "ace2_distance")

time_dists <- cophenetic(time_tree) %>% 
  data.frame() %>% 
  rownames_to_column("from_species") %>% 
  pivot_longer(-.data$from_species, names_to = "to_species", values_to = "time_distance")

combined_dists <- ace2_dists %>% 
  full_join(time_dists, by = c("from_species", "to_species"))

dist_plot <- ggplot(combined_dists, aes(x = time_distance, y = ace2_distance)) +
  geom_point(colour = "grey70", size = 0.8) +
  geom_smooth(method = "lm", colour = TRENDLINE_COLOUR) +
  labs(x = "Phylogenetic distance (timetree)",
       y = "Phylogenetic distance (ACE2 phylogeny)")

cat(sprintf("\nSpearman correlation in phylogenetic distances: %3.3f\n", 
            cor(combined_dists$ace2_distance, combined_dists$time_distance,
                method = "spearman")))


# ---- Tanglegram ----------------------------------------------------------------------------------
ace2_dendro <- ladder(ace2_dendro)
time_dendro <- ladder(time_dendro)


aligned_dendrograms <- intersect_trees(ace2_dendro, time_dendro) %>% 
  untangle(method = "step2side")

# Plot
infection_status <- infection_data$infected
names(infection_status) <- infection_data$species

colour_vec <- labels(aligned_dendrograms[[1]]) %>% 
  infection_status[.] %>% 
  INFECTION_STATUS_COLOURS[.] %>% 
  if_else(condition = is.na(.), MISSING_DATA_COLOUR, .)


pdf("output/plots/phylogeny_congruence.pdf", width = 7, height = 4)
slot_width <- 0.5*1/3  # tanglegram takes 3 slots (dend, connectors, dend)
layout(matrix(1:4, nrow = 1), widths = c(slot_width, slot_width, slot_width, 0.5))

tanglegram(aligned_dendrograms, 
           main_left = "ACE2 phylogeny\n(midpoint rooted)",
           main_right = "Timetree",
           color_lines = colour_vec,
           faster = TRUE, 
           lwd = 1, lab.cex = 1e-5, margin_inner = 0.4, 
           cex_main = 0.8,
           just_one = FALSE)

# ggplot alongside this:
vp <- viewport(height = unit(1, "npc"), width = unit(0.5, "npc"), 
               just=c("left","top"), 
               y=1, x=0.5)
print(dist_plot, vp = vp)

dev.off()
