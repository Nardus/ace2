# Cluster ACE2 positions by correlation (across all available sequences)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(parallel)
  library(recipes)
  library(apcluster)
  library(ggplot2)
  library(cowplot)
  
  source("scripts/utils/aa_distance_utils.R")
  source("scripts/utils/feature_calc_utils.R")
  source("scripts/utils/plot_utils.R")
  source("scripts/plotting/plotting_constants.R")
})

N_CORES <- 20  # Number of parallel threads used

# ---- Get clusters --------------------------------------------------------------------------------
# Load data
variable_sites <- readRDS("data/calculated/features_variable_sites.rds")

# Calculate amino acid distances at all sites
# - A similar operation is used during training, but here we are:
#     - using ALL data, not just the training data
#     - using calculating distance from the population consensus, not the consensus among
#       positives only
metadata <- data.frame(ace2_accession = variable_sites$ace2_accession,
                       label = "True", 
                       stringsAsFactors = FALSE)

site_distances <- get_consensus_dist(variable_sites, variable_sites, metadata) %>% 
  select(-.data$ace2_accession)

# Remove uninformative sites
# - These will cause errors in cor() below
# - These sites won't make it into the final model anyway
site_distances <- recipe(site_distances) %>% 
  step_zv(everything()) %>% 
  prep() %>% 
  bake(NULL)

# Get pairwise correlations while ignoring missing values:
get_cor <- function(var1, var2, data = site_distances) {
  res <- cor(data[[var1]], data[[var2]], method = "spearman", use = "pairwise.complete.obs")
  
  data.frame(var1 = var1, var2 = var2, spearman_cor = res)
}

all_combos <- expand_grid(var1 = colnames(site_distances),
                          var2 = colnames(site_distances))

correlations <- mcmapply(get_cor,
                         var1 = all_combos[["var1"]],
                         var2 = all_combos[["var2"]],
                         SIMPLIFY = FALSE,
                         mc.cores = N_CORES) %>% 
  bind_rows() %>% 
  pivot_wider(names_from = .data$var2, values_from = .data$spearman_cor)

# Cluster
cor_mat <- correlations %>% 
  select(-.data$var1) %>% 
  as.matrix()

rownames(cor_mat) <- correlations$var1

clusters <- apcluster(s = cor_mat, details = TRUE)

# Extract clustered site names
clustered_locations <- sapply(1:length(clusters@clusters), 
                              function(i) data.frame(cluster = i,
                                                     exemplar = names(clusters@exemplars)[i],
                                                     feature = names(clusters@clusters[[i]])),
                              simplify = FALSE) %>% 
  bind_rows() %>% 
  add_readable_feature_names()  # add corrected feature positions


# ---- Output --------------------------------------------------------------------------------------
dir.create("output/plots/intermediates/", recursive = TRUE)

# Save cluster info
clustered_locations %>% 
  select(.data$cluster, .data$exemplar, .data$feature_position, .data$feature_position_corrected) %>% 
  write_rds("output/plots/intermediates/feature_clusters.rds")
