# Check correlation among sites

library(dplyr)
library(tidyr)
library(stringr)
library(apcluster)
library(caret)

source("scripts/utils/plot_utils.R")


full_model <- readRDS("output/all_data/infection/all_features/training_results.rds")

# - Using dist_ variables for this, as the only continuous representation of sites
continuous_data <- full_model$trainingData %>% 
  select(-species, -ace2_accession, -evidence_level, -label) %>% 
  select(starts_with("dist_"))

uniniformative <- nearZeroVar(continuous_data)
continuous_data <- continuous_data[, -uniniformative]


# Get distances while ignoring missing values:
get_cor <- function(var1, var2) {
  res <- cor(continuous_data[[var1]], continuous_data[[var2]], method = "spearman", use = "pairwise.complete.obs")
  
  data.frame(var1 = var1, var2 = var2, spearman_cor = res)
}

all_combos <- expand_grid(var1 = colnames(continuous_data),
                          var2 = colnames(continuous_data))

correlations <- mapply(get_cor,
                       var1 = all_combos[["var1"]],
                       var2 = all_combos[["var2"]],
                       SIMPLIFY = FALSE) %>% 
  bind_rows() %>% 
  pivot_wider(names_from = .data$var2, values_from = .data$spearman_cor)

cor_mat <- correlations %>% 
  select(-.data$var1) %>% 
  as.matrix()

rownames(cor_mat) <- correlations$var1


clusters <- apcluster(s = cor_mat, details = TRUE)


clustered_locations <- sapply(1:length(clusters@clusters), 
                              function(i) data.frame(cluster = i,
                                                     exemplar = names(clusters@exemplars)[i],
                                                     feature = names(clusters@clusters[[i]])),
                              simplify = FALSE) %>% 
  bind_rows() %>% 
  add_readable_feature_names()  # add corrected feature positions



# Plot:
# Binding sites taken from human ACE2 genbank entry (NP_001358344.1)
s_binding_sites <- tribble(
  ~start_pos, ~stop_pos, ~name,
  30,         41,        "ECO:0000269|PubMed:15791205 1",
  82,         84,        "ECO:0000269|PubMed:15791205 2",
  353,        357,       "ECO:0000269|PubMed:15791205 3"
)

all_s_binding_inds <- mapply(seq,
                             from = s_binding_sites$start_pos, 
                             to = s_binding_sites$stop_pos) %>% 
  unlist()


cluster_binding <- clustered_locations %>% 
  group_by(.data$cluster) %>% 
  summarise(includes_binding_site = any(.data$feature_position_corrected %in% all_s_binding_inds))

clustered_locations <- clustered_locations %>% 
  left_join(cluster_binding, by = "cluster")


ggplot(clustered_locations) +
  geom_rect(aes(xmin = start_pos, xmax = stop_pos, ymin = -Inf, ymax = Inf), 
            fill = "grey80", data = s_binding_sites) +
  geom_tile(aes(x = feature_position_corrected, y = factor(cluster), fill = includes_binding_site)) +
  labs(x = NULL, y = "Number of replicates") +
  PLOT_THEME +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 5.5, r = 5.5, b = 0, l = 5.5))



# Save cluster info
clustered_locations %>% 
  select(.data$cluster, .data$exemplar, .data$feature_position_corrected) %>% 
  write_rds("output/plots/feature_clusters_infection.rds")
