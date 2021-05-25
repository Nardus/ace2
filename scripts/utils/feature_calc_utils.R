## Functions to calculate features specific to each training set
# - These features rely on a summary of known hosts, and hence must be calculated
#   from training species only to avoid a data leak

require(dplyr)
require(tidyr)


#' Find the most common character.
#' If there is a tie, samples randomly.
#' 
most_common_value <- function(observed) {
  observed <- na.omit(observed)
  
  if (!is.factor(observed)) 
    observed <- factor(observed)
  
  counts <- tabulate(observed)
  res <- levels(observed)[counts == max(counts)]
  
  if (length(res) > 1)
    return(sample(res, size = 1))
  
  res
}


# ---- Distance to closest infectable species -----------------------------------------------------
#' Get the closest infectable species (excluding self) for each evidence level, e.g. how far away
#' is the current ACE2 sequence from that of a host with level-1 evidence that it can be infected?
#' 
#' @param pairwise_dist_data: data frame containing pairwise distances in long format, filtered
#'                            to include only distances to training-set hosts
#' @param metadata: meta-data listing accessions along with a "label" column specifying
#'                  the response being predicted (with character values "True" and 
#'                  "False" indicating positives and negatives, respectively)
#' 
get_dist_to_closest_positive <- function(pairwise_dist_data, metadata) {
  dist_data <- pairwise_dist_data %>% 
    left_join(metadata, by = c("other_seq" = "ace2_accession")) %>% 
    
    filter(.data$label == "True") %>%                                     # Only positive neighbours
    filter(.data$ace2_accession != .data$other_seq)                       # Exclude self
  
  overall_dist <- dist_data %>% 
    group_by(.data$ace2_accession) %>% 
    summarise(closest_positive_overall = min(.data$distance), 
              .groups = "drop")
  
  level_dists <- dist_data %>% 
    group_by(.data$ace2_accession, .data$evidence_level) %>% 
    summarise(closest_positive = min(.data$distance), 
              .groups = "drop") %>% 
    pivot_wider(id_cols = .data$ace2_accession, 
                names_from = "evidence_level", 
                names_prefix = "closest_positive_l",
                values_from = "closest_positive")
  
  full_join(overall_dist, level_dists, by = "ace2_accession")
}


# ---- Variable sites (distance to consensus) -------------------------------------------
#' Get most common amino acid found among positives (for a given binary reponse being 
#' predicted, e.g. susceptible hosts) at each site, then calculate distance to this 
#' amino acid.
#' 
#' @param variable_sites_training: data frame of variable sites observed in training data
#' @param variable_sites_all: data frame of variable sites for all hosts
#' @param metadata: meta-data listing accessions, along with a "label" column specifying
#'                  the response being predicted (with character values "True" and 
#'                  "False" indicating positives and negatives, respectively)
#' 
get_consensus_dist <- function(variable_sites_training, variable_sites_all, metadata) {
  consensus_vals <- variable_sites_training %>% 
    left_join(metadata, by = "ace2_accession") %>% 
    filter(.data$label == "True") %>% 
    ungroup() %>% 
    summarise(across(starts_with("variable_site"), most_common_value)) %>% 
    pivot_longer(everything(), names_to = "site", values_to = "most_common")
  
  variable_sites_all %>% 
    pivot_longer(starts_with("variable_site"), names_to = "site", values_to = "observed") %>% 
    left_join(consensus_vals, by = "site") %>% 
    group_by(.data$ace2_accession, .data$site) %>% 
    summarise(distance = get_site_dist_vectorised(.data$most_common, .data$observed, 
                                                  type = "grantham", ignore_na = TRUE),
              .groups = "drop") %>% 
    pivot_wider(id_cols = .data$ace2_accession, names_from = "site", values_from = "distance", 
                names_prefix = "dist_")
}
