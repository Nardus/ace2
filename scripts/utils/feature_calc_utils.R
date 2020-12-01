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
#' @param infection_data: meta-data listing accessions and whether or not the host is susceptible
#' 
get_dist_to_closest_infectable <- function(pairwise_dist_data, infection_data) {
  pairwise_dist_data %>% 
    left_join(infection_data, by = c("other_seq" = "ace2_accession")) %>% 
    
    filter(.data$infected) %>%                                            # Only infected neighbours
    filter(.data$ace2_accession != .data$other_seq) %>%                   # Exclude self
    
    group_by(.data$ace2_accession, .data$evidence_level) %>% 
    summarise(closest_infectable = min(.data$distance), 
              .groups = "drop") %>% 
    pivot_wider(id_cols = .data$ace2_accession, 
                names_from = "evidence_level", 
                names_prefix = "closest_infectable_l",
                values_from = "closest_infectable")
}


# ---- Variable sites (distance to consensus) -------------------------------------------
#' Get most common amino acid found among susceptibles at each site, then calculate 
#' distance to this amino acid.
#' 
#' @param variable_sites_training: data frame of variable sites observed in training data
#' @param variable_sites_all: data frame of variable sites for all hosts
#' @param infection_data: meta-data listing accessions and whether or not the host is susceptible
#' 
get_consensus_dist <- function(variable_sites_training, variable_sites_all, infection_data) {
  consensus_vals <- variable_sites_training %>% 
    left_join(infection_data, by = "ace2_accession") %>% 
    filter(.data$infected) %>% 
    ungroup() %>% 
    summarise(across(starts_with("variable_site"), most_common_value)) %>% 
    pivot_longer(everything(), names_to = "site", values_to = "most_common")
  
  variable_sites_all %>% 
    pivot_longer(starts_with("variable_site"), names_to = "site", values_to = "observed") %>% 
    left_join(consensus_vals, by = "site") %>% 
    group_by(.data$ace2_accession, .data$site) %>% 
    summarise(distance = get_site_dist(.data$most_common, .data$observed, 
                                       type = "grantham", ignore_na = TRUE),
              .groups = "drop") %>% 
    pivot_wider(id_cols = .data$ace2_accession, names_from = "site", values_from = "distance", 
                names_prefix = "dist_")
}