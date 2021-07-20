## Recipe steps to calculate training-set specific features on the fly
#  - This is needed because some variables depend on the set of viruses considered "known" 
#    at given point in time, which for model training and evaluation is the training set
#  - See feature_calc_utils.R for the actual feature calculations

#' Calculate distaince features specific to a give training set
#' To use, prepare the recipe with the training data, then bake with each unique dataset
#' 
#' `options` should be a list containing the following named data frames:
#'   - all_pairwise_dists: data frame containing pairwise distances in long format
#'   - all_variable_sites: data frame of variable sites for all hosts
#'   - metadata:           meta-data listing accessions along with a "label" column specifying
#'                         the response being predicted (with character values "True" and 
#'                         "False" indicating positives and negatives, respectively)
step_distance_features <- function(
  recipe, 
  ..., 
  role = NA, 
  trained = FALSE, 
  train_accessions = NULL,
  options = list(all_pairwise_dists = NULL, 
                 all_variable_sites = NULL, 
                 metadata = NULL),
  skip = FALSE,
  id = rand_id("distance_features")
) {
  terms <- ellipse_check(...) 
  
  add_step(
    recipe, 
    step_distance_features_new(
      terms = terms, 
      trained = trained,
      role = role, 
      train_accessions = train_accessions,
      options = options,
      skip = skip,
      id = id
    )
  )
}


step_distance_features_new <- 
  function(terms, role, trained, train_accessions, options, skip, id) {
    step(
      subclass = "distance_features", 
      terms = terms,
      role = role,
      trained = trained,
      train_accessions = train_accessions,
      options = options,
      skip = skip,
      id = id
    )
  }


prep.step_distance_features <- function(x, training, info = NULL, ...) {
  # Record accessions from the training set
  train_accessions = unique(training$ace2_accession)
  
  # Use the constructor to return the updated object, setting `trained` to TRUE:
  step_distance_features_new(
    terms = x$terms, 
    trained = TRUE,
    role = x$role, 
    train_accessions = train_accessions,
    options = x$options,
    skip = x$skip,
    id = x$id
  )
}


bake.step_distance_features <- function(object, new_data, ...) {
  # Use training accessions recorded by prep() step to adjust calculations to this training set
  pairwise_dists_training <- object$options$all_pairwise_dists %>% 
    filter(.data$other_seq %in% object$train_accessions)
  
  variable_sites_training <- object$options$all_variable_sites %>% 
    filter(.data$ace2_accession %in% object$train_accessions)
  
  # Calculate new features
  closest_positive <- get_dist_to_closest_positive(pairwise_dist_data = pairwise_dists_training, 
                                                   metadata = object$options$metadata)
  
  consensus_dists <- get_consensus_dist(variable_sites_training = variable_sites_training, 
                                        variable_sites_all = object$options$all_variable_sites, 
                                        metadata = object$options$metadata)
  
  # Attach to input data and return
  new_data %>% 
    left_join(closest_positive, by = "ace2_accession") %>% 
    left_join(consensus_dists, by = "ace2_accession") %>% 
    as_tibble()
}


# Globals which should be exported to run above steps in parallel:
recipe_globals <- c("GRANTHAM_MAT", "BLOSUM_MAT", "WAG_MAT", "EXPERIMENTAL_MAT",
                    ".check_alignment", 
                    "most_common_value",
                    "get_site_dist",
                    "get_site_dist_vectorised",
                    "get_dist_to_closest_positive",
                    "get_consensus_dist",
                    "step_distance_features_new",
                    "prep.step_distance_features",
                    "bake.step_distance_features")
