## Utility functions for calculating amino-acid distances
require(dplyr)
require(Matrix)

# ---- Load raw distance matrices -----------------------------------------------------------------
# Grantham, derived from the AAindex database 
#   (https://www.genome.jp/dbget-bin/www_bget?aaindex:GRAR740104)
grantham_mat <- read.csv("data/internal/distance_matrices/grantham.csv")
rownames(grantham_mat) <- colnames(grantham_mat)
grantham_mat <- as.matrix(grantham_mat)
grantham_mat <- forceSymmetric(grantham_mat, uplo = "L")
grantham_mat <- as.matrix(grantham_mat)

# BLOSUM100, downloaded from ftp://ftp.ncbi.nih.gov/blast/matrices/:
blosum_mat <- read.table("data/internal/distance_matrices/BLOSUM100", check.names = FALSE)
blosum_mat <- as.matrix(blosum_mat)

# WAG, derived from https://www.ebi.ac.uk/goldman-srv/WAG/:
wag_mat <- read.csv("data/internal/distance_matrices/wag.csv")
rownames(wag_mat) <- colnames(wag_mat)
wag_mat <- as.matrix(wag_mat)
diag(wag_mat) <- 0
wag_mat <- forceSymmetric(wag_mat, uplo = "L")
wag_mat <- as.matrix(wag_mat)


# EX, derived from PMID:PMC1449787's supplement:
experimental_mat <- read.csv("data/internal/distance_matrices/experimental_exchangeability.csv", 
                             row.names = 1, na.strings = ".")
experimental_mat <- as.matrix(experimental_mat)
diag(experimental_mat) <- 1
experimental_mat <- 1 - experimental_mat # convert to distance (original measures "exchangeability")



# ---- Functions ----------------------------------------------------------------------------------
#' Retrieve a distance score for a single site.
#' 
#' @param aa1: the wild-type or origin amino acid
#' @param aa2: the novel/mutated amino acid
#' @param type: type of distance (see Details)
#' @param symmetric: should scores be symmetric (see details)?
#' @param ignore_gaps: should gaps be treated as missing data?
#' @output double
#' 
#' @details
#' Four distance metrics are supported:
#' - `grantham`: Physicochemical distance of Grantham, 1974 (PMID: 4843792), which is based on 
#'               molecular volume, polarity, and composition
#' - `blosum`: BLOSUM100, giving an alignment score (or odds of relatedness) based on substitution 
#'             probabilities among HIGHLY conserved protein homologs, while accounting for how 
#'             common different amino acids are (PMID: 1438297)
#' - `wag`: Evolutionary distance of Whelan and Goldman (PMID: 11319253), based on empirically 
#'          observed substitution rates (WAG)
#' - `experimental`: The experimental exchangeability (EX) distance of Yampolsky & Stoltzfus, 2005 
#'                   (PMID: PMC1449787)
#' 
#' By default, a symmetric score is returned, by averaging over scores in both directions when 
#' needed (i.e., the mean of d(aa1, aa2) and d(aa2, aa1) is returned). Note that of the currently 
#' implemented distance types, only `experimental` is actually asymmetric.
#' 
#' If `ignore_gaps=TRUE`, gaps are (silently) ignored and an NA returned. Otherwise, the distance
#' matrix will be checked for a gap character and a warning issued if it is not present (in
#' which case the returned result is again NA).
#' 
get_site_dist <- function(aa1, aa2, type = c("grantham", "blosum", "wag", "experimental"), 
                          symmetric = TRUE, ignore_gaps = TRUE) {
  if (!nchar(aa1) == 1 & nchar(aa2) == 1)
    stop("Single amino acid character expected")
  
  type <- match.arg(type)
  
  if (ignore_gaps & "-" %in% c(aa1, aa2))
    return(NA_real_)
  
  distance_matrix <- switch(type,
                            grantham = grantham_mat,
                            blosum = blosum_mat,
                            wag = wag_mat,
                            experimental = experimental_mat)
  
  if (!aa1 %in% colnames(distance_matrix) | !aa2 %in% rownames(distance_matrix)) {
    warning("Amino acid not found in distance matrix, returing NA")
    return(NA_real_)
  }
  
  if (symmetric & type == "experimental") {
    forward <- distance_matrix[aa1, aa2]
    reverse <- distance_matrix[aa2, aa1]
    return(mean(c(forward, reverse)))
  }
    
  as.double(distance_matrix[aa1, aa2])
}


#' Get the distance between two sequences.
#' 
#' @param from: the wild-type or origin sequence, as a character vector
#' @param to: the novel/mutated sequence, as a character vector
#' @param type: type of distance (see `get_site_dist`)
#' @param symmetric: should scores be symmetric (see `get_site_dist`)?
#' @output double
#' 
get_sequence_distance <- function(from, to, type = c("grantham", "blosum", "wag", "experimental"),
                                  symmetric = TRUE) {
  site_distances <- mapply(get_site_dist,
         aa1 = from,
         aa2 = to,
         MoreArgs = list(type = type,
                         symmetric = symmetric))
  
  sum(site_distances, na.rm = TRUE)
}


#' Calculate symmetric pairwise distances between all sequences in an alignment.
#' 
#' @param alignment: list of character vectors representing an alignment (as produced by 
#'                   e.g. `seqinr::read_fasta()`)
#' @param type: type of distance (see `get_site_dist`)
#' @output a symmetric n x n distance matrix
#' 
#' @details 
#'
get_all_distances <- function(alignment, type = c("grantham", "blosum", "wag", "experimental")) {
  all_combinations <- expand.grid(names(alignment), names(alignment), stringsAsFactors = FALSE)
  
  dist_vec <- mapply(get_sequence_distance, 
                     from = alignment[all_combinations[, 1]],
                     to = alignment[all_combinations[, 2]],
                     MoreArgs = list(type = type,
                                     symmetric = TRUE),
                     USE.NAMES = FALSE)
  
  dist_mat <- matrix(dist_vec, nrow = length(alignment))
  colnames(dist_mat) <- rownames(dist_mat) <- names(alignment)
  
  dist_mat
}
