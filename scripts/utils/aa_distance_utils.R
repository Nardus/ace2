## Utility functions for calculating amino-acid distances
require(Matrix)
require(parallel)

# ---- Load raw distance matrices -----------------------------------------------------------------
# Grantham, derived from the AAindex database 
#   (https://www.genome.jp/dbget-bin/www_bget?aaindex:GRAR740104)
GRANTHAM_MAT <- read.csv("data/internal/distance_matrices/grantham.csv")
rownames(GRANTHAM_MAT) <- colnames(GRANTHAM_MAT)
GRANTHAM_MAT <- as.matrix(GRANTHAM_MAT)
GRANTHAM_MAT <- forceSymmetric(GRANTHAM_MAT, uplo = "L")
GRANTHAM_MAT <- as.matrix(GRANTHAM_MAT)

# BLOSUM100, downloaded from ftp://ftp.ncbi.nih.gov/blast/matrices/:
BLOSUM_MAT <- read.table("data/internal/distance_matrices/BLOSUM100", check.names = FALSE)
BLOSUM_MAT <- as.matrix(BLOSUM_MAT)

# WAG, derived from https://www.ebi.ac.uk/goldman-srv/WAG/:
WAG_MAT <- read.csv("data/internal/distance_matrices/wag.csv")
rownames(WAG_MAT) <- colnames(WAG_MAT)
WAG_MAT <- as.matrix(WAG_MAT)
diag(WAG_MAT) <- 0
WAG_MAT <- forceSymmetric(WAG_MAT, uplo = "L")
WAG_MAT <- as.matrix(WAG_MAT)


# EX, derived from PMID:PMC1449787's supplement:
EXPERIMENTAL_MAT <- read.csv("data/internal/distance_matrices/experimental_exchangeability.csv", 
                             row.names = 1, na.strings = ".")
EXPERIMENTAL_MAT <- as.matrix(EXPERIMENTAL_MAT)
diag(EXPERIMENTAL_MAT) <- 1
EXPERIMENTAL_MAT <- 1 - EXPERIMENTAL_MAT # convert to distance (original measures "exchangeability")



# ---- Functions ----------------------------------------------------------------------------------
#' INTERNAL: Check for and warn about unknown characters
.check_alignment <- function(alignment, type = c("grantham", "blosum", "wag", "experimental"), 
                             ignore_na = FALSE,
                             matrices = list(grantham = GRANTHAM_MAT,
                                             blosum = BLOSUM_MAT,
                                             wag = WAG_MAT,
                                             experimental = EXPERIMENTAL_MAT)) {
  type <- match.arg(type)
  distance_matrix <- matrices[[type]]
  alignment_chars <- unique(unlist(alignment))
  alignment_chars <- alignment_chars[alignment_chars != "-"]
  
  if (ignore_na)
    alignment_chars <- alignment_chars[!is.na(alignment_chars)]
  
  
  unknown_chars <- alignment_chars[!alignment_chars %in% colnames(distance_matrix)]
  
  if (length(unknown_chars) > 0)
    warning(sprintf("Amino acid(s) not found in %s distance matrix will be ignored: %s", 
                    type, 
                    paste(unknown_chars, collapse = ", ")))
}


#' Retrieve a distance score for a single site.
#' 
#' @param aa1: the wild-type or origin amino acid
#' @param aa2: the novel/mutated amino acid
#' @param type: type of distance (see Details)
#' @param symmetric: should scores be symmetric (see details)?
#' @param ignore_gaps: should gaps be treated as missing data?
#' @param check: should input amino acids be checked for validity?
#' @param ignore_na: when checking amino acids, should NA's generate a warning?
#' @param matrices: list of lookup matrixes specifying distances between all bases
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
                          symmetric = TRUE, ignore_gaps = TRUE, check = TRUE, ignore_na = FALSE,
                          matrices = list(grantham = GRANTHAM_MAT,
                                          blosum = BLOSUM_MAT,
                                          wag = WAG_MAT,
                                          experimental = EXPERIMENTAL_MAT)) {
  if (!nchar(aa1) == 1 & nchar(aa2) == 1)
    stop("Single amino acid character expected")
  
  type <- match.arg(type)
  distance_matrix <- matrices[[type]]
  
  if (ignore_gaps & "-" %in% c(aa1, aa2))
    return(NA_real_)
  
  if (!aa1 %in% colnames(distance_matrix) | !aa2 %in% rownames(distance_matrix)) {
    if (check)
      .check_alignment(c(aa1, aa2), type = type, matrices = matrices, ignore_na = ignore_na)
    
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
#' @param check: should input amino acids be checked for validity?
#' @output double
#' 
get_sequence_distance <- function(from, to, type = c("grantham", "blosum", "wag", "experimental"),
                                  symmetric = TRUE, check = TRUE) {
  
  if (check)
    .check_alignment(list(from, to), type = type)
  
  site_distances <- mapply(get_site_dist,
         aa1 = from,
         aa2 = to,
         MoreArgs = list(type = type,
                         symmetric = symmetric,
                         check = FALSE))
  
  sum(site_distances, na.rm = TRUE)
}


#' Calculate symmetric pairwise distances between all sequences in an alignment.
#' 
#' @param alignment: list of character vectors representing an alignment (as produced by 
#'                   e.g. `seqinr::read_fasta()`)
#' @param type: type of distance (see `get_site_dist`)
#' @param cores: number of parallel cores to use
#' @param check: should input amino acids be checked for validity?
#' @output a symmetric n x n distance matrix
#' 
#' @details 
#'
get_all_distances <- function(alignment, type = c("grantham", "blosum", "wag", "experimental"), 
                              cores = 1, check = TRUE) {
  if (check)
    .check_alignment(alignment, type = type)
  
  all_combinations <- expand.grid(names(alignment), names(alignment), stringsAsFactors = FALSE)
  
  dist_vec <- mcmapply(get_sequence_distance, 
                       from = alignment[all_combinations[, 1]],
                       to = alignment[all_combinations[, 2]],
                       MoreArgs = list(type = type,
                                       symmetric = TRUE,
                                       check = FALSE),
                       USE.NAMES = FALSE, 
                       mc.cores = cores)
  
  dist_mat <- matrix(dist_vec, nrow = length(alignment))
  colnames(dist_mat) <- rownames(dist_mat) <- names(alignment)
  
  dist_mat
}
