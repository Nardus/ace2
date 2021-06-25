## Functions for dealing with AAindex (https://www.genome.jp/aaindex/) files
# - These files contain physicochemical properties for individual amino acids

#' Read an AAindex file and return a lookup-table of amino acid properties
#' 
#' Note that this is only compatible with AAindex1 files, not AAindex2, etc
#' 
#' @param path: path to an AAindex file
#' 
read_aaindex <- function(path) {
  # Order is always the same:
  aa_order <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
                "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  
  # Read / parse
  lines <- readLines(path)
  data_start <- which(startsWith(lines, "I")) + 1
  data_lines <- lines[data_start:(data_start+1)] 
  
  values <- data_lines %>% 
    trimws() %>% 
    strsplit("[[:space:]]+") %>%  # Number of spaces used as separator varies
    unlist() %>% 
    as.numeric()
  
  # Return NA when for gap characters:
  values <- c(values, NA)
  names(values) <- c(aa_order, "-")
  
  values
}


#' Read multiple AAindex files and return a data frame
#' 
#' @param index_names: vector of index file names (minus the extension)
#' @param base_path: path to the folder containing AAindex files to read
#' 
read_aaindex_multi <- function(index_names, base_path) {
  paths <- paste0(base_path, "/", index_names, ".txt")
  
  properties <- lapply(paths, read_aaindex) %>% 
    do.call(cbind, .)
  
  colnames(properties) <- index_names
  
  properties %>% 
    as.data.frame() %>% 
    rownames_to_column("AA")
}
