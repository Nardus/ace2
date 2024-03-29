## Build human-readable supplementary table describing:
#   - Training data & accessions
#   - Shedding data

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(readr)
  library(readxl)
  library(bibtex)
  library(writexl)
})

# Raw data
infection_metadata <- read_excel("data/internal/infection_data.xlsx")
ace2_metadata <- read_csv("data/internal/ace2_accessions.csv")

bib <- read.bib("data/internal/data_citations.bib")

# Cleaned data
infection_data <- readRDS("data/calculated/cleaned_infection_data.rds")
shedding_data <- readRDS("data/calculated/cleaned_shedding_data.rds")


# ---- References ----------------------------------------------------------------------------------
# Make human-readable references (as a seperate worksheet in the final xlsx file)
extract_citation <- function(x) {
    # Build an author-year citation from a bibentry record
    authors <- x$author$family
    authors <- sapply(authors, function(x) paste(x, collapse = " ")) # Join multi-part surnames
    
    citation_authors <- if_else(length(authors) > 2,
                                paste(authors[[1]], "et al."),
                                paste(authors, collapse = " & "))
                                
    citation_year <- x$year
    
    paste(citation_authors, citation_year)
}

clean_reference <- function(x) {
    # Create a human-readable, plain text reference from a bibentry record
    x = format(x, style = "text") %>% 
        str_remove_all("[“”*\\\\]") %>% 
        str_replace_all("_([[:alnum:]&\\s]+)_", "\\1") %>% 
        str_replace_all("\\s-([[:alnum:],\\s]+)-\\s", "\\1") %>% 
        str_replace_all("\\n", " ") %>% 
        str_replace(" \\(URL: .+\\)\\.$", ".")
    
    tibble(reference = x)
}

stopifnot(all(infection_metadata$citation_key %in% names(bib)))

full_citation <- sapply(bib, extract_citation)

reference_df <- sapply(bib, clean_reference, USE.NAMES = TRUE, simplify = FALSE) %>% 
    bind_rows(.id = "citation_key") %>% 
    filter(.data$citation_key %in% infection_metadata$citation_key) %>% 
    mutate(citation = full_citation[.data$citation_key]) %>% 
    select(.data$citation, .data$reference) %>% 
    arrange(.data$citation)


# ---- Infection data ------------------------------------------------------------------------------
# Add formatted citations
collapsed_citations <- infection_metadata %>% 
    arrange(.data$citation_key) %>%
    mutate(citation = full_citation[.data$citation_key]) %>% 
    group_by(.data$species) %>% 
    summarise(citation = paste(citation, collapse = "; "),
              .groups = "drop")

infection_data <- infection_data %>% 
    left_join(collapsed_citations, by = "species")


# Add accession replacement information
stopifnot(all(infection_data$ace2_accession %in% ace2_metadata$ace2_accession))

replacement_data <- ace2_metadata %>% 
    mutate(ace2_subsituted = !.data$ace2_species == "this_species",
           ace2_subsituted_species = if_else(.data$ace2_subsituted, .data$ace2_species, ""),
           ace2_subsituted = if_else(.data$ace2_subsituted, "True", "False")) %>% 
    select(.data$species, .data$ace2_accession, 
           .data$ace2_subsituted, .data$ace2_subsituted_species)

infection_data <- infection_data %>% 
    left_join(replacement_data, by = c("species", "ace2_accession"))


# Human-readable evidence levels
format_evidence_levels <- function(x) {
    x %>% 
        str_replace("1", "Observed infection") %>% 
        str_replace("2", "Experimental infection") %>%
        str_replace("3", "Cell culture") %>%
        str_replace("4", "Cell culture (het-ACE2)") %>% 
        str_replace_all(",", "; ")
}

infection_data <- infection_data %>% 
    mutate(all_evidence_true = format_evidence_levels(.data$all_evidence_true),
           all_evidence_false = format_evidence_levels(.data$all_evidence_false),
           evidence_level = as.character(.data$evidence_level),
           evidence_level = format_evidence_levels(.data$evidence_level),
           viruses_true = str_replace_all(.data$viruses_true, ",", "; "),
           viruses_false = str_replace_all(.data$viruses_false, ",", "; "))


# Clean column names
infection_data <- infection_data %>% 
    select(.data$species, .data$infected, .data$viruses_true, .data$viruses_false,
           evidence_true = .data$all_evidence_true, 
           evidence_false = .data$all_evidence_false,
           .data$citation,
           .data$ace2_accession, .data$ace2_subsituted, .data$ace2_subsituted_species)


# Check all joins worked
stopifnot(nrow(infection_data) == n_distinct(infection_data$species))


# ---- Shedding data -------------------------------------------------------------------------------
shedding_data <- shedding_data %>% 
    left_join(collapsed_citations, by = "species") %>% 
    left_join(replacement_data, by = c("species", "ace2_accession")) %>% 
    
    mutate(all_evidence_true = format_evidence_levels(.data$all_evidence_true),
           all_evidence_false = format_evidence_levels(.data$all_evidence_false),
           evidence_level = as.character(.data$evidence_level),
           evidence_level = format_evidence_levels(.data$evidence_level),
           viruses_true = str_replace_all(.data$viruses_true, ",", "; "),
           viruses_false = str_replace_all(.data$viruses_false, ",", "; ")) %>% 
    
    select(.data$species, .data$shedding, .data$viruses_true, .data$viruses_false,
           evidence_true = .data$all_evidence_true, 
           evidence_false = .data$all_evidence_false,
           .data$citation,
           .data$ace2_accession, .data$ace2_subsituted, .data$ace2_subsituted_species)

# Check all joins worked
stopifnot(nrow(shedding_data) == n_distinct(shedding_data$species))


# ---- Output --------------------------------------------------------------------------------------
dir.create("output/si_tables/")

# Training data
training_data <- list(infection = infection_data,
                      shedding = shedding_data,
                      references = reference_df)

write_xlsx(training_data, "output/si_tables/supplement_training_data.xlsx")
