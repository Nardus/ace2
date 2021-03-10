## Summarise recorded meta-analysis data to one record per species

library(dplyr)
library(tidyr)
library(readxl)
library(readr)


# ---- Load ---------------------------------------------------------------------------------------
metadata <- read_excel("data/internal/infection_data.xlsx")

ace2_accessions <- read_excel("data/internal/ace2_accessions.xlsx")


# ---- Check data ---------------------------------------------------------------------------------
with(metadata, {
  stopifnot(all(experiment_type %in% c("cell", "animal")))
  stopifnot(all(infection_type %in% c("natural", "experimental", "modified cell", "whole cell")))
  
  stopifnot(all(infection_type[experiment_type == "animal"] %in% c("natural", "experimental")))
  stopifnot(all(infection_type[experiment_type == "cell"] %in% c("whole cell", "modified cell")))
  
  stopifnot(all(virus_name %in% c("SARS-CoV-2", "SARS-CoV-1", 
                                  "SARS-like CoV", "ACE2-utilizing SARS-like CoV")))
  
  stopifnot(all(infected %in% c(0, 1, NA_integer_)))
  stopifnot(all(shedding_virus %in% c(0, 1, NA_integer_)))
  stopifnot(all(transmission %in% c(0, 1, NA_integer_)))
  
  stopifnot(all(infected[infection_type == "natural"] == 1))
  
  stopifnot(all(species %in% ace2_accessions$species))
})


with(ace2_accessions, {
  stopifnot(all(species %in% metadata$species))
  stopifnot(length(unique(species)) == length(species)) # Duplicates will cause issues in join below
})



# ---- Relevant columns ---------------------------------------------------------------------------
# Note: from here onwards, shedding is formally defined as "shedding infectious virus"; detection
# of RNA at external sites is ignored
ace2_accessions <- ace2_accessions %>% 
  select(.data$species, .data$ace2_accession)

# Columns kept in all datasets below:
base_columns <- c("species", "experiment_type", "infection_type", "virus_name")

metadata <- metadata %>% 
  rename(shedding = .data$shedding_virus) %>% 
  select(all_of(c(base_columns, "infected", "shedding", "transmission"))) %>% 
  mutate(across(all_of(c("infected", "shedding", "transmission")), as.logical))


# ---- Evidence levels ----------------------------------------------------------------------------
# Evidence quality has 4 levels, with level 1 being most reliable:
#     1. Natural infection observed 
#     2. Experimental infection attempted (whole animal)
#     3. Experimental infection attempted (regular cell line)
#     4. Experimental infection attempted (cell line engineered to express ACE2 of a given species)
metadata <- metadata %>% 
  mutate(evidence_level = case_when(.data$experiment_type == "animal" & .data$infection_type == "natural" ~ 1L,
                                    .data$experiment_type == "animal" & .data$infection_type == "experimental" ~ 2L,
                                    .data$experiment_type == "cell" & .data$infection_type == "whole cell" ~ 3L,
                                    .data$experiment_type == "cell" & .data$infection_type == "modified cell" ~ 4L,
                                    TRUE ~ NA_integer_))

base_columns <- c(base_columns, "evidence_level")


# ---- Infectivity --------------------------------------------------------------------------------
# Summarise infection evidence:
# - Conflicting evidence from less reliable evidence levels are ignored if evidence is already 
#   known from a more reliable level. Similarly, conflicting evidence from the same level will be 
#   ignored, so we are simply asking "has it _ever_ been reported to be infected?"
# - Which virus was involved is currently ignored, but this is recorded for potential downstream
#   use
#
# - Animals recorded as shedding or transmitting must have been infected

# Summarise evidence levels supporting the fact that response_var takes on a given focal_value:
summarise_evidence <- function(evidence_levels, response_var, focal_value) {
   evidence_levels[response_var == focal_value] %>% 
    unique() %>% 
    sort() %>% 
    paste(collapse = ",")
}

infection_data <- metadata %>% 
  mutate(infected = case_when(is.na(.data$infected) & .data$shedding ~ TRUE,
                              is.na(.data$infected) & .data$transmission ~ TRUE,
                              TRUE ~ .data$infected)) %>% 
  
  select(all_of(base_columns), .data$infected) %>% 
  distinct() %>% 
  
  group_by(.data$species) %>% 
  summarise(all_evidence_true = summarise_evidence(.data$evidence_level, .data$infected, TRUE),
            all_evidence_false = summarise_evidence(.data$evidence_level, .data$infected, FALSE),
            viruses_true = summarise_evidence(.data$virus_name, .data$infected, TRUE),
            viruses_false = summarise_evidence(.data$virus_name, .data$infected, FALSE),
            evidence_level = if_else(any(.data$infected),
                                     min(c(100L, .data$evidence_level[.data$infected])),  # 100L prevents warnings when all(.data$infected) == FALSE
                                     min(.data$evidence_level)),
            infected = any(.data$infected),
            .groups = "drop")


# ---- Shedding -----------------------------------------------------------------------------------
# Same procedure as for infectivity data - asking "has shedding _ever_ been observed, and at what
# evidence level?". 
# 
# - transmission means shedding must have happened
# - animals recorded as not infective should not provide negative data here (otherwise, infection
#   and shedding become confounded)
shedding_data <- metadata %>% 
  mutate(shedding = case_when(!is.na(.data$shedding) & !is.na(.data$transmission) ~ .data$shedding | .data$transmission,
                              !is.na(.data$shedding) ~ .data$shedding,
                              !is.na(.data$transmission) ~ .data$transmission,
                              TRUE ~ NA),
         shedding = if_else(!.data$infected & !.data$shedding, NA, .data$shedding)) %>% 
  
  select(all_of(base_columns), .data$shedding) %>% 
  filter(!is.na(.data$shedding)) %>% 
  distinct() %>% 
  group_by(.data$species) %>% 
  summarise(all_evidence_true = summarise_evidence(.data$evidence_level, .data$shedding, TRUE),
            all_evidence_false = summarise_evidence(.data$evidence_level, .data$shedding, FALSE),
            viruses_true = summarise_evidence(.data$virus_name, .data$shedding, TRUE),
            viruses_false = summarise_evidence(.data$virus_name, .data$shedding, FALSE),
            evidence_level = if_else(any(.data$shedding),
                                     min(c(100L, .data$evidence_level[.data$shedding])),
                                     min(.data$evidence_level)),
            shedding = any(.data$shedding),
            .groups = "drop")


# ---- Convert column types -----------------------------------------------------------------------
# Response variables converted to character, because values should be usable as valid column names
infection_data <- infection_data %>% 
  mutate(infected = if_else(.data$infected, "True", "False"))

shedding_data <- shedding_data %>% 
  mutate(shedding = if_else(.data$shedding, "True", "False"))


# ---- Add ACE2 accessions ------------------------------------------------------------------------
ace2_accessions <- ace2_accessions %>% 
  filter(!is.na(.data$ace2_accession))

infection_data <- inner_join(infection_data, ace2_accessions, by = "species")
shedding_data <- inner_join(shedding_data, ace2_accessions, by = "species")


# ---- Output -------------------------------------------------------------------------------------
write_rds(infection_data, "data/calculated/cleaned_infection_data.rds")
write_rds(shedding_data, "data/calculated/cleaned_shedding_data.rds")
