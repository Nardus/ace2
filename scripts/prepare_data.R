## Summarise recorded meta-analysis data to one record per species

library(dplyr)
library(tidyr)
library(readr)


# ---- Load ---------------------------------------------------------------------------------------
metadata <- read_csv("data/internal/ace2_metadata.csv", 
                     na = c("", "?", "NA"),
                     col_types = cols(.default = col_character(),
                                      `Outcome of infection [1=infection, 0=no infection]` = col_double(),
                                      `Shedding [0=no, 1=PCR,isolation,or evidence of natural or experimental transmission, ?=unknown or not discussed]` = col_double(),
                                      `Transmission [1=yes, 0=no, ?=unknown]` = col_double(),
                                      `Delta G (kcal/mol)` = col_double(),
                                      `Binding Affinity (M)` = col_double()))


# ---- Check data ---------------------------------------------------------------------------------
check_acc <- metadata %>% 
  distinct(.data$`Host genus`, .data$`Host species`, .data$`ACE2 sequence`) %>% 
  add_count(.data$`ACE2 sequence`) %>% 
  filter(n > 1)

if (nrow(check_acc) > 1)
  warning(length(unique(check_acc$`ACE2 sequence`)), 
          " ACE2 accession numbers are duplicated across >1 species")

check_spp <- metadata %>% 
  distinct(.data$`Host genus`, .data$`Host species`, .data$`ACE2 sequence`) %>% 
  add_count(.data$`Host genus`, .data$`Host species`) %>% 
  filter(n > 1)

if (nrow(check_spp) > 1)
  warning(n_distinct(check_spp$`Host genus`, check_spp$`Host species`), 
          " species linked to different ACE2 accesion numbers in different rows")


# ---- Relevant columns ---------------------------------------------------------------------------
metadata <- metadata %>% 
  rename(study_type = .data$`Animal or cell study`,
         infection_type = .data$`Infection type [natural or experimental inoculation]`,
         ace2_accession = .data$`ACE2 sequence`) %>% 
  mutate(species = paste(.data$`Host genus`, .data$`Host species`),
         infected = .data$`Outcome of infection [1=infection, 0=no infection]` == 1,
         shedding = .data$`Shedding [0=no, 1=PCR,isolation,or evidence of natural or experimental transmission, ?=unknown or not discussed]` == 1,
         transmission = .data$`Transmission [1=yes, 0=no, ?=unknown]` == 1)

# Columns kept in all datasets below:
base_columns <- c("species", "study_type", "infection_type", "ace2_accession")


# ---- Evidence levels ----------------------------------------------------------------------------
# Evidence quality has 4 levels, with level 1 being most reliable:
#     1. Natural infection observed 
#     2. Experimental infection attempted (whole animal)
#     3. Experimental infection attempted (primary cells)
#     4. (TODO) Experimental infection attempted (cell line engineered to express ACE2 of a given species)
warning("TODO: cells engineered to express ACE2 are not currently identifiable")

metadata <- metadata %>% 
  mutate(evidence_level = case_when(.data$study_type == "animal" & .data$infection_type == "natural" ~ 1L,
                                    .data$study_type == "animal" & .data$infection_type == "experimental" ~ 2L,
                                    .data$study_type == "cell" ~ 3L,
                                    TRUE ~ NA_integer_))

stopifnot(all(metadata$infected[metadata$infection_type == "natural"] == TRUE))
base_columns <- c(base_columns, "evidence_level")


# ---- Infectivity --------------------------------------------------------------------------------
# Summarise infection evidence:
# - Conflicting evidence from less reliable evidence levels are ignored if evidence is already 
#   known from a more reliable level. Similarly, conflicting evidence from the same level will be 
#   ignored, so we are simply asking "has it _ever_ been reported to be infected?"
infection_data <- metadata %>% 
  select(all_of(base_columns), .data$infected) %>% 
  distinct() %>% 
  
  group_by(.data$species, .data$ace2_accession) %>% 
  summarise(evidence_level = if_else(any(.data$infected),
                                     min(c(100L, .data$evidence_level[.data$infected])),  # 100L prevents warnings when all(.data$infected) == FALSE
                                     min(.data$evidence_level)),
            infected = any(.data$infected),
            .groups = "drop")

stopifnot(length(unique(infection_data$species)) == nrow(infection_data))  # Expect one ace2_accession per species


# ---- Shedding -----------------------------------------------------------------------------------
# Same procedure as for infectivity data - asking "has shedding _ever_ been observed, and at what
# evidence level?". However, transmission also means shedding must have happened, so evidence of
# either is counted as evidence that shedding occurred.
shedding_data <- metadata %>% 
  mutate(shedding = case_when(!is.na(.data$shedding) & !is.na(.data$transmission) ~ .data$shedding | .data$transmission,
                              !is.na(.data$shedding) ~ .data$shedding,
                              !is.na(.data$transmission) ~ .data$transmission,
                              TRUE ~ NA)) %>% 
  
  select(all_of(base_columns), .data$shedding) %>% 
  filter(!is.na(.data$shedding)) %>% 
  distinct() %>% 
  group_by(.data$species, .data$ace2_accession) %>% 
  summarise(evidence_level = if_else(any(.data$shedding),
                                     min(c(100L, .data$evidence_level[.data$shedding])),
                                     min(.data$shedding)),
            shedding = any(.data$shedding),
            .groups = "drop")

stopifnot(length(unique(shedding_data$species)) == nrow(shedding_data))


# ---- Transmission -------------------------------------------------------------------------------
transmission_data <- metadata %>% 
  select(all_of(base_columns), .data$transmission) %>% 
  filter(!is.na(.data$transmission)) %>% 
  distinct() %>% 
  group_by(.data$species, .data$ace2_accession) %>% 
  summarise(evidence_level = if_else(any(.data$transmission),
                                     min(c(100L, .data$evidence_level[.data$transmission])),
                                     min(.data$transmission)),
            transmission = any(.data$transmission),
            .groups = "drop")



# ---- Convert column types -----------------------------------------------------------------------
# Response variables converted to character, becuase values should be usable as valid column names
infection_data <- infection_data %>% 
  mutate(infected = if_else(.data$infected, "True", "False"))

shedding_data <- shedding_data %>% 
  mutate(shedding = if_else(.data$shedding, "True", "False"))

transmission_data <- transmission_data %>% 
  mutate(transmission = if_else(.data$transmission, "True", "False"))


# ---- Output -------------------------------------------------------------------------------------
write_rds(infection_data, "data/calculated/cleaned_infection_data.rds")
write_rds(shedding_data, "data/calculated/cleaned_shedding_data.rds")
write_rds(transmission_data, "data/calculated/cleaned_transmission_data.rds")