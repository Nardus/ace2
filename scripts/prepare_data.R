## Summarise recorded meta-analysis data to one record per species

library(dplyr)
library(tidyr)
library(readr)


# ---- Load ---------------------------------------------------------------------------------------
metadata <- read_csv("data/internal/ace2_metadata.csv", 
                     col_types = cols(.default = col_character(),
                                      `Outcome of infection [1=infection, 0=no infection]` = col_double(),
                                      `Grantham distance to human` = col_double(),
                                      `Grantham distance to closest related species` = col_double(),
                                      `Grantham distance to closest related species not including self` = col_double(),
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
infection_data <- metadata %>% 
  distinct(.data$`Host genus`, .data$`Host species`, .data$`Animal or cell study`, 
           .data$`ACE2 sequence`, .data$`Infection type [natural or experimental inoculation]`, 
           .data$`Outcome of infection [1=infection, 0=no infection]`) %>% 
  mutate(species = paste(.data$`Host genus`, .data$`Host species`),
         infected = .data$`Outcome of infection [1=infection, 0=no infection]` == 1) %>% 
  select(.data$species,
         study_type = .data$`Animal or cell study`,
         infection_type = .data$`Infection type [natural or experimental inoculation]`,
         ace2_accession = .data$`ACE2 sequence`,
         .data$infected)


# ---- Infectivity --------------------------------------------------------------------------------
# Summarise infection evidence:
# - Evidence quality has xx levels, with level one being most reliable. Conflicting evidence from
#   later levels are ignored if evidence is already known from a more reliable level. Similarly,
#   conflicting evidence from the same level is ignored, so we are simply asking "has it _ever_ 
#   been reported to be infected?"
# - Levels are:
#     1. Natural infection observed 
#     2. Experimental infection attempted (whole animal)
#     3. Experimental infection attempted (primary cells)
#     4. (TODO) Experimental infection attempted (cell line engineered to express ACE2 of a given species)

warning("TODO: cells engineered to express ACE2 are not currently identifiable")

stopifnot(all(infection_data$infected[infection_data$infection_type == "natural"] == 1))

infection_data <- infection_data %>% 
  mutate(evidence_level = case_when(.data$study_type == "animal" & .data$infection_type == "natural" ~ 1L,
                                    .data$study_type == "animal" & .data$infection_type == "experimental" ~ 2L,
                                    .data$study_type == "cell" ~ 3L,
                                    TRUE ~ NA_integer_)) %>%  
  group_by(.data$species, .data$ace2_accession) %>% 
  summarise(evidence_level = if_else(any(.data$infected),
                                     min(c(100L, .data$evidence_level[.data$infected])),  # 100L prevents warnings when all(.data$infected) == FALSE
                                     min(.data$evidence_level)),
            infected = any(.data$infected),
            .groups = "drop")

stopifnot(length(unique(infection_data$species)) == nrow(infection_data))  # Expect one ace2_accession per species


# ---- Output -------------------------------------------------------------------------------------
write_rds(infection_data, "data/calculated/cleaned_infection_data.rds")