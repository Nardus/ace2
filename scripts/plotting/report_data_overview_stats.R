## General data overview stats to mention in text

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(readxl)
  library(ggplot2)
  
  source("scripts/plotting/plotting_constants.R")
})


infection_data <- readRDS("data/calculated/cleaned_infection_data.rds")
shedding_data <- readRDS("data/calculated/cleaned_shedding_data.rds")
shedding_raw <- read_excel("data/internal/infection_data.xlsx")

# ---- Infection data ------------------------------------------------------------------------------
# Viruses
# - Direct observations
l1_data <- infection_data %>% 
  filter(.data$evidence_level == 1)

cov1 <- grepl("SARS-CoV-1", l1_data$viruses_true, fixed = TRUE)
cov2 <- grepl("SARS-CoV-2", l1_data$viruses_true, fixed = TRUE)
cov_multi <- grepl(",", l1_data$viruses_true)

sprintf("%i natural infections involve SARS-CoV-2", sum(cov2))
sprintf("%i of these also infected by other sarbecoviruses, primarily SARS-CoV (N=%i)", 
        sum(cov2 & cov_multi),
        sum(cov2 & cov1))
cat("\n")

# - Experimental infections
l2_data <- infection_data %>% 
  filter(.data$evidence_level == 2)

cov2_pos <- grepl("SARS-CoV-2", l2_data$viruses_true, fixed = TRUE)
cov2_neg <- grepl("SARS-CoV-2", l2_data$viruses_false, fixed = TRUE)
cov_multi_pos <- grepl(",", l2_data$viruses_true)
cov_multi_all <- cov_multi_pos | grepl(",", l2_data$viruses_false)

sprintf("%i positive experimental infections involve SARS-CoV-2", sum(cov2_pos))
sprintf("%i negative experimental infections involve SARS-CoV-2", sum(cov2_neg))
sprintf("%i of %i hosts tested for >1 virus was susceptible for >1 virus", sum(cov_multi_pos), sum(cov_multi_all))


# Evidence levels
print("Infection data by evidence level")
infection_data %>% 
  mutate(evidence_label = EVIDENCE_LABELS[.data$evidence_level]) %>% 
  group_by(.data$evidence_level, .data$evidence_label) %>% 
  summarise(N = n(),
            Positive = sum(.data$infected == "True"),
            Negative = sum(.data$infected == "False"),
            .groups = "keep") %>% 
  mutate(Proportion = .data$Positive/.data$N,
         CI_lower = binom.test(.data$Positive, .data$N)$conf.int[1],
         CI_upper = binom.test(.data$Positive, .data$N)$conf.int[2]) %>% 
  print()


# ---- Shedding data ------------------------------------------------------------------------------
# RNA
rna_shedding <- shedding_raw %>% 
  group_by(.data$species) %>% 
  summarise(data_reported = any(!is.na(.data$shedding_rna)),
            shedding_rna = any(.data$shedding_rna == 1, na.rm = TRUE))

sprintf("RNA shedding detected in %i of %i species for which data are available",
        sum(rna_shedding$shedding_rna),
        sum(rna_shedding$data_reported))

# Evidence levels (infectious virus shedding)
print("Shedding data by evidence level (infectious virus only)")
shedding_data %>% 
  mutate(evidence_label = EVIDENCE_LABELS[.data$evidence_level]) %>% 
  group_by(.data$evidence_level, .data$evidence_label) %>% 
  summarise(N = n(),
            Positive = sum(.data$shedding == "True"),
            Negative = sum(.data$shedding == "False"),
            .groups = "keep") %>% 
  mutate(Proportion = .data$Positive/.data$N,
         CI_lower = binom.test(.data$Positive, .data$N)$conf.int[1],
         CI_upper = binom.test(.data$Positive, .data$N)$conf.int[2]) %>% 
  print()
