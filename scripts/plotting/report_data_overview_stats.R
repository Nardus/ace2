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
raw_data <- read_excel("data/internal/infection_data.xlsx")

taxonomy <- readRDS("data/calculated/taxonomy.rds")
mammals <- taxonomy %>% 
  filter(.data$class == "Mammalia") %>% 
  pull(.data$internal_name)

# ---- Infection data ------------------------------------------------------------------------------
# Overview
cat("\n\n")
sprintf("Infection data: %i susceptible, %i not susceptible",
        sum(infection_data$infected == "True"),
        sum(infection_data$infected == "False"))

sprintf("%.2f%% of records involve SARS-CoV-2 in some way (might not be maximum evidence level)",
        sum(str_detect(infection_data$viruses_true, "SARS-CoV-2") |
              str_detect(infection_data$viruses_false, "SARS-CoV-2")) / nrow(infection_data) * 100)

# Percentages
total <- binom.test(sum(infection_data$infected == "True"), nrow(infection_data))
mammal_only <- binom.test(sum(infection_data$infected == "True" & infection_data$species %in% mammals),
                          sum(infection_data$species %in% mammals))

cat("\n")
sprintf("In total %0.2f%% susceptible, (CI: %0.2f - %0.2f)",
        total$estimate * 100,
        total$conf.int[1] * 100,
        total$conf.int[2] * 100)

sprintf("Mammals only: %0.2f%% susceptible, (CI: %0.2f - %0.2f)",
        mammal_only$estimate * 100,
        mammal_only$conf.int[1] * 100,
        mammal_only$conf.int[2] * 100)

cat("\n")

# Viruses
# - Direct observations
l1_data <- infection_data %>% 
  filter(.data$evidence_level == 1)

cov1 <- grepl("SARS-CoV-1", l1_data$viruses_true, fixed = TRUE)
cov2 <- grepl("SARS-CoV-2", l1_data$viruses_true, fixed = TRUE)
cov_multi <- grepl(",", l1_data$viruses_true)

sprintf("%i observations of natural infections", nrow(l1_data))
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

sprintf("%i observations from experimental infections", nrow(l2_data))
sprintf("%i positive experimental infections involve SARS-CoV-2", sum(cov2_pos))
sprintf("%i negative experimental infections involve SARS-CoV-2", sum(cov2_neg))
sprintf("%i of %i hosts tested for >1 virus was susceptible for >1 virus", sum(cov_multi_pos), sum(cov_multi_all))


# Evidence levels
print("Infection data by evidence level (all)")
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
  ungroup() %>%
  select(-.data$evidence_level) %>% 
  print()

cat("\n\n(mammals only)")
infection_data %>% 
  filter(species %in% mammals) %>% 
  mutate(evidence_label = EVIDENCE_LABELS[.data$evidence_level]) %>% 
  group_by(.data$evidence_level, .data$evidence_label) %>% 
  summarise(N = n(),
            Positive = sum(.data$infected == "True"),
            Negative = sum(.data$infected == "False"),
            .groups = "keep") %>% 
  mutate(Proportion = .data$Positive/.data$N,
         CI_lower = binom.test(.data$Positive, .data$N)$conf.int[1],
         CI_upper = binom.test(.data$Positive, .data$N)$conf.int[2]) %>% 
  ungroup() %>%
  select(-.data$evidence_level) %>% 
  print()


# Hosts infected by ACE2-utilizing sarbecoviruses
cat("\n\nACE2-utilization:\n")
ace2_viruses <- c("SARS-CoV-2", "SARS-CoV-1", "ACE2-utilizing SARS-like CoV")

cat("Overall:\n")
ace2_used <- raw_data %>%
  group_by(.data$species) %>%
  summarise(ACE2_used = any(.data$virus_name %in% ace2_viruses))

sprintf("%i of %i all hosts (%.3f%%) are known to use ACE2",
        sum(ace2_used$ACE2_used),
        nrow(ace2_used),
        sum(ace2_used$ACE2_used) / nrow(ace2_used) * 100)

cat("Susceptible:\n")
susceptible_spp <- infection_data %>%
  filter(.data$infected == "True") %>%
  pull(.data$species)

ace2_used <- raw_data %>%
  filter(.data$species %in% susceptible_spp) %>%
  group_by(.data$species) %>%
  summarise(ACE2_used = any(.data$virus_name %in% ace2_viruses))

sprintf("%i of %i susceptible hosts (%.3f%%) are known to use ACE2",
        sum(ace2_used$ACE2_used),
        nrow(ace2_used),
        sum(ace2_used$ACE2_used) / nrow(ace2_used) * 100)


# SARS-CoV-1 vs. 2
cat("\n\nHost range similarity (all evidence levels):\n")
both_tested <- infection_data %>% 
  filter((grepl("SARS-CoV-1", .data$viruses_true, fixed = TRUE) | 
            grepl("SARS-CoV-1", .data$viruses_false, fixed = TRUE)) &
           (grepl("SARS-CoV-2", .data$viruses_true, fixed = TRUE) | 
              grepl("SARS-CoV-2", .data$viruses_false, fixed = TRUE)))

cov1 <- grepl("SARS-CoV-1", both_tested$viruses_true, fixed = TRUE)
cov2 <- grepl("SARS-CoV-2", both_tested$viruses_true, fixed = TRUE)

sprintf("%i species susceptible to both SARS-CoV-1 and SARS-CoV-2, out of %i tested for both", 
        sum(cov2 & cov1),
        nrow(both_tested))

cat("\n\n")

# ---- Shedding data ------------------------------------------------------------------------------
# RNA
rna_shedding <- raw_data %>% 
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
