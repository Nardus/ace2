# Plot range maps for terrestrial wildlife predictions

suppressPackageStartupMessages({
  library(ape)
  library(sf)
  library(raster)
  library(fasterize)
  library(dplyr)
  library(tibble)
  library(readr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(colorspace)
  library(cowplot)
  library(parallel)
  
  source("scripts/utils/plot_utils.R")
  source("scripts/plotting/plotting_constants.R")
})

RESOLUTION <- 1/6


# ---- Data ----------------------------------------------------------------------------------------
ensemble_predictions <- read_rds("output/all_data/infection/ensemble_all_features_phylogeny/holdout_predictions.rds")
phylogeny_predictions <- read_rds("output/all_data/infection/phylogeny/holdout_predictions.rds")

infection_data <- read_rds("data/calculated/cleaned_infection_data.rds")
taxonomy <- read_rds("data/calculated/taxonomy.rds")

mammal_tree <- read.tree("data/internal/timetree_mammalia.nwk")
bird_tree <- read.tree("data/internal/timetree_aves.nwk")

continent_outlines <- st_read("data/external/iucn_base/", "Land_Masses_and_Ocean_Islands")
iucn_range_terrestrial <- st_read("data/iucn_range_maps/", "MAMMALS_TERRESTRIAL_ONLY")
iucn_range_freshwater <- st_read("data/iucn_range_maps/", "MAMMALS_FRESHWATER") # Semi-aquatic
iucn_range_marineter <- st_read("data/iucn_range_maps/", "MAMMALS_MARINE_AND_TERRESTRIAL") # Semi-aquatic (marine)

iucn_ranges <- iucn_range_terrestrial
excluded_ranges <- rbind(iucn_range_freshwater, iucn_range_marineter)


# ---- Fix taxonomy --------------------------------------------------------------------------------
# Only plotting land-based, wild mammals:
domestic_species <- c("Bos taurus", "Felis catus", "Equus caballus", "Cavia porcellus",
                      "Cricetulus griseus", "Ovis aries", "Capra hircus", "Bubalus bubalis",
                      "Vicugna pacos", "Equus przewalskii", "Camelus bactrianus",
                      "Equus asinus", "Bos indicus", "Bos indicus x Bos taurus",
                      "Camelus dromedarius", "Canis familiaris")

mammals <- mammal_tree$tip.label %>% 
  str_replace("_", " ")

plot_species <- data.frame(species = c(ensemble_predictions$species, 
                                       phylogeny_predictions$species)) %>% 
  distinct() %>% 
  filter(.data$species %in% mammals) %>% 
  filter(.data$species != "Homo sapiens") %>% 
  filter(!.data$species %in% domestic_species) %>% 
  filter(!.data$species %in% excluded_ranges$binomial) %>% 
  pull(.data$species)


missing_spp <- plot_species[!plot_species %in% iucn_ranges$binomial]

if(length(missing_spp) != 0)
  warning("Some species do not have range data: ", paste(missing_spp, collapse = ", "))

ensemble_predictions <- ensemble_predictions %>% 
  filter(.data$species %in% plot_species)

phylogeny_predictions <- phylogeny_predictions %>% 
  filter(.data$species %in% plot_species)


# ---- Base plot -----------------------------------------------------------------------------------
blank_raster <- raster(resolution = RESOLUTION, crs = st_crs(iucn_ranges))  

# Basic plot settings always the same:
base_plot <- ggplot() +
  labs(x = NULL, y = NULL) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 7))

# When counting, no data in raster means 0 (when on land)
base_plot_counts <- base_plot +
  geom_sf(aes(fill = 0), size = 0.05, data = continent_outlines) +
  coord_sf(ylim = c(-62.5, 90), expand = FALSE)

# When calculating frequencies, all range data are retained in the denominator, so
# areas with no species are missing data
base_plot_freq <- base_plot +
  geom_sf(fill = "white", size = 0.05, data = continent_outlines) +
  coord_sf(ylim = c(-62.5, 90), expand = FALSE)


# ---- Total susceptible -----------------------------------------------------------------------------
get_predicted_raster <- function(predictions, base_raster = blank_raster, range_data = iucn_ranges) {
  # Get a raster counting the number of species predicted as susceptible in each cell
  susceptible_ranges <- range_data %>% 
    filter(.data$binomial %in% predictions$species[predictions$predicted_label == "True"])
  
  fasterize(susceptible_ranges, base_raster, fun = "sum")
}

plot_count_raster <- function(raster_obj, label = "Susceptible\nspecies", base = base_plot_counts, 
                              guide = TRUE, limits = NULL) {
  rasterdf <- as.data.frame(raster_obj, xy = TRUE)
  
  p <- base +  
    geom_raster(aes(x = x, y = y, fill = layer), data = rasterdf) +
    scale_fill_viridis_c(breaks = breaks_pretty(), direction = 1, na.value = NA, limits = limits) +
    labs(fill = label)
  
  if (!guide) {
    p <- p + guides(fill = "none")
  }
  
  p
}

susceptible_raster_ensemble <- get_predicted_raster(ensemble_predictions)
susceptible_raster_phylogeny <- get_predicted_raster(phylogeny_predictions)

p_count_ensemble <- plot_count_raster(susceptible_raster_ensemble) +
  ggtitle("ACE2 / phylogeny ensemble")

p_count_phylogeny <- plot_count_raster(susceptible_raster_phylogeny) +
  ggtitle("Phylogeny-only")


# ---- Observed / expected -------------------------------------------------------------------------
get_obs_exp_raster <- function(predictions, susceptible_raster, susceptible_frequency,
                             base_raster = blank_raster, range_data = iucn_ranges) {
  # Expected
  expected_ranges <- range_data %>% 
    filter(.data$binomial %in% predictions$species)
  
  species_count <- fasterize(expected_ranges, base_raster, fun = "sum") # Number of species available in each raster cell
  
  species_present <- fasterize(expected_ranges, base_raster, fun = "any") # 1 if at least one species with data present in a given raster cell
  expected_frequency <- species_present * susceptible_frequency
  
  # Observed
  observed_frequency <- susceptible_raster / species_count
  
  # Proportion
  observed_frequency / expected_frequency
}

plot_oe_raster <- function(raster_obj, label = "Observed/\nexpected", base = base_plot_freq, 
                              guide = TRUE, limits = NULL) {
  rasterdf <- as.data.frame(raster_obj, xy = TRUE)
  
  p <- base +  
    geom_raster(aes(x = x, y = y, fill = layer), data = rasterdf) +
    scale_fill_continuous_diverging("Blue-Red 2", breaks = breaks_pretty(), na.value = NA, 
                                    limits = limits, mid = 1) +
    labs(fill = label)
  
  if (!guide) {
    p <- p + guides(fill = "none")
  }
  
  p
}

# Expected frequency: based on observations among mammals
mammal_data <- infection_data %>% 
  rename(internal_name = .data$species) %>% 
  left_join(taxonomy, by = "internal_name") %>% 
  filter(.data$class == "Mammalia")

base_frequency <- sum(mammal_data$infected == "True")/nrow(mammal_data)

oe_ensemble <- get_obs_exp_raster(ensemble_predictions, susceptible_raster_ensemble, 
                                  susceptible_frequency = base_frequency)

oe_phylogeny <- get_obs_exp_raster(phylogeny_predictions, susceptible_raster_phylogeny,
                                   susceptible_frequency = base_frequency)


# Plot - using same legend scale for both
shared_min_oe <- min(c(oe_ensemble@data@values,
                       oe_phylogeny@data@values),
                     na.rm = TRUE)

shared_max_oe <- max(c(oe_ensemble@data@values,
                       oe_phylogeny@data@values),
                     na.rm = TRUE)

p_obs_ensemble <- plot_oe_raster(oe_ensemble, 
                                 limits = c(shared_min_oe, shared_max_oe))

p_obs_phylogeny <- plot_oe_raster(oe_phylogeny, 
                                  limits = c(shared_min_oe, shared_max_oe))


# ---- Output --------------------------------------------------------------------------------------
combined_plot <- plot_grid(p_count_ensemble, p_count_phylogeny,
                           p_obs_ensemble, p_obs_phylogeny,
                           nrow = 2, rel_heights = c(1.11, 1),
                           align = "v", axis = "lrtb",
                           labels = c("A", "B", "C", "D"),
                           vjust = c(2.5, 2.5, 1.5, 1.5),
                           label_colour = c("grey20", "white", "grey20", "grey20"))

ggsave2("output/plots/prediction_maps.png", combined_plot,
        width = 7, height = 2.74)


# ---- Values in text ------------------------------------------------------------------------------
birds <- bird_tree$tip.label %>% 
  str_replace("_", " ")

# Observed data
bird_data <- infection_data %>% 
  rename(internal_name = .data$species) %>% 
  left_join(taxonomy, by = "internal_name") %>% 
  filter(.data$class == "Aves")

cat(sprintf("\n%3.3f%% of mammals amd %3.3f%% of birds in current data are susceptible (N =%3i, %3i)\n", 
            sum(mammal_data$infected == "True") / nrow(mammal_data) * 100,
            sum(bird_data$infected == "True") / nrow(bird_data) * 100,
            nrow(mammal_data),
            nrow(bird_data)))

# Predictions
ensemble_predictions <- read_rds("output/all_data/infection/ensemble_all_features_phylogeny/holdout_predictions.rds")
phylogeny_predictions <- read_rds("output/all_data/infection/phylogeny/holdout_predictions.rds")

get_stats <- function(predictions, class_label) {
  pred_positive <- sum(predictions$predicted_label == "True" & predictions$class == class_label)
  n <- sum(predictions$class == class_label)
  
  list(
    freq = pred_positive/n,
    n = n
  )
}

ensemble <- ensemble_predictions %>% 
  rename(internal_name = .data$species) %>% 
  left_join(taxonomy, by = "internal_name") %>% 
  filter(.data$class %in% c("Mammalia", "Aves"))

phylo <- phylogeny_predictions %>% 
  filter(.data$species %in% c(mammals, birds)) %>% 
  mutate(class = if_else(.data$species %in% mammals, "Mammalia", "Aves"))


ensemble_mammal <- get_stats(ensemble, "Mammalia")
ensemble_bird <- get_stats(ensemble, "Aves")

cat(sprintf("Ensemble model predicts %3.3f%% of mammals and %3.3f%% of birds to be susceptible (N = %3i, %3i)\n",
            ensemble_mammal$freq * 100, 
            ensemble_bird$freq * 100,
            ensemble_mammal$n,
            ensemble_bird$n))

phylo_mammal <- get_stats(phylo, "Mammalia")
phylo_bird <- get_stats(phylo, "Aves")

cat(sprintf("Phylogeny model predicts %3.3f%% of mammals and %3.3f%% of birds to be susceptible (N = %3i, %3i)\n",
            phylo_mammal$freq * 100, 
            phylo_bird$freq * 100,
            phylo_mammal$n,
            phylo_bird$n))


# ---- Susceptible species in hotspots -------------------------------------------------------------
get_species_max <- function(species, ranges, raster) {
  # Find the maximum value in 'raster' within a given species' range
  sp_range <- ranges %>% 
    filter(.data$binomial == species)
  
  sp_raster <- mask(raster, sp_range)
  
  maxValue(sp_raster)
}

find_hotspot_species <- function(oe_raster, predictions, cutoff, all_ranges = iucn_ranges) {
  valid_ranges <- all_ranges %>% 
    filter(.data$binomial %in% predictions$species)
  
  all_species <- unique(valid_ranges$binomial)
  
  range_max <- mclapply(all_species, get_species_max,
                        ranges = valid_ranges,
                        raster = oe_raster,
                        mc.cores = 8)
  
  range_max <- unlist(range_max)
  range_max[is.na(range_max)] <- -999  # NA's caused by very small ranges (below resolution of raster, but we only want the main species visible on the map)
  
  all_species[range_max > cutoff]
}

hotspot_spp_ensemble <- find_hotspot_species(oe_ensemble, ensemble_predictions, cutoff = 1.1)
hotspot_spp_phylogeny <- find_hotspot_species(oe_phylogeny, phylogeny_predictions, cutoff = 1.1)

hotspot_susceptibles_ensemble <- ensemble_predictions %>% 
  filter(.data$predicted_label == "True") %>% 
  filter(.data$species %in% hotspot_spp_ensemble) %>% 
  pull(.data$species)

hotspot_susceptibles_phylogeny <- phylogeny_predictions %>% 
  filter(.data$predicted_label == "True") %>% 
  filter(.data$species %in% hotspot_spp_phylogeny) %>% 
  pull(.data$species)

# Plot these (for internal use)
plot_ranges <- function(all_species, susceptibles, out_folder,
                        all_ranges = iucn_ranges, base_plot = base_plot_freq) {
  for (spp in all_species) {
    rnge <- all_ranges %>% 
      filter(.data$binomial == spp)
    
    label <- if_else(spp %in% susceptibles,
                     paste(spp, "(susceptible)"),
                     paste(spp, "(not susceptible)"))
    
    p <- base_plot + 
      geom_sf(colour = "red", data = rnge) + 
      ggtitle(label)
    
    out_name <- paste(spp, ".png") %>% 
      paste(out_folder, ., sep = "/")
      
    ggsave2(out_name, p, width = 7, height = 3.2)
  }
}

dir.create("output/plots/hotspot_species/ensemble", recursive = TRUE)
dir.create("output/plots/hotspot_species/phylogeny")

plot_ranges(hotspot_spp_ensemble, hotspot_susceptibles_ensemble,
            out_folder = "output/plots/hotspot_species/ensemble")


plot_ranges(hotspot_spp_phylogeny, hotspot_susceptibles_phylogeny,
            out_folder = "output/plots/hotspot_species/phylogeny")
