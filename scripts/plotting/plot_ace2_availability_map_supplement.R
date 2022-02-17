## Plot host ranges of all mammals with available ACE2 sequences
#   - Supplement to "plot_holdout_maps.R"

suppressPackageStartupMessages({
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
  
  source("scripts/utils/plot_utils.R")
  source("scripts/plotting/plotting_constants.R")
})

RESOLUTION <- 1/6


# ---- Data ----------------------------------------------------------------------------------------
ensemble_predictions <- read_rds("output/all_data/infection/ensemble_all_features_phylogeny/holdout_predictions.rds")

taxonomy <- read_rds("data/calculated/taxonomy.rds")

continent_outlines <- st_read("data/external/iucn_base/", "Land_Masses_and_Ocean_Islands")
iucn_range_terrestrial <- st_read("data/iucn_range_maps/", "MAMMALS_TERRESTRIAL_ONLY")
iucn_range_freshwater <- st_read("data/iucn_range_maps/", "MAMMALS_FRESHWATER") # Semi-aquatic

iucn_ranges <- rbind(iucn_range_terrestrial, iucn_range_freshwater)


# ---- Fix taxonomy --------------------------------------------------------------------------------
# Only plotting land-based, wild mammals:
domestic_species <- c("Bos taurus", "Felis catus", "Equus caballus", "Cavia porcellus",
                      "Cricetulus griseus", "Ovis aries", "Capra hircus", "Bubalus bubalis",
                      "Vicugna pacos", "Equus przewalskii", "Camelus bactrianus",
                      "Equus asinus", "Bos indicus", "Bos indicus x Bos taurus",
                      "Camelus dromedarius", "Canis familiaris")
marine_species <- c("Orcinus orca", "Tursiops truncatus", 
                    "Physeter catodon", "Balaenoptera acutorostrata",
                    "Delphinapterus leucas", "Lagenorhynchus obliquidens", "Monodon monoceros",
                    "Globicephala melas", "Phocoena sinus", "Balaenoptera musculus")
semimarine_species <- c("Odobenus rosmarus", "Leptonychotes weddellii", "Ursus maritimus", 
                        "Neomonachus schauinslandi", "Enhydra lutris", "Callorhinus ursinus",
                        "Zalophus californianus", "Eumetopias jubatus", "Mirounga leonina",
                        "Halichoerus grypus") 

mammals <- taxonomy %>% 
  filter(.data$class == "Mammalia") %>% 
  pull(.data$species)

plot_species <- data.frame(species = ensemble_predictions$species) %>% 
  distinct() %>% 
  filter(.data$species %in% mammals) %>% 
  filter(.data$species != "Homo sapiens") %>% 
  filter(!.data$species %in% domestic_species) %>% 
  filter(!.data$species %in% marine_species) %>% 
  filter(!.data$species %in% semimarine_species) %>% 
  pull(.data$species)

missing_spp <- plot_species[!plot_species %in% iucn_ranges$binomial]

if(length(missing_spp) != 0)
  warning("Some species do not have range data: ", paste(missing_spp, collapse = ", "))


# ---- Base plot -----------------------------------------------------------------------------------
blank_raster <- raster(resolution = RESOLUTION, crs = st_crs(iucn_ranges))  

# When counting, no data in raster means 0 (when on land)
base_plot_counts <- ggplot() +
  geom_sf(aes(fill = 0), size = 0.05, data = continent_outlines) +
  coord_sf(expand = FALSE) +
  scale_x_continuous(breaks = seq(-180, 180, by = 20)) +
  scale_y_continuous(limits = c(-62.5, 90),
                     breaks = seq(-60, 80, by = 20)) +
  labs(x = NULL, y = NULL)


# ---- Sequences available -------------------------------------------------------------------------
available_ranges <- iucn_ranges %>% 
  filter(.data$binomial %in% plot_species)

available_raster <- fasterize(available_ranges, blank_raster, fun = "sum")
rasterdf <- as.data.frame(available_raster, xy = TRUE)

p <- base_plot_counts +  
  geom_raster(aes(x = x, y = y, fill = layer), data = rasterdf) +
  scale_fill_viridis_c(breaks = breaks_pretty(), direction = 1, na.value = NA) +
  labs(fill = "Species with\nACE2 sequences")


# ---- Save ----------------------------------------------------------------------------------------
ggsave2("output/plots/ace2_availability_map_supplement.png", p,
        width = 7, height = 3.2)