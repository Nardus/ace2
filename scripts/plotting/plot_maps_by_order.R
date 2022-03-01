## Plot range maps for terrestrial wildlife predictions, separated by taxonomic order
#     - Supplement to "plot_holdout_maps.R"


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
phylogeny_predictions <- read_rds("output/all_data/infection/phylogeny/holdout_predictions.rds")

infection_data <- read_rds("data/calculated/cleaned_infection_data.rds")
taxonomy <- read_rds("data/calculated/taxonomy.rds")

continent_outlines <- st_read("data/external/iucn_base/", "Land_Masses_and_Ocean_Islands")
iucn_ranges <- st_read("data/iucn_range_maps/", "MAMMALS_TERRESTRIAL_ONLY")


# ---- Fix taxonomy --------------------------------------------------------------------------------
# Only plotting land-based, wild mammals:
domestic_species <- c("Bos taurus", "Felis catus", "Equus caballus", "Cavia porcellus",
                      "Cricetulus griseus", "Ovis aries", "Capra hircus", "Bubalus bubalis",
                      "Vicugna pacos", "Equus przewalskii", "Camelus bactrianus",
                      "Equus asinus", "Bos indicus", "Bos indicus x Bos taurus",
                      "Camelus dromedarius", "Canis familiaris")

plot_species <- data.frame(species = phylogeny_predictions$species) %>% 
  distinct() %>% 
  filter(.data$species != "Homo sapiens") %>% 
  filter(!.data$species %in% domestic_species) %>% 
  filter(.data$species %in% iucn_ranges$binomial) %>% 
  pull(.data$species)

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

# Display no data as grey (on land)
base_plot_prop <- base_plot +
  geom_sf(fill = "grey20", size = 0.05, data = continent_outlines) +
  coord_sf(ylim = c(-62.5, 90), expand = FALSE)


# ---- Proportion of each order susceptible --------------------------------------------------------
get_order_raster <- function(predictions, taxonomic_order, base_raster = blank_raster, range_data = iucn_ranges) {
  # Get a raster counting the proportion of species from a given order predicted as susceptible 
  # in each cell
  
  # Count the number of species from this order in each cell
  order_preds <- predictions %>% 
    filter(.data$order == taxonomic_order)
  
  order_ranges <- range_data %>% 
    filter(.data$binomial %in% order_preds$species)
    
  all_species_raster <- fasterize(order_ranges, base_raster, fun = "sum")
  
  # Count susceptibles
  susceptible_species <- order_preds %>% 
    filter(.data$predicted_label == "True") %>%
    pull(.data$species)
  
  susceptible_ranges <- range_data %>% 
    filter(.data$binomial %in% susceptible_species)
  
  susceptible_raster <- fasterize(susceptible_ranges, base_raster, fun = "sum")
  
  # Return proportion
  susceptible_raster/all_species_raster
}

predictions <- phylogeny_predictions %>% 
    rename(internal_name = .data$species) %>% 
    left_join(taxonomy, by = "internal_name") %>% 
    filter(!is.na(.data$order))

# Only plot orders with at least one susceptible species
tax_orders <- predictions %>% 
    filter(.data$predicted_label == "True") %>% 
    distinct(.data$order) %>%
    arrange(.data$order) %>% 
    pull(.data$order)

rasters <- sapply(tax_orders, get_order_raster, predictions = predictions, 
                  simplify = FALSE, USE.NAMES = TRUE)


# ---- Plot ----------------------------------------------------------------------------------------
plot_order_raster <- function(raster_obj, taxonomic_order, label = "Proportion\nsusceptble", 
                              base = base_plot_prop, guide = FALSE, limits = c(0, 1)) {
  rasterdf <- as.data.frame(raster_obj, xy = TRUE)
  
  p <- base +  
    geom_raster(aes(x = x, y = y, fill = layer), data = rasterdf) +
    scale_fill_viridis_c(breaks = breaks_pretty(), direction = 1, na.value = NA, limits = limits) +
    labs(title = taxonomic_order, fill = label)
  
  if (!guide) {
    p <- p + guides(fill = "none")
  }
  
  p
}

plots <- mapply(plot_order_raster,
                raster_obj = rasters, 
                taxonomic_order = names(rasters),
                SIMPLIFY = FALSE, 
                USE.NAMES = TRUE)

# Combine
legend_plot <- plot_order_raster(rasters$Chiroptera, "Chiroptera", guide = TRUE)
legend <- get_legend(legend_plot)

plots <- c(plots, list(legend))

final_plot <- plot_grid(plotlist = plots, 
                        ncol = 3,
                        align = "v", axis = "lrtb")

ggsave2("output/plots/prediction_maps_by_order.png", final_plot, width = 6, height = 6)
