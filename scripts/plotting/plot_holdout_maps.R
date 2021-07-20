# Plot range maps for terrestrial wildlife predictions

suppressPackageStartupMessages({
  library(sf)
  library(raster)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(taxize)
  library(ggplot2)
  library(scales)
  library(cowplot)
  library(parallel)
})

source("scripts/utils/plot_utils.R")


# ---- Data ----------------------------------------------------------------------------------------
holdout_metadata <- read_csv("data/internal/NCBI_ACE2_orthologs.csv")

continent_outlines <- st_read("data/external/iucn_base/", "Land_Masses_and_Ocean_Islands")
iucn_range_terrestrial <- st_read("data/iucn_range_maps/", "MAMMALS_TERRESTRIAL_ONLY")
iucn_range_freshwater <- st_read("data/iucn_range_maps/", "MAMMALS_FRESHWATER") # Semi-aquatic
iucn_range_marineter <- st_read("data/iucn_range_maps/", "MAMMALS_MARINE_AND_TERRESTRIAL") # Semi-aquatic (marine)

iucn_ranges <- rbind(iucn_range_terrestrial, iucn_range_freshwater, iucn_range_marineter)


# ---- Fix taxonomy --------------------------------------------------------------------------------
taxonomy_table <- classification(unique(holdout_metadata$`Scientific name`), db = "ncbi") 

taxonomy_table <- taxonomy_table %>% 
  lapply(data.frame) %>% 
  bind_rows(.id = "lookup_species") %>% 
  select(-.data$id) %>% 
  filter(!.data$rank %in% c("no rank", "clade")) %>% 
  pivot_wider(names_from = .data$rank, values_from = .data$name)

# Only plotting land-based, wild mammals:
domestic_species <- c("Bos taurus", "Felis catus", "Equus caballus", "Cavia porcellus",
                      "Cricetulus griseus", "Ovis aries", "Capra hircus", "Bubalus bubalis",
                      "Vicugna pacos", "Equus przewalskii", "Camelus bactrianus",
                      "Equus asinus", "Bos indicus", "Bos indicus x Bos taurus",
                      "Camelus dromedarius")
marine_species <- c("Orcinus orca", "Tursiops truncatus", 
                    "Physeter catodon", "Balaenoptera acutorostrata",
                    "Delphinapterus leucas", "Lagenorhynchus obliquidens", "Monodon monoceros",
                    "Globicephala melas", "Phocoena sinus", "Balaenoptera musculus")
semimarine_species <- c("Odobenus rosmarus", "Leptonychotes weddellii", "Ursus maritimus", 
                        "Neomonachus schauinslandi", "Enhydra lutris", "Callorhinus ursinus",
                        "Zalophus californianus", "Eumetopias jubatus", "Mirounga leonina",
                        "Halichoerus grypus") # These not currently excluded (seals, etc.)

taxonomy_table <- taxonomy_table %>% 
  filter(class == "Mammalia") %>% 
  filter(species != "Homo sapiens") %>% 
  filter(!species %in% domestic_species) %>% 
  filter(!species %in% marine_species)

# Fix known issues / synonyms:
# - Based on homotypic synonyms listed in NCBI taxonomy
taxonomy_table <- taxonomy_table %>% 
  mutate(species = case_when(species == "Tupaia chinensis" ~ "Tupaia belangeri",
                             species == "Nannospalax galili" ~ "Nannospalax ehrenbergi",
                             species == "Grammomys surdaster" ~ "Grammomys dolichurus",
                             TRUE ~ species))

missing_spp <- taxonomy_table$species[!taxonomy_table$species %in% iucn_ranges$binomial]
stopifnot(length(missing_spp) == 0)


# ---- Species with ACE2 sequences available -------------------------------------------------------
all_holdout_ranges <- iucn_ranges %>% 
  filter(binomial %in% taxonomy_table$species)

blank_raster <- raster(resolution = 0.5, crs = st_crs(iucn_ranges))  
all_holdout_raster <- rasterize(all_holdout_ranges, blank_raster, 
                                field = "binomial", fun = n_distinct)

all_holdout_rasterdf <- as.data.frame(all_holdout_raster, xy = TRUE)

ggplot() +
  geom_sf(aes(fill = 0), size = 0.05, data = continent_outlines) +  # No data in raster means 0 (when on land)
  geom_raster(aes(x = x, y = y, fill = layer), data = all_holdout_rasterdf) +
  scale_fill_viridis_c(breaks = breaks_pretty(), limits = c(0, NA),
                       direction = 1, na.value = NA) +
  scale_y_continuous(limits = c(-62.5, 90),
                     breaks = seq(-60, 80, by = 20)) +
  coord_sf(expand = FALSE) +
  labs(x = NULL, y = NULL, fill = "Number\nof species") +
  theme_bw()


# ---- Predicted hosts -----------------------------------------------------------------------------

# TODO: left join predictions to iucn_ranges, then filter?


# ---- Predictions by host order -------------------------------------------------------------------
rasterize_by_order <- function(tax_order, 
                               all_ranges = all_holdout_ranges, 
                               target_raster = blank_raster) {
  order_ranges <- all_ranges %>% 
    filter(.data$order_ == str_to_upper(tax_order))
  
  filled_raster <- rasterize(order_ranges, target_raster, 
                             field = "binomial", fun = n_distinct)
  
  as.data.frame(filled_raster, xy = TRUE) %>% 
    mutate(order = tax_order)
}

MIN_SIZE <- 5

ranges_by_order <- all_holdout_ranges %>% 
  st_drop_geometry() %>% 
  group_by(.data$order_) %>% 
  summarise(n = n_distinct(.data$binomial), 
            .groups = "drop") %>% 
  filter(.data$n > MIN_SIZE) %>% 
  pull(.data$order_) %>% 
  mclapply(rasterize_by_order, mc.cores = 16) %>% 
  bind_rows() %>% 
  mutate(order = str_to_sentence(.data$order))


# TODO: filter predicted hosts only

ggplot() +
  geom_sf(aes(fill = 0), size = 0.05, data = continent_outlines) +  # No data in raster means 0 (when on land)
  geom_raster(aes(x = x, y = y, fill = layer), data = ranges_by_order) +
  facet_wrap(vars(order)) +
  scale_fill_viridis_c(breaks = breaks_pretty(), limits = c(0, NA),
                       direction = 1, na.value = NA) +
  scale_y_continuous(limits = c(-62.5, 90),
                     breaks = seq(-60, 80, by = 20)) +
  coord_sf(expand = FALSE) +
  labs(x = NULL, y = NULL, fill = "Number\nof species") +
  theme_bw()
