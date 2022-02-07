## Constants used across all plots

# Site info
source("scripts/utils/site_info.R")

# Theme
PLOT_THEME <- theme_bw(base_size = 7) +
  theme(strip.text = element_text(margin = margin(2.5, 2.5, 2.5, 2.5)),
        legend.key.height = unit(0.6, 'lines'),
        legend.key.width = unit(0.6, 'lines'),
        legend.margin = margin(t = 2.5, r = 5.5, b = 2.5, l = 0),
        panel.grid = element_blank())
# plot.background = element_blank()

theme_set(PLOT_THEME)


# Symbols
# - These can be combined to illustrate multiple viruses in a given host
VIRUS_SHAPES <- c("SARS-CoV" = 24, "SARS-CoV-2" = 25, "Other sarbecovirus" = 21)


# Colours
TRENDLINE_COLOUR <- "#009988"
INFECTION_STATUS_COLOURS <- c("True" = "#EE7733", "False" = "#0077BB")  # Paul Tol's "vibrant" palette
MISSING_DATA_COLOUR <- "#BBBBBB"


EVIDENCE_LABELS <- c("1" = "Observed infection",
                     "2" = "Experimental infection",
                     "3" = "Cell culture",
                     "4" = "Cell culture (het-ACE2)")
