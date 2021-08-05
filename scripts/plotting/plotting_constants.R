# Constants used across all plots

PLOT_THEME <- theme_bw(base_size = 7) +
  theme(strip.text = element_text(margin = margin(2.5, 2.5, 2.5, 2.5)),
        legend.key.height = unit(0.6, 'lines'),
        legend.key.width = unit(0.6, 'lines'),
        legend.margin = margin(t = 2.5, r = 5.5, b = 2.5, l = 0),
        panel.grid = element_blank(),
        plot.background = element_blank())

theme_set(PLOT_THEME)


# Symbols
# - These can be combined to illustrate multiple viruses in a given host
VIRUS_SHAPES <- c("SARS-CoV" = 24, "SARS-CoV-2" = 25, "Other sarbecovirus" = 21)
