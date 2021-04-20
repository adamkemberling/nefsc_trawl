##### 
# Size Spectrum App Tools  ####
# Goal: Create functions to quickly estimate and plot 
# size spectrum fits for varying time windows and different groups
# for the maine and nh trawl survey as well as the NEFSC data.

#### Code Status: 12/30/2020
# Needs to be reviewed since changes upstream to data builds occurred


####  Packages  ####

library(here)
library(gmRi)
library(janitor)
library(sizeSpectra)
library(patchwork)
library(tidyverse)

####  Support Functions  ####
source(here("R/support/sizeSpectra_support.R"))



####  Data  ####

# Load the Catch of each station, with the length and weight bins for each species
# there has been no biomass cutoff at this stage
nefsc_dbin <- read_csv(here("data/NEFSC/nefsc_ss_bins.csv"), col_types = cols())
menh_dbin  <- read_csv(here("data/MENH/menh_ss_bins.csv"), col_types = cols())






####  Biomass Cutoff Controls  ####

# Set the mass cutoff in grams
mass_cutoff <- 20

# Filter biomass by cutoff
nefsc_dbin <- filter(nefsc_dbin, wmin >= mass_cutoff)
menh_dbin <- filter(menh_dbin, wmin >= mass_cutoff)



####  Grouping Variable Options  ####

# Create Grouping Var
menh_dbin <- menh_dbin %>% mutate(group_var = str_c(Year, season, sep = "_"))
nefsc_dbin <- nefsc_dbin %>% mutate(group_var = str_c(Year, season, sep = "_"))



####  Time Window  ####

# So here you probably wan't to just duplicate the dataset
# Set A will have a custom year-range, default the full range
# Set B will have the grouping variables


####  Slope Estimation  ####

# Using NEFSC as tester
dataBin <- nefsc_dbin

# Map through instead of looping
grouped_size_spectra <- dataBin %>% 
  split(.$group_var) %>% 
  imap_dfr(group_mle_calc) %>% 
  mutate(
    stdErr = (abs(confMin - b) + abs(confMax - b)) / (2 * 1.96),
    Year = str_sub(group_var, 1, 4),
    Year = as.numeric(as.character(Year)),
    season = str_sub(group_var, 6, -1),
    season = factor(season, levels = c("SPRING", "SUMMER", "FALL", "WINTER")),
    C = (b != -1 ) * (b + 1) / ( xmax^(b + 1) - xmin^(b + 1) ) + (b == -1) * 1 / ( log(xmax) - log(xmin)))





####  Plotting Annual / Group Slopes  ####
# Better plot
grouped_size_spectra %>% 
  filter(season %in% c("SPRING", "FALL")) %>% 
  ggplot(aes(Year, b, color = season)) +
  geom_line(aes(group = season)) +
  geom_segment(aes(x = Year, xend = Year, y = confMin, yend = confMax)) +
  geom_point(aes(y = b)) + 
  geom_smooth(formula = y ~ x,
              method = "lm") +
  labs(x = NULL,
       y = "Size Spectrum Slope (b)") +
  facet_wrap(~season, ncol = 1)




####  Plotting Curves  ####

# Data inputs
dbin_group_list <- dataBin %>% 
  split(.$group_var) %>% 
  map(isd_plot_prep)

# global limits for plot
xlim.global <- c( min(dataBin$wmin), max(dataBin$wmax) )

# Vector of years, for naming
group_names <- names(dbin_group_list)
names(group_names) <- names(dbin_group_list)

# Size Spectra Coefficients for each group
ss_results_list  <- grouped_size_spectra %>% split(.$group_var)      

# Verify names match
group_names == names(ss_results_list)

# Map through the groups for each plot
seasonal_isd <- map2(dbin_group_list, ss_results_list, .f = ggplot_isd) %>% setNames(group_names)

# Take a peak
seasonal_isd[["2005_SPRING"]]["stacked"]

# Put them in a list by year

