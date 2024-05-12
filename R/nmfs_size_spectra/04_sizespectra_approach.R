####  NEFSC Trawl Survey Size Spectrum  ####
# Quick and dirty NEFSC trawl data run for quick display with app
# update to use wigley growth coefficients 10/27/2020
# Options to use abundances from survey directly, 
# or their aadjustments to total stratum areas


# Output files:
### size spectrum prepared data is exported to data/NEFSC/nefsc_ss_bins.csv


####_____________________####
####  Packages  ####
library(targets)
library(here)
library(gmRi)
library(janitor)
library(sizeSpectra)
library(patchwork)
library(tidyverse)

####  Support Functions  ####
source(here("R/support/sizeSpectra_support.R"))
source(here("R/support/nefsc_ss_build_nodrop.R"))

# data paths
nmfs_path  <- shared.path("unix", "RES_Data", "NMFS_Trawl")

####  Data  ####

# Use targets to load data
tar_load(nefsc_stratified)

# stop pretending that you prefer est_year to year
nefsc_stratified <- nefsc_stratified %>% 
  mutate(year = factor(est_year)) 




####_____________________####
####  Exploratory Plots  ####

# Exploratory plot of individual body mass by species
nefsc_stratified %>% 
  ggplot( aes( y = fct_reorder(comname, ind_weight_kg, .fun = mean, .desc = FALSE), 
               x = ind_weight_kg)) + 
  geom_boxplot() + 
  labs(x = "Individual Bodymass (kg)", y = "Common Name")


# Plot of size ranges by year, weighted by numlen_adj
nefsc_stratified %>% 
  ggplot(aes(ind_weight_kg, fct_rev(year), weight = numlen_adj)) +
  geom_boxplot(outlier.shape = NA) +
  scale_x_continuous(limits = quantile(nefsc_stratified$ind_weight_kg, c(0, 0.99))) +
  labs(x = "Individual Bodymass (kg)", y = "")




####  SizeSpectra Setup  ####

# Starting point for sizeSpectra steps
dataOrig <- nefsc_stratified %>% filter(is.na(ind_weight_kg) == FALSE)


# Keep desired columns, name them for the vignette
# NOTE: names do not reflect NEFSC survey design: CPUE_bio_per_hour not biomass per hour


data <- dataOrig %>%  
  as_tibble() %>% 
  mutate(season = str_to_title(season),
         season = factor(season, levels = c("Spring", "Fall"))) %>%
  rename(
    Year       = est_year,                     # year
    samp_month = est_month,
    samp_dat   = est_day,
    vessel     = svvessel,
    SpecCode   = comname,              # common name
    LngtClass  = length,              # length bin
    Number     = numlen_adj,             # adjusted number at length
    LWa        = a,                         # length/weight param
    LWb        = b) %>%                     # length/weight param
  mutate(
    bodyMass           = ind_weight_kg * 1000,      # weight in grams of individual
    Biomass            = Number * bodyMass,         # number at size times the mass
    Stratified_Biomass = expanded_abund_s * bodyMass,    # number scaled to stratum, time the mass of ind.
    lw_group           = str_c(SpecCode, season, catchsex)) %>% 
  arrange(Year, season, SpecCode, LngtClass)







####__  Calculate upper bin weights  ####
# So for each measurement increment there is a minimum 
# and maximum weight for what a fish could weigh

# Pretty sure since the growth coefficients,
# stratified numbers, and number at length are all together we can 
# just use mutate here...
dataBin <- data %>% 
  mutate(
    LngtMax         = LngtClass + 1,
    Ln_wmax         = (ln_a + LWb * log(LngtMax)),
    wmin            = bodyMass,                     # minimum weight of individual
    wmax            = exp(Ln_wmax) * 1000,          # max weight of individual
    wmin_sum        = Number * wmin,                # wmin * number caught actual
    wmax_sum        = Number * wmax,                # wmax * number caught actual  
    wmin_area_strat = expanded_abund_s * wmin,           # wmin * stratified abundance
    wmax_area_strat = expanded_abund_s * wmax            # wmax * stratified abundance
  )





#### __ Export NEFSC dataBin  ####
# write_csv(dataBin, here("data/NEFSC/nefsc_ss_bins.csv"))
tar_load(nefsc_databin)
glimpse(nefsc_databin )

# name change
dataBin <- nefsc_databin




####_____________________####

####  Set Bodymass Cutoff and Groups  ####

# Set bodymass lower limit
# subjective
mass_cutoff <- 0 #grams

# Filter for lower end of gear selectivity
dbin_trunc <- filter(dataBin, wmin >= mass_cutoff)

# Create Grouping Var
dbin_trunc <- dbin_trunc %>% 
  mutate(group_var = str_c(Year, season, sep = "_"))


####  1. Survey Abundance Curves  ####

# Get MLE parameters for each group
mle_bins <- dbin_trunc %>% 
  split(.$group_var) %>% 
  imap_dfr(group_mle_calc)  


# Make formatting changes for labels
mle_seasons <- mle_bins %>% 
  mutate(Year = str_sub(group_var, 1, 4),
         Year = as.numeric(as.character(Year)),
         season = str_sub(group_var, 6, -1),
         season = factor(season, levels = c("Spring", "Fall")),
         stdErr = (abs(confMin - b) + abs(confMax - b)) / (2 * 1.96),
         C = (b != -1 ) * (b + 1) / ( xmax^(b + 1) - xmin^(b + 1) ) + (b == -1) * 1 / ( log(xmax) - log(xmin))
         )


# Plot the slope parameter
mle_seasons %>% 
  filter(season %in% c("Spring", "Fall")) %>% 
  ggplot(aes(Year, b, color = season)) +
  geom_line(aes(group = season)) +
  geom_segment(aes(x = Year, xend = Year, y = confMin, yend = confMax)) +
  geom_point(aes(y = b)) + 
  geom_smooth(formula = y ~ x,
              method = "lm") +
  labs(x = NULL,
       y = "Size Spectrum Slope (b)",
       caption = paste0("Minimum Bodymass Cutoff: ", mass_cutoff)) +
  facet_wrap(~season, ncol = 1)



#####  Annual SS Curves  ####
# using ggplot

# change to dumb name used in package
MLEbins.res <- mle_seasons

# Data inputs

# 1. Raw data, prepped for plotting size bins
dataRecommend.isd <- dbin_trunc %>% 
  select(group_var, Year, season, wmin, wmax, Number) # Recommended data for ISD plots

# Take each group and generate the plots
data.year.list <- dataRecommend.isd %>% 
  split(.$group_var) %>% 
  map(isd_plot_prep)


# Global limits for plot
xlim.global <- c( min(dataRecommend.isd$wmin), max(dataRecommend.isd$wmax) )

# Vector of years, for naming
group_names <- sort(unique(dbin_trunc$group_var))


# 2. MLE coefficients that match them
MLEbins.res.list  <- MLEbins.res %>% split(.$group_var)    

# Loop through years for plots
seasonal_isd <- map2(data.year.list, MLEbins.res.list, .f = ggplot_isd) %>% 
  setNames(group_names)



# Take a peak
seasonal_isd$`2005_Spring`$stacked | seasonal_isd$`2005_Fall`$stacked
seasonal_isd$`1975_Spring`$stacked | seasonal_isd$`1975_Fall`$stacked





#### 2.  Stratified Abundance Curves  ####


# Prep data using stratified abundance:
mle_stratified <- dbin_trunc %>% 
  filter(expanded_abund_s != 0) %>% 
  split(.$group_var) %>% 
  imap_dfr(strat_abund_mle_calc)  


mle_strat_abund <- mle_stratified %>% 
  mutate(Year = str_sub(group_var, 1, 4),
         Year = as.numeric(as.character(Year)),
         season = str_sub(group_var, 6, -1),
         season = factor(season, levels = c("Spring", "Fall")))


# Plot the slope parameter
mle_strat_abund %>% 
  filter(season %in% c("Spring", "Fall")) %>% 
  ggplot(aes(Year, b, color = season)) +
  geom_line(aes(group = season)) +
  geom_segment(aes(x = Year, xend = Year, y = confMin, yend = confMax)) +
  geom_point(aes(y = b)) + 
  geom_smooth(formula = y ~ x,
              method = "lm") +
  labs(x = NULL,
       y = "Size Spectrum Slope (b)",
       caption = paste0("Minimum Bodymass Cutoff: ", mass_cutoff),
       subtitle = "Stratified Biomass") +
  facet_wrap(~season, ncol = 1)




##### Annual SS Curves  ####
# Data inputs:

# Raw data to plot
strat_isd_data  <- dbin_trunc %>% 
  select(group_var, Year, season, vessel, wmin, wmax, wmin_area_strat, wmax_area_strat, expanded_abund_s) %>% 
  filter(expanded_abund_s > 0)

# Prepped for plotting the size bins
strat_year_list <- strat_isd_data %>% 
  split(.$group_var) %>% 
  map(strat_isd_prep)

# MLE results for fit, split into groups to match the data
strat_res_list  <- mle_strat_abund %>% 
  split(.$group_var) 


# Global limits for plots to compare across seasons
xlim.global <- c( min(strat_isd_data$wmin), 
                  max(strat_isd_data$wmax) )

# Vector of years, for naming
group_names <- sort(unique(dbin_trunc$group_var))

# Loop through years for plots
seasonal_strat_isd <- map2(strat_year_list, strat_res_list, .f = ggplot_strat_isd) %>% 
  setNames(group_names)

# Take a peak
seasonal_strat_isd$`2005_Spring`$obs_y
seasonal_strat_isd$`2005_Spring`$stacked | seasonal_strat_isd$`2005_Fall`$stacked
