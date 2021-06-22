####
# Digging into size spectra results from nmfs trawl data
# data processed in groups to get group ss slopes


####  Packages  ####
library(targets)
library(here)
library(sf)
library(tidyverse)
library(gmRi)
library(patchwork)

# Support functions
source(here("R/support/sizeSpectra_support.R"))
source(here("R/support/temp_support.R"))

####  Data  ####
tar_load(nefsc_1g)                # Input data for nmfs ss, truncated
# tar_load(nmfs_group_mle_ss)       # results for survey abundance
# tar_load(nmfs_stratified_mle_ss)  # results for stratified abundance
# tar_load(nmfs_log10_slopes)       # Output from manual bins
tar_load(size_spectrum_indices)   # SS and manual bins together
tar_load(regional_oisst)          # OISST for all the regions in one


# quick format
size_spectrum_indices <- size_spectrum_indices %>% 
  mutate(survey_area = factor(survey_area, 
                              levels = c("GoM", "GB", "SNE", "MAB")),
         yr = as.numeric(as.character(Year)))

# # Code to run slopes for a custom group - decades
# decade_slopes <- nefsc_1g %>% 
#   mutate(decade = floor_decade(Year)) %>% 
#   group_ss_slope_estimate(min_weight_g = 50, 
#                           abundance_vals = "observed",
#                           .group_cols = c("Year", "decade"))






####  Hypotheses  ####

# Hypothesis 1
"
H1: Warming enhances seasonal stratification, leading to enhanced phytoplankton production in 
the winter, a more diverse zooplankton community, and changes in higher trophic levels...

...Rigorously comparing the decadal changes in the plankton and fish community and comparing the 
communities in the Gulf of Maine and mid-Atlantic is a major goal of our project...
"


# Hypothesis 2
"
H2: Warming alters the community through the direct influence of temperature on metabolism, 
growth, and population productivity...

...We expect that changes in growth rates, especially in immature fish, will be the first 
response in the community to increases in temperature and may be especially prominent during 
the heatwaves. After several years of warm conditions, the size of adults will decline...

...As warming occurs, predators like sharks may become more important to ecosystem dynamics in 
northern regions and change ecosystem structure and functioning (Sagarese et al. 2010). Thus, 
we expect that food-webs in the northwest Atlantic will shift from ones where a few species
such as Calanus finmarchicus, herring, and cod dominate the flow of biomass to higher trophic 
levels to food webs where there are more pathways and where biomass is spread more evenly...
"


# Hypothesis 3
"
H3: Warming causes marine communities to be more diverse leading to changes in trophic pathways 
and ecosystem stability, but abrupt temperature changes degrade ecosystem function...

...We expect to see a clear gradient between our focal regions, with ecosystem properties of 
Gulf of Maine gradually approaching those of the historic mid-Atlantic as temperatures increase. 
We also expect that a large disturbance like a marine heatwave will lead to quantifiable changes 
in how biomass flows through the food web, even if the biomass of most species remains relatively
unchanged. For example, the 2012 heatwave lengthened the period of summer-like conditions by more
than a month. This could tilt the flow of energy in the ecosystem toward predators that migrate 
into the region during the warm season (Pershing et al. 2015).
"


####  EDA  ####


#####  Decade Check  ####
# How do the different size spectrum groups play out across decades?

# Pull the group ID for the slopes grouped only on decades
decade_slopes <- size_spectrum_indices %>% filter(str_detect(`group ID`, "decade"))

slope_timeline <- decade_slopes %>% 
  ggplot(aes(decade, b, color = survey_area)) +
  geom_line(aes(group = survey_area)) +
  geom_pointrange(aes(ymin = confMin, ymax = confMax),
                  alpha = 0.6) +
  scale_x_discrete(labels = function(x) paste0(x, "'s")) +
  scale_color_gmri() +
  labs(x = "Decade", 
       y = "Size Spectrum Slope (b)", 
       color = "",
       subtitle = "Results from Area-Stratified Abundance")


# How do they track with temperature, (lagged or otherwise)
temp_decades <- regional_oisst %>% 
  group_by(decade, survey_area) %>% 
  summarise(across(c(sst, sst_clim, sst_anom), mean)) %>% 
  ungroup() 

temp_timeline <- ggplot(data = temp_decades, aes(decade, sst_anom, color = survey_area)) +
  geom_line(aes(group = survey_area)) +
  geom_point(alpha = 0.6) +
  scale_x_discrete(labels = function(x) paste0(x, "'s")) +
  scale_color_gmri() +
  labs(x = "Decade", 
       color = "",
       y = expression("Temperature Anomaly "~~degree~C), 
       subtitle = "Anomalies from 1982-2011 Climatology")


# slopes over temperatures
slope_timeline / temp_timeline

  


#####  Regional ISD  ####
# ISD = Individual Size Distribution

# The following code is used to plot the Individual Size Distributions (size spectra)
# curves for each region

# it takes forever. and is only useful to verify they don't look horrible


# Pull Data for different regions for ISD or binned slope comparison
region_groups <- nefsc_1g %>% 
  mutate(group_var = survey_area) %>% 
  split(.$group_var)



# Prep Georges Bank for plotting
GB <- region_groups$GB

# 1. Take each group and run totals for how fishes fit into each discrete size bin
gb_isd_prep <- isd_plot_prep(biomass_data = GB, 
                             stratified_abundance = T, 
                             min_weight_g = 1)




# 2. Pull out MLE coefficients that match the groups

# LME Slopes
gb_lme_results <- size_spectrum_indices %>% 
  filter(`group ID` == "only regions",
         survey_area == "GB")

# 3. Set global limits for plots
xlim_global <- c( 1, max(nefsc_1g$wmax) )

# Plot the individual size distributions
gb_isd <- ggplot_isd(isd_data_prepped =  gb_isd_prep, 
                     mle_results = gb_lme_results, 
                     abundance_vals = "stratified", 
                     plot_rects = F, 
                     show_pl_fit = T, 
                     xlim_global = xlim_global, 
                     group_name = "Georges Bank - All Data")

# # Plot one
gb_isd$obs_y
gb_isd$log10_y




#####  Binned Slope Plot  ####

# LME Slopes & log10 slopes
gb_lme_results <- size_spectrum_indices %>% 
  filter(`group ID` == "only regions",
         survey_area == "GB")



# Set up manual log10 bins

# assign the bins to the data
GB_l10 <- region_groups$GB %>% 
  assign_log10_bins()



# plotting manual bins
plot_log10_ss(GB_l10)




####  Single Year Comparison  ####

# Georges bank 2010 as tester
gb_2010 <- nefsc_1g %>% 
  filter(survey_area == "GB",
         Year == 2010) %>% 
  tibble()

# results from both methods
gb_2010_results <- size_spectrum_indices %>% filter(`group ID` == "single years * region") %>% 
  filter(survey_area == "GB",
         Year == "2010")

# ISD plot
gb2010_isd_prep <- isd_plot_prep(biomass_data = gb_2010, stratified_abundance = T, min_weight_g = 1)
gb2010_isd <- ggplot_isd(isd_data_prepped = gb2010_isd_prep, 
                         mle_results = gb_2010_results, 
                         abundance_vals = "stratified", 
                         plot_rects = TRUE, 
                         show_pl_fit = T, group_name = "Georges Bank - 2010")
# un-transformed y
gb2010_isd$obs_y 
gb2010_isd$log10_y 

# Binned version
gb2010_assigned <- assign_log10_bins(gb_2010)
plot_log10_ss(l10_assigned = gb2010_assigned)

 
####  Annual Patterns  ####

##### Yr * Area * Season  ####
# How do the different size spectrum groups play out across decades?

# Pull the group ID for the slopes grouped only on decades
yr_area_season_slopes <- size_spectrum_indices %>% 
  filter(`group ID`== "single years * season * region") %>% 
  mutate(season = fct_rev(season),
         survey_area = factor(survey_area, 
                              levels = c("GoM", "GB", "SNE", "MAB")),
         yr = as.numeric(as.character(Year)))


lme_method <- yr_area_season_slopes %>% 
  ggplot(aes(yr, b, color = survey_area)) +
  geom_line(aes(group = survey_area)) +
  geom_point(alpha = 0.6) +
  facet_grid(season~survey_area) +
  scale_color_gmri() +
  labs(x = "", 
       y = "Size Spectrum Slope (b)", 
       color = "",
       subtitle = "Results from Area-Stratified Abundance - Edwards Method")



l10_method <- yr_area_season_slopes %>% 
  ggplot(aes(yr, l10_slope_strat, color = survey_area)) +
  geom_line(aes(group = survey_area)) +
  geom_point(alpha = 0.6) +
  facet_grid(season~survey_area) +
  scale_color_gmri() +
  labs(x = "", 
       y = "Size Spectrum Slope (b)", 
       color = "",
       subtitle = "Results from Area-Stratified Abundance - log10 bins")

# methods comparison
lme_method / l10_method

##### Yr * Area  ####
# How do the different size spectrum groups play out across decades?

# Pull the group ID for the slopes grouped only on decades
yr_area_slopes <- size_spectrum_indices %>% 
  filter(`group ID`== "single years * region")
  


lme_method  <- yr_area_slopes %>% 
  ggplot(aes(yr, b, color = survey_area)) +
  geom_line(aes(group = survey_area)) +
  geom_point(alpha = 0.6) +
  facet_grid(~survey_area) +
  scale_color_gmri() +
  labs(x = "", 
       y = "Size Spectrum Slope (b)", 
       color = "",
       subtitle = "Results from Area-Stratified Abundance - Edwards Method")



l10_method <- yr_area_slopes %>% 
  ggplot(aes(yr, l10_slope_strat, color = survey_area)) +
  geom_line(aes(group = survey_area)) +
  geom_point(alpha = 0.6) +
  facet_grid(~survey_area) +
  scale_color_gmri() +
  labs(x = "", 
       y = "Size Spectrum Slope (b)", 
       color = "",
       subtitle = "Results from Area-Stratified Abundance - log10 bins")


# methods comparison
lme_method / l10_method


####  Seasonal Differences  ####

season_slopes <- size_spectrum_indices %>% 
  filter(`group ID` == "single years * season ")

season_slopes %>% 
  ggplot(aes(yr, b, color = season)) +
  geom_line(aes(group = season)) +
  geom_point(alpha = 0.6) +
  facet_grid(~season) +
  scale_color_gmri() +
  labs(x = "", 
       y = "Size Spectrum Slope (b)", 
       color = "",
       subtitle = "Overall Seasonal Trends - Edwards Method")

season_slopes %>% 
  ggplot(aes(yr, l10_slope_strat, color = season)) +
  geom_line(aes(group = season)) +
  geom_point(alpha = 0.6) +
  facet_grid(~season) +
  scale_color_gmri() +
  labs(x = "", 
       y = "Size Spectrum Slope (b)", 
       color = "",
       subtitle = "Overall Seasonal Trends - log10 bins")


