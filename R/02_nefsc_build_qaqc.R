####  Verification that NOAA Data isn't changing on us ###
#### look at gross metrics like stratified mean abundance etc.


####  Comparing Groundfish Data from 2020 to 2016  ####



####____________________####
###__ Packages  ####
library(targets)
library(here)
library(janitor)
library(gmRi)
library(patchwork)
library(rnaturalearth)
library(sf)
library(tidyverse)



# Load the build code and stratification function
box_paths  <- research_access_paths(os.use = "unix")
mills_path <- box_paths$mills
res_path   <- box_paths$res
nmfs_path  <- shared.path("unix", "RES_Data", "NMFS_Trawl")


####__  Size Spectra Builds  ####

# Load data build function
source(here("R/01_nefsc_ss_build_nodrop.R"))



####  Comparing the results of original build vs. no column drops
# # Run cleanup functions (drop or not drop) for comparison
# source(here("R/01_old_nefsc_ss_build.R"))
# survdat_21 <- survdat_prep(survdat = survdat_raw) %>%
#   mutate(survdat_source = "survdat_2021")
# 
# # Run the no drop cleanup
# allcols_21 <- survdat_prep_nodrop(survdat = survdat_raw)
# 
# # leftovers
# lefties <- anti_join(survdat_21, allcols_21)
# 
# # export for lindsay 3/15/2021
# #write_csv(survdat_21, here("data/march_2021_survdat_filtered.csv"))
# #write_csv(allcols_21, here("data/march_2021_survdat_allcols.csv"))
####


####__  Loading Clean Data  ####

# Loading from targets:
tar_load(nefsc_stratified)
weights_21 <- nefsc_stratified


####__  Export SS Prepped Length-Weight Data  ####


# # Load cleaned data, and drop species with off coefficients
# weights_20 <- survdat_prep(survdat = survdat_raw) %>% 
#   add_lw_info(cutoff = T) %>% 
#   add_area_stratification(include_epu = F)

# # export for speedy recovery
# write_csv(weights_16, here::here("data/ss_prepped_data/survdat_2016_ss.csv"))
# write_csv(weights_19, here::here("data/ss_prepped_data/survdat_2019_ss.csv"))
# write_csv(weights_20, here::here("data/ss_prepped_data/survdat_2020_ss.csv"))
# write_csv(weights_21, here::here("data/ss_prepped_data/survdat_2021_ss.csv"))


####____________________####
####  1. Comparing Annual Differences  ####

# run summaries
# summ_16 <- ss_annual_summary(weights_16, include_epu = F) %>% mutate(source = "2016 survdat")
# summ_19 <- ss_annual_summary(weights_19, include_epu = F) %>% mutate(source = "2019 survdat")
# summ_20 <- ss_annual_summary(weights_20, include_epu = F) %>% mutate(source = "2020 survdat")
  summ_21 <- ss_annual_summary(weights_21, include_epu = F) %>% mutate(source = "2021 survdat")

# build aggregates
# summs <- bind_rows(list(summ_16, summ_19, summ_20))
# summs <- bind_rows(list(summ_19, summ_20))

# reduce to single table
summs <- bind_rows(summ_21)



# Total Fish Caught - Actual in survey
summs %>% 
  ggplot(aes(est_year, total_survey_abund)) +
  geom_line() +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(x = "", y = "Survey Total Abundance")

# Total Biomass - From Survey
summs %>% 
  ggplot() +
  geom_line(aes(est_year, fscs_biomass_kg, color = "L-W Regression Biomass")) +
  geom_line(aes(est_year, lw_biomass_kg, color = "Measured Biomass")) +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(x = "", y = "Total Survey Biomass (kg)")


# Total Fish Caught - Projected
summs %>% 
  ggplot() +
  geom_line(aes(est_year, strat_abundance_s)) +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(x = "", y = "Effort & Area Expanded Total Abundance")



# Total Projected Biomass - whole survey
summs  %>% 
  ggplot() +
    geom_line(aes(est_year, lw_strat_biomass_s/ 1000000)) +
    geom_vline(xintercept = 2008.5, linetype = 2, alpha = 0.5) +
    scale_y_continuous(labels = scales::comma_format()) +
    labs(x = "", y = "Projected Total Biomass (million kg)\n L-W Derived")



# Total Projected Biomass - whole survey
summs %>% 
  ggplot() +
  geom_line(aes(est_year, fscs_strat_biomass_s/ 1000000)) +
  geom_vline(xintercept = 2008.5, linetype = 2, alpha = 0.5) +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(x = "", y = "Projected Total Biomass (million kg)\n FSCS Derived")





####  2. Regional Differences  ####

# run summaries
# area_summ_16 <- ss_regional_differences(weights_16) %>% mutate(source = "2016 survdat")
# area_summ_19 <- ss_regional_differences(weights_19) %>% mutate(source = "2019 survdat")
# area_summ_20 <- ss_regional_summary(weights_20) %>% mutate(source = "2020 survdat")
area_summ_21 <- ss_regional_summary(weights_21, epu_regions = F) %>% mutate(source = "2021 survdat")
# reg_summs <- bind_rows(list(area_summ_16, area_summ_19, area_summ_20))
reg_summs <- bind_rows(area_summ_21)

# Total Biomass
p1 <- reg_summs %>% 
  ggplot(aes(est_year, lw_biomass_kg)) +
  geom_line(show.legend = F) +
  scale_y_continuous(labels = scales::comma_format()) +
  facet_wrap(~survey_area, ncol = 1, scales = "free") +
  labs(x = "", y = "Total Survey Biomass (kg)\n(L-W Regressions)")

# Total Biomass - FSCS
p2 <- reg_summs %>% 
  ggplot(aes(est_year, fscs_biomass_kg)) +
  geom_line() +
  scale_y_continuous(labels = scales::comma_format()) +
  facet_wrap(~survey_area, ncol = 1 , scales = "free") +
  labs(x = "", y = "Total Survey Biomass (kg)\n(FSCS Haul Weights)")

# effort
p3 <- reg_summs %>% 
  ggplot() +
  geom_line( aes(est_year, lw_strat_biomass_s / 1000000), show.legend = F) +
  scale_y_continuous(labels = scales::comma_format()) +
  facet_wrap(~survey_area, ncol = 1, scales = "free") +
  labs(x = "", y = "Stratified Biomass (million kg)\nL-W")

# Species 
p4 <- reg_summs %>% 
  ggplot() +
  geom_line(aes(est_year, fscs_strat_biomass_s/ 1000000), show.legend = F) +
  scale_y_continuous(labels = scales::comma_format()) +
  facet_wrap(~survey_area, ncol = 1, scales = "free") +
  labs(x = "", y = "Stratified Biomass (million kg)\nFSCS")

p1 + p2 + p3 + p4















