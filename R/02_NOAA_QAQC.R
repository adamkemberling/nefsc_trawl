####  Verification that NOAA Data isn't changing on us ###
#### look at gross metrics like stratified mean abundance etc.


####  Comparing Groundfish Data from 2020 to 2016  ####



####____________________####
###__ Packages  ####
library(here)
library(janitor)
library(gmRi)
library(patchwork)
library(tidyverse)



# Load the build code and stratification function
mills_path <- shared.path(os.use = "unix", group = "Mills Lab", folder = NULL)
res_path <- shared.path(os.use = "unix", group = "RES Data", folder = NULL)



####__  Data  ####


# 2019 survdat data
load(str_c(mills_path, "Data/Survdat_Nye_allseason_2019.RData"))
trawl_19 <- survdat %>% clean_names()
rm(survdat)


# # 2020 Raw data
# trawl_20 <- read_csv(here("data/NEFSC/2020survdat_nye.csv"), 
#                      guess_max = 1e6, 
#                      col_types = cols()) %>% clean_names()







####__  Size Spectra Builds  ####
source(here("R/01_nefsc_ss_build.R"))
weights_16 <- load_2016_ss_data() 
weights_19 <- load_ss_data(survdat = trawl_19, survdat_source = "2019")
weights_20 <- load_ss_data(survdat_source = "2020")


# Here are the main columns for weights, stratified weights, etc
weights_20 %>% select(id, comname, biom_adj, numlen_adj, ind_weight_kg, sum_weight_kg, 
                      strat_wt_bio_fscs, strat_wt_bio_lw, fscs_strat_bio, lw_strat_bio)




####__  Export SS Prepped Data  ####

# # export for speedy recovery
# write_csv(weights_16, here::here("data/ss_prepped_data/survdat_2016_ss.csv"))
# write_csv(weights_19, here::here("data/ss_prepped_data/survdat_2019_ss.csv"))
# write_csv(weights_20, here::here("data/ss_prepped_data/survdat_2020_ss.csv"))


####____________________####
####  1. Comparing Annual Differences  ####

# run summaries
summ_16 <- ss_annual_summary(weights_16) %>% mutate(source = "2016")
summ_19 <- ss_annual_summary(weights_19) %>% mutate(source = "2019")
summ_20 <- ss_annual_summary(weights_20) %>% mutate(source = "2020")
summs <- bind_rows(list(summ_16, summ_19, summ_20))


# Total Biomass - L-W
summs %>% 
  ggplot(aes(est_year, lw_biomass_kg, color = source)) +
  geom_line() +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(x = "", y = "Total Biomass (kg) \n (L-W Regressions)")

# Total Biomass - FSCS
summs %>% 
  ggplot(aes(est_year, fscs_biomass_kg, color = source)) +
  geom_line() +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(x = "", y = "Total Biomass (kg) \n (FSCS Haul Weights)")

# Number of Species
summs %>% 
  ggplot(aes(est_year, n_species, color = source)) +
  geom_line() +
  labs(x = "", y = "Number of Species")


# stratum weighted catch /  nautical mile - length weight derived
summs %>% 
  ggplot(aes(est_year, lw_strat_biomass, color = source)) +
    geom_line() +
    scale_y_continuous(labels = scales::comma_format()) +
    labs(x = "", y = "Effort & Area Stratified Biomass \n L-W Derived")


# stratum weighted catch /  nautical mile - length weight derived
summs %>% 
  ggplot(aes(est_year, fscs_strat_biomass, color = source)) +
  geom_line() +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(x = "", y = "Effort & Area Stratified Biomass \n FSCS Derived")


####  2. Regional Differences  ####

# run summaries
summ_16 <- ss_regional_differences(weights_16) %>% mutate(source = "2016")
summ_19 <- ss_regional_differences(weights_19) %>% mutate(source = "2019")
summ_20 <- ss_regional_differences(weights_20) %>% mutate(source = "2020")
reg_summs <- bind_rows(list(summ_16, summ_19, summ_20))


# Total Biomass
p1 <- reg_summs %>% 
  ggplot(aes(est_year, lw_biomass_kg, color = source)) +
  geom_line(show.legend = F) +
  scale_y_continuous(labels = scales::comma_format()) +
  facet_wrap(~survey_area, ncol = 1, scales = "free") +
  labs(x = "", y = "Total Biomass \n (L-W Regressions)")

# Total Biomass - FSCS
p2 <- reg_summs %>% 
  ggplot(aes(est_year, fscs_biomass_kg, color = source)) +
  geom_line() +
  scale_y_continuous(labels = scales::comma_format()) +
  facet_wrap(~survey_area, ncol = 1 , scales = "free") +
  labs(x = "", y = "Total Biomass \n (FSCS Haul Weights)")

# effort
p3 <- reg_summs %>% 
  ggplot(aes(est_year, lw_strat_biomass, color = source)) +
  geom_line(show.legend = F) +
  scale_y_continuous(labels = scales::comma_format()) +
  facet_wrap(~survey_area, ncol = 1, scales = "free") +
  labs(x = "", y = "Stratified Abundance - L-W")

# Species 
p4 <- reg_summs %>% 
  ggplot(aes(est_year, fscs_strat_biomass, color = source)) +
  geom_line(show.legend = F) +
  scale_y_continuous(labels = scales::comma_format()) +
  facet_wrap(~survey_area, ncol = 1, scales = "free") +
  labs(x = "", y = "Stratified Abundance - FSCS")

p1 + p2 + p3 + p4
