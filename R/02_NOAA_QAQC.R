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
# weights_16 <- load_2016_ss_data() 
# weights_19 <- load_ss_data(survdat = trawl_19, survdat_source = "2019")
weights_20 <- load_ss_data(survdat_source = "2020")


# Here are the main columns for weights, stratified weights, etc
weights_20 %>% select(
  est_year, stratum, epu, id, comname, biom_adj, # station and catch info
  numlen_adj, ind_weight_kg, sum_weight_kg,      # number at length, individual and group weights
  abund_tow_s, abund_tow_epu,                    # abundance / number of tows in that strata that year
  biom_tow_s, biom_tow_epu,                      # biomass / number of tows in that strata that year
  expanded_abund_s, expanded_abund_epu,          # projected abundance, weighted by stratum cpue, expanded to total area
  expanded_biom_s, expanded_biom_epu,            # projected biomass, weighted by stratum cpue, expanded to total area
  expanded_lwbio_s, expanded_lwbio_epu)          # projected biomass as l-w biomass of individuals * stratified abundance


# column groups, by what they are and what stratification they are

# annual summaries, rough approach
test_summaries <- weights_20 %>% 
  group_by(est_year) %>% 
  summarise(across(c(numlen_adj, sum_weight_kg, biom_adj, abund_tow_s:expanded_lwbio_epu), ~ sum(.x, na.rm = T)),
            .groups = "keep") %>% 
  pivot_longer(names_to = "Metric", values_to = "total", cols = c(numlen_adj:expanded_lwbio_epu)) %>% 
  mutate(
    strat_level = case_when(
      str_detect(Metric, "_epu") ~ "epu",
      str_detect(Metric, "sum_weight_kg|biom_adj|numlen_adj") ~ "un-adjusted",
      str_detect(Metric, "_s") ~ "nmfs_strata"),
    type = case_when(
      str_detect(Metric, "kg|biom_adj")     ~ "survey biomass",
      str_detect(Metric, "tow")             ~ "catch / tow",
      str_detect(Metric, "wt")              ~ "weighted cpue",
      str_detect(Metric, "expanded")        ~ "expanded total",
      str_detect(Metric, "numlen_adj|abund") ~ "abundance"),
    `measurement source` = case_when(
      str_detect(Metric, "lwbio|sum_weight_kg")  ~ "LW Regression",
      str_detect(Metric, "biom|biom_adj")        ~ "Shipboard Biomass",
      str_detect(Metric, "abund|numlen_adj")     ~ "Abundance"),
    vessel = ifelse(est_year > 2008, "HB", "AL")
  )


test_summaries %>% 
  split(.$vessel) %>%
  map(function(x){
    ggplot(x, aes(est_year, total)) +
    geom_line(aes(color = `measurement source`, linetype = )) +
    labs(x = "", y = "Annual Total") +
    facet_grid(type~strat_level, scales = "free")
  })



####__  Export SS Prepped Data  ####

# # export for speedy recovery
# write_csv(weights_16, here::here("data/ss_prepped_data/survdat_2016_ss.csv"))
# write_csv(weights_19, here::here("data/ss_prepped_data/survdat_2019_ss.csv"))
# write_csv(weights_20, here::here("data/ss_prepped_data/survdat_2020_ss.csv"))


####____________________####
####  1. Comparing Annual Differences  ####

# run summaries
#summ_16 <- ss_annual_summary(weights_16) %>% mutate(source = "2016 survdat")
#summ_19 <- ss_annual_summary(weights_19) %>% mutate(source = "2019 survdat")
summ_20 <- ss_annual_summary(weights_20) %>% mutate(source = "2020 survdat")
# summs <- bind_rows(list(summ_16, summ_19, summ_20))
# summs <- bind_rows(list(summ_19, summ_20))
summs <- bind_rows(summ_20)


# Number of Species
summs %>% 
  ggplot(aes(est_year, n_species, color = source)) +
  geom_line() +
  labs(x = "", y = "Number of Species")

# Total Fish Caught - Actual in survey
summs %>% 
  ggplot(aes(est_year, total_survey_abund, color = source)) +
  geom_line() +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(x = "", y = "Survey Total Abundance")

# Total Biomass - From Survey
summs %>% 
  ggplot() +
  geom_line(aes(est_year, lw_biomass_kg, color = "L-W Regression Biomass")) +
  geom_line(aes(est_year, fscs_biomass_kg, color = "Measured Biomass")) +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(x = "", y = "Total Survey Biomass (kg)")


# Total Fish Caught - Projected
summs %>% 
  ggplot() +
  geom_line(aes(est_year, strat_abundance_s, color = "Stratum Effort/Area Weighted")) +
  geom_line(aes(est_year, strat_abundance_epu, color = "EPU Effort/Area Weighted")) +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(x = "", y = "Effort & Area Expanded Total Abundance")



# Total Projected Biomass - whole survey
summs  %>% 
  ggplot() +
    geom_line(aes(est_year, lw_strat_biomass_s/ 1000000, color = "Stratum Weighted")) +
    geom_line(aes(est_year, lw_strat_biomass_epu/ 1000000, color = "EPU Weighted")) +
    facet_wrap(~`Research Vessel`, scales = "free") +
    scale_y_continuous(labels = scales::comma_format()) +
    labs(x = "", y = "Projected Total Biomass (million kg)\n L-W Derived")



# Total Projected Biomass - whole survey
summs %>% 
  ggplot() +
  geom_line(aes(est_year, fscs_strat_biomass_s/ 1000000, color = "Stratum Weighted")) +
  geom_line(aes(est_year, fscs_strat_biomass_epu/ 1000000, color = "EPU Weighted")) +
  facet_wrap(~`Research Vessel`, scales = "free") +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(x = "", y = "Projected Total Biomass (million kg)\n FSCS Derived")





####  2. Regional Differences  ####

# run summaries
# summ_16 <- ss_regional_differences(weights_16) %>% mutate(source = "2016 survdat")
# summ_19 <- ss_regional_differences(weights_19) %>% mutate(source = "2019 survdat")
summ_20 <- ss_regional_summary(weights_20) %>% mutate(source = "2020 survdat")
# reg_summs <- bind_rows(list(summ_16, summ_19, summ_20))
reg_summs <- bind_rows(summ_20)

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
  geom_line( aes(est_year, lw_strat_biomass_s / 1000000, color = "Stratum"), show.legend = F) +
  geom_line( aes(est_year, lw_strat_biomass_epu/ 1000000, color = "EPU"), show.legend = F) +
  scale_y_continuous(labels = scales::comma_format()) +
  facet_wrap(~survey_area, ncol = 1, scales = "free") +
  labs(x = "", y = "Stratified Biomass (million kg)\nL-W")

# Species 
p4 <- reg_summs %>% 
  ggplot() +
  geom_line(aes(est_year, fscs_strat_biomass_s/ 1000000, color = "Stratum"), show.legend = F) +
  geom_line(aes(est_year, fscs_strat_biomass_epu/ 1000000, color = "EPU"), show.legend = F) +
  scale_y_continuous(labels = scales::comma_format()) +
  facet_wrap(~survey_area, ncol = 1, scales = "free") +
  labs(x = "", y = "Stratified Biomass (million kg)\nFSCS")

p1 + p2 + p3 + p4
