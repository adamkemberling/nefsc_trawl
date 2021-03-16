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
box_paths <- research_access_paths(os.use = "unix")
mills_path <- box_paths$mills
res_path <- box_paths$res
nmfs_path  <- shared.path("unix", "RES_Data", "NMFS_Trawl")


####__  Data  ####


####__  Size Spectra Builds  ####

# Load data build function
source(here("R/01_nefsc_ss_build.R"))


# Data we just received in 2021 with errors located and corrected
load(paste0(nmfs_path, "NEFSC_BTS_all_seasons_03032021.RData"))
survdat_21 <- survey$survdat %>% clean_names()


# Run cleanup
survdat_21 <- survdat_prep(survdat = survdat_21) %>% 
  mutate(survdat_source = "survdat_2021")



# export for lindsay 3/15/2021
#write_csv(survdat_21, here("data/march_2021_survdat_filtered.csv"))




# Load cleaned data, and drop species with off coefficients
weights_20 <- survdat_prep(survdat_source = "2020") %>% 
  add_lw_info(cutoff = T) %>% 
  add_area_stratification(include_epu = F)



# Here are the main columns for weights, stratified weights, etc
weights_20 %>% select(
  est_year, stratum, epu, id, comname, biom_adj, # station and catch info
  numlen_adj, ind_weight_kg, sum_weight_kg,      # number at length, individual and group weights
  abund_tow_s, #abund_tow_epu,                    # abundance / number of tows in that strata that year
  biom_tow_s, #biom_tow_epu,                      # biomass / number of tows in that strata that year
  expanded_abund_s, #expanded_abund_epu,         # projected abundance, weighted by stratum cpue, expanded to total area
  expanded_biom_s, #expanded_biom_epu,           # projected biomass, weighted by stratum cpue, expanded to total area
  expanded_lwbio_s#, expanded_lwbio_epu          # projected biomass as l-w biomass of individuals * stratified abundance
)



####__  Export SS Prepped Data  ####

# # export for speedy recovery
# write_csv(weights_16, here::here("data/ss_prepped_data/survdat_2016_ss.csv"))
# write_csv(weights_19, here::here("data/ss_prepped_data/survdat_2019_ss.csv"))
# write_csv(weights_20, here::here("data/ss_prepped_data/survdat_2020_ss.csv"))


####____________________####
####  1. Comparing Annual Differences  ####

# run summaries
# summ_16 <- agg_species_metrics(weights_16, est_year, include_epu = F) %>% mutate(source = "2016 survdat")
# summ_19 <- agg_species_metrics(weights_19, est_year, include_epu = F) %>% mutate(source = "2019 survdat")
summ_20 <- agg_strat_metrics(weights_20, est_year, area_stratification = "stratum") %>% mutate(source = "2020 survdat")
# summs <- bind_rows(list(summ_16, summ_19, summ_20))
# summs <- bind_rows(list(summ_19, summ_20))
summs <- bind_rows(summ_20)



# Total Fish Caught - Actual in survey
summs %>% 
  ggplot(aes(est_year, sum_abund)) +
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
  geom_line(aes(est_year, strat_abundance_s, color = "Stratum Effort/Area Weighted"), show.legend = F) +
  #geom_line(aes(est_year, strat_abundance_epu, color = "EPU Effort/Area Weighted")) +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(x = "", y = "Effort & Area Expanded Total Abundance")



# Total Projected Biomass - whole survey
summs  %>% 
  ggplot() +
    geom_line(aes(est_year, lw_strat_biomass_s/ 1000000, color = "Stratum Weighted"), show.legend = F) +
    #geom_line(aes(est_year, lw_strat_biomass_epu/ 1000000, color = "EPU Weighted")) +
    facet_wrap(~`Research Vessel`, scales = "free") +
    scale_y_continuous(labels = scales::comma_format()) +
    labs(x = "", y = "Projected Total Biomass (million kg)\n L-W Derived")



# Total Projected Biomass - whole survey
summs %>% 
  ggplot() +
  geom_line(aes(est_year, fscs_strat_biomass_s/ 1000000, color = "Stratum Weighted"), show.legend = F) +
  # geom_line(aes(est_year, fscs_strat_biomass_epu/ 1000000, color = "EPU Weighted")) +
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












####  Mapping the trawl areas, exists elsewhere but needed for something  ####

# Load polygons for mapping
library(rnaturalearth)
library(sf)
box_paths <- research_access_paths()
res_path <- box_paths$res
new_england <- ne_states("united states of america") %>% st_as_sf(crs = 4326) 
canada <- ne_states("canada") %>% st_as_sf(crs = 4326)
trawl_strata <- read_sf(str_c(res_path, "Shapefiles/BottomTrawlStrata/BTS_Strata.shp"))
trawl_strata <- trawl_strata %>%  clean_names()

# Stratum Key for filtering specific areas
strata_key <- list(
  "Georges Bank"          = as.character(13:23),
  "Gulf of Maine"         = as.character(24:40),
  "Southern New England"  = str_pad(as.character(1:12), width = 2, pad = "0", side = "left"),
  "Mid-Atlantic Bight"    = as.character(61:76))

# Add labels to the data
trawl_strata <- trawl_strata %>%
  mutate(
    strat_num = str_sub(strata, 2, 3),
    survey_area =  case_when(
      strat_num %in% strata_key$`Georges Bank`         ~ "Georges Bank",
      strat_num %in% strata_key$`Gulf of Maine`        ~ "Gulf of Maine",
      strat_num %in% strata_key$`Southern New England` ~ "Southern New England",
      strat_num %in% strata_key$`Mid-Atlantic Bight`   ~ "Mid-Atlantic Bight",
      TRUE                                             ~ "not found"))


# Optional, Use strata_select to pull the strata we want individually
strata_select <- c(
  strata_key$`Georges Bank`,
  strata_key$`Gulf of Maine`,
  strata_key$`Southern New England`,
  strata_key$`Mid-Atlantic Bight`)


# Filtering with strata_select
trawl_strata <- trawl_strata %>% 
  filter(
    strata >= 01010,
    strata <= 01760,
    strata != 1310,
    strata != 1320,
    strata != 1330,
    strata != 1350,
    strata != 1410,
    strata != 1420,
    strata != 1490,
    strat_num %in% strata_select)


ggplot() +
  geom_sf(data = trawl_strata, aes(fill = survey_area)) +
  geom_sf(data = new_england) +
  geom_sf(data = canada) +
  coord_sf(xlim = c(-76, -66), ylim = c(35,46)) +
  theme(legend.title = element_blank())



# Load and reshape the Gulf of Maine timeseries for Jason Johnston

gom_temps <- read_csv(str_c(box_paths$okn, "oisst/regional_timeseries/nmfs_trawl_regions/OISSTv2_anom_gulf_of_maine.csv"))

# get yearxmonth means
monthlyavg <- gom_temps %>% 
  mutate(year = lubridate::year(time), 
         month = lubridate::month(time)) %>% 
  group_by(year, month) %>% 
  summarise(
    mean_sst = mean(sst, na.rm = T),
    historic_avg = mean(sst_clim, na.rm = T),
    temp_anomaly = mean_sst - historic_avg,
    climate_ref_period = "1982-2011",
    data_source = "NOAA OISST",
    region = "Gulf of Maine",
    region_def = "NMFS Trawl Strata 24-40")

# # Export csv for Jason to use
# write_csv(monthlyavg, here("data/GoM_monthly_oisst.csv"))


ggplot() +
  geom_sf(data = filter(trawl_strata, survey_area == "Gulf of Maine"), aes(fill = survey_area)) +
  geom_sf(data = new_england) +
  geom_sf(data = canada) +
  coord_sf(xlim = c(-71, -65.5), ylim = c(41,46)) +
  theme(legend.title = element_blank())
