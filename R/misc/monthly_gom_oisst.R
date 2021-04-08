# Creating Monthly OISST Averages for Gulf of Maine
# this was done as a one-off task for Jason Johnston to give 
# his class temperature data to work with

# Packages




# Processing

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