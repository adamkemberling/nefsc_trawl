# Pulling data for a single species, keeping all stations



# Packages
library(targets)
library(gmRi)
library(tidyverse)

# Load data - 
# after general cleanup steps, before length-weight or area stratification steps
tar_load(nefsc_clean)


# All Unique stations
stations <- nefsc_clean %>% 
  distinct(id,station, cruise6, svvessel, est_towdate, est_year, season, 
           stratum, tow, decdeg_beglon, decdeg_beglat, 
           avgdepth, surftemp, surfsalin, bottemp, botsalin)



# Pull out (menhaden) catch:
menhad <- nefsc_clean %>% 
  filter(comname == "menhaden") %>% 
  select(id, comname, catchsex, abundance, biomass, length, numlen, numlen_adj)


# Re-join stations to get stations where they weren't caught
menhad_all <- left_join(stations, menhad, by = "id") %>% 
  mutate(comname    = "menhaden",
         abundance  = ifelse(is.na(abundance), 0, abundance),
         biomass    = ifelse(is.na(biomass), 0, biomass),
         numlen     = ifelse(is.na(numlen), 0, numlen),
         numlen_adj = ifelse(is.na(numlen_adj), 0, numlen_adj)) %>% 
  select(id, everything())


# Save out
write_csv(menhad_all, here::here("data/single_species/survdat_menhaden_2019.csv"))

# Column explanation
data_dict <- data.frame(column_name = names(menhad_all), description = "")

# Write them up
data_dict %>% 
  mutate(
    description = case_when(
      column_name == "id" ~ "Unique identifier for tow",
      column_name == "cruise6" ~ "Six digit identifier of cruise",
      column_name == "station" ~ "Station identifier within cruise",
      column_name == "stratum" ~ "Survey stratum number",
      column_name == "tow" ~ "Tow number",
      column_name == "svvessel" ~ "Survey vessel",
      column_name == "est_year" ~ "Year",
      column_name == "season" ~ "Season",
      column_name == "decdeg_beglat" ~ "Latitude of starting position (decimal degrees)",
      column_name == "decdeg_beglon" ~ "Longitude of starting position (decimal degrees)",
      column_name == "est_towdate" ~ "Datetime for start of the station sampling",
      column_name == "avgdepth" ~ "Average depth",
      column_name == "surfsalin" ~ "Surface salinity",
      column_name == "surftemp" ~ "Surface temperature",
      column_name == "bottemp" ~ "Bottom temperature",
      column_name == "botsalin" ~ "Bottom salinity",
      column_name == "comname" ~ "Common name",
      column_name == "catchsex" ~ "Species sex code",
      column_name == "abundance" ~ "Total abundance for species and that sex, of any size caught at station",
      column_name == "biomass" ~ "Total biomass for species and that sex, of any size caught at station",
      column_name == "length" ~ "Length (cm) of a subset of total abundance",
      column_name == "numlen" ~ "Number of individuals for species and sex, that were of the length recorded",
      column_name == "numlen_adj" ~ "Adjusted numlen value, adjusted to ensure that the sum of numlen equals abundance column",
    )
  )
