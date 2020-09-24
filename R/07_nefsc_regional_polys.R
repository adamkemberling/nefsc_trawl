# Goal: Use the stratum lists to create polygons for calculating oisst anomalies from
# do the anomalies using python and Jupyter Notebooks
# Keep Everything together in this repository
# Python navigation to this project from terminal
# cd -> this repo
# conda activate py36
# OISST processing scripts can be found at Documents/oisst_mainstays

####  Packages  ####
library(here)
library(gmRi)
library(sf)
library(janitor)
library(sizeSpectra)
library(patchwork)
library(tidyverse)
library(rnaturalearth)


####  Support Functions  ####
source(here("R/support/sizeSpectra_support.R"))

# Set Theme
theme_set(theme_minimal())


# File Paths
mills_path <- shared.path(os.use = "unix", group = "Mills Lab", folder = NULL)
res_path <- shared.path(os.use = "unix", group = "RES Data", folder = NULL)


####  Data  ####

# Stratum Shapefiles
survey_strata <- read_sf(str_c(res_path, "Shapefiles/BottomTrawlStrata/BTS_Strata.shp"))  %>% 
  clean_names() %>% 
  filter(strata >= 01010 ,
         strata <= 01760,
         strata != 1310,
         strata != 1320,
         strata != 1330,
         strata != 1350,
         strata != 1410,
         strata != 1420,
         strata != 1490) 


# Key to which strata = which regions
strata_key <- list(
  "Georges Bank"          = as.character(13:23),
  "Gulf of Maine"         = as.character(24:40),
  "Southern New England"  = str_pad(as.character(1:12), width = 2, pad = "0", side = "left"),
  "Mid-Atlantic Bight"    = as.character(61:76))


# Assign Areas
strata_assigned <- survey_strata %>% 
  mutate(
    strata = str_pad(strata, width = 5, pad = "0", side = "left"),
    area = case_when(
    str_sub(strata, 3, 4) %in% strata_key$`Georges Bank` ~ "Georges Bank",
    str_sub(strata, 3, 4) %in% strata_key$`Gulf of Maine` ~ "Gulf of Maine",
    str_sub(strata, 3, 4) %in% strata_key$`Southern New England` ~ "Southern New England",
    str_sub(strata, 3, 4) %in% strata_key$`Mid-Atlantic Bight` ~ "Mid-Atlantic Bight",
    TRUE ~ "Outside Major Study Areas"
  )) 

# Plot to check
new_england <- ne_states("united states of america") %>% st_as_sf(crs = 4326) 

ggplot() +
  geom_sf(data = new_england) +
  geom_sf(data = strata_assigned, aes(fill = area)) +
  coord_sf(xlim = c(-79, -65.5), ylim = c(32, 44.5), expand = FALSE)

# Export Single Areas
gb <- strata_assigned %>% filter(area == "Georges Bank")
gm <- strata_assigned %>% filter(area == "Gulf of Maine")
mab <- strata_assigned %>% filter(area == "Mid-Atlantic Bight")
sne <- strata_assigned %>% filter(area == "Southern New England")


# Export the collection, try and do the regional timelines by area
st_write(strata_assigned, here("data/polys/trawl_regions_collection.geojson"))
st_write(filter(strata_assigned, area == "Gulf of Maine"), here("data/polys/nmfs_trawl_gulf_of_maine.geojson"))
st_write(filter(strata_assigned, area == "Georges Bank"), here("data/polys/nmfs_trawl_georges_bank.geojson"))
st_write(filter(strata_assigned, area == "Southern New England"), here("data/polys/nmfs_trawl_southern_new_england.geojson"))
st_write(filter(strata_assigned, area == "Mid-Atlantic Bight"), here("data/polys/nmfs_trawl_mid_atlantic_bight.geojson"))





