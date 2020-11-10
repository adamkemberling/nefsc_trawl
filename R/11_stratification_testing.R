####  Testing Sean Luceys Stratification Functions  ####
#### 11/4/2020
####


####  Load Packages  ####
library(janitor)
library(ecodata)
library(here)
library(gmRi)
library(rgdal)
library(sf)
library(sp)
library(data.table)
library(tidyverse)


####  Resource Paths  ####
mills_path <- shared.path(os.use = "unix", group = "Mills Lab", folder = NULL)
res_path <- shared.path(os.use = "unix", group = "RES Data", folder = NULL)


####  Load Data  ####
source(here("R/kathy_ss_code/slucey_survdat_functions.R"))  # Sean Lucey Survdat functions
load(here("data/NEFSC/Survdat_Nye_Aug 2020.RData"))         # Survdat Nye 2020
epu_sf <- ecodata::epu_sf                                   # EPU Shapefiles

# Trawl Survey Strata we typically use
trawl_strata <- survey_strata <- read_sf(str_c(res_path, "Shapefiles/BottomTrawlStrata/BTS_Strata.shp"))  %>% 
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

# Stratum area values used currently
strat_area_csv <- read_csv(str_c(mills_path, "Projects/NSF_CAccel/Data/strata area.csv"), 
                           col_types = cols()) %>% 
  filter(stratum >= 01010 ,
         stratum <= 01760,
         stratum != 1310,
         stratum != 1320,
         stratum != 1330,
         stratum != 1350,
         stratum != 1410,
         stratum != 1420,
         stratum != 1490) 









####  1. Getting Area from Strata  ####

trawl_sp     <- sf::as_Spatial(trawl_strata)
epu_sp       <- sf::as_Spatial(epu_sf)
strata_areas <- getarea(trawl_sp, strat.col = "strata") %>% setNames(c("stratum", "lucey_area"))
epu_areas    <- getarea(epu_sp, strat.col = "EPU")


# Compare to csv for stratum areas
area_test <- full_join(strata_areas, strat_area_csv, cy = c("stratum"))
conv_nm_to_km <- 1.852
conv_sqnm_to_sqkm <- 3.4299

test <- area_test %>% 
  select(-c(type, stratum_num)) %>% 
  mutate(
  nm_to_km         = area * conv_nm_to_km,         # if value were nautical miles
  sqnm_to_sqkm     = area * conv_sqnm_to_sqkm,     # if csv were in square nautical miles
  off_mark         = lucey_area - sqnm_to_sqkm) # how far the ratio is from actual conversion

range(test$off_mark, na.rm = T)
hist(test$off_mark, main = "Sean Lucey Calculation (sq-km) - Spreadsheet area (sq-nm converted to sq-km)")





####  2. Poststrat  ####
survdat_poststrat <- survdat %>% 
  mutate(
    LAT = DECDEG_BEGLAT,
    LON = DECDEG_BEGLON) %>% 
  poststrat(survdat = ., stratum = epu_sp, strata.col = "EPU") %>% 
  rename(EPU = newstrata,
         YEAR = EST_YEAR)


####  3. Stratprep  ####
survdat_prepped <- stratprep(survdat = survdat_poststrat, areas = epu_areas, strat.col = "EPU", area.col = "Area")


####  4. stratmean  ####
survdat_stratmean <- stratmean(survdat_prepped)



####  5. Swept area  ####
survdat_swept_area <- sweptarea(survdat = survdat_prepped, stratmean = survdat_stratmean, strat.col = "EPU", area.col = "Area", group.col = "SVSPP")



# # Export Lucey Stratified Abundances
# write_csv(survdat_swept_area, here("data/slucey/slucey_survdat_EPU.csv"))


strat_data <- survdat_swept_area %>% 
  rename_all(tolower)


ann_summs <- strat_data %>% 
  group_by(year) %>% 
  summarise(total_abundance = sum(tot.abundance),
            total_biomass   = sum(tot.biomass),
            strat_abund     = sum(strat.abund),
            strat_biomass   = sum(strat.biomass))

ggplot(ann_summs, aes(year, total_abundance)) +
  geom_line() +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(x = "", y = "Total Abundance - area swept")

ggplot(ann_summs, aes(year, total_biomass)) +
  geom_line() +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(x = "", y = "Total Biomass - area swept")


