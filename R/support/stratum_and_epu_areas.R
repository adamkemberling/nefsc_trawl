# Getting Areas in km2 for EPU's and Trawl Stratum 
# Used for stratified area estimations


####____________________####
#### Prep EPU polygon areas  ####
# This only needs to be done once, but before the data prep function
# Also found on Box/RES_Data/NMFS_trawl/slucey_functions/


# source(here("R/kathy_ss_code/slucey_survdat_functions.R"))
# 
# # packages just for this:
# library(sf)
# library(sp)
# library(data.table)
# 
# # paths:
# res_path   <- shared.path(os.use = "unix", group = "RES Data", folder = NULL)
# mills_path <- shared.path(os.use = "unix", group = "Mills Lab", folder = "")
# 
# # Load EPU Shapefiles, get area in sq km
# epu_sf    <- ecodata::epu_sf
# epu_sp    <- sf::as_Spatial(epu_sf)
# epu_areas <- getarea(epu_sp, strat.col = "EPU") %>% setNames(c("epu", "epu_area_km2"))
# write_csv(epu_areas, str_c(mills_path, "Projects/NSF_CAccel/Data/EPU_areas_km2.csv"))
# 
# # Load strata shapefiles, get area, also assign their EPU
# 
# trawl_strata <- read_sf(str_c(res_path, "Shapefiles/BottomTrawlStrata/BTS_Strata.shp"))
# trawl_strata <- trawl_strata %>%  clean_names()
# trawl_sp     <- sf::as_Spatial(trawl_strata)
# # get area
# strata_areas <- getarea(trawl_sp, strat.col = "strata") %>% setNames(c("stratum", "s_area_km2"))
# write_csv(strata_areas, str_c(mills_path, "Projects/NSF_CAccel/Data/strata_areas_km2.csv"))
