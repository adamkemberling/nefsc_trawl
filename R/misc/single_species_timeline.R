
# Single Species Abundance/Biomass and Stratified Abundance/Biomass
# Load target pipeline to display single species timelines



####. Packages. ####
library(targets)
library(here)
library(sf)
library(gmRi)
library(patchwork)
library(tidyverse)
library(scales)



####. Data. ####
tar_load("survdat_clean")

# Set Common Name
species_name <- "american lobster"

# Set seasons to include
season_opts <- c("Spring", "Fall")

# Set Theme
theme_set(theme_gmri())


####. Station Level Aggregation. ####

# biomass and abundance numbers should come from station specific totals
station_totals <- survdat_clean %>% 
  filter(season %in% season_opts) %>% 
  distinct(est_year, survey_area, stratum, tow, est_towdate, season, comname, catchsex, .keep_all = T) %>%
  group_by(est_year, survey_area, stratum, tow, est_towdate, season, 
           avgdepth, bottemp, decdeg_beglat, decdeg_beglon, comname) %>% 
  summarise(
      station_abundance = sum(abundance, na.rm = T),
      station_biomass_kg = sum(biomass_kg, na.rm = T), 
      .groups = "drop")

# ^ All species present to preserve stations that may not have hit target species
single_species_raw <- filter(station_totals, comname == "american lobster")

# Now we can total up for the year
single_species_annual <- single_species_raw %>% 
  group_by(comname, est_year) %>% 
  summarise(
    total_abundance = sum(station_abundance, na.rm = T),
    total_biomass_kg = sum(station_biomass_kg, na.rm = T),
    .groups = "drop") %>% 
  pivot_longer(cols = c(total_abundance, total_biomass_kg),
               names_to = "metric",
               values_to = "annual_total")
  
# Plot the Annual Totals
annual_total_plot <- single_species_annual  %>% 
  ggplot(aes(est_year, annual_total)) +
  geom_line() +
  facet_grid(.~metric) +
  scale_y_continuous(labels = label_comma()) +
  labs(x = "Year", y = "Annual Total", title = species_name)


####. Area Stratified. ####


area_stratify_station_totals <- function(survdat_clean, box_location = "cloudstorage"){
    
    # https://noaa-edab.github.io/survdat/articles/calc_strat_mean.html
    
    # Switch for mac users with different box storage location
    path_fun <- boxpath_switch(box_location = box_location)
    
    
    ####  1. Import supplemental files  ####
    nmfs_path <- path_fun(box_group = "RES_Data", subfolder = "NMFS_trawl")
    
    # Stratum Area Information
    stratum_area_path <- stringr::str_c(nmfs_path, "Metadata/strata_areas_km2.csv")
    stratum_area      <- readr::read_csv(stratum_area_path, col_types = readr::cols())
    stratum_area      <- dplyr::mutate(stratum_area, stratum = as.character(stratum))
    
    
    
    ####  2. Set Constants:  ####
    
    # Area covered by an albatross standard tow in km2
    alb_tow_km2 <- 0.0384
    
    # catchability coefficient - ideally should change for species guilds
    q <- 1
    
    
    ####  Get Annual Stratum Effort, and Area Ratios
    # number of tows in each stratum by year
    # area of a stratum relative to total area of all stratum sampled that year
    
    
    ####  3. Stratum Area & Effort Ratios  ####
    
    # Merge in the area of strata in km2
    # (excludes ones we do not care about via left join)
    survdat_clean <- dplyr::left_join(survdat_clean, stratum_area, by = "stratum")
    survdat_clean <- dplyr::arrange(survdat_clean, survdat_clean, id)
    
    
    # Get Total area of all strata sampled in each year
    total_stratum_areas <- dplyr::group_by(survdat_clean, est_year)
    total_stratum_areas <- dplyr::distinct(total_stratum_areas, stratum, .keep_all = T)
    total_stratum_areas <- dplyr::summarise(total_stratum_areas,
                                            tot_s_area =  sum(s_area_km2, na.rm = T),
                                            .groups = "keep")
    total_stratum_areas <- dplyr::ungroup(total_stratum_areas)
    
    
    # Calculate strata area relative to total area i.e. stratio or stratum weights
    survdat_clean <- dplyr::left_join(survdat_clean, total_stratum_areas, by = "est_year")
    survdat_clean <- dplyr::mutate(survdat_clean, st_ratio = s_area_km2 / tot_s_area)
    
    
    # We have total areas, now we want effort within each
    # Number of unique tows per stratum, within each season
    yr_strat_effort <- dplyr::group_by(survdat_clean, est_year, season, stratum)
    yr_strat_effort <- dplyr::summarise(yr_strat_effort, strat_ntows = dplyr::n_distinct(id), .groups = "keep")
    yr_strat_effort <- dplyr::ungroup(yr_strat_effort)
    
    
    
    # Add those yearly effort counts back for later
    # (area stratified abundance)
    survdat_clean <- dplyr::left_join(survdat_clean, yr_strat_effort, by = c("est_year", "season", "stratum"))
    
    
    
    
    ####  4. Derived Stratum Area Estimates ####
    
    # a. Catch / tow, for that year & season
    survdat_clean <-  dplyr::mutate(
      survdat_clean,
                                      
      # Abundance
      abund_tow   = numlen_adj / strat_ntows,
      
      # All size biomass
      # Biomass is repeated across length classes at each station by species
      # the number of length classes is tallied where the conversion factor is done
      biom_per_lclass = (biomass_kg / n_len_class),
      biom_tow = biom_per_lclass / strat_ntows)
    
    
    # b. Stratified Mean Catch Rates
    survdat_clean <-  dplyr::mutate(
      survdat_clean,
                                      
      # Stratified mean abundance CPUE, weighted by the stratum areas
      strat_mean_abund = abund_tow * st_ratio,
      
      # Stratified mean BIOMASS CPUE
      strat_mean_bio = biom_tow * st_ratio)
    
    
    # c. Stratified Total Abundance/Biomass
    # convert from catch rate by area swept to total catch for entire stratum
    survdat_clean <-  dplyr::mutate(
      survdat_clean,
                                      
      # Total Abundance
      strat_total_abund = round((strat_mean_abund * tot_s_area / alb_tow_km2) / q),
      
      # Total BIOMASS from the biomass of all lengths
      strat_total_bio = (strat_mean_bio * tot_s_area / alb_tow_km2) / q)
  
}





####  Add Area Stratifications  ####
clean_stratified <- survdat_clean %>% 
  filter(season %in% season_opts,
         est_year <= 2019) %>% 
  area_stratify_station_totals(box_location = "cloudstorage")





#### Add the Lobster Management Zones  ####

# So from the ASMFC 2020 PDF, there is some modified spatial areas that go along with lobster management breaks

# We will add these in by stratum
# Shapefiles for the fisheries stat zones
# Resource Path
library(sf)
res_path <- gmRi::cs_path("res")
stat_zones <- read_sf(str_c(res_path, "Shapefiles/Statistical_Areas/Statistical_Areas_2010_withNames.shp"))

# Get the paths to the shapefiles used to mask
trawl_paths <- gmRi::get_timeseries_paths(region_group = "nmfs_trawl_regions", box_location = "cloudstorage")

# Polygons for each region
stat_zones <- read_sf(trawl_paths[["inuse_strata"]][["shape_path"]]) 

# Make a list of zones to roughly match the survey areas:
lobster_zones <- list(
  "Georges Bank / Gulf of Maine" = c(
    as.character(13:23),  # GB
    as.character(24:40), # GoM
    str_pad(as.character(9:12), width = 2, pad = "0", side = "left")),
  "Southern New England / Mid-Atlantic"  = 
    c(stringr::str_pad(as.character(1:8),width = 2, pad = "0", side = "left"),
      as.character(61:76))
)


# Trim out what we don't need and label them
lob_stat_zones <- stat_zones %>% 
  mutate(
    lobster_area = case_when(
      strat_num %in% lobster_zones$"Georges Bank / Gulf of Maine" ~ "Georges Bank / Gulf of Maine",
      strat_num %in% lobster_zones$"Southern New England / Mid-Atlantic" ~ "Southern New England / Mid-Atlantic")
    ) %>% 
  filter(survey_area %in% c(
    "Georges Bank", "Gulf of Maine", "Southern New England", "Mid-Atlantic Bight"))


# Map Check:
ggplot() +
  #geom_sf(data = us_poly) +
  #geom_sf(data = canada) +
  geom_sf(data = lob_stat_zones, aes(fill = lobster_area), alpha = 0.8) +
  coord_sf(xlim = c(-76.4, -64.4), ylim = c(35, 45.5), expand = F) +
  scale_fill_gmri() +
  theme_bw() +
  map_theme(legend.position = c(0.6, 0.125), 
            legend.background = element_rect(fill = "white"),
            plot.caption = element_text(hjust = 0)) +
  guides(fill = guide_legend(title = "", nrow = 2)) 



####  Aggregate Single Species Stratified Totals  ####


# Aggregate the stratified totals using the per-lenclass values
# Divide by a million
single_species_stratified <- clean_stratified %>%
  filter(comname == species_name) %>% 
  mutate(
    lobster_areas = case_when(
      strat_num %in% lobster_zones$"Georges Bank / Gulf of Maine" ~ "Georges Bank / Gulf of Maine",
      strat_num %in% lobster_zones$"Southern New England / Mid-Atlantic" ~ "Southern New England / Mid-Atlantic")
    ) %>% 
  group_by(comname, est_year, lobster_areas) %>% 
  summarise(
     area_stratified_abundance = sum(strat_total_abund, na.rm = T)/1e6,
     area_stratified_biomass = sum(strat_total_bio, na.rm = T)/1e6,
     .groups = "drop")  %>% 
  pivot_longer(cols = starts_with("area_stratified"),
               names_to = "metric",
               values_to = "annual_total")

# Plot the Annual Totals
(annual_stratified_plot <- single_species_stratified  %>% 
  ggplot(aes(est_year, annual_total)) +
  geom_line() +
  facet_grid(lobster_areas~metric, scales = "free") +
  scale_y_continuous(labels = label_comma(suffix = "m")) +
  labs(x = "Year", y = "Annual Total", title = species_name))
