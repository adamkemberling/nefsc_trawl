####
#### NEFSC Trawl Data - Size Spectra Build
#### 9/28/2020
####
#### Objective: 
#### Load 2020 "survdat" data from Janet Nye, perform all filtering and adjustment for Size Spectra Analyses
####




#### Packages  ####
library(janitor)
library(magrittr)
library(here)
library(gmRi)
library(tidyverse)



####____________________####

#### 2019-2020 Size Spectrum Build  ####
load_ss_data <- function(survdat = NULL, survdat_source = "2020"){

  ####  Paths  ####
  mills_path <- shared.path(group = "Mills Lab", folder = "")
  nsf_path <- shared.path(group = "NSF OKN", folder = "")
  res_path   <- shared.path(os.use = "unix", group = "RES Data", folder = NULL)
  
  # Support Functions
  `%not in%` <-  purrr::negate(`%in%`) # to filter out ones that have matches
  
  
  ####  Load Data  ####
  
  # If providing a starting point pass it in:
  survdat <- survdat 
  
  # 2020 groundfish survey data
  if(survdat_source == "2020"){
    load(here("data/NEFSC/Survdat_Nye_Aug 2020.RData"))
    survdat <- clean_names(survdat)
  } 
  
  
  # Survey Stratumm Areas
  stratum_area <- read_csv(str_c(mills_path, "Projects/NSF_CAccel/Data/strata area.csv"), 
                           col_types = cols())
  
  
  
  ####  2020 Survey Data Prep  ####
  
  
  ####__ 1. Column Changes  ####
  trawldat <- survdat %>% 
    mutate(comname = tolower(comname),
           id = format(id, scientific = FALSE),
           biomass = ifelse(is.na(biomass) == TRUE & abundance > 0, 0.0001, biomass),
           biomass = ifelse(biomass == 0 & abundance > 0, 0.01, biomass),
           abundance = ifelse(is.na(abundance) == TRUE & biomass > 0, 1, abundance),
           abundance = ifelse(abundance == 0 & biomass > 0, 1, abundance),
           tstrat = str_sub(stratum, 2, 3)) 
  
  
  ####__ 2. Column Selection  ####

  trawldat <- trawldat %>%
    select(id, est_year, est_month, est_day, season, stratum, tstrat, decdeg_beglat, decdeg_beglon,
           svspp, comname, catchsex, biomass, avgdepth, abundance, length, numlen) 
  
  
  
  ####__ 3. Filtering  ####
  trawldat <- trawldat %>%
    filter(
      # Eliminate Candian Strata and Not in Use Strata
      stratum >= 01010,
      stratum <= 01760,
      stratum != 1310,
      stratum != 1320,
      stratum != 1330,
      stratum != 1350,
      stratum != 1410,
      stratum != 1420,
      stratum != 1490,
      # Filter to just Spring and Fall
      season %in% c("SPRING", "FALL"),
      est_year >= 1970,
      est_year < 2020,
      # Drop na Biomass and Abundance Records
      !is.na(biomass),
      !is.na(abundance))
  
  
  
  ####__ 4. Stratum Area Ratios  ####
  trawldat <- left_join(trawldat, stratum_area, by = "stratum")
  trawldat <- arrange(trawldat, id)
  
  
  # Compute Stratum Area Ratio & Adjusted Abundances
  total_stratum_area <- distinct(stratum_area, stratum, .keep_all = T) %$% sum(area, na.rm = T)
  trawldat <- trawldat %>% 
    mutate(stratio = area / total_stratum_area,
           biom_adj = ifelse(biomass == 0 & abundance > 0, 0.0001, biomass),
           abund_adj = ifelse(abundance == 0 & biomass > 0, 1, abundance),
           indwt = biom_adj / abund_adj)
  
  
  
  ###__ 5. Explicit Species Exclusion  ####
  trawldat <- trawldat %>% 
    filter(
      # Exclude the Shrimps
      svspp %not in% c(285:299, 305, 306, 307, 316, 323, 910:915, 955:961),
      # Exclude the unidentified fish
      svspp %not in% c(0, 978, 979, 980, 998)
    )
  
  
  ####__ 6. Spatial Filtering  ####
  
  # If we know which strata represent which area we can use this:
  # Stratum Key for filtering specific areas
  strata_key <- list(
    "Georges Bank"          = as.character(13:23),
    "Gulf of Maine"         = as.character(24:40),
    "Southern New England"  = str_pad(as.character(1:12), width = 2, pad = "0", side = "left"),
    "Mid-Atlantic Bight"    = as.character(61:76))
  
  # Add labels to the data
  trawldat <- trawldat %>% 
    mutate(area =  case_when(
      tstrat %in% strata_key$`Georges Bank`         ~ "GB", 
      tstrat %in% strata_key$`Gulf of Maine`        ~ "GoM",
      tstrat %in% strata_key$`Southern New England` ~ "SNE",
      tstrat %in% strata_key$`Mid-Atlantic Bight`   ~ "MAB",
      TRUE                                          ~ "not found"))
  
  # Use strata_select to pull the strata we want
  strata_select <- c(
    strata_key$`Georges Bank`, 
    strata_key$`Gulf of Maine`,
    strata_key$`Southern New England`,
    strata_key$`Mid-Atlantic Bight`)
  
  
  # Filtering
  trawldat <- trawldat %>% filter(tstrat %in% strata_select)
  
  
  
  
  ####__ 7. Adjusted NumLength  ####
  conv_factor <- trawldat %>% 
    group_by(id, comname, catchsex, numlen, abundance) %>% 
    summarise(abundance_raw = sum(numlen),         # Total abundance by numlen aka how many measured at that size
              #abundance = mean(abundance),         # abundance is number measured of that species, same across the groups so mean pulls one value
              convers =  abundance/abundance_raw)  # convers is the difference between abundance and the number measured
  
  
  # Merge back and convert the numlen field
  trawldat <- trawldat %>% 
    left_join(conv_factor, by = c("id", "comname", "catchsex", "numlen", "abundance")) %>% 
    mutate(numlen_adj = numlen * convers) 
  
  # remove conversion factor from environment
  rm(conv_factor, strata_key, strata_select, stratum_area)
  
  
  
  #### Length to Bodymass Conversions  ####
  
  
  
  ####__ 1. Distinct Station & Species Length Info   ####
  
  # For each station we need unique combinations of
  # station_id, species, catchsex, length, adjusted_numlen
  # Record of unique station catches: # rows for each species * sex * length
  trawl_lens <- trawldat %>% 
    filter(is.na(length) == FALSE,
           is.na(numlen) == FALSE,
           numlen > 0) %>% 
    distinct(id, svspp, comname, catchsex, length, numlen, numlen_adj, biom_adj) %>% 
    mutate(svspp = as.character(svspp),
           svspp = str_pad(svspp, 3, "left", "0"))
    
  
  
  # Pull distinct records of the stations themselves and metadata that match
  trawl_stations <- trawldat %>% 
    select(id, est_year, est_month,  est_day, season, stratum, stratum_num, area, decdeg_beglat, decdeg_beglon, 
           type, avgdepth, stratio) %>% 
    distinct() 
  
  # recombine with the distinct station info
  trawl_spectra <- trawl_stations %>% 
    left_join(trawl_lens, by = "id")
  
  
  
  ####__ 2. Combine with Growth Coefficients  ####
  
  #Now we want to use the lw_combined here instead of just the fishbase lengths
  lw_combined <- read_csv(here::here("data/biomass_key_combined.csv"),
                          col_types = cols()) %>% 
    mutate(svspp = as.character(svspp),
           svspp = str_pad(svspp, 3, "left", "0"))
  
  # Do a priority pass with the filter(lw_combined, source == "wigley)
  # merge on comname, season, and catchsex
  w_trimmed <- filter(lw_combined, source == "wigley") %>% 
    select(source, season, svspp, comname, scientific_name, spec_class, 
           hare_group, fishery, catchsex, a, b, ln_a)
  
  
  # Do a second pass with the filter(lw_combined, source == "fishbase")
  # merge on common names only
  fb_trimmed <- filter(lw_combined, source == "fishbase") %>% 
    select(source, svspp, comname, scientific_name, spec_class, 
           hare_group, fishery, a, b, ln_a)
  
  # First Pass - Wigley
  pass_1 <- trawl_spectra %>% 
    # Join just by svspp to account for name changes
    select(-comname) %>% 
    inner_join(w_trimmed)
  
  
  # Second Pass - Fishbase, for the stragglers if any
  # currently not joining well because lack of svspp and comname agreement 
  pass_2 <- trawl_spectra %>% 
    filter(comname %not in% w_trimmed$comname) %>% 
    inner_join(fb_trimmed)
  
  
  
  # Join them with bind rows (implicitly drops things that don't have growth coefs)
  trawl_weights <- bind_rows(pass_1, pass_2) %>% 
    arrange(est_year, season) %>% 
    mutate(b = as.numeric(b),
           a = as.numeric(a),
           a = ifelse(is.na(a) & !is.na(ln_a), exp(ln_a), a),
           ln_a = ifelse(is.na(ln_a), log(a), ln_a),          # log of a used if ln_a is isn't already there (some fish just had ln_a reported)
           llen = log(length),
           ind_log_wt = ln_a + (b * llen),
           ind_weight_kg = exp(ind_log_wt),                    # weight of an individual in size class
           sum_weight_kg = ind_weight_kg * numlen_adj) %>%      # Individual weight * adjusted numlen
    drop_na(ind_weight_kg) %>% 
    select(-ind_log_wt, -llen)
  
  
  # clean up environment
  rm(pass_1, pass_2, fb_trimmed, w_trimmed)
  
  
  
  
  ####__ 3. Use Coefficients to Re-calculate Biomass  ####
  
  # calculate total biomass again using weights from key 
  # make a key for the length weight coefficient sources
  trawl_weights <- trawl_weights %>%  
    arrange(est_year, season, comname, length) %>% 
    mutate(lw_group = str_c(comname, season, catchsex)) %>% 
    as_tibble()
  
  ####__ 4. Area Weighted Biomass  ####
  trawl_weights <- trawl_weights %>% 
    mutate(
      strat_wt_bio_fscs = stratio * biom_adj,
      strat_wt_bio_lw = stratio * sum_weight_kg
    )
  
  return(trawl_weights)
}



####____________________####

####____________________####
####  2016 Size Spectrum Build  ####

load_2016_ss_data <- function(){
  
  
  ####  Paths  ####
  mills_path <- shared.path(group = "Mills Lab", folder = "")
  nsf_path <- shared.path(group = "NSF OKN", folder = "")
  res_path   <- shared.path(os.use = "unix", group = "RES Data", folder = NULL)
  
  
  #  Support Functions  
  `%not in%` <-  purrr::negate(`%in%`) 
  
  
  ####  Load Data  ####
  load(str_c(mills_path, "Projects/WARMEM/Old survey data/Survdat.RData"))
  survdat <- clean_names(survdat) %>% 
    mutate(svspp = as.character(svspp),
           svspp = str_pad(svspp, 3, "left", "0"))
  
  
  # Extra prep for 2016 data, getting species names via svspp codes
  spp_classes <- read_csv(here("data/kmills/sppclass.csv"),
                          col_types = cols()) %>% 
    clean_names() %>% 
    mutate(common_name = str_to_lower(common_name),
           scientific_name = str_to_lower(scientific_name)) %>% 
    distinct(svspp, comname = common_name, scientific_name)
  
  
  # Add the common names over and format for rest of build
  survdat <- left_join(survdat, clean_names(spp_classes)) %>% 
    drop_na(comname) %>% 
    mutate(cruise6 = str_pad(cruise6, 6, "left", "0"),
           station = str_pad(station, 3, "left", "0"),
           stratum = str_pad(stratum, 4, "left", "0"),
           id = str_c(cruise6, station, stratum)) %>% 
    select(id, est_year = year, station, stratum, svvessel, season, lat, lon, depth, 
           surftemp, surfsalin, bottemp, botsalin, svspp, comname, scientific_name, everything()) %>% 
    mutate(comname = str_to_lower(comname),
           scientific_name = str_to_lower(scientific_name))
  
  
  
  
  ####  Data Prep  ####
  
  
  
  
  
  ####__ 1. Column Changes  ####
  trawldat <- survdat %>% 
    mutate(comname = tolower(comname),
           id = format(id, scientific = FALSE),
           biomass = ifelse(is.na(biomass) == TRUE & abundance > 0, 0.01, biomass),
           biomass = ifelse(biomass == 0 & abundance > 0, 0.01, biomass),
           abundance = ifelse(is.na(abundance) == TRUE & biomass > 0, 1, abundance),
           abundance = ifelse(abundance == 0 & biomass > 0, 1, abundance),
           tstrat = str_sub(stratum, 2, 3)) 
  
  rm(survdat)
  
  ####__ 2. Column Selection  ####
  trawldat <- trawldat %>%
    select(id, est_year,season, stratum, tstrat, decdeg_beglat = lat, decdeg_beglon = lon,
           svspp, comname, catchsex, biomass, avgdepth = depth, abundance, length, numlen) 
  
  
  
  ####__ 3. Filtering  ####
  trawldat <- trawldat %>%
    filter(
      # Eliminate Candian Strata and Not in Use Strata
      stratum >= 01010,
      stratum <= 01760,
      stratum != 1310,
      stratum != 1320,
      stratum != 1330,
      stratum != 1350,
      stratum != 1410,
      stratum != 1420,
      stratum != 1490,
      # Filter to just Spring and Fall
      season %in% c("SPRING", "FALL"),
      est_year >= 1970,
      est_year < 2020,
      # Drop na Biomass and Abundance Records
      !is.na(biomass),
      !is.na(abundance))
  
  
  
  ####__ 4. Stratum Area Ratios  ####
  # Survey Stratumm Areas
  stratum_area <- read_csv(str_c(mills_path, "Projects/NSF_CAccel/Data/strata area.csv"), 
                           col_types = cols()) %>% 
    mutate(stratum = as.character(stratum))
  
  trawldat <- left_join(trawldat, stratum_area, by = "stratum")
  trawldat <- arrange(trawldat, id)
  
  
  # Compute Stratum Area Ratio & Adjusted Abundances
  total_stratum_area <- distinct(stratum_area, stratum, .keep_all = T) %$% sum(area, na.rm = T)
  trawldat <- trawldat %>% 
    mutate(stratio = area / total_stratum_area,
           biom_adj = ifelse(biomass == 0 & abundance > 0, 0.0001, biomass),
           abund_adj = ifelse(abundance == 0 & biomass > 0, 1, abundance),
           indwt = biom_adj / abund_adj)
  
  
  
  ###__ 5. Explicit Species Exclusion  ####
  trawldat <- trawldat %>% 
    filter(
      # Exclude the Shrimps
      svspp %not in% c(285:299, 305, 306, 307, 316, 323, 910:915, 955:961),
      # Exclude the unidentified fish
      svspp %not in% c(0, 978, 979, 980, 998)
    )
  
  
  ####__ 6. Spatial Filtering  ####
  
  # If we know which strata represent which area we can use this:
  # Stratum Key for filtering specific areas
  strata_key <- list(
    "Georges Bank"          = as.character(13:23),
    "Gulf of Maine"         = as.character(24:40),
    "Southern New England"  = str_pad(as.character(1:12), width = 2, pad = "0", side = "left"),
    "Mid-Atlantic Bight"    = as.character(61:76))
  
  # Add labels to the data
  trawldat <- trawldat %>% 
    mutate(area =  case_when(
      tstrat %in% strata_key$`Georges Bank`         ~ "GB", 
      tstrat %in% strata_key$`Gulf of Maine`        ~ "GoM",
      tstrat %in% strata_key$`Southern New England` ~ "SNE",
      tstrat %in% strata_key$`Mid-Atlantic Bight`   ~ "MAB",
      TRUE                                          ~ "not found"))
  
  # Use strata_select to pull the strata we want
  strata_select <- c(
    strata_key$`Georges Bank`, 
    strata_key$`Gulf of Maine`,
    strata_key$`Southern New England`,
    strata_key$`Mid-Atlantic Bight`)
  
  
  # Filtering
  trawldat <- trawldat %>% filter(tstrat %in% strata_select)
  
  
  
  
  ####__ 7. Adjusted NumLength  ####
  conv_factor <- trawldat %>% 
    group_by(id, comname, catchsex, numlen, abundance) %>% 
    summarise(abundance_raw = sum(numlen),         # Total abundance by numlen
              #abundance = mean(abundance),         # abundance is the same across the groups so mean pulls one value
              convers =  abundance/abundance_raw)  # convers is the difference between abundance and the number measured
  
  
  # Merge back and convert the numlen field
  trawldat <- trawldat %>% 
    left_join(conv_factor, by = c("id", "comname", "catchsex", "numlen", "abundance")) %>% 
    mutate(numlen_adj = numlen * convers) 
  
  # remove conversion factor from environment
  rm(conv_factor, strata_key, strata_select, stratum_area)
  
  
  
  
  #### Length to Bodymass Conversions  ####
  
  
  
  ####__ 1. Distinct Station & Species Length Info   ####
  
  
  
  #-#
  trawl_lens <- trawldat %>% 
    filter(is.na(numlen) == FALSE,
           numlen > 0) %>% 
    distinct(id, svspp, comname, catchsex, length, numlen, numlen_adj, biom_adj) %>% 
    mutate(svspp = as.character(svspp),
           svspp = str_pad(svspp, 3, "left", "0"))
  
  
  
  # Pull distinct records of the stations themselves and metadata that match
  trawl_stations <- trawldat %>% 
    select(id, est_year, season, stratum, stratum_num, area, decdeg_beglat, decdeg_beglon, 
           type, avgdepth, stratio) %>% 
    distinct() 
  
  # recombine with the distinct station info
  trawl_spectra <- trawl_stations %>% 
    left_join(trawl_lens, by = "id")
  
  
  
  ####__ 2. Combine with Growth Coefficients  ####
  
  #Now we want to use the lw_combined here instead of just the fishbase lengths
  lw_combined <- read_csv(here::here("data/biomass_key_combined.csv"),
                          col_types = cols()) %>% 
    mutate(svspp = as.character(svspp),
           svspp = str_pad(svspp, 3, "left", "0"))
  
  # Do a priority pass with the filter(lw_combined, source == "wigley)
  # merge on comname, season, and catchsex
  w_trimmed <- filter(lw_combined, source == "wigley") %>% 
    select(source, season, svspp, comname, scientific_name, spec_class, 
           hare_group, fishery, catchsex, a, b, ln_a)
  
  
  # Do a second pass with the filter(lw_combined, source == "fishbase")
  # merge on common names only
  fb_trimmed <- filter(lw_combined, source == "fishbase") %>% 
    select(source, svspp, comname, scientific_name, spec_class, 
           hare_group, fishery, a, b, ln_a)
  #-#
  
  
  
  # First Pass - Wigley
  pass_1 <- trawl_spectra %>% 
    # Join just by svspp to account for name changes
    select(-comname) %>% 
    inner_join(w_trimmed)
  
  
  # Second Pass - Fishbase, for the stragglers if any
  # currently not joining well because lack of svspp and comname agreement 
  pass_2 <- trawl_spectra %>% 
    filter(comname %not in% w_trimmed$comname) %>% 
    inner_join(fb_trimmed)
  

  
  
  
  # Join them with bind rows (implicitly drops things that don't have growth coefs)
  trawl_weights <- bind_rows(pass_1, pass_2) %>% 
    arrange(est_year, season) %>% 
    mutate(b = as.numeric(b),
           a = as.numeric(a),
           a = ifelse(is.na(a) & !is.na(ln_a), exp(ln_a), a),
           ln_a = ifelse(is.na(ln_a), log(a), ln_a),
           ind_log_wt = ln_a + (b * log(length)),
           ind_weight_kg = exp(ind_log_wt),                    # weight of an individual in size class
           sum_weight_kg = ind_weight_kg * numlen_adj) %>%      # Individual weight * adjusted numlen
    drop_na(ind_weight_kg) %>% 
    select(-ind_log_wt)
  
  # clean up environment
  rm(pass_1, pass_2, fb_trimmed, w_trimmed)
  
  
  
  ####__ 3. Use Coefficients to Re-calculate Biomass  ####
  
  # calculate total biomass again using weights from key 
  # make a key for the length weight coefficient sources
  trawl_weights <- trawl_weights %>%  
    arrange(est_year, season, comname, length) %>% 
    mutate(lw_group = str_c(comname, season, catchsex)) %>% 
    as_tibble()
  
  
  rm(trawl_spectra, trawl_stations, trawl_lens)
  
  
  
  ####__ 4. Area Weighted Biomass  ####
  trawl_weights <- trawl_weights %>% 
    mutate(
      strat_wt_bio_fscs = stratio * biom_adj,
      strat_wt_bio_lw = stratio * sum_weight_kg
    )
  
  
  # Pass out 2016 weights
  return(trawl_weights)
  
  
  
}



####____________________####
####____________________####
####  Stratified Catch/Biomass  ####
####  AND  ####
###  EDA Functions ####

#### Function to get stratified means at various group levels
# stratified = metric / num tows
ss_stratification <- function(.data = data, ...){
  
  # number of effort (tows) in the groups
  num_tows <- .data %>% 
    group_by(...) %>% 
    summarise(tow_count = n())
  
  # Add effort to the data
  data_stratified <- .data %>% 
    left_join(num_tows)
  
  # now get the stratified values at each fish growth group
  data_stratified <- data_stratified  %>% 
    group_by(..., comname, lw_group) %>% 
    summarise(
      ## Things that are the same across groups ##
      a = mean(a),                   
      b = mean(b),
      ln_a = mean(ln_a),
      ind_weight_kg = mean(ind_weight_kg),
      ## Things that we want to total for the stratum  ##
      numlen_adj = sum(numlen_adj),       # Total count of individuals
      sum_weight_kg = sum(sum_weight_kg),   # Total Biomass of species at length for group strata
      mbiom_per_tow = sum_weight_kg / mean(tow_count),  # total weight averaged by n-tows
      wtmean_per_tow = mbiom_per_tow * mean(stratio))  # relative to stratum area
  
  # Recreate the catchsex bit
  data_stratified <- data_stratified %>% 
    mutate(catchsex = case_when(
      str_detect(lw_group, "0") ~ "0",
      str_detect(lw_group, "1") ~ "1")) %>% 
    ungroup()
  
  # Return the data aggreagted to group levels
  return(data_stratified)
  
}

# # Test it out
# # Just years and stratum
# strat_test <- ss_stratification(data, est_year, stratum_num)
# glimpse(strat_test)
# 
# # Year as well as season and stratum
# t2 <- ss_stratification(data, est_year, season, stratum_num)
# glimpse(t2)
# 
# # Individual tows
# tow_id_test <- ss_stratification()




####__####



# Return an annual summary
ss_annual_summary <- function(survey_data) {
  
  # Summarize on a measurement length level
  summ_data <- survey_data %>% 
    group_by(est_year) %>% 
    summarise(
      area = "All Areas",
      season = "Spring + Fall",
      n_stations = n_distinct(id),
      n_species = n_distinct(comname),
      lw_biomass_kg = sum(sum_weight_kg),
      lw_biomass_per_station = lw_biomass_kg / n_stations,
      total_fish = sum(numlen_adj),
      mean_ind_length = weighted.mean(length, numlen_adj),
      mean_ind_bodymass = weighted.mean(ind_weight_kg, numlen_adj),
      lw_strat_biomass = sum(strat_wt_bio_lw)
    )
  
  # Biomass from aggregate catch by species at station, 
  #these are distinct to each tow but repeated for different sizes of the same species
  agg_biomass <- survey_data %>% 
    group_by(est_year, id, comname,  biom_adj) %>% 
    summarise(
      haul_biomass = biom_adj,
      haul_strat_biomass = mean(strat_wt_bio_fscs, na.rm = T)) %>% 
    ungroup() %>% 
    group_by(est_year) %>% 
    summarise(
      fscs_biomass = sum(haul_biomass),
      fscs_strat_biomass = sum(haul_strat_biomass))
  
  # rejoin
  summ_data <- left_join(summ_data, agg_biomass, by = c("est_year")) %>% 
    mutate(
      fscs_biom_per_station = fscs_biomass / n_stations,
      mean_fscs_bio = fscs_biomass / total_fish)
  return(summ_data)
  
}




####__####



# Include region as a group
ss_regional_differences <- function(survey_data) {
  
  # Summarize on a measurement length level
  summ_data <- survey_data %>% 
    group_by(est_year, area) %>% 
    summarise(
      season = "Spring + Fall",
      n_stations = n_distinct(id),
      n_species = n_distinct(comname),
      lw_biomass_kg = sum(sum_weight_kg),
      lw_biomass_per_station = lw_biomass_kg / n_stations,
      total_fish = sum(numlen_adj),
      mean_ind_length = weighted.mean(length, numlen_adj),
      mean_ind_bodymass = weighted.mean(ind_weight_kg, numlen_adj),
      lw_strat_biomass = sum(strat_wt_bio_lw)
  )
  
  # Biomass from aggregate catch by species at station, 
  #these are distinct to each tow but repeated for different sizes of the same species
  agg_biomass <- survey_data %>% 
    group_by(est_year, area, id, comname, biom_adj) %>% 
    summarise(
      haul_biomass = biom_adj,
      haul_strat_biomass = mean(strat_wt_bio_fscs, na.rm = T)) %>% 
    ungroup() %>% 
    group_by(est_year, area) %>% 
    summarise(
      fscs_biomass = sum(haul_biomass),
      fscs_strat_biomass = sum(haul_strat_biomass))
  
  # rejoin
  summ_data <- left_join(summ_data, agg_biomass, by = c("est_year", "area")) %>% 
    mutate(
      fscs_biom_per_station = fscs_biomass / n_stations,
      mean_fscs_bio = fscs_biomass / total_fish)
  return(summ_data)
  
}



####__####

# Include season as a group
ss_seasonal_summary <- function(survey_data){
  
  # Summarize on a measurement length level
  summ_data <- survey_data %>% 
    group_by(est_year, season) %>% 
    summarise(
      area = "All Areas",
      n_stations = n_distinct(id),
      n_species = n_distinct(comname),
      lw_biomass_kg = sum(sum_weight_kg),
      lw_biomass_per_station = lw_biomass_kg / n_stations,
      total_fish = sum(numlen_adj),
      mean_ind_length = weighted.mean(length, numlen_adj),
      mean_ind_bodymass = weighted.mean(ind_weight_kg, numlen_adj),
      lw_strat_biomass = sum(strat_wt_bio_lw)
    ) 
  
  # Biomass from aggregate catch by species at station, 
  #these are distinct to each tow but repeated for different sizes of the same species
  agg_biomass <- survey_data %>% 
    group_by(est_year, season, id, comname, biom_adj) %>% 
    summarise(
      haul_biomass = biom_adj,
      haul_strat_biomass = mean(strat_wt_bio_fscs, na.rm = T)) %>% 
    ungroup() %>% 
    group_by(est_year, season) %>% 
    summarise(
      fscs_biomass = sum(haul_biomass),
      fscs_strat_biomass = sum(haul_strat_biomass))
  
  # rejoin
  summ_data <- left_join(summ_data, agg_biomass, by = c("est_year", "season")) %>% 
    mutate(
      fscs_biom_per_station = fscs_biomass / n_stations,
      mean_fscs_bio = fscs_biomass / total_fish)
  return(summ_data)
  
}




####  Testing merge by svspp  ####

# look for #-# for insert

# For each station we need unique combinations of
# station_id, species, catchsex, length, adjusted_numlen
# Record of unique station catches: # rows for each species * sex * length
# trawl_lens <- trawldat %>% 
#   filter(is.na(numlen) == FALSE,
#          numlen > 0) %>% 
#   distinct(id, comname, catchsex, length, numlen, numlen_adj, biom_adj) 
# 
# 
# # Pull distinct records of the stations themselves and metadata that match
# trawl_stations <- trawldat %>% 
#   select(id, est_year, season, stratum, stratum_num, area, decdeg_beglat, decdeg_beglon, 
#          type, avgdepth, stratio) %>% 
#   distinct() 
# 
# 
# # drop unnecessary space
# rm(trawldat)
# 
# 
# # recombine with the distinct station info
# trawl_spectra <- trawl_stations %>% 
#   left_join(trawl_lens, by = "id")
# 
# 
# 
# 
# ####__ 2. Combine with Growth Coefficients  ####
# 
# #Now we want to use the lw_combined here instead of just the fishbase lengths
# lw_combined <- read_csv(here::here("data/biomass_key_combined.csv"),
#                         col_types = cols())
# 
# # Do a priority pass with the filter(lw_combined, source == "wigley)
# # merge on comname, season, and catchsex
# w_trimmed <- filter(lw_combined, source == "wigley") %>% 
#   select(source, season, comname, scientific_name, spec_class, 
#          hare_group, fishery, catchsex, a, b, ln_a)
# 
# 
# # Do a second pass with the filter(lw_combined, source == "fishbase")
# # merge on common names only
# fb_trimmed <- filter(lw_combined, source == "fishbase") %>% 
#   select(source, comname, scientific_name, spec_class, hare_group, fishery, a, b, ln_a)
