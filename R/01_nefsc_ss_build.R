####
#### NEFSC Trawl Data - Size Spectra Build
#### 9/28/2020
####
#### Objective: 
#### Load 2020 "survdat" data from Janet Nye, perform all filtering and adjustment for Size Spectra Analyses
####




####  Load Packages  ####
library(janitor)
library(magrittr)
library(here)
library(gmRi)
library(sf)
library(tidyverse)



####____________________####
#### Prepping EPU area file  ####
# This only needs to be done once
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
# trawl_strata <- survey_strata <- read_sf(str_c(res_path, "Shapefiles/BottomTrawlStrata/BTS_Strata.shp"))
# trawl_strata <- trawl_strata %>%  clean_names()
# trawl_sp     <- sf::as_Spatial(trawl_strata)
# # get area
# strata_areas <- getarea(trawl_sp, strat.col = "strata") %>% setNames(c("stratum", "s_area_km2"))
# write_csv(strata_areas, str_c(mills_path, "Projects/NSF_CAccel/Data/strata_areas_km2.csv"))


####____________________####

#### 2019-2020 Size Spectrum Build  ####
#' @title  Survdat Size Spectrum Build
#'
#' @description Processing function to prepare survdat data for size spectra analyses. 
#' 
#'
#' @param survdat optional candidate dataframe in the R environment to run through size spectra build.
#' @param survdat_source String indicating which survdat file to load from box
#'
#' @return Returns a dataframe filtered and tidy-ed for size spectrum analysis.
#' @export
#'
#' @examples
load_ss_data <- function(survdat = NULL, survdat_source = "2020"){

  ####  Paths  
  mills_path <- shared.path(os.use = "unix", group = "Mills Lab", folder = "")
  nsf_path   <- shared.path(os.use = "unix", group = "NSF OKN", folder = "")
  res_path   <- shared.path(os.use = "unix", group = "RES Data", folder = NULL)
  
  
  ####  Load Data  ####
  
  # If providing a starting point pass it in:
  if(exists("survdat")){ survdat <- survdat }
  
  # 2020 groundfish survey data
  if(survdat_source == "2020"){
    load(here("data/NEFSC/Survdat_Nye_Aug 2020.RData"))
    survdat <- clean_names(survdat)
  } 
  
  
  # # Load Survey Stratum Area Info - old and in nautical miles (with errors)
  # stratum_area <- read_csv(str_c(mills_path, "Projects/NSF_CAccel/Data/strata area.csv"), col_types = cols())
  
  
  # Load fresh stratum areas
  stratum_area <- read_csv(str_c(mills_path, "Projects/NSF_CAccel/Data/strata_areas_km2.csv"), col_types = cols())
  
  # EPU are info: source slucey_survdat_functions.R and {ecodata}
  epu_areas <- read_csv(str_c(mills_path, "Projects/NSF_CAccel/Data/EPU_areas_km2.csv"), col_types = cols())
  
  
  
  ####__ 1. Column Changes  ####
  trawldat <- survdat %>% 
    mutate(
      comname   = tolower(comname),
      id        = format(id, scientific = FALSE),
      biomass   = ifelse(is.na(biomass) == TRUE & abundance > 0, 0.0001, biomass),
      abundance = ifelse(is.na(abundance) == TRUE & biomass > 0, 1, abundance),
      #stratum number excluding leading and trailing codes for inshore/offshore
      strat_num    = str_sub(stratum, 2, 3)) %>%  
    mutate(biom_adj  = ifelse(biomass == 0 & abundance > 0, 0.0001, biomass), .after = biomass) %>% 
    mutate(abund_adj = ifelse(abundance == 0 & biomass > 0, 1, abundance), .after = abundance)
  
  ####__ 2. Column Selection ####

  # currently a light-weight group of columns, 
  # leaves behind CTD and shipboard instrument details.
  # Favors the larger categorical group metadata
  trawldat <- trawldat %>%
    select(
      id, 
      cruise6,
      station,
      est_year, 
      est_month, 
      est_day, 
      svvessel, 
      season, 
      stratum, 
      strat_num, 
      decdeg_beglat, 
      decdeg_beglon,
      avgdepth, 
      svspp, 
      comname, 
      catchsex, 
      length, 
      numlen,
      abundance, 
      abund_adj,
      biomass, 
      biom_adj
    ) 
  
  
  
  ####__ 3. Data Filtering  ####
  # 1. Strata
  # 2. Seasons
  # 3. Year limits
  # 4. Vessels
  # 5. Species Exclusion
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
      # Only the Albatross and Henry Bigelow
      svvessel %in% c("AL", "HB"),
      est_year >= 1970,
      est_year < 2020,
      # Drop na Biomass and Abundance Records
      !is.na(biomass),
      !is.na(abundance),
      # Exclude the Skrimps
      svspp %not in% c(285:299, 305, 306, 307, 316, 323, 910:915, 955:961),
      # Exclude the unidentified fish
      svspp %not in% c(0, 978, 979, 980, 998)
    )
  
  
  # ####__ 6. Spatial Filtering  ####
  
  
  #### EPU assignment for survdat stations - function from Sean Luceys "RSurvey" repo
  source(here("R/kathy_ss_code/slucey_survdat_functions.R"))
  epu_sf <- ecodata::epu_sf
  epu_sp <- as_Spatial(epu_sf)
  
  
  # assign EPU labels
  trawldat <- trawldat %>% 
    rename(CRUISE6 = cruise6,
           STATION = station,
           STRATUM = stratum,
           LAT     = decdeg_beglat,
           LON     = decdeg_beglon) %>% 
    poststrat(survdat = ., stratum = epu_sp, strata.col = "EPU") %>% 
    rename(epu           = newstrata,
           cruise6       = CRUISE6,
           station       = STATION,
           stratum       = STRATUM,
           decdeg_beglat = LAT,
           decdeg_beglon = LON)
    
  
  # Stratum Key for filtering specific areas
  strata_key <- list(
    "Georges Bank"          = as.character(13:23),
    "Gulf of Maine"         = as.character(24:40),
    "Southern New England"  = str_pad(as.character(1:12), width = 2, pad = "0", side = "left"),
    "Mid-Atlantic Bight"    = as.character(61:76))

  # Add labels to the data
  trawldat <- trawldat %>%
    mutate(
      survey_area =  case_when(
      strat_num %in% strata_key$`Georges Bank`         ~ "GB",
      strat_num %in% strata_key$`Gulf of Maine`        ~ "GoM",
      strat_num %in% strata_key$`Southern New England` ~ "SNE",
      strat_num %in% strata_key$`Mid-Atlantic Bight`   ~ "MAB",
      TRUE                                             ~ "not found"))

  
  # Optional, Use strata_select to pull the strata we want individually
  strata_select <- c(
    strata_key$`Georges Bank`,
    strata_key$`Gulf of Maine`,
    strata_key$`Southern New England`,
    strata_key$`Mid-Atlantic Bight`)


  # Filtering with strata_select
  trawldat <- trawldat %>% filter(strat_num %in% strata_select)
  
  
  
  
  ####__ 4. Stratum Area/Effort Ratios  ####
  
  
  
  # Join to the files containing area of each stratum, epu
  trawldat <- trawldat %>% 
    left_join(stratum_area, by = "stratum") %>% 
    left_join(epu_areas, by = "epu") %>% 
    arrange(trawldat, id)
  
  
  
  
  # Get Total stratum area of all strata 
  # (excludes ones we do not care about via left join)
  total_stratum_areas <- trawldat %>% 
    group_by(est_year) %>% 
    distinct(stratum, .keep_all = T) %>%  
    summarise(tot_s_area =  sum(s_area_km2, na.rm = T),
              .groups = "keep") 
  
  total_epu_areas <- trawldat %>% 
    group_by(est_year) %>% 
    distinct(epu, .keep_all = T) %>%  
    summarise(tot_epu_area =  sum(epu_area_km2, na.rm = T), 
              .groups = "keep")
  
  
  # Calculate strata area relative to total area i.e. stratio or stratum weights
  trawldat <- trawldat %>% 
    left_join(total_stratum_areas, by = "est_year") %>% 
    left_join(total_epu_areas, by = "est_year") %>% 
    mutate(st_ratio   = s_area_km2 / tot_s_area,
           epu_ratio  = epu_area_km2 / tot_epu_area) 
  
  
  # Get number of unique tows per stratum
  yr_strat_effort <- trawldat %>% 
    group_by(est_year, stratum) %>% 
    summarise(strat_ntows = n_distinct(id), 
              .groups = "keep")
  
  
  yr_epu_effort <-  trawldat %>% 
    group_by(est_year, epu) %>% 
    summarise(epu_ntows = n_distinct(id), 
              .groups = "keep")
  
  
  # Add those yearly effort counts back for later
  trawldat <- trawldat %>% 
    left_join(yr_strat_effort, by = c("est_year", "stratum")) %>% 
    left_join(yr_epu_effort, by = c("est_year", "epu"))
  
  



  ####__ 7. Adjusted NumLength  ####
  # scales difference between the sum(numlen) and the reported abundance
  # convers is the difference between abundance and the number measured
  conv_factor <- trawldat %>%
    group_by(id, comname, catchsex, numlen, abundance) %>%
    summarise(
      abund_raw = sum(numlen),               
      convers   =  abund_adj / abund_raw,
      .groups   = "keep")  


  # Merge back and convert the numlen field
  trawldat <- trawldat %>%
    left_join(conv_factor, by = c("id", "comname", "catchsex", "numlen", "abundance")) %>%
    mutate(numlen_adj = numlen * convers, .after = numlen) %>% 
    select(-c(abund_raw, convers))

  # remove conversion factor from environment
  rm(conv_factor, strata_key, strata_select, stratum_area, epu_areas, epu_sf, epu_sp)
  
  
  
  #### Length to Bodymass Conversions  ####
  
  
  
  ####__ 1. Distinct Station & Species Length Info   ####
  
  # For each station we need unique combinations of
  # station_id, species, catchsex, length, adjusted_numlen
  # Record of unique station catches: # rows for each species * sex * length
  trawl_lens <- trawldat %>% 
    filter(is.na(length) == FALSE,
           is.na(numlen) == FALSE,
           numlen_adj > 0) %>% 
    distinct(id, svspp, comname, catchsex, abundance, length,  numlen, numlen_adj, biom_adj) %>% 
    mutate(svspp = as.character(svspp),
           svspp = str_pad(svspp, 3, "left", "0"))
    
  
  
  # Pull distinct records of the stations themselves and metadata that match
  # these are susceptiblt to upstream changes
  trawl_stations <- trawldat %>% 
    select(
      # Tow Identification details
      id, est_year, est_month,  est_day, svvessel, season,
      # physical location details
      decdeg_beglat, decdeg_beglon, avgdepth,
      # NMFS/NEFSC Survey Stratum
      stratum, strat_num, s_area_km2, st_ratio, strat_ntows, tot_s_area,
      # Aggregate Regions/EPU's
      survey_area, epu, epu_area_km2, epu_ratio, epu_ntows, tot_epu_area) %>% 
    distinct() 
  
  
  # recombine with the distinct station info
  trawl_spectra <- trawl_stations %>% 
    left_join(trawl_lens, by = "id")
  
  
  
  ####__ 2. Combine with Growth Coefficients  ####
  
  # This table is a combined table of wigley and fishbase L-W coefficients
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
    mutate(
      b             = as.numeric(b),
      a             = as.numeric(a),
      a             = ifelse(is.na(a) & !is.na(ln_a), exp(ln_a), a),
      ln_a          = ifelse(is.na(ln_a), log(a), ln_a),  # log of a used if ln_a is isn't already there (some fish just had ln_a reported)
      llen          = log(length),
      ind_log_wt    = ln_a + (b * llen),
      ind_weight_kg = exp(ind_log_wt),                    # weight of an individual in size class
      sum_weight_kg = ind_weight_kg * numlen_adj) %>%     # Individual weight * adjusted numlen
    drop_na(ind_weight_kg) %>% 
    select(-ind_log_wt, -llen)
  
  
  # clean up environment
  rm(pass_1, pass_2, fb_trimmed, w_trimmed)
  
  
  
  
  ####__ 3. Use Coefficients to Re-calculate Biomass  ####
  
  # calculate total biomass again using weights from key 
  # make a key for the length weight coefficient sources
  trawl_weights <- trawl_weights %>%  
    arrange(est_year, season, comname, length) %>% 
    mutate(lw_group = str_c(comname, season, catchsex)) 
  
  
  
  
  ####__ 4. Area Weighted Biomass  ####

  # 1. Get the average abundance and biomass by effort in that stratum by year
  # 2. Multiply that mean catch by the area ratio
  
  
  # Good to go up until here:
  # Issue is that we have 2+ scales going on
  # The individual lengths and their length specific bodymass
  # and then the species aggregate biomasses
  # and then the effort by stratum each year to scale with
  
  # The effort (tows) below is from that year and in that stratum
 alb_tow_km2 <- 0.0384 # average area covered by an albatross standard tow in km2
 q <- 1                # catchability coefficient, ideally should shift for different species guilds
  trawl_weights <- trawl_weights  %>% 
    mutate(
      
      # 1. Calculate abundance per tow
      # abundance / ntows for the year in that strata/epu
      abund_tow_s   = numlen_adj / strat_ntows,     
      abund_tow_epu = numlen_adj / epu_ntows,       
      
      # Calculate mean biomass/tow for the "biomass" column
      biom_tow_s    = biom_adj / strat_ntows,
      biom_tow_epu  = biom_adj / epu_ntows,
      
      
      # 2. Calculate stratified mean abundances and fscs biomasses, weighted by the stratum areas
      wt_abund_s = abund_tow_s * st_ratio,                       
      wt_abund_epu   = abund_tow_epu * epu_ratio,
      
      # again with the biomass column from fscs
      wt_biom_s   = biom_tow_s * st_ratio,                       
      wt_biom_epu = biom_tow_epu * epu_ratio,
      
      # Variance - to account for zero catch
      # Not sure how to do this here, number of zeros ddepends on stratification level
      
      # convert from catch rate by area swept to total catch for entire stratum
      # So catch/tow times the total area, divided by how many tows would cover that area
      expanded_abund_s   = round((wt_abund_s * tot_s_area / alb_tow_km2) / q), 
      expanded_abund_epu = round((wt_abund_epu * tot_epu_area/ alb_tow_km2) / q),
      
      # Get fscs biomass from the weighted biomass
      expanded_biom_s   = round((wt_biom_s * tot_s_area / alb_tow_km2) / q), 
      expanded_biom_epu = round((wt_biom_epu * tot_epu_area/ alb_tow_km2) / q),
      
      # Get Biomass from Projected abundances: Biomass = abundance * weight
      expanded_lwbio_s    = sum_weight_kg * expanded_abund_s,
      expanded_lwbio_epu  = sum_weight_kg * expanded_abund_epu
      
    ) 
  
  
  # Remove instances where there were fish that were measured but not weighed
  trawl_weights <- trawl_weights %>% 
    filter(expanded_lwbio_s != -Inf)
  
  return(trawl_weights)
}



####____________________####
####____________________####




####  Aggregate Summary Functions  ####

# For testing
# dat_2020 <- load_ss_data(survdat_source = "2020")




#### 1. Stratified Species Summary  ####
# Stratified Mean = metric / num tows
# Area Stratified Mean = stratified mean * stratum area ratio 
# stratified abundance = strat_mean * total area / area of a tow

#' @title Stratified Means by Species from Size Spectrum Data
#'
#' @description Takes one or more inputs for stratification level, then aggregates abundances and catch rates
#' to the stratum weighted and the epu stratified level. These levels are built into the base area weights
#' and are fixed before this point.
#'
#' @param .data Input data, generated by load_ss_data() function.
#' @param ... Factor columns indicating what you wish to group_by for aggregated metrics ex. season, stratum, epu
#'
#' @return species_data_stratified, a table with aggregate numbers, weights, and cpue by grouping stratum 
#' @export
#'
#' @examples
agg_species_metrics <- function(.data = data, ...){
  
  
  # number of effort (tows) in the groups
  num_tows <- .data %>% 
    group_by(...) %>% 
    summarise(n_stations = n_distinct(id), 
              n_species = n_distinct(comname),
              .groups = "keep")
  
  # Add effort back to the data
  data_stratified <- .data %>% 
    left_join(num_tows)
  
  ## Things that we want to total for the stratum that vary by length increments  ##
  species_data_stratified <- data_stratified  %>% 
    group_by(..., comname) %>% 
    summarise(
      n_tows                   = max(n_stations),
      n_species                = max(n_species),                           # Distinct species
      lw_biomass_kg            = sum(sum_weight_kg),                       # sum weight across species
      fscs_biomass_kg          = sum(biom_adj),                            # total fscs weight
      abundance                = sum(numlen_adj),
      abund_tow                = abundance / n_tows,
      lwbio_tow                = lw_biomass_kg / n_tows ,                  # total weight / tows
      biom_tow                 = fscs_biomass_kg / n_tows,                 # total fscs weight / tows
      sum_survey_abund         = sum(numlen_adj),                             # total number across species
      strat_abundance_s        = sum(expanded_abund_s),                          # total abundance, stratified by nmfs strata
      strat_abundance_epu      = sum(expanded_abund_epu),                        # total abundance, stratified by epu
      mean_ind_length          = weighted.mean(length, numlen_adj),           # average ind length
      mean_ind_bodymass_lw     = weighted.mean(ind_weight_kg, numlen_adj),    # average ind weight
      mean_ind_length_lw       = weighted.mean(length, numlen_adj),           # average ind length
      lw_strat_biomass_s       = sum(expanded_lwbio_s, na.rm = T),              # total biomass across all areas,  weighted by stratum
      lw_strat_biomass_epu     = sum(expanded_lwbio_epu, na.rm = T),            # total biomass across all areas,  weighted by epu
      fscs_strat_biomass_s     = sum(expanded_biom_s),                       # total biomass across all areas,  weighted by stratum
      fscs_strat_biomass_epu   = sum(expanded_biom_epu)                     # total biomass across all areas, weighted by epu
      ,
      
      # #### NOTE:  Would want to do species-level variance here ####
      # # Could do variance at this point, capture the zero catch tows etc.
      # n_zero   = n_tows - n_distinct(id), #works because we are at new group level
      # # biomass
      # zero_var_b = n_zero * (0 - biom_tow)^2,
      # vari_b   = (sum_weight_kg - biom_tow)^2,     # should this be the raw values
      # sh_2_b   = (zero_var_b + sum(vari_b)) / (n_tows - 1),
      # sh_2_b   = ifelse(is.nan(sh_2_b), 0, sh_2_b),
      # #abundance
      # zero_var_a = n_zero * (0 - abund_tow)^2,
      # vari_a   = (numlen_adj - abund_tow)^2,
      # sh_2_a   = (zero_var_a + sum(vari_a)) / (n_tows - 1),
      # sh_2_a   = ifelse(is.nan(sh_2_a), 0, sh_2_a),
      
      
      
      
      .groups = "keep")
  
  
  
  
  
  # Return the data aggreagted to group levels
  return(species_data_stratified)
  
}

# # Test it out
# agg_species_metrics(dat_2020, est_year, epu) %>% glimpse()




####  2. Aggregate Stratified Means and Totals  ####
# This function can be used as a double-check for whether or not the abundances/biomasses
# are equivelant when you do them as individuals vs at an aggregate
# also lets you aggregate the effort and cpue for specific stratifications
# Drops the species specific angle so any cpue or variance is at a group aggregate

#' @title Stratified Mean Biomass and Abundance
#'
#' @description Report aggregate stratified means and totals from survey data.
#' This function follows steps found in the original build, but differs in that 
#' any weighted cpue or biomass adjustment occurs at an aggregate level also 
#' adds in a variance component to capture variability.
#' 
#'
#' @param survdat_clean 
#' @param ... 
#' @param area_stratification 
#'
#' @return
#' @export
#'
#' @examples
agg_strat_metrics <- function(survdat_clean,
                              #aggregation_level = c("year"),
                              ...,
                              area_stratification = "epu"){
  
  # subset columns to sppeed everything up:
  focus_data <- survdat_clean %>% 
    select(id, est_year, est_month, svvessel, season, stratum, s_area_km2, epu, epu_area_km2, 
           comname, catchsex, length, numlen_adj, sum_weight_kg, biom_adj)
  
  
  
  # catch-all for spellings
  if(tolower(area_stratification) %in% c("epu", "stratum", "strata")){
    post_strat <- ifelse(tolower(area_stratification) == "epu", TRUE, FALSE)
    strat_col  <- ifelse(post_strat, sym("epu"), sym("stratum"))
    area_col   <- ifelse(post_strat, sym("epu_area_km2"), sym("s_area_km2"))
  } else {return(message("area stratification not accepted"))}
  
  
  
  #### effort and stratification weights by the chosen aggregation level
  # Get the aggregate total areas
  total_areas <- focus_data %>% 
    group_by( ... ) %>%
    distinct( {{strat_col}}, .keep_all = T) %>%
    summarise(tot_s_area =  sum({{area_col}}, na.rm = T), .groups = "keep")
  # return(total_areas)
  
  # Join those to focus data, get stratum area ratios
  focus_data <- left_join(focus_data, total_areas) %>%
    mutate(stratum_area_ratio = {{area_col}} / tot_s_area)
  # return(focus_data)
  
  
  # get the aggregate effort in each stratification
  q <- 1
  tow_area <- 0.0384
  aggregate_outcomes <- focus_data %>%
    group_by( ..., {{strat_col}}) %>%
    summarise(strat_ntows   = n_distinct(id),
              sum_abund     = sum(numlen_adj),
              sum_lwbio     = sum(sum_weight_kg),
              sum_biom      = sum(biom_adj),
              abund_tow     = sum_abund / strat_ntows,
              lwbio_tow     = sum_lwbio / strat_ntows,
              biom_tow      = sum_biom  / strat_ntows,
              wt_abund_tow  = abund_tow * mean(stratum_area_ratio),
              wt_lwbio_tow  = lwbio_tow * mean(stratum_area_ratio),
              wt_biom_tow   = biom_tow * mean(stratum_area_ratio),
              strat_abund   = wt_abund_tow * (mean(tot_s_area) / tow_area) / q,
              strat_lwbio   = wt_abund_tow * (mean(tot_s_area) / tow_area) / q,
              strat_biom    = wt_abund_tow * (mean(tot_s_area) / tow_area) / q)
  
  
  # # Calculate Variance (depends on stratification structure)
  # if(post_strat == TRUE){
  #   aggregate_outcomes <- aggregate_outcomes %>% 
  #     mutate(
  #       
  #     )
  # }

  return(aggregate_outcomes)
  
  
  
  
  
}


# Testing:
# weights_20 <- load_ss_data(survdat_source = "2020")

# Testing ...
# agg_strat_metrics(survdat_clean = weights_20, est_year, season, area_stratification = "epu")
# agg_strat_metrics(survdat_clean = weights_20, est_year, svvessel, season, area_stratification = "stratum")



#### 3. Annual Summary  ####
# Return an annual summary
ss_annual_summary <- function(survey_data) {

  
  # Summarize on an individual measurement length level
  summ_data <- survey_data %>% 
    group_by(est_year) %>% 
    summarise(
      area                     = "All Areas",
      season                   = "Spring + Fall",
      n_stations               = n_distinct(id),                              # Distinct tows
      n_species                = n_distinct(comname),                         # Distinct species
      lw_biomass_kg            = sum(sum_weight_kg),                          # sum weight across species
      fscs_biomass_kg          = sum(biom_adj),
      lw_biomass_per_station   = lw_biomass_kg / n_stations,                  # total weight / tows
      fscs_biomass_per_station = fscs_biomass_kg /n_stations,
      total_survey_abund       = sum(numlen_adj),                             # total number across species
      strat_abundance_s        = sum(expanded_abund_s),
      strat_abundance_epu      = sum(expanded_abund_epu),
      mean_ind_length          = weighted.mean(length, numlen_adj),           # average ind length
      mean_ind_bodymass_lw     = weighted.mean(ind_weight_kg, numlen_adj),    # average ind weight
      mean_ind_length_lw       = weighted.mean(length, numlen_adj),           # average ind length
      lw_strat_biomass_s       = sum(expanded_lwbio_s, na.rm = T),              # total biomass across all areas,  weighted by stratum
      lw_strat_biomass_epu     = sum(expanded_lwbio_epu, na.rm = T),            # total biomass across all areas,  weighted by epu
      fscs_strat_biomass_s     = sum(expanded_biom_s),                       # total biomass across all areas,  weighted by stratum
      fscs_strat_biomass_epu   = sum(expanded_biom_epu),                     # total biomass across all areas, weighted by epu
      .groups = "keep"
    )  %>% 
    mutate(`Research Vessel` = ifelse(est_year > 2008, "HB", "AL"))
  
  return(summ_data)
  
}

# test: ss_annual_summary(dat_2020)


#### 4. Regional Summary####



# Include region as a group
ss_regional_summary <- function(survey_data, epu = FALSE) {
  
  if(epu == FALSE){
  summ_data <- survey_data %>% 
    group_by(est_year, survey_area) 
  } else if(epu == TRUE){
    summ_data <- survey_data %>% 
      group_by(est_year, epu) 
  }
  
  summ_data <- summ_data %>% 
    summarise(
      season                   = "Spring + Fall",
      n_stations               = n_distinct(id),                              # Distinct tows
      n_species                = n_distinct(comname),                         # Distinct species
      lw_biomass_kg            = sum(sum_weight_kg),                          # sum weight across species
      fscs_biomass_kg          = sum(biom_adj),
      lw_biomass_per_station   = lw_biomass_kg / n_stations,                  # total weight / tows
      fscs_biomass_per_station = fscs_biomass_kg /n_stations,
      total_survey_abund       = sum(numlen_adj),                             # total number across species
      strat_abundance_s        = sum(expanded_abund_s),
      strat_abundance_epu      = sum(expanded_abund_epu),
      mean_ind_length          = weighted.mean(length, numlen_adj),           # average ind length
      mean_ind_bodymass_lw     = weighted.mean(ind_weight_kg, numlen_adj),    # average ind weight
      mean_ind_length_lw       = weighted.mean(length, numlen_adj),           # average ind length
      lw_strat_biomass_s       = sum(expanded_lwbio_s, na.rm = T),              # total biomass across all areas,  weighted by stratum
      lw_strat_biomass_epu     = sum(expanded_lwbio_epu, na.rm = T),            # total biomass across all areas,  weighted by epu
      fscs_strat_biomass_s     = sum(expanded_biom_s),                       # total biomass across all areas,  weighted by stratum
      fscs_strat_biomass_epu   = sum(expanded_biom_epu),                     # total biomass across all areas, weighted by epu
      .groups = "keep"
    )  %>% 
    mutate(`Research Vessel` = ifelse(est_year > 2008, "HB", "AL"))
  
  
return(summ_data)
  
}



#### 5. Seasonal Summary  ####

# Include season as a group
ss_seasonal_summary <- function(survey_data){
  
  summ_data <- survey_data %>%
    group_by(est_year, season) %>%
    summarise(
      area                     = "All Areas",
      n_stations               = n_distinct(id),                              # Distinct tows
      n_species                = n_distinct(comname),                         # Distinct species
      lw_biomass_kg            = sum(sum_weight_kg),                          # sum weight across species
      fscs_biomass_kg          = sum(biom_adj),
      lw_biomass_per_station   = lw_biomass_kg / n_stations,                  # total weight / tows
      fscs_biomass_per_station = fscs_biomass_kg /n_stations,
      total_survey_abund       = sum(numlen_adj),                             # total number across species
      strat_abundance_s        = sum(expanded_abund_s),
      strat_abundance_epu      = sum(expanded_abund_epu),
      mean_ind_length          = weighted.mean(length, numlen_adj),           # average ind length
      mean_ind_bodymass_lw     = weighted.mean(ind_weight_kg, numlen_adj),    # average ind weight
      mean_ind_length_lw       = weighted.mean(length, numlen_adj),           # average ind length
      lw_strat_biomass_s       = sum(expanded_lwbio_s, na.rm = T),              # total biomass across all areas,  weighted by stratum
      lw_strat_biomass_epu     = sum(expanded_lwbio_epu, na.rm = T),            # total biomass across all areas,  weighted by epu
      fscs_strat_biomass_s     = sum(expanded_biom_s),                       # total biomass across all areas,  weighted by stratum
      fscs_strat_biomass_epu   = sum(expanded_biom_epu),                     # total biomass across all areas, weighted by epu
      .groups = "keep"
    )  %>% 
    mutate(`Research Vessel` = ifelse(est_year > 2008, "HB", "AL"))
 
  return(summ_data)
  
}









####____________________####
####____________________####

# ####  2016 Size Spectrum Build  ####
# 
# # needs to be brought up to speed with the 2020 build for stratified abundances and those changes
# #### CHANGE 1: Better Area Measures  ####
# #### CHANGE X : Can this move to data filtering and stratum effort  ####
# ####  CHANGE 2: EPU stratification  ####
# ####  CHANGE 3: Annual Total Stratum Area  ####
# ####  CHANGE 4 : Swept Area Abundances  ####
# 
# load_2016_ss_data <- function(){
#   
#   
#   ####  Paths  ####
#   mills_path <- shared.path(group = "Mills Lab", folder = "")
#   nsf_path <- shared.path(group = "NSF OKN", folder = "")
#   res_path   <- shared.path(os.use = "unix", group = "RES Data", folder = NULL)
#   
#   
#   ####  Load Data  ####
#   load(str_c(mills_path, "Projects/WARMEM/Old survey data/Survdat.RData"))
#   survdat <- clean_names(survdat) %>% 
#     mutate(svspp = as.character(svspp),
#            svspp = str_pad(svspp, 3, "left", "0"))
#   
#   
#   # Extra prep for 2016 data, getting species names via svspp codes
#   spp_classes <- read_csv(here("data/kmills/sppclass.csv"),
#                           col_types = cols()) %>% 
#     clean_names() %>% 
#     mutate(common_name = str_to_lower(common_name),
#            scientific_name = str_to_lower(scientific_name)) %>% 
#     distinct(svspp, comname = common_name, scientific_name)
#   
#   
#   # Add the common names over and format for rest of build
#   survdat <- left_join(survdat, clean_names(spp_classes)) %>% 
#     drop_na(comname) %>% 
#     mutate(cruise6 = str_pad(cruise6, 6, "left", "0"),
#            station = str_pad(station, 3, "left", "0"),
#            stratum = str_pad(stratum, 4, "left", "0"),
#            id = str_c(cruise6, station, stratum)) %>% 
#     select(id, est_year = year, station, stratum, svvessel, season, lat, lon, depth, 
#            surftemp, surfsalin, bottemp, botsalin, svspp, comname, scientific_name, everything()) %>% 
#     mutate(comname = str_to_lower(comname),
#            scientific_name = str_to_lower(scientific_name))
#   
#   
#   
#   
#   ####  Data Prep  ####
#   
#   
#   
#   
#   
#   ####__ 1. Column Changes  ####
#   trawldat <- survdat %>% 
#     mutate(comname = tolower(comname),
#            id = format(id, scientific = FALSE),
#            biomass = ifelse(is.na(biomass) == TRUE & abundance > 0, 0.01, biomass),
#            biomass = ifelse(biomass == 0 & abundance > 0, 0.01, biomass),
#            abundance = ifelse(is.na(abundance) == TRUE & biomass > 0, 1, abundance),
#            abundance = ifelse(abundance == 0 & biomass > 0, 1, abundance),
#            tstrat = str_sub(stratum, 2, 3)) 
#   
#   rm(survdat)
#   
#   ####__ 2. Column Selection  ####
#   trawldat <- trawldat %>%
#     select(id, est_year, svvessel, season, stratum, tstrat, decdeg_beglat = lat, decdeg_beglon = lon,
#            svspp, comname, catchsex, biomass, avgdepth = depth, abundance, length, numlen) 
#   
#   
#   
#   ####__ 3. Data Filtering  ####
#   trawldat <- trawldat %>%
#     filter(
#       # Eliminate Candian Strata and Not in Use Strata
#       stratum >= 01010,
#       stratum <= 01760,
#       stratum != 1310,
#       stratum != 1320,
#       stratum != 1330,
#       stratum != 1350,
#       stratum != 1410,
#       stratum != 1420,
#       stratum != 1490,
#       # Filter to just Spring and Fall
#       season %in% c("SPRING", "FALL"),
#       # Only the Albatross and Henry Bigelow
#       svvessel %in% c("AL", "HB"),
#       est_year >= 1970,
#       est_year < 2020,
#       # Drop na Biomass and Abundance Records
#       !is.na(biomass),
#       !is.na(abundance))
#   
#   
#   
#   ####__ 4. Import Stratum Area Ratios  ####
#   # Survey Stratumm Areas
#   stratum_area <- read_csv(str_c(mills_path, "Projects/NSF_CAccel/Data/strata area.csv"), 
#                            col_types = cols()) %>% 
#     mutate(stratum = as.character(stratum))
#   
#   trawldat <- left_join(trawldat, stratum_area, by = "stratum")
#   trawldat <- arrange(trawldat, id)
#   
#   
#   # Get Total stratum area of all strata (excluding ones we do not care about via left join)
#   total_stratum_area <- trawldat %>% 
#     distinct(stratum, .keep_all = T) %$% 
#     sum(area, na.rm = T)
#   
#   
#   # Calculate strata area relative to total area i.e. stratio
#   trawldat <- trawldat %>% 
#     mutate(stratio = area / total_stratum_area,
#            biom_adj = ifelse(biomass == 0 & abundance > 0, 0.0001, biomass),
#            abund_adj = ifelse(abundance == 0 & biomass > 0, 1, abundance))
#   
#   
#   # Get number of unique tows per stratum
#   yr_strat_effort <- trawldat %>% 
#     group_by(est_year, stratum) %>% 
#     summarise(strat_effort = n_distinct(id))
#   
#   
#   # Add those yearly effort counts back for later
#   trawldat <- trawldat %>% 
#     left_join(yr_strat_effort, by = c("est_year", "stratum"))
#   
#   
#   ###__ 5. Explicit Species Exclusion  ####
#   trawldat <- trawldat %>% 
#     filter(
#       # Exclude the Shrimps
#       svspp %not in% c(285:299, 305, 306, 307, 316, 323, 910:915, 955:961),
#       # Exclude the unidentified fish
#       svspp %not in% c(0, 978, 979, 980, 998)
#     )
#   
#   
#   ####__ 6. Spatial Filtering  ####
#   
#   # If we know which strata represent which area we can use this:
#   # Stratum Key for filtering specific areas
#   strata_key <- list(
#     "Georges Bank"          = as.character(13:23),
#     "Gulf of Maine"         = as.character(24:40),
#     "Southern New England"  = str_pad(as.character(1:12), width = 2, pad = "0", side = "left"),
#     "Mid-Atlantic Bight"    = as.character(61:76))
#   
#   # Add labels to the data
#   trawldat <- trawldat %>% 
#     mutate(survey_area =  case_when(
#       tstrat %in% strata_key$`Georges Bank`         ~ "GB", 
#       tstrat %in% strata_key$`Gulf of Maine`        ~ "GoM",
#       tstrat %in% strata_key$`Southern New England` ~ "SNE",
#       tstrat %in% strata_key$`Mid-Atlantic Bight`   ~ "MAB",
#       TRUE                                          ~ "not found"))
#   
#   # Use strata_select to pull the strata we want
#   strata_select <- c(
#     strata_key$`Georges Bank`, 
#     strata_key$`Gulf of Maine`,
#     strata_key$`Southern New England`,
#     strata_key$`Mid-Atlantic Bight`)
#   
#   
#   # Filtering
#   trawldat <- trawldat %>% filter(tstrat %in% strata_select)
#   
#   
#   
#   
#   ####__ 7. Adjusted NumLength  ####
#   conv_factor <- trawldat %>% 
#     group_by(id, comname, catchsex, numlen, abundance) %>% 
#     summarise(abundance_raw = sum(numlen),         # Total abundance by numlen
#               #abundance = mean(abundance),         # abundance is the same across the groups so mean pulls one value
#               convers =  abundance/abundance_raw)  # convers is the difference between abundance and the number measured
#   
#   
#   # Merge back and convert the numlen field
#   trawldat <- trawldat %>% 
#     left_join(conv_factor, by = c("id", "comname", "catchsex", "numlen", "abundance")) %>% 
#     mutate(numlen_adj = numlen * convers) 
#   
#   # remove conversion factor from environment
#   rm(conv_factor, strata_key, strata_select, stratum_area)
#   
#   
#   
#   
#   #### Length to Bodymass Conversions  ####
#   
#   
#   
#   ####__ 1. Distinct Station & Species Length Info   ####
#   
#   
#   
#   #-#
#   trawl_lens <- trawldat %>% 
#     filter(is.na(numlen) == FALSE,
#            numlen > 0) %>% 
#     distinct(id, svspp, comname, catchsex, length, numlen, numlen_adj, biom_adj) %>% 
#     mutate(svspp = as.character(svspp),
#            svspp = str_pad(svspp, 3, "left", "0"))
#   
#   
#   
#   # Pull distinct records of the stations themselves and metadata that match
#   trawl_stations <- trawldat %>% 
#     select(id, est_year, season, stratum, stratum_num, survey_area, decdeg_beglat, decdeg_beglon, 
#            type, avgdepth, stratio, strat_effort) %>% 
#     distinct() 
#   
#   # recombine with the distinct station info
#   trawl_spectra <- trawl_stations %>% 
#     left_join(trawl_lens, by = "id")
#   
#   
#   
#   ####__ 2. Combine with Growth Coefficients  ####
#   
#   #Now we want to use the lw_combined here instead of just the fishbase lengths
#   lw_combined <- read_csv(here::here("data/biomass_key_combined.csv"),
#                           col_types = cols()) %>% 
#     mutate(svspp = as.character(svspp),
#            svspp = str_pad(svspp, 3, "left", "0"))
#   
#   # Do a priority pass with the filter(lw_combined, source == "wigley)
#   # merge on comname, season, and catchsex
#   w_trimmed <- filter(lw_combined, source == "wigley") %>% 
#     select(source, season, svspp, comname, scientific_name, spec_class, 
#            hare_group, fishery, catchsex, a, b, ln_a)
#   
#   
#   # Do a second pass with the filter(lw_combined, source == "fishbase")
#   # merge on common names only
#   fb_trimmed <- filter(lw_combined, source == "fishbase") %>% 
#     select(source, svspp, comname, scientific_name, spec_class, 
#            hare_group, fishery, a, b, ln_a)
#   #-#
#   
#   
#   
#   # First Pass - Wigley
#   pass_1 <- trawl_spectra %>% 
#     # Join just by svspp to account for name changes
#     select(-comname) %>% 
#     inner_join(w_trimmed)
#   
#   
#   # Second Pass - Fishbase, for the stragglers if any
#   # currently not joining well because lack of svspp and comname agreement 
#   pass_2 <- trawl_spectra %>% 
#     filter(comname %not in% w_trimmed$comname) %>% 
#     inner_join(fb_trimmed)
#   
# 
#   
#   
#   
#   # Join them with bind rows (implicitly drops things that don't have growth coefs)
#   trawl_weights <- bind_rows(pass_1, pass_2) %>% 
#     arrange(est_year, season) %>% 
#     mutate(b = as.numeric(b),
#            a = as.numeric(a),
#            a = ifelse(is.na(a) & !is.na(ln_a), exp(ln_a), a),
#            ln_a = ifelse(is.na(ln_a), log(a), ln_a),
#            ind_log_wt = ln_a + (b * log(length)),
#            ind_weight_kg = exp(ind_log_wt),                    # weight of an individual in size class
#            sum_weight_kg = ind_weight_kg * numlen_adj) %>%      # Individual weight * adjusted numlen
#     drop_na(ind_weight_kg) %>% 
#     select(-ind_log_wt)
#   
#   # clean up environment
#   rm(pass_1, pass_2, fb_trimmed, w_trimmed)
#   
#   
#   
#   ####__ 3. Use Coefficients to Re-calculate Biomass  ####
#   
#   # calculate total biomass again using weights from key 
#   # make a key for the length weight coefficient sources
#   trawl_weights <- trawl_weights %>%  
#     arrange(est_year, season, comname, length) %>% 
#     mutate(lw_group = str_c(comname, season, catchsex)) %>% 
#     as_tibble()
#   
#   
#   rm(trawl_spectra, trawl_stations, trawl_lens)
#   
#   
#   
#   ####__ 4. Area Weighted Biomass  ####
#   # the effort (tows) below is from that year and in that stratum
#   #the  question is whether you can apply stratio here and aggregate up or not
#   trawl_weights <- trawl_weights  %>% 
#     mutate(
#       # 1. Stratum Area Weighted Abundance
#       mean_abund        = numlen_adj / strat_effort,                # how many / number of tows
#       wt_mean_abund     = mean_abund * stratio,                     # adjust by relative stratum size
#       strat_abund       = (wt_mean_abund / .01) * total_stratum_area, # convert to nautical miles, multiply by total area
#       
#       # 2. Stratum Weighted Biomass
#       fscs_strat_bio = biom_adj * strat_abund,     # stratified biomass using fscs biomass
#       lw_strat_bio   = sum_weight_kg * strat_abund # stratified biomass using L-W coefficients
#       
#       # # 1. Stratum Weighted biomass cpue
#       # strat_wt_bio_fscs = stratio * (biom_adj / strat_effort),
#       # strat_wt_bio_lw   = stratio * (sum_weight_kg /strat_effort),
#       # # 2. Stratum Weighted Biomass in kg
#       # fscs_strat_bio = log10((strat_wt_bio_fscs / .01) * total_stratum_area),
#       # lw_strat_bio   = log10((strat_wt_bio_lw / .01) * total_stratum_area)
#       
#       
#     )
#   
#   
#   # Remove instances where there were fish that were measured but not weighed
#   trawl_weights <- trawl_weights %>% 
#     filter(lw_strat_bio != -Inf)
#   
#   # Pass out 2016 weights
#   return(trawl_weights)
#   
#   
#   
# }



####____________________####
####____________________####
