####
#### NEFSC Trawl Data - Size Spectra Build
#### 3/24/2021
####
#### Objective: 
#### Load and clean any survdat pull, without loss of columns. Additional steps to add LW information
#### And to get area expanded (stratified) abundances and biomasses




####  _____Load Packages____  ####
library(janitor)
library(magrittr)
library(here)
library(gmRi)
library(sf)
library(data.table)
library(tidyverse)




####____________________####
#### 1. Survdat Data Prep No Column Drop  ####
#' @title  Load survdat file with standard data filters, keep all columns
#' 
#'
#' @description Processing function to prepare survdat data for size spectra analyses. 
#' Options to select various survdat pulls, or provide your own available.
#' 
#'
#' @param survdat optional starting dataframe in the R environment to run through size spectra build.
#' @param survdat_source String indicating which survdat file to load from box
#'
#' @return Returns a dataframe filtered and tidy-ed for size spectrum analysis.
#' @export
#'
#' @examples
survdat_prep_nodrop <- function(survdat = NULL, survdat_source = "2020"){
  
  ####  Resource Paths  
  box_paths   <- research_access_paths(os.use = "unix")
  mills_path  <- box_paths$mills
  res_path    <- box_paths$res
  
  
  ####  Import supplemental data  ####
  
  # # Load Survey Stratum Area Info - old and in nautical miles (with errors)
  # stratum_area <- read_csv(str_c(caccel_path, "strata area.csv"), col_types = cols())
  
  # Load fresh stratum areas
  stratum_area <- read_csv(str_c(res_path, "NMFS_trawl/Metadata/strata_areas_km2.csv"), col_types = cols())  %>% 
    mutate(stratum = as.character(stratum))
  
  # EPU are info: source slucey_survdat_functions.R and {ecodata}
  epu_areas <- read_csv(str_c(res_path, "NMFS_trawl/Metadata/EPU_areas_km2.csv"), col_types = cols())
  
  
  
  
  ####  Import SURVDAT Data  ####
  
  # Testing:
  #survdat_source <- "bigelow"     ; survdat <- NULL
  #survdat_source <- "most recent" ; survdat <- NULL
  
  # convenience change to make it lowercase
  survdat_source <- tolower(survdat_source)
  
  
  # Build Paths to survdat for standard options
  survdat_path <- switch(EXPR = survdat_source,
                         "2016"        = paste0(mills_path, "Projects/WARMEM/Old survey data/Survdat_Nye2016.RData"),
                         "2019"        = paste0(res_path,   "NMFS_trawl/prior_survdat_data/Survdat_Nye_allseason.RData"),
                         "2020"        = paste0(res_path,   "NMFS_trawl/prior_survdat_data/Survdat_Nye_Aug 2020.RData"),
                         "2021"        = paste0(res_path,   "NMFS_trawl/prior_survdat_data/survdat_slucey_01152021.RData"),
                         "bigelow"     = paste0(res_path,   "NMFS_trawl/2021_survdat/survdat_Bigelow_slucey_01152021.RData"),
                         "most recent" = paste0(res_path,   "NMFS_trawl/2021_survdat/NEFSC_BTS_all_seasons_03032021.RData"), 
                         "bio"         = paste0(res_path,   "NMFS_trawl/2021_survdat/NEFSC_BTS_2021_bio_03192021.RData") )
  
  
  
  # If providing a starting point for survdat pass it in:
  if(is.null(survdat) == FALSE){ 
    trawldat <- survdat %>% clean_names() 
  } else if(is.null(survdat) == TRUE){
    
    # If not then load using the correct path
    load(survdat_path)
    
    
    # Bigelow data doesn't load in as "survdat"
    if(survdat_source == "bigelow"){
      survdat <- survdat.big
      rm(survdat.big)}
    
    # Most recent pulls load a list containing survdat
    if(survdat_source %in% c("bio", "most recent")){
      survdat <- survey$survdat }
    
    # clean names up for convenience
    trawldat <- survdat %>% clean_names() 
  }
  
  # remove survdat once the data is in
  rm(survdat)
  
  
  
  ####_ 1.  Special Steps for Different SURVDAT versions  ####
  
  ####__ a. Missing data flags  ####
  
  # Flags for missing columns that need to be merged in or built
  has_comname  <- "comname" %in% names(trawldat)
  has_id_col   <- "id" %in% names(trawldat)
  has_towdate  <- "est_towdate" %in% names(trawldat)
  has_month    <- "est_month" %in% names(trawldat)
  
  # Flags for renaming or subsetting the data due to presence/absence of columns
  has_year      <- "est_year" %in% names(trawldat)
  has_catchsex  <- "catchsex" %in% names(trawldat)
  has_decdeg    <- "decdeg_beglat" %in% names(trawldat)
  has_avg_depth <- "avgdepth" %in% names(trawldat)
  
  
  
  
  
  
  
  ####__ b. Missing comname  ####
  
  # Use SVSPP to get common names for species
  if(has_comname == FALSE){
    message("no comnames found, merging records in with NMFS_trawl/spp_keys/sppclass.csv")
    # Load sppclass codes and common names
    spp_classes <- read_csv(
      paste0(res_path, "NMFS_trawl/spp_keys/sppclass.csv"), 
      col_types = cols()) %>%
      clean_names() %>%
      mutate(comname         = str_to_lower(common_name),
             scientific_name = str_to_lower(scientific_name)) %>%
      distinct(svspp, comname, scientific_name)
    
    
    # Add the common names over and format for rest of build
    trawldat <- mutate(trawldat, svspp = str_pad(svspp, 3, "left", "0")) %>% 
      left_join(spp_classes, by = "svspp") 
    
  }
  
  
  ####__ c. Missing ID  ####
  if(has_id_col == FALSE) {
    message("creating station id from cruise-station-stratum fields")
    trawldat <- trawldat %>%
      
      # Build ID column
      mutate(cruise6 = str_pad(cruise6, 6, "left", "0"),
             station = str_pad(station, 3, "left", "0"),
             stratum = str_pad(stratum, 4, "left", "0"),
             id = str_c(cruise6, station, stratum))}
  
  
  ####__ d. Field renaming  ####
  
  # Rename select columns for consistency
  if(has_year == FALSE)      {
    message("renaming year column to est_year")
    trawldat <- rename(trawldat, est_year = year) }
  if(has_decdeg == FALSE) {
    message("renaming lat column to decdeg_beglat")
    trawldat <- rename(trawldat, decdeg_beglat = lat) }
  if(has_decdeg == FALSE) {
    message("renaming lon column to decdeg_beglon")
    trawldat <- rename(trawldat, decdeg_beglon = lon) }
  if(has_avg_depth == FALSE)      {
    message("renaming depth column to avgdepth")
    trawldat <- rename(trawldat, avgdepth = depth) }
  
  
  
  ####____ d. build date structure for quick grab of date components
  if(has_towdate == TRUE) {
    message("building month/day columns from est_towdate")
    trawldat <- mutate(trawldat,
                       est_month = str_sub(est_towdate, 6,7),
                       est_month = as.numeric(est_month),
                       est_day   = str_sub(est_towdate, -2, -1),
                       est_day   = as.numeric(est_day))}
  
  
  
  ####_ 2. Column Changes  ####
  trawldat <- trawldat %>% 
    mutate(
      
      # Text Formatting 
      comname = tolower(comname),
      id      = format(id, scientific = FALSE),
      svspp   = as.character(svspp),
      svspp   = str_pad(svspp, 3, "left", "0"),
      season  = tolower(season), 
      
      # Stratum number, 
      # exclude leading and trailing codes for inshore/offshore, 
      # used for matching to stratum areas
      strat_num = str_sub(stratum, 2, 3)) %>%  
    
    # Replace NA's where there is some biomass/abundance
    mutate(biom_adj  = ifelse(biomass == 0 & abundance > 0, 0.0001, biomass), 
           .after = biomass) %>% 
    mutate(abund_adj = ifelse(abundance == 0 & biomass > 0, 1, abundance), 
           .after = abundance)
  
  
  
  

  
  ####_ 3. Row Filtering  ####
  
  # Things filtered:
  # 1. Strata
  # 2. Seasons
  # 3. Year limits
  # 4. Vessels
  # 5. Species Exclusion
  trawldat <- trawldat %>%
    filter(
      
      # Eliminate Canadian Strata and Strata No longer in Use 
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
      season %in% c("spring", "fall"),
      
      # Only the Albatross and Henry Bigelow
      #svvessel %in% c("AL", "HB"),
      
      # Filter years
      est_year >= 1970,
      est_year < 2020,
      
      # Drop NA Biomass and Abundance Records
      !is.na(biomass),
      !is.na(abundance),
      
      # Exclude the Skrimps
      svspp %not in% c(285:299, 305, 306, 307, 316, 323, 910:915, 955:961),
      
      # Exclude the unidentified fish
      svspp %not in% c(0, 978, 979, 980, 998)
    )
  
  
  
  
  
  
  
  
  ####_ 4. Spatial Filtering  ####
  
  # This section assigns EPU's via overlay using Sean Lucey's code,
  # And also merges with stratum area information,
  # these are used to relate catch/effort to physical areas in km squared
  
  #### EPU assignment for survdat stations - function from Sean Luceys "RSurvey" repo
  
  # Area stratification code from Sean Lucey's Repo, renamed to signify what it is
  source(paste0(res_path, "NMFS_trawl/slucey_functions/slucey_survdat_functions.R"))
  
  
  # EPU shapefiles can be loaded from ecodata package
  epu_sf <- ecodata::epu_sf
  epu_sp <- suppressWarnings(as_Spatial(epu_sf))
  
  
  # Rename columns to match expected formats for Sean Lucey's poststrat()
  trawldat <- trawldat %>% 
    rename(CRUISE6 = cruise6,
           STATION = station,
           STRATUM = stratum,
           LAT     = decdeg_beglat,
           LON     = decdeg_beglon)
  
  # Post stratify the station positions using EPU polygons
  trawldat <- trawldat %>% 
    poststrat(survdat = ., stratum = epu_sp, strata.col = "EPU") %>% 
    rename(epu           = newstrata,
           cruise6       = CRUISE6,
           station       = STATION,
           stratum       = STRATUM,
           decdeg_beglat = LAT,
           decdeg_beglon = LON)
  
  
  # Stratum Area Key for which stratum correspond to larger regions we use
  strata_key <- list(
    "Georges Bank"          = as.character(13:23),
    "Gulf of Maine"         = as.character(24:40),
    "Southern New England"  = str_pad(as.character(1:12), width = 2, pad = "0", side = "left"),
    "Mid-Atlantic Bight"    = as.character(61:76))
  
  
  # Add the labels to the data
  trawldat <- trawldat %>%
    mutate(
      survey_area =  case_when(
        strat_num %in% strata_key$`Georges Bank`         ~ "GB",
        strat_num %in% strata_key$`Gulf of Maine`        ~ "GoM",
        strat_num %in% strata_key$`Southern New England` ~ "SNE",
        strat_num %in% strata_key$`Mid-Atlantic Bight`   ~ "MAB",
        TRUE                                             ~ "stratum not in key"))
  
  
  # Optional, Use strata_select to pull the strata we want individually
  strata_select <- c(strata_key$`Georges Bank`, strata_key$`Gulf of Maine`,
                     strata_key$`Southern New England`, strata_key$`Mid-Atlantic Bight`)
  
  
  # Filtering with strata_select
  trawldat <- trawldat %>% filter(strat_num %in% strata_select) %>% 
    mutate(stratum = as.character(stratum))
  
  
  
  ####_ 5. Stratum Area/Effort Ratios  ####
  # Stratum area ratio is the ratio between the area of the select 
  # stratum to the total area of all stratum sampled that year
  
  # Join to the files containing area of each stratum, epu
  trawldat <- trawldat %>% 
    left_join(stratum_area, by = "stratum") %>% 
    left_join(epu_areas, by = "epu") %>% 
    arrange(trawldat, id)
  
  
  
  
  # Get Total stratum area of all strata sampled in a year
  # (excludes ones we do not care about via left join)
  total_stratum_areas <- trawldat %>% 
    group_by(est_year) %>% 
    distinct(stratum, .keep_all = T) %>%  
    summarise(tot_s_area =  sum(s_area_km2, na.rm = T),
              .groups = "keep") %>% ungroup()
  
  # Get total area of EPU's sampled in a year
  total_epu_areas <- trawldat %>% 
    group_by(est_year) %>% 
    distinct(epu, .keep_all = T) %>%  
    summarise(tot_epu_area =  sum(epu_area_km2, na.rm = T), 
              .groups = "keep") %>% ungroup()
  
  
  # Calculate strata area relative to total area i.e. stratio or stratum weights
  trawldat <- trawldat %>% 
    left_join(total_stratum_areas, by = "est_year") %>% 
    left_join(total_epu_areas, by = "est_year") %>% 
    mutate(st_ratio   = s_area_km2 / tot_s_area,
           epu_ratio  = epu_area_km2 / tot_epu_area) 
  
  # We have to total areas, now we want effort in each
  # Number of unique tows per stratum
  yr_strat_effort <- trawldat %>% 
    group_by(est_year, stratum) %>% 
    summarise(strat_ntows = n_distinct(id), 
              .groups = "keep") %>% ungroup()
  
  # Number of unique tows per EPU
  yr_epu_effort <-  trawldat %>% 
    group_by(est_year, epu) %>% 
    summarise(epu_ntows = n_distinct(id), 
              .groups = "keep") %>% ungroup()
  
  
  # Add those yearly effort counts back for later 
  #(area stratified abundance)
  trawldat <- trawldat %>% 
    left_join(yr_strat_effort, by = c("est_year", "stratum")) %>% 
    left_join(yr_epu_effort, by = c("est_year", "epu"))
  
  
  
  
  
  ####_ 6. Adjusting NumLength to Match Abundance  ####
  # Sometimes there are more/less measured than initially tallied* in abundance
  
  
  # If catchsex is not a column then total abundance is assumed pooled
  if(has_catchsex == TRUE){
    abundance_groups <- c("id", "comname", "catchsex")
  } else {
    message("catchsex column not found, ignoring sex for numlen adjustments")
    abundance_groups <- c("id", "comname")}
  
  
  
  # Get the abundance value for each sex 
  # arrived at by summing across each length
  abundance_check <- trawldat %>%
    group_by(!!!syms(abundance_groups), abundance) %>%
    summarise(
      abund_actual = sum(numlen),               
      n_len_class  = n_distinct(length),
      .groups      = "keep") %>% ungroup()
  
  # Get the ratio between the original abundance column 
  # and the sum of numlen we just grabbed
  conv_factor <- trawldat %>% 
    #distinct(id, comname, catchsex, length, abund_adj) %>% 
    distinct(!!!syms(abundance_groups), length, abund_adj) %>% 
    #inner_join(abundance_check, by = c("id", "comname", "catchsex")) %>% 
    inner_join(abundance_check, by = abundance_groups) %>% 
    mutate(convers = abund_adj / abund_actual)
  
  
  
  # Merge back and convert the numlen field
  # original numlen * conversion factor = numlength adjusted
  survdat_processed <- trawldat %>%
    #left_join(conv_factor, by = c("id", "comname", "catchsex", "length", "abundance", "abund_adj")) %>%
    left_join(conv_factor, by = c(abundance_groups, "length", "abundance", "abund_adj")) %>%
    mutate(numlen_adj = numlen * convers, .after = numlen) %>% 
    select(-c(abund_actual, convers))
  
  # remove conversion factor from environment
  rm(abundance_check, conv_factor, strata_key, strata_select, stratum_area, epu_areas, epu_sf, epu_sp)
  
  
  
  
  
  ####_ 7. Distinct Station & Species Length Info   ####
  
  # For each station we need unique combinations of
  # station_id, species, catchsex, length, adjusted_numlen
  # to capture what and how many of each length fish is caught
  
  # Record of unique station catches: 
  # One row for every species * sex * length, combination in the data
  trawl_lens <- survdat_processed %>% 
    filter(is.na(length) == FALSE,
           is.na(numlen) == FALSE,
           numlen_adj > 0) 
  
  
  # Do we want to just keep all the station info here as well?
  # question to answer is whether any other columns repeat,
  # or if these are the only ones
  trawl_spectra <- trawl_lens %>%
    distinct(id, svspp, comname, catchsex, abundance, n_len_class,
             length, numlen, numlen_adj, biom_adj, .keep_all = TRUE)
  
  
  
  
  
  # Reorder the main columns
  
  
  
  
  
  # Return the dataframe
  # Row for each length class of every species caught
  return(trawl_spectra)
  
}




####____________________####
#### 2. Length to Bodymass Conversions  ####

# testing data
# survdat_clean <- survdat_prep(survdat_source = "2020")

#' @title Add species length weight information, calculate expected biomass
#'
#' @param survdat_clean Survdat data, after usual preparations are completed.
#' These include removal of old strata, labeling of areas of interest, and inclusion
#' of the annual effort in each.
#'
#' @return
#' @export
#'
#' @examples
add_lw_info <- function(survdat_clean, cutoff = FALSE){
  
  
  ####__ 1. Match Species with Growth Coefficients  ####
  
  # This table is a combined table of wigley and fishbase L-W coefficients
  lw_combined <- read_csv(paste0(res_path, "NMFS_trawl/length_weight_keys/fishbase_wigley_combined_key.csv"), col_types = cols()) %>% 
    mutate(svspp = str_pad(svspp, 3, "left", "0"),
           season = tolower(season))
  
  
  # Do a priority pass with the filter(lw_combined, source == "wigley)
  # merge on comname, season, and catchsex
  wigley_coefficients <- filter(lw_combined, source == "wigley") %>% 
    select(source, season, svspp, comname, scientific_name, spec_class, 
           hare_group, fishery, catchsex, a, b, ln_a)
  
  
  # Do a second pass with the filter(lw_combined, source == "fishbase")
  # merge on common names only
  fishbase_coefficients <- filter(lw_combined, source == "fishbase") %>% 
    select(source, -svspp, comname, scientific_name, spec_class, 
           hare_group, fishery, a, b, ln_a)  
  
  
  # First Pass - Wigley
  # Join just by svspp to account for name changes
  pass_1 <- survdat_clean %>% 
    select(-comname) %>% 
    inner_join(wigley_coefficients)
  
  
  # Second Pass - Fishbase, for the stragglers if any
  pass_2 <- survdat_clean %>% 
    filter(comname %not in% wigley_coefficients$comname,
           svspp %not in% pass_1$svspp) %>% 
    inner_join(fishbase_coefficients) 
  
  
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
  rm(pass_1, pass_2, fishbase_coefficients, wigley_coefficients)
  
  
  
  
  ####__ 2. Use Coefficients to Re-calculate Biomass  ####
  
  # calculate total biomass again using weights from key 
  # make a key for the length weight coefficient sources
  survdat_weights <- trawl_weights %>%  
    arrange(est_year, season, comname, length) %>% 
    mutate(lw_group = str_c(comname, season, catchsex)) 
  
  
  
  
  
  
  ####__ 3. Drop comnames that don't align well with BIOMASS
  
  # these species were dropped at 50% mismatch threshold
  # code: 02_survdat_stratification_validation
  # list updated : 3/24/2021
  cutoff_25 <- c(
    "acadian redfish"          ,      "alewife"                 ,
    "american plaice"          ,      "american shad"           ,
    "atlantic angel shark"     ,      "atlantic cod"            ,
    "atlantic croaker"         ,      "atlantic halibut"        ,
    "atlantic herring"         ,      "atlantic mackerel"       ,
    "atlantic sharpnose shark" ,      "atlantic sturgeon"       ,
    "atlantic wolffish"        ,      "barndoor skate"          ,
    "black sea bass"           ,      "bluefish"                ,
    "bluntnose stingray"       ,      "buckler dory"            ,
    "bullnose ray"             ,      "butterfish"              ,
    "chain dogfish"            ,      "clearnose skate"         ,
    "cownose ray"              ,      "cunner"                  ,
    "cusk"                     ,      "fawn cusk-eel"           ,
    "fourspot flounder"        ,      "haddock"                 ,
    "little skate"             ,      "longhorn sculpin"        ,
    "northern kingfish"        ,      "northern searobin"       ,
    "ocean pout"               ,      "offshore hake"           ,
    "pollock"                  ,      "red hake"                ,
    "rosette skate"            ,      "roughtail stingray"      ,
    "round herring"            ,      "sandbar shark"           ,
    "sea raven"                ,      "silver hake"             ,
    "smooth butterfly ray"     ,      "smooth dogfish"          ,
    "smooth skate"             ,      "southern stingray"       ,
    "spanish mackerel"         ,      "spiny butterfly ray"     ,
    "spiny dogfish"            ,      "spot"                    ,
    "spotted hake"             ,      "striped bass"            ,
    "summer flounder"          ,      "tautog"                  ,
    "thorny skate"             ,      "weakfish"                ,
    "white hake"               ,      "windowpane flounder"     ,
    "winter flounder"          ,      "winter skate"            ,
    "witch flounder"           ,      "yellowtail flounder" 
  )
  
  # Filter to use species that meet cutoff criteria
  # source: 02_survdat_stratification_validation
  if(cutoff == TRUE){
    survdat_weights <- survdat_weights %>% filter(comname %in% cutoff_25)
  }
  
  
  
  
  return(survdat_weights)
}





####____________________####
####  3. Area Stratified Abundance / Biomass  ####


# Add area stratified biomass function
#' @title Add Survey Area Stratified Abundances and Biomasses
#' 
#' @description Take the survdat data paired with length weight relationships and
#' return estimates of area stratified catch rates and their expected abundances and
#' biomasses when applied to the total areas of stratum.
#'
#' @param survdat_weights Input dataframe, produced by add_lw_info
#' @param include_epu Flag for calculating the EPU rates in addition to the stratum regions we use.
#'
#' @return
#' @export
#'
#' @examples
add_area_stratification <- function(survdat_weights, include_epu = F){ 
  
  # Constants:
  # average area covered by an albatross standard tow in km2
  alb_tow_km2 <- 0.0384 
  
  # catchability coefficient, ideally should shift for different species guilds
  q <- 1                
  
  
  # Derived Estimates:
  survdat_weights <- survdat_weights  %>% 
    mutate(
      # Abundance per tow
      # abundance / ntows for the year within that strata/epu
      abund_tow_s   = numlen_adj / strat_ntows,    
      
      # Biomass is repeated across length classes at each station by species
      # the number of length classes is tallied where the conversion factor is done
      biom_per_lclass = (biom_adj / n_len_class),
      
      # Mean biomass/tow for the "biomass" column
      biom_tow_s      = biom_per_lclass / strat_ntows,
      
      # Stratified mean abundance, weighted by the stratum areas
      wt_abund_s     = abund_tow_s * st_ratio, 
      
      # Stratified mean BIOMASS
      wt_biom_s   = biom_tow_s * st_ratio,
      
      # convert from catch rate by area swept to total catch for entire stratum
      # So catch/tow times the total area, divided by how many tows would cover that area
      expanded_abund_s   = round((wt_abund_s * tot_s_area / alb_tow_km2) / q),
      
      # Total BIOMASS from the weighted biomass
      expanded_biom_s   = round((wt_biom_s * tot_s_area / alb_tow_km2) / q), 
      
      # Total lw-biomass from Projected abundances: Biomass = abundance * lw_weight
      expanded_lwbio_s    = sum_weight_kg * expanded_abund_s
      
    ) 
  
  ####  Optional - weighted by EPU areas
  if(include_epu == TRUE){
    survdat_weights <- survdat_weights  %>% 
      mutate(
        # Abundance per tow
        abund_tow_epu = numlen_adj / epu_ntows,       
        # Mean biomass/tow for the BIOMASS column
        biom_tow_epu    = biom_per_lclass / epu_ntows,
        # Stratified mean abundances, weighted by the stratum areas
        wt_abund_epu   = abund_tow_epu * epu_ratio,
        # Stratified mean BIOMASS
        wt_biom_epu = biom_tow_epu * epu_ratio,
        # Total catch for entire stratum
        expanded_abund_epu = round((wt_abund_epu * tot_epu_area/ alb_tow_km2) / q),
        # Total Biomass from the weighted biomass
        expanded_biom_epu = round((wt_biom_epu * tot_epu_area/ alb_tow_km2) / q),
        # LW Biomass from Expanded abundances: Biomass = abundance * lw_weight
        expanded_lwbio_epu  = sum_weight_kg * expanded_abund_epu)} 
  
  
  # Remove instances where there were fish that were measured but not weighed
  survdat_weights <- survdat_weights %>% 
    filter(expanded_lwbio_s != -Inf)
  
  return(survdat_weights)
}

# survdat_clean <- survdat_prep(survdat_source = "2020")
# weights_20 <- add_lw_info(survdat_clean = survdat_clean)
# strat_biomass <- add_area_stratification(survdat_weights = weights_20, include_epu = F)




####____________________####

####  Taking Out Individual Components  ####


# Grab all station location and other abiotic details
# with and without filters
pull_station_details <- function(survdat_raw, filter = TRUE){
  
  # filter if you want
  if(filter == TRUE){
    survdat_raw <- survdat_raw %>% 
      filter()
  }
  
  # If not, just pull distinct values
  
  
  
}

# Pull station totals for things that do not change across a
# single station for each species i.e. station totals for a species for either sex
pull_station_totals <- function(){
  
  survdat_weights %>% 
    distinct(id, station, svvessel, season, 
             svspp, comname, catchsex, 
             abundance, biom_adj, sum_weight_kg,
             expanded_abund_s, expanded_biom_s, expanded_lwbio_s)
  
}




# Pull the length specific totals
# i.e. how many of each length are caught, and what that biomass should be
pull_numlen_details <- function(){
  survdat_weights %>% 
    distinct(id, station, svvessel, season, 
             svspp, comname, catchsex, numlen_adj, lw_biom)
}




####____________________####

####  Aggregation Summary Functions  ####

# For testing
# dat_2020 <- strat_biomass


####____________________####


#### 1. Group Stratified Species Summaries  ####

#' @title Stratified Means by Species from Size Spectrum Data
#'
#' @description Takes one or more inputs for stratification level to use in addition to comname.
#' Aggregates abundances and catch rates will be totaled for desired stratification
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
agg_species_metrics <- function(.data = data, ..., include_epu = FALSE){
  
  
  # Identify number of effort (tows) in the groups
  num_tows <- .data %>% 
    group_by(...) %>% 
    summarise(n_stations = n_distinct(id), 
              n_species = n_distinct(comname),
              .groups = "keep")
  
  # Add effort back to the data
  data_stratified <- .data %>% 
    left_join(num_tows)
  
  
  ## Group Aggregate Metrics ##
  species_data_stratified <- data_stratified  %>% 
    group_by(..., comname) %>% 
    summarise(
      n_tows                   = max(n_stations),                           # Distinct stations
      n_species                = max(n_species),                            # Distinct species
      lw_biomass_kg            = sum(sum_weight_kg),                        # Sum lw weight across species
      fscs_biomass_kg          = sum(biom_per_lclass),                      # Total BIOMASS, split evenly across length classes
      abundance                = sum(numlen_adj),                           # Abundance
      abund_tow                = abundance / n_tows,                        # Abundance / tow
      lwbio_tow                = lw_biomass_kg / n_tows ,                   # total weight / tows
      biom_tow                 = fscs_biomass_kg / n_tows,                  # total fscs weight / tows
      mean_ind_length          = weighted.mean(length, numlen_adj),         # average ind length
      mean_ind_bodymass_lw     = weighted.mean(ind_weight_kg, numlen_adj),  # average ind weight
      mean_ind_length_lw       = weighted.mean(length, numlen_adj),         # average ind length
      strat_abundance_s        = sum(expanded_abund_s),                     # total abundance, stratified by nmfs strata
      lw_strat_biomass_s       = sum(expanded_lwbio_s, na.rm = T),          # total lw biomass across all areas,  weighted by stratum
      fscs_strat_biomass_s     = sum(expanded_biom_s),                      # total biomass across all areas, weighted by epu
      .groups = "keep") %>% 
    ungroup()
  
  # Optional step to include EPU area aggregates
  if(include_epu == T){
    epu_strat <- data_stratified  %>% 
      group_by(..., comname) %>% 
      summarise(
        strat_abundance_epu      = sum(expanded_abund_epu),               # total abundance, stratified by epu
        lw_strat_biomass_epu     = sum(expanded_lwbio_epu, na.rm = T),    # total lw biomass across all areas, weighted by epu
        fscs_strat_biomass_epu   = sum(expanded_biom_epu),                # total biomass across all areas, weighted by epu
        .groups = "keep") %>% 
      ungroup()
    
    # Add to stratum area aggregates
    species_data_stratified <- bind_cols(species_data_stratified, epu_strat)
    
  }
  
  
  # Return the data aggreagted to group levels
  return(species_data_stratified)
  
}

# # Test it out
# agg_species_metrics(dat_2020, est_year, epu) %>% glimpse()




####  2. Aggregate Stratified Means and Totals  ####


# This function can be used as a double-check for whether or not the abundances/biomasses
# are equivalent when you do them as individuals vs at an aggregate
# also lets you aggregate the effort and cpue for specific stratification
# Drops the species specific angle so any cpue or variance is at a group aggregate

#' @title Stratified Mean Biomass and Abundance by Stratum or EPU by year
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
agg_strat_metrics <- function(.data,
                              ...,
                              area_stratification = "stratum"){
  
  # subset columns to speed everything up:
  focus_data <- .data %>% 
    select(id, est_year, est_month, svvessel, season, stratum, s_area_km2, epu, epu_area_km2, 
           comname, catchsex, length, numlen_adj, sum_weight_kg, biom_adj, biom_per_lclass)
  
  
  
  # catch-all for spellings, determine whether to use EPU area or stratum area
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
              sum_biom      = sum(biom_per_lclass),
              abund_tow     = sum_abund / strat_ntows,
              lwbio_tow     = sum_lwbio / strat_ntows,
              biom_tow      = sum_biom  / strat_ntows,
              wt_abund_tow  = abund_tow * mean(stratum_area_ratio),
              wt_lwbio_tow  = lwbio_tow * mean(stratum_area_ratio),
              wt_biom_tow   = biom_tow * mean(stratum_area_ratio),
              strat_abund   = wt_abund_tow * (mean(tot_s_area) / tow_area) / q,
              strat_lwbio   = wt_abund_tow * (mean(tot_s_area) / tow_area) / q,
              strat_biom    = wt_abund_tow * (mean(tot_s_area) / tow_area) / q) %>% 
    ungroup()
  
  
  
  return(aggregate_outcomes)
  
}


# Testing:
# weights_20 <- load_ss_data(survdat_source = "2020")

# Testing ...
# agg_strat_metrics(survdat_clean = weights_20, est_year, season, area_stratification = "epu")
# agg_strat_metrics(survdat_clean = weights_20, est_year, svvessel, season, area_stratification = "stratum")



#### 3. Annual Summary  ####
# Return an annual summary
ss_annual_summary <- function(survey_data, include_epu = F) {
  
  
  # Summarize on an individual measurement length level
  summ_data <- survey_data %>% 
    group_by(est_year) %>% 
    summarise(
      area                     = "All Areas",
      season                   = "Spring + Fall",
      n_stations               = n_distinct(id),                              # Distinct tows
      n_species                = n_distinct(comname),                         # Distinct species
      lw_biomass_kg            = sum(sum_weight_kg),                          # sum weight across species
      fscs_biomass_kg          = sum(biom_per_lclass),
      lw_biomass_per_station   = lw_biomass_kg / n_stations,                  # total weight / tows
      fscs_biomass_per_station = fscs_biomass_kg /n_stations,
      total_survey_abund       = sum(numlen_adj),                             # total number across species
      mean_ind_length          = weighted.mean(length, numlen_adj),           # average ind length
      mean_ind_bodymass_lw     = weighted.mean(ind_weight_kg, numlen_adj),    # average ind weight
      mean_ind_length_lw       = weighted.mean(length, numlen_adj),           # average ind length
      strat_abundance_s        = sum(expanded_abund_s),
      lw_strat_biomass_s       = sum(expanded_lwbio_s, na.rm = T),            # total biomass across all areas,  weighted by stratum
      fscs_strat_biomass_s     = sum(expanded_biom_s),                        # total biomass across all areas,  weighted by stratum
      .groups = "keep"
    )  %>% 
    mutate(`Research Vessel` = ifelse(est_year > 2008, "HB", "AL"))
  
  
  
  # Optional step to include EPU area aggregates
  if(include_epu == T){
    epu_strat <- survey_data  %>% 
      group_by(est_year) %>% 
      summarise(
        strat_abundance_epu      = sum(expanded_abund_epu),               # total abundance, stratified by epu
        lw_strat_biomass_epu     = sum(expanded_lwbio_epu, na.rm = T),    # total lw biomass across all areas, weighted by epu
        fscs_strat_biomass_epu   = sum(expanded_biom_epu),                # total biomass across all areas, weighted by epu
        .groups = "keep") %>% 
      ungroup()
    
    # Add to stratum area aggregates
    summ_data <- bind_cols(summ_data, epu_strat)
    
  }
  
  return(summ_data)
  
}

# test: ss_annual_summary(dat_2020)


#### 4. Regional Summary####



# Include region as a group
ss_regional_summary <- function(survey_data, epu_regions = F, include_epu = FALSE) {
  
  # toggle whether to use regions or epu's
  if(epu_regions == FALSE){
    summ_data <- survey_data %>% 
      group_by(est_year, survey_area) 
  } else if(epu_regions == TRUE){
    summ_data <- survey_data %>% 
      group_by(est_year, epu) 
  }
  
  # Get group metrics
  summ_data <- summ_data %>% 
    summarise(
      season                   = "Spring + Fall",
      n_stations               = n_distinct(id),                              # Distinct tows
      n_species                = n_distinct(comname),                         # Distinct species
      lw_biomass_kg            = sum(sum_weight_kg),                          # sum weight across species
      fscs_biomass_kg          = sum(biom_per_lclass),
      lw_biomass_per_station   = lw_biomass_kg / n_stations,                  # total weight / tows
      fscs_biomass_per_station = fscs_biomass_kg /n_stations,
      total_survey_abund       = sum(numlen_adj),                             # total number across species
      mean_ind_length          = weighted.mean(length, numlen_adj),           # average ind length
      mean_ind_bodymass_lw     = weighted.mean(ind_weight_kg, numlen_adj),    # average ind weight
      mean_ind_length_lw       = weighted.mean(length, numlen_adj),           # average ind length
      strat_abundance_s        = sum(expanded_abund_s),
      lw_strat_biomass_s       = sum(expanded_lwbio_s, na.rm = T),            # total biomass across all areas,  weighted by stratum
      fscs_strat_biomass_s     = sum(expanded_biom_s),                        # total biomass across all areas,  weighted by stratum
      .groups = "keep")  %>% 
    ungroup()
  
  
  # Optional step to include EPU area aggregates
  if(include_epu == T){
    epu_strat <- summ_data  %>% 
      group_by(est_year) %>% 
      summarise(
        strat_abundance_epu      = sum(expanded_abund_epu),               # total abundance, stratified by epu
        lw_strat_biomass_epu     = sum(expanded_lwbio_epu, na.rm = T),    # total lw biomass across all areas, weighted by epu
        fscs_strat_biomass_epu   = sum(expanded_biom_epu),                # total biomass across all areas, weighted by epu
        .groups = "keep") %>% 
      ungroup()
    
    # Add to stratum area aggregates
    summ_data <- bind_cols(summ_data, epu_strat)
    
  }
  
  return(summ_data)
  
  
}



#### 5. Seasonal Summary  ####

# Include season as a group
ss_seasonal_summary <- function(survey_data, include_epu = F){
  
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
      mean_ind_length          = weighted.mean(length, numlen_adj),           # average ind length
      mean_ind_bodymass_lw     = weighted.mean(ind_weight_kg, numlen_adj),    # average ind weight
      mean_ind_length_lw       = weighted.mean(length, numlen_adj),           # average ind length
      lw_strat_biomass_s       = sum(expanded_lwbio_s, na.rm = T),              # total biomass across all areas,  weighted by stratum
      fscs_strat_biomass_s     = sum(expanded_biom_s),                       # total biomass across all areas,  weighted by stratum
      .groups = "keep"
    )  %>% 
    mutate(`Research Vessel` = ifelse(est_year > 2008, "HB", "AL"))
  
  # Optional step to include EPU area aggregates
  if(include_epu == T){
    epu_strat <- survey_data  %>% 
      group_by(est_year, season) %>% 
      summarise(
        strat_abundance_epu      = sum(expanded_abund_epu),               # total abundance, stratified by epu
        lw_strat_biomass_epu     = sum(expanded_lwbio_epu, na.rm = T),    # total lw biomass across all areas, weighted by epu
        fscs_strat_biomass_epu   = sum(expanded_biom_epu),                # total biomass across all areas, weighted by epu
        .groups = "keep") %>% 
      ungroup()
    
    # Add to stratum area aggregates
    summ_data <- bind_cols(summ_data, epu_strat)
    
  }
  
  return(summ_data)
  
}