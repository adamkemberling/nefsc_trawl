####
#### NEFSC Trawl Data - Size Spectra Build
#### 9/28/2020
####
#### Objective: 
#### Load 2020 "survdat" data from Janet Nye, perform all filtering and adjustment for Size Spectra Analyses
####

load_ss_data <- function(){
  #### Packages  ####
  library(janitor)
  library(tidyverse)
  library(magrittr)
  library(here)
  library(gmRi)
  
  
  ####  Paths  ####
  mills_path <- shared.path(group = "Mills Lab", folder = "")
  nsf_path <- shared.path(group = "NSF OKN", folder = "")
  res_path   <- shared.path(os.use = "unix", group = "RES Data", folder = NULL)
  
  
  
  ####  Support Functions  ####
  source(here("R/support/sizeSpectra_support.R"))
  `%notin%` <-  purrr::negate(`%in%`) # to filter out ones that have matches
  
  
  
  ####  Data  ####
  
  # 2020 groundfish survey data
  load(here("data/NEFSC/Survdat_Nye_Aug 2020.RData"))
  survdat <- clean_names(survdat)
  
  
  # Survey Stratumm Areas
  stratum_area <- read_csv(str_c(mills_path, "Projects/NSF_CAccel/Data/strata area.csv"), 
                           col_types = cols())
  
  
  
  # Fishbase Growth Coefficients
  nefsc_lw <- read_csv(here("data/NEFSC/nefsc_lw_key_filled.csv"),
                       guess_max = 1e3,
                       col_types = cols())
  
  # Wrigley Paper, Load length weight coefficients
  load(paste0(res_path, "/NMFS_trawl/lwreg.Rdata"))
  
  
  # Species class/groups based on life history
  spp_classes <- read_csv(here("data/kmills/sppclass.csv"),
                          col_types = cols())
  
  ####____________________####
  ####  Survey Data Prep  ####
  
  
  ####__ 1. Column Changes  ####
  trawldat <- survdat %>% 
    mutate(comname = tolower(comname),
           id = format(id, scientific = FALSE),
           biomass = ifelse(is.na(biomass) == TRUE & abundance > 0, 0.01, biomass),
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
           biom.adj = ifelse(biomass == 0 & abundance > 0, 0.0001, biomass),
           abund.adj = ifelse(abundance == 0 & biomass > 0, 1, abundance),
           indwt = biom.adj / abund.adj)
  
  
  
  ###__ 5. Explicit Species Exclusion  ####
  trawldat <- trawldat %>% 
    filter(
      # Exclude the Shrimps
      svspp %notin% c(285:299, 305, 306, 307, 316, 323, 910:915, 955:961),
      # Exclude the unidentified fish
      svspp %notin% c(0, 978, 979, 980, 998)
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
              abundance = mean(abundance),         # abundance is the same across the groups so mean pulls one value
              convers =  abundance/abundance_raw)  # convers is the difference between abundance and the number measured
  
  
  # Merge back and convert the numlen field
  trawldat <- trawldat %>% 
    left_join(conv_factor, by = c("id", "comname", "catchsex", "numlen", "abundance")) %>% 
    mutate(numlen_adj = numlen * convers,
           numlen_adj = round(numlen_adj, 2)) 
  
  # remove conversion factor from environment
  rm(conv_factor, strata_key, strata_select, stratum_area)
  
  
  ####____________________####
  ####  Growth Coefficients Prep  ####
  
  ####__ 1.  Fishbase Coefficients  ####
  # Fill in NA growth coefficients using related species (typically same genera)
  nefsc_lw <- nefsc_lw %>% 
    mutate(a = as.numeric(a),
           b = as.numeric(b),
           ln_a = log(a),
           related_ln_a = log(related_a),
           ln_a = ifelse(is.na(ln_a), related_ln_a, ln_a),
           a = ifelse(is.na(a), exp(ln_a), a))
  
  
  
  ####__ 2. Wrigley Paper Coefficients  ####
  lwreg <- lwreg %>% clean_names()  
  lwreg <- lwreg %>%  
    mutate(comname = str_to_lower(common_name), .before = scientific_name,
           common_name = NULL,
           scientific_name = str_to_lower(scientific_name),
           svspp = str_pad(svspp, width = 3, side = "left", pad = "0")) %>% 
    mutate(ln_a = lna, 
           a = exp(ln_a),
           lna = NULL) %>% 
    select(svspp, comname, scientific_name, season, catchsex, a, ln_a, everything())
  
  
  
  ####__ 3. Classes & Economic Groups  ####
  # Species Life History & Fisheries Relevance
  spp_classes <- spp_classes %>% 
    clean_names() %>% 
    mutate(common_name = str_to_lower(common_name),
           scientific_name = str_to_lower(scientific_name)) %>% 
    pivot_longer(gf:inv, names_to = "spec_class", values_to = "flag") %>% 
    pivot_longer(com:nc, names_to = "fishery", values_to = "flag2") %>% 
    filter(is.na(flag) == FALSE,
           is.na(flag2) == FALSE) %>% 
    select(svspp, comname = common_name, everything(), -flag, -flag2)
  
  
  # Merge the growth coefficients with the species classes
  wigley_lw <- left_join(lwreg, spp_classes, by = c("svspp", "comname", "scientific_name"))
  
  
  ####__ 4.  Combining Growth Coefficient Sources  ####
  
  
  # 1. use these growth coefficients to get biomasses
  # 2. Potentially merge in my fishbase ones so that there is just one master
  
  wigley_lw <- wigley_lw %>% mutate(source = "wigley") %>% select(source, everything())
  nefsc_lw  <- nefsc_lw %>% mutate(source = "fishbase") %>% select(source, everything())
  lw_combined <- full_join(nefsc_lw, wigley_lw, by = c("source", "comname", "b", "a", "ln_a")) 
  
  # Fill in gaps for things that should be consistent
  # Lot of matching by common name and pulling the first value where it matches to fill the NA values
  lw_combined <- lw_combined %>% 
    arrange(comname) %>% 
    mutate(
      hare_group = ifelse(is.na(hare_group),
                          nefsc_lw[ match(comname,  nefsc_lw$comname)[1], "hare_group"][[1]],
                          hare_group),
      svspp = ifelse(is.na(svspp),
                     wigley_lw[ match(comname, wigley_lw$comname)[1], "svspp"][[1]],
                     svspp),
      scientific_name = ifelse(is.na(scientific_name),
                               wigley_lw[ match(comname, wigley_lw$comname)[1], "scientific_name"][[1]],
                               scientific_name),
      fishery = ifelse(is.na(fishery),
                       wigley_lw[ match(comname, wigley_lw$comname)[1], "fishery"][[1]],
                       fishery),
      spec_class = ifelse(is.na(spec_class),
                          wigley_lw[ match(comname, wigley_lw$comname)[1], "spec_class"][[1]],
                          spec_class))
  
  # Grab and order columns we need
  lw_combined <- lw_combined %>% 
    select(source, svspp, comname, scientific_name, spec_class, hare_group, fishery, season, catchsex,
           b, a, ln_a, units, related_species, related_a, related_b, related_ln_a,
           wigleymin, wigleymax, wigleyvalue)
  
  
  # Remove the lw building components
  rm(wigley_lw, nefsc_lw, spp_classes, lwreg)
  
  # Proceed to Size Spectra Steps:
  
  ####____________________####
  #### Length to Bodymass Conversions  ####
  
  
  
  ####__ 1. Distinct Station & Species Length Info   ####
  
  # For each station we need unique combinations of
  # station_id, species, catchsex, length, adjusted_numlen
  # Record of unique station catches: # rows for each species * sex * length
  trawl_lens <- trawldat %>% 
    filter(is.na(numlen) == FALSE,
           numlen > 0) %>% 
    distinct(id, comname, catchsex, length, numlen, numlen_adj) 
  
  
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
  
  # Do a priority pass with the filter(lw_combined, source == "wigley)
  # merge on comname, season, and catchsex
  w_trimmed <- filter(lw_combined, source == "wigley") %>% 
    select(source, season, comname, scientific_name, spec_class, 
           hare_group, fishery, catchsex, a, b, ln_a)
  
  
  # Do a second pass with the filter(lw_combined, source == "fishbase")
  # merge on common names only
  fb_trimmed <- filter(lw_combined, source == "fishbase") %>% 
    select(source, comname, scientific_name, spec_class, hare_group, fishery, a, b, ln_a)
  
  # First Pass - Wigley
  pass_1 <- trawl_spectra %>% 
    inner_join(w_trimmed)
  
  
  # Second Pass - Fishbase, for the stragglers if any
  pass_2 <- trawl_spectra %>% 
    filter(comname %notin% w_trimmed$comname) %>% 
    inner_join(fb_trimmed)
  
  
  
  # Join them with bind rows (implicitly drops things that don't have growth coefs)
  trawl_weights <- bind_rows(pass_1, pass_2) %>% 
    arrange(est_year, season) %>% 
    mutate(b = as.numeric(b),
           a = as.numeric(a),
           a = ifelse(is.na(a) & !is.na(ln_a), exp(ln_a), a),
           ln_a = ifelse(is.na(ln_a), log(a), ln_a),
           #ln_weight = (ln_a + b * log(length)),
           #single_weight_g = exp(ln_weight),                 # weight of individual in size class
           single_weight_g = exp(ln_a + b * log(length)),                 # weight of individual in size class
           sum_weight_g = single_weight_g * numlen_adj) %>%  # Individual weight * adjusted numlen
    drop_na(single_weight_g)
  
  
  # clean up environment
  rm(pass_1, pass_2, fb_trimmed, w_trimmed)
  
  
  
  
  
  
  
  ####__ 3. Use Coefficients to Re-calculate Biomass  ####
  
  # calculate total biomass again using weights from key 
  # make a key for the length weight coefficient sources
  trawl_weights <- trawl_weights %>%  
    arrange(est_year, season, comname, length) %>% 
    mutate(total_biomass = numlen_adj * single_weight_g,
           lw_group = str_c(comname, season, catchsex))
  
  
  
  return(trawl_weights)
}

####____________________####
####  Stratified Catch/Biomass  ####



####__ 1. Year & Stratum Counts  ####





#### Function to get stratified means at various group levels

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
      single_weight_g = mean(single_weight_g),
      ## Things that we want to total for the stratum  ##
      numlen_adj = sum(numlen_adj),       # Total count of individuals
      sum_weight_g = sum(sum_weight_g),   # Total Biomass of species at length for group strata
      mbiom_per_tow = sum_weight_g / mean(tow_count),  # total weight averaged by n-tows
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




