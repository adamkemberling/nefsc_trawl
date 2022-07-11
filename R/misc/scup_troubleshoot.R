# Scup test
# Why isn't it making it through add_lw_info:


# load libraries
#library(targets)
library(tidyverse)
library(gmRi)


# Load the data before this step:
survdat_clean_test <- gmri_survdat_prep(survdat = NULL, survdat_source = "most recent", "mojave")
survdat_clean <- survdat_clean_test %>% filter(comname == "scup")

# # Why is this Null?!
# survdat_clean <- tar_load(survdat_clean)




# Debugging:
# See what could be happening that lets scup data only pass through for some early years
# but not others. If systemic it likely affects other species:







#####  Add_lw_info  ####

#### 1. Match Species to LW Coefficients  ####

# Switch for mojave/Monterey users or other mac versions with CloudStorage folder
path_fun <- os_fun_switch(mac_os = "mojave")

####  Resource Paths
nmfs_path   <- path_fun(box_group = "RES_Data", subfolder = "NMFS_trawl")
lw_key_path <- paste0(nmfs_path, "length_weight_keys/fishbase_wigley_combined_key.csv")

# This table is a combined table of wigley and fishbase L-W coefficients
lw_combined <- readr::read_csv(lw_key_path, col_types = readr::cols())
lw_combined <- dplyr::mutate(lw_combined,
                             svspp = stringr::str_pad(svspp, 3, "left", "0"),
                             season = stringr::str_to_title(season))


# Do a priority pass with the dplyr::filter(lw_combined, source == "wigley)
# merge on comname, season, and catchsex
wigley_coefficients <- dplyr::filter(lw_combined, source == "wigley")
wigley_coefficients <- dplyr::select(wigley_coefficients,
                                     source, season, svspp, comname, scientific_name, spec_class,
                                     hare_group, fishery, catchsex, a, b, ln_a)


# Do a second pass with the dplyr::filter(lw_combined, source == "fishbase")
# merge on common names only
fishbase_coefficients <- dplyr::filter(lw_combined, source == "fishbase")
fishbase_coefficients <- dplyr::select(fishbase_coefficients,
                                       source, -svspp, comname, scientific_name, spec_class,
                                       hare_group, fishery, a, b, ln_a)



# Mismatched svspp codes
wigley_lookup <- function(x){unique(wigley_coefficients$svspp[which(wigley_coefficients$comname == x)])}
survdat_clean <- dplyr::mutate(
  .data = survdat_clean,
  svspp = if_else(comname == "scup", wigley_lookup("scup"), svspp))



# First Pass - Wigley
# Join just by svspp to account for name changes
pass_1 <- dplyr::select(survdat_clean, -comname)
pass_1 <- dplyr::inner_join(pass_1, wigley_coefficients)


# Second Pass - Fishbase, for the stragglers if any
pass_2 <- dplyr::filter(survdat_clean,
                        comname %not in% wigley_coefficients$comname,
                        svspp %not in% pass_1$svspp)
pass_2 <- dplyr::inner_join(pass_2, fishbase_coefficients)


# Join them with bind rows (implicitly drops things that don't have growth coefs)
trawl_weights <- dplyr::bind_rows(pass_1, pass_2)
trawl_weights <- dplyr::arrange(trawl_weights, est_year, season)
trawl_weights <- dplyr::mutate(trawl_weights,
                               b             = as.numeric(b),
                               a             = as.numeric(a),
                               a             = ifelse(is.na(a) & !is.na(ln_a), exp(ln_a), a),
                               ln_a          = ifelse(is.na(ln_a), log(a), ln_a),  # log of a used if ln_a is isn't already there (some fish just had ln_a reported)
                               llen          = log(length_cm),
                               ind_log_wt    = ln_a + (b * llen),
                               ind_weight_kg = exp(ind_log_wt),                # weight of an individual in size class
                               sum_weight_kg = ind_weight_kg * numlen_adj)     # Individual weight * adjusted numlen
trawl_weights <- tidyr::drop_na(trawl_weights, ind_weight_kg)
trawl_weights <- dplyr::select(trawl_weights, -ind_log_wt, -llen)


# clean up environment
rm(pass_1, pass_2, fishbase_coefficients, wigley_coefficients)




####  2. Use Coefficients to Re-calculate Biomass  ####

# calculate total biomass again using weights from key
# make a key for the length weight coefficient sources
survdat_weights <- dplyr::arrange(trawl_weights, est_year, season, comname, length_cm)
survdat_weights <- dplyr::mutate(survdat_weights, lw_group = stringr::str_c(comname, season, catchsex))






####  3. Filter Bad Fits  ####
# Drop comnames that don't align well with BIOMASS

# Length weight biomasses were checked against survdat$biomass using data prior to the vessel change
# 15% difference in either direction were flagged for removal
# code: github.com/adamkemberling/nefsc_trawl/R/qa_qc_reports/stratification_validation
# list updated : 8/27/2021
#
cutoff_15 <- c(
  "acadian redfish", "american plaice",
  "american shad",
  "atlantic angel shark",
  "atlantic cod",
  "atlantic croaker",
  "atlantic halibut",
  "atlantic herring",
  "atlantic mackerel",
  "atlantic sharpnose shark",
  "atlantic spadefish",
  "atlantic sturgeon",
  "atlantic thread herring",
  "atlantic wolffish",
  "barndoor skate",
  "black sea bass",
  "blackbelly rosefish",
  "blueback herring",
  "bluefish",
  "buckler dory",
  "bullnose ray",
  "butterfish",
  "chain dogfish",
  "clearnose skate",
  "cownose ray",
  "cunner",
  "cusk",
  "fawn cusk-eel",
  "fourspot flounder",
  "goosefish",
  "greater amberjack",
  "haddock",
  "little skate",
  "longhorn sculpin",
  "northern kingfish",
  "northern searobin",
  "ocean pout",
  "offshore hake",
  "pollock",
  "red hake",
  "rosette skate",
  "roughtail stingray",
  "round herring",
  "sand tiger",
  "sandbar shark",
  "scup",
  "sea raven",
  "silver hake",
  "smooth butterfly ray",
  "smooth dogfish",
  "smooth skate",
  "southern kingfish",
  "spanish mackerel",
  "spanish sardine",
  "spiny butterfly ray",
  "spiny dogfish",
  "spot",
  "spotted hake",
  "striped bass",
  "summer flounder",
  "thorny skate",
  "weakfish",
  "white hake",
  "windowpane flounder",
  "winter flounder",
  "winter skate",
  "witch flounder",
  "yellowtail flounder")

# these species were dropped, for the record:
cutoff_15_dropped <-c(
  "alewife", "atlantic torpedo", "bluntnose stingray",
  "southern stingray", "tautog", "vermillion snapper")

# Filter to use species that meet cutoff criteria
# source: 02_survdat_stratification_validation
if(cutoff == TRUE){
  survdat_weights <- dplyr::filter(survdat_weights, comname %in% cutoff_15)
}


