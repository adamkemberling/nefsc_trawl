# Title: Maine Neh Hampshire Trawl Cleanup
# 4/8/2021
# Adam A. Kemberling



####  Load Packages  ####
library(gmRi)
library(janitor)
library(tidyverse)






####  Loading Data  ####
load_menh_data <- function(data_option = c("length frequencies", "biological"), os.use = "unix"){
  
  
  # Paths
  res_path <- shared.path(os.use = os.use, group = "res data", folder = "")
  data_path <- switch (tolower(data_option),
    "length frequencies" = paste0(res_path, "Maine_NH_Trawl/shiny_length_freq.csv"),
    "biological"         = paste0(res_path, "Maine_NH_Trawl/full_me_dmr_biological.csv"),
    "expcatch"           = paste0(res_path, "Maine_NH_Trawl/full_me_dmr_expcatch.csv"),
    "lobster"            = paste0(res_path, "Maine_NH_Trawl/full_me_dmr_lobster_LF.csv")
  )
  
  
  # Length Frequencies from Data Portal 
  menh_data <- read_csv(data_path, guess_max = 1e5, col_types = cols()) %>% 
    clean_names()
  return(menh_data)
}


####  Cleanup Functions  ####

# Match length frequencies to lw relationship coefficients
add_lw_to_menh <- function(menh_length_frequencies = menh_lens, os.use = "unix"){
  
  # Paths
  nmfs_path <- shared.path(os.use = os.use, group = "res data", folder = "NMFS_trawl")
  
  # Format columns
  menh_data <- menh_length_frequencies %>% 
    rename(est_year = year,
           comname = common_name) %>% 
    mutate(catchsex = case_when(
      sex == "Male"    ~ 1,
      sex == "Female"  ~ 2,
      sex == "Unknown" ~ 0),
      comname = tolower(comname))
  
  
  # clean up names to be consistent with NMFS Common names
  # common pattern of fish adjective, example: cod atlantic
  # this needs to be fixed
  
  
  
  
  
 
  # Length Weight Lookup using Wigley and Fishbase
  lw_combined <- read_csv(paste0(nmfs_path, "length_weight_keys/fishbase_wigley_combined_key.csv"), col_types = cols()) %>% 
    mutate(svspp = str_pad(svspp, 3, "left", "0"),
           season = tolower(season)) %>% 
    mutate(season = str_to_title(season))

  

  ####__ 1. Pair LW Coefficients to Maine New Hampshire Data  ####

  
  
  # Do a priority pass with the filter(lw_combined, source == "wigley)
  # merge on comname, season, and catchsex
  wigley_coefficients <- filter(lw_combined, source == "wigley") %>% 
    select(source, season, svspp, comname, scientific_name, spec_class, 
           hare_group, fishery, catchsex, a, b, ln_a)
  
  
  # Do a second pass with the filter(lw_combined, source == "fishbase")
  # merge on common names only
  fishbase_coefficients <- filter(lw_combined, source == "fishbase") %>% 
    select(source, svspp, comname, scientific_name, spec_class, 
           hare_group, fishery, a, b, ln_a)  
  
  # First Pass - Wigley
  # Join just by svspp to account for name changes
  pass_1 <- menh_data %>% 
    inner_join(wigley_coefficients)
  
  # Second Pass - Fishbase, for the stragglers if any
  pass_2 <- menh_data %>% 
    filter(comname %not in% wigley_coefficients$comname) %>% 
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
      sum_weight_kg = ind_weight_kg * frequency) %>%     # Individual weight * adjusted numlen
    drop_na(ind_weight_kg) %>% 
    select(-ind_log_wt, -llen) %>%  
    arrange(est_year, season, comname, length) %>% 
    mutate(lw_group = str_c(comname, season, catchsex)) 
  
  
  # clean up environment
  rm(pass_1, pass_2, fishbase_coefficients, wigley_coefficients)
  
 
  # Get Length-Weight Biomasses
  return(trawl_weights)
  
  
}




