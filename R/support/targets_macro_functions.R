# Targets Workflow Macro-Functions
# Goal: Combine multiple discrete processing steps into singular functions:


# Packages
library(ggpmisc)
suppressWarnings(library(scales))
suppressWarnings(library(patchwork))
library(sizeSpectra)
library(tidyverse)


# Building Block Functions
source(here::here("R/support/sizeSpectra_support.R"))


# Macro 1: survdat catch and survdat bio
# Load, add_lw, supplement functional groups
import_and_tidy_catch <- function(box_location){
  catch_data <- gmri_survdat_prep(survdat = NULL, 
                    survdat_source = "most recent", 
                    box_location = box_location) %>% 
    add_lw_info(cutoff = T, 
                box_location = box_location) %>% 
    add_area_stratification(include_epu = F, 
                            box_location = box_location) %>% 
    fill_func_groups() 
  return(catch_data)
}
    

# Same steps (minus area stratification) for biological data
import_and_tidy_bio <- function(box_location){
  bio_data <- gmri_survdat_prep(survdat = NULL, 
                    survdat_source = "bio", 
                    box_location = box_location) %>% 
    add_lw_info(cutoff = T, 
                box_location = box_location)%>% 
    mutate(Year = est_year,
           season = str_to_title(season)) %>% 
    fill_func_groups()  
  return(bio_data)
}





# Macro 2: Prepare Size Spectra
# area stratification, size truncation, set units, 
# prepare upper/lower weight limits for 1cm length increments
size_spectrum_prep <- function(
    catch_data, 
    min_weight_g = 1, 
    max_weight_g = 10^5,
    bin_increment = 1){
  spectra_input <- prep_sizeSpectra_data(lw_trawl_data = catch_data) %>% 
    min_weight_cutoff(catch_lw = ., min_weight_g = min_weight_g) %>% 
    max_weight_cutoff(catch_lw = ., max_weight_g = max_weight_g) %>% 
    #size_bin_formatting(catch_1g = .) %>%  # purely aesthetic binning
    assign_log10_bins(., l10_increment = bin_increment)
  return(spectra_input)
}







# Macro 2: Prepare Size Spectra
# area stratification, size truncation, set units, 
# prepare upper/lower weight limits for 1cm length increments
LBNbiom_prep <- function(
    catch_data, 
    min_weight_g = 1, 
    max_weight_g = 2^13,
    bin_increment = 1){
  spectra_input <- prep_sizeSpectra_data(lw_trawl_data = catch_data) %>% 
    min_weight_cutoff(catch_lw = ., min_weight_g = min_weight_g) %>% 
    max_weight_cutoff(catch_lw = ., max_weight_g = max_weight_g) %>% 
    assign_log2_bins(., log2_increment = bin_increment)
  return(spectra_input)
}




# Macro 3: Estimate Spectra for Groups
# Split on the proper group columns
# Run the size spectrum function
# Add in the rest of the groups not explicitly used as "all"
# Compile them together in one table
process_size_spectra <- function(){}
  
  


# Macro 4: Process Regional SST
process_regional_sst <- function(box_location){
  
  
  # Load regional timeseries:
  gom_oisst <- oisst_access_timeseries(
      region_family = "nmfs trawl regions", 
      poly_name = "gulf of maine", 
      box_location = box_location )
  gb_oisst <-  oisst_access_timeseries(
      region_family = "nmfs trawl regions", 
      poly_name = "georges bank", 
      box_location = box_location )
  mab_oisst <- oisst_access_timeseries(
      region_family = "nmfs trawl regions", 
      poly_name = "mid atlantic bight", 
      box_location = box_location )
  sne_oisst <- oisst_access_timeseries(
      region_family = "nmfs trawl regions", 
      poly_name = "southern new england", 
      box_location = box_location )
  inuse_strata_oisst <- oisst_access_timeseries(
      region_family = "nmfs trawl regions", 
      poly_name = "inuse strata", 
      box_location = box_location )
  
  # Make every daily timeseries into a yearly one
  gom_yrly <- make_yearly(gom_oisst)
  gb_yrly <- make_yearly(gb_oisst)
  mab_yrly <- make_yearly(mab_oisst)
  sne_yrly <- make_yearly(sne_oisst)
  all_yrly <- make_yearly(inuse_strata_oisst)
  
  # Put them together for plotting/merging
  regional_oisst <- bind_rows(
      list(
        "GoM" = gom_yrly,
        "GB"  = gb_yrly,
        "MAB" = mab_yrly,
        "SNE" = sne_yrly,
        "all" = all_yrly
      ), .id = "survey_area")
  
  
  # Return the regional sst yearly averages
  return(regional_oisst)
  
}
