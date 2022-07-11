#### NMFS + Maine NH Trawl Survey Size Spectrum Analysis
# Targets Workflow


####  Load Packages
suppressPackageStartupMessages(suppressWarnings(suppressMessages(library(readr))))
suppressPackageStartupMessages(suppressWarnings(suppressMessages(library(tarchetypes))))
suppressPackageStartupMessages(suppressWarnings(suppressMessages(library(here))))
suppressPackageStartupMessages(suppressWarnings(suppressMessages(library(janitor))))
suppressPackageStartupMessages(suppressWarnings(suppressMessages(library(sf))))
suppressPackageStartupMessages(suppressWarnings(suppressMessages(library(tidyverse))))
suppressPackageStartupMessages(suppressWarnings(suppressMessages(library(gmRi))))
suppressPackageStartupMessages(suppressWarnings(suppressMessages(library(targets))))


####  Build code and stratification functions  ####

# Survdat cleanup functions
source("~/Documents/Repositories/gmRi/R/nefsc_groundfish_access.R")

# Size Spectra Build Functions
# source(here("R/support/maine_nh_trawl_build.R"))
source(here("R/support/sizeSpectra_support.R"))
source(here("R/support/temp_support.R"))




####_____________________________####

####__  Groundfish Data Preparation  __####

# Define target pipeline: Outlines high-level steps of the analysis
# Format is just a list of all the targets
# Order is not important, package sorts out connections for everything
list(
  
  ##### 1. NMFS Data Import  ####
  
  
  #####__ a. Full Survdat  ####
  # Preparing Survdat Data
  tar_target(targets_os, 
             command = "mojave"),
  tar_target(
    name = survdat_clean,
    command = gmri_survdat_prep(survdat = NULL, 
                                survdat_source = "most recent", 
                                mac_os = targets_os)),
  tar_target(
    name = survdat_lw,
    command = add_lw_info(survdat_clean, 
                          cutoff = T, 
                          mac_os = targets_os) ),
  tar_target(
    name = nefsc_stratified,
    command = add_area_stratification(survdat_lw, 
                                      include_epu = F, 
                                      mac_os = targets_os) ),
  
  #####__ b. Biological Data  ####
  # survdat biological data - for actual length relationships
  tar_target(
    name = survdat_biological,
    command = gmri_survdat_prep(survdat_source = "bio", 
                                mac_os = targets_os) ),
  tar_target(
    name = survdat_bio_lw,
    command = add_lw_info(survdat_biological, cutoff = T, 
                          mac_os = targets_os) %>% 
      mutate(Year = est_year,
             season = str_to_title(season))),
  

  #####  2. Size Spectrum Prep  #####
  
  # Prep the minimum size (wmin) and-
  # the maximum size (wmax) in grams for all the data
  tar_target(
    name = wmin_grams,
    command = prep_sizeSpectra_data(lw_trawl_data = nefsc_stratified)),
  
  # Apply a minimum size cutoff
  tar_target(
    name = nefsc_1g,
    command = min_weight_cutoff(nefsc_lw = wmin_grams, min_weight_g = 1)
  ),
  
  
  # Format the group labels and add discrete size groups for length/width
  tar_target(
    name = nefsc_1g_labelled,
    command = size_bin_formatting(nefsc_1g)
  ),
  
  

  
  
  ##### 4. log10 SS Slopes  ####
  
  # Assing the bin structure to the lw data
  tar_target(nefsc_1g_binned,
             assign_log10_bins(nefsc_1g_labelled)),
  
  

  
  # Run the different groupings through the slope estimation
  tar_target(nmfs_log10_slopes,
             log10_ss_all_groups(wmin_grams = nefsc_1g_binned,
                                 min_weight_g = 1)),
  
  
  
  
  ##### 5. sizeSpectra Exponents  ####
  
  
  
  # Run the mle calculation on stratified abundance
  tar_target(
    name = strat_total_mle_results,
    command = ss_slopes_all_groups(nefsc_1g_labelled, 
                                   min_weight_g = 1, 
                                   abundance_vals = "stratified")),
  
  # Join the MLE and Binned Resukts into a table
  tar_target(size_spectrum_indices,
             full_join(strat_total_mle_results, nmfs_log10_slopes,
                       by = c("group ID", "group_var", "Year", "season", "survey_area", "decade"))),
  
  
  
  
  
  ##### 6. Spectra - Species Sensitivity  ####
  
  # This section will repeat the estimation of log10 ss slopes,
  # repeating the steps with an omitted species to see the influence each has on 
  # the size spectrum estimate
  
  # # Create parallel groups that do not contain a species at each iteration
  # get the l1- slope info
  tar_target(species_ommission_dat,
             species_omit_spectra(start_dat = nefsc_1g_binned)),

  
  
  # Join the ommission data to the all species slopes to get the changes
  tar_target(year_region_only,
             filter(nmfs_log10_slopes, `group ID` == "single years * region")),
  
  
  # Code to run for this target
  tar_target(species_sensitivity_shifts,
             species_omit_changes(spec_omit_results = species_ommission_dat,
                                  all_spec_results = year_region_only)),

  
  
  ####_____________________________####
  
  ####__  Physical Drivers  ####
  
  #### 1. Temperature Data  ####
  tar_target(
    gom_oisst,
    oisst_access_timeseries(
      region_family = "nmfs trawl regions", 
      poly_name = "gulf of maine", 
      mac_os = targets_os )),
  tar_target(
    gb_oisst, 
    oisst_access_timeseries(
      region_family = "nmfs trawl regions", 
      poly_name = "georges bank", 
      mac_os = targets_os )),
  tar_target(
    mab_oisst,
    oisst_access_timeseries(
      region_family = "nmfs trawl regions", 
      poly_name = "mid atlantic bight", 
      mac_os = targets_os )),
  tar_target(
    sne_oisst,
    oisst_access_timeseries(
      region_family = "nmfs trawl regions", 
      poly_name = "southern new england", 
      mac_os = targets_os )),
  tar_target(
    inuse_strata_oisst,
    oisst_access_timeseries(
      region_family = "nmfs trawl regions", 
      poly_name = "inuse strata", 
      mac_os = targets_os )),
  
  # Make every daily timeseries into a yearly one
  tar_target(
    gom_yrly,
    make_yearly(gom_oisst) ),
  tar_target(
    gb_yrly,
    make_yearly(gb_oisst) ),
  tar_target(
    mab_yrly,
    make_yearly(mab_oisst) ),
  tar_target(
    sne_yrly,
    make_yearly(sne_oisst) ),
  tar_target(
    all_yrly,
    make_yearly(inuse_strata_oisst) ),
  
  # Put them together for plotting/merging
  tar_target(
    regional_oisst,
    bind_rows(
      list(
        "GoM" = gom_yrly,
        "GB"  = gb_yrly,
        "MAB" = mab_yrly,
        "SNE" = sne_yrly,
        "all" = all_yrly
      ), .id = "survey_area"
    )
  ),
  
  
  
  
  
  ####_____________________________####
  ####__ Body Size Changes  __####
  
  ##### 1. Mean Size Change  ####
  
  # "
  # Time series of mean length and weight will be constructed for the fish 
  # community using length-weight conversions based on observed lengths:
  # Year, season, taxonomic group, functional group, economic status
  # "
  
  # NOTE: not using BIO data because it is a subset of lengths
  # # Run the weighted mean sizes
  # tar_target(
  #   mean_individual_sizes,
  #   group_size_metrics(size_data = survdat_bio_lw, 
  #                      .group_cols = c("Year", "survey_area", "season", "spec_class", "fishery"))),
  # 
  
  # NOTE: Use all the data since all fish are measured, and we have LW weights
  # Run the main suite of groupings
  tar_target(
    mean_sizes_ss_groups,
    mean_sizes_all_groups(size_data = rename(as.data.frame(nefsc_stratified), Year = est_year),
                          min_weight_g = 0, 
                          abund_vals = "stratified")),
 
  # Run the size change for each species across years using stratified abundances
  tar_target(
    annual_individual_sizes,
    group_size_metrics(size_data = rename(as.data.frame(nefsc_stratified), Year = est_year),
                       .group_cols = c("comname", "Year", "season"),
                       abund_vals = "stratified")),
  
  
  
  #####  2. Size at Age Characteristics  ####
  
  # starts with "survdat_biological"
  # ref code: size_at_age_exploration.Rmd
  
  tar_target(vonbert_species_bio,
             select_vonbert_species(survdat_biological, rank_cutoff = 17))
  
  
  # tar_target(vonbert_growth_coef,
  #            estimate_vonbert_coef(vonbert_species_bio))
  # 
  
  
  
  
  
  
  ####_____________________________####
  ####__  Summary Tables  __####
  
  # ##### 1. Table of indices  ####
  # tar_target(
  #   group_community_indicators,
  #   full_join(mean_sizes_ss_groups, strat_total_mle_results)
  
  
  
  ##### 11. Chronological Clustering  ##
  # "
  # Based on these indicators, we will identify several multi-year “stanzas” 
  # with contrasting ecosystem conditions. These will be identified objectively
  # using chronological clustering. From prior work, we expect that the stanzas 
  # will roughly correspond to the 1980s, 1990s, 2000s, and 2010s 
  # (Pershing et al. 2010, Perretti et al. 2017). For simplicity, we will refer 
  # to these periods in our plan below, but the exact year ranges may change 
  # based on the analysis.
  # "
  # 
  
  
  
  

  
  
  
  
  
  
  
) 
####__  Close Pipeline  __####




#### Visualization of workflow
# targets::tar_visnetwork()

####  Processing the pipeline
# targets::tar_make()







####______####
####  Additional non-necessary steps  ####




