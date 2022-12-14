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


####  Targets Workflow Function Libraries  ####

# 1.
# NEFSC Survey Data "Survdat" cleanup functions (from gmRi package)
source("~/Documents/Repositories/gmRi/R/nefsc_groundfish_access.R")

# 2.
# Size Spectra Preparation and Analysis Functions
# source(here("R/support/maine_nh_trawl_build.R"))
source(here::here("R/support/sizeSpectra_support.R"))

# 3. 
# OISST Loading and Reshaping
source(here::here("R/support/temp_support.R"))

# 4.
# Macro-Functions for Condensing the visible pipeline
source(here("R/support/targets_macro_functions.R"))



####_____________________________####

####____  Groundfish Data Preparation:  __####

####_____________________________####

# Define target pipeline: Outlines high-level steps of the analysis
# Format is just a list of all the targets
# Order is not important, package sorts out connections for everything
list(
  
  
  
  ####  Configure Global Options:  ####
  
  
  # These are linked to the remaining analysis steps
  # with exception to analysis of average length/width, those are not filtered
  tar_target(analysis_options, 
             command = list(
               # Controls the data that reaches the analysis, and is used for context
               # Also Sets min/max for ISD exponent estimates
               min_input_weight_g = 2^0,
               max_input_weight_g = 10^4, # To pick a reasonable sounding limit
               # max_input_weight_g = 2^13, # To match log2 bins
               
               # Set/enforce the bin structure used for spectra analysis
               # These enforce what bins go to binned size spectra analysis
               # min and max set the left limits of the bin range: 
               # i.e. max_l10_bin of 3 = 14^3 to 10^4
               l10_bin_width = 1,
               min_l10_bin = 0,
               max_l10_bin = 3,
               
               # log2 bin limits - LBNbiom method
               log2_bin_width = 1,
               min_log2_bin = 0,
               max_log2_bin = 12
               )
             ),
  
  
  
  ##### 1. NMFS Data Import  ####
  
  
  # Pointer to data on Box,
  # Can be changed to trigger a full workflow reset
  tar_target(boxdata_location, 
             command = "cloudstorage"),
  
  
  
  ###### a. Catch Data  ####
  
  
  
  
  # Perform standard cleanup without LW and stratification for species that drop
  # This target is used to investigate data immediately after the base cleanup
  tar_target(
    name = survdat_clean,
    command = gmri_survdat_prep(survdat = NULL,
                                survdat_source = "most recent",
                                box_location = boxdata_location)),
  
  
  # Run all the import and tidying for catch data to use for the pipeline
  tar_target(catch_complete,
             command = import_and_tidy_catch(box_location = boxdata_location)),
  
  ###### b. Biological Data  ####
  
  # Import and tidy the biological dataset
  tar_target(bio_complete,
             command = import_and_tidy_bio(box_location = boxdata_location)),
  
  # survdat biological data - for actual length relationships
  tar_target(
    name = survdat_biological,
    command = gmri_survdat_prep(survdat_source = "bio",
                                box_location = boxdata_location) ),
  

  #####  2. Size Spectrum Prep  #####
  

  
  # # Prepare the data for use in size spectrum analysis:
  # # Set weight column units to grams
  # # impose a max/min weight
  # # set up the size bins - for plotting
  # # set a min/max weight a fish could be within 1cm growth increments for lme method
  # tar_target(
  #   catch_1g_labelled,
  #   command = size_spectrum_prep(
  #     catch_data = catch_complete, 
  #     min_weight_g = analysis_options[["min_input_weight_g"]], 
  #     max_weight_g = analysis_options[["max_input_weight_g"]], 
  #     bin_increment = analysis_options[["l10_bin_width"]])),
  

  
  # Assign the log size bins and apply min/max cutoffs
  tar_target(
    catch_log2_labelled,
    command = LBNbiom_prep(
      catch_data = catch_complete, 
      min_weight_g = analysis_options[["min_input_weight_g"]], # 2^0, 
      max_weight_g = analysis_options[["max_input_weight_g"]], # 2^13, 
      bin_increment = analysis_options[["log2_bin_width"]])),
  
  ##### 4a. log10 SS Slopes  ####
  

  # # Run the different groupings through the slope estimation
  # # lower end is >=, upper end is <
  # tar_target(warmem_log10_slopes,
  #            warmem_l10_estimates(
  #              # The input data for all the group estimates, labelled for 
  #              # both the binned and ISD analyses
  #              wmin_grams = catch_1g_labelled,
  #              
  #              # These two filter the input data
  #              min_weight_g = analysis_options[["min_input_weight_g"]],
  #              max_weight_g = analysis_options[["max_input_weight_g"]], 
  #              
  #              # These two enforce the bins used to estimate the spectra
  #              min_l10_bin = analysis_options[["min_l10_bin"]], 
  #              max_l10_bin = analysis_options[["max_l10_bin"]], 
  #              bin_increment = analysis_options[["l10_bin_width"]])),
  
  
  
  
  
  ##### 4b. LBNbiom Method Spectra  ####
  # See Edwards et al. 2016
  # Method LBNbiom
  
  
  
  
  # Run the spectrum fitting using the log2 binning and LBNbiom methodology
  tar_target(warmem_log2_slopes,
             warmem_log2_estimates(
               # The input data for all the group estimates, labelled for 
               # both the binned and ISD analyses
               wmin_grams = catch_log2_labelled,
               
               # These two filter the input data
               min_weight_g = analysis_options[["min_input_weight_g"]],
               max_weight_g = analysis_options[["max_input_weight_g"]], 
               
               # These two enforce the bins used to estimate the spectra
               min_log2_bin = analysis_options[["min_log2_bin"]],
               max_log2_bin = analysis_options[["min_log2_bin"]], 
               bin_increment = analysis_options[["log2_bin_width"]])),
  
  
  
  
  ##### 5. sizeSpectra Exponents  ####
  
  # These should match the max size available to the size spectra bins
  # or be close to it
  # its obvious from the composition qmd data those largest sizes just aren't sampled in many groups
  # This can be achieved with the max size filtering, or with a control in this function
  
  # Run the individual size distribution calculation on stratified abundances
  tar_target(
    name = warmem_isd_results,
    command = warmem_isd_estimates(
      wmin_grams = catch_log2_labelled, 
      min_weight_g = analysis_options[["min_input_weight_g"]],
      max_weight_g = analysis_options[["max_input_weight_g"]],
      isd_xmin = analysis_options[["min_input_weight_g"]],
      isd_xmax = analysis_options[["max_input_weight_g"]],
      abundance_vals = "stratified")),
  
  
  
  # Join the MLE and Binned Results into a table
  tar_target(size_spectrum_indices,
             full_join(warmem_isd_results, warmem_log2_slopes,
                       by = c("group ID", "group_var", "Year", "season", "survey_area", "decade"))),
  

  
  
  ####_____________________________####
  
  ####______  Physical Drivers  _______####
  
  # ##### 1. Temperature Data  ####
  tar_target(
    name = regional_oisst,
    command = process_regional_sst(box_location = boxdata_location)
  ),
  
  
  
  
  
  ####_____________________________####
  ####_____ Body Size Changes  ______####
  
  ##### 1. Mean Size Change  ####
  
  # "
  # Time series of mean length and weight will be constructed for the fish 
  # community using length-weight conversions based on observed lengths:
  # Year, season, taxonomic group, functional group, economic status
  # "
  
  # NOTE: Use all the data since all fish are measured, and we have LW weights
  # Run the main suite of groupings
  tar_target(
    mean_sizes_ss_groups,
    mean_sizes_all_groups(
      size_data = rename(as.data.frame(catch_complete), Year = est_year),
      min_weight_g = 0, 
      abund_vals = "numlen_adj")
      #abund_vals = "stratified")
    ),
 
  # Run the size change for each species across years using stratified abundances
  tar_target(
    annual_individual_sizes,
    group_size_metrics(
      size_data = rename(as.data.frame(catch_complete), Year = est_year),
      .group_cols = c("comname", "Year", "season"),
      abund_vals = "numlen_adj")
      #abund_vals = "stratified")
    ),
  
  
  
  #####  2. Size at Age Characteristics  ####
  
  # starts with "survdat_biological"
  # ref code: size_at_age_exploration.Rmd
  
  tar_target(vonbert_species_bio,
             select_vonbert_species(bio_complete, rank_cutoff = 17))
  
  
  # tar_target(vonbert_growth_coef,
  #            estimate_vonbert_coef(vonbert_species_bio))
  # 
  
  
  
  
  
  

  
  
  
  
  
  
  
) ####_____________________________####
####______  END PIPELINE  __________####




#### Visualization of workflow
# targets::tar_visnetwork()

####  Processing the pipeline
# targets::tar_make()











####_____________________________####
####__ Additional phased-out steps  __####



#### Condensed into Macro Steps:  ####



##### 1. Catch Data  ####

# # Perform standard cleanup
# tar_target(
#   name = survdat_clean,
#   command = gmri_survdat_prep(survdat = NULL,
#                               survdat_source = "most recent",
#                               box_location = boxdata_location)),
# 
# # Assign Length-Weight Relationships
# tar_target(
#   name = survdat_lw,
#   command = add_lw_info(survdat_clean, 
#                         cutoff = T, 
#                         box_location = boxdata_location) ),
# 
# # Perform Area Stratification
# tar_target(
#   name = catch_stratified,
#   command = add_area_stratification(survdat_lw, 
#                                     include_epu = F, 
#                                     box_location = boxdata_location) ),
# 
# # Fill in Gaps in functional groups
# tar_target(
#   name = catch_complete,
#   command = fill_func_groups(catch_stratified)
# ),

##### 2. Bio Data  ####


# # survdat biological data - for actual length relationships
# tar_target(
#   name = survdat_biological,
#   command = gmri_survdat_prep(survdat_source = "bio", 
#                               box_location = boxdata_location) ),
# tar_target(
#   name = survdat_bio_lw,
#   command = add_lw_info(survdat_biological, cutoff = T, 
#                         box_location = boxdata_location) %>% 
#     mutate(Year = est_year,
#            season = str_to_title(season))),

##### 3.  Spectra Prep  ####
# # Prep the minimum size (wmin) and-
# # the maximum size (wmax) in grams for all the data based on 
# # length in cm, and length+1 in cm using LW relationships
# tar_target(
#   name = wmin_grams,
#   command = prep_sizeSpectra_data(lw_trawl_data = catch_complete)),
# 
# # Apply a minimum size cutoff
# tar_target(
#   name = catch_mincut,
#   command = min_weight_cutoff(catch_lw = wmin_grams, min_weight_g = 1)),
# 
# 
# # Apply a maximum weight cutoff
# tar_target(
#   name = catch_maxcut,
#   command = max_weight_cutoff(catch_lw = catch_mincut, max_weight_g = 10^5)),
# 
# 
# # Format the group labels and add discrete size groups for length/width
# tar_target(
#   name = catch_1g_labelled,
#   command = size_bin_formatting(catch_maxcut)),
# 
# # Assign the bin structure to the lw data 
# tar_target(catch_1g_binned,
#            assign_log10_bins(catch_1g_labelled)),





##### 4. Temperature Data  ####
# tar_target(
#   gom_oisst,
#   oisst_access_timeseries(
#     region_family = "nmfs trawl regions", 
#     poly_name = "gulf of maine", 
#     box_location = boxdata_location )),
# tar_target(
#   gb_oisst, 
#   oisst_access_timeseries(
#     region_family = "nmfs trawl regions", 
#     poly_name = "georges bank", 
#     box_location = boxdata_location )),
# tar_target(
#   mab_oisst,
#   oisst_access_timeseries(
#     region_family = "nmfs trawl regions", 
#     poly_name = "mid atlantic bight", 
#     box_location = boxdata_location )),
# tar_target(
#   sne_oisst,
#   oisst_access_timeseries(
#     region_family = "nmfs trawl regions", 
#     poly_name = "southern new england", 
#     box_location = boxdata_location )),
# tar_target(
#   inuse_strata_oisst,
#   oisst_access_timeseries(
#     region_family = "nmfs trawl regions", 
#     poly_name = "inuse strata", 
#     box_location = boxdata_location )),
# 
# # Make every daily timeseries into a yearly one
# tar_target(
#   gom_yrly,
#   make_yearly(gom_oisst) ),
# tar_target(
#   gb_yrly,
#   make_yearly(gb_oisst) ),
# tar_target(
#   mab_yrly,
#   make_yearly(mab_oisst) ),
# tar_target(
#   sne_yrly,
#   make_yearly(sne_oisst) ),
# tar_target(
#   all_yrly,
#   make_yearly(inuse_strata_oisst) ),
# 
# # Put them together for plotting/merging
# tar_target(
#   regional_oisst,
#   bind_rows(
#     list(
#       "GoM" = gom_yrly,
#       "GB"  = gb_yrly,
#       "MAB" = mab_yrly,
#       "SNE" = sne_yrly,
#       "all" = all_yrly
#     ), .id = "survey_area"
#   )
# ),





# ##### 6. Spectra - Species Sensitivity  ####
# 
# # This section will repeat the estimation of log10 ss slopes,
# # repeating the steps with an omitted species to see the influence each has on
# # the size spectrum estimate
# 
# # # Create parallel groups that do not contain a species at each iteration
# # get the l1- slope info
# tar_target(species_ommission_dat,
#            species_omit_spectra(
#              start_dat = catch_1g_labelled,
#              .group_cols = c("Year", "survey_area"),
#              # To filter the input data
#              min_weight_g = analysis_options[["min_input_weight_g"]],
#              # These two enforce the bins used to estimate the spectra
#              min_l10_bin = analysis_options[["min_l10_bin"]], 
#              max_l10_bin = analysis_options[["max_l10_bin"]], 
#              bin_increment = analysis_options[["l10_bin_width"]])),
# 
# 
# 
# # Join the ommission data to the all species slopes to get the changes
# tar_target(year_region_only,
#            filter(warmem_log10_slopes, `group ID` == "single years * region")),
# 
# 
# # Code to run for this target
# tar_target(species_sensitivity_shifts,
#            species_omit_changes(spec_omit_results = species_ommission_dat,
#                                 all_spec_results = year_region_only)),

####____________________________####


####__  Summary Tables  ####

# ##### 1. Table of indices  ####
# tar_target(
#   group_community_indicators,
#   full_join(mean_sizes_ss_groups, warmem_isd_results)



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

