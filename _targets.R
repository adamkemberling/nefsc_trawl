#### NMFS + Maine NH Trawl Survey Size Spectrum Analysis
# Targets Workflow


####  Load Packages
library(targets)
library(tarchetypes)
library(here)
library(janitor)
library(gmRi)
library(patchwork)
library(rnaturalearth)
library(sf)
library(tidyverse)

####  Resource Paths  
oisst_path <- box_path(box_group = "RES_Data", subfolder = "/OISST/oisst_mainstays")

####  Build code and stratification functions  ####

# Size Spectra Build Functions
# source(here("R/support/nefsc_ss_build_nodrop.R"))
# source(here("R/support/maine_nh_trawl_build.R"))
source(here("R/support/sizeSpectra_support.R"))
source(here("R/support/temp_support.R"))




####__  Targets Pipeline  __####

# Define target pipeline: Outlines high-level steps of the analysis
# Format is just a list of all the targets
# Order is not important, package sorts out connections for everything
list(
  
  ##### 1. NMFS Data Import  ####
  
  
  #####__ a. Full Survdat  ####
  # Preparing Survdat Data
  tar_target(
    name = nefsc_clean,
    command = gmri_survdat_prep(survdat_source = "most recent") ),
  tar_target(
    name = nefsc_survdat_lw,
    command = add_lw_info(nefsc_clean, cutoff = T) ),
  tar_target(
    name = nefsc_stratified,
    command = add_area_stratification(nefsc_survdat_lw, include_epu = F) ),
  
  #####__ b. Biological Data  ####
  # survdat biological data - for actual length relationships
  tar_target(
    name = nefsc_biological,
    command = gmri_survdat_prep(survdat_source = "bio") ),
  tar_target(
    name = nefsc_bio_lw,
    command = add_lw_info(nefsc_biological, cutoff = T) %>% 
      mutate(Year = est_year,
             season = str_to_title(season))),
  
  # tar_target(
  #   name = bio_stratified,
  #   command = add_area_stratification(nefsc_bio_lw, include_epu = F) %>% 
  #     rename(Year = est_year) %>% 
  #     mutate(season = str_to_title(season))),

  #####  2. Prep Size Spectrum Groups  #####
  
  # Prep wmin and wmax for all the data
  tar_target(
    name = wmin_grams,
    command = prep_wmin_wmax(lw_trawl_data = nefsc_stratified)),
  tar_target(
    name = nefsc_1g,
    command = min_weight_cutoff(nefsc_lw = wmin_grams, min_weight_g = 1)
  ),
  
  
  
  
  
  
  ##### 4. Manual log10 Size Spectra  ####
  
  # tar_target(nefsc_l10_assigned,
  #            assign_log10_bins(wmin_grams = nefsc_1g)),
  tar_target(nmfs_log10_slopes,
             log10_ss_all_groups(wmin_grams = nefsc_1g,
                                 min_weight_g = 1)),
  
  
  
  
  ##### 5. Results for WARMEM Groups  ####
  
  # Run the mle calculation on survey abundance
  tar_target(
    name = nmfs_group_mle_ss,
    command = ss_slopes_all_groups(nefsc_1g, 
                                  min_weight_g = 1, 
                                  abundance_vals = "observed")),
  # Run the mle calculation on stratified abundance
  tar_target(
    name = nmfs_stratified_mle_ss,
    command = ss_slopes_all_groups(nefsc_1g, 
                                  min_weight_g = 1, 
                                  abundance_vals = "stratified")),
  tar_target(size_spectrum_indices,
             full_join(nmfs_stratified_mle_ss, nmfs_log10_slopes)),
  
  
  
  
  
  
  ##### 6. Temperature Data  ####
  tar_target(
    gom_oisst,
    oisst_access_timeseries(oisst_path = oisst_path, 
                            region_family = "nmfs trawl regions", 
                            poly_name = "gulf of maine") ),
  tar_target(
    gb_oisst, 
    oisst_access_timeseries(oisst_path = oisst_path, 
                            region_family = "nmfs trawl regions", 
                            poly_name = "georges bank") ),
  tar_target(
    mab_oisst,
    oisst_access_timeseries(oisst_path = oisst_path, 
                            region_family = "nmfs trawl regions", 
                            poly_name = "mid atlantic bight") ),
  tar_target(
    sne_oisst,
    oisst_access_timeseries(oisst_path = oisst_path, 
                            region_family = "nmfs trawl regions", 
                            poly_name = "southern new england") ),
  tar_target(
    inuse_strata_oisst,
    oisst_access_timeseries(oisst_path = oisst_path, 
                            region_family = "nmfs trawl regions", 
                            poly_name = "inuse strata") ),
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
  
  
  
  
  
  ##### 7. Mean Sizes  ####
  
  # "
  # Time series of mean length and weight will be constructed for the fish 
  # community using length-weight conversions based on observed lengths:
  # Year, season, taxonomic group, functional group, economic status
  # "
  
  # Run the weighted mean sizes
  tar_target(
    mean_individual_sizes,
    group_size_metrics(nefsc_bio = nefsc_bio_lw, 
                       .group_cols = c("Year", "survey_area", "season", "spec_class", "fishery"))),
  
  tar_target(
    mean_sizes_ss_groups,
    mean_sizes_all_groups(nefsc_bio = nefsc_bio_lw, 
                          min_weight_g = 0))
  
  
  # ##### 8. Assemble Table of indices  ####
  # tar_target(
  #   group_community_indicators,
  #   full_join(mean_sizes_ss_groups, nmfs_stratified_mle_ss)
  
  
  
  ##### 9. Chronological Clustering  ####
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


##### 3. Process LME Group SS  #### 


# # Individual Years
# tar_target(
#   name = stratified_annual_ss,
#   command = group_mle_slope_estimate(wmin_grams = nefsc_1g,
#                                     min_weight_g = 1, 
#                                     abundance_vals = "stratified",
#                                     .group_cols = c("Year"))),
# 
# # Specific Areas
# tar_target(
#   name = stratified_area_ss,
#   command = group_mle_slope_estimate(wmin_grams = nefsc_1g,
#                                     min_weight_g = 1, 
#                                     abundance_vals = "stratified",
#                                     .group_cols = c("survey_area"))),
# # Seasons
# tar_target(
#   name = stratified_season_ss,
#   command = group_mle_slope_estimate(wmin_grams = nefsc_1g,
#                                     min_weight_g = 1, 
#                                     abundance_vals = "stratified",
#                                     .group_cols = c("season"))),
# # Years and seasons
# tar_target(
#   name = stratified_yr_season_ss,
#   command = group_mle_slope_estimate(wmin_grams = nefsc_1g,
#                                     min_weight_g = 1,
#                                     abundance_vals = "stratified",
#                                     .group_cols = c("Year", "season"))),
# # Years and areas
# tar_target(
#   name = stratified_yr_area_ss,
#   command = group_mle_slope_estimate(wmin_grams = nefsc_1g,
#                                     min_weight_g = 1,
#                                     abundance_vals = "stratified",
#                                     .group_cols = c("Year", "survey_area"))),
# # Year Season Area
# tar_target(
#   name = stratified_decade_area_ss,
#   command = group_mle_slope_estimate(wmin_grams = nefsc_1g,
#                                     min_weight_g = 1,
#                                     abundance_vals = "stratified",
#                                     .group_cols = c("Year", "survey_area", "season"))),


# # Preparing Maine + NH Survey Data
# # name matching super broken 4/8/2021
# tar_target(
#   name = menh_lens,
#   command = load_menh_data(data_option = "length frequencies") ),
# tar_target(
#   name = menh_weights,
#   command = add_lw_to_menh(menh_length_frequencies = menh_lens) ),
# 
# # menh
# tar_target(
#   name = menh_filtered_lens,
#   command = min_length_cutoff(trawl_lens = menh_weights, cutoff_cm = 1) ),
# tar_target(
#   name = menh_databin,
#   command = prep_wmin_wmax(lw_trawl_data = menh_filtered_lens) )


# # Global plot limits for individual size distributions
# tar_target(name = isd_plot_lims,
#            command = c( 0, max(nefsc_1g$wmax) )),


# Size Spectrum Group Comparisons
# tar_target(factor_groupings,
#            command = list(
#              "All Data"            = c("All Data"),
#              "Annual slopes"       = c("Year"),
#              "Years and seasons"   = c("Year", "season"),
#              "Years and regions"   = c("Year", "survey_area"),
#              "Seasons Overall"     = c("season"),
#              "Seasons and Regions" = c("season", "survey_area"),
#              "Regions Overall"     = c("survey_area"),
#              "Years, Regions, and Seasons" = c("Year", "season", "survey_area") )),
# 
# tar_target(group_sizespectra_slopes, w + x, pattern = map(w, x)),



