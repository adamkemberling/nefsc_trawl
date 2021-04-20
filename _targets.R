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


# # Load the build code and stratification function
# box_paths  <- research_access_paths(os.use = "unix")
# mills_path <- box_paths$mills
# res_path   <- box_paths$res
# nmfs_path  <- shared.path("unix", "RES_Data", "NMFS_Trawl")


# Size Spectra Build Functions
source(here("R/support/nefsc_ss_build_nodrop.R"))
source(here("R/support/maine_nh_trawl_build.R"))
source(here("R/support/sizeSpectra_support.R"))


####  Resource Paths  
box_paths   <- research_access_paths(os.use = "unix")
mills_path  <- box_paths$mills
res_path    <- box_paths$res


####  Targets Pipeline  ####

# Define target pipeline: Outlines high-level steps of the analysis
# Format is just a list of all the targets
# Order is not important, package sorts out connections for everything
list(
  
  #### CMIP Import  ####
  
  # Preparing Survdat Data
  tar_target(
    name = nefsc_clean,
    command = survdat_prep_nodrop(survdat_source = "most recent") ),
  tar_target(
    name = nefsc_survdat_lw,
    command = add_lw_info(nefsc_clean, cutoff = T) ),
  tar_target(
    name = nefsc_stratified,
    command = add_area_stratification(nefsc_survdat_lw, include_epu = F) ),
  
  # survdat biological data
  tar_target(
    name = nefsc_biological,
    command = survdat_prep_nodrop(survdat_source = "bio") ),


  ####  Prepping Size Spectrum Groups  #####
  
  # nmfs size spectra
  tar_target(
    name = nefsc_filtered_lens,
    command = min_length_cutoff(trawl_lens = nefsc_stratified, cutoff_cm = 1) ),
  tar_target(
    name = nefsc_databin,
    command = prep_wmin_wmax(lw_trawl_data = nefsc_filtered_lens)
  )
  
  
  
  # # Preparing Maine + NH Survey Data
  # # name matching super broken 4/8/2021
  # tar_target(
  #   name = menh_lens,
  #   command = load_menh_data(data_option = "length frequencies", os.use = "unix") ),
  # tar_target(
  #   name = menh_weights,
  #   command = add_lw_to_menh(menh_length_frequencies = menh_lens, os.use = "unix")
  # ),
  # 
  # # menh
  # tar_target(
  #   name = menh_filtered_lens,
  #   command = min_length_cutoff(trawl_lens = menh_weights, cutoff_cm = 1) ),
  # tar_target(
  #   name = menh_databin,
  #   command = prep_wmin_wmax(lw_trawl_data = menh_filtered_lens) )
  
)




#### Visualization of workflow
# targets::tar_visnetwork()

####  Processing the pipeline
# targets::tar_make()


