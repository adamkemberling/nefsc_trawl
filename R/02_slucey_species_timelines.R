# Single Species Dig into comparing survdat pulls
# Goal: Take Sean's Data
# Create some display functions for annual abundance, biomass timelines
# Quick Checks of What abundance or biomass look like across sources:


####  Load Packages  ####
library(here)
library(janitor)
library(gmRi)
library(patchwork)
library(tidyverse)


# cleanup function
source(here("R/01_nefsc_ss_build_nodrop.R"))

# Load the build code and stratification function
box_paths  <- research_access_paths(os.use = "unix")
mills_path <- box_paths$mills
res_path   <- box_paths$res
nmfs_path  <- shared.path("unix", "RES_Data", "NMFS_Trawl")

#### Set theme  ####
theme_set(theme_bw() + theme(axis.text.y = element_text(size = 11)))

# function to check abundance/biomass of raw data
check_raw <- function(raw_survdat, spec_name){
  raw_survdat %>% 
    filter(comname == toupper(spec_name)) %>% 
    distinct(station, year, svvessel, season, comname, abundance, biomass) %>% 
    group_by(year) %>% 
    summarise(n_stations = n(),
              tot_abund = sum(abundance, na.rm = T),
              tot_bio = sum(biomass, na.rm = T)) %>% 
    ungroup()
}





####  Load Data  ####

# Species of highest concern
species_check <- read_csv(here("data/andrew_species/Assesmentfishspecies.csv"), 
                          col_types = cols()) %>% 
  clean_names() %>% 
  mutate(svspp = str_pad(svspp, 3, pad = "0", side = "left"),
         comname = tolower(comname),
         species = str_to_title(species)) %>% 
  arrange(svspp)



# 2019 Data used in 2020 paper
load(paste0(nmfs_path, "Survdat_Nye_allseason.RData"))
survdat_nye_raw  <- survdat %>% clean_names()


# 2020 data received last august
load(paste0(nmfs_path, "Survdat_Nye_Aug 2020.RData"))
survdat_20_raw <- clean_names(survdat) 


# Data we just received in 2021 with errors located and corrected w/ additional columns
load(paste0(nmfs_path, "2021_survdat/NEFSC_BTS_all_seasons_03032021.RData"))
survdat_21_raw <- survey$survdat %>% clean_names()

# 2021 with fixes, but also bio changes
load(paste0(nmfs_path, "2021_survdat/NEFSC_BTS_2021_bio_03192021.RData"))
survdat_21_bio_raw <- survey$survdat %>% clean_names()


# Remove base files
rm(survdat, survey)

####  Clean up  ####

survdat_nye <- survdat_prep_nodrop(survdat = survdat_nye_raw) %>% 
  mutate(survdat_source = "survdat_nye")

survdat_20 <- survdat_prep_nodrop(survdat = survdat_20_raw) %>% 
  mutate(survdat_source = "survdat_nye2020")

survdat_21 <- survdat_prep_nodrop(survdat = survdat_21_raw) %>% 
  mutate(survdat_source = "survdat_2021") 

survdat_21_bio <- survdat_prep_nodrop(survdat = survdat_21_bio_raw) %>% 
  mutate(survdat_source = "survdat_2021_bio") 





####  Plotting Functions  ####

annual_plots <- function(survdat1, survdat2, spec_name, source1, source2){
  
  
  # Make a key of all years and the data source to merge into and preserve NA's
  yr_df <- data.frame(
    "est_year" = sort(unique( c(unique(survdat1$est_year), unique(survdat2$est_year) ))))  %>% 
      mutate(s1 = source1, 
             s2 = source2) %>% 
      pivot_longer(names_to = "short_source", values_to = "Data Source", cols = c("s1", "s2")) 
  
  
  
  # get annual summary
  comp_data <- map(list(survdat1, survdat2), function(sdat){ 
    
    
    # Get Annual summaries
    comp_data <- sdat %>% 
      filter(comname == spec_name) %>% 
      distinct(id, est_year, season, svvessel, 
               comname, svspp, catchsex, abundance, biom_adj) %>% 
      group_by(est_year) %>% 
      summarise(
        n_stations_present = n(),
        tot_bio = sum(biom_adj, na.rm = T),
        tot_abund = sum(abundance, na.rm = T),
        .groups = "keep") %>% 
      ungroup()
  }) %>% setNames(c(source1, source2)) %>% 
    bind_rows(.id = "Data Source")
  
  # Merge to key to preserve gaps
  comp_data <- left_join(yr_df, comp_data, by = c("Data Source", "est_year")) %>% 
    mutate(yr_numeric = as.numeric(as.character(est_year)))
      
  # abundance plot
  abund_plot <- ggplot() +
    geom_point(data = comp_data, aes(x = yr_numeric, y = tot_abund, color = `Data Source`), show.legend = F) +
    geom_point(data = comp_data, aes(x = yr_numeric, y = tot_abund, color = `Data Source`), show.legend = F) +
    geom_line(data = comp_data,  aes(x = yr_numeric, y = tot_abund, color = `Data Source`), show.legend = F) +
    geom_line(data = comp_data,  aes(x = yr_numeric, y = tot_abund, color = `Data Source`), show.legend = F) +
    scale_color_gmri() +
    labs(x = "", y = "Annual Abundance", color = "", title = str_to_title(spec_name)) +
    scale_y_continuous(labels = scales::comma_format()) +
    theme(axis.text.x = element_blank())
  
  # biomass plot
  bio_plot <- ggplot() +
    geom_point(data = comp_data, aes(x = yr_numeric, y = tot_bio, color = `Data Source`)) +
    geom_point(data = comp_data, aes(x = yr_numeric, y = tot_bio, color = `Data Source`)) +
    geom_line(data = comp_data,  aes(x = yr_numeric, y = tot_bio, color = `Data Source`)) +
    geom_line(data = comp_data,  aes(x = yr_numeric, y = tot_bio, color = `Data Source`)) +
    scale_color_gmri() +
    labs(x = "", y = "Annual Biomass", color = "") +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    scale_y_continuous(labels = scales::comma_format()) +
    theme(legend.position = "bottom")
  
  return(abund_plot / bio_plot)
      
  
}







####  Digging Into Plots  ####



# Manual Check
annual_plots(survdat1 = survdat_nye, 
             survdat2 = survdat_21, 
             spec_name = "longfin squid", 
             source1 = "2019", 
             source2 = "2021 Mar 3rd")

annual_plots(survdat1 = survdat_20, 
             survdat2 = survdat_21_bio, 
             spec_name = "longfin squid", 
             source1 = "2020", 
             source2 = "2021  BIO")




#### Checking them all  ####

# Set list of species
andrew_species <- species_check$comname %>% setNames(species_check$comname)
andrew_species <- sort(andrew_species)

# loop through
all_check <- map(andrew_species, function(andrew_spec){
  annual_plots(survdat1 = survdat_nye, 
               survdat2 = survdat_21, 
               spec_name = andrew_spec, 
               source1 = "2019", 
               source2 = "2021") }) 


# check important ones
all_check$`acadian redfish`
all_check$`atlantic sturgeon`
all_check$`atlantic cod`
all_check$`haddock`
all_check$`butterfish`
all_check$`alewife`
all_check$`longfin squid`
all_check$`striped bass`
all_check$`spiny dogfish`
all_check$`blue crab`
all_check$


# are we losing it duting filtering?
check_raw(survdat_21_raw, "spiny dogfish")
check_raw(survdat_21_bio_raw, "spiny dogfish")
check_raw(survdat_21_raw, "longfin squid")
check_raw(survdat_21_bio_raw, "longfin squid")



####  Checking Tables

# Annual tables - stations by stratum season vessel etc.
annual_tables <- function(survdat1, survdat2, source1, source2, spec_name){
  
  # slightly different group_by here
  comp_data <- map(list(survdat1, survdat2), function(sdat){ 
    sdat %>% 
      filter(comname == spec_name) %>% 
      distinct(id, est_year, season, svvessel, 
               comname, svspp, catchsex, abundance, biom_adj) %>% 
      group_by(est_year, season, svvessel) %>% 
      summarise(
        n_stations_present = n(),
        tot_bio = sum(biom_adj, na.rm = T),
        tot_abund = sum(abundance, na.rm = T),
        .groups = "keep") %>% 
      ungroup()
  }) %>% setNames(c(source1, source2)) %>% 
    bind_rows(.id = "Data Source") 
  
  
  # Expand grid to get at NA's
  comp_data <- comp_data %>% 
    expand(`Data Source`, est_year, season, svvessel) %>% 
    left_join(comp_data) %>% 
    mutate(yr_numeric = as.numeric(as.character(est_year)))
 
  
  return(comp_data)
  
}


####  Dig Into Tables  ####


# Plotting how many stations the species was caught at
species_presence <- function(my_species){
  annual_tables(survdat1 = survdat_nye, 
                survdat2 = survdat_21_bio, 
                source1 = "2019", 
                source2 = "2021", 
                spec_name = my_species)  %>% 
    arrange(est_year, season, svvessel) %>% 
    ggplot(aes(yr_numeric, n_stations_present, color = `Data Source`)) +
      geom_point() +
      geom_line() +
      facet_grid(season~svvessel) + 
      labs(x = "", y = "Stations Present", title = str_to_title(my_species)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
            legend.position = "bottom")
}


species_presence("spiny dogfish")
