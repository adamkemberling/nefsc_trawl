####  Landings Data EDA  ####
# 6/7/2022

####  Packages  ####
library(readxl)
library(here)
library(tidyverse)
library(gmRi)
library(janitor)


####  Data  ####

# Organized by port
by_port <- read_xlsx(
  path = here("data/landings_data/KMills_landings by area 1964-2021 - FINFISH ONLY_MAY 2022.xlsx"), 
  sheet = 1, 
  skip = 0) %>% 
  clean_names()

# Organized by stat zone
by_szone <- read_xlsx(
  path = here("data/landings_data/KMills_landings by area 1964-2021 - FINFISH ONLY_MAY 2022.xlsx"), 
  sheet = 2, 
  skip = 0) %>% 
  clean_names()

# List of species
spec_list <- read_xlsx(
  path = here("data/landings_data/KMills_landings by area 1964-2021 - FINFISH ONLY_MAY 2022.xlsx"), 
  sheet = 3, 
  skip = 0) %>% 
  clean_names()


####  Data Exploration  ####


# 1. Aggregate Statzones by the region definitions of the project

# Define the regions
# Stratum Area Key for which stratum correspond to larger regions we use
strata_key <- list(
  "Georges Bank"          = as.character(13:23),
  "Gulf of Maine"         = as.character(24:40),
  "Southern New England"  = stringr::str_pad(
    as.character(1:12),
    width = 2, pad = "0", side = "left"),
  "Mid-Atlantic Bight"    = as.character(61:76))


# Add the labels to the data
trawldat <- dplyr::mutate(
  by_szone,
  survey_area =  dplyr::case_when(
    stat_area %in% strata_key$`Georges Bank`         ~ "GB",
    stat_area %in% strata_key$`Gulf of Maine`        ~ "GoM",
    stat_area %in% strata_key$`Southern New England` ~ "SNE",
    stat_area %in% strata_key$`Mid-Atlantic Bight`   ~ "MAB",
    TRUE                                             ~ "stratum not in key"))


# Join to landings
by_szone %>% 
  group_by(stat_area)
  
  
#### Aggregate by Area: ####
# Aggregate the landings by survey regions:

