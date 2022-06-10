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


# Join to landings
by_szone %>% 
  group_by(stat_area)
  
  
  
# Aggregate the landings