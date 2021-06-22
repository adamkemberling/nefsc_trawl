####  Size Structure EDA  ####

# Goal of this script is to dig into 
# results from mean length/weight targets


####  Packages  ####
library(targets)
library(here)
library(sf)
library(gmRi)
library(patchwork)
library(rioja)
library(vegan)
library(tidyverse)

# Support functions
source(here("R/support/sizeSpectra_support.R"))
# source(here("R/support/temp_support.R"))

####  Data  ####
tar_load("nefsc_bio_lw")
tar_load("mean_individual_sizes")
tar_load("mean_sizes_ss_groups")


# rename them to shorter names - cleanupp for plotting
nefsc_bio_lw <- nefsc_bio_lw %>% 
  filter(is.na(indwt) == FALSE) %>% 
  mutate(yr = as.numeric(as.character(Year)),
         season = fct_rev(season))


mean_sizes_taxons <- mean_individual_sizes  %>% 
  mutate(yr = as.numeric(as.character(Year)),
         season = fct_rev(season))
mean_size_ss_groups  <- mean_sizes_ss_groups  %>% 
  mutate(yr = as.numeric(as.character(Year)),
         season = fct_rev(season)) 

####  EDA  ####

# Data peak
glimpse(mean_sizes_ss_groups)
glimpse(mean_sizes_taxons)
glimpse(nefsc_bio_lw)

# Targeted Inspection Columns
target_cols <- c("Year", "comname", "spec_class", "length", "biomass", "indwt", "numlen")



#####  Lengths  ####

# mean sizes were grouped on:
# Year, survey_area, season, spec_class, fishery
gb_com_data <- mean_sizes_taxons %>% 
  filter(survey_area == "GB",
         fishery == "com")

# Data that goes with it
gb_bio <- tibble(nefsc_bio_lw) %>% 
  filter(fishery == "com",
         survey_area == "GB")


# Georges Bank Lengths
gb_com_data %>% 
  ggplot(aes(yr, mean_len_cm, color = spec_class)) +
  geom_line() +
  geom_point() +
  facet_grid(spec_class~season) +
  labs(y = "Avg Length (cm)", 
       x = "",
       subtitle = "GB - Commercial Species Only")

# odd lengths of demersals in 1991
gb_bio %>% 
  filter(Year == 1992,
         spec_class == "dem") %>% 
  select(all_of(target_cols)) %>% 
  arrange(desc(length))



#####  Weights  ####


gb_com_data %>% 
  ggplot(aes(yr, mean_wt_kg, color = spec_class)) +
  geom_line() +
  geom_point() +
  facet_grid(spec_class~season,
             scales = "free") +
  labs(y = "Avg Individual Bodymass (kg)", 
       x = "",
       subtitle = "GB - Commercial Species Only")







