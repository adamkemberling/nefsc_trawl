####  NEFSC Trawl Data Access  ####



####  Packages  ####
library(here)
library(gmRi)
library(tidyverse)


####  Data  ####

# From Janet Nye
mills_path <- shared.path(group = "Mills Lab", folder = "")
load(str_c(mills_path, "Data/Survdat_Nye_allseason.RData"))

# Gonna save it locally for convenience/laziness
#write_csv(survdat, here("data/NEFSC/2019survdat_nye.csv"))

# 2020 edition
load(here("data/NEFSC/Survdat_Nye_Aug 2020.RData"))

# Gonna save it locally for convenience/laziness
# write_csv(survdat, here("data/NEFSC/2020survdat_nye.csv"))
