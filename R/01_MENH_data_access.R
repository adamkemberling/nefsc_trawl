####  ME/NH Trawl Data Access  ####



####  Packages  ####
library(here)
library(gmRi)
library(tidyverse)


####  Data  ####

# Downloaded from Shiny Data Access Portal
# https://mainedmr.shinyapps.io/MaineDMR_Trawl_Survey_Portal/
menh_len <- read_csv(here("data/MENH/shiny_length_freq.csv"), guess_max = 1e6)
                              
# Downloaded via correspondence
bio         <- read_csv(here("data/MENH/full_me_dmr_biological.csv"))
expcatch    <- read_csv(here("data/MENH/full_me_dmr_expcatch.csv"))
len_freq    <- read_csv(here("data/MENH/full_me_dmr_lengthfreq.csv"))
lobster_len <- read_csv(here("data/MENH/full_me_dmr_lobster_LF.csv"))

# Length Weight Relationships via fishbase, compiled by Ben Resek, REU
lw_coef <- read_csv("data/MENH/listfishusingfishbase.csv" )