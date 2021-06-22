##### 
# 4/28/2021
# Paralell steps to running size spectra following mizer methods


####  Packages  ####

library(here)
library(targets)
library(gmRi)
library(janitor)
library(mizer)
library(sizeSpectra)
library(patchwork)
library(tidyverse)

####  Support Functions  ####
source(here("R/support/sizeSpectra_support.R"))



####  Data  ####

# Load the data from targets
tar_load(nefsc_databin)

# name change
dataBin <- nefsc_databin


####  Binning  ####

"
 It was used in the R pack-age mizer 
(Scott, Blanchard & Andersen 2014), which simu-lates the potential consequences 
of various ﬁshing patternsusing an approach based on the McKendrick–von
Foerster equation and calculates resulting size spectra. The user speciﬁesthe 
number of bins, and the lower bounds of the lowest andhighest bins.


In the R package mizer (Scott et al., 2014) the user specifies bins by giving the lower
bounds of the smallest and largest bins (min w and max w) and also the number of bins
(no w). The lower bounds of the bins are then calculated as:

w < −10^(seq(from = log10(min w), to = log10(max w), length.out = no w))
"



# So using that definition, here is how mizer would be specified
min_w <- min(dataBin$wmin) #smallest bin
max_w <- max(dataBin$wmax) #largest bin

# Assuming minimum and maximum weights determine start and end bins:
bin1     <- log10(min_w)
bin_last <- log10(max_w)

# need to solve an equation for constant b
# the bin width that satisfies the equal width requirements
# issue is that it depends on k... and I don't know how to get around that
binwidth_fun <- function(b, k){
  b^(k-2) * (2*b - 1) - (max_w/min_w)
}



# Thanks to Andrew edwards for this
# code in eightMethods()

# number of bins A. Edwards Uses
no_w <- 8

# How Beta is derived
mizer_beta <- nlm(LBmizbinsFun, 2, xmin = min_w, xmax = max_w, k = num_binz)

# setting bin number from weight
bin_weight  <- function(w){w < -10^(seq(from = log10(min_w), to = log10(max_w), length.out = no_w))}

# Setting Bins
# check mizerbins.R for more detailed methods
