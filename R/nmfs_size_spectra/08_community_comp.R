#### Redundancy Analysis
"
Community composition: Changes in community composition will be tracked through time using
redundancy analysis (RDA), a multivariate ordination method that is robust to zero-inflated
data and that can be used to test how explanatory variables influences species assemblages
(Legendre and Legendre 2006). 

Before applying the distance-based RDA, the species matrix 
will be Hellinger transformed (Legendre and Gallagher 2001) to meet assumptions of normality.
The RDA creates principal components that are linear combinations of the species data, thereby
reducing assemblage data to only a few axes. These axes will be used to represent changes in
zooplankton and fish community composition over time. The RDA will be performed using the 
R vegan package (Oksanen et al. 2018) for zooplankton and fish separately.

Resources:
https://mb3is.megx.net/gustame/constrained-analyses/rda
"

####  Packages  ####
library(vegan)
library(gmRi)
library(targets)
library(tidyverse)


####  Data  ####

# Want to load in the abundances as a matrix here:
# want there to be all species not just the length-weight species
tar_load(nefsc_databin)



#  Pivot abundances
glimpse(nefsc_databin)

# So we have five levels of abundances to choose from...:
# 1. survey catch at length: numlen_adj
# 2. survey catch of all lengths: abundance
# 3. catch per tow at length: abund_tow_s
# 4. catch/tow at length, weighted by stratum area ratio: wt_abund_s
# 5. total expected catch across whole strata (at length): expanded_abund_s

# Hellinger transformation