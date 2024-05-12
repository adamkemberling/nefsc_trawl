#### SizeSpectra MLE Methodology Comparisons
# Date: 5/16/2024
# Status: fine, but strange behavior around weighting biomass based on counts/cpue
# Looking into/reviewing whether the best "negLL.fn" was being used to estimate exponents of size spectra
# Needed to be confident for methods reporting

# I had initially been using the MLE bins methodology: negLL.fn = negLL.PLB.bins.species
# But I believe that the growth rate differences are negligible between species at a 1cm measurement scale
# Meaning that we probably could simplify the code down






####  Packages   ####
library(sizeSpectra)
library(gmRi)
library(targets)
library(tidyverse)




####  Starting Points  ####

tar_load(catch_log2_labelled)


# Goal, swap group_isd_calc
# for a function that applies estimation using non-binned mle estimation


# Current Setup
# group_isd_calc() uses the negLL.PLB.bins.species
# This approach takes a dataframe containing:
# species name/code
# wmin & wmax for the size(s) of the fish
# number of counts for each species x bin





# Did alot of prep to have information for each method in the data
# can simplify a lot here
# Follow this documentation from Andrew Edwards
# https://htmlpreview.github.io/?https://raw.githubusercontent.com/andrew-edwards/sizeSpectra/master/doc/MEPS_IBTS_1.html
glimpse(catch_log2_labelled)


# Right now we have all lengths x species x tows as records
# we can aggregate to whatever group level we're using and just use the stratified abundances or CPUES as weights


# This is where group column choices come into play:
# And also the weighting to use
# Function Options
.group_cols <- c("Year", "survey_area")
abundance_vals <- "strat_mean"
min_mass_g <- 4 # From blanchard et al 2005

# Decide which abundance/cpue value to weight with using "Number"
.abund_col = switch(
  EXPR = abundance_vals,
  "observed"   = sym("numlen_adj"),
  "strat_mean" = sym("strat_mean_abund_s"),
  "stratified" = sym("strat_total_abund_s"))

# Now aggregate by body-weight within the groups we are measuring spectra for:
.agg_cols <- c(.group_cols, "comname", "ind_weight_g")
ss_input <- catch_log2_labelled %>% 
  group_by(!!!syms(.agg_cols)) %>% 
  summarise(
    Number = sum(!!.abund_col),
    LWa = unique(a), 
    LWb = unique(b),
    .groups = "drop")  

# Drop columns we don't need, rename to match edwards code
ss_input <- ss_input %>% 
  unite(col = "group_var", {{.group_cols}}, sep = "_", remove = FALSE, na.rm = FALSE) %>% 
  select(
    group_var, 
    Number, 
    LWa, 
    LWb, 
    bodyMass = ind_weight_g)

# Now filter the range of sizes we want to include
ss_input <- ss_input %>% 
  filter(bodyMass > min_mass_g)



# Set the PLB Constants using All Group Data:
# Set using all-group min/max
x <- rep(ss_input$bodyMass, ss_input$Number)
sum.log.x <- sum( log(x) )
isd_xmin = min(x)
isd_xmax = max(x)

# Set starting point for non-linear fit solver
# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
PL.bMLE = 1/( log(min(ss_input$bodyMass)) - sum.log.x/length(x)) - 1


# ----- Everything below should be done on the individual groups


# Proceed with one group of the list for speed:
ss_input_i <- split(ss_input, f = "group_var")[[1]]



####  calcLike Function ####


# Now Estimate b using negLL.PLB
n_i <- sum(ss_input_i$Number) # How many total

#### Weighting Trouble - Memory Exhausting  ####
# Weight the vector by repeating bodymass x abundance
# For strat abundance its too much
x_i <- 
# If we use stratified mean cpue then we can get relative differences but its strange
x_i <- rep(ss_input_i$bodyMass, (ss_input_i$Number/ min(ss_input_i$Number)))
sum_log_x_i <- sum( log(x_i) )

# How is this working without giving it any bodymass data
# Its not changing if I supply an x vector
calcLike(
  negLL.fn = negLL.PLB, 
  p = PL.bMLE, 
  x = x_i,
  xmin = isd_xmin, 
  xmax = isd_xmax,
  n = n_i,
  sumlogx = sum_log_x_i)

# -----------------------------

#### NLM Directly  ####

# How its done in eightmethods:
# NOTES: Performs unexpectedly when given xmin/xmax values that aren't specific to group


# Use the group subset bodymass vector as "x"
x <- ss_input_i$bodyMass # Raw
x <- x_i # Weighted by strat abundance etc

# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# as a starting point for nlm for MLE of b for PLB model.
PL.bMLE <- 1/( log(min(x)) - sum(log(x)) / length(x)) - 1

PLB.minLL =  nlm(
  f = negLL.PLB, 
  p = PL.bMLE, 
  x = x, 
  n = length(x),
  xmin = min(x), 
  xmax = max(x), 
  sumlogx = sum(log(x)))

# Check code, should be 1 or 2 to pass
PLB.minLL$code

# Pull slope estimate
PLB.bMLE = PLB.minLL$estimate
PLB.bMLE

# Confidence intervals
# 95% confidence intervals for MLE method.

PLB.minNegLL = PLB.minLL$minimum

# Values of b to test to obtain confidence interval. For the real movement data
#  sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
#  symmetric interval here.

bvec = seq(PLB.bMLE - 0.75, PLB.bMLE + 0.75, 0.001)  # If make 0.0001 then do
# get an interval for raw 1980 data

PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
for(i in 1:length(bvec)){
  PLB.LLvals[i] = negLL.PLB(
    bvec[i], 
    x = x, 
    n = length(x), 
    xmin = xmin,
    xmax = xmax, 
    sumlogx = sum.log.x)
}
critVal = PLB.minNegLL  + qchisq(0.95,1)/2
# 1 degree of freedom, Hilborn and Mangel (1997) p162.
bIn95 = bvec[ PLB.LLvals < critVal ]

# b values in 95% confidence interval
PLB.MLE.bConf = c(min(bIn95), max(bIn95))
if(PLB.MLE.bConf[1] == min(bvec) | PLB.MLE.bConf[2] == max(bvec)){ 
  plot(bvec, PLB.LLvals)
  abline(h = critVal, col="red")
  stop("Need to make bvec larger - see R window")   # Could automate
}


