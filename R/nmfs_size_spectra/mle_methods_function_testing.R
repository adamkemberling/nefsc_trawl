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

# Set theme
theme_set(
  theme_gmri() + 
    theme(
      text = element_text(family = "Avenir", size = 11),
      plot.margin = margin(5,5,5,5), 
      legend.position = "bottom", 
      legend.direction = "horizontal", 
      legend.title = element_text(size = 12, face = "bold")) )


####  Starting Points  ####

# Load the data we are feeding in to the analysis
tar_load(catch_log2_labelled)





####  Objectives  ####
# Goal, swap function in group_isd_calc
# for a function that applies estimation using non-binned mle estimation method


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




####__________________####
####  Write Function  ####


# # These were above
# .group_cols <- c("Year", "survey_area")
# abundance_vals <- "observed"
# min_mass_g <- 4 # From blanchard et al 2005





#' @title {MLE Size Spectra Estimation}
#'
#' @param ss_input Dataframe containing a column of abundance and a column of sizes
#' @param grouping_vars string identifiers of columns to group_by prior to analysis
#' @param abundance_vals string indicating column of abundances
#' @param weight_vals string indicating column with individual weight data
#' @param isd_xmin lower limit for size distribution fitting
#' @param isd_xmax upper limit for size distribution fitting
#'
#' @return
#' @export
#'
#' @examples
group_mle_estimates <-  function(
    ss_input, 
    grouping_vars, 
    abundance_vals = "numlen_adj",
    weight_vals = "ind_weight_g",
    isd_xmin = NULL,
    isd_xmax = NULL){
  
  
  # 1. toggle which abundance/weight columns to use with switch
  .abund_col <- sym(abundance_vals)
  .weight_col <- sym(weight_vals)
  .group_cols <- grouping_vars
  .agg_cols <- c(.group_cols, "comname", weight_vals)
  
  # 2. Select only the columns we need:
  # Now aggregate by body-weight within the groups we are measuring spectra for:
  ss_input <- ss_input %>% 
    group_by(!!!syms(.agg_cols)) %>% 
    summarise(
      Number = sum(!!.abund_col),
      LWa = unique(a), 
      LWb = unique(b),
      .groups = "drop")  
  
  
  # Create a group column that we can loop through
  # Drop columns we don't need, rename to match edwards code
  # edwards calls this df databinforlike
  ss_input <- ss_input %>% 
    unite(col = "group_var", {{.group_cols}}, sep = "-", remove = FALSE, na.rm = FALSE) %>% 
    select(
      group_var, 
      Number, 
      bodyMass = !!.weight_col)
  
  
  
  #------ Set all-group constants
  # set these one time
  
  # 3. Set the power-law limits:
  # set left bound & right bounds, grams
  if(is.null(isd_xmin)){ isd_xmin <- min(ss_input$bodyMass, na.rm = T)}
  if(is.null(isd_xmax)){ isd_xmax <- max(ss_input$bodyMass, na.rm = T)}
  
  
  # Now filter the range of sizes we want to include
  ss_input <- ss_input %>% 
    filter(bodyMass > isd_xmin,
           bodyMass < isd_xmax)
  
  
  # Get many total individuals
  n <- sum(ceiling(ss_input$Number) )

  # Vector of body weights, repeat bodymass x abundance
  x <- rep(ss_input$bodyMass, ceiling(ss_input$Number))
  
  # Sum( log(bodymass) )
  sum_log_x <- sum( log(x) )
  
  # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
  PL.bMLE  <- 1/( log(min(x)) - sum_log_x / length(x)) - 1
  
  
  #------ Loop through Groups  ------
  group_results_df <- ss_input %>% 
    split(.$group_var) %>% 
    map_df(
      function(ss_input_i){
        
        # Get the subgroup inputs:
        # Vector of body weights, repeat bodymass x abundance
        x_i <- rep(ss_input_i$bodyMass, ceiling(ss_input_i$Number))
        
        # Get many total individuals
        n_i <- sum(ceiling(ss_input_i$Number) )
        
        # Sum( log(bodymass) )
        sum_log_xi <- sum( log(x_i)) 
    
        # Do the estimate for the group
        group_est <- calcLike(
          negLL.fn = negLL.PLB,  # MLE function, takes a vector of weights
          p        = PL.bMLE, 
          x        = x_i,
          xmin     = isd_xmin, 
          xmax     = isd_xmax,
          n        = n_i,
          sumlogx  = sum_log_xi)
        
        # Put it into a df
        mle_group_results <- data.frame(
          xmin_fit = isd_xmin,
          xmax_fit = isd_xmax,
          xmin_actual = min(ss_input_i$bodyMass),
          xmax_actual = max(ss_input_i$bodyMass),
          n = n_i,
          b = group_est$MLE,
          confMin = group_est$conf[1],
          confMax = group_est$conf[2])
        
        
        # Process C and standard error
        mle_group_results <- mle_group_results %>% 
          mutate(
            stdErr = (abs(confMin - b) + abs(confMax - b)) / (2 * 1.96),
            C = (b != -1 ) * (b + 1) / ( xmax_fit^(b + 1) - xmin_fit^(b + 1) ) + (b == -1) * 1 / ( log(xmax_fit) - log(xmin_fit)))
        return(mle_group_results)
    
    }, .id = "group_var") %>% 
    separate(group_var, sep = "-", into = grouping_vars)
  
  # Spit it out
  return(group_results_df)
  
}



# Does it work?c- yes, cool
regional_mle_spectra <- catch_log2_labelled %>% 
  drop_na(numlen_adj, ind_weight_g) %>% 
  group_mle_estimates(
    ss_input = .,
    grouping_vars = c("Year", "survey_area"), 
    abundance_vals = "numlen_adj", 
    weight_vals = "ind_weight_g", 
    isd_xmin = 4, 
    isd_xmax = 10000)




# Plot
regional_mle_spectra %>% 
  mutate(yr_num = as.numeric(as.character(Year)),
         survey_area = factor(survey_area, levels = c("GoM", "GB", "SNE", "MAB" ))) %>% 
ggplot(aes(yr_num, b, color = survey_area)) +
  geom_line(linewidth = 1) + 
  scale_color_gmri() +
  facet_wrap(~survey_area)


# What does the seasonal view look like
seasonal_mle_spectra <- catch_log2_labelled %>% 
  drop_na(numlen_adj, ind_weight_g) %>% 
  group_mle_estimates(
    ss_input = .,
    grouping_vars = c("Year", "survey_area", "season"), 
    abundance_vals = "numlen_adj", 
    weight_vals = "ind_weight_g", 
    isd_xmin = 4, 
    isd_xmax = 10000)


# Plot it - regret looking into it\
seasonal_mle_spectra %>% 
  mutate(
    yr_num = as.numeric(as.character(Year)),
    survey_area = factor(survey_area, levels = c("GoM", "GB", "SNE", "MAB" )),
    season = factor(season, levels = c("Spring", "Fall"))) %>% 
  ggplot(aes(yr_num, b, color = season)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm") +
  scale_color_gmri() +
  scale_fill_gmri() +
  facet_wrap(~survey_area, ncol = 2)



####_____________________####
# Everything below was done when prototyping, then shifted below the function for testing

####  Input/Function Controls  ####

# Right now we have all lengths x species x tows as records
# we can aggregate to whatever group level we're using and just use the stratified abundances or CPUES as weights


# This is where group column choices come into play:
# And also the weighting to use
# Function Options
.group_cols <- c("Year", "survey_area")
abundance_vals <- "observed"
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



####  Set Global Controls  ####

# For consistency across groups/times
# set a starting point for estimation using:
# total abundance
# global min/max bounds



# Get the sizes as one long vector by expanding out abundances
x <- rep(ss_input$bodyMass, ceiling(ss_input$Number))

# Sum of log(biomass)
sum.log.x <- sum( log(x) )

# Set bounds using all-group min/max
isd_xmin = min(x)
isd_xmax = max(x)


# Set starting point for non-linear fit solver
# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
PL.bMLE <- 1/( log(min(x)) - sum.log.x / length(x)) - 1



#### Perform Sub-Group Estimations  ####

# ----- Everything below should be done on the individual groups


# Proceed with one group of the list for speed:
ss_input_i <- split(ss_input, f = "group_var")[[1]]


####  Using the calcLike Function ####

# Estimate b using negLL.PLB


# How many total
n_i <- sum(ceiling(ss_input_i$Number) )

# Vector of body weights
# Expand the weights vector by repeating bodymass x abundance
x_i <- rep(ss_input_i$bodyMass, ceiling(ss_input_i$Number))

# Sum( log(bodymass) )
sum_log_x_i <- sum( log(x_i) )


# How is this working without giving it any bodymass data
# Its not changing if I supply an x vector
calcLike(
  negLL.fn = negLL.PLB,  # MLE function, takes a vector of weights
  p = PL.bMLE, 
  x = x_i,
  xmin = isd_xmin, 
  xmax = isd_xmax,
  n = n_i,
  sumlogx = sum_log_x_i)




#### Calling NLM Directly  ####

# This is how its done in eightmethods:'
# Gets the same result
# NOTES: Performs unexpectedly when given xmin/xmax values that aren't specific to group


# Use the group subset bodymass vector as "x"
x_i <- x_i

# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# as a starting point for nlm for MLE of b for PLB model.
PL.bMLE <- 1/( log(min(x_i)) - sum(log(x_i)) / length(x)) - 1

PLB.minLL =  nlm(
  f = negLL.PLB, 
  p = PL.bMLE, 
  x = x_i, 
  n = length(x_i),
  xmin = min(x_i), 
  xmax = max(x_i), 
  sumlogx = sum(log(x_i)))

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
    xmin = isd_xmin,
    xmax = isd_xmax, 
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






####  Using negLL.PLB.bins.species  ####
