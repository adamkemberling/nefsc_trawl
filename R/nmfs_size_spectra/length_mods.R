# Load raw catch
library(tidyverse)
library(gmRi)
library(sizeSpectra)

####  Load Data  ####

# Just read it in, no* species filtering
catch_basic <- gmri_survdat_prep(
  survdat = NULL, 
  survdat_source = "most recent", 
  box_location = "cloudstorage")

# Subset to species from wigley
catch_wigley <- add_lw_info(catch_basic, cutoff = F, box_location = "cloudstorage")


# Subset to finfish
ecodata::species_groupings

# --- quick approach, remove things measured in mm

# So if the species is a decapod or shellfish then they estimate length to the nearest mm. 
# Otherwise length is to the nearest cm. So this is important to distinguish for lobster vs. finfish predators.
crusts <- catch_basic %>%
  mutate(length_char = as.character(length_cm)) %>%
  filter(grepl("[.]", length_char)) %>%
  distinct(comname) %>% 
  pull(comname)
crusts

catch_basic %>% 
  filter(comname %not in% crusts) %>% 
  distinct(comname) %>% 
  pull()

# thorough -- join taxa info, subset finfish families



# --- Drop duplicate length*species*


#----- Run length spectra
finfish_catch <- filter(catch_basic, comname %not in% crusts)

# test one group
test_g <- finfish_catch %>% 
  filter(est_year == 2016, season == "Fall", survey_area == "GoM")

 

####  Test Group - mle method  ####


# summarise numbers at length, by length
ss_input_summs <- test_g %>% 
  group_by(est_year, season, survey_area, length_cm) %>% 
  summarise(
    total_number_at_length = sum(numlen_adj),
    .groups = "drop") 




# Get the subgroup inputs:

# Vector of body weights, repeat length_cm x abundance
x_i <- rep(ss_input_summs$length_cm, ceiling(ss_input_summs$total_number_at_length))

# Get many total individuals
n_i <- sum(ceiling(ss_input_summs$total_number_at_length) )

# Sum( log(length_cm) )
sum_log_xi <- sum( log(x_i)) 

# Analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
PL.bMLE  <- 1/( log(min(x_i)) - sum_log_xi / length(x_i)) - 1

# mmin/max
max_i <- max(ss_input_summs$length_cm)


# Do the exponent estimation for group i
isd_results <- group_est <- calcLike(
  negLL.fn = negLL.PLB,  # MLE function, takes a vector of lengths
  p        = PL.bMLE,    # Starting value
  x        = x_i,        # Vector
  xmin     = 1,          # min
  xmax     = max_i,      # max
  n        = n_i,        # total number
  sumlogx  = sum_log_xi) # 



####  Test Group - mleBins method  ####
 #-------- bin method

# Next Step
# Use the non-integer count version - mlebins

bin_summary <- test_g %>% 
  group_by(est_year, season, survey_area, length_cm, comname) %>% 
  summarise(
    total_number_at_length = sum(numlen_adj),
    .groups = "drop")  %>% 
  mutate(wmin = length_cm,
         wmax = length_cm + 1,
         Number = total_number_at_length)


# just rename it for Edwards
dataBinForLike <- bin_summary


# 3. Set the power-law limits:
# n, xmin, xmax

# Total individuals
n    <- sum(dataBinForLike$Number)

# left bound, grams

# Set manually
isd_xmin <- 1
isd_xmax <- NULL

# Set based on data
if(is.null(isd_xmax)){ isd_xmax <- max(dataBinForLike$wmax, na.rm = T)}
if(is.null(isd_xmin)){isd_xmin <- min(dataBinForLike$wmin, na.rm = T)}


# Calculate b
bin_results <- calcLike(
  negLL.fn          = negLL.PLB.bins.species,
  p                 = -1.9,
  suppress.warnings = TRUE,
  dataBinForLike    = dataBinForLike,
  n                 = n,
  xmin              = isd_xmin,
  xmax              = isd_xmax)


# Compare the two outputs
isd_results$MLE
bin_results$MLE







# ----- Make it a function  ------------


# Code to do mle for one group of data:



#' @title {MLE Size Spectra Estimation}
#'
#' @param ss_input Dataframe containing a column of abundance and a column of sizes
#' @param grouping_vars string identifiers of columns to group_by prior to analysis
#' @param abundance_vals string indicating column of abundances
#' @param size_vals string indicating column with individual length/weight data
#' @param isd_xmin lower limit for size distribution fitting
#' @param isd_xmax upper limit for size distribution fitting
#' @param global_min T/F Enforce a minimum size across all groups
#' @param global_max T/F Enforce a maximum size across all groups
#'
#' @return
#' @export
#'
#' @examples
group_binspecies_mle <-  function(
    ss_input, 
    grouping_vars, 
    abundance_vals = "numlen_adj",
    size_vals = "ind_weight_g",
    isd_xmin = NULL,
    isd_xmax = NULL,
    global_min = TRUE,
    global_max = TRUE,
    bin_width = 1){
  
  # 1. toggle which abundance/weight columns to use with switch
  .abund_col  <- sym(abundance_vals)
  .size_col   <- sym(size_vals)
  .group_cols <- grouping_vars
  .agg_cols   <- c(.group_cols, "comname", size_vals)
  
  
  # 2. Select only the columns we need:
  # Aggregate numbers by size 
  # within the groups we are measuring:
  ss_input_summs <- ss_input %>% 
    group_by(!!!syms(.agg_cols)) %>% 
    summarise(
      Number = sum(!!.abund_col),
      .groups = "drop")
  
  
  # Create a group column that we can map() through
  # Drop columns we don't need, rename to match edwards code 
  # edwards calls this df databinforlike in eightMethods()
  mle_input <- ss_input_summs  %>% 
    unite(col = "group_var", {{.group_cols}}, sep = "-", remove = FALSE, na.rm = FALSE)  %>% 
    select(
      group_var, 
      Number, 
      wmin = !!.size_col)  %>% 
    # Set up the min/max for the bins
    mutate(wmax = wmin + bin_width)

  #------ Optional - Set Global Constants
  
  # Power-law limits:
  # set left bound & right bounds, grams
  if(global_min == TRUE){
    if(is.null(isd_xmin)){ isd_xmin <- min(mle_input$wmin, na.rm = T)}
  }
  if(global_max == TRUE){
    if(is.null(isd_xmax)){ isd_xmax <- max(mle_input$wmax, na.rm = T)}
  }
  
  
  #------ Loop through Groups 
  group_results_df <- mle_input %>% 
    split(.$group_var) %>% 
    map_df(
      function(ss_input_i){
        
        # 3. Tune subgroup the power-law limits:
        # set left bound & right bounds, grams/cm
        min_i <- min(ss_input_i$wmin, na.rm = T)
        max_i <- max(ss_input_i$wmax, na.rm = T)
        if(global_min == FALSE){
          if(is.null(isd_xmin)){ isd_xmin <- min_i}
        }
        if(global_max == FALSE){
          if(is.null(isd_xmax)){ isd_xmax <- max_i}
        }
        
        # Filter the range of sizes we want to include
        ss_input_i <- ss_input_i %>% 
          filter(wmin >= isd_xmin,
                 wmax <= isd_xmax)
        
        # Total individuals in subgroup
        n_i <- sum(ceiling(ss_input_i$Number) )
        
        # Do the exponent estimation for group i
        group_est <- calcLike(
          negLL.fn          = negLL.PLB.bins.species,
          p                 = -1.9,
          suppress.warnings = TRUE,
          dataBinForLike    = ss_input_i,
          n                 = n_i,
          xmin              = isd_xmin,
          xmax              = isd_xmax)
        
        # Put it into a dataframe to rejoin neatly
        mle_group_results <- data.frame(
          xmin_fit = isd_xmin,
          xmax_fit = isd_xmax,
          xmin_actual = min_i,
          xmax_actual = max_i,
          n = n_i,
          b = group_est$MLE,
          confMin = group_est$conf[1],
          confMax = group_est$conf[2])
        
        
        # Process C and standard error
        mle_group_results <- mle_group_results %>% 
          mutate(
            stdErr = (abs(confMin - b) + abs(confMax - b)) / (2 * 1.96),
            C = (b != -1 ) * (b + 1) / ( xmax_fit^(b + 1) - xmin_fit^(b + 1) ) + (b == -1) * 1 / ( log(xmax_fit) - log(xmin_fit)))
        return(mle_group_results)}, 
      
      # Name the column with group ID
      .id = "group_var") %>% 
    
    # Decompose the group ID back into original columns
    separate(
      group_var, 
      sep = "-", 
      into = grouping_vars, 
      remove = F)
  
  
  # Spit it out
  return(group_results_df)
  
}




# Test it

# All fiunfish
length_binspectra <- group_binspecies_mle(
  ss_input = finfish_catch, 
  grouping_vars = c("est_year", "season", "survey_area"), 
  abundance_vals = "numlen_adj", 
  size_vals = "length_cm", 
  isd_xmin = 1, 
  isd_xmax = NULL, 
  global_min = TRUE, 
  global_max = FALSE, 
  bin_width = 1)

# Does the single group match the values run together
test_g %>% distinct(est_year, season, survey_area)
bin_results$MLE
length_binspectra %>% filter(est_year == 2016, season == "Fall", survey_area == "GoM")
-0.4435251


# Plot it
library(tidyquant)
length_binspectra %>% 
  filter(season %in% c("Spring", "Fall")) %>% 
  mutate(yr_num = as.numeric(as.character(est_year))) %>% 
  ggplot(aes(yr_num, b, color = season)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_ma(aes(linetype = "5-Year Moving Average"),n = 5, ma_fun = SMA) +
  geom_smooth(method = "lm", linewidth = 1, se = F, 
              aes(linetype = "Regression Fit")) +
  facet_wrap(~survey_area) +
  scale_color_gmri() +
  labs(title = "Length Spectra Slope - MLE Bins Method",
       subtitle = "Enforced xmin = 1, xmax = max(length_cm + 1)",
       y = "b",
       x = "Year",
       color = "Season")




# # Just the Wigley species
# length_binspectra_wigley <- group_binspecies_mle(
#     ss_input = catch_wigley, 
#     grouping_vars = c("est_year", "season", "survey_area"), 
#     abundance_vals = "numlen_adj", 
#     size_vals = "length_cm", 
#     isd_xmin = 1, 
#     isd_xmax = NULL, 
#     global_min = TRUE, 
#     global_max = FALSE, 
#     bin_width = 1)




####  Covariate Data  ####

# levels for faceting areas
area_levels <- c("GoM", "GB", "SNE", "MAB")
area_levels_long <- c("Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight")

# table to join for swapping shorthand for long-hand names
area_df <- data.frame(
  area = c("Scotian Shelf", "Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight", "All"),
  survey_area = c("SS", "GoM", "GB", "SNE", "MAB", "Northeast Shelf"),
  area_titles = c("Scotian Shelf", "Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight", "Northeast Shelf"))


# Load raw landings values
landings_raw <- read_csv(here::here("data/unscaled_spectra_predictor_dataset_wide.csv")) %>% 
  select(year, survey_area, landings_raw = landings)


# Read in du pontavice bottom temperatures
bot_temps <- read_csv(str_c(cs_path("res", "Du_Pontavice_Combined_BT/RegionalTimeseries"), "trawl_region_bottom_temps.csv"))







####  LMER With Covariates  ####



# Put it all together and clarify the column names
seasonal_df <- length_binspectra %>% 
  left_join(area_df) %>% 
  select(-survey_area) %>% 
  rename(year = est_year,
         survey_area = area)  %>% 
  mutate(year = as.numeric(as.character(year))) %>% 
  left_join(landings_raw) %>% 
  left_join(bot_temps) %>% 
  mutate(
    survey_area = factor(
      survey_area, 
      levels = area_levels_long),
    region = str_replace_all(survey_area, "-| ", "_"),
    region = factor(
      region, 
      levels = c("Gulf_of_Maine", "Georges_Bank", 
                 "Mid_Atlantic_Bight", "Southern_New_England")))





library(lme4)
library(lmerTest)
library(emmeans)
library(merTools)


glimpse(seasonal_df)


#lm first

# Drop NA
model_df <- drop_na(seasonal_df, landings_raw, bot_temp, b) %>% 
  filter(season %in% c("Spring", "Fall"))

# remove 

# Proportion/residuals of landings not explained by bot temp
landings_temp_lm <- lm(landings_raw ~ bot_temp, data = model_df)
model_df$land_resid <- resid(landings_temp_lm)




# lmer
elmer <- lmerTest::lmer(
  formula = b ~ scale(bot_temp) * region * season + scale(land_resid) * region + (1 | year),
  data = model_df)

summary(elmer)


# ggpredictions
library(ggeffects)
simple_preds <- as.data.frame(
  ggpredict(elmer, ~ bot_temp + season + region) )

# Plot over observed data
simple_preds %>% 
  mutate(
    season = factor(group, levels = c("Spring", "Fall")),
    region = factor(facet, levels = str_replace_all(area_levels_long, " |-", "_"))) %>% 
  ggplot() +
  geom_ribbon(aes(x, ymin = conf.low, ymax = conf.high, group = season), alpha = 0.1) +
  geom_line(
    aes(x, predicted, color = season, linetype = "Predicted"), 
    linewidth = 1) +
  facet_wrap(~region) +
  geom_point(
    data = model_df, 
    aes(bot_temp, b, color = season),
    alpha = 0.6) +
  geom_line(
    data = model_df,
    aes(bot_temp, b, color = season, linetype = "Observed"),
    alpha = 0.3, linewidth = 0.5) +
  scale_linetype_manual(values = c(3, 1)) +
  labs(y = "b")



# not pairwise
slope_phoc <- emtrends(
  object = elmer, 
  specs =  ~ region * season,
  var = "bot_temp",
  adjust = "sidak")



slope_phoc %>% 
  as_tibble() %>% 
  mutate(flag_alpha = ifelse(lower.CL > 0 | upper.CL <0, 1, 0.3)) %>% 
  ggplot(aes(region, bot_temp.trend, ymin = lower.CL, ymax = upper.CL)) +
  geom_hline(yintercept = 0, linetype = 1, color = "gray30", linewidth = 1) +
  geom_pointrange(aes(alpha = I(flag_alpha)), color = gmri_cols("blue"), size = 1) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  facet_wrap(~fct_rev(season), ncol = 2) +
  labs(
    title = "Region * Seasonal Slopes",
    x = NULL,
    y = "yr_num Trend")






# Landings residuals
# not pairwise
slope_phoc <- emtrends(
  object = elmer, 
  specs =  ~ region,
  var = "land_resid",
  adjust = "sidak")

slope_phoc %>% 
  as_tibble() %>% 
  mutate(flag_alpha = ifelse(lower.CL > 0 | upper.CL <0, 1, 0.3)) %>% 
  ggplot(aes(region, land_resid.trend, ymin = lower.CL, ymax = upper.CL)) +
  geom_hline(yintercept = 0, linetype = 1, color = "gray30", linewidth = 1) +
  geom_pointrange(aes(alpha = I(flag_alpha)), color = gmri_cols("blue"), size = 1) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  labs(
    title = "Region Slopes",
    x = NULL,
    y = "yr_num Trend")
