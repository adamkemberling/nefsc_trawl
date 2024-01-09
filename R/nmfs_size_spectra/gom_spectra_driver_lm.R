####  Size Spectra Linear Model Selection  
# Taking this out of manuscript file to dig into outputs

####  Packages  ####
library(EnvCpt)
library(targets)
library(here)
library(gmRi)
library(patchwork)
library(scales)
library(emmeans)
library(tidyverse)


# Package Conflicts
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")

# Support functions
source(here("R/support/sizeSpectra_support.R"))

# Resource Path
res_path <- gmRi::cs_path("res")

# Set a seed for reproducing stochastic elements
set.seed(123)




####  Abundance Data  ####

# 1. Biological data used as input
tar_load(catch_log2_labelled)  

# rename and format
# Add the area titles
catch_size_bins <- catch_log2_labelled #%>% fill_func_groups()




#### 2.   Spectra Slope Estimates  ####
tar_load(size_spectrum_indices)


# Grab SS Groups we care about
region_indices <- size_spectrum_indices  %>% 
  filter(`group ID` == "single years * region") %>% 
  mutate(
    yr = as.numeric(as.character(Year)),
    survey_area = factor(survey_area, levels = c("GoM", "GB", "SNE", "MAB")))




# # Export for Bart
# write_csv(region_indices, here::here("data/multispecies_spectra_results.R"))


####  Spectra Driver Datasets  ####


# This is the unscaled predictor dataset



# these are the environmental drivers that are hypothesized to matter
# scaled by mean+sd within region
drivers_scaled <- read_csv(here::here("data/env_drivers_scaled.R"))




####  Lagging Predictors  ####

driver_ccf_prep <- drivers_scaled %>% 
  select(-All_sst_anom) %>% 
  #rownames_to_column(var = "year") %>% 
  pivot_longer(names_to = "spectra_param", 
               values_to = "spectra_values", 
               cols = ends_with("slope") | ends_with("int")) %>% 
  pivot_longer(names_to = "driver_var", 
               values_to = "driver_values", 
               cols = -matches("spectra|year")) %>% 
  mutate(
    # C. These are the areas associated with the Spectra Features
    spectra_area = case_when(
      str_detect(spectra_param, "All")          ~ "All",
      str_detect(spectra_param, "Georges")      ~ "Georges Bank",
      str_detect(spectra_param, "Gulf")         ~ "Gulf of Maine",
      str_detect(spectra_param, "Southern")     ~ "Southern New England",
      str_detect(spectra_param, "Mid-Atlantic") ~ "Mid-Atlantic Bight"),
    spectra_area = factor(spectra_area, 
                          levels = c("All", "Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight")),
    
    # D. These are the areas of the drivers
    driver_area = case_when(
      str_detect(driver_var, "All")          ~ "All",
      str_detect(driver_var, "Georges")      ~ "Georges Bank",
      str_detect(driver_var, "Gulf")         ~ "Gulf of Maine",
      str_detect(driver_var, "Southern")     ~ "Southern New England",
      str_detect(driver_var, "Mid-Atlantic") ~ "Mid-Atlantic Bight"),
    driver_area = factor(driver_area, 
                         levels = c("All", "Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight"))) 



# Filter it out so there is only cases where the driver area matches
driver_ccf_prep <- driver_ccf_prep %>% 
  filter(driver_area == spectra_area | driver_area == "All") %>% 
  arrange(year, driver_var)



# Function to grab the correlation data and lag data
get_ccf_vector <- function(x,y){
  
  # Run the ccf
  ccf_data <- ccf(x,y,plot= F , lag.max = 15)
  
  # Get the signif:
  # https://www.squaregoldfish.co.uk/programming/r_acf_significance.md/
  # Not 100% sure n is the same for ccf as it is for acf, but...
  significance_level <- qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(x)))
  
  data.frame(
    "acf" = ccf_data$acf,
    "lag" = ccf_data$lag,
    "sigpos" = significance_level,
    "signeg" = significance_level*-1
  )
}


# Walk through each xvariable, and see how it correlates with each yvar
# split on xvar
# then split on yvars
# make sure order is good
ccf_relationships <- driver_ccf_prep %>% 
  drop_na() %>% 
  split(.$spectra_param) %>% 
  map_dfr(function(x_param){
    x_param %>% 
      split(.$driver_var) %>% 
      map_dfr(function(driver_y_data){
        ccf_df <- get_ccf_vector(
          x = driver_y_data$spectra_values,
          y = driver_y_data$driver_values
        )
      }, .id = "driver_var")
  }, .id = "spectra_param")




# Build back out the labels for plotting
ccf_plot_data <- ccf_relationships %>% 
  mutate(
    # Flag what the driver type was
    driver_type = case_when(
      str_detect(driver_var, "landings") ~ "Commercial Landings",
      str_detect(driver_var, "sst") ~ "SST",
      str_detect(driver_var, "index") ~ "GSI",
      str_detect(driver_var, "small") ~ "ZP-S",
      str_detect(driver_var, "large") ~ "ZP-L",
      TRUE ~ "Missed Something"),
    param_feature = case_when(
      # Flag what the Size Distribution Parameter was
      str_detect(spectra_param, "int") ~ "Spectra Intercept",
      str_detect(spectra_param, "isd") ~ "ISD Exponent",
      str_detect(spectra_param, "slope") ~ "Spectra Slope",
      TRUE ~ "Missed Something"),
    spectra_region = case_when(
      # Flag what region the driver was coming from
      str_detect(spectra_param, "All") ~ "All",
      str_detect(spectra_param, "Georges") ~ "Georges Bank",
      str_detect(spectra_param, "Gulf") ~ "Gulf of Maine",
      str_detect(spectra_param, "Southern") ~ "Southern New England",
      str_detect(spectra_param, "Mid-Atlantic") ~ "Mid-Atlantic Bight"),
    # Flag when it crosses threshold
    sig_flag = ifelse(acf < signeg | acf > sigpos, T, F),
    # Set Factor Levels
    spectra_region = factor(spectra_region, levels = c("All", "Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight")))


# Limit to one Response:
ccf_plot_data <- filter(ccf_plot_data, param_feature == "ISD Exponent") %>% 
  filter(spectra_region != "All")


# Plot them
ccf_plot_data %>% 
  filter(spectra_region == "Gulf of Maine") %>% 
  ggplot(aes(
    x = lag, 
    y = acf, 
    group = driver_var, 
    color = driver_type, 
    fill = driver_type)) +
  geom_col(alpha = 0.65) +
  geom_col(
    data = filter(ccf_plot_data, sig_flag, spectra_region == "Gulf of Maine"),
    aes(lag, acf, fill = driver_type,  group = driver_var), color = "black") +
  geom_text(
    data = filter(ccf_plot_data, sig_flag, spectra_region == "Gulf of Maine"),
    aes(lag, y = 0, label = lag, vjust = ifelse(acf < 0, -0.5,1.5)), color = "black", size = 2.5) +
  geom_line(aes(x = lag, y = sigpos), linetype = 2, color = "gray25") +
  geom_line(aes(x = lag, y = signeg), linetype = 2, color = "gray25") +
  geom_hline(yintercept = 0, color = "gray25", linewidth = 1) +
  geom_vline(xintercept = 0, color = "black", linewidth = 1) +
  scale_color_gmri() +
  scale_fill_gmri() +
  facet_grid(spectra_region~driver_type) +
  scale_x_continuous(limits = c(-10, 10), breaks = seq(from = -15, to = 12, by = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 8)) +
  labs(x = "Driver Lag", y = "Correlation (CCF)", 
       fill = "Driver Type:", color = "Driver Type:",
       title = "Auto-correlation Function:\nSpectrum Slope ~ Predictors",
       caption = "All predictors scaled over year range of 1982-2019. Lags exceeding significance threshold outlined and numbered.")





####  Prepare Regression Dataframe  ####


# Use scaled drivers for the regressions
# Need to drop zooplankton because they are non-unique since we split the MAB EPU and repeated it
regression_df <- drivers_scaled %>% 
  select(!ends_with(c("zp_small","zp_large")))

# Build back out the metadata:
# Identify what region is associated with both: spectrum values & driver values
regression_df <- regression_df %>% 
  pivot_longer(names_to = "spectra_param", 
               values_to = "spectra_values", 
               cols = ends_with("slope") | ends_with("int")) %>% 
  pivot_longer(names_to = "driver_var", 
               values_to = "driver_values", 
               cols = -matches("spectra|year")) %>% 
  mutate(
    year = as.numeric(year),
    # A. Flag what the driver type was
    driver_type = case_when(
      str_detect(driver_var, "landings") ~ "Commercial Landings",
      str_detect(driver_var, "sst") ~ "SST",
      str_detect(driver_var, "index") ~ "GSI",
      TRUE ~ "Missed Something"),
    # B. Flag what the spectrum feature was
    param_feature = case_when(
      # Flag what the Size Distribution Parameter was
      str_detect(spectra_param, "int") ~ "Spectra Intercept",
      str_detect(spectra_param, "isd") ~ "ISD Exponent",
      str_detect(spectra_param, "slope") ~ "Spectra Slope",
      TRUE ~ "Missed Something"),
    # C. These are the areas associated with the Spectra Features
    spectra_area = case_when(
      str_detect(spectra_param, "All") ~ "All",
      str_detect(spectra_param, "Georges") ~ "Georges Bank",
      str_detect(spectra_param, "Gulf") ~ "Gulf of Maine",
      str_detect(spectra_param, "Southern") ~ "Southern New England",
      str_detect(spectra_param, "Mid-Atlantic") ~ "Mid-Atlantic Bight"),
    spectra_area = factor(spectra_area, levels = c("All", "Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight")),
    
    # D. These are the areas of the drivers
    driver_area = case_when(
      str_detect(driver_var, "All") ~ "All",
      str_detect(driver_var, "Georges") ~ "Georges Bank",
      str_detect(driver_var, "Gulf") ~ "Gulf of Maine",
      str_detect(driver_var, "Southern") ~ "Southern New England",
      str_detect(driver_var, "Mid-Atlantic") ~ "Mid-Atlantic Bight"),
    driver_area = factor(driver_area, levels = c("All", "Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight")),
    
    # E. Add Decade in
    decade = floor_decade(year),
    #
    #regime = ifelse(year < 2010, "regime_1", "regime_2" )
  ) 


# Filter it so predictors are only within matching regions, 
# or just overarching (GSI)
regression_df <- regression_df %>% 
  filter(driver_area == spectra_area | driver_type == "GSI") %>% 
  arrange(year, driver_var)




# Last wrangling:
# columns for: year, slope, isd, spectra region, drver_region, GSI, sst_anom, landings
# GSI is so annoying here, just rejoin it without the area column

# Join the climate modes together
gsi <- ecodata::gsi  %>% 
  mutate(
    Time = str_pad(as.character(Time), width = 7, side = "right", pad = "0"),
    year = as.numeric(str_sub(Time, 1, 4)),
    month = str_sub(Time, -2, -1),
    Time = as.Date(str_c(year, month, "01", sep = "-")))
clim_drivers <- gsi %>% filter(EPU == "All")


# Put climate drivers on annual schedule
clim_idx <- clim_drivers %>% 
  group_by(year, area = EPU, var = Var) %>% 
  summarise(value = mean(Value, na.rm = T),
            .groups = "drop")

# Join stupid GSI back in...
regression_df <- regression_df %>% 
  select(-c(spectra_param, driver_var)) %>% 
  pivot_wider(names_from = "driver_type", values_from = "driver_values") %>% 
  pivot_wider(names_from = "param_feature", values_from = "spectra_values") %>% 
  select(-GSI) %>% 
  filter(driver_area != "All") %>% 
  left_join( select(clim_idx, year, GSI = value), by = c("year") ) %>% 
  mutate(spectra_area = fct_drop(spectra_area))




# Load original data
zp_sli <- ecodata::zoo_sli_anom %>% filter(EPU != "SS")

# Reformat to match the other indices
zp_index <- zp_sli %>% 
  #EPUS combine SNE and MAB, we can repeat the values here
  bind_rows(filter(zp_sli, EPU == "MAB") %>% mutate(EPU = "SNE")) %>% 
  mutate(
    survey_area = case_when(
      EPU == "GB" ~ "Georges Bank",
      EPU == "GOM" ~ "Gulf of Maine",
      EPU == "MAB" ~ "Mid-Atlantic Bight",
      EPU == "SNE" ~ "Southern New England"),
    Units = str_c("zp_", Var)) %>% 
  select(year = Time, survey_area, Value, Units) %>% 
  pivot_wider(names_from = Units, values_from = Value)



# Add zooplankton here so it can be checked in the regressions
# Zooplankton values are duplicated in MAB and SNE b/c its the same EPU
regression_df <- left_join(
  regression_df, 
  zp_index %>% rename(spectra_area = survey_area), 
  by = join_by(year, spectra_area))







####  Breakpoint Structure  ####

# Does it make sense to do this?


# Pull Regions to test independently:
regression_df <- regression_df %>% 
  rename(
    landings = `Commercial Landings`,
    b = `ISD Exponent`)

# Gulf of Maine
gom_df <- filter(regression_df, spectra_area == "Gulf of Maine") %>% 
  drop_na()



# Check each for different changepoint structures
gom_cpt <- envcpt(gom_df$b, verbose = F)



# Pulling out details of best supported model
cpt_res <- list(
  "gom" =  list("cpt_res" = gom_cpt)#,
  # "gb"  =  list("cpt_res" = gb_cpt),
  # "sne" =  list("cpt_res" = sne_cpt),
  # "mab" =  list("cpt_res" = mab_cpt)
  ) %>% map(function(x){
    best <- names(which.min(AIC(x$cpt_res)))
    best_summ <- x$cpt_res[[best]]
    x$best <- best
    x$best_summ <- best_summ
    return(x)
  })


# # Plot GOM, the only one with changepoints
plot(gom_cpt, main = "Gulf of Maine, Break Point Analysis")



# Put the rest of the results into a delta AIC table:
bind_rows(map(cpt_res, ~tibble("Most Supported Model" = .x[["best"]])), .id = "Region") %>% knitr::kable()






####  Gulf of Maine  ####


gom_df <- mutate(gom_df, 
                 regime = case_when(
                   year < 1999 ~ "yrs_1982-1998",
                   year < 2007 ~ "yrs_1999_2006",
                   year >= 2007 ~ "yrs_2007_2019"))


# Add the autoregressive predictors neatly/manually
gom_df_lag <- gom_df %>% 
  mutate(
    blag1   = dplyr::lag(b, 1),
    blag2   = dplyr::lag(b, 2),
    zpslag1 = dplyr::lag(zp_small, 1),
    zpslag2 = dplyr::lag(zp_small, 2)) %>% 
  drop_na(blag2) %>% 
  mutate(regime = factor(regime, levels = c("yrs_1982-1998", "yrs_1999_2006", "yrs_2007_2019")))


# kitchen sink
gom_ar_models = lm(
  b ~ 
    regime * SST + 
    regime * GSI + 
    regime * landings + 
    regime * zp_small + 
    regime * zp_large + 
    regime * blag1 + 
    regime * blag2, 
  data = gom_df_lag)



# Model selection on models with lags at 1+2 years
library(MuMIn)
library(broom)

# change na. action for dredge
options(na.action = "na.fail")
results_gom <- dredge(gom_ar_models)


# Best autoregressive model
best_gom <- get.models(results_gom, subset = delta == 0)[[1]]



####  "Best" Model Exploration ####

# What do the partial dependence plots look like:
# they look like trash because the interactions are NA


# Create plot-data data frame with mean values of the other predictors
# And the different factor levels too
# langings ranges from the observed minimum to maximum


# Set up the daata to plot the range of predicted values for landings across regimes
plotdata <- expand.grid(
  landings = seq(min(gom_df_lag$landings), max(gom_df_lag$landings), by = .01),
  regime = unique(gom_df_lag$regime),
  blag2 = mean(gom_df_lag$blag2, na.rm = T),
  zp_large = mean(gom_df_lag$zp_large, na.rm = T),
  zp_small = mean(gom_df_lag$zp_small, na.rm = T),
  SST = mean(gom_df_lag$SST, na.rm = T)
)



# Restrict the ranges to the ranges seen within each regime?
# Makes it easier to see the actual ranges we observed
plotdata <- plotdata %>% 
  split(.$regime) %>% 
  imap_dfr(function(x,y){
    observed_range <- filter(gom_df_lag, regime == y) %>% pull(landings)
    filter(x, landings >= min(observed_range, na.rm = T), landings <= max(observed_range, na.rm = T))
  })


# Add Best model predictions
bmod_se <- predict(best_gom, plotdata, se.fit = T)[["se.fit"]]
plotdata <- plotdata %>% 
  mutate(bmod_preds = predict(best_gom, plotdata),
         landings_se = bmod_se,
         landings_hi = bmod_preds + 1.96*landings_se,
         landings_lo = bmod_preds - 1.96*landings_se)

# Plotting Differential Impact of landings
ggplot() +
  geom_ribbon(
    data = plotdata,
    aes(ymin = landings_lo, ymax = landings_hi, x = landings, fill = regime), alpha = 0.5)+
  geom_line(
    data = plotdata,
    aes(landings, bmod_preds, color = regime), linewidth = 1) +
  facet_wrap(~regime, scales = "free") +
  geom_point(data = gom_df_lag, aes(landings, b, color = regime), size = 2) +
  scale_color_gmri() +
  scale_fill_gmri() +
  labs(x = "Landings (z)", y = "Multi-Species Spectrum Slope", caption = "Points = Observed Values")








# Do it again for large zooplankton


# Set up the daata to plot the range of predicted values for landings across regimes
plotdata <- expand.grid(
  landings = mean(gom_df_lag$landings, na.rm = T),
  regime = unique(gom_df_lag$regime, na.rm = T),
  blag2 = mean(gom_df_lag$blag2, na.rm = T),
  zp_large = seq(min(gom_df_lag$zp_large), max(gom_df_lag$zp_large), by = .01),
  zp_small = mean(gom_df_lag$zp_small, na.rm = T),
  SST = mean(gom_df_lag$SST, na.rm = T)
)



# Restrict the ranges to the ranges seen within each regime?
# Makes it easier to see the actual ranges we observed
plotdata <- plotdata %>% 
  split(.$regime) %>% 
  imap_dfr(function(x,y){
    observed_range <- filter(gom_df_lag, regime == y) %>% pull(zp_large)
    filter(x, zp_large >= min(observed_range, na.rm = T), zp_large <= max(observed_range, na.rm = T))
  })


# Add Best model predictions
bmod_se <- predict(best_gom, plotdata, se.fit = T)[["se.fit"]]
plotdata <- plotdata %>% 
  mutate(bmod_preds = predict(best_gom, plotdata),
         zp_large_se = bmod_se,
         zp_large_hi = bmod_preds + 1.96*zp_large_se,
         zp_large_lo = bmod_preds - 1.96*zp_large_se)


# Plotting Differential Impact of zp_large
ggplot() +
  geom_ribbon(
    data = plotdata,
    aes(ymin = zp_large_lo, ymax = zp_large_hi, x = zp_large, fill = regime), alpha = 0.5)+
  geom_line(
    data = plotdata,
    aes(zp_large, bmod_preds, color = regime), linewidth = 1) +
  facet_wrap(~regime, scales = "free") +
  geom_point(data = gom_df_lag, aes(zp_large, b, color = regime), size = 2) +
  scale_color_gmri() +
  scale_fill_gmri() +
  labs(x = "Large Zooplankton Index (z)", y = "Multi-Species Spectrum Slope", caption = "Points = Observed Values",
       title = "Really Terrible Fit for Large Zooplankton and Regime Interaction",
       subtitle = "*At Mean Values for Other Predictors")



#### EMMEANS Interactions Check  ####
# Use emmeans to speed life up:

# Partial Dependence plots, way less work...
emmip(best_gom, regime ~ landings, 
      at = list(
        landings = gom_df_lag$landings))

# Zooplankton Large
emmip(best_gom, regime ~ zp_large, 
      at = list(
        zp_large = gom_df_lag$zp_large
      ))




####  Top Performing Model Weights/PRedictions  ####

# Top performers
# Pull details for models with delta AIC < 2
top_mods_gom <- get.models(results_gom, subset = delta < 2)

# Peformance, Gets the AIC/AICc order 
top_perf <- map_dfr(top_mods_gom, function(x){
  broom::glance(x) %>% 
    mutate(AICc = AICc(x))}, 
  .id = "id") %>% 
  arrange(AICc)

# Calculate delta AICc
top_perf$Delta <- top_perf$AICc - min(top_perf$AICc)

# Give the ranking based on that order
top_perf <- top_perf %>% 
  arrange(Delta) %>% 
  mutate(`Model ID` = row_number(),
         `Model ID` = str_c("Model ", `Model ID`))



# Retrieve Model Coefficients from "top" models
# join in the model number based on AIC order
top_coef <- map_dfr(top_mods_gom, tidy, .id = "id") %>% 
  left_join(select(top_perf, id, `Model ID`)) %>% 
  select(-id)


# Reshape Performance Specs to bind_row later
top_perf <- select(top_perf, `Model ID`, r2 = r.squared, Delta) %>% 
  pivot_longer(cols = -1, names_to = "term", values_to = "estimate")




###  Model Weights  ####

# Sum of weights - predictor importance
sw(top_mods_gom)

# Model averaging
gom.ests <- model.avg(top_mods_gom, revised.var = TRUE)
gom.ests


# Model preds for each
model.preds <- sapply(setNames(top_mods_gom, c("mod 1", "mod 2", "mod 3")), predict, newdata = gom_df_lag)

# Weighted average prediction
out.put <- model.sel(top_mods_gom$`10459`, top_mods_gom$`2073`, top_mods_gom$`10491`)
mod.ave.preds <- model.preds %*% Weights(out.put)

gom_df_lag %>% cbind(mod.ave.preds) %>% 
ggplot() +
  geom_line(aes(year, b, color = "Actual"), linewidth = 1) +
  geom_line(aes(year, mod.ave.preds, color = "Model-Averaged Predicition"), linewidth = 1) +
  labs(title = "Model Averaged Prediction")
  



#### Model Averaged Partial Dependence Plots  ####

# More interesting application: everything but one predictor set to mean value

# langings ranges from the observed minimum to maximum
landings <- seq(min(gom_df_lag$landings), max(gom_df_lag$landings), by = .01)

# Create plot-data data frame with mean values of the other predictors
# And the different factor levels too
# Set up the daata to plot the range of predicted values for landings across regimes
plotdata <- expand.grid(
  landings = seq(min(gom_df_lag$landings), max(gom_df_lag$landings), by = .01),
  regime = unique(gom_df_lag$regime),
  blag2 = mean(gom_df_lag$blag2, na.rm = T),
  zp_large = mean(gom_df_lag$zp_large, na.rm = T),
  zp_small = mean(gom_df_lag$zp_small, na.rm = T),
  SST = mean(gom_df_lag$SST, na.rm = T)
)



# Restrict the ranges to the ranges seen within each regime?
# Makes it easier to see the actual ranges we observed
plotdata <- plotdata %>% 
  split(.$regime) %>% 
  imap_dfr(function(x,y){
    observed_range <- filter(gom_df_lag, regime == y) %>% pull(landings)
    filter(x, landings >= min(observed_range, na.rm = T), landings <= max(observed_range, na.rm = T))
  })


# Predict response for the plot data with each model
model.preds = sapply(top_mods_gom, predict, newdata = plotdata)

# Weight the prediction from each model by its AIC weight
# and sum (matrix multiplication)
mod.ave4plot <- model.preds %*% Weights(out.put)

# plot the model averaged predicted densities vs elevation

#plot(mod.ave4plot ~ plotdata$landings, type = 'l', xlab="Elevation (m)", ylab="Model averaged predicted density")

plotdata %>% 
  mutate(avg_pred = mod.ave4plot[,1]) %>% 
  ggplot() +
  geom_line(aes(landings, avg_pred, color = regime, linetype = "Model Averaged Prediction")) +
  geom_point(data = gom_df_lag, aes(landings, b, color = regime), size = 2) +
  facet_wrap(~regime, nrow = 1, scales = "free") +
  scale_color_gmri() +
  scale_fill_gmri() +
  labs(
    y = "Multi-Species Spectrum Slope", x = "Landings Index (z)", color = "Regime")




# What about with confidence intervals?
# documentation: https://cran.r-project.org/web/packages/MuMIn/MuMIn.pdf
avg_preds <- predict(model.avg(top_mods_gom), plotdata, se.fit = T)
plotdata %>% 
  mutate(avg_pred = avg_preds$fit,
         pred_se = avg_preds$se.fit,
         pred_hi = avg_pred + 1.96*pred_se,
         pred_lo = avg_pred - 1.96*pred_se) %>% 
  ggplot() +
  geom_ribbon(aes(landings, ymin = pred_lo,  ymax = pred_hi, fill = regime), alpha = 0.3) +
  geom_line(aes(landings, avg_pred, color = regime, linetype = "Model Averaged Prediction")) +
  geom_point(data = gom_df_lag, aes(landings, b, color = regime), size = 2) +
  scale_color_gmri() +
  scale_fill_gmri() +
  facet_wrap(~regime, scales = "free") +
  labs(
    y = "Multi-Species Spectrum Slope", x = "Landings Index (z)", color = "Regime", fill = "Regime",
    title = "AIC Top Model Ensemble Prediction: Model-Averaging Prediction")
