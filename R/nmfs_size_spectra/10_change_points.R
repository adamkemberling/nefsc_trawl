# Change point analyses
# Goal: test methodologies for change point or break point analyses as
# alternatives to chronological clustering

# Options: bcp, mcp

####  Packages  ####
library(targets)
library(here)
library(gmRi)
library(patchwork)
library(rioja)
library(vegan)
library(tidyverse)
library(bcp)
library(mcp)


# ggplot theme
theme_set(theme_minimal() + theme(legend.position = "bottom"))


####  Size Spectra Data  ####

# OISST for all the regions
withr::with_dir(rprojroot::find_root('_targets.R'), tar_load(regional_oisst))    

# SS and manual bins together
withr::with_dir(rprojroot::find_root('_targets.R'), tar_load(size_spectrum_indices))   
size_spectrum_indices <- size_spectrum_indices  


# Filter just regions and years
year_region_slopes <- size_spectrum_indices %>% 
  filter(`group ID` == "single years * region") %>% 
  drop_na() # drop the "all years"

# Set factor order
year_region_slopes <- year_region_slopes%>% 
  mutate(season = fct_rev(season),
         survey_area = factor(survey_area, 
                              levels = c("GoM", "GB", "SNE", "MAB")),
         yr = as.numeric(as.character(Year))) 






# Gulf of Maine Tester
gom_test <- filter(year_region_slopes, survey_area == "GoM")
gom_sst <- filter(regional_oisst, survey_area == "GoM")





####  bcp package  ####
# https://rdrr.io/cran/bcp/man/bcp.html#heading-8




# Testing bcp
bcp_gom <- bcp(y = gom_test$b, burnin = 10000, mcmc = 10000, return.mcmc = T)
plot(bcp_gom)

# Looks like there's one area about a third of the way in:
bcp_sum <- as.data.frame(summary(bcp_gom))

# Let's filter the data frame and identify the year:
bcp_sum$id <- gom_test$yr
(sel <- bcp_sum[which(bcp_gom$posterior.prob > 0.3), ])

# plot with breakpoint
ggplot(gom_test, aes(yr, b)) + 
  geom_line() + 
  geom_point() + 
  geom_vline(data = sel, aes(xintercept = id), linetype = 2, color = "gray50") +
  ggrepel::geom_label_repel(data = sel, aes(x = id, y = -1.05, label = round(Probability, 2)))+
  labs(x = "Year", y = "sizeSpectra Slope Exponent (b)", subtitle = "GoM Bayesian Change Points")


#### BCP Meta Functions  ####


# run bcp and flag changepoints using threshold
flag_bcp <- function(response_var,  predictor_var = NULL, id_var, burnin = 10000, mcmc = 10000, return.mcmc = T, threshold_prob = 0.3){
  
  bcp_x <- bcp(y = response_var, x = predictor_var, burnin = burnin, mcmc = mcmc, return.mcmc = return.mcmc)
  sink("/dev/null")
  bcp_sum <- as.data.frame(summary(bcp_x))
  sink()
  bcp_sum$id <- id_var
  flag_dat <- bcp_sum[which(bcp_x$posterior.prob > threshold_prob), ]
  flag_dat <- as.data.frame(flag_dat)
  return(flag_dat)
  
}

# plot bcp breaks
plot_bcp_breaks <- function(bcp_data, bcp_breaks, region_label){
  # plot with breakpoint
  ggplot(bcp_data, aes(yr, b)) + 
    geom_line() + 
    geom_point() + 
    geom_vline(data = bcp_breaks, aes(xintercept = id), linetype = 2, color = "gray50") +
    ggrepel::geom_label_repel(data = bcp_breaks, 
                              aes(x = id, y = -1.05, label = str_c(round(Probability, 2), "%"))) +
    labs(x = "Year", y = "sizeSpectra Slope Exponent (b)", 
         subtitle = str_c(region_label, " - Bayesian Change Points"))
}

#### Run Different Areas  ####

# GoM
gom_bcp <- flag_bcp(response_var = gom_test$b, id_var = gom_test$yr, threshold_prob = 0.3)
plot_bcp_breaks(gom_test, gom_bcp, "Gulf of Maine")

# GB
gb_test <- filter(year_region_slopes, survey_area == "GB")
gb_bcp <- flag_bcp(response_var = gb_test$b, id_var = gb_test$yr, threshold_prob = 0.3)
plot_bcp_breaks(gb_test, gb_bcp, "Georges Bank")

#SNE
sne_test <- filter(year_region_slopes, survey_area == "GB")
sne_bcp <- flag_bcp(response_var = sne_test$b, id_var = sne_test$yr, threshold_prob = 0.3)
plot_bcp_breaks(sne_test, sne_bcp, "Southern NE")

# MAB
mab_test <- filter(year_region_slopes, survey_area == "MAB")
mab_bcp <- flag_bcp(response_var = mab_test$b, id_var = mab_test$yr, threshold_prob = 0.3)
plot_bcp_breaks(mab_test, mab_bcp, "Mid-Atlantic Bight")







####  Temperature  ####

# Plot sst breaks
plot_sst_breaks <- function(region){
  # filter data
  reg_sst <- filter(regional_oisst, survey_area == region, yr %in% c(1982:2019))
  
  # run changepoints
  reg_sst_bcp <- flag_bcp(reg_sst$sst_anom, id_var = reg_sst$yr, threshold_prob = 0.3)
  
  # build plot
  p1 <- ggplot(reg_sst, aes(yr, sst_anom)) +
    geom_line() +
    geom_point() 
  
  if(nrow(reg_sst_bcp > 0)){
    p1 <- p1 +
      geom_vline(data = reg_sst_bcp, aes(xintercept = id), linetype = 2, color = "gray50") +
      ggrepel::geom_label_repel(data = reg_sst_bcp, 
                                aes(x = id, y = -1.05, label = str_c(round(Probability, 2), "%")))}
   
  
    p1 + labs(x = "Year", 
              y = expression("Temperature Anomaly"~degree*C), 
              subtitle = str_c(region, " SST", " - Bayesian Change Points"))
}

# Gulf of Maine
plot_sst_breaks("GoM")

# Georges Bank
plot_sst_breaks("GB")

# southern New England
plot_sst_breaks("SNE")

# MAB
plot_sst_breaks("MAB")






####  bcp with covariates  ####

# The bcp package gives the option of using external variables in the changepoint analysis
# under the assumption that a linear model is appropriate

# # Tweak for the inclusion of sst
# bcp_with_sst <- function(region, threshold){
#   
#   sst_dat <- filter(regional_oisst, survey_area == region, yr %in% c(1982:2019))
#   bcp_dat <- filter(year_region_slopes, survey_area == region, yr %in% c(1982:2019)) %>% 
#     left_join(sst_dat, by = c("yr"))
#   
#   lm_bcp <- flag_bcp(response_var = bcp_dat$b, predictor_var = bcp_dat$sst_anom,  id_var = bcp_dat$yr, threshold_prob = threshold)
#   return(lm_bcp)
#   
# }


# # run regions
# gom_sst_bcp <- bcp_with_sst("GoM", 0.3)
# gb_sst_bcp <- bcp_with_sst("GB", 0.3)
# mab_sst_bcp <- bcp_with_sst("MAB", 0.3)
# sne_sst_bcp <- bcp_with_sst("SNE", 0.3)



# Plot both sst and b
plot_bcp_with_sst <- function(region, threshold){
  
  # Filter sst data
  reg_sst <- filter(regional_oisst, survey_area == region, yr %in% c(1982:2019))
  
  # Grab sizeSpectrra results
  bcp_dat <- filter(year_region_slopes, survey_area == region, yr %in% c(1982:2019)) %>% 
    left_join(reg_sst, by = c("yr"))
  
  
  reg_sst_bcp <- flag_bcp(response_var = bcp_dat$b, predictor_var = bcp_dat$sst_anom,  id_var = bcp_dat$yr, threshold_prob = threshold)
  
  
  
  
  # Build temperature plot
  p1 <- ggplot(reg_sst, aes(yr, sst_anom)) +
    geom_line() +
    geom_point() 
  
  if(nrow(reg_sst_bcp > 0)){
    p1 <- p1 +
      geom_vline(data = reg_sst_bcp, aes(xintercept = id), linetype = 2, color = "gray50")}
  
    p1 <- p1 + labs(x = "Year", 
                    y = expression("Temperature Anomaly"~degree*C),
                    subtitle = "Linear Predictor: sst_anom")
  
    
    
    # Build sizespectra plot with breakpoint
    p2 <- ggplot(bcp_dat, aes(yr, b)) + 
      geom_line() + 
      geom_point() 
    
    if(nrow(reg_sst_bcp > 0)){
    p2 <- p2 + 
      geom_vline(data = reg_sst_bcp, aes(xintercept = id), linetype = 2, color = "gray50") +
      ggrepel::geom_label_repel(data = reg_sst_bcp, 
                                aes(x = id, y = -1.05, label = str_c(round(Probability, 2), "%")))}
    p2 <- p2 +
      labs(x = "Year", y = "sizeSpectra Slope Exponent (b)", 
           subtitle = str_c(region, " | b ~ sst_anom - Bayesian Change Points"))
    
    return(p2 / p1)
    
}




#plot sst and the 
plot_bcp_with_sst("GoM", 0.1)
plot_bcp_with_sst("GB", 0.1)
plot_bcp_with_sst("MAB", 0.1)
plot_bcp_with_sst("SNE", 0.1)


####  mcp  ####

# # demo example has varying intercepts and slopes
# mcp_example("demo")
# 
# # Define the model
# model = list(
#   b ~ 1,    # Plateau (Int_1)
#   ~ 0 + yr, # joined slope (time_2) at cp_1
#   ~ 1 + yr  # disjoined slope (int_3, time_3) at cp_2
# )
# 
# # Run the mcp model
# gom_mcp <- mcp(model = model, data = gom_test)
# 
# # Summary
# summary(gom_mcp)
# 
# # Plot it
# plot_pars(gom_mcp, regex_pars = "cp_")
# 
# 
# 
# 
# # Null Model
# # Define the model
# model_null <- list(
#   b ~ 1 + yr,  # intercept (int_1) and slope (time_1)
#   ~ 1 + yr            # disjoined slope (int_2, time_1)
# )
# 
# # Fit it
# fit_null = mcp(model_null, gom_test)
# 
# 
# # Model Comparison
# gom_mcp$loo <- loo::loo(gom_mcp)
# fit_null$loo <- loo::loo(fit_null) 
