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
size_spectrum_indices <- size_spectrum_indices  %>% 
  mutate(season = fct_rev(season),
         survey_area = factor(survey_area, 
                              levels = c("GoM", "GB", "SNE", "MAB")),
         yr = as.numeric(as.character(Year)))


# Filter just regions and years
year_region_slopes <- size_spectrum_indices %>% 
  filter(`group ID` == "single years * region")



# Gulf of Maine Tester
gom_test <- filter(year_region_slopes, survey_area == "GoM")






####  bcp package  ####

# Testing bcp
bcp_gom <- bcp(y = gom_test$b, burnin = 10000, mcmc = 10000, return.mcmc = T)
plot(bcp_gom)

# Looks like theres one area about a third of the way in:
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
flag_bcp <- function(response_var, id_var, burnin = 10000, mcmc = 10000, return.mcmc = T, threshold_prob = 0.3){
  
  bcp_x <- bcp(y = response_var, burnin = burnin, mcmc = mcmc, return.mcmc = return.mcmc)
  bcp_sum <- as.data.frame(summary(bcp_x))
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
regional_oisst %>% distinct(survey_area)

# Gulf of Maine
gom_sst <- filter(regional_oisst, survey_area == "GoM")
gom_sst_bcp <- flag_bcp(gom_sst$sst_anom, gom_sst$yr, threshold_prob = 0.3)
ggplot(gom_sst, aes(yr, sst_anom)) +
  geom_line() +
  geom_point() +
  geom_vline(data = gom_sst_bcp, aes(xintercept = id), linetype = 2, color = "gray50") +
  ggrepel::geom_label_repel(data = gom_sst_bcp, 
                            aes(x = id, y = -1.05, label = str_c(round(Probability, 2), "%"))) +
  labs(x = "Year", y = "sizeSpectra Slope Exponent (b)", 
       subtitle = str_c("GOM SST", " - Bayesian Change Points"))

# Georges Bank
gb_sst <- filter(regional_oisst, survey_area == "GB")
gb_sst_bcp <- flag_bcp(gb_sst$sst_anom, gb_sst$yr, threshold_prob = 0.3)
ggplot(gom_sst, aes(yr, sst_anom)) +
  geom_line() +
  geom_point() +
  geom_vline(data = gb_sst_bcp, aes(xintercept = id), linetype = 2, color = "gray50") +
  ggrepel::geom_label_repel(data = gb_sst_bcp, 
                            aes(x = id, y = -1.05, label = str_c(round(Probability, 2), "%"))) +
  labs(x = "Year", y = "sizeSpectra Slope Exponent (b)", 
       subtitle = str_c("GB SST", " - Bayesian Change Points"))


####  mcp  ####

# demo example has varying intercepts and slopes
mcp_example("demo")

# Define the model
model = list(
  b ~ 1,    # Plateau (Int_1)
  ~ 0 + yr, # joined slope (time_2) at cp_1
  ~ 1 + yr  # disjoined slope (int_3, time_3) at cp_2
)

# Run the mcp model
gom_mcp <- mcp(model = model, data = gom_test)

# Summary
summary(gom_mcp)

# Plot it
plot_pars(gom_mcp, regex_pars = "cp_")




# Null Model
# Define the model
model_null <- list(
  b ~ 1 + yr,  # intercept (int_1) and slope (time_1)
  ~ 1 + yr            # disjoined slope (int_2, time_1)
)

# Fit it
fit_null = mcp(model_null, gom_test)


# Model Comparison
gom_mcp$loo <- loo::loo(gom_mcp)
fit_null$loo <- loo::loo(fit_null) 
