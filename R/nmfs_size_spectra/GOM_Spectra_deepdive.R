##### 9/12/2023
##### Digging in to GOM Spectra Composition
##### Goals

##### 1. Try to find any signs of gaps in the size spectra
##### 2. Highlight where a major species or functional group is occupying a gap
##### 3. Look at residuals from MLE method predicted abundance
##### 4. Look at large fish and small fish indices as supportive evidence




####  Libraries  ####


####  Packages  ####
library(targets)
library(here)
library(gmRi)
library(patchwork)
library(gt)
library(tidyverse)
library(scales)
library(gtsummary)
library(dismo)
library(rcartocolor)

# Package Conflicts
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")

# Support functions
source(here("R/support/sizeSpectra_support.R"))

# Resource Path
res_path <- gmRi::cs_path("res")

# GGplot theme
theme_set(theme_gmri() + theme(
  panel.grid.major = element_line(linewidth = 0.5, linetype = 1, color = "gray"),
  panel.grid.major.x = element_line(linewidth = 0.5, linetype = 3, color = "gray"),
  panel.grid.minor = element_line(linewidth = 0.5, linetype = 3, color = "gray90"),
  axis.line = element_line(color = "black"),
  axis.line.y = element_line(color = "black"), 
  strip.background = element_rect(colour = gmri_cols("teal")),
  legend.position = "bottom"))


####  Abundance Data  ####

# 1. Biological data used as input
tar_load(catch_log2_labelled)  

# rename and format
# Add the area titles
catch_size_bins <- catch_log2_labelled %>% 
  fill_func_groups()



####  Spectra Estimates  ####

tar_load(size_spectrum_indices)


# Grab SS Groups we care about
region_indices <- size_spectrum_indices  %>% 
  filter(`group ID` == "single years * region") %>% 
  mutate(yr = as.numeric(as.character(Year)),
         survey_area = factor(survey_area, levels = c("GoM", "GB", "SNE", "MAB")),
         sig_fit = ifelse(log2_sig_strat < 0.05, "Significant", "Non-Significant"))



####  Plotting Spectra & Residuals  ####



####  Preparing the Fitted Curve  
# Prepare inputs directly:
yr_op <- 1978
reg <- "GoM"
mle_results <- region_indices %>% filter(yr == yr_op, survey_area == reg)
actual_vals <- catch_size_bins %>% filter(Year == yr_op, survey_area == reg) 


# Power law parameters and summary details for the group of data:
b.MLE           <- mle_results$b
total_abundance <- mle_results$n
b.confMin       <- mle_results$confMin
b.confMax       <- mle_results$confMax


# Create range of x values from the group to get power law predictions
# PLB = bounded power-law
# min and max weights for predictions
# xmin  <- mle_results$xmin
# xmax  <- mle_results$xmax
xmin  <- 1
xmax  <- 10^4

# Create x values (individual bodymass) to predict across
# break up the Xlim into pieces between min and max
x.PLB <- seq(from = xmin, 
             to   = xmax,
             by = 0.1)
             #length.out = 2000)   

# get the length of that vector
x.PLB.length <- length(x.PLB)  

# # remove last entry, add an entry .9999 of the way there, and cap it with the last entry wtf
# x.PLB <- c(
#   x.PLB[-x.PLB.length], 
#   0.99999 * x.PLB[x.PLB.length], 
#   x.PLB[x.PLB.length])


# Y values for plot limits/bounds/predictions from bounded power law pdf
y.PLB         <- (1 - sizeSpectra::pPLB(x = x.PLB, b = b.MLE, xmin = min(x.PLB), xmax = max(x.PLB))) * total_abundance
y.PLB.confMin <- (1 - sizeSpectra::pPLB(x = x.PLB, b = b.confMin, xmin = min(x.PLB), xmax = max(x.PLB))) * total_abundance
y.PLB.confMax <- (1 - sizeSpectra::pPLB(x = x.PLB, b = b.confMax, xmin = min(x.PLB), xmax = max(x.PLB))) * total_abundance


# Put it in a df to make it easier
PLB_df <- data.frame(
  x.PLB   = x.PLB,
  y.PLB   = y.PLB,
  confMin = y.PLB.confMin,
  confMax = y.PLB.confMax) %>% 
  mutate()



# 2. solve the power law function again in mutate for residuals
# need to do it twice for wmin_g and wmax_g

# Aggregate across overlapping size ranges
p_prepped <- actual_vals %>% 
  isd_plot_prep(stratified_abundance = T, min_weight_g = 1) %>% 
  left_join(actual_vals) %>% 
  select(comname, hare_group,  wmin_g, wmax_g, lowCount, highCount, countGTEwmin)


# Just do it in mutate
p_prepped <- p_prepped %>% 
  mutate(
    # Estimated Abundance
    wmin_estimate = (1 - sizeSpectra::pPLB(x = wmin_g, b = b.MLE, xmin = xmin, xmax = xmax)) * total_abundance,
    wmax_estimate = (1 - sizeSpectra::pPLB(x = wmax_g, b = b.MLE, xmin = xmin, xmax = xmax)) * total_abundance,
    # Residuals
    b_resid = countGTEwmin - ((wmin_estimate + wmax_estimate)/2)) 


#####  Plotting Raw  ####

# Plotting the pieces
ggplot() +
  geom_segment(
    data = p_prepped,
    aes(x = wmin_g, 
        xend = wmax_g, 
        y = countGTEwmin, 
        yend = countGTEwmin),
    color = gmri_cols("blue"),
    alpha = 0.5,
    linewidth = 0.8) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)), limits = c(10^6, 10^9)) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)), limits = c(xmin, xmax)) +
  labs(x = "Weight (g)", y = "Abundance", title = "Actual Vals")

#####  Plotting Functional Groups  ####


# Plotting the pieces
ggplot() +
  geom_line(data = PLB_df, aes(x.PLB, y.PLB), color = gmri_cols("orange"), linewidth = 1) +
  geom_segment(
    data = p_prepped %>% filter(lowCount >0),
    aes(x = wmin_g, 
        xend = wmax_g, 
        y = countGTEwmin, 
        yend = countGTEwmin),
    color = gmri_cols("blue"),
    alpha = 0.5,
    linewidth = 0.8) +
  facet_wrap(~hare_group) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)), limits = c(10^6, 10^9)) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)), limits = c(xmin, xmax)) +
  labs(x = "Weight (g)", y = "Abundance", title = "Actual Vals + Fit")

                


#####  Plotting Residuals  ####


# Plotting the pieces
ggplot() +
  geom_hline(yintercept = 0, color = "black", linewidth = 1) +
  geom_segment(
    data = p_prepped %>% filter(lowCount >0),
    aes(x = wmin_g, 
        xend = wmax_g, 
        y = b_resid, 
        yend = b_resid,
        color = hare_group),
    alpha = 0.8,
    linewidth = 0.8) +
  facet_wrap(~hare_group) +
  scale_color_gmri() +
  scale_y_continuous(labels = comma_format(scale = 1/1e6, suffix = "M")) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)), limits = c(xmin, xmax)) +
  labs(x = "Weight (g)", y = "Residual Abundance", title = "Size Spectra Estimate Residuals")








##### Specific Regimes  ####

####  Regime 1: 1982-1998  ####


# Prepare inputs directly:
yr_op <- c(1982:1998)
reg <- "GoM"


# Data for the spectra
actual_vals <- catch_size_bins %>% 
  filter(
    Year %in% yr_op, 
    survey_area == reg) 


# Get the slope over those years
mle_results <- actual_vals  %>% 
  group_isd_estimation(
    min_weight_g = 1, 
    max_weight_g = 10000,
    isd_xmin = 1,
    isd_xmax = 10000,
    abundance_vals = "stratified",
    .group_cols = c("survey_area")) 


# Power law parameters and summary details for the group of data:
b.MLE           <- mle_results$b
total_abundance <- mle_results$n
b.confMin       <- mle_results$confMin
b.confMax       <- mle_results$confMax


# Create range of x values from the group to get power law predictions
# min and max weights for power law
xmin  <- 1
xmax  <- 10^4

# Create x values (individual bodymass) to predict across
# break up the Xlim into pieces between min and max
x.PLB <- seq(
  from = xmin, 
  to   = xmax,
  by = 0.1)  

# get the length of that vector
x.PLB.length <- length(x.PLB)  

# Y values for plot limits/bounds/predictions from bounded power law pdf
y.PLB         <- (1 - sizeSpectra::pPLB(x = x.PLB, b = b.MLE, xmin = min(x.PLB), xmax = max(x.PLB))) * total_abundance
y.PLB.confMin <- (1 - sizeSpectra::pPLB(x = x.PLB, b = b.confMin, xmin = min(x.PLB), xmax = max(x.PLB))) * total_abundance
y.PLB.confMax <- (1 - sizeSpectra::pPLB(x = x.PLB, b = b.confMax, xmin = min(x.PLB), xmax = max(x.PLB))) * total_abundance


# Put it in a df to make it easier
PLB_df <- data.frame(
  x.PLB   = x.PLB,
  y.PLB   = y.PLB,
  confMin = y.PLB.confMin,
  confMax = y.PLB.confMax) %>% 
  mutate()



# 2. solve the power law function again in mutate for residuals
# need to do it twice for wmin_g and wmax_g

# Aggregate across overlapping size ranges
p_prepped <- actual_vals %>% 
  isd_plot_prep(stratified_abundance = T, min_weight_g = 1) %>% 
  left_join(actual_vals) %>% 
  select(comname, hare_group,  wmin_g, wmax_g, lowCount, highCount, countGTEwmin) %>% 
  # Get the residuals
  mutate(
    # Estimated Abundance
    wmin_estimate = (1 - sizeSpectra::pPLB(x = wmin_g, b = b.MLE, xmin = xmin, xmax = xmax)) * total_abundance,
    wmax_estimate = (1 - sizeSpectra::pPLB(x = wmax_g, b = b.MLE, xmin = xmin, xmax = xmax)) * total_abundance,
    # Residuals
    b_resid = countGTEwmin - ((wmin_estimate + wmax_estimate)/2)) 



#####  Plotting Raw  ####

# Plotting the pieces
r1_1 <- ggplot() +
  geom_segment(
    data = p_prepped,
    aes(x = wmin_g, 
        xend = wmax_g, 
        y = countGTEwmin, 
        yend = countGTEwmin),
    color = gmri_cols("blue"),
    alpha = 0.5,
    linewidth = 0.8) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)), limits = c(xmin, xmax)) +
  labs(x = "Weight (g)", y = "Abundance", title = "Regime 1: Actual Vals")

#####  Plotting Functional Groups  ####


# Plotting the pieces
r1_2 <- ggplot() +
  geom_line(data = PLB_df, aes(x.PLB, y.PLB), color = gmri_cols("orange"), linewidth = 1) +
  geom_segment(
    data = p_prepped %>% filter(lowCount >0),
    aes(x = wmin_g, 
        xend = wmax_g, 
        y = countGTEwmin, 
        yend = countGTEwmin, 
        color = hare_group),
    #color = gmri_cols("blue"),
    alpha = 0.5,
    linewidth = 0.8) +
  #facet_wrap(~hare_group) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)), limits = c(xmin, xmax)) +
  labs(x = "Weight (g)", y = "Abundance", title = "Regime 1: Actual Vals + Fit")




#####  Plotting Residuals  ####


# Plotting the pieces
r1_3 <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linewidth = 1) +
  geom_segment(
    data = p_prepped %>% filter(lowCount >0),
    aes(x = wmin_g, 
        xend = wmax_g, 
        y = b_resid, 
        yend = b_resid,
        color = hare_group),
    alpha = 0.8,
    linewidth = 0.8) +
  #facet_wrap(~hare_group) +
  scale_color_gmri() +
  scale_y_continuous(labels = comma_format(scale = 1/1e6, suffix = "M")) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)), limits = c(xmin, xmax)) +
  labs(x = "Weight (g)", y = "Residual Abundance", title = "Regime 1: Size Spectra Estimate Residuals")








####  Regime 2: 1999-2006  ####

# Prepare inputs directly:
yr_op <- c(1999:2006)
reg <- "GoM"


# Data for the spectra
actual_vals <- catch_size_bins %>% 
  filter(
    Year %in% yr_op, 
    survey_area == reg) 


# Get the slope over those years
mle_results <- actual_vals  %>% 
  group_isd_estimation(
    min_weight_g = 1, 
    max_weight_g = 10000,
    isd_xmin = 1,
    isd_xmax = 10000,
    abundance_vals = "stratified",
    .group_cols = c("survey_area")) 


# Power law parameters and summary details for the group of data:
b.MLE           <- mle_results$b
total_abundance <- mle_results$n
b.confMin       <- mle_results$confMin
b.confMax       <- mle_results$confMax


# Create range of x values from the group to get power law predictions
# min and max weights for power law
xmin  <- 1
xmax  <- 10^4

# Create x values (individual bodymass) to predict across
# break up the Xlim into pieces between min and max
x.PLB <- seq(
  from = xmin, 
  to   = xmax,
  by = 0.1)  

# get the length of that vector
x.PLB.length <- length(x.PLB)  

# Y values for plot limits/bounds/predictions from bounded power law pdf
y.PLB         <- (1 - sizeSpectra::pPLB(x = x.PLB, b = b.MLE, xmin = min(x.PLB), xmax = max(x.PLB))) * total_abundance
y.PLB.confMin <- (1 - sizeSpectra::pPLB(x = x.PLB, b = b.confMin, xmin = min(x.PLB), xmax = max(x.PLB))) * total_abundance
y.PLB.confMax <- (1 - sizeSpectra::pPLB(x = x.PLB, b = b.confMax, xmin = min(x.PLB), xmax = max(x.PLB))) * total_abundance


# Put it in a df to make it easier
PLB_df <- data.frame(
  x.PLB   = x.PLB,
  y.PLB   = y.PLB,
  confMin = y.PLB.confMin,
  confMax = y.PLB.confMax) %>% 
  mutate()



# 2. solve the power law function again in mutate for residuals
# need to do it twice for wmin_g and wmax_g

# Aggregate across overlapping size ranges
p_prepped <- actual_vals %>% 
  isd_plot_prep(stratified_abundance = T, min_weight_g = 1) %>% 
  left_join(actual_vals) %>% 
  select(comname, hare_group,  wmin_g, wmax_g, lowCount, highCount, countGTEwmin) %>% 
  # Get the residuals
  mutate(
    # Estimated Abundance
    wmin_estimate = (1 - sizeSpectra::pPLB(x = wmin_g, b = b.MLE, xmin = xmin, xmax = xmax)) * total_abundance,
    wmax_estimate = (1 - sizeSpectra::pPLB(x = wmax_g, b = b.MLE, xmin = xmin, xmax = xmax)) * total_abundance,
    # Residuals
    b_resid = countGTEwmin - ((wmin_estimate + wmax_estimate)/2)) 



#####  Plotting Raw  ####

# Plotting the pieces
r2_1 <- ggplot() +
  geom_segment(
    data = p_prepped,
    aes(x = wmin_g, 
        xend = wmax_g, 
        y = countGTEwmin, 
        yend = countGTEwmin),
    color = gmri_cols("blue"),
    alpha = 0.5,
    linewidth = 0.8) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)), limits = c(xmin, xmax)) +
  labs(x = "Weight (g)", y = "Abundance", title = "Regime 2: Actual Vals")

#####  Plotting Functional Groups  ####


# Plotting the pieces
r2_2 <- ggplot() +
  geom_line(data = PLB_df, aes(x.PLB, y.PLB), color = gmri_cols("orange"), linewidth = 1) +
  geom_segment(
    data = p_prepped %>% filter(lowCount >0),
    aes(x = wmin_g, 
        xend = wmax_g, 
        y = countGTEwmin, 
        yend = countGTEwmin, 
        color = hare_group),
    #color = gmri_cols("blue"),
    alpha = 0.5,
    linewidth = 0.8) +
  #facet_wrap(~hare_group) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)), limits = c(xmin, xmax)) +
  labs(x = "Weight (g)", y = "Abundance", title = "Regime 2: Actual Vals + Fit")




#####  Plotting Residuals  ####


# Plotting the pieces
r2_3 <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linewidth = 1) +
  geom_segment(
    data = p_prepped %>% filter(lowCount >0),
    aes(x = wmin_g, 
        xend = wmax_g, 
        y = b_resid, 
        yend = b_resid,
        color = hare_group),
    alpha = 0.8,
    linewidth = 0.8) +
  #facet_wrap(~hare_group) +
  scale_color_gmri() +
  scale_y_continuous(labels = comma_format(scale = 1/1e6, suffix = "M")) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)), limits = c(xmin, xmax)) +
  labs(x = "Weight (g)", y = "Residual Abundance", title = "Regime 2: Size Spectra Estimate Residuals")






####  Regime 3: 2007-2019 ####







# Prepare inputs directly:
yr_op <- c(2007:2019)
reg <- "GoM"


# Data for the spectra
actual_vals <- catch_size_bins %>% 
  filter(
    Year %in% yr_op, 
    survey_area == reg) 


# Get the slope over those years
mle_results <- actual_vals  %>% 
  group_isd_estimation(
    min_weight_g = 1, 
    max_weight_g = 10000,
    isd_xmin = 1,
    isd_xmax = 10000,
    abundance_vals = "stratified",
    .group_cols = c("survey_area")) 


# Power law parameters and summary details for the group of data:
b.MLE           <- mle_results$b
total_abundance <- mle_results$n
b.confMin       <- mle_results$confMin
b.confMax       <- mle_results$confMax


# Create range of x values from the group to get power law predictions
# min and max weights for power law
xmin  <- 1
xmax  <- 10^4

# Create x values (individual bodymass) to predict across
# break up the Xlim into pieces between min and max
x.PLB <- seq(
  from = xmin, 
  to   = xmax,
  by = 0.1)  

# get the length of that vector
x.PLB.length <- length(x.PLB)  

# Y values for plot limits/bounds/predictions from bounded power law pdf
y.PLB         <- (1 - sizeSpectra::pPLB(x = x.PLB, b = b.MLE, xmin = min(x.PLB), xmax = max(x.PLB))) * total_abundance
y.PLB.confMin <- (1 - sizeSpectra::pPLB(x = x.PLB, b = b.confMin, xmin = min(x.PLB), xmax = max(x.PLB))) * total_abundance
y.PLB.confMax <- (1 - sizeSpectra::pPLB(x = x.PLB, b = b.confMax, xmin = min(x.PLB), xmax = max(x.PLB))) * total_abundance


# Put it in a df to make it easier
PLB_df <- data.frame(
  x.PLB   = x.PLB,
  y.PLB   = y.PLB,
  confMin = y.PLB.confMin,
  confMax = y.PLB.confMax) %>% 
  mutate()



# 2. solve the power law function again in mutate for residuals
# need to do it twice for wmin_g and wmax_g

# Aggregate across overlapping size ranges
p_prepped <- actual_vals %>% 
  isd_plot_prep(stratified_abundance = T, min_weight_g = 1) %>% 
  left_join(actual_vals) %>% 
  select(comname, hare_group,  wmin_g, wmax_g, lowCount, highCount, countGTEwmin) %>% 
  # Get the residuals
  mutate(
    # Estimated Abundance
    wmin_estimate = (1 - sizeSpectra::pPLB(x = wmin_g, b = b.MLE, xmin = xmin, xmax = xmax)) * total_abundance,
    wmax_estimate = (1 - sizeSpectra::pPLB(x = wmax_g, b = b.MLE, xmin = xmin, xmax = xmax)) * total_abundance,
    # Residuals
    b_resid = countGTEwmin - ((wmin_estimate + wmax_estimate)/2)) 



#####  Plotting Raw  ####

# Plotting the pieces
r3_1 <- ggplot() +
  geom_segment(
    data = p_prepped,
    aes(x = wmin_g, 
        xend = wmax_g, 
        y = countGTEwmin, 
        yend = countGTEwmin),
    color = gmri_cols("blue"),
    alpha = 0.5,
    linewidth = 0.8) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)), limits = c(xmin, xmax)) +
  labs(x = "Weight (g)", y = "Abundance", title = "Regime 3: Actual Vals")

#####  Plotting Functional Groups  ####


# Plotting the pieces
r3_2 <- ggplot() +
  geom_line(data = PLB_df, aes(x.PLB, y.PLB), color = gmri_cols("orange"), linewidth = 1) +
  geom_segment(
    data = p_prepped %>% filter(lowCount >0),
    aes(x = wmin_g, 
        xend = wmax_g, 
        y = countGTEwmin, 
        yend = countGTEwmin, 
        color = hare_group),
    #color = gmri_cols("blue"),
    alpha = 0.5,
    linewidth = 0.8) +
  #facet_wrap(~hare_group) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)), limits = c(xmin, xmax)) +
  labs(x = "Weight (g)", y = "Abundance", title = "Regime 3: Actual Vals + Fit")




#####  Plotting Residuals  ####


# Plotting the pieces
r3_3 <- ggplot() +
  geom_hline(yintercept = 0, color = "black", linewidth = 1) +
  geom_segment(
    data = p_prepped %>% filter(lowCount >0),
    aes(x = wmin_g, 
        xend = wmax_g, 
        y = b_resid, 
        yend = b_resid,
        color = hare_group),
    alpha = 0.8,
    linewidth = 0.8) +
  #facet_wrap(~hare_group) +
  scale_color_gmri() +
  scale_y_continuous(labels = comma_format(scale = 1/1e6, suffix = "M")) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)), limits = c(xmin, xmax)) +
  labs(x = "Weight (g)", y = "Residual Abundance", title = "Regime 3: Size Spectra Estimate Residuals")







#### Comparisons
library(patchwork)
r1_1 / r2_1 / r3_1
r1_2 / r2_2 / r3_2
r1_3 / r2_3 / r3_3


#### The 1kg Hole - Single Species Checks  ####


label_l10_bins <- function(catch_size_bins){
  catch_size_bins <- catch_size_bins %>% 
    mutate(
      weight_bin = case_when(
        ind_weight_kg <= 0.001 ~ "0 - 1g",
        ind_weight_kg <= 0.010 ~ "1 - 10g",
        ind_weight_kg <= 0.100 ~ "10 - 100g",
        ind_weight_kg <= 1.000 ~ "100g - 1kg",
        ind_weight_kg <= 10.000 ~ "1kg - 10kg",
        ind_weight_kg >= 10.00 ~ ">= 10kg" ),
      weight_bin = factor(weight_bin, levels = c(
        "0 - 1g", "1 - 10g", "10 - 100g", "100g - 1kg",
        "1kg - 10kg", ">= 10kg")))
  
  return(catch_size_bins)
}



# Pick a Species to plot
spec_choice <- "atlantic herring"
spec_choice <- "atlantic cod"
spec_choice <- "haddock"
spec_choice <- "spiny dogfish"

catch_size_bins <- catch_size_bins %>% label_l10_bins()


catch_size_bins %>% 
  filter(comname == spec_choice) %>% 
  filter(survey_area == reg) %>% 
  ggplot(aes(Year, wmin_g)) +
  geom_point(
    aes(color = weight_bin, 
        size = strat_total_lwbio_s/1e3), 
    position = position_jitter(width = 0.15, height = 0), 
    alpha = 0.1) +
  scale_color_carto_d(palette = "Burg") +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  labs(title = str_c(str_to_sentence(spec_choice), " Size Size Structure"),
       y = "Individual Weight (g)",
       subtitle = reg,
       color = "Size Class", size = "Area-Adjusted Catch (thousands)")




####  Digging Into Size Classes  ####

##### Large Fish Index  ####


# Label the regime and drop early years
regime_testing <- catch_size_bins %>% 
  filter(Year >= 1982,
         #comname == spec_choice,
         survey_area == "GoM") %>% 
  mutate(
    regime = case_when(
      Year <= 2000 ~ "1982-1999",
      Year > 2006 ~ "2000-2006",
      TRUE ~ "2007-2019"))

# Get 95th percentile as threshold
DescTools::Quantile(regime_testing$length_cm, 
                    weights = regime_testing$strat_total_abund_s, 
                    probs = c(0.05, 0.5, 0.95))

regime_testing <- regime_testing %>% 
  mutate(large = ifelse(length_cm > 71, T, F))



# Get total bio
total_bio <-  regime_testing %>% 
  group_by(survey_area, Year, regime) %>% 
  summarise(total_bio = sum(strat_total_lwbio_s),
            .groups = "drop")
    
# Get large fish bio
lf_bio <- regime_testing %>% 
  group_by(survey_area, Year, regime) %>% 
  summarise(large_bio = sum(strat_total_lwbio_s*large),
            .groups = "drop")

# Join and divide
lfi <- left_join(total_bio, lf_bio) %>% 
  mutate(LFI = large_bio / total_bio)


# Plot the LFI
gom_lfi <- ggplot(lfi, aes(Year, LFI, fill = regime)) +
  geom_col() +
  scale_fill_gmri() +
  scale_y_continuous(limits = c(0,1)) +
  labs(
    x = "Year",
    y = "Large Fish Index",
    title = "Gulf of Maine Large Fish Index",
    fill = "Breakpoint Regimes"
  )

# Plot the size Spectra
gom_spectra <- region_indices %>% 
  filter(survey_area == "GoM") %>% 
  mutate(Year = as.numeric(Year),
         regime = case_when(
           Year <= 2000 ~ "1982-1999",
           Year > 2006 ~ "2000-2006",
           TRUE ~ "2007-2019")) %>% 
  filter(Year >1981) %>% 
  ggplot(aes(Year, b)) +
  geom_path(linewidth = 0.5, color = "darkgray") +
  geom_path(aes(color = regime), linewidth = 0.8) +
  scale_color_gmri() +
  geom_point(aes(color = regime), size = 2) +
  labs(
    y = "Spectra Slope (b)",
    x = "Year",
    title = "Gulf of Maine Size Spectra",
    color = "Breakpoint Regimes")

gom_spectra / gom_lfi

##### Large Species Indicator  ####

# Do we have L infinity for large species indicator?
# Shephard et al 2012
# Can just check by species if any fish ever is over 1.38kg or 85cm
ever_ls <- catch_size_bins %>% 
  split(.$comname) %>% 
  imap_dfr(function(x, y){
    #gets_big = ifelse(any(x$ind_weight_kg > 1.38), T, F)
    gets_big = ifelse(any(x$length_cm > 85), T, F)
    data.frame(
      "comname" = y, 
      large_species = gets_big)})


# Get large fish bio
ls_bio <- regime_testing %>% 
  left_join(ever_ls) %>% 
  group_by(survey_area, Year, regime) %>% 
  summarise(large_species_bio = sum(strat_total_lwbio_s*large_species),
            .groups = "drop")

# Join and divide
lsi <- left_join(total_bio, ls_bio) %>% 
  mutate(LSI = large_species_bio / total_bio)


# Plot the LSI
gom_lsi <- ggplot(lsi, aes(Year, LSI, fill = regime)) +
  geom_col() +
  scale_fill_gmri() +
  scale_y_continuous(limits = c(0,1)) +
  labs(
    x = "Year",
    y = "Large Species Index",
    title = "Gulf of Maine Large Species Index",
    fill = "Breakpoint Regimes"
  )



gom_spectra / gom_lfi / gom_lsi




####  Normalized  ###


gom_spectra_norm <- region_indices %>% 
  filter(survey_area == "GoM") %>% 
  mutate(Year = as.numeric(Year),
         regime = case_when(
           Year <= 2000 ~ "1982-1999",
           Year > 2006 ~ "2000-2006",
           TRUE ~ "2007-2019")) %>% 
  filter(Year >1981) %>% 
  mutate(b_norm = (b - mean(.$b))/ mean(.$b)) %>% 
  ggplot(aes(Year, b_norm)) +
  geom_path(linewidth = 0.5, color = "darkgray") +
  geom_path(aes(color = regime), linewidth = 0.8) +
  scale_color_gmri() +
  geom_point(aes(color = regime), size = 2) +
  labs(
    y = "Spectra Slope (b) Index (z-score)",
    x = "Year",
    title = "Gulf of Maine Size Spectra",
    color = "Breakpoint Regimes")
           




####   Log2 bin Composition Structure  ####
# If we look at the Gulf of Maine,
# and use bins to simplify what we are looking at,
# what would we see animating through the compositions of the bins

####  REMINDER: Can't Add on Log scale, i.e. you can't display stacked bars this way

# Prepare inputs directly:
yr_op <- c(1982:2019)
reg <- "GoM"


# Data for the spectra
actual_vals <- catch_size_bins %>% 
  filter(
    Year %in% yr_op, 
    survey_area == reg) 


# Pull a year
actual_vals %>% 
  mutate(log2_bins = as.numeric(as.character(log2_bins))) %>% 
  filter(Year == 2018) %>% 
  ggplot(aes(x = left_lim, y = stratified_sum_bodymass, fill = hare_group)) +
  geom_col(position = "stack") 




# Debugging how to get groupings to work here...
# Problem is the group aggregation, and preserving the bins using 

aggregate_log2_bins <- function(
    log2_assigned, 
    min_log2_bin = 0, 
    max_log2_bin = 13, 
    bin_increment = 1,
    grouping_vars){
  
  
  # Full Possible Bin Structure
  # Fills in any gaps
  log2_bin_structure <- define_log2_bins(
    log2_min       = min_log2_bin, 
    log2_max       = max_log2_bin, 
    log2_increment = bin_increment)
  
  
  # Capture all the group levels with a cheeky expand()
  #return(grouping_vars)
  if(length(grouping_vars) > 0){
    log2_bin_structure <- log2_bin_structure %>% 
      tidyr::expand(left_lim, distinct(log2_assigned, !!!syms(grouping_vars))) %>% 
      full_join(log2_bin_structure)
  }
  #return(log2_bin_structure)
  
  
  # Get bin breaks
  log2_breaks <- sort(unique(c(log2_bin_structure$left_lim, log2_bin_structure$right_lim)))
  
  
  # Get Totals for bodymass and abundances
  log2_aggregates <- log2_assigned %>% 
    group_by(left_lim, !!!syms(grouping_vars)) %>% 
    summarise(
      observed_abundance   = sum(numlen_adj, na.rm = T),
      observed_weight_g    = sum(wmin_g, na.rm = T),
      stratified_abundance = sum(strat_total_abund_s, na.rm = T),
      stratified_weight_g  = sum(wmin_area_strat, na.rm = T),
      .groups = "drop")
  
  
  # Join back in what the limits and labels are
  # The defined bins and their labels enforce the size limits
  log2_prepped <- left_join(
    x = log2_bin_structure, 
    y = log2_aggregates)
  
  
  #### Fill Gaps with Zero's?? 
  # This ensures that any size bin that is intended to be in use is actually used
  log2_prepped <- log2_prepped %>% 
    mutate(
      observed_abundance   = ifelse(is.na(observed_abundance), 1, observed_abundance),
      observed_weight_g    = ifelse(is.na(observed_weight_g), 1, observed_weight_g),
      stratified_abundance = ifelse(is.na(stratified_abundance), 1, stratified_abundance),
      stratified_weight_g  = ifelse(is.na(stratified_weight_g), 1, stratified_weight_g))
  
  
  #### normalize abundances using the bin widths
  log2_prepped <- log2_prepped %>% 
    mutate(
      normalized_abund = observed_abundance / bin_width,
      norm_strat_abund = stratified_abundance / bin_width,
      # # Remove Bins Where Normalized Biomass < 0? No!
      # normalized_abund = ifelse(normalized_abund < 2^0, NA, normalized_abund),
      # norm_strat_abund = ifelse(norm_strat_abund < 2^0, NA, norm_strat_abund)
    )
  
  # Add de-normalized abundances (abundance * bin midpoint)
  log2_prepped <- log2_prepped %>% 
    mutate(
      denorm_abund = normalized_abund * bin_midpoint,
      denorm_strat_abund = norm_strat_abund * bin_midpoint)
  
  # Return the aggregations
  return(log2_prepped)
  
}



# Testing
aggregate_log2_bins(
  log2_assigned = actual_vals,
  min_log2_bin = 0, 
  max_log2_bin = 13, 
  bin_increment = 1, grouping_vars = "hare_group")

aggregate_log2_bins(
  log2_assigned = actual_vals,
  min_log2_bin = 0, 
  max_log2_bin = 13, 
  bin_increment = 1, grouping_vars = "hare_group") %>% 
  ggplot(aes(left_lim, norm_strat_abund)) +
  geom_col(aes(fill = hare_group), position = "fill")  #+
  #scale_y_log10(labels = trans_format("log10", math_format(10^.x)))



#' @title Plot normalized abundance of Size Spectra
#' 
#' @description Single panel plot of the distribution of abundance on log size bins.
#'
#' @param dat_l2_labelled Catch data prepped using assign_log2_bins
#' @param stratified Stratified or survey abundances
#'
#' @return
#' @export
#'
#' @examples
plot_grouped_log2_ss <- function(dat_l2_labelled, stratified = TRUE, grouping_var){
  
  # Get totals for each bin:
  bins_aggregated <- aggregate_log2_bins(dat_l2_labelled, grouping_vars = "right_lim") %>% drop_na()
  groups_aggregated <- aggregate_log2_bins(dat_l2_labelled, grouping_vars = grouping_var)
  
  #### Plots Correcting for the bin widths
  norm_strat_abund_plot <- bins_aggregated  %>% 
    ggplot(aes(left_lim, norm_strat_abund)) +
    geom_col(data = groups_aggregated, aes(fill = !!sym(grouping_var)), 
             position = "stack", color = "transparent") +
    scale_y_log10(
      labels = trans_format("log2", math_format(2^.x))) + 
    labs(x = "Log2(bodymass)", y = "Abundances Density (normalized)")
  
  
  # Stratified Abundances
  p1 <- norm_strat_abund_plot +
    geom_smooth(
      formula = y ~ x,
      method = "lm",
      color = gmri_cols("orange")) +
    stat_poly_eq(
      formula = y ~ x,
      aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
      label.y = 1.1, 
      parse = TRUE) +
    labs(caption = "Abundance normalized by bin-width.")
  
  # Return plot
  return(p1)
}


# Test that...
plot_grouped_log2_ss(dat_l2_labelled = actual_vals, stratified = T, grouping_var = "hare_group")
