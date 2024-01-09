####  Single Species Spectrum Checking
# Check fishery targeted species against landings trends

# Grab Cod

# Get Cod only spectra timeseries

# Plot that, intercept and slope




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
library(readxl)

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

# 1. trawl data used as input
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




####  Log2 Bin Cod Spectra  ####

# Filter to just cod
just_cod <- filter(catch_log2_labelled, 
                   comname == "atlantic cod",
                   survey_area == "GoM")




# Use this function on each year for Gulf of Maine
cod_spectra <- group_log2_spectra(
  wmin_grams = just_cod, 
  min_log2_bin = 0, 
  max_log2_bin = 10, 
  min_weight_g = 1, 
  bin_increment = 1, 
  .group_cols = c("Year", "survey_area"))



#####  log2 Spectra Timeseries  ####

# Plot slope
cod_spectra %>% 
  mutate(Year = as.numeric(Year)) %>% 
  ggplot(aes(Year, log2_slope_strat)) +
  geom_point() +
  geom_line() +
  facet_wrap(~survey_area)

# Plot intercept
cod_spectra %>% 
  mutate(Year = as.numeric(Year)) %>% 
  ggplot(aes(Year, log2_int_strat)) +
  geom_point() +
  geom_line() +
  facet_wrap(~survey_area)



#####  Single Year Bins  ####


# Plot one
plot_log2_ss <- function(dat_l2_labelled, stratified = TRUE, remove_empty_bins = FALSE){
  
  # Get totals for each bin:
  bins_aggregated <- aggregate_log2_bins(dat_l2_labelled)
  
  if(remove_empty_bins){bins_aggregated <- filter(bins_aggregated, norm_strat_abund>1)}
  
  #### Plots Correcting for the bin widths
  norm_strat_abund_plot <- bins_aggregated  %>% 
    ggplot(aes(left_lim, norm_strat_abund)) +
    geom_col(color = gmri_cols("green"), fill = gmri_cols("green")) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) + 
    scale_x_continuous(labels = math_format(2^.x)) + 
    labs(x = "Bodyweight (g)", y = "Abundances Density (normalized)")
  
  
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


# Single year plot
plot_log2_ss(dat_l2_labelled = filter(just_cod, Year == 2000), stratified = T) +
  labs(title = "Atlantic Cod Size Spectrum, 2000, Gulf of Maine",
       x = "Bodyweight (g)")

# Try it again with TRUE NA's & less than 1 abundance as NA

plot_log2_ss(dat_l2_labelled = filter(just_cod, Year == 2000), stratified = T, remove_empty_bins = T) +
  labs(title = "Atlantic Cod Size Spectrum, 2000, Gulf of Maine",
       x = "Bodyweight (g)")

# All years
plot_log2_ss(dat_l2_labelled = just_cod, stratified = T) +
  labs(title = "Atlantic Cod Size Spectrum, All-Years 1970-2019, Gulf of Maine",
       x = "Bodyweight (g)")


### MLE Cod Spectra  ####


# Get the slope over those years
year_test <- just_cod %>% 
  filter(Year == 2000)

# use MLE method
mle_results <- year_test  %>% 
  group_isd_estimation(
    min_weight_g = 1, 
    max_weight_g = 10000,
    isd_xmin = 1,
    isd_xmax = 10000,
    abundance_vals = "stratified",
    .group_cols = c("survey_area")) 

# Prepare pieces for Comparable MLE plot
# Power law parameters and summary details for the group of data:
b.MLE           <- mle_results$b
total_abundance <- mle_results$n
b.confMin       <- mle_results$confMin
b.confMax       <- mle_results$confMax


# Create range of x values from the group to get power law predictions
# min and max weights for power law
xmin  <- 1
xmax  <- 10^4


# Aggregate across overlapping size ranges to get bar height for observations
p_prepped <- year_test %>% 
  isd_plot_prep(stratified_abundance = T, min_weight_g = 1) %>% 
  left_join(year_test) %>% 
  select(comname, hare_group,  wmin_g, wmax_g, lowCount, highCount, countGTEwmin)

# Make prediction over the range
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


# Plotting the pieces
ggplot() +
  geom_ribbon(data = PLB_df, aes(x = x.PLB,ymin = confMin, ymax = confMax), fill = "gray90") +
  geom_line(data = PLB_df, aes(x = x.PLB, y = y.PLB), color = "orange", linewidth = 1) +
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
  annotate(x = 10^1, y = 10^3, label = )
  labs(x = "Bodyweight (g)", y = "Abundance", 
       title = "Atlantic Cod Size Distribution, 2000, Gulf of Maine")


  
  
# Do all years for Gulf of Maine Cod using MLE,
# should be less susceptible to these biases of absent size classes

  
cod_mle_results <- just_cod  %>% 
  group_isd_estimation(
    min_weight_g = 1, 
    max_weight_g = 10000,
    isd_xmin = 1,
    isd_xmax = 10000,
    abundance_vals = "stratified",
    .group_cols = c("Year", "survey_area"))  %>% 
  mutate(year = as.numeric(Year))

  

# Plot the timeseries
cod_mle_results %>% 
  ggplot(aes(year, b)) +
  geom_line() +
  geom_point() +
  labs(x = "Year", y = "Exponent of Individual Size Distribution (b)",
       title = "Atlantic Cod Size Distribution Exponent, 1970-2019, Gulf of Maine")







####  Cod Spectra and Fishing  ####

# Landings of finfish* sheet 5
res_path <- cs_path("res")

landings <- read_xlsx(
  path = str_c(res_path, "GARFO_landings/KMills_landings by area 1964-2021_JUN 2022.xlsx"), sheet = 5) %>% 
  rename_all(tolower)



# 1. Aggregate Statzones by the region definitions of the project

# Make a list of zones to roughly match the survey areas:
fish_zones <- list(
  "Gulf of Maine" = c(511:515, 464, 465),
  "Georges Bank" = c(521, 522, 525, 561, 562),
  "Southern New England" = c(611, 612, 613, 616, 526, 537, 538, 539),
  "Mid-Atlantic Bight" = c(614:615, 621, 622, 625, 626, 631, 632))


# Join to landings
# Add the labels into the landings data and remove what we don't need there:
area_levels_long <- c("Northeast Shelf", "Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight")
landings <- landings %>% 
  mutate(
    survey_area = case_when(
      `stat area` %in% fish_zones$"Gulf of Maine" ~ "Gulf of Maine",
      `stat area` %in% fish_zones$"Georges Bank" ~ "Georges Bank",
      `stat area` %in% fish_zones$"Southern New England" ~ "Southern New England",
      `stat area` %in% fish_zones$"Mid-Atlantic Bight" ~ "Mid-Atlantic Bight")) %>% 
  filter(survey_area %in% c("Georges Bank", "Gulf of Maine", "Southern New England", "Mid-Atlantic Bight")) %>% 
  mutate(survey_area = factor(survey_area, area_levels_long)) %>% 
  filter(between(year, 1960, 2019),
         !str_detect(sppname, "CONFIDENTIAL")) %>% 
  mutate(decade = floor_decade(year),
         sppname = str_to_title(sppname))




#####  Landings Breakdown by Groundfish  ####
# Goal:
# Separate landings into groups that might better explain in a mechanistic sense
# how they would impact a trawl-sampled community


# Maybe a plot is better:
landings_yrly <- landings %>% 
  expand(year, survey_area, sppname) %>% 
  left_join(landings) %>% 
  filter(survey_area!= "Northeast Shelf") %>% 
  group_by(survey_area, year,  sppname) %>% 
  summarise(avg_landings_lb = mean(`landed lbs`, na.rm = T),
            across(.cols = c(landings_lb = `landed lbs`, 
                             value = value), 
                   .fns = sum, 
                   .names = "total_{.col}"),
            .groups = "drop") %>% 
  mutate(across(4:6, ~ifelse(is.na(.x), 0, .x)))



# Assign groundfish label (Jonathan Labaree, market based)
gf_species <- tolower(
  c("cod, atlantic", 
    "Flounder, American Plaice", 
    "Flounder, Winter", 
    "Flounder, witch",
    "flounder, yellowtail", 
    "haddock", 
    "hake, white", 
    "hake, red", 
    "hake, silver", 
    "halibut, atlantic", 
    "pollock", 
    "redfish, acadian"))


# What does that give us if the groundfish are labeled?
landings_yrly <- mutate(
  landings_yrly, is_gf = ifelse(tolower(sppname) %in% gf_species, "groundfish", "other"))

# Gulf of Maine only
gom_landings <- filter(landings_yrly, survey_area == "Gulf of Maine")



# Plots
gom_landings %>% 
  ggplot() +
  geom_line(aes(year, total_landings_lb, group = sppname, color = is_gf), show.legend = F) +
  #facet_grid(survey_area~is_gf)
  facet_grid(is_gf~survey_area) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(x = "Year", y = "Single Species Yearly Total Landings (lb.)")



# Join landings to spectra
gom_fishing_cod <- gom_landings %>% 
  #select(-c(avg_landings_lb, total_value)) %>% 
  pivot_wider(values_from = "total_landings_lb", names_from = "is_gf") %>% 
  group_by(year) %>% 
  summarise(gf_landings = sum(groundfish, na.rm = T),
            other_landings = sum(other, na.rm = T)) %>% 
  left_join(cod_mle_results, join_by(year))




# Are groundfish landings correlated?
gom_fishing_cod %>% 
  ggplot(aes(gf_landings, b)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_poly_eq(
    formula = y ~ x,
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
    label.y = 1.1, 
    parse = TRUE) +
  labs(x = "Groundfish Species Landings (lb.)", y = "Cod Size Exponent of Distribution (b)")


# What about lags
# Run the ccf
ccf(
  drop_na(gom_fishing_cod) %>% pull(gf_landings), 
  drop_na(gom_fishing_cod) %>% pull(b), 
  plot= T , lag.max = 10, 
  main = "Groundfish Landings Lagged Impact on Spectrum Slope") 


# Are all landings correlated?
gom_fishing_cod %>% 
  ggplot(aes(other_landings, b)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_poly_eq(
    formula = y ~ x,
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
    label.y = 1.1, 
    parse = TRUE)  +
  labs(x = "Non-Groundfish Species Landings (lb.)", y = "Cod Size Exponent of Distribution (b)")

ccf(
  drop_na(gom_fishing_cod) %>% pull(other_landings), 
  drop_na(gom_fishing_cod) %>% pull(b), 
  plot= T , lag.max = 10, 
  main = "Non-Groundfish Landings Lagged Impact on Spectrum Slope") 





# Only species in the trawl survey?
all_species <- distinct(catch_log2_labelled,comname)
gom_landings$sppname[!landings$sppname %in% all_species]












# #####______________####
# ####  Playing with Log 10 Bins  ####
# 
# # Re-label the bin information
# just_cod_10 <- just_cod %>% 
#   select(-c(left_lim, right_lim, bin_label, bin_midpoint, log2_bins, log2_weight, bin_width)) %>% 
#   assign_log10_bins(l10_increment = 1)
# 
# 
# 
# # Use this function on each year for Gulf of Maine
# cod_l10_spectra <- group_l10_spectra(
#   wmin_grams = just_cod_10, 
#   min_l10_bin = 0, 
#   max_l10_bin = 4, 
#   min_weight_g = 1, 
#   bin_increment = 1, 
#   .group_cols = c("Year", "survey_area"))
# 
# 
# 
# ####  Spectra Timeseries  ####
# 
# # Plot slope
# cod_l10_spectra %>% 
#   mutate(Year = as.numeric(Year)) %>% 
#   ggplot(aes(Year, l10_slope_strat)) +
#   geom_point() +
#   geom_line() +
#   facet_wrap(~survey_area)
# 
# # Plot intercept
# cod_l10_spectra %>% 
#   mutate(Year = as.numeric(Year)) %>% 
#   ggplot(aes(Year, l10_int_strat)) +
#   geom_point() +
#   geom_line() +
#   facet_wrap(~survey_area)
# 
# 
# # Single Year
# plot_log10_ss(l10_assigned = filter(just_cod_10, Year == 2000)) + 
#   labs(title = "Atlantic Cod Size Spectrum, 2000, Gulf of Maine")
# 
# # Single Year
# plot_log10_ss(l10_assigned = just_cod_10) + 
#   labs(title = "Atlantic Cod Size Spectrum, All-Years, Gulf of Maine")
