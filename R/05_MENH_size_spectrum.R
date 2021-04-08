####  ME/NH Trawl Survey Size Spectrum  ####


####_____________________####
####  Packages  ####
library(here)
library(gmRi)
library(janitor)
library(sizeSpectra)
library(patchwork)
library(tidyverse)

####  Support Functions  ####
source(here("R/support/sizeSpectra_support.R"))
source(here("R/support/maine_nh_trawl_build.R"))

# Paths
res_path <- shared.path(os.use = "unix", group = "res data", folder = "")

####  Data  ####

# https://mainedmr.shinyapps.io/MaineDMR_Trawl_Survey_Portal/
# menh_len <- read_csv(here("data/MENH/shiny_length_freq.csv"), guess_max = 1e6, col_types = cols()) %>% clean_names()

# Load length frequencies using build code
menh_len <- load_menh_data("length frequencies")








####_____________________####
####  Setup  ####

###__  Calculating Weights  ####

# Using Ben's Coefficients

# # Length Weight Relationships via fishbase, compiled by Ben Resek, REU
# lw_coef <- read_csv("data/MENH/listfishusingfishbase.csv", col_types = cols()) %>% 
#   clean_names() %>% select(1:5) %>% 
#   arrange(common_name)
# 
# menh_weights <- left_join(menh_len, lw_coef, by = "common_name") %>% 
#   select(common_name, ln_a, a, b, everything())  %>% 
#   mutate(ln_weight = (ln_a + b * log(length)),
#          weight_kg = exp(ln_weight),
#          freq_weigth = weight_kg * frequency)

# Using Wigley & Fishbase
menh_weights <- add_lw_to_menh(menh_length_frequencies = menh_len, os.use = "unix")


####__  Checking Weights  ####

# filter to only use cm
menh_weights <- menh_weights %>% 
  filter(unit_of_length == "CM",
         comname != "sea urchins green",
         ind_weight_kg > 0) %>% 
  mutate(comname = fct_drop(comname))

# test species
test_species <- c("alewife", "atlantic herring", "atlantic cod", "alligatorfish", "smelt rainbow")
menh_weights %>% 
  filter(comname %in% test_species) %>% 
  ggplot(aes(length, ind_weight_kg)) +
  geom_point() +
  facet_wrap(~comname)

# Visualizing average weights by species
menh_weights %>% 
  ggplot( aes(y = fct_reorder(comname, ind_weight_kg, .fun = mean, .desc = TRUE), 
              x = ind_weight_kg)) + 
  geom_boxplot() + 
  labs(x = "Mean Weight (kg)", y = "Common Name")


# Got some gram measurements baked in:
# Species that are in grams
gram_species <- c(
  "Monkfish",
  "Lumpfish",
  "Squid Short-Finned",
  "Sculpin Shorthorn",
  "Crab Jonah",
  "Wrymouth",
  "Fourbeard Rockling",
  "Crab Atlantic Rock",
  "Blenny Snake",
  "Smelt Rainbow",
  "Sand Lance American",
  "Menhaden Atlantic",
  "Silverside Atlantic",
  "Alligatorfish"
)


# Change those species to kg to be consistent
menh_weights <- menh_weights %>% 
  mutate(weight_kg = ifelse(comname %in% gram_species, weight_kg/1000, weight_kg),
         freq_weight = weight_kg * frequency)

# check plot again
menh_weights %>% 
  ggplot( aes( y = fct_reorder(comname, weight_kg, .fun = mean, .desc = TRUE), x = weight_kg)) + 
  geom_boxplot() + 
  labs(x = "Mean Weight (kg)", y = "Common Name")


####__  SizeSpectra Setup  ####

# Starting point for sizeSpectra steps
dataOrig <- menh_weights %>% filter(is.na(weight_kg) == FALSE)


# Keep desired columns, name them for the vignette
data <- dataOrig %>%  
  mutate(
    a = exp(ln_a),
    weight_g = weight_kg * 1000) %>% 
  select(
    Year = year,                           # year
    SpecCode = comname,                # common name
    LngtClass = length,                    # length bin
    Number = frequency,                    # CPUE in n/effort
    LWa = a,                               # length/weight param
    LWb = b,                               # length/weight param
    bodyMass = weight_g,                   # weight of an individual from that size class
    CPUE_bio_per_hour = freq_weight        # CPUE in biomass
  ) %>% arrange(Year, SpecCode, LngtClass)



####__  Aggregate into Groups Year/Area etc.  ####
data <- data %>% 
  group_by(Year, SpecCode, LngtClass) %>% 
  summarise(
    Number = sum(Number), #/grouping_variable (numAreas),
    LWa = unique(LWa),
    LWb = unique(LWb),
    bodyMass = unique(bodyMass)) %>% 
  ungroup()


# Total number of fish 
total_fish <- sum(data$Number)

# Biomass data
data_biomass <- data %>% 
  mutate(Biomass = Number * bodyMass)



####__  Track Unique Species and Length Classes  ####
dataSumm <- data %>% 
  group_by(Year) %>% 
  summarise(uniqLngtClass = length(unique(LngtClass)),
            uniqSpec = length(unique(SpecCode)))
ggplot(dataSumm) +
  geom_line(aes(Year, uniqLngtClass, color = "Unique Length Classes")) +
  geom_line(aes(Year, uniqSpec, color = "Unique Species")) +
  labs(y = "Count", x = "Null")



####__  Size Spectra L/W Bins  ####

species_splits <- data_biomass %>% 
  split(.$SpecCode)


# Get a key for each species and size class
# assumes 1cm bins for all species as written
# returns what the wmin and wmax is based on length weight relationship
data_bin_key <- species_splits %>% 
  map_dfr(make_bin_key)

# Check Alewife, a species that was fine
data_bin_key %>% filter(SpecCode == "Alewife")  %>% arrange(LngtMin) %>% head()

# Add the bins back into the original and clean up
dataBin <- data_biomass %>% 
  select(Year, SpecCode, LngtMin = LngtClass, Number, Biomass) %>% 
  left_join(data_bin_key, by = c("SpecCode", "LngtMin")) 

# Create Grouping Var
dataBin <- dataBin %>% 
  mutate(group_var = Year)


####__  Set Bodymass Cutoffs  ####

# Set bodymass lower limit
# Filter for lower end of gear selectivity
mass_cutoff <- 400 #grams
dataBin <- filter(dataBin, wmin >= mass_cutoff)




####_____________________####
####  Estimation  ####



####__  MLEbins Size Spectra  ####

# Map through instead of looping
# Returns Power Law Coefficients for the Each group
mle_bins <- dataBin %>% 
  split(.$group_var) %>% 
  imap_dfr(group_mle_calc) 




# Need the standard error for weighted linear regression,
#  see eightMethods.count() for details:
mle_bins <- mle_bins %>% 
  mutate(stdErr = (abs(confMin-b) + abs(confMax - b)) / (2 * 1.96),
         Year = group_var,
         Year = as.numeric(as.character(Year)))


# change name
mle_years <- mle_bins

####__  Plotting Size Spectra  ####

# Better plot
mle_years %>% 
  ggplot(aes(Year, b)) +
  geom_line(group = 1) +
  geom_segment(aes(x = Year, xend = Year, y = confMin, yend = confMax)) +
  geom_point(aes(y = b)) + 
  geom_smooth(formula = y ~ x,
              method = "lm") +
  labs(x = NULL,
       y = "Size Spectrum Slope (b)")


####__  Summary Table  ####

MLEbins.res <- mle_years %>%  
  mutate(C = (b != -1 ) * (b + 1) / ( xmax^(b + 1) - xmin^(b + 1) ) + (b == -1) * 1 / ( log(xmax) - log(xmin)))



####_____________________####
####  Plotting w/  ggplot  ####


# Data inputs
dataRecommend.isd <- dataBin %>% select(group_var, Year, wmin, wmax, Number) # Recommended data for ISD plots
MLEbins.res.list  <- MLEbins.res %>% split(.$Year)                    # MLE coefficients that match them
data.year.list <- dataRecommend.isd %>% 
  split(.$group_var) %>% 
  map(isd_plot_prep)

# global limits for plot
xlim.global <- c( min(dataRecommend.isd$wmin), max(dataRecommend.isd$wmax) )

# Vector of years, for naming
group_names <- sort(unique(dataBin$group_var))


# Loop through years for plots
annual_isd <- map2(.x = data.year.list,
                   .y =  MLEbins.res.list, 
                   .f = ggplot_isd, 
                   show_pl_fit = TRUE,
                   plot_rects = FALSE) %>% 
  setNames(group_names)

# Take a peak
annual_isd$`2000`$stacked |annual_isd$`2005`$stacked




####_____________________####
####_____________________####
####_____________________####
#### Year x Season Slopes ####


####__ Setup  ####

# Same starting Data
dataOrig <- menh_weights %>% filter(is.na(weight_kg) == FALSE)

#number of areas/stratum/seasons that we need to aggregate on later
numAreas  <- length(unique(dataOrig$region))


# Keep desired columns, name them for the vignette
data <- dataOrig %>%  
  mutate(
    a = exp(ln_a),
    weight_g = weight_kg * 1000) %>% 
  select(
    Year = year,                     # year
    season = season,                 # season
    SpecCode = comname,          # common name
    LngtClass = length,              # length bin
    Number = frequency,              # CPUE in n/effort
    LWa = a,                         # length/weight param
    LWb = b,                         # length/weight param
    bodyMass = weight_g,             # weight of an individual from that size class
    CPUE_bio_per_hour = freq_weight  # CPUE in biomass
  ) %>% arrange(Year, SpecCode, LngtClass)

####__  Aggregate into Groups Year/Area etc.  ####
data <- data %>% 
  group_by(Year, season, SpecCode, LngtClass) %>% 
  summarise(
    #Divide Number / grouping_variable (numAreas),
    Number = sum(Number), 
    LWa = unique(LWa),
    LWb = unique(LWb),
    bodyMass = unique(bodyMass)) %>% 
  ungroup()



# Total number of fish 
total_fish <- sum(data$Number)

# Biomass data
data_biomass <- data %>% 
  mutate(Biomass = Number * bodyMass)



####__  Track Unique Species and Length Classes  ####
dataSumm <- data %>% 
  group_by(Year, season) %>% 
  summarise(uniqLngtClass = length(unique(LngtClass)),
            uniqSpec = length(unique(SpecCode))) %>% 
  ungroup()

ggplot(dataSumm) +
  geom_line(aes(Year, uniqLngtClass, color = "Unique Length Classes", linetype = season)) +
  geom_line(aes(Year, uniqSpec, color = "Unique Species", linetype = season)) +
  labs(y = "Count", x = NULL)



####__  Size Spectra Bins  ####

species_splits <- data_biomass %>% 
  split(.$SpecCode)


# Get a key for each species and size class
# assumes 1cm bins for all species as written
data_bin_key <- species_splits %>% 
  map_dfr(make_bin_key)


# Add the bins back into the original and clean up
dataBin <- data_biomass %>% 
  select(Year, season, SpecCode, LngtMin = LngtClass, Number, Biomass) %>% 
  left_join(data_bin_key, by = c("SpecCode", "LngtMin")) 





####  Export MENH dataBin  ####
write_csv(dataBin, here("data/MENH/menh_ss_bins.csv"))


####__  Set Bodymass Cutoff and Groups  ####

# Set bodymass lower limit
# Filter for lower end of gear selectivity
mass_cutoff <- 400 #grams
dataBin <- filter(dataBin, wmin >= mass_cutoff)

# Create Grouping Var
dataBin <- dataBin %>% 
  mutate(group_var = str_c(Year, season, sep = "_"))




####_____________________####
####  Estimation  ####

# Map through instead of looping
mle_bins <- dataBin %>% 
  split(.$group_var) %>% 
  imap_dfr(group_mle_calc)  

mle_bins <- mle_bins %>% 
  mutate(stdErr = (abs(confMin - b) + abs(confMax - b)) / (2 * 1.96),
         Year = str_sub(group_var, 1, 4),
         Year = as.numeric(as.character(Year)),
         season = str_sub(group_var, 6, -1))

# change dumb name
mle_seasons <- mle_bins


# Better plot
mle_seasons %>% 
  ggplot(aes(Year, b, color = season)) +
  geom_line(aes(group = season)) +
  geom_segment(aes(x = Year, xend = Year, y = confMin, yend = confMax)) +
  geom_point(aes(y = b)) + 
  geom_smooth(formula = y ~ x,
              method = "lm") +
  labs(x = NULL,
       y = "Size Spectrum Slope (b)")


####__  Summary Table  ####
MLEbins.res <- mle_seasons %>%  
  mutate(C = (b != -1 ) * (b + 1) / ( xmax^(b + 1) - xmin^(b + 1) ) + (b == -1) * 1 / ( log(xmax) - log(xmin)))

####_____________________####
####  Plotting Annual Curves, with ggplot  ####


# Data inputs
dataRecommend.isd <- dataBin %>% select(group_var, Year, season, wmin, wmax, Number) # Recommended data for ISD plots
MLEbins.res.list  <- MLEbins.res %>% split(.$group_var)                    # MLE coefficients that match them
data.year.list <- dataRecommend.isd %>% 
  split(.$group_var) %>% 
  map(isd_plot_prep)

# global limits for plot
xlim.global <- c( min(dataRecommend.isd$wmin), max(dataRecommend.isd$wmax) )

# Vector of years, for naming
group_names <- sort(unique(dataBin$group_var))


# Loop through years for plots
seasonal_isd <- map2(data.year.list, MLEbins.res.list, .f = ggplot_isd) %>% setNames(group_names)

# Take a peak
seasonal_isd$`2005_Spring`$stacked | seasonal_isd$`2005_Fall` $stacked


seasonal_isd


map(seasonal_isd, function(x){ x[[3]]})
