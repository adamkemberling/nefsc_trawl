####  NEFSC Trawl Survey Size Spectrum  ####
# Quick and dirty NEFSC trawl data run for quick display with app

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


####  Data  ####

# NEFSC Trawl Data
nefsc_trawl <- read_csv(here("data/NEFSC/2019survdat_nye.csv"), guess_max = 1e5, col_types = cols()) %>% clean_names() 

nefsc_trawl <- nefsc_trawl %>% 
    mutate(comname = tolower(comname),
           id = as.character(id))


# # Length Weight Relationships via fishbase, compiled by Ben Resek, REU
# lw_coef <- read_csv("data/MENH/listfishusingfishbase.csv", col_types = cols()) %>% 
#   clean_names() %>% select(1:5) %>% arrange(common_name) %>% 
#   mutate(comname = tolower(common_name),
#          common_name = NULL,
#          num_obs = NULL)

# List of distinct common names that were gonna need to operate with
# # Merge with Ben's
# nefsc_list <- nefsc_trawl %>% 
#   distinct(comname) %>% 
#   arrange(comname) %>% 
#   left_join(lw_coef) %>% 
#   mutate(hare_group = NA,
#          units = NA)
# 
# # Export as Master
# write_csv(nefsc_list, here("data/NEFSC/nefsc_lw_key.csv"))




# Load the lw key with coefficients from fishbase
nefsc_lw <- read_csv(here("data/NEFSC/nefsc_lw_key_filled.csv"),
                     guess_max = 1e3,
                     col_types = cols())
  
# Fill in NA growth coefficients for related species
nefsc_lw <- nefsc_lw %>% 
  mutate(a = as.numeric(a),
         b = as.numeric(b),
         ln_a = log(a),
         related_ln_a = log(related_a),
         ln_a = ifelse(is.na(ln_a), related_ln_a, ln_a))














####_____________________####
####  Setup  ####


nefsc_trawl %>% glimpse()

nefsc <- nefsc_trawl %>% 
  select(c(id,
           est_year:svvessel, 
           decdeg_beglat, 
           decdeg_beglon, 
           comname:area, 
           station:numlen))


glimpse(nefsc)


####  Individual Lengths  ####
nefsc_lens <- nefsc %>% 
  filter(is.na(numlen) == FALSE,
         numlen > 0) %>% 
  #mutate(id = str_c(station, comname, length, sep = "_")) %>% 
  #group_by(id, numlen) %>% 
  distinct(id, comname, length, numlen) #%>% uncount(weights =  numlen)


# recombine with the distinct station info
nefsc_spectra <- nefsc %>% 
  select(-c(est_julian_day, comname, numlen, biomass, abundance, length, numlen)) %>% 
  distinct() %>% 
  left_join(nefsc_lens, by = "id")



# Combine with the length weight coefficients
####  Calculate Weights  ####
nefsc_weights <- nefsc_spectra %>% 
  left_join(nefsc_lw[,1:5], by = "comname")  %>% 
  mutate(b = as.numeric(b),
         a = as.numeric(a),
         ln_a = log(a),
         ln_weight = (ln_a + b * log(length)),
         weight_g = exp(ln_weight),
         freq_weight = weight_g * numlen) %>% 
  drop_na(weight_g)


# Check weights
nefsc_weights %>% 
  group_by(comname) %>% 
  summarise(mean_weight =  mean(weight_g)) %>% 
  arrange(desc(mean_weight)) 


nefsc_weights %>% 
  ggplot( aes( y = fct_reorder(comname, weight_g, .fun = mean, .desc = TRUE), x = weight_g)) + 
  geom_boxplot() + 
  labs(x = "Mean Weight (g)", y = "Common Name")





####__  SizeSpectra Setup  ####

# Starting point for sizeSpectra steps
dataOrig <- nefsc_weights %>% filter(is.na(weight_g) == FALSE)


# Keep desired columns, name them for the vignette
data <- dataOrig %>%  
  mutate(season = factor(season, levels = c("SPRING", "SUMMER", "FALL", "WINTER"))) %>%
  select(
    Year = est_year,                     # year
    samp_month = est_month,
    samp_dat = est_day,
    samp_time = est_time,
    vessel = svvessel,
    season,
    SpecCode = comname,          # common name
    hare_group,
    LngtClass = length,              # length bin
    Number = numlen,              # CPUE in n/effort
    LWa = a,                         # length/weight param
    LWb = b,                         # length/weight param
    bodyMass = weight_g,             # weight of an individual from that size class
    CPUE_bio_per_hour = freq_weight  # CPUE in biomass
  ) %>% arrange(Year, season, SpecCode, LngtClass)




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
paste0("Total Number of Fish in Analysis: ", sum(data$Number))

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
  geom_line(aes(Year, uniqLngtClass, color = "Unique Length Classes")) +
  geom_line(aes(Year, uniqSpec, color = "Unique Species")) +
  labs(y = "Count", x = NULL) +
  facet_wrap(~season, ncol = 1)





####__  Size Spectra Bins  ####

species_splits <- data_biomass %>% 
  split(.$SpecCode)


# Get a key for each species and size class
# assumes 1cm bins for all species as written
data_bin_key <- species_splits %>% 
  map_dfr(function(species_df){
    
    #pull the distinct length bins
    species_df <- species_df %>% 
      distinct(LngtClass, .keep_all = T) %>% 
      arrange(LngtClass)
    
    
    # Add the max length for the bin, and its weight
    binned_df <- species_df %>% 
      mutate(
        LngtMax = LngtClass + 1, 
        wmax    = exp(log(LWa) + LWb * log(LngtMax)))  %>%
      select(SpecCode, LWa, LWb, LngtMin = LngtClass, wmin = bodyMass, LngtMax, wmax,
             -c(Number, Biomass))
    
    # return the clean data
    return(binned_df)})


# Add the bins back into the original and clean up
dataBin <- data_biomass %>% 
  select(Year, season, SpecCode, LngtMin = LngtClass, Number, Biomass) %>% 
  left_join(data_bin_key, by = c("SpecCode", "LngtMin")) 





####  Export NEFSC dataBin  ####
write_csv(dataBin, here("data/NEFSC/nefsc_ss_bins.csv"))


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
         season = str_sub(group_var, 6, -1),
         season = factor(season, levels = c("SPRING", "SUMMER", "FALL", "WINTER")))


# change dumb name
mle_seasons <- mle_bins


# Better plot
mle_seasons %>% 
  filter(season %in% c("SPRING", "FALL")) %>% 
  ggplot(aes(Year, b, color = season)) +
  geom_line(aes(group = season)) +
  geom_segment(aes(x = Year, xend = Year, y = confMin, yend = confMax)) +
  geom_point(aes(y = b)) + 
  geom_smooth(formula = y ~ x,
              method = "lm") +
  labs(x = NULL,
       y = "Size Spectrum Slope (b)") +
  facet_wrap(~season, ncol = 1)



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
seasonal_isd$`2005_SPRING`$stacked | seasonal_isd$`2005_FALL`$stacked
seasonal_isd$`1975_SPRING`$stacked | seasonal_isd$`1975_FALL`$stacked


seasonal_isd


map(seasonal_isd, function(x){ x[[3]]})
