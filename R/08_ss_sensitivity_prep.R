####

# Goal: attribute ocean region to nefsc trawl size spectra
# Do areas of inherently different temperature reflect differences in slope/shape
# Is this holding up, does GOM resemble mid atlantic now
# Want to end up with a chronological record of what the slope is for every year/season/region and how that region did for the entire year as well


# Output files:
# dataBin & results from 5g cutoff are saved for use in 09_ss_sensitivity


####  Packages  ####
library(here)
library(gmRi)
library(sf)
library(janitor)
library(sizeSpectra)
library(patchwork)
library(tidyverse)



####  Support Functions  ####
source(here("R/support/sizeSpectra_support.R"))

group_mle_plot <- function(mle_res){
  plot <- mle_res %>% 
    ggplot(aes(Year, b, color = area, shape = season)) +
    geom_pointrange(aes(x = Year, y = b, ymin = confMin, ymax = confMax)) +
    labs(x = "Year",
         y = "Size Spectrum Slope (b)") 
  
  return(plot)
}

# File Paths
mills_path <- shared.path(os.use = "unix", group = "Mills Lab", folder = NULL)
res_path   <- shared.path(os.use = "unix", group = "RES Data", folder = NULL)





####_____________________________________________________####
####  Load Prepped Data  ####
nefsc_weights <- read_csv(here::here("data/ss_prepped_data/survdat_2020_ss.csv"),
                          col_types = cols(),
                          guess_max = 1e5)


# Check weights - seems like some are in grams
nefsc_weights %>% group_by(comname) %>% summarise(mean_weight =  mean(ind_weight_kg)) %>% arrange(desc(mean_weight)) 
nefsc_weights %>% group_by(comname) %>% summarise(mean_weight =  mean(ind_weight_kg)) %>% arrange(mean_weight)


#plot them all
weights_summ <- nefsc_weights %>% 
  group_by(comname) %>% 
  summarise(mean_weight = mean(ind_weight_kg, na.rm = T),
            bin = case_when(
              mean_weight <= .25 ~ "0 - .25",
              mean_weight <= 5 ~ ".25 - .5",
              mean_weight <= 1 ~ ".5 - 1",
              mean_weight <= 5 ~ "1 - 5",
              mean_weight <= 10 ~ "5 - 10",
              mean_weight <= 50 ~ "10 - 50",
              TRUE ~ "50+")) %>% 
  left_join(nefsc_weights, by = "comname") 

weights_summ %>% 
  #filter(weight_g > 5) %>% 
  ggplot( aes( y = fct_reorder(comname, ind_weight_kg, .fun = mean, .desc = TRUE), x = ind_weight_kg)) + 
  geom_boxplot(outlier.alpha = 0.2, outlier.size = 0.5, outlier.shape = 3) + 
  labs(x = "Mean Individual Weight (kg)", y = "Common Name") +
  facet_wrap(~bin, scales = "free", ncol = 2) +
  theme(axis.text = element_text(size = 6))







 ####_____________________________________________________####
 
 ####  Size Spectra Prep  ####
 
 
 # Starting point for sizeSpectra steps
 dataOrig <- nefsc_weights %>% filter(is.na(ind_weight_kg) == FALSE)
 
 
 # Keep desired columns, name them for the vignette
 data <- dataOrig %>%  
   mutate(season = factor(season, levels = c("SPRING", "SUMMER", "FALL", "WINTER"))) %>%
   rename(
     Year = est_year,                 # year
     SpecCode = comname,              # common name
     LngtClass = length,              # length bin
     Number = numlen_adj,             # CPUE in n/effort
     LWa = a,                         # length/weight param
     LWb = b,                         # length/weight param
     #bodyMass = ind_weight_g,             # weight of an individual from that size class
     #CPUE_bio_per_hour = freq_weight  # CPUE in biomass
   ) %>% arrange(Year, season, SpecCode, LngtClass) %>% 
   mutate(
     # calculate total biomass, make a key for the length weight coefficient sources
     bodyMass = ind_weight_kg * 1000,            # weight of an individual from that size class
     #CPUE_bio_per_hour = sum_weight_kg * 1000,  # unnecessary
     Biomass = Number * bodyMass,
     lw_group = str_c(SpecCode, season, catchsex)
   )
 

 
 
 
 
 
 
 ####__  SS Species Size Bins  ####
 
 # Split for every species  still??
 # at this point you want to break the data into the groups 
 # that reflect their unique length weight coefficients
 species_splits <- data %>% 
   split(.$SpecCode)
 
 
 # Get a key for each species and size class, with the relevant season and sex
 # assumes 1cm bins for all species as written
 data_bin_key <- species_splits %>% 
   map_dfr(function(species_split){
     
     # pull the distinct length bins
     # and the related lw_groups
     species_split <- species_split %>% 
       distinct(LngtClass, lw_group, .keep_all = T) %>% 
       arrange(LngtClass)
     
     
     # Add the max length for the bin, and its weight
     binned_df <- species_split %>% 
       mutate(
         LngtMax = LngtClass + 1, 
         ln_wmax = (ln_a + LWb * log(LngtMax)),
         wmax    = exp(ln_wmax) * 1000)  %>%
       select(lw_group, 
              SpecCode, 
              season, 
              catchsex, 
              LWa, 
              ln_a, 
              LWb, 
              LngtMin = LngtClass, 
              wmin = bodyMass, 
              LngtMax, 
              wmax,
              -c(Number, Biomass))
     
     # return the clean data
     return(binned_df)
})
 

 
 
 
 # Add the bins back into the original and clean up
 dataBin <- data %>% 
   select(Year, season, area, catchsex, lw_group, SpecCode, spec_class, LngtMin = LngtClass, Number, Biomass) %>% 
   left_join(data_bin_key, by = c("SpecCode", "lw_group", "season", "catchsex", "LngtMin")) 

 
 
 # # Export DataBin
 # write_csv(dataBin, here("data/NEFSC/nefsc_databin_allsizes.csv"))
 
 
 ####__  Set Bodymass Cutoff and Groups  ####
 
 # Set bodymass lower limit
 # Filter for lower end of gear selectivity
 mass_cutoff <- 5 #grams
 dataBin <- filter(dataBin, wmin >= mass_cutoff) %>% 
   filter(season %in% c("SPRING", "FALL")) %>% 
   mutate(season = forcats::fct_drop(season))
 

 
 
 ####_____________________________________________________#### 
 ####  Size Spectra Group Calculations  ####
 
 #####__ 1.  All years, every region  ####
 g1 <- dataBin %>% 
   mutate(group_level = "all_data") %>% 
   split(.$group_level) 
 
 # get SS results
 g1_res <- g1 %>% 
   imap_dfr(group_mle_calc) %>% 
   mutate(stdErr = (abs(confMin - b) + abs(confMax - b)) / (2 * 1.96),
          Year = "all",
          season = "all",
          area = "all",
          C = (b != -1 ) * (b + 1) / ( xmax^(b + 1) - xmin^(b + 1) ) + (b == -1) * 1 / ( log(xmax) - log(xmin)))
 
 
 # Better plot function for group comparisons

 


# plot group comparisons
group_mle_plot(g1_res)



#####__ 2. All Years, each season  ####

# get SS results
g2_res <- dataBin  %>% 
  mutate(group_level = season) %>% 
  split(.$group_level) %>% 
  imap_dfr(group_mle_calc) %>% 
  mutate(stdErr = (abs(confMin - b) + abs(confMax - b)) / (2 * 1.96),
         Year = "all",
         season = group_var,
         area = "all",
         C = (b != -1 ) * (b + 1) / ( xmax^(b + 1) - xmin^(b + 1) ) + (b == -1) * 1 / ( log(xmax) - log(xmin)))

 

# plot group comparisons
group_mle_plot(g2_res)


 

#####__ 3. All Years, regions  ####
 
# get SS results
g3_res <- dataBin  %>% 
  mutate(group_level = area) %>% 
  split(.$group_level) %>% 
  imap_dfr(group_mle_calc, vecDiff = 2) %>% 
  mutate(stdErr = (abs(confMin - b) + abs(confMax - b)) / (2 * 1.96),
         Year = "all",
         season = "all",
         area = group_var,
         C = (b != -1 ) * (b + 1) / ( xmax^(b + 1) - xmin^(b + 1) ) + (b == -1) * 1 / ( log(xmax) - log(xmin)))



# plot group comparisons
group_mle_plot(g3_res)




#####__ 3. All Years, seasons * regions  ####

# get SS results
g4_res <- dataBin  %>% 
  mutate(group_level = str_c(season, area)) %>% 
  split(.$group_level) %>% 
  imap_dfr(.f = group_mle_calc) %>% 
  mutate(stdErr = (abs(confMin - b) + abs(confMax - b)) / (2 * 1.96),
         Year = "all",
         season = case_when(
           str_detect(group_var, "FALL") ~ "FALL",
           str_detect(group_var, "SPRING") ~ "SPRING"),
         area = case_when(
           str_detect(group_var, "GoM") ~ "GoM",
           str_detect(group_var, "SNE") ~ "SNE",
           str_detect(group_var, "MAB") ~ "MAB",
           str_detect(group_var, "GB") ~ "GB"),
         C = (b != -1 ) * (b + 1) / ( xmax^(b + 1) - xmin^(b + 1) ) + (b == -1) * 1 / ( log(xmax) - log(xmin)))


# plot group comparisons
group_mle_plot(g4_res)

 
 #####__ 4. Every year, entire survey  ####

# get SS results
g5_res <- dataBin  %>% 
  mutate(group_level = Year) %>% 
  split(.$group_level) %>% 
  imap_dfr(.f = group_mle_calc) %>% 
  mutate(stdErr = (abs(confMin - b) + abs(confMax - b)) / (2 * 1.96),
         Year = group_var,
         season = "all",
         area = "all",
         C = (b != -1 ) * (b + 1) / ( xmax^(b + 1) - xmin^(b + 1) ) + (b == -1) * 1 / ( log(xmax) - log(xmin)))



# plot group comparisons
group_mle_plot(g5_res)
 
 #####__ 5. every year, every region  ####
 
# get SS results
g6_res <- dataBin  %>% 
  mutate(group_level = str_c(Year, area)) %>% 
  split(.$group_level) %>% 
  imap_dfr(.f = group_mle_calc) %>% 
  mutate(stdErr = (abs(confMin - b) + abs(confMax - b)) / (2 * 1.96),
         Year = str_sub(group_var, 1, 4),
         season = "all",
         area = str_sub(group_var, 5, -1),
         C = (b != -1 ) * (b + 1) / ( xmax^(b + 1) - xmin^(b + 1) ) + (b == -1) * 1 / ( log(xmax) - log(xmin)))
 

# plot group comparisons
group_mle_plot(g6_res)




 #####__ 6. every year, only seasons  ####

# get SS results
g7_res <- dataBin  %>% 
  mutate(group_level = str_c(Year, season)) %>% 
  split(.$group_level) %>% 
  imap_dfr(.f = group_mle_calc, vecDiff = 2) %>% 
  mutate(stdErr = (abs(confMin - b) + abs(confMax - b)) / (2 * 1.96),
         Year = str_sub(group_var, 1, 4),
         season = str_sub(group_var, 5, -1),
         area = "all",
         C = (b != -1 ) * (b + 1) / ( xmax^(b + 1) - xmin^(b + 1) ) + (b == -1) * 1 / ( log(xmax) - log(xmin)))


# plot group comparisons
group_mle_plot(g7_res)
 



 #####__ 7. every year, region * season  ####


# get SS results
g8_res <- dataBin  %>% 
  mutate(group_level = str_c(Year, season, area)) %>% 
  split(.$group_level) %>% 
  imap_dfr(.f = group_mle_calc, vecDiff = 2) %>% 
  mutate(stdErr = (abs(confMin - b) + abs(confMax - b)) / (2 * 1.96),
         Year = str_sub(group_var, 1, 4),
         season = case_when(
           str_detect(group_var, "FALL") ~ "FALL",
           str_detect(group_var, "SPRING") ~ "SPRING"),
         area = case_when(
           str_detect(group_var, "GoM") ~ "GoM",
           str_detect(group_var, "SNE") ~ "SNE",
           str_detect(group_var, "MAB") ~ "MAB",
           str_detect(group_var, "GB") ~ "GB"),
         C = (b != -1 ) * (b + 1) / ( xmax^(b + 1) - xmin^(b + 1) ) + (b == -1) * 1 / ( log(xmax) - log(xmin)))


# plot group comparisons
group_mle_plot(g8_res) +
  facet_wrap(~ area) +
  theme(axis.text.x = element_text(angle = 90, size = 6, vjust = 0.5)) +
  xlab("")




# super table
table_complete <- bind_rows(
  list(
    "Overall"                        = g1_res,
    "only seasons"                   = g2_res,
    "only regions"                   = g3_res,
    "region * seasons"               = g4_res,
    "single years"                   = g5_res,
    "single years * region"          = g6_res,
    "single years * season "         = g7_res,
    "single years * season * region" = g8_res
), .id = "group ID")



# Export table
write_csv(table_complete, path = here("data/size_spectra_results/nefsc_5grams.csv"))

# super plot
table_complete %>% 
  mutate(Year = as.numeric(as.character(Year))) %>% 
  group_mle_plot() + 
  facet_grid(season ~ area) + 
  geom_smooth(method = "gam", show.legend = FALSE, se = FALSE, 
              formula = y ~ s(x, bs = "cs")) +
  geom_hline(data = filter(table_complete, Year == "all"), 
             aes(yintercept = b), linetype = 2, alpha = 0.7, size = 1) +
  theme(axis.text.x = element_text(angle = 90, size= 6, vjust = 0.5)) +
  labs(x = "") + 
  guides(shape = "none",
         color = "none")




table_complete %>% 
  mutate(Year = as.numeric(as.character(Year))) %>% 
  group_mle_plot() + 
  facet_grid(season ~ area) + 
  geom_smooth(method = "gam", show.legend = FALSE, se = FALSE, 
              formula = y ~ s(x, bs = "cs")) +
  geom_hline(data = filter(table_complete, Year == "all"), 
             aes(yintercept = b), linetype = 2, alpha = 0.7, size = 1) +
  theme(axis.text.x = element_text(angle = 90, size= 6, vjust = 0.5)) +
  labs(x = "") + 
  guides(shape = "none",
         color = "none")












####_____________________________________________________#### 

####  Vessel Change  ####
dataBin <- data %>% 
  select(Year, svvessel, season, area, catchsex, lw_group, SpecCode, spec_class, LngtMin = LngtClass, Number, Biomass) %>% 
  left_join(data_bin_key, by = c("SpecCode", "lw_group", "season", "catchsex", "LngtMin")) 



# # Export DataBin
# write_csv(dataBin, here("data/NEFSC/nefsc_databin_allsizes.csv"))


####__  Set Bodymass Cutoff and Groups  ####

# Set bodymass lower limit
# Filter for lower end of gear selectivity
mass_cutoff <- 5 #grams
dataBin <- filter(dataBin, wmin >= mass_cutoff) %>% 
  filter(season %in% c("SPRING", "FALL")) %>% 
  mutate(season = forcats::fct_drop(season))



# get SS results - vessel, season, area
ves_res_1 <- dataBin  %>% 
  mutate(group_level = str_c(svvessel, area, season)) %>% 
  split(.$group_level) %>% 
  imap_dfr(.f = group_mle_calc, vecDiff = 2) %>% 
  mutate(stdErr = (abs(confMin - b) + abs(confMax - b)) / (2 * 1.96),
         Year = "all",
         season = case_when(
           str_detect(group_var, "FALL") ~ "FALL",
           str_detect(group_var, "SPRING") ~ "SPRING"),
         area = case_when(
           str_detect(group_var, "GoM") ~ "GoM",
           str_detect(group_var, "SNE") ~ "SNE",
           str_detect(group_var, "MAB") ~ "MAB",
           str_detect(group_var, "GB") ~ "GB"),
         svvessel = str_sub(group_var, 1, 2),
         C = (b != -1 ) * (b + 1) / ( xmax^(b + 1) - xmin^(b + 1) ) + (b == -1) * 1 / ( log(xmax) - log(xmin)))


# plot group comparisons
group_mle_plot(ves_res_1) +
  facet_grid(svvessel ~ area) +
  theme(axis.text.x = element_text(angle = 90, size = 6, vjust = 0.5)) +
  xlab("")



# get SS results - year, vessel, season, area
ves_res_2 <- dataBin  %>% 
  mutate(group_level = str_c(Year, season, area, svvessel)) %>% 
  split(.$group_level) %>% 
  imap_dfr(.f = group_mle_calc, vecDiff = 2) %>% 
  mutate(stdErr = (abs(confMin - b) + abs(confMax - b)) / (2 * 1.96),
         Year = str_sub(group_var, 1, 4),
         season = case_when(
           str_detect(group_var, "FALL") ~ "FALL",
           str_detect(group_var, "SPRING") ~ "SPRING"),
         area = case_when(
           str_detect(group_var, "GoM") ~ "GoM",
           str_detect(group_var, "SNE") ~ "SNE",
           str_detect(group_var, "MAB") ~ "MAB",
           str_detect(group_var, "GB") ~ "GB"),
         svvessel = str_sub(group_var, 5, 6),
         C = (b != -1 ) * (b + 1) / ( xmax^(b + 1) - xmin^(b + 1) ) + (b == -1) * 1 / ( log(xmax) - log(xmin)))


# plot group comparisons
group_mle_plot(ves_res_2) +
  facet_grid(vessel ~ area) +
  theme(axis.text.x = element_text(angle = 90, size = 6, vjust = 0.5)) +
  xlab("")









