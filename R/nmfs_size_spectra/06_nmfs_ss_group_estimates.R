####

# Goal: attribute ocean region to nefsc trawl size spectra
# Do areas of inherently different temperature reflect differences in slope/shape
# Is this holding up, does GOM resemble mid atlantic now
# Want to end up with a chronological record of what the slope is for every year/season/region and how that region did for the entire year as well
#
# Results from MLE estimation of size spectrum slope is done for each factor group split.
# To produce yearly, seasonally, slope measures across regions

# Output files:
# dataBin & results from 5g cutoff are saved for use in 09_ss_sensitivity


# STATUS 6/8/2020: Moved into targets workflow as a function


####  Packages  ####
library(here)
library(gmRi)
library(sf)
library(janitor)
library(sizeSpectra)
library(patchwork)
library(tidyverse)
library(targets)



####  Support Functions  ####
source(here("R/support/sizeSpectra_support.R"))



# File Paths
mills_path <- shared.path(os.use = "unix", group = "Mills Lab", folder = NULL)
res_path   <- shared.path(os.use = "unix", group = "RES Data", folder = NULL)





####_____________________________________________________####
####  Load Prepped Data  ####

# # source file : 02_NOAA_QAQC.R
# nefsc_weights <- read_csv(here::here("data/ss_prepped_data/survdat_2020_ss.csv"),
#                           col_types = cols(),
#                           guess_max = 1e5)

tar_load("nefsc_stratified")
nefsc_weights <- nefsc_stratified

# Check weights 
nefsc_weights %>% group_by(comname) %>% summarise(mean_weight =  mean(ind_weight_kg)) %>% arrange(desc(mean_weight)) 
nefsc_weights %>% group_by(comname) %>% summarise(mean_weight =  mean(ind_weight_kg)) %>% arrange(mean_weight)


#plot them all
weights_summ <- nefsc_weights %>% 
  group_by(comname) %>% 
  summarise(mean_weight = mean(ind_weight_kg, na.rm = T),
            bin = case_when(
              mean_weight <= .25 ~ "0 - 0.25",
              mean_weight <= 5 ~ "0.25 - 0.5",
              mean_weight <= 1 ~ "0.5 - 1",
              mean_weight <= 5 ~ "1 - 5",
              mean_weight <= 10 ~ "5 - 10",
              mean_weight <= 50 ~ "10 - 50",
              TRUE ~ "50+")) %>% 
  left_join(nefsc_weights, by = "comname") 



weights_summ %>% 
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
   mutate(season = str_to_title(season),
          season = factor(season, levels = c("Spring", "Fall"))) %>%
   rename(
     Year = est_year,                 # year
     SpecCode = comname,              # common name
     LngtClass = length,              # length bin
     Number = numlen_adj,             # CPUE in n/effort
     LWa = a,                         # length/weight param
     LWb = b) %>%                          # length/weight param
   mutate(
     # calculate total biomass, make a key for the length weight coefficient sources
     bodyMass = ind_weight_kg * 1000,            # weight of an individual from that size class
     Biomass = Number * bodyMass,                # number at size times the mass
     Stratified_Biomass = expanded_abund_s * bodyMass, # number scaled to stratum, time the mass of ind.
     lw_group = str_c(SpecCode, season, catchsex)) %>% 
   arrange(Year, season, SpecCode, LngtClass)
 


 # Pretty sure since the growth coefficients,
 # stratified numbers, and number at length are all together we can just use mutate here...
 dataBin <- data %>% 
   mutate( 
     LngtMax         = LngtClass + 1,
     Ln_wmax         = (ln_a + LWb * log(LngtMax)),
     wmin            = bodyMass,                     # minimum weight of individual
     wmax            = exp(Ln_wmax) * 1000,          # max weight of individual
     wmin_sum        = Number * wmin,                # wmin * number caught actual
     wmax_sum        = Number * wmax                 # wmax * number caught actual  
    ) %>% 
   filter(expanded_abund_s != 0)
 
 
 # # Export DataBin
 # write_csv(dataBin, here("data/NEFSC/nefsc_databin_allsizes.csv"))
 # 
 
 
 
 ####__  Set Bodymass Cutoff and Groups  ####
 
 # Set bodymass lower limit
 # Filter for lower end of gear selectivity
 mass_cutoff <- 1 #grams
 dbin_truncated <- filter(dataBin, wmin >= mass_cutoff) 
 

 
 
 ####_____________________________________________________#### 
 ####__  Survey Abundance Group Comparisons  __####
 
 #####__ 1.  All years, every region  ####
 g1 <- dbin_truncated %>% 
   mutate(group_level = "all_data") %>% 
   split(.$group_level) 
 
 
 # get SS results
 g1_res <- g1 %>% 
   imap_dfr(group_mle_calc) %>% 
   mutate(Year = "all",
          season = "all",
          area = "all")
 
 
# plot group comparisons
group_mle_plot(g1_res)



#####__ 2. All Years, each season  ####

# get SS results
g2_res <- dbin_truncated  %>% 
  mutate(group_level = season) %>% 
  split(.$group_level) %>% 
  imap_dfr(group_mle_calc) %>% 
  mutate(Year = "all",
         season = group_var,
         area = "all")

# plot group comparisons
group_mle_plot(g2_res)


 

#####__ 3. All Years, regions  ####
 
# get SS results
g3_res <- dbin_truncated  %>% 
  mutate(group_level = survey_area) %>% 
  split(.$group_level) %>% 
  imap_dfr(group_mle_calc, vecDiff = 2) %>% 
  mutate(Year = "all",
         season = "all",
         area = group_var)



# plot group comparisons
group_mle_plot(g3_res)




#####__ 3. All Years, seasons * regions  ####

# get SS results
g4_res <- dbin_truncated  %>% 
  mutate(group_level = str_c(season, survey_area)) %>% 
  split(.$group_level) %>% 
  imap_dfr(.f = group_mle_calc) %>% 
  mutate(Year = "all",
         season = case_when(
           str_detect(group_var, "Fall") ~ "Fall",
           str_detect(group_var, "Spring") ~ "Spring"),
         area = case_when(
           str_detect(group_var, "GoM") ~ "GoM",
           str_detect(group_var, "SNE") ~ "SNE",
           str_detect(group_var, "MAB") ~ "MAB",
           str_detect(group_var, "GB") ~ "GB"))


# plot group comparisons
group_mle_plot(g4_res)

 
 #####__ 4. Every year, entire survey  ####

# get SS results
g5_res <- dbin_truncated  %>% 
  mutate(group_level = Year) %>% 
  split(.$group_level) %>% 
  imap_dfr(.f = group_mle_calc) %>% 
  mutate(Year = group_var,
         season = "all",
         area = "all")



# plot group comparisons
group_mle_plot(g5_res)
 
 #####__ 5. every year, every region  ####
 
# get SS results
g6_res <- dbin_truncated  %>% 
  mutate(group_level = str_c(Year, survey_area)) %>% 
  split(.$group_level) %>% 
  imap_dfr(.f = group_mle_calc) %>% 
  mutate(Year = str_sub(group_var, 1, 4),
         season = "all",
         area = str_sub(group_var, 5, -1))
 

# plot group comparisons
group_mle_plot(g6_res)




 #####__ 6. every year, only seasons  ####

# get SS results
g7_res <- dbin_truncated  %>% 
  mutate(group_level = str_c(Year, season)) %>% 
  split(.$group_level) %>% 
  imap_dfr(.f = group_mle_calc, vecDiff = 2) %>% 
  mutate(Year = str_sub(group_var, 1, 4),
         season = str_sub(group_var, 5, -1),
         area = "all")


# plot group comparisons
group_mle_plot(g7_res)
 



 #####__ 7. every year, region * season  ####


# get SS results
g8_res <- dbin_truncated  %>% 
  mutate(group_level = str_c(Year, season, survey_area)) %>% 
  split(.$group_level) %>% 
  imap_dfr(.f = group_mle_calc, vecDiff = 2) %>% 
  mutate(Year = str_sub(group_var, 1, 4),
         season = case_when(
           str_detect(group_var, "Fall") ~ "Fall",
           str_detect(group_var, "Spring") ~ "Spring"),
         area = case_when(
           str_detect(group_var, "GoM") ~ "GoM",
           str_detect(group_var, "SNE") ~ "SNE",
           str_detect(group_var, "MAB") ~ "MAB",
           str_detect(group_var, "GB") ~ "GB"))


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
    "single years * season * region" = g8_res), 
  .id = "group ID")




# # plotting all the changes
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


# # Export table
# write_csv(table_complete, path = here("data/size_spectra_results/nefsc_5grams.csv"))
















####_____________________________________________________#### 

####  Vessel Change Tracking  ####

# Set bodymass lower limit
# Filter for lower end of gear selectivity
mass_cutoff <- 5 #grams
dbin_vessels <- filter(dataBin, wmin >= mass_cutoff) %>% 
  filter(season %in% c("Spring", "Fall")) %>% 
  mutate(season = forcats::fct_drop(season))



#####__ 1. Vessel Area Season  ####
ves_res_1 <- dbin_vessels  %>% 
  mutate(group_level = str_c(svvessel, survey_area, season)) %>% 
  split(.$group_level) %>% 
  imap_dfr(.f = group_mle_calc) %>% 
  mutate(Year = "all",
         season = case_when(
           str_detect(group_var, "Fall") ~ "Fall",
           str_detect(group_var, "Spring") ~ "Spring"),
         area = case_when(
           str_detect(group_var, "GoM") ~ "GoM",
           str_detect(group_var, "SNE") ~ "SNE",
           str_detect(group_var, "MAB") ~ "MAB",
           str_detect(group_var, "GB") ~ "GB"),
         svvessel = str_sub(group_var, 1, 2))


# plot group comparisons
group_mle_plot(ves_res_1) +
  facet_grid(area ~ svvessel) +
  theme(axis.text.x = element_text(angle = 90, size = 6, vjust = 0.5)) +
  xlab("")



####__ 2. Year Vessel Season Area  ####
ves_res_2 <-  dbin_vessels  %>% 
  mutate(group_level = str_c(Year, season, survey_area, svvessel)) %>% 
  split(.$group_level) %>% 
  imap_dfr(.f = group_mle_calc) %>% 
  mutate(Year = str_sub(group_var, 1, 4),
         season = case_when(
           str_detect(group_var, "Fall") ~ "Fall",
           str_detect(group_var, "Spring") ~ "Spring"),
         area = case_when(
           str_detect(group_var, "GoM") ~ "GoM",
           str_detect(group_var, "SNE") ~ "SNE",
           str_detect(group_var, "MAB") ~ "MAB",
           str_detect(group_var, "GB") ~ "GB"),
         svvessel = case_when(
           str_detect(group_var, "AL") ~ "AL",
           str_detect(group_var, "HB", ) ~ "HB"))


# plot group comparisons
group_mle_plot(ves_res_2) +
  facet_grid(area ~ svvessel, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, size = 6, vjust = 0.5)) +
  xlab("")

# Different style
ggplot(ves_res_2, aes(y = b, x = svvessel, fill = area)) +
  geom_violin() +
  geom_vline(xintercept = 1.5) +
  theme_bw()

ggplot(ves_res_2, aes(y = C, x = svvessel, fill = area)) +
  geom_violin() +
  geom_vline(xintercept = 1.5) +
  theme_bw()

####__ 3. Only Vessel  ####
ves_res_3 <-  dbin_vessels  %>% 
  mutate(group_level = svvessel) %>% 
  split(.$group_level) %>% 
  imap_dfr(.f = group_mle_calc) %>% 
  mutate(Year = "all",
         season = "all",
         area = "all",
         svvessel = case_when(
           str_detect(group_var, "AL") ~ "AL",
           str_detect(group_var, "HB", ) ~ "HB"))


# plot
group_mle_plot(ves_res_3) +
  facet_wrap(~svvessel)



# Combine the vessel comparison data for export
vessel_results <- bind_rows(
  list(
    "all_years_all_factors"     = ves_res_1,
    "single_years_all_factors"  = ves_res_2,
    "vessel_only"               = ves_res_3),
  .id = "group ID"
)


# # Export table
# write_csv(vessel_results, path = here("data/size_spectra_results/nefsc_5grams_vessel.csv"))





####_____________________________________________________#### 
####_____________________________________________________#### 
####__  Stratified Abundances  __####


#####__ 1.  All years, every region  ####

# get SS results
strat_g1_res <- dbin_truncated %>% 
  mutate(group_level = "all_data") %>% 
  split(.$group_level)  %>% 
  imap_dfr(strat_abund_mle_calc) %>% 
  mutate(Year = "all",
         season = "all",
         area = "all")


# plot group comparisons
group_mle_plot(strat_g1_res)



#####__ 2. All Years, each season  ####

# get SS results
strat_g2_res <- dbin_truncated  %>% 
  mutate(group_level = season) %>% 
  split(.$group_level) %>% 
  imap_dfr(strat_abund_mle_calc) %>% 
  mutate(Year = "all",
         season = group_var,
         area = "all")

# plot group comparisons
group_mle_plot(strat_g2_res)




#####__ 3. All Years, regions  ####

# get SS results
strat_g3_res <- dbin_truncated  %>% 
  mutate(group_level = survey_area) %>% 
  split(.$group_level) %>% 
  imap_dfr(strat_abund_mle_calc, vecDiff = 2) %>% 
  mutate(Year = "all",
         season = "all",
         area = group_var)



# plot group comparisons
group_mle_plot(strat_g3_res)




#####__ 3. All Years, seasons * regions  ####

# get SS results
strat_g4_res <- dbin_truncated  %>% 
  mutate(group_level = str_c(season, survey_area)) %>% 
  split(.$group_level) %>% 
  imap_dfr(.f = strat_abund_mle_calc) %>% 
  mutate(Year = "all",
         season = case_when(
           str_detect(group_var, "Fall") ~ "Fall",
           str_detect(group_var, "Spring") ~ "Spring"),
         area = case_when(
           str_detect(group_var, "GoM") ~ "GoM",
           str_detect(group_var, "SNE") ~ "SNE",
           str_detect(group_var, "MAB") ~ "MAB",
           str_detect(group_var, "GB") ~ "GB"))


# plot group comparisons
group_mle_plot(strat_g4_res)


#####__ 4. Every year, entire survey  ####

# get SS results
strat_g5_res <- dbin_truncated  %>% 
  mutate(group_level = Year) %>% 
  split(.$group_level) %>% 
  imap_dfr(.f = strat_abund_mle_calc) %>% 
  mutate(Year = group_var,
         season = "all",
         area = "all")



# plot group comparisons
group_mle_plot(strat_g5_res)

#####__ 5. every year, every region  ####

# get SS results
strat_g6_res <- dbin_truncated  %>% 
  mutate(group_level = str_c(Year, survey_area)) %>% 
  split(.$group_level) %>% 
  imap_dfr(.f = strat_abund_mle_calc) %>% 
  mutate(Year = str_sub(group_var, 1, 4),
         season = "all",
         area = str_sub(group_var, 5, -1))


# plot group comparisons
group_mle_plot(strat_g6_res)




#####__ 6. every year, only seasons  ####

# get SS results
strat_g7_res <- dbin_truncated  %>% 
  mutate(group_level = str_c(Year, season)) %>% 
  split(.$group_level) %>% 
  imap_dfr(.f = strat_abund_mle_calc, vecDiff = 2) %>% 
  mutate(Year = str_sub(group_var, 1, 4),
         season = str_sub(group_var, 5, -1),
         area = "all")


# plot group comparisons
group_mle_plot(strat_g7_res)




#####__ 7. every year, region * season  ####


# get SS results
strat_g8_res <- dbin_truncated  %>% 
  mutate(group_level = str_c(Year, season, survey_area)) %>% 
  split(.$group_level) %>% 
  imap_dfr(.f = strat_abund_mle_calc, vecDiff = 2) %>% 
  mutate(Year = str_sub(group_var, 1, 4),
         season = case_when(
           str_detect(group_var, "Fall") ~ "Fall",
           str_detect(group_var, "Spring") ~ "Spring"),
         area = case_when(
           str_detect(group_var, "GoM") ~ "GoM",
           str_detect(group_var, "SNE") ~ "SNE",
           str_detect(group_var, "MAB") ~ "MAB",
           str_detect(group_var, "GB") ~ "GB"))


# plot group comparisons
group_mle_plot(strat_g8_res) +
  facet_wrap(~ area) +
  theme(axis.text.x = element_text(angle = 90, size = 6, vjust = 0.5)) +
  xlab("")




# super table
strat_table_complete <- bind_rows(
  list(
    "Overall"                        = strat_g1_res,
    "only seasons"                   = strat_g2_res,
    "only regions"                   = strat_g3_res,
    "region * seasons"               = strat_g4_res,
    "single years"                   = strat_g5_res,
    "single years * region"          = strat_g6_res,
    "single years * season "         = strat_g7_res,
    "single years * season * region" = strat_g8_res), 
  .id = "group ID")




# # plotting all the changes
strat_table_complete %>% 
  mutate(Year = as.numeric(as.character(Year))) %>% 
  group_mle_plot() + 
  facet_grid(season ~ area) + 
  geom_smooth(method = "gam", show.legend = FALSE, se = FALSE, 
              formula = y ~ s(x, bs = "cs")) +
  geom_hline(data = filter(strat_table_complete, Year == "all"), 
             aes(yintercept = b), linetype = 2, alpha = 0.7, size = 1) +
  theme(axis.text.x = element_text(angle = 90, size= 6, vjust = 0.5)) +
  labs(x = "") + 
  guides(shape = "none",
         color = "none")


# # Export table
# write_csv(strat_table_complete, path = here("data/size_spectra_results/nefsc_5grams_strat.csv"))



