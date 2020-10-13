####

# Goal: attribute ocean region to nefsc trawl size spectra
# Do areas of inherently different temperature reflect differences in slope/shape
# Is this holding up, does GOM resemble mid atlantic now
# Want to end up with a chronological record of what the slope is for every year/season/region and how that region did for the entire year as well


# Output files:
# dataBin & results from 5g cutoff are saved for use in 06_ss_sensitivity


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
`%notin%` <-  purrr::negate(`%in%`) # to filter out ones that have matches


# File Paths
mills_path <- shared.path(os.use = "unix", group = "Mills Lab", folder = NULL)
res_path   <- shared.path(os.use = "unix", group = "RES Data", folder = NULL)



####_____________________________________________________####
####  X MOVED to 01_nefsc_ss_build.R ####
# #### X  Data  ####
# 
# #####X __ a. NEFSC bottom trawl  ####
# nefsc_master <- read_csv(here("data/NEFSC/2020survdat_nye.csv"), 
#                         col_types = cols(), 
#                         guess_max = 1e6) %>% 
#   clean_names()
# 
# # Which months do the seasons cover
# nefsc_master %>% filter(season %in% c("SPRING", "FALL")) %>%  count(season, est_month)
# nefsc_master %>% filter(season %in% c("SPRING", "FALL")) %>%  count(est_year)
# 
# 
# # Drop unhelpful columnns to free up space
# 
# # Cleanup Code from Kathy
# nefsc_master <- nefsc_master %>% 
#   mutate(comname = tolower(comname),
#          id = format(id, scientific = FALSE)) %>%
#   select(id, est_year, season, stratum, decdeg_beglat, decdeg_beglon,
#          svspp, comname, catchsex, biomass, avgdepth, abundance, length, numlen) %>%
#   filter(stratum >= 01010,
#          stratum <= 01760,
#          stratum != 1310,
#          stratum != 1320,
#          stratum != 1330,
#          stratum != 1350,
#          stratum != 1410,
#          stratum != 1420,
#          stratum != 1490,
#          season %in% c("SPRING", "FALL"),
#          est_year >= 1970,
#          est_year < 2020) %>%
#   mutate(biomass = ifelse(is.na(biomass) == TRUE & abundance > 0, 0.01, biomass),
#          biomass = ifelse(biomass == 0 & abundance > 0, 0.01, biomass),
#          abundance = ifelse(is.na(abundance) == TRUE & biomass > 0, 1, abundance),
#          abundance = ifelse(abundance == 0 & biomass > 0, 1, abundance),
#          tstrat = str_sub(stratum, 2, 3)) %>%  
#   filter(!is.na(biomass),
#          !is.na(abundance))
# 
# 
# ####X __ b. L/W coefficients  ####
# 
# # Load the lw key with coefficients from fishbase
# nefsc_lw <- read_csv(here("data/NEFSC/nefsc_lw_key_filled.csv"),
#                      guess_max = 1e3,
#                      col_types = cols())
# 
# # Fill in NA growth coefficients using related species
# nefsc_lw <- nefsc_lw %>% 
#   mutate(a = as.numeric(a),
#          b = as.numeric(b),
#          ln_a = log(a),
#          related_ln_a = log(related_a),
#          ln_a = ifelse(is.na(ln_a), related_ln_a, ln_a),
#          a = ifelse(is.na(a), exp(ln_a), a))
# 
# 
# 
# # Wrigley Paper, length weight coefficients
# load("/Users/akemberling/Box/RES Data/NMFS_trawl/lwreg.Rdata")
# lwreg <- lwreg %>% clean_names()  
# lwreg <- lwreg %>%  
#   mutate(comname = str_to_lower(common_name), .before = scientific_name,
#          common_name = NULL,
#          scientific_name = str_to_lower(scientific_name),
#          svspp = str_pad(svspp, width = 3, side = "left", pad = "0")) %>% 
#   mutate(ln_a = lna, 
#          a = exp(ln_a),
#          lna = NULL,
#          .before = b)
# 
# 
# 
# # Species class/groups based on life history
# spp_classes <- read_csv(here("data/kmills/sppclass.csv"),
#                         col_types = cols())
# spp_classes <- spp_classes %>% 
#   clean_names() %>% 
#   mutate(common_name = str_to_lower(common_name),
#          scientific_name = str_to_lower(scientific_name)) %>% 
#   pivot_longer(gf:inv, names_to = "spec_class", values_to = "flag") %>% 
#   pivot_longer(com:nc, names_to = "fishery", values_to = "flag2") %>% 
#   filter(is.na(flag) == FALSE,
#          is.na(flag2) == FALSE) %>% 
#   select(svspp, comname = common_name, everything(), -flag, -flag2)
#  
# 
# # Merge the growth coefficients with the species classes
# wigley_lw <- left_join(lwreg, spp_classes, by = c("svspp", "comname", "scientific_name"))
# 
# 
# 
# ####X __ c.  Combining L/W  Sources  ####
# 
# # 1. use these growth coefficients to get biomasses
# # 2. Potentially merge in my fishbase ones so that there is just one master
# 
# wigley_lw <- wigley_lw %>% mutate(source = "wigley", .before = svspp)
# nefsc_lw  <- nefsc_lw %>% mutate(source = "fishbase", .before = comname) 
# lw_combined <- full_join(nefsc_lw, wigley_lw) 
# 
# # Fill in gaps for things that should be consistent
# lw_combined <- lw_combined %>% 
#   arrange(comname) %>% 
#   mutate(
#     hare_group = ifelse(is.na(hare_group),
#                              nefsc_lw[ match(comname,  nefsc_lw$comname)[1], "hare_group"][[1]],
#                              hare_group),
#          svspp = ifelse(is.na(svspp),
#                               wigley_lw[ match(comname, wigley_lw$comname)[1], "svspp"][[1]],
#                               svspp),
#          scientific_name = ifelse(is.na(scientific_name),
#                                   wigley_lw[ match(comname, wigley_lw$comname)[1], "scientific_name"][[1]],
#                                   scientific_name),
#          fishery = ifelse(is.na(fishery),
#                           wigley_lw[ match(comname, wigley_lw$comname)[1], "fishery"][[1]],
#                           fishery),
#          spec_class = ifelse(is.na(spec_class),
#                              wigley_lw[ match(comname, wigley_lw$comname)[1], "spec_class"][[1]],
#                              spec_class)) %>% 
#   select(source, svspp, comname, scientific_name, spec_class, hare_group, fishery, season, catchsex,
#          b, a, ln_a, units, related_species, related_a, related_b, related_ln_a,
#          wigleymin, wigleymax, wigleyvalue)
#   
# 
# 
# # 3. From there we can do the rest of the size spectra steps
# # continue to get individual length bins
# 
# 
# 
# 
# 
# 
# 
# 
# ####  X Data  Prep  ####
# 
# 
# #####X __ 1.  Stratum-Area Key  ####
# 
# 
# # If we know which strata go where we can use this:
# # Stratum Key for filtering specific areas
# strata_key <- list(
#   "Georges Bank"          = as.character(13:23),
#   "Gulf of Maine"         = as.character(24:40),
#   "Southern New England"  = as.character(1:10),
#   "Mid-Atlantic Bight"    = as.character(61:76))
# 
# # Use strata_select to pull the strata we want
# strata_select <- c(
#   strata_key$`Georges Bank`, 
#   strata_key$`Gulf of Maine`,
#   strata_key$`Southern New England`,
#   strata_key$`Mid-Atlantic Bight`)
# 
# 
# 
# 
# 
# ####X __ 2. Temporal Filtering  ####
# #Filtering
# nefsc <- nefsc_master %>% 
#   mutate(area =  case_when(
#            tstrat %in% strata_key$`Georges Bank`         ~ "GB", 
#            tstrat %in% strata_key$`Gulf of Maine`        ~ "GoM",
#            tstrat %in% strata_key$`Southern New England` ~ "SNE",
#            tstrat %in% strata_key$`Mid-Atlantic Bight`   ~ "MAB",
#            TRUE                                          ~ "not found")) %>% 
#   filter(#est_year %in% c(1982:2019),
#          season   %in% c("SPRING", "FALL"),
#          tstrat   %in% strata_select)
# 
# 
# # Check what is going on with numlen and abundance, make conversion factor
# conv_factor <- nefsc %>% 
#   group_by(id, comname, catchsex, numlen, abundance) %>% 
#   summarise(abundance_raw = sum(numlen),
#             abundance = mean(abundance), 
#             convers =  abundance/abundance_raw) 
# 
# 
# # Merge back and convert the numlen field
# nefsc <- nefsc %>% 
#   left_join(conv_factor, by = c("id", "comname", "catchsex", "numlen", "abundance")) %>% 
#   mutate(numlen_adj = numlen * convers,
#          numlen_adj = round(numlen_adj, 2)) 
# 
# 
# 
# glimpse(nefsc)
# 
# 
# 
# 
# 
# ####X __ 3. Get Individual Lengths  ####
# 
# # Record of unique station catches: # rows for each species * sex * length
# nefsc_lens <- nefsc %>% 
#   filter(is.na(numlen) == FALSE,
#          numlen > 0) %>% 
#   distinct(id, comname, catchsex, length, numlen_adj) 
# 
# 
# # recombine with the distinct station info
# nefsc_spectra <- nefsc %>% 
#   # drop the columns with individual catch info, we want all distinct records merging in ok
#   select(-c(catchsex, svspp, avgdepth, comname, numlen_adj, biomass, abundance, length, numlen, abundance_raw, convers)) %>% 
#   distinct() %>% 
#   left_join(nefsc_lens, by = "id")
# 
# 
# glimpse(nefsc_spectra)
# 
# 
# 
# ####X __ 4.  Calculate Weights/Biomass  ####
# 
# 
# 
# #Now we want to use the lw_combined here instead of just the fishbase lengths
# 
# # Do a priority pass with the filter(lw_combined, source == "wigley)
# # merge on comname, season, and catchsex
# w_trimmed <- filter(lw_combined, source == "wigley") %>% 
#   select(source, season, comname, scientific_name, spec_class, hare_group, catchsex, a, b, ln_a)
# 
# 
# # Do a second pass with the filter(lw_combined, source == "fishbase")
# # merge on common names only
# fb_trimmed <- filter(lw_combined, source == "fishbase") %>% 
#   select(source, comname, scientific_name, spec_class, hare_group, a, b, ln_a)
# 
# # First Pass - Wigley
# pass_1 <- nefsc_spectra %>% 
#   inner_join(w_trimmed)
# 
# 
# # Second Pass - Fishbase, for the stragglers if any
# pass_2 <- nefsc_spectra %>% 
#   filter(comname %notin% w_trimmed$comname) %>% 
#   inner_join(fb_trimmed)
# 
# 
# 
# # Join them with bind rows (implicitly drops things that don't have growth coefs)
# nefsc_weights <- bind_rows(pass_1, pass_2) %>% 
#   arrange(est_year, season) %>% 
#   mutate(b = as.numeric(b),
#          a = as.numeric(a),
#          a = ifelse(is.na(a) & !is.na(ln_a), exp(ln_a), a),
#          ln_a = ifelse(is.na(ln_a), log(a), ln_a),
#          ln_weight = (ln_a + b * log(length)),
#          weight_g = exp(ln_weight),
#          freq_weight = weight_g * numlen_adj) %>% 
#   drop_na(weight_g)
# 
# 
# # Export the key you are using
# # write_csv(nefsc_weights, here("data/size_spectra_results/lw_kew_combined.csv"))


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
     bodyMass = weight_g,             # weight of an individual from that size class
     CPUE_bio_per_hour = freq_weight  # CPUE in biomass
   ) %>% arrange(Year, season, SpecCode, LngtClass)
 
 # calculate total biomass, make a key for the length weight coefficient sources
 data <- data %>% 
   mutate(Biomass = Number * bodyMass,
          lw_group = str_c(SpecCode, season, catchsex))
 
 
 
 
 
 
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
         wmax    = exp(ln_wmax))  %>%
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
group_mle_plot <- function(mle_res){
plot <- mle_res %>% 
   ggplot(aes(Year, b, color = area, shape = season)) +
   geom_pointrange(aes(x = Year, y = b, ymin = confMin, ymax = confMax)) +
   labs(x = "Year",
        y = "Size Spectrum Slope (b)") 

return(plot)
}
 


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
