# Date 3/21/2023
# Subject: Stratified Abundance Changes, decadal variability paper:
# Premise: McCalls basin hypothesis suggests that species with shrinking/growing
# populations will differently track a favorable environment
# Goal is to see which species fall into either category
# Range of specific interest is 2000-2010-2019

# Packages
library(gmRi)
library(targets)
library(tidyverse)
library(patchwork)


# Load the clean survdat data with targets
tar_load("survdat_clean")



#----------  Stratum Area Details. ----------

# Need the number of all tows in each stratum, the area of each stratum:

####  1. Import supplemental files  ####
nmfs_path <- cs_path(box_group = "RES_Data", subfolder = "NMFS_trawl")

# Stratum Area Information File
stratum_area_path <- stringr::str_c(nmfs_path, "Metadata/strata_areas_km2.csv")
stratum_area      <- readr::read_csv(stratum_area_path, col_types = readr::cols())
stratum_area      <- dplyr::mutate(stratum_area, stratum = as.character(stratum))



####  2.  Set Constants:  ####

# Area covered by an albatross standard tow in km2
alb_tow_km2 <- 0.0384

# catchability coefficient - ideally should change for species guilds or functional groups. 
q <- 1


####  3. Stratum Area & Effort Ratios  ####
#  Get Annual Stratum Effort, and Area Ratios
# number of tows in each stratum by year
# area of a stratum relative to total area of all stratum sampled that year


# Merge in the area of strata in km2
# (excludes ones we do not care about via left join)
survdat_areas <- dplyr::left_join(survdat_clean, stratum_area, by = "stratum") %>% 
  arrange(id)


# Get Total area of all strata sampled in each year:
total_stratum_areas <- dplyr::group_by(survdat_areas, est_year) %>% 
  distinct(stratum, .keep_all = T) %>% 
  summarise(tot_s_area =  sum(s_area_km2, na.rm = T), .groups = "drop")


# Calculate individual strata area relative to total area that year
# i.e. stratio or stratum weights
survdat_areas <- dplyr::left_join(survdat_areas, total_stratum_areas, by = "est_year")
survdat_areas <- dplyr::mutate(survdat_areas, st_ratio = s_area_km2 / tot_s_area)


# We have total areas, now we want effort within each
# Number of unique tows per stratum, within each season
yr_strat_effort <- dplyr::group_by(survdat_areas, est_year, season, stratum) %>% 
  summarise(strat_ntows = dplyr::n_distinct(id), .groups = "drop")



# Add those yearly effort counts back for later
# (area stratified abundance)
survdat_areas <- dplyr::left_join(survdat_areas, yr_strat_effort, 
                                by = c("est_year", "season", "stratum"))





#--------- Getting Species Specific Abundance Totals --------

# Abundance and biomass across all sizes is measured once for each species
# Need tos trip out the repeated rows to avoid inflating the numbers

species_station_totals <- survdat_areas %>% 
  distinct(id, est_year, season, stratum, comname, abundance, biomass_kg, strat_ntows, st_ratio, tot_s_area)







#--------- Strata specific Abundance Densities. -------------------


# a. Catch / tow, for that year & season
stratified_abundance <-  species_station_totals %>% 
  dplyr::mutate(strata_abundance_cpue = abundance / strat_ntows)


# b. Stratified Mean Catch Rates
# Stratified mean abundance CPUE, weighted by the stratum areas
stratified_abundance <-  stratified_abundance %>% 
  dplyr::mutate(strat_mean_abund_s = strata_abundance_cpue * st_ratio)


# c. Stratified Totals
# convert from catch rate by area swept to total catch for entire stratum
# Depends on catchability (q) at this step, and the albatross area-towed
stratified_abundance <-  stratified_abundance %>% 
  dplyr::mutate(
    # Total Abundance
    strat_total_abund_s = round((strat_mean_abund_s * tot_s_area / alb_tow_km2) / q))






# ---------- Species Filtering. -------------


# Pick the different species out that we are using

species <- c(
  "acadian redfish",
  "alewife",
  "american lobster",
  "american plaice",
  "american shad",
  "atlantic cod",
  "atlantic hagfish",
  "atlantic herring",
  "atlantic mackerel",
  "atlantic rock crab",
  "black sea bass",
  "blackbelly rosefish",
  "butterfish",
  "chain dogfish",
  "cusk",
  "fawn cusk-eel",
  "fourbeard rockling",
  "fourspot flounder",
  "goosefish",
  "gulf stream flounder",
  "haddock",
  "jonah crab",
  "little skate",
  "longfin squid",
  "longhorn sculpin",
  "northern sand lance",
  "northern searobin",
  "northern shortfin squid",
  "ocean pout",
  "offshore hake",
  "pollock",
  "red hake",
  "rosette skate",
  "scup",
  "sea raven",
  "sea scallop",
  "silver hake",
  "smooth dogfish",
  "smooth skate",
  "spiny dogfish",
  "spotted hake",
  "summer flounder",
  "thorny skate",
  "white hake",
  "windowpane",
  "winter flounder",
  "winter skate",
  "witch flounder",
  "yellowtail flounder"
) %>% sort()


# Filter to just those species:
decadal_dat <- stratified_abundance %>% 
  filter(comname %in% tolower(species))







# ---------- Annual Aggregations. ------------


# Do we take the mean across area-weighted catch rates? 
# or the total abundances estimated for the entire area?
annual_abundance_summary <- decadal_dat %>% 
  group_by(year = est_year, comname) %>% 
  summarise(
    estimated_abundance = sum(strat_total_abund_s, na.rm = T),
    area_wtd_density = mean(strat_mean_abund_s, na.rm = T),
    .groups = "drop")



# Catch Density:
annual_abundance_summary %>% 
  filter(comname == "atlantic cod") %>% 
  ggplot(aes(year, area_wtd_density)) +
  geom_line() +
  geom_point() +
  theme_gmri() +
  labs(y = "Mean area-weighted CPUE", title = "Atlantic Cod")

# Stratified Total Abundance Estimate
annual_abundance_summary %>% 
  filter(comname == "atlantic cod") %>% 
  ggplot(aes(year, estimated_abundance)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(y = "Total Area-Stratified Abundance", title = "Atlantic Cod")




# Catch Density:
annual_abundance_summary %>% 
  filter(comname == "atlantic cod") %>% 
  ggplot(aes(year, area_wtd_density)) +
  geom_line() +
  geom_point() +
  theme_gmri() +
  labs(y = "Mean area-weighted CPUE", title = "Atlantic Cod")



#------------  Design the Plot --------

# Do a trendline of the last 20 years
# Color based on direction and significance

plot_strat_abund <- function(species_x){
  one_species <- annual_abundance_summary %>% filter(comname == tolower(species_x)) %>% 
    mutate(comname = toupper(comname), 
           abund_mill = estimated_abundance/1e6)
  # lm_coef <- lm(estimated_abundance ~ year, data = one_species) %>% coef() %>% as.numeric()
  # direction_col <- ifelse(lm_coef[2] > 0, gmri_cols("blue"), gmri_cols("orange"))
  ggplot(one_species, aes(year, abund_mill)) +
    geom_line(color = gmri_cols("light gray"), linewidth = 0.5) +
    geom_point(color = "black", size = 0.5) +
    #geom_smooth(method = "lm", color = direction_col, formula = y~x, se = F, linewidth = 1) +
    #facet_wrap(~comname) +
    theme_gmri() +
    #scale_y_continuous(labels = scales::comma_format()) +
  #scale_y_continuous(labels = scales::number_format(suffix = "M")) +
    scale_y_continuous(labels = scales::label_comma(suffix = "M")) +
    labs(y = "Abundance", x = "Year", title = toupper(species_x))
  
  
}




# Run it for all of them:
abundance_figs <- map(species, plot_strat_abund) %>% 
  setNames(species)





# ----------- Saving -------------





# Path to decadal folder:
decadal_folder <- cs_path(box_group = "mills", subfolder = "Projects/Decadal Variability/Graphics/Strat_Abund")


# Pages of 24:
species_abund_1 <- wrap_plots(abundance_figs[1:24], ncol = 4, nrow = 6, widths = 3, heights = 3)
species_abund_2 <- wrap_plots(abundance_figs[25:48], ncol = 4, nrow = 6, widths = 3, heights = 3)
species_abund_3 <- wrap_plots(abundance_figs[49], ncol = 4, nrow = 6, widths = 3, heights = 3)


ggsave(str_c(decadal_folder, "strat_abundance_p1.pdf"), species_abund_1, height = 15, width = 12.5, units ="in")
ggsave(str_c(decadal_folder, "strat_abundance_p2.pdf"), species_abund_2, height = 15, width = 12.5, units ="in")
ggsave(str_c(decadal_folder, "strat_abundance_p3.pdf"), species_abund_3, height = 15, width = 12.5, units ="in")


# Pages of 9 (3x3), 5 pages



# # 3 x 3
# # Using patchwork?
# species_abund_1 <- wrap_plots(abundance_figs[1:24], ncol = 4, nrow = 6, widths = 3, heights = 3)
# species_abund_2 <- wrap_plots(abundance_figs[10:18], ncol = 3, nrow = 3, widths = 3, heights = 3)
# species_abund_3 <- wrap_plots(abundance_figs[19:27], ncol = 3, nrow = 3, widths = 3, heights = 3)
# species_abund_4 <- wrap_plots(abundance_figs[28:36], ncol = 3, nrow = 3, widths = 3, heights = 3)
# species_abund_5 <- wrap_plots(abundance_figs[37:45], ncol = 3, nrow = 3, widths = 3, heights = 3)
# 
# # Save pages
# ggsave(str_c(decadal_folder, "strat_abundance_p1.pdf"), species_abund_1, height = 15, width = 12.5, units ="in")
# ggsave(str_c(decadal_folder, "strat_abundance_p2.pdf"), species_abund_2, height = 15, width = 12.5, units ="in")
# ggsave(str_c(decadal_folder, "strat_abundance_p3.pdf"), species_abund_3, height = 15, width = 12.5, units ="in")
# ggsave(str_c(decadal_folder, "strat_abundance_p4.pdf"), species_abund_4, height = 15, width = 12.5, units ="in")
# ggsave(str_c(decadal_folder, "strat_abundance_p5.pdf"), species_abund_5, height = 15, width = 12.5, units ="in")






# # Group the data based on pre and post 2010
# one_species %>% 
#   mutate(comparison_period = ifelse(year <2010, "period_1", "period_2")) %>% 
#   t.test(estimated_abundance ~ comparison_period, data = .) %>% 
#   broom::tidy() %>% 
#   mutate(outcome = case_when(
#     p.value <= 0.05 & estimate1 > estimate2 ~ "pop. decrease",
#     p.value <= 0.05 & estimate1 < estimate2 ~ "pop. increase",
#     p.value > 0.05 & estimate1 > estimate2 ~ "no change"),
#     perc_change = ((estimate2-estimate1)/estimate1)*100) %>% 
#   select(estimate1, estimate2, p.value, outcome, perc_change)
# 
# 
# # Stratified Total Abundance Estimate
# annual_abundance_summary %>% 
#   filter(comname == species_x) %>% 
#   ggplot(aes(year, area_wtd_density)) +
#   geom_line() +
#   geom_point() +
#   scale_y_continuous(labels = scales::comma_format()) +
#   labs(y = "Estimated Area-Stratified Abundance",
#        title = species_x)



