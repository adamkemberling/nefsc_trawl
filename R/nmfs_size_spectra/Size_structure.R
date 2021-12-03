

####  Size Structure of Haddock, Cod and Species Guilds
# uses weighted average of biomass / stratified total abundance

####  Load Packages  ####
library(targets)
library(rnaturalearth)
library(gmRi)
library(patchwork)
library(tidyverse)
library(ggforce)
library(ggridges)
library(scales)
library(reldist)

# ggplot theme
theme_set(theme_bw() + theme(plot.title = element_text(hjust = 0.05)))


####  Load Data  ####


# Load the area-stratified biomass/abundances that used species cutoffs
withr::with_dir(rprojroot::find_root('_targets.R'), 
                tar_load(nefsc_stratified))


# Do some text formatting
nefsc_stratified <- nefsc_stratified %>% 
  mutate(
    id = as.character(id),
    spec_class = case_when(
      spec_class == "dem" ~ "Demersal",
      spec_class == "gf" ~ "Groundfish",
      spec_class == "pel" ~ "Pelagic",
      spec_class == "el" ~ "Elasmobranch",
      TRUE ~ "Unclassified"),
    season = factor(season, levels = c("Spring", "Fall")),
    Year = fct_rev(factor(est_year)),
    decade = floor_decade(est_year),
    decade = fct_rev(decade))


####______________________####


####  Plot Functions  ####


# Plot stratified abundance and changes in size
plot_species_dist <- function(comname_choice){
  
  species_dat <- filter(nefsc_stratified, comname == comname_choice)
  
  # Get yearly totals
  species_summs <- species_dat %>% 
    group_by(Year, length_cm) %>% 
    summarise(n_lengths = n(),
              survey_catch = sum(numlen_adj),
              survey_bio = sum(biom_per_lclass),
              strat_catch = sum(strat_total_abund_s),
              strat_bio = sum(strat_total_lwbio_s),
              .groups = "drop") %>% 
    mutate(yr = as.numeric(as.character(Year)))
  
  # Abundances
  p1 <- species_summs  %>% 
    group_by(yr) %>% 
    summarise(strat_catch = sum(strat_catch)/1e6, .groups = "drop") %>% 
    ggplot(aes(x = yr, y = strat_catch)) +
    geom_segment(aes(xend = yr, yend = 0)) + 
    geom_point(aes(color = "strat_abund"), show.legend = FALSE) + 
    scale_color_gmri() +
    scale_y_continuous(labels = comma_format()) +
    labs(y = "Stratified Abundance\n(millions)", title = comname_choice) +
    theme(axis.title.x = element_blank())
  
  # Biomass
  p2 <- species_summs  %>% 
    group_by(yr) %>% 
    summarise(strat_bio = sum(strat_bio)/1e6, .groups = "drop") %>% 
    ggplot(aes(x = yr, y = strat_bio)) +
    geom_segment(aes(xend = yr, yend = 0)) + 
    geom_point(aes(color = "strat_bio"), show.legend = FALSE) + 
    scale_color_gmri() +
    scale_y_continuous(labels = comma_format()) +
    labs(y = "Stratified Biomass\n(million kg)") +
    theme(axis.title.x = element_blank())
  
  # Mean and 95% size
  size_summaries <- species_summs %>% 
    group_by(yr) %>% 
    summarise(mean_len = weighted.mean(length_cm, strat_catch),
              max_len = max(length_cm), 
              len_95 = wtd.quantile(length_cm, q = 0.95, weight = strat_catch, na.rm = T),
              len_75 = wtd.quantile(length_cm, q = 0.75, weight = strat_catch, na.rm = T),
              len_25 = wtd.quantile(length_cm, q = 0.25, weight = strat_catch, na.rm = T),
              len_05 = wtd.quantile(length_cm, q = 0.05, weight = strat_catch, na.rm = T),
              .groups = "drop")
  
  # scaling and displaying distribution
  p3 <- ggplot(size_summaries) +
    geom_point(aes(x = yr, y = mean_len, color = "Mean Length (cm)")) +
    geom_smooth(aes(yr, max_len, color = "Max Length (cm)"), formula = y ~ s(x, bs = "cs", k = 14), method = "gam", se = F) +
    geom_point(aes(x = yr, y = max_len, color = "Max Length (cm)")) +
    geom_smooth(aes(yr, mean_len, color = "Mean Length (cm)"), formula = y ~ s(x, bs = "cs", k = 14), method = "gam", se = F) +
    geom_point(aes(x = yr, y = len_05, color = "05pct. Length (cm)")) +
    geom_smooth(aes(yr, len_05, color = "05pct. Length (cm)"), formula = y ~ s(x, bs = "cs", k = 14), method = "gam", se = F) +
    geom_point(aes(x = yr, y = len_25, color = "25pct. Length (cm)")) +
    geom_smooth(aes(yr, len_25, color = "25pct. Length (cm)"), formula = y ~ s(x, bs = "cs", k = 14), method = "gam", se = F) +
    geom_point(aes(x = yr, y = len_75, color = "75pct. Length (cm)")) +
    geom_smooth(aes(yr, len_75, color = "75pct. Length (cm)"), formula = y ~ s(x, bs = "cs", k = 14), method = "gam", se = F) +
    geom_point(aes(x = yr, y = len_95, color = "95pct. Length (cm)")) +
    geom_smooth(aes(yr, len_95, color = "95pct. Length (cm)"), formula = y ~ s(x, bs = "cs", k = 14), method = "gam", se = F) +
    labs(x = "Year", y = "Length (cm)", color = "") +
    scale_color_gmri() +
    theme(legend.position = "bottom")
  
  
  p_out <- (p1 / p2 / p3) 
  return(p_out)
}


# Test
plot_species_dist("atlantic cod")
plot_species_dist("haddock")
plot_species_dist("spiny dogfish")
plot_species_dist("winter skate")
plot_species_dist("acadian redfish")






####  Species Groups  ####

# function to plot a species guild
plot_guild_dist <- function(guild){
  
  
  guild_dat <- filter(nefsc_stratified, spec_class == guild)
  
  # Get yearly totals
  guild_summs <- guild_dat %>% 
    group_by(Year, length_cm) %>% 
    summarise(n_lengths = n(),
              survey_catch = sum(numlen_adj),
              survey_bio = sum(biom_per_lclass),
              strat_catch = sum(strat_total_abund_s),
              strat_bio = sum(strat_total_lwbio_s),
              .groups = "drop") %>% 
    mutate(yr = as.numeric(as.character(Year)))
  
  # No Seasons
  p1 <- guild_summs  %>% 
    group_by(yr) %>% 
    summarise(strat_catch = sum(strat_catch)/1e6, .groups = "drop") %>% 
    ggplot(aes(x = yr, y = strat_catch)) +
    geom_segment(aes(xend = yr, yend = 0)) + 
    geom_point(aes(color = "strat_abund"), show.legend = FALSE) + 
    scale_color_gmri() +
    scale_y_continuous(labels = comma_format()) +
    labs(y = "Stratified Abundance\n(millions)", title = guild) +
    theme(axis.title.x = element_blank())
  
  # Biomass
  p2 <- guild_summs  %>% 
    group_by(yr) %>% 
    summarise(strat_bio = sum(strat_bio)/1e6, .groups = "drop") %>% 
    ggplot(aes(x = yr, y = strat_bio)) +
    geom_segment(aes(xend = yr, yend = 0)) + 
    geom_point(aes(color = "strat_bio"), show.legend = FALSE) + 
    scale_color_gmri() +
    scale_y_continuous(labels = comma_format()) +
    labs(y = "Stratified Biomass\n(million kg)") +
    theme(axis.title.x = element_blank())
  
  # Mean and 95% size
  size_summaries <- guild_summs %>% 
    group_by(yr) %>% 
    summarise(mean_len = weighted.mean(length_cm, strat_catch),
              max_len = max(length_cm), 
              len_95 = wtd.quantile(length_cm, q = 0.95, weight = strat_catch, na.rm = T),
              len_75 = wtd.quantile(length_cm, q = 0.75, weight = strat_catch, na.rm = T),
              len_25 = wtd.quantile(length_cm, q = 0.25, weight = strat_catch, na.rm = T),
              len_05 = wtd.quantile(length_cm, q = 0.05, weight = strat_catch, na.rm = T),
              .groups = "drop")
  
  # scaling and displaying distribution
  p3 <- ggplot(size_summaries) +
    geom_point(aes(x = yr, y = mean_len, color = "Mean Length (cm)")) +
    geom_smooth(aes(yr, mean_len, color = "Mean Length (cm)"), formula = y ~ s(x, bs = "cs", k = 14), method = "gam", se = F) +
    # geom_smooth(aes(yr, max_len, color = "Max Length (cm)"), formula = y ~ s(x, bs = "cs", k = 14), method = "gam", se = F) +
    # geom_point(aes(x = yr, y = max_len, color = "Max Length (cm)")) +
    geom_point(aes(x = yr, y = len_05, color = "05pct. Length (cm)")) +
    geom_smooth(aes(yr, len_05, color = "05pct. Length (cm)"), formula = y ~ s(x, bs = "cs", k = 14), method = "gam", se = F) +
    geom_point(aes(x = yr, y = len_25, color = "25pct. Length (cm)")) +
    geom_smooth(aes(yr, len_25, color = "25pct. Length (cm)"), formula = y ~ s(x, bs = "cs", k = 14), method = "gam", se = F) +
    geom_point(aes(x = yr, y = len_75, color = "75pct. Length (cm)")) +
    geom_smooth(aes(yr, len_75, color = "75pct. Length (cm)"), formula = y ~ s(x, bs = "cs", k = 14), method = "gam", se = F) +
    geom_point(aes(x = yr, y = len_95, color = "95pct. Length (cm)")) +
    geom_smooth(aes(yr, len_95, color = "95pct. Length (cm)"), formula = y ~ s(x, bs = "cs", k = 14), method = "gam", se = F) +
    labs(x = "Year", y = "Length (cm)", color = "") +
    scale_color_gmri() +
    guides(color = guide_legend(nrow = 2)) +
    theme(legend.position = "bottom")
  
  p_out <- (p1 / p2 / p3) 
  return(p_out)
}



####  Species Guild Plots  ####


#####  Groundfish  ####
plot_guild_dist("Groundfish")


#####  Elasmobranchs  ####
plot_guild_dist("Elasmobranch")



#####  Pelagics  ####
plot_guild_dist("Pelagic")


#####  Demersals  ####
plot_guild_dist("Demersal")

####  Single Species Ridgeline Plots  ####


#####  Cod    ####

# filter cod
cod_dat <- filter(nefsc_stratified, comname == "atlantic cod")

# Summarize
cod_summs <- cod_dat %>% 
  group_by(Year, season, length_cm) %>% 
  summarise(n_lengths = n(),
            survey_catch = sum(numlen_adj),
            strat_catch = sum(strat_total_abund_s),
            .groups = "drop")

(p1 | p2) + plot_annotation(title = "cod")

# Season Facets
cod_summs  %>% 
  ggplot(aes(x = length_cm, y = Year, height = strat_catch)) +
  geom_density_ridges(stat = "identity", scale=1.25) +
  facet_wrap(~season) +
  labs(x = "Length (cm)", y = "Year")

# scaling and displaying distribution
ggplot(cod_summs, aes(x = length_cm, y = Year))+
  geom_density_ridges(alpha = 0.7, scale = 1.25, rel_min_height = 0.02, quantile_lines = TRUE, quantiles = 2)+
  facet_wrap(~season) +
  labs(x = "Length (cm)", y = "Year")

#####  Haddock   ####

hadd_dat <- filter(nefsc_stratified, comname == "haddock")

hadd_summs <- hadd_dat %>% 
  group_by(Year, season, length_cm) %>% 
  summarise(n_lengths = n(),
            survey_catch = sum(numlen_adj),
            strat_catch = sum(strat_total_abund_s),
            .groups = "drop")

hadd_summs  %>% 
  ggplot(aes(x = length_cm, y = Year, height = strat_catch)) +
  geom_density_ridges(stat = "identity", scale = 1.25) +
  facet_wrap(~season) +
  labs(x = "Length (cm)", y = "Year")

# scaling and displaying distribution
ggplot(hadd_summs, aes(x = length_cm, y = Year))+
  geom_density_ridges(alpha = 0.7, scale = 1.25, rel_min_height = 0.02, quantile_lines = TRUE, quantiles = 2)+
  facet_wrap(~season) +
  labs(x = "Length (cm)", y = "Year")

