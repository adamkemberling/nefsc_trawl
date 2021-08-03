

####  Size Structure of Haddock, Cod and Species Guilds


####  Load Packages  ####
library(targets)
library(rnaturalearth)
library(gmRi)
library(patchwork)
library(tidyverse)
library(ggforce)
library(ggridges)

theme_set(theme_bw())


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
      TRUE ~ "Unclassified"
    ),
    season = factor(season, levels = c("Spring", "Fall")),
    Year = fct_rev(factor(est_year)),
    decade = floor_decade(est_year),
    decade = fct_rev(decade))



####  Cod  ####

cod_dat <- filter(nefsc_stratified, comname == "atlantic cod")


cod_summs <- cod_dat %>% 
  group_by(Year, season, length_cm) %>% 
  summarise(n_lengths = n(),
            survey_catch = sum(numlen_adj),
            strat_catch = sum(strat_total_abund_s),
            .groups = "drop")

# No Seasons
cod_summs  %>% 
  ggplot(aes(x = length_cm, y = Year, height = strat_catch)) +
  geom_density_ridges(stat = "identity", scale=1.25) +
  # facet_wrap(~season) +
  labs(x = "Length (cm)", y = "Year")

# scaling and displaying distribution
ggplot(hadd_summs, aes(x = length_cm, y = Year))+
  geom_density_ridges(alpha = 0.7, scale = 1.25, rel_min_height = 0.02, quantile_lines = TRUE, quantiles = 2)+
  # facet_wrap(~season) +
  labs(x = "Length (cm)", y = "Year")


# Season Facets
cod_summs  %>% 
  ggplot(aes(x = length_cm, y = Year, height = strat_catch)) +
  geom_density_ridges(stat = "identity", scale=1.25) +
  facet_wrap(~season) +
  labs(x = "Length (cm)", y = "Year")

# scaling and displaying distribution
ggplot(hadd_summs, aes(x = length_cm, y = Year))+
  geom_density_ridges(alpha = 0.7, scale = 1.25, rel_min_height = 0.02, quantile_lines = TRUE, quantiles = 2)+
  facet_wrap(~season) +
  labs(x = "Length (cm)", y = "Year")

####  Haddock  ####

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


####  Groundfish  ####



####  Elasmobranchs  ####




####  Pelagics  ####



####  Demersals  ####



