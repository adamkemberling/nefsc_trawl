# calanus vs cpr community PCA's
library(tidyverse)
library(gmRi)
library(ecodata)




# From CPR Analysis:
# Uses 7 focal taxa, 2 calanus stages
cpr_indices <- read_csv("~/Documents/Repositories/continuous_plankton_recorder/results_data/cpr_focal_pca_timeseries_period_1961-2017.csv")


# Do some Reshaping
cpr_indices <- cpr_indices %>% 
  select(c(year, 2:3)) %>% 
  rename(
    `Principal Component 1` = "First Mode",
    `Principal Component 2` = "Second Mode") %>% 
  pivot_longer(names_to = "pca_mode", 
               values_to = "pca_loading", 
               cols = c(2:3))



cpr_idx <- cpr_indices %>% 
  mutate(area = "Gulf of Maine",
         var = str_c("cpr_pca_", str_sub(pca_mode, -1, -1))) %>% 
  group_by(year, area, var) %>% 
  summarise(value = mean(pca_loading),
            .groups = "drop")


# Zooplankton densities from ECOMON?
# tidy the calanus stage abundance data
cal <- ecodata::calanus_stage %>% 
  rename(year = Time,
         var = Var,
         area = EPU) %>% 
  group_by(var, area) %>% 
  mutate(value = as.numeric(scale(Value))) %>% 
  ungroup() %>% 
  select(-Value)


# Put them together with cpr indices to plot
cal_test <- left_join(
  cal, 
  pivot_wider(cpr_idx, values_from = value, names_from = var) %>% select(-area))


# Plot them somehow
ggplot(cal_test, aes(year, value, group = var)) +
  geom_col(fill = "gray60") +
  geom_line(aes(year, cpr_pca_1, color = "PCA 1")) +
  geom_line(aes(year, cpr_pca_2, color = "PCA 2")) +
  scale_color_gmri() +
  facet_grid(area~var)


# Ignore years, correlations
ggplot(cal_test, aes(cpr_pca_1, value)) +
  geom_point() +
  geom_smooth(formula = y~x, method = "lm") +
  facet_grid(area~var, scales = "free") +
  labs(y = "Z-Scaled Density")

# PCA 2 (the  "calanus" mode)
ggplot(cal_test, aes(cpr_pca_2, value)) +
  geom_point() +
  geom_smooth(formula = y~x, method = "lm") +
  facet_grid(area~var, scales = "free") +
  labs(y = "Z-Scaled Density")
