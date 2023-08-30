# Area Stratification, Effort Check


# Effort should be allocated based on area of strata



# It is "possible" to change effort to better sample strata with high variance,
# or to not succesfully sample it and have lower effort
# If the effort isn't proportional to area, then
# metrics of habitat use can be biased



# PAckages
library(targets)
library(tidyverse)
library(scales)
library(gmRi)
library(gt)


####  Goal  ####
# targets data - match the decadal variability inputs
# should have stratum area in there
tar_load(catch_complete)


# Save for Kathy, she can read in here with read_csv()
# write_csv(catch_complete, here::here("data/nmfs_survdat_clean.csv"))
glimpse(catch_complete)


# Is sampling proportional to area?
# catch_complete$strat_ntows is number of tows in a strata within a season
# also have area in km, total area of strata sampled that season, and the ratio between those two
# that ratio should be ~~ 
strata_effort <- distinct(catch_complete, stratum, est_year, season, strat_ntows, s_area_km2, tot_s_area, st_ratio)


# Get the total seasonal effort:
strata_effort <- strata_effort %>% 
  group_by(est_year, season) %>% 
  summarise(season_total_tows = sum(strat_ntows, na.rm = T),
            .groups = "drop") %>% 
  left_join(strata_effort, join_by(est_year, season)) %>% 
  mutate(strata_effort_ratio = strat_ntows/season_total_tows)


# Is it different in some years?
ggplot(strata_effort, aes(est_year, strata_effort_ratio, group = stratum)) +
  geom_line() +
  facet_grid(stratum~season, scales = "free") +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(x = "Year", y = "Effort in Strata (ntows) / total effort (ntows all strata)    |    Done each season independently")



# Scale it or something, maybe just make a table...
strata_effort <- strata_effort %>% 
  group_by(stratum, season) %>% 
  summarise(avg_effort_proportion = mean(strata_effort_ratio, na.rm = T),
            effort_proportion_sd = sd(strata_effort_ratio, na.rm = T),
            .groups = "drop") %>% 
  left_join(strata_effort) 


strata_effort %>% 
  mutate(
    percent_diff = ((strata_effort_ratio - avg_effort_proportion) / avg_effort_proportion)*100,
    diff_z =  (strata_effort_ratio - avg_effort_proportion) / effort_proportion_sd ) %>% 
  ggplot(aes(stratum, percent_diff)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Stratum", y = "Percent Difference from All-Year Average", 
       title = "Difference in Relative Sampling Effort by Year",
       subtitle = "Number of tows in each strata / tows across all strata")



strata_effort %>% 
  distinct(est_year, season, stratum, strata_effort_ratio, st_ratio) %>%
  mutate(after_2010 = ifelse(est_year >2009, "Post-2009", "Pre-2010")) %>% 
  ggplot() +
    geom_point(aes(st_ratio, strata_effort_ratio, color = after_2010)) +
    geom_abline(intercept = 0, slope = 1) +
    scale_color_manual(values = c("orange", "gray")) +
    labs(
      x = "Stratum's Percentage of Total Area",
      y = "Stratam's Percentage of Total Effort")
  


# ?debug
strata_effort %>% 
  distinct(stratum, avg_effort_proportion) %>% 
  mutate(stratum = fct_reorder(stratum, avg_effort_proportion)) %>% 
  ggplot(aes(y = stratum, avg_effort_proportion)) +
  geom_col() +
  scale_x_continuous(labels = label_percent())
  


strata_effort %>% 
  distinct(stratum, avg_effort_proportion) %>% 
  mutate(stratum = fct_reorder(stratum, avg_effort_proportion)) %>% 
  ggplot(aes(y = stratum, avg_effort_proportion)) +
  geom_col() +
  scale_x_continuous(labels = label_percent())



####  T-testing  ####
# Is the fraction of the effort before/after 2010

test_data <- strata_effort %>% 
  filter(est_year %in% c(1970:2019)) %>% 
  mutate(decade = ifelse(est_year < 2010, 
                         "1970-2009", 
                         "2010-2019")) %>% 
  distinct(est_year, season, stratum, strata_effort_ratio, decade)


# Do it for each strata independently
effort_differences <- test_data %>% 
  split(.$stratum) %>% 
  map_dfr(function(x){
    t.test(strata_effort_ratio ~ decade, data = x) %>% 
      tidy() %>% 
      select(
        effort_fraction_1970to2009 = estimate1, 
        effort_fraction_2010to2019 = estimate2, 
        method, 
        p.value) %>% 
      mutate(different = ifelse(p.value <= 0.05, T, F))
  }, .id = "stratum")


# How many ?
effort_differences %>% filter(different)  %>% 
  mutate(across(where(is.numeric), \(x) round(x, 3))) %>% 
  arrange(stratum) %>% 
  gt::gt()

# What percent
mean(effort_differences$different)  
sum(effort_differences$different)
nrow(effort_differences)

