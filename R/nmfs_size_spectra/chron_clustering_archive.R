####  

# Chronological Clustering Code
# Removed from size spectrum biomass exploration
# isolated here for reference when we decided to explore changepoints:




####  Packages  ####
library(targets)
library(here)
library(sf)
library(gmRi)
library(patchwork)
library(rioja)
library(vegan)
library(gt)
library(knitr)
library(tidyverse)

# Support functions
source(here("R/support/sizeSpectra_support.R"))

theme_set(theme_minimal() + 
            theme(legend.position = "bottom", 
                  plot.title = element_text(hjust = 0.5)))



####  Data  ####

# OISST for all the regions
withr::with_dir(rprojroot::find_root('_targets.R'), 
                tar_load(regional_oisst))     

# SS and manual bins together
withr::with_dir(rprojroot::find_root('_targets.R'), 
                tar_load(size_spectrum_indices))   

# Format Columns
size_spectrum_indices <- size_spectrum_indices  %>% 
  mutate(season = fct_rev(season),
         survey_area = factor(survey_area, 
                              levels = c("GoM", "GB", "SNE", "MAB")),
         yr = as.numeric(as.character(Year)))

# 1. Biological data used as input
withr::with_dir(rprojroot::find_root('_targets.R'), 
                tar_load(nefsc_1g))            

# quick little format
nefsc_1g <- nefsc_1g %>% 
  mutate(fishery = case_when(
    fishery == "com" ~ "Commercially Targeted",
    fishery == "nc" ~ "Not Commercially Targeted",
    TRUE ~ "Not Labelled"))


#### Prep  ####

# Pull the group ID for the slopes grouped on year, season, and region
warmem_group_slopes <- size_spectrum_indices %>% 
  filter(`group ID`== "single years * season * region")


# Or just regions and years
year_region_slopes <- size_spectrum_indices %>% 
  filter(`group ID` == "single years * region")


####  Exploration  ####

# Just years and regions
(ss_patterns <- year_region_slopes %>% 
   ggplot(aes(yr, b, color = survey_area)) +
   geom_line(aes(group = survey_area), linetype = 3) +
   geom_point(alpha = 0.6) +
   geom_smooth(method = "lm", formula = y ~ x, 
               alpha = 0.2) +
   facet_grid(survey_area~.) +
   scale_color_gmri() +
   labs(x = "", 
        y = "Size Spectrum Slope (b)", 
        color = "",
        title = "Results from Area-Stratified Abundance - Edwards Method"))


####  Clustering  ####

#### Chronological Clustering

# Prep and Clustering Code# 
  
  
# Reshaping
# Goal: Row for Years, columns for each different group
cluster_dat <- year_region_slopes %>% 
  select(Year, survey_area, b) 

# Pivot wider
ss_dat <- cluster_dat %>% 
  pivot_wider(names_from = c(survey_area), 
              values_from = b, 
              names_sep = "_") %>% 
  column_to_rownames(var = "Year")





# Function to run chronological clustering
run_chron_clust <- function(in_dat){
  
  # Get Euclidean Distances
  eucdist <- vegdist(in_dat,
                     method = "euclidean",
                     binary = FALSE,
                     diag = FALSE,
                     upper = FALSE,
                     na.rm = TRUE)
  
  # Perform Chronological Clustering on distances
  cl <- chclust(eucdist, method = "coniss")
  
  # Return list
  return(list(eucdist = eucdist, cl = cl))
}






# Run for chronological clustering for each area
regions <- setNames(names(ss_dat), names(ss_dat))
region_clusters <- map(regions, function(region_name){
  ss_chron <- run_chron_clust(in_dat = select(ss_dat, one_of(region_name)))
  return(ss_chron)
})





# Broomstick Plot# 
  
  
# broken stick model
par(mfrow = c(2,2))
imap(region_clusters, function(region_clust, region_id){
  vegan::bstick(region_clust$cl, plot=T)
  title(main = region_id)
})

par(mfrow = c(1,1))


# Dendrogram# 
  
  
par(mfrow = c(2,2))
imap(region_clusters, function(region_clust, region_id){
  plot(region_clust$cl, 
       hang = -0.1, 
       axes = FALSE, 
       cex = 0.8, 
       horiz = F)
  axis(side = 2, cex.axis = 1)
  title(main = paste0("sizeSpectra Clustering - ", region_id), cex = 1)
  mtext(side = 2, line = 2.3, "Sum of squares", cex = 1, las = 0)
})
par(mfrow = c(1,1))




#### Slopes with Breakpoints




# put breakpoints in table
ss_breaks <- list(
  "GB"  = data.frame(year = c(2008, 2003, 2002, 1980, 1982) + 0.5),
  "GoM" = data.frame(year = c(1984, 1999, 1997, 2008, 2015) + 0.5),
  "SNE" = data.frame(year = c(1992, 1994, 1990, 1987, 2012) + 0.5),
  "MAB" = data.frame(year = c(2005, 2013, 2012, 2007, 2001) + 0.5)
  
) %>%  map_dfr(., ~ .x, .id = "survey_area")  %>% 
  mutate(survey_area = factor(survey_area, levels = c("GoM", "GB", "SNE", "MAB")))


# original code
ss_patterns +
  geom_vline(data = ss_breaks, aes(xintercept = year), color = "gray40", linetype = 2, size = 0.7)






### Year, Region, and Season {.tabset}


#### Long-Term Patterns



# Plot the sizeSpectra slopes - year*region*season
(ss_patterns <- warmem_group_slopes %>% 
    ggplot(aes(yr, b, color = survey_area)) +
    geom_line(aes(group = survey_area), linetype = 3) +
    geom_point(alpha = 0.6) +
    geom_smooth(formula = y ~ x, method = "lm", alpha = 0.2) +
    facet_grid(survey_area~season) +
    scale_color_gmri() +
    labs(x = "", 
         y = "Size Spectrum Slope (b)", 
         color = "",
         title = "Results from Area-Stratified Abundance - Edwards Method") 
)







#### Chronological Clustering

# Prep and Clustering Code# 
  
  
# Reshaping
# Goal: Row for Years, columns for each different group
cluster_dat <- warmem_group_slopes %>% 
  select(Year, season, survey_area, b) 

# Pivot wider
ss_dat <- cluster_dat %>% 
  pivot_wider(names_from = c(season, survey_area), 
              values_from = b, 
              names_sep = "_") %>% 
  column_to_rownames(var = "Year")





# Function to run chronological clustering
run_chron_clust <- function(in_dat){
  
  # Get Euclidean Distances
  eucdist <- vegdist(in_dat,
                     method = "euclidean",
                     binary = FALSE,
                     diag = FALSE,
                     upper = FALSE,
                     na.rm = TRUE)
  
  # Perform Chronological Clustering on distances
  cl <- chclust(eucdist, method = "coniss")
  
  # Return list
  return(list(eucdist = eucdist, cl = cl))
}


# Run for sizeSpectra Results
ss_chron <- run_chron_clust(in_dat = ss_dat)





# Broomstick Plot# 
  
  
# broken stick model
vegan::bstick(ss_chron$cl, plot=T)


# Dendrogram# 
  
  

plot(ss_chron$cl, 
     hang = -0.1, 
     axes = FALSE, 
     cex = 0.8, 
     horiz = F)
axis(side = 2, cex.axis = 1)
title("sizeSpectra Clustering Metrics", cex = 1)
mtext(side = 2, line = 2.3, "Sum of squares", cex = 1, las = 0)



#### Slopes with Breakpoints



ss_patterns +
  geom_vline(aes(xintercept = 1987.5), color = "gray40", linetype = 2, size = 0.7) +
  geom_vline(aes(xintercept = 2008.5), color = "gray40", linetype = 2, size = 0.7) +
  geom_vline(aes(xintercept = 2004.5), color = "gray40", linetype = 2, size = 0.7) +
  geom_vline(aes(xintercept = 2001.5), color = "gray40", linetype = 2, size = 0.7) +
  geom_vline(aes(xintercept = 1973.5), color = "gray40", linetype = 2, size = 0.7) +
  geom_vline(aes(xintercept = 1974.5), color = "gray40", linetype = 2, size = 0.7) 





## Log10 Binning {.tabset}


### Slope {.tabset}

#### Long-Term Patterns



# plot trends of log10 slopes
(l10_patterns <- warmem_group_slopes %>% 
    ggplot(aes(yr, l10_slope_strat, color = survey_area)) +
    geom_line(aes(group = survey_area), linetype = 3) +
    geom_point(alpha = 0.6) +
    geom_smooth(formula = y ~ x, method = "lm", alpha = 0.2) +
    facet_grid(survey_area~season) +
    scale_color_gmri() +
    labs(x = "", 
         y = "Size Spectrum Slope (b)", 
         color = "",
         title = "Results from Area-Stratified Abundance - log10 bins"))





#### Chronological Clustering

# Prep and Clustering Code# 
  
  
# Reshaping
# Goal: Row for Years, columns for each different group
l10_cluster_dat <- warmem_group_slopes %>% 
  select(Year, season, survey_area, l10_slope_strat) 

# Pivot wider
l10_dat <- l10_cluster_dat %>% 
  pivot_wider(names_from = c(season, survey_area), 
              values_from = l10_slope_strat, 
              names_sep = "_") %>% 
  column_to_rownames(var = "Year")

# Run clustering
l10_clust <- run_chron_clust(in_dat = l10_dat)



# Broomstick Plot# 
  
  
# broken stick model
vegan::bstick(l10_clust$cl, plot=T)


# Dendrogram# 
  
  

plot(l10_clust$cl, 
     hang = -0.1, 
     axes = FALSE, 
     cex = 0.8, 
     horiz = F)
axis(side = 2, cex.axis = 1.3)
title("Log10 Size Spectrum Clustering Metrics", cex = 1)
mtext(side = 2, line = 2.3, "Sum of squares", cex = 1, las = 0)




#### Slopes with Breakpoints



l10_patterns +
  geom_vline(aes(xintercept = 1987.5), color = "gray40", linetype = 2, size = 0.7) +
  geom_vline(aes(xintercept = 1996.5), color = "gray40", linetype = 2, size = 0.7) +
  geom_vline(aes(xintercept = 2008.5), color = "gray40", linetype = 2, size = 0.7) 




### Intercept {.tabset}

#### Long-Term Patterns



# plot trends of log10 slopes
(int_patterns <- warmem_group_slopes %>% 
    ggplot(aes(yr, l10_int_strat, color = survey_area)) +
    geom_line(aes(group = survey_area), linetype = 3) +
    geom_point(alpha = 0.6) +
    geom_smooth(formula = y ~ x, method = "lm", alpha = 0.2) +
    facet_grid(survey_area~season) +
    scale_color_gmri() +
    labs(x = "", 
         y = "Size Spectrum Intercept", 
         color = "",
         title = "Results from Area-Stratified Abundance - log10 bins")
)



#### Chronological Clustering

# Prep and Clustering Code# 
  
  
# Reshaping
# Goal: Row for Years, columns for each different group
l10_intercept_dat <- warmem_group_slopes %>% 
  select(Year, season, survey_area, l10_int_strat) 

# Pivot wider
intercept_dat <- l10_intercept_dat %>% 
  pivot_wider(names_from = c(season, survey_area), 
              values_from = l10_int_strat, 
              names_sep = "_") %>% 
  column_to_rownames(var = "Year")

# Run Clustering
l10_ints <- run_chron_clust(in_dat = intercept_dat)



# Broomstick Plot# 
  
  
# broken stick model
vegan::bstick(l10_ints$cl, plot=T)


# Dendrogram# 
  
  

plot(l10_ints$cl, 
     hang = -0.1, 
     axes = FALSE, 
     cex = 0.8, 
     horiz = F)
axis(side = 2, cex.axis = 1.3)
title("Log10 Bins Intercept - Clustering Metrics", cex = 1)
mtext(side = 2, line = 2.3, "Sum of squares", cex = 1, las = 0)




#### Intercepts with Breakpoints



# Re-plot patterns with breaks
int_patterns +
  geom_vline(aes(xintercept = 1987.5), color = "gray50", linetype = 2, size = 0.7) +
  geom_vline(aes(xintercept = 2008.5), color = "gray50", linetype = 2, size = 0.7) +
  geom_vline(aes(xintercept = 1995.5), color = "gray50", linetype = 2, size = 0.7) +
  geom_vline(aes(xintercept = 1996.5), color = "gray50", linetype = 2, size = 0.7) +
  labs(subtitle = "Chronological Clustering Breakpoints")


### Fit Statidtics (adjusted r-squared) 

# For each of the slope/intercepts derived using a linear model and binned data I also pulled the adjusted r-square to get a sense of whether or not certain groups had poor fits that should be investigated.


# Reshaping
# Goal: Row for Years, columns for each different group
l10_fit_dat <- warmem_group_slopes %>% 
  select(Year, season, survey_area, l10_slope_strat, l10_int_strat, l10_rsqr_strat) %>% 
  mutate(yr = as.numeric(Year))


l10_fit_dat %>% 
  ggplot(aes(yr, l10_rsqr_strat, color = survey_area)) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = 0, ymax = 0.2, 
            fill = "gray80", color = "transparent") +
  geom_hline(yintercept = 0.2, linetype = 2, size = 0.5, color = "gray50") +
  geom_line(aes(group = survey_area), linetype = 3) +
  geom_point(alpha = 0.6) +
  geom_smooth(formula = y ~ x, method = "lm", alpha = 0.2) +
  facet_grid(survey_area~season) +
  scale_color_gmri() +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  labs(x = "", 
       y = "Linear Model R-Square", 
       color = "",
       title = "Results from Area-Stratified Abundance - log10 bins")








##  Sea Surface Temperature {.tabset}

### Long-Term Patterns



# plot trends of log10 slopes
(sst_patterns <- regional_oisst %>% 
    ggplot(aes(yr, sst_anom, color = survey_area)) +
    geom_line(aes(group = survey_area), linetype = 3) +
    geom_point(alpha = 0.6) +
    geom_smooth(formula = y ~ x, method = "lm", alpha = 0.2) +
    facet_grid(survey_area~.) +
    scale_color_gmri() +
    labs(x = "", 
         y = "Mean Temperature Anomaly", 
         color = "",
         title = "")
)




### Chronological Clustering

# Prep and Clustering Code# 
  
  
# Pivot OISST Wider
oisst_dat <- regional_oisst %>% 
  select(yr, survey_area, sst_anom) %>% 
  pivot_wider(names_from = survey_area, values_from = sst_anom) %>% 
  column_to_rownames(var = "yr")

# Run clustering
sst_clust <- run_chron_clust(oisst_dat)





# Broomstick Plot# 
  
  
# broken stick model
vegan::bstick(sst_clust$cl, plot=T)


# Dendrogram# 
  
  

plot(sst_clust$cl, 
     hang = -0.1, 
     axes = FALSE, 
     cex = 0.8, 
     horiz = F)
axis(side = 2, cex.axis = 1.3)
title("SST Anomalies - Clustering Metrics", cex = 1)
mtext(side = 2, line = 2.3, "Sum of squares", cex = 1, las = 0)




### Anomalies with Breakpoints



# Re-plot patterns with breaks
sst_patterns +
  geom_vline(aes(xintercept = 2011.5), color = "gray50", linetype = 2, size = 0.7) +
  geom_vline(aes(xintercept = 1998.5), color = "gray50", linetype = 2, size = 0.7) +
  geom_vline(aes(xintercept = 2002.5), color = "gray50", linetype = 2, size = 0.7) +
  geom_vline(aes(xintercept = 2005.5), color = "gray50", linetype = 2, size = 0.7) 


# Biomass Patterns with Breakpoints {.tabset}

# This section seeks to track different aspects of the biomass distribution during each of the breaks determined from the break point analyses.

# Primarily
# 1. How is biomass changing for the different functional groups   
# 2. What is the relative balance of biomass across functional groups, and within them, size bins




# Make groups based on breaks
edwards_breaks <- list("1970-1982" = 1970:1982,
                       "1983-2003" = 1983:2003,
                       "2004-2008" = 2004:2008,
                       "2009-2016" = 2009:2016,
                       "2017-2020" = 2017:2020)



temp_breaks <- list("1981-1998" = 1981:1998, 
                    "1999-2002" = 1999:2002, 
                    "2003-2005" = 2003:2005,
                    "2006-2011" = 2006:2011, 
                    "2012-2021" = 2012:2021)


####  Digging into Actual Biomass data  ####


# Cut up discrete length and weight bins
nefsc_size_bins <- nefsc_1g %>% 
  mutate(
    length_bin = case_when(
      length_cm <= 5   ~ "0 - 5cm",
      length_cm <= 10  ~ "5 - 10cm",
      length_cm <= 25  ~ "10 - 25cm",
      length_cm <= 50  ~ "25 - 50cm",
      length_cm <= 75  ~ "50 - 75cm",
      length_cm <= 100 ~ "75 - 100cm",
      length_cm >= 100 ~ ">= 100cm"),
    length_bin = factor(length_bin, levels = c(
      "0 - 5cm", "5 - 10cm", "10 - 25cm",
      "25 - 50cm", "50 - 75cm", "75 - 100cm", ">= 100cm")))


# Weight bins
nefsc_size_bins <- nefsc_size_bins %>% 
  mutate(
    weight_bin = case_when(
      ind_weight_kg <= 0.001 ~ "0 - 1g",
      ind_weight_kg <= 0.005 ~ "1 - 5g",
      ind_weight_kg <= 0.010 ~ "5 - 10g",
      ind_weight_kg <= 0.050 ~ "10 - 50g",
      ind_weight_kg <= 0.100 ~ "50 - 100g",
      ind_weight_kg <= 0.500 ~ "100 - 500g",
      ind_weight_kg <= 1.000 ~ ".5 - 1kg",
      ind_weight_kg <= 5.000 ~ "1 - 2kg",
      ind_weight_kg <= 5.000 ~ "2 - 5kg",
      ind_weight_kg <= 10.00 ~ "5 - 10kg",
      ind_weight_kg >= 10.00 ~ ">= 10kg" ),
    weight_bin = factor(weight_bin, levels = c(
      "0 - 1g", "1 - 5g", "5 - 10g", "10 - 50g",
      "50 - 100g", "100 - 500g", ".5 - 1kg",
      "1 - 2kg", "2 - 5kg", "5 - 10kg", ">= 10kg")))


# Rename the functional groups
nefsc_size_bins <- nefsc_size_bins %>% 
  mutate(
    spec_class = case_when(
      spec_class == "gf"  ~ "Groundfish",
      spec_class == "pel" ~ "Pelagic",
      spec_class == "dem" ~ "Demersal",
      spec_class == "el"  ~ "Elasmobranch",
      spec_class == "<NA>" ~ "NA"))

# Make regions go N->S
nefsc_size_bins <- nefsc_size_bins %>% 
  mutate(survey_area = factor(survey_area, levels = c("GoM", "GB", "SNE", "MAB")),
         season = factor(season, levels = c("Spring", "Fall")))



## label the edwards bins
nefsc_size_bins <- imap_dfr(edwards_breaks, function(x,y){
  nefsc_size_bins %>% 
    filter(Year %in% x) %>% 
    mutate(e_break = y,
           e_break = factor(e_break, levels = names(edwards_breaks)))
}) %>% arrange(Year)

# label the temp bins
nefsc_size_bins <- imap_dfr(temp_breaks, function(x,y){
  nefsc_size_bins %>% 
    filter(Year %in% x) %>% 
    mutate(t_break = y,
           t_break = factor(t_break, levels = names(temp_breaks)))
}) %>% arrange(Year)






# Track how things are changing across slope exponent breaks
functional_group_splits <- nefsc_size_bins %>% 
  split(.$e_break) %>% 
  imap(function(e_data, label){
    
    # 1. How is biomass trending within each functional group
    
    # Get that yearly summary
    yearly_data <- e_data %>% 
      group_by(Year, spec_class) %>% 
      summarise(
        total_strat_biom = sum(strat_total_lwbio_s), 
        avg_ind_mass = weighted.mean(ind_weight_kg, numlen),
        avg_ind_len = weighted.mean(length_cm, numlen),
        .groups = "drop")
    
    # run linear models
    biom_trend_lm   <- broom::tidy(lm(total_strat_biom ~ Year, 
                                      data = yearly_data))
    length_trend_lm <- broom::tidy(lm(avg_ind_len ~ Year, 
                                      data = yearly_data))
    imass_trend_lm <- broom::tidy(lm(avg_ind_mass ~ Year, 
                                     data = yearly_data))
    
    # report direction
    biom_slope     <- biom_trend_lm[[2,"estimate"]]
    len_slope      <- length_trend_lm[[2,"estimate"]]
    imass_slope      <- imass_trend_lm[[2,"estimate"]]
    biom_direction <- ifelse(biom_slope > 0, "increase", "decrease")
    len_direction  <- ifelse(len_slope > 0, "increase", "decrease")
    imass_direction  <- ifelse(imass_slope > 0, "increase", "decrease")
    
    # report significance
    biom_sig <- ifelse(biom_trend_lm[[2,"p.value"]] < 0.05, "significant", "not")
    len_sig  <- ifelse(length_trend_lm[[2,"p.value"]] < 0.05, "significant", "not")
    imass_sig  <- ifelse(imass_trend_lm[[2,"p.value"]] < 0.05, "significant", "not")
    
    # make table
    trend_table <- data.frame("ebreak"       = rep(label, 3),
                              "data_source"  = c("stratified_biomass", "avg_length", "avg_mass"),
                              "trend"        = c(biom_direction, len_direction, imass_direction),
                              "rate"         = c(biom_slope, len_slope, imass_slope),
                              "significance" = c(biom_sig, len_sig, imass_sig))
    
    # 2. Relative Biomass
    total_biom <- sum(e_data$strat_total_lwbio_s)
    relative_biom <- e_data %>% 
      group_by(spec_class) %>% 
      summarise(n_years = n_distinct(Year),
                group_biom_total = sum(strat_total_lwbio_s),
                avg_ann_biom = group_biom_total / n_years,
                biom_frac = group_biom_total / total_biom,
                avg_ind_mass = weighted.mean(ind_weight_kg, numlen),
                avg_ind_len = weighted.mean(length_cm, numlen),
                .groups = "drop") %>% 
      mutate(ebreak = label)
    
    # return tables
    return(list(#"biom_lm" = biom_trend_lm,
      #"len_lm" = length_trend_lm,
      "trends" = trend_table, 
      "rel_bio" = relative_biom))
    
  })


# Pull out trends as a table
fgroup_trends <- map_dfr(functional_group_splits, ~ .x[["trends"]])


# Pull out relative biomass as a table
fgroup_rel_bio <- map_dfr(functional_group_splits, ~ .x[["rel_bio"]])


## Significant Trends

# For each breakpoint period a linear regression was used to determine whether there was a net increase/decrease within that period. The results from all significant regressions (none) are reported below.


fgroup_trends %>% filter(significance == "significant")


## Relative Biomass

# The following figure tracks the relative distribution of catch across the different functional groups.



fgroup_rel_bio %>% 
  ggplot(aes(ebreak, y = biom_frac, fill = spec_class)) +
  geom_col(position = "dodge") +
  scale_fill_gmri() +
  labs(x = "Stanzas from ISD Exponent Breakpoint Analysis", 
       y = "Relative Proportion of Stratified Biomass",
       fill = "")


## Average Annual Biomass

# This figure visualizes the average total biomass (per year) for each functional group within the different size spectrum stanzas.



fgroup_rel_bio %>% 
  ggplot(aes(ebreak, y = avg_ann_biom/1e6, fill = spec_class)) +
  geom_col(position = "dodge") +
  scale_fill_gmri() +
  labs(x = "Stanzas from ISD Exponent Breakpoint Analysis", 
       y = "Avg. Annual Biomass (million kg)",
       fill = "")



## Average Individual Mass

# This figure visualizes the average individual bodymass for each functional group within the different size spectrum stanzas.

fgroup_rel_bio %>% 
  ggplot(aes(ebreak, y = avg_ind_mass, fill = spec_class)) +
  geom_col(position = "dodge") +
  scale_fill_gmri() +
  labs(x = "Stanzas from ISD Exponent Breakpoint Analysis", 
       y = "Average Individual Bodymass (kg)",
       fill = "")



## Average Individual Length

# This figure does the same, but for individual length, not bodymass.


fgroup_rel_bio %>% 
  ggplot(aes(ebreak, y = avg_ind_len, fill = spec_class)) +
  geom_col(position = "dodge") +
  scale_fill_gmri() +
  labs(x = "Stanzas from ISD Exponent Breakpoint Analysis", 
       y = "Average Individual Length (cm)",
       fill = "")




## Slope Check

# After seeing no large difference in the things above that we expected to see, I plotted the distribution of size spectrum slopes for each stanza. Lot of within gorup variation, not much across group variation.


warmem_group_slopes <- imap_dfr(edwards_breaks, function(x,y){
  warmem_group_slopes %>% 
    filter(Year %in% x) %>% 
    mutate(e_break = y,
           e_break = factor(e_break, levels = names(edwards_breaks)))
}) %>% arrange(Year)


warmem_group_slopes %>% 
  ggplot(aes(x = e_break, y = b, color = e_break)) +
  scale_color_gmri() +
  geom_boxplot() +
  labs(x = "Stanzas from ISD Exponent Breakpoint Analysis", 
       y = "Individual Size Distribution Exponent",
       color = "")






# # Needed a place to put this testing code
# 
# col_select <- function(data_in, abundance_vals = "stratified"){
#   
#  # toggle which abundance calue to use with switch
#   abundance_col <- switch(abundance_vals,
#     "observed"   = sym("numlen"),
#     "strat_mean" = sym("strat_mean_abund_s"),
#     "stratified" = sym("strat_total_abund_s"))
#   
#   
#   # Select just the columns we need:
#   # Rename them to match the code for sizeSpectra
#   data_out <- data_in %>% 
#     dplyr::select(
#       SpecCode = comname,
#       wmin = wmin_g,
#       wmax = wmax_g,
#       #Number = one_of(abundance_col))
#       Number = !!abundance_col)
#   
#   return(data_out)
# 
# }
# 
# # data to test with
# data_test <- nefsc_1g %>% slice(1:20)
# 
# 
# 
# col_select(data_in = data_test, abundance_vals = "stratified")
# col_select(data_in = data_test, abundance_vals = "observed")
# col_select(data_in = data_test, abundance_vals = "strat_mean")




