####  Decadal Variability  ####
#### Length at Age and Weight at Age Figures  
# Isolated from original nefsc_trawl repo for code availability with publication

# What Figures:

# Percent Change in Length at age: Figure 5
# Percent Change in Weight at Age: S8
# Average Length at Age: S9
# Average Weight at Age: S10
# Decadal Differences in Length at Age: Figure S11
# Decadal Differences in Weight at Age: Figure S12


####  Packages  ####

{
  library(car)
  library(scales)
  library(patchwork)
  library(directlabels)
  library(tidyverse)
  library(gmRi)
  library(targets)
  library(broom)
  library(grid)
  library(gridExtra)
}

# Path to decadal folder:
decadal_hires_folder <- cs_path(
  box_group = "mills", 
  subfolder = "Projects/Decadal Variability/Revisions/Growth/")


# Set seed
set.seed(3905)

# Celsius symbol
deg_c <- "\u00b0C"

# Where to put the figures
decadal_hires_folder <- cs_path("mills", "Projects/Decadal Variability/publication_figures/")

# Relative Path to Data
data_save_path <- here::here("data/decadal_variability_data/")

# Make them again with theme tweaks:
# From Carly
theme_ices <- function(...){
  theme_gmri() +
    theme(plot.title    = element_text(size = 8),
          axis.title    = element_text(size = 7),
          axis.text     = element_text(size = 6),
          panel.grid    = element_line(linewidth = 0.2),
          axis.line.x   = element_line(linewidth = 0.1),
          axis.ticks.x  = element_line(linewidth = 0.1),
          plot.margin   = margin(t = 8, b = 4, r = 8, l = 4),
          ...)
}

# Set it
theme_set(theme_ices())


####  Data  ####


# Access Biological Data with {targets}
withr::with_dir(rprojroot::find_root('_targets.R'), 
                tar_load(survdat_biological))

# Save this starting point:
write_csv(survdat_biological, str_c(data_save_path, "survdat_biological.csv"))


##### Control Options  ####

# Seasons to include in center of biomass
season_filter <- list(
  "cob" = c("Spring", "Fall"),
  "age" = c("Spring", "Fall"))


# Species to display with their own panels
vb_species <- sort(c(
  #"acadian redfish",
  "american plaice",
  "atlantic cod",
  "atlantic herring",
  "atlantic mackerel",
  "black sea bass",
  #"bluefish",
  "butterfish",
  "haddock", 
  "pollock",
  #"red hake",
  "scup",
  "silver hake",
  #"summer flounder",
  "white hake",
  #"windowpane",
  "winter flounder",
  "witch flounder",
  "yellowtail flounder"
))


# Get the groups we want to compare
year_groups <- list("1970-2009" = c(1970:2009),
                    "2010-2019" = c(2010:2019))




# Set the Age ranges for each species:
ages_to_test <- list(
  # "acadian redfish"     = c(1:15),
  "american plaice"     = c(1:9),
  "atlantic cod"        = c(1:7),
  "atlantic herring"    = c(1:8),
  "atlantic mackerel"   = c(1:5),
  "black sea bass"      = c(1:5),
  "butterfish"          = c(1:4),
  "haddock"             = c(1:10),
  "pollock"             = c(1:8),
  # "red hake"            = c(1:12),
  "scup"                = c(1:4),
  "silver hake"         = c(1:5),
  # "summer flounder"     = c(1:12),
  "white hake"          = c(1:7),
  # "windowpane"          = c(1:8),
  "winter flounder"     = c(1:8),
  "witch flounder"      = c(1:8),
  "yellowtail flounder" = c(1:6)
)

####  Data Prep  ####

##### 1. Bio Data  ####

# 1. Prepare Bio Data for Age analysis for regimes:

# Supplement the data
# 1. Add regime and decade labels
# 2. Add condition Facotr: Fulton's condition factor
# Drop NA age data, restrict to both summer and fall
bio_data <- survdat_biological %>%
  filter(is.na(age) == FALSE,
         season %in% season_filter$age) %>% 
  as.data.frame() %>% 
  mutate(
    yearclass = est_year - (age-1),
    regime = ifelse(est_year < 2010, "Early Regime", "Warm Regime"),
    decade = floor_decade(est_year),
    # Calculate Fultons condition factor:
    fK = (100 * indwt*1000) / length_cm^3) 


##### 2.  Growth Data Prep  #### 

# Filter to just the data for the two regimes
# Prep data separately for table
age_truncated <- bio_data %>% 
  filter(
    est_year >= 1970, est_year <= 2019,
    comname %in% vb_species,
    age > 0) %>% 
  complete(comname, est_year) %>% 
  mutate(regime = ifelse(est_year < 2010, "Early Regime", "Warm Regime"),
         age = ifelse(age>10, ">10", age)) %>% 
  mutate(age = fct_drop(age),
         age_tally = ifelse(!is.na(age), 1, 0))

# tidy up
rm(survdat_biological)
rm(bio_data)


## Apply an Outlier catcher here
# There are year1 herring over 10g that are obviously errors, they can be removed here
# And similarly for other problematic data points
age_truncated <- age_truncated %>% 
  mutate(indwt = ifelse(
    comname == "atlantic herring" & age == 1 & indwt > 10, 
    indwt*.001, 
    indwt))



# Get data split into the groups of years being tested
period_groups <- map(year_groups, function(period){
  # filter years:
  period_dat <- filter(age_truncated, est_year %in% period)
  return(period_dat)})


####  Analysis Functions  ####

# This function is looped through once for each ageclass


# Function to compare a single ageclass:


#' @title Age specific t-test 
#'
#' @param age_x Age class to filter out and test within
#' @param species_data dataframe containing the ageclass and dependent variable information
#' @param dependent_var String indicating dependent variable column being tested
#'
#' @return
#' @export
#'
#' @examples
ageclass_ttest <-  function(age_x, species_data = both_periods_df, dependent_var = {{test_var}}){
  
  # Subset the age class we're testing
  age_subset <- dplyr::filter(species_data, age == age_x)
  
  # Set formula to change dynamically to function input:
  test_formula <- formula(str_c(dependent_var, " ~ comparison_period"))
  
  # Run the comparison of means
  # Run Students/Welch's (depending on sample sizes)
  tidy_t <-  t.test(test_formula, data = age_subset) %>% 
    tidy() %>% 
    select(estimate1, estimate2, method, p.value)
  
  # Run a Bartletts
  tidy_b <- bartlett.test(test_formula, data = age_subset) %>% 
    tidy() %>% 
    select(variance_test = method, 
           bartlett_p =  p.value)
  
  # Make an output table
  bind_cols(
    list(
      data.frame(
        "comparison_var" = dependent_var,
        "statistical_test" = "two sample t.test"),
      tidy_t,
      tidy_b))}



# Function to take each species and compare some variable across comparison periods
# Reports the mean across both periods, and t-test results

#' @title Welch's T-Tests for Single Variable
#' 
#' @description This is a wrapper function used to handle different behaviors for
#' ageclass_ttest and manage outputs so that they can be stored in a single table.
#' 
#'
#' @param comname_x String indicating the common name of the species
#' @param group_1 String indicating time period 1, corresponds to a dataframe stored in list: "period_groups"
#' @param group_2 String indicating time period 2, corresponds to a dataframe stored in list: "period_groups"
#' @param test_var String indicating the dataframe variable to test
#'
#' @return
#' @export
#'
#' @examples
species_ageclass_ttests <- function(
    comname_x, 
    group_1 = "1970-2009", 
    group_2 = "2010-2019", 
    test_var = "length_cm"){
  
  
  
  # 1. Make names for the two groups/regimes
  regime_names <- c(group_1, group_2)
  
  
  # 2. Grab their data and put into a table, set the levels so group_1 is the first factor level
  both_periods_df <- bind_rows(
    list(filter(period_groups[[group_1]], comname == comname_x) %>% 
           drop_na(age, {{test_var}}),
         filter(period_groups[[group_2]], comname == comname_x) %>% 
           drop_na(age, {{test_var}})) %>% 
      setNames(regime_names),
    .id = "comparison_period") %>% 
    mutate(comparison_period = factor(comparison_period, levels = regime_names))
  
  
  # If the variable is indwt start it at 1992
  if(test_var == "indwt"){
    both_periods_df <- filter(both_periods_df, est_year >= 1992)}
  
  
  
  # 2. Run the ttest comparison for each ageclass:
  
  # Map the ages for that species
  ageclass_comparisons <- map_dfr(
    ages_to_test[[comname_x]], 
    possibly(
      .f = ~ageclass_ttest(.x, species_data = both_periods_df, dependent_var = test_var),
      
      # Make a fail state that lets some ages through
      otherwise = 
        data.frame(
          "comparison_var" = test_var,
          "statistical_test" = "two sample t.test", 
          "estimate1" = NA,    
          "estimate2" = NA,    
          "method" = NA,           
          "p.value" = NA,          
          "variance_test" = NA,    
          "bartlett_p" = NA)),
    .id = "age") 
  
  
  # 3. Do minor renaming and adjusting
  ageclass_comparisons <- ageclass_comparisons %>% 
    mutate(age = as.numeric(as.character(age))) %>% 
    rename({{group_1}} := estimate1) %>% 
    rename({{group_2}} := estimate2)
  
  
  # Get the number of values for each age-class in each comparison period
  ages_tested <- both_periods_df %>% 
    filter(age %in% ages_to_test[[comname_x]]) %>% 
    group_by(comparison_period, age) %>% 
    summarise(
      n_obs = n(),
      n_years = n_distinct(est_year),
      .groups = "drop") %>% 
    mutate(
      obs_per_yr = n_obs / n_years,
      age = as.numeric(as.character(age)))
  
  # Reshape it so it can join to the comparison results:
  species_results <- ages_tested %>% 
    pivot_wider(names_from = comparison_period, values_from = c(n_obs, n_years, obs_per_yr)) %>% 
    left_join(ageclass_comparisons, by = "age") %>% 
    arrange(age)
  
  
  # Add what year the earliest measurement came from:
  earliest_measure <- both_periods_df %>% pull(est_year) %>% min()
  species_results <- mutate(species_results, earliest_measure = earliest_measure, .after = "comparison_var")
  
  return(species_results)
}





####  Size at Age Tests  ####



##### a. Length at Age T Tests  ####
laa_1970_start <- map_dfr(
  # Species to loop through
  .x = setNames(vb_species, vb_species), 
  
  # Handle errors/failures
  .f = possibly(
    
    # Function to loop through for species
    .f = ~ species_ageclass_ttests(
      comname_x = .x, 
      group_1 = "1970-2009", 
      group_2 = "2010-2019", 
      test_var = "length_cm"),
    
    # Use otherwise to fill data gaps
    otherwise = data.frame(
      "age" = NA,
      "n_obs_1970-2009" = NA,
      "n_obs_2010-2019" = NA,
      "n_years_1970-2009" = NA,    
      "n_years_2010-2019" = NA,    
      "obs_per_yr_1970-2009" = NA,
      "obs_per_yr_2010-2019" = NA, 
      "comparison_var" = NA,
      "earliest_measure" = NA,
      "statistical_test" = "two sample t.test", 
      "1970-2009" = NA,    
      "2010_2019" = NA,   
      "method" = NA,         
      "p.value" = NA,          
      "variance_test" = NA,    
      "bartlett_p" = NA)), 
  .id = "comname")



##### b. Weight at Age T Tests  ####
waa_1970_start <- map_dfr(
  .x = setNames(vb_species, vb_species), 
  .f = possibly(
    
    .f = ~ species_ageclass_ttests(
      comname_x = .x, 
      group_1 = "1970-2009", 
      group_2 = "2010-2019", 
      test_var = "indwt"),
    
    otherwise = data.frame(
      "age" = NA,
      "n_obs_1970-2009" = NA,
      "n_obs_2010-2019" = NA,
      "n_years_1970-2009" = NA,    
      "n_years_2010-2019" = NA,    
      "obs_per_yr_1970-2009" = NA,
      "obs_per_yr_2010-2019" = NA, 
      "comparison_var" = NA,
      "earliest_measure" = NA,
      "statistical_test" = "two sample t.test", 
      "1970-2009" = NA,    
      "2010_2019" = NA,   
      "method" = NA,         
      "p.value" = NA,          
      "variance_test" = NA,    
      "bartlett_p" = NA)), 
  .id = "comname")








##### c. Percent Change Tidy  ####


# Only using one comparison period for all results now
# make a "late period" for the 2010-2019 also, then rename early & late period means

# Pull the comparisons for 1970-2009 vs. 2010-2019
perc_change_results <- bind_rows(
  list(
    "length-at-age-1970-2009"    = laa_1970_start,  
    "weight-at-age-1970-2009"    = waa_1970_start), 
  .id = "comparison_variable") %>% 
  mutate(early_period_coverage = str_c(earliest_measure, "-2009"),
         late_period_coverage = "2010-2019") %>% 
  select(
    comparison_variable,
    comparison_var,
    early_period_coverage, 
    late_period_coverage, 
    comname,
    age,
    early_period_obs_per_yr = `obs_per_yr_1970-2009`,
    late_period_obs_per_yr = `obs_per_yr_2010-2019`,
    early_period_mean = `1970-2009`,
    late_period_mean = `2010-2019`,
    p.value,
    bartlett_p) %>% 
  mutate(percent_change = ((late_period_mean-early_period_mean)/early_period_mean)*100, 
         .after = "late_period_mean") %>% 
  mutate(across(where(is.numeric), \(x) round(x, 3)))


# Data Reshaping for Barplots of Percent Change

# Set palette colors into the dataframe
palette = c(gmri_cols("light gray"), "#0098B8", "#EACC1B")

# Throw the percent change stuff in the data prep:
# Don't need to nest the data, can just use groups
barplot_prep <- perc_change_results  %>% 
  select(comname, comparison_var, age, percent_change, early_period_mean, late_period_mean, p.value, percent_change) %>% 
  drop_na() %>% 
  group_by(comname, comparison_var) %>% 
  arrange(age) %>% 
  mutate(
    # Make "wins" list column for binary outcome plots
    color = case_when(
      p.value > 0.05 ~ palette[1],
      p.value < 0.05 & early_period_mean > late_period_mean ~ palette[3],
      p.value < 0.05 & early_period_mean < late_period_mean ~ palette[2])) %>% 
  ungroup() 




#### Size at Age in Time  ####

# Get the average ageclass size for all species:
size_at_ages_all <- age_truncated %>% 
  group_by(comname, year = est_year, age) %>% 
  summarise(
    avg_len = mean(length_cm, na.rm =T),
    avg_wt = mean(indwt, na.rm =T),
    avg_condition = mean(fK, na.rm =T),
    .groups = "drop") %>% 
  complete(comname, year, age) %>% 
  mutate(age = factor(age, levels = c(1:10, ">10"))) %>% 
  arrange(comname, year, age)







####  Plotting Functions  ####




# Make a Legend for Size at Age

# Tall
(saa_legend <- tibble(
  label = c("Significant Increase", "No Change", "Significant Decrease"),
  fill  = c(palette[[2]], palette[[1]], palette[[3]]),
  xmin  = c(0.5, 0.5, 0.5),
  xmax = c(1,1,1),
  ymin  = c(1.5, 0.8, 0.1),
  ymax  = c(2, 1.3, 0.6),
  lab_x = c(1.75, 1.75, 1.75),
  lab_y = c(1.75, 1.05, 0.35)) %>% 
   ggplot() +
   geom_rect(aes(xmin =  xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = I(fill))) +
   geom_text(aes(x = xmin, y = lab_y, label = label), hjust = 1.1,  size = 3) +
   #labs(title = "Decadal Variability in Size at Age Characteristics:") +
   scale_x_continuous(expand = expansion(add=c(5,0.5))) +
   scale_y_continuous(expand = expansion(add=c(0.7,0.7))) +
   theme_void() +
   theme(plot.title = element_text(size = 7, face = "bold")))



# Function to plot the Size at Age differences
mean_saa_plots <- function(comname_x, test_var = avg_len, var_label = "Average Length (cm)", is_wt = F){
  
  # Data for lines
  size_at_age_data <- size_at_ages_all %>% 
    filter(comname == comname_x,
           age %in% ages_to_test[[comname_x]]) 
  
  
  # Set a more accurate left-lim for species that were aged starting later
  true_left <- size_at_age_data %>% 
    drop_na({{test_var}}) %>% 
    pull(year) %>% 
    min()
  
  # If weights start all at 1992
  if(is_wt){
    size_at_age_data <- size_at_age_data %>% filter(year >= 1992)
    true_left <- 1992
  }
  
  
  # Data for period averages
  period_list <- list("period_1" = c(1970:2009), "period_3" = c(2010:2019))
  
  # Data to get averages
  species_dat <- age_truncated %>% 
    filter(comname == comname_x, 
           age %in% ages_to_test[[comname_x]])
  
  
  # Get the means over that range, make a dataframe of x, xend, y, yend for segments
  period_means <- imap_dfr(period_list, function(period_years, period_name){
    species_dat %>% 
      filter(est_year %in% period_years) %>% 
      group_by(comname, age) %>% 
      summarise(
        avg_len       = mean(length_cm, na.rm =T),
        avg_wt        = mean(indwt, na.rm =T),
        avg_condition = mean(fK, na.rm =T),
        .groups = "drop") %>% 
      mutate(
        x = ifelse(
          period_name == "period_1", 
          true_left, min(period_years)), 
        xend = max(period_years))
  }, .id = "period") %>% 
    mutate(age_alpha = as.numeric(as.character(age)),
           age_alpha = ifelse(is.na(age_alpha), max(age_alpha, na.rm = T) + 1, age_alpha))
  
  #return(period_means)
  
  # Make the plot with means for periods
  size_at_age_data %>% 
    ggplot(aes(year, {{test_var}})) +
    geom_line(aes(group = age), color = gmri_cols("light gray"), linewidth = 0.5) +
    geom_point(aes(group = age), color = "black", size = 0.5) +
    geom_segment(
      data = period_means,
      aes(x = x, xend = xend, y = {{test_var}}, yend = {{test_var}}, group = age, color = period),
      linewidth = .75, show.legend = FALSE) +
    scale_color_manual(values = c(gmri_cols("blue"), gmri_cols("orange"))) +
    scale_x_continuous(limits = c(true_left, 2019)) +
    labs(
      title = str_to_sentence(comname_x), 
      y = NULL, 
      x = ""
      # y = var_label, 
      # x = "Year"
    )
  
}


# Display the lines for size at age starting at a 0 and extending out

decadal_growth_lines <- function(comname_x, option = "length"){
  
  # Data to get averages
  species_dat <- filter(age_truncated, 
                        comname == comname_x,
                        age %in% ages_to_test[[comname_x]])
  
  #If weight start all at 1992
  if(option == "weight"){
    species_dat <- filter(species_dat, est_year >= 1992)
  }
  
  # Get the averages for each year, give labels for display groups
  age_lines <- species_dat %>% 
    filter(yearclass >= 1980) %>% 
    group_by(comname, birthyear = yearclass,  age) %>% 
    summarise(
      avg_len = mean(length_cm, na.rm =T),
      avg_wt  = mean(indwt, na.rm =T),
      avg_condition = mean(fK, na.rm =T),
      .groups = "drop") %>% 
    mutate(birthdecade =floor_decade(birthyear),
           age = factor(age, levels = c(1:10, ">10")))
  
  #return(age_lines)
  
  # Make decadal averages
  decade_lines <- species_dat %>%
    filter(yearclass>= 1980) %>% 
    mutate(birthdecade = floor_decade(yearclass)) %>% 
    group_by(comname, birthdecade,  age) %>% 
    summarise(
      avg_len = mean(length_cm, na.rm =T),
      avg_wt = mean(indwt, na.rm =T),
      avg_condition = mean(fK, na.rm =T),
      .groups = "drop")
  
  #return(decade_lines)
  
  
  # Average Length
  len_plot <- ggplot(age_lines, aes(age, avg_len)) +
    geom_line(aes(group = birthyear, color = birthdecade), alpha = 0.2, linewidth = 0.5) +
    scale_color_gmri() +
    geom_line(data = decade_lines, aes(group = birthdecade, color = birthdecade), linewidth = 1) +
    geom_dl(data = decade_lines, 
            aes(group = birthdecade, color = birthdecade, label = birthdecade, fontface = "bold"), 
            method = list("last.points", cex = 0.5)) +
    scale_x_discrete(expand = expansion(mult = c(0,.4))) +
    labs(title = str_to_sentence(comname_x), y = "Average Length (cm)", x = "Age", color = "Decade of Birth")
  
  
  # Average Weight
  wt_plot <- ggplot(age_lines, aes(age, avg_wt)) +
    geom_line(aes(group = birthyear, color = birthdecade), alpha = 0.2, linewidth = 0.5) +
    scale_color_gmri() +
    geom_line(data = decade_lines, aes(group = birthdecade, color = birthdecade), linewidth = 1) +
    geom_dl(
      data = decade_lines, 
      aes(group = birthdecade, color = birthdecade, label = birthdecade, fontface = "bold"), 
      method = list("last.points", cex = 0.5)) +
    scale_x_discrete(expand = expansion(mult = c(0,.4))) +
    labs(title = str_to_sentence(comname_x), y = "Average Weight (kg)", x = "Age", color = "Decade of Birth")
  
  
  
  if(option == "length"){return(len_plot)}
  if(option == "weight"){return(wt_plot)}
  
  if(option == "both"){
    both_p <- (len_plot / (wt_plot + labs(title = NULL))) + plot_layout(guides = "collect")
    return(both_p)
  }
}


####  Generate Plots  ####

##### 1. Length at Age Barplot: F5  ####

# Just length
percent_change_laa <- barplot_prep %>% 
  split(.$comname) %>% 
  imap(function(x,y){
    
    # Ages
    len_dat <- x %>% filter(comparison_var == "length_cm")
    max_age <- max(len_dat$age)
    
    # Make length plot
    len_plot <- len_dat %>% 
      ggplot() +
      geom_col(aes(x = age, y = percent_change, fill = I(color))) +
      #geom_line(aes(x = age, y = 0), linewidth = 0.5, color = "darkgray") +
      labs(
        x = "Age", 
        y = "% Change") +
      scale_x_continuous(
        limits = c(0.5, 10.5),
        breaks  = c(1:max_age),
        expand = expansion(add = c(0,0))) +
      scale_y_continuous(
        limits = c(-20, 20),
        breaks = pretty_breaks(),
        labels = scales::label_percent(scale = 1), 
        oob = scales::squish) +
      labs(x = "Age",
           y = NULL,
           title = str_to_sentence(y))
    return(len_plot)
  })


# Add the legend as its own panel
percent_change_laa[[15]] <- saa_legend




# Arrange them all into one page of figures
laa_changes_all <- gridExtra::marrangeGrob(
  percent_change_laa, 
  layout_matrix = matrix(1:16,  nrow = 4, ncol=4, byrow=TRUE), 
  top=NULL,
  left = textGrob(
    expression(bold("Length at Age Change (%)")), rot = 90,
    gp = gpar(col = "black", fontsize = 8))
)


##### 2. Weight at Age Barplot: S8  ####

# Do the steps again for the Weight at Age Plots
percent_change_waa <- barplot_prep %>% 
  split(.$comname) %>% 
  imap(function(x,y){
    
    # Make length plot
    wt_dat <- x %>% filter(comparison_var == "indwt")
    max_age <- max(wt_dat$age)
    wt_plot <- wt_dat  %>% 
      ggplot() +
      geom_col(aes(x = age, y = percent_change, fill = I(color))) +
      #geom_line(aes(x = age, y = 0), linewidth = 0.5, color = "darkgray") +
      labs(
        x = "Age", 
        y = "% Change") +
      scale_x_continuous(
        limits = c(0.5, 10.5),
        breaks  = c(1:max_age),
        expand = expansion(add = c(0,0))) +
      scale_y_continuous(
        limits = c(-50, 50),
        breaks = pretty_breaks(),
        labels = scales::label_percent(scale = 1), 
        oob = scales::squish) +
      labs(x = "Age",
           y = NULL,
           title = str_to_sentence(y))
    
    return(wt_plot)
  })

# Add the legend as its own panel
percent_change_waa[[15]] <- saa_legend


# Arrange them all into one page
waa_changes_all <- gridExtra::marrangeGrob(
  percent_change_waa, 
  layout_matrix = matrix(
    1:16,  
    nrow = 4, 
    ncol=4, 
    byrow=TRUE), 
  top=NULL,
  left = textGrob(
    expression(bold("Weight at Age Change (%)")), 
    rot = 90,
    gp = gpar(col = "black", fontsize = 8)))




##### 3. Size at Age: S9 & S10  ####

#------ Length at age across decades

# Generate Species Specific Length at Age Figures
laa_figs <- map(
  .x = vb_species, 
  .f = mean_saa_plots, 
  test_var = avg_len, 
  var_label = "Individual Length (cm)") %>% 
  setNames(vb_species)


# Do as one page
laa_figs <- gridExtra::marrangeGrob(
  laa_figs, 
  layout_matrix = matrix(1:16,  nrow = 4, ncol=4, byrow=TRUE), 
  top=NULL,
  left = textGrob(
    expression(bold("Individual Length (cm)")), rot = 90,
    gp = gpar(col = "black", fontsize = 8)))





# Generate Species Specific Weight at Age Figures
waa_figs <- map(
  .x = vb_species, 
  .f = mean_saa_plots, 
  test_var = avg_wt,
  var_label = "Individual Weight (kg)", 
  is_wt = T) %>% 
  setNames(vb_species)


# Using arrange
# Do as one page
waa_figs <- gridExtra::marrangeGrob(
  waa_figs, 
  layout_matrix = matrix(1:16,  nrow = 4, ncol=4, byrow=TRUE), 
  top=NULL,
  left = textGrob(
    expression(bold("Individual Weight (kg)")), rot = 90,
    gp = gpar(col = "black", fontsize = 8)))




##### 4. Decadal Growth Change: S11 & S12  ####

#------ length at age decades

# Export as pdf
decadal_laa_figs <- map(
  vb_species, 
  decadal_growth_lines,
  option = "length") %>% 
  setNames(vb_species) %>% 
  map(~.x + labs(y = NULL))

# Make single page
decadal_species_laa <- patchworkGrob(wrap_plots(
  decadal_laa_figs[1:14], ncol = 4, nrow = 4, widths = 3, heights = 3, guides = "collect"))

# Use grid.arrange to make a shared y axis label
decadal_species_laa <- gridExtra::grid.arrange(
  decadal_species_laa, 
  left = textGrob(
    expression(bold("Average Length (cm)")), rot = 90,
    gp = gpar(col = "black", fontsize = 8))
)


#------ Weight at age decades


# Export them as one figure
decadal_waa_figs <- map(
  .x = vb_species,
  .f =  decadal_growth_lines,
  option = "weight") %>% 
  setNames(vb_species) %>% 
  map(~.x + labs(y = NULL))


# Assemble them into one plot
decadal_species_waa <- patchworkGrob(wrap_plots(
  decadal_waa_figs[1:14], ncol = 4, nrow = 4, widths = 3, heights = 3, guides = "collect"))


# Add a shared y axis
decadal_species_waa <- gridExtra::grid.arrange(
  decadal_species_waa, 
  left = textGrob(
    expression(bold("Average Weight (kg)")), rot = 90,
    gp = gpar(col = "black", fontsize = 8))
)




####  Saving Figures  #####

# Height per Row: 45mm
# Width for 2 columns: 170mm


# % Change Barplots - length at age F5
ggsave(
  str_c(decadal_hires_folder, "Figure_5.pdf"),
  laa_changes_all, 
  height = 180, width = 170, units ="mm", 
  dpi = 500)


# % Change Barplots - Weight at Age S8
ggsave(
  str_c(decadal_hires_folder, "Figure_S8.pdf"),
  waa_changes_all, 
  height = 180, width = 170, units ="mm", 
  dpi = 500)


# Timeseries - Length at Age S9
ggsave(
  filename = str_c(decadal_hires_folder, "Figure_S9.pdf"),
  laa_figs, 
  height = 180, width = 170, units ="mm", 
  dpi = 500)

# Timeseries - Weight at Age S10
ggsave(
  filename = str_c(decadal_hires_folder, "Figure_S10.pdf"),
  waa_figs, 
  height = 180, width = 170, units ="mm", 
  dpi = 500)


# SAA Lines - Length S11
ggsave(
  filename = str_c(decadal_hires_folder, "Figure_S11.pdf"), 
  decadal_species_laa, 
  height = 180, width = 170, units ="mm", 
  dpi = 500)


# SAA Lines - Weight S12
ggsave(
  filename = str_c(decadal_hires_folder, "Figure_S12.pdf"), 
  decadal_species_waa, 
  height = 180, width = 170, units ="mm", 
  dpi = 500)





##### a. Export Size at Age Tables  ####
# These have been exported as a csv file

# Reviewers wanted a table - Length at Ages
laa_table_byyr <- size_at_ages_all %>% 
  drop_na(avg_len) %>% 
  select(-c(avg_condition, avg_wt)) %>% 
  pivot_wider(names_from = "age", values_from = "avg_len") %>% 
  relocate(`>10`, .after = last_col())

# Reviewers wanted a table - Weight at Ages
waa_table_byyr <- size_at_ages_all %>% 
  drop_na(avg_wt) %>% 
  select(-c(avg_condition, avg_len)) %>% 
  pivot_wider(names_from = "age", values_from = "avg_wt") %>% 
  relocate(`>10`, .after = last_col())

# Relative Path to Files
data_save_path <- here::here("data/decadal_variability_data/")

write_csv(laa_table_byyr, str_c(data_save_path, "Annual_mean_length_at_age.csv"))
write_csv(waa_table_byyr, str_c(data_save_path, "Annual_mean_weight_at_age.csv"))
