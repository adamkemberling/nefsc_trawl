#### Size at Age Testing
#### 11/11/2021



####  Packages  ####
library(FSAdata)
library(FSA)
library(car)
library(tidyverse)
library(gmRi)
library(ggridges)


# set theme
theme_set(theme_minimal() + 
            theme(legend.position = "bottom", 
                  plot.title = element_text(hjust = 0.5),
                  axis.line.x = element_line(),
                  axis.ticks.x = element_line(), 
                  axis.text = element_text(size = 11),
                  strip.text = element_text(color = "white", 
                                            face = "bold",
                                            size = 11),
                  strip.background = element_rect(
                    color = "white", 
                    fill = "#36454F", 
                    size = 1, 
                    linetype="solid"), 
                  legend.background = element_rect(fill = "transparent", 
                                                   color = "black"))
          )



#####  Load Data  ####

# Access Biological Data with {gmRi}
nmfs_bio <- gmRi::gmri_survdat_prep(survdat_source = "bio")

# Pull Cod
cod <- filter(nmfs_bio, comname == "atlantic cod")




####  EDA  ####
cod %>% ggplot(aes(age, length_cm)) +
  geom_point(alpha = 0.1, color = "gray60") +
  facet_grid(season~svvessel)


# pull relevant columns for length-age and weight-age relationship
cod_prep <- cod %>% 
  select(year = est_year, tow_id = id, est_towdate, survey_area, comname, 
         indid, length_cm, indwt, age, sex, maturity)


#### VBert Example  ####

# Check the model structure/code
( vb <- vbFuns(param="Typical") )

# Get reasonable starting points:
cod_starts <- vbStarts(length_cm ~ age, data = cod_prep)


# Use nls to estimate VBGF parameters using starting points
cod_fit <- nls(length_cm ~ vb(age, Linf, K, t0), data = cod_prep, start = cod_starts)

# Access parameters with coef()
coef(cod_fit)


# Use boot() to get bootstrapped confidence intervals
# cod_boot1 <- Boot(cod_fit)  # Be patient! Be aware of some non-convergence
# confint(cod_boot1)



####  Plotting Predicted Values  ####

# Get prediction and confidence intervals over range of ages
ages <- seq(0, 16, by = 0.2)
preds1 <- data.frame(ages,
                     predict(cod_fit, data.frame(age = ages))#,
                     #confint(cod_boot1)
                     )
names(preds1) <- c("age","fit"#,"LCI","UCI"
                   )



# Restrict predictions to only ages found in the data
agesum <- group_by(cod_prep, sex) %>%
  summarize(minage = min(age, na.rm = T),
            maxage = max(age, na.rm = T))

# filter
preds2 <- filter(preds1, 
                 age >= min(agesum$minage), 
                 age <= max(agesum$maxage))


# Plot the overall age-length relationship
ggplot() + 
  geom_point(data = filter(cod_prep, year == 2019), aes(y = length_cm, x = age), size=2, alpha=0.1) +
  geom_line(data=preds2, aes(x=age, y = fit), size = 1.5, linetype = 1, color = gmri_cols("orange")) +
  labs(x = "age", y = "Total Length (cm)", subtitle = "Atlantic Cod - All Data")


# Coefficients
coef(cod_fit)

# Get Sizes at age as a table
preds2 %>% filter(age %in% c(0:15))
preds2 %>% filter(age %% 1 == 0) # integers only




#####___________########
####  Von Bert Workflow  ####


####  Data Prep  #### 
# Setting up for group analyses

# Rank species by how many measurements there are
species_abunds <- nmfs_bio %>% count(comname) %>% arrange(desc(n)) # ordered by number measured

# species to run it on:
vonbert_species <- sort(c("atlantic cod", "haddock", "summer flounder", "winter flounder", "acadian redfish",
                     "silver hake", "red hake", "american plaice", "yellowtail flounder",
                     "atlantic herring", "white hake", "butterfish","witch flounder", 
                     "atlantic mackerel", "windowpane"))
names(vonbert_species) <- vonbert_species

seasons <- c("Summer", "Fall")


# Does not work for:
# skates, dogfish, fourspot flounder, goosefish


#### 5 year increments



##### pull species and make group id for later
species_data <- map(vonbert_species, function(species_name){
  
  yr5_breaks <- seq(1970, 2020, by = 5)
  yr5_ends <- seq(1974, 2024, by = 5)
  yr5_labels <- str_c(yr5_breaks, "-", yr5_ends)
  yr5_labels <- yr5_labels[1:(length(yr5_labels)-1)]
  
  nmfs_bio %>%
    filter(comname == species_name,
           season %in% seasons,
           is.na(age) == FALSE) %>% 
    mutate(five_yr_group = cut(est_year, 
                         breaks = yr5_breaks, 
                         labels = yr5_labels,
                         include.lowest = TRUE,
                         ordered_result = TRUE),
           decade = floor_decade(est_year),
           group_id = str_c(est_year, comname, sep = "-"))  %>% 
    as.data.frame()
})






####  Processing Von Bert Coefficients  ####


# Function to Pull Vonbert Coefficients
vonbert_coef <- function(length_age_dat, start_points, min_obs = 20){
  
  # Von Bert Fitting Using: {FSA}
  # Don't run on fewer than a minimum number of obervations
  num_aged <- nrow(length_age_dat)
  if(num_aged < min_obs){
    na_df <- data.frame("Linf" = NA, "K" = NA, "t0" = NA, "n_aged" = num_aged)
    return(na_df)
  }
  
  # Use nls to estimate VBGF parameters using starting points
  vbert_fit <- nls(length_cm ~ vb(age, Linf, K, t0), 
                   data = length_age_dat, 
                   start = start_points)
  
  # Access parameters with coef()
  vbert_coef <- as.data.frame(t(coef(vbert_fit))) %>% 
    mutate(n_aged = num_aged)
  return(vbert_coef)
  
  
}



# Function for running that for a  single Species
testing_species <- function(comname, split_col, min_obs){
  
  # Get starting points from all data
  test_starts <-  vbStarts(length_cm ~ age, data = species_data[[comname]])
  
  # Get coefficients for groups
  species_data[[comname]] %>% 
    split(.[split_col]) %>% 
    map_dfr(vonbert_coef, test_starts, min_obs = min_obs, .id = split_col) %>% 
    mutate(comname = comname)
}


# Just run through them individually
testing_species(vonbert_species[1], "five_yr_group", min_obs = 25)
testing_species(vonbert_species[2], "five_yr_group", min_obs = 25)
testing_species(vonbert_species[3], "five_yr_group", min_obs = 25)
testing_species(vonbert_species[4], "five_yr_group", min_obs = 25)
testing_species(vonbert_species[5], "five_yr_group", min_obs = 25)
testing_species(vonbert_species[6], "five_yr_group", min_obs = 25)
testing_species(vonbert_species[7], "five_yr_group", min_obs = 25)
testing_species(vonbert_species[8], "five_yr_group", min_obs = 25) 
testing_species(vonbert_species[9], "five_yr_group", min_obs = 25)
testing_species(vonbert_species[10], "five_yr_group", min_obs = 25)
testing_species(vonbert_species[11], "five_yr_group", min_obs = 25)
testing_species(vonbert_species[12], "five_yr_group", min_obs = 25)
testing_species(vonbert_species[13], "five_yr_group", min_obs = 25)
testing_species(vonbert_species[14], "five_yr_group", min_obs = 25)
testing_species(vonbert_species[15], "five_yr_group", min_obs = 25)




# Running everything?
species_coef <- vonbert_species %>% 
  map_dfr(testing_species, split_col =  "five_yr_group", min_obs = 25)


####  Plotting Coefficients  ####

# pick species and coefficients
plot_vonbert_coef <- function(comnames, x_col, coef_name){
  
  x_col_sym <- sym(x_col)
  coef_sym <- sym(coef_name)
  coef_label <- switch (
    coef_name,
    Linf = "Asymptotic Max Length (cm)",
    K = "Body Growth Coefficient (K)",
    t0 = "t0",
    n_aged = "Number Aged"
  )
  species_coef %>% 
    filter(comname %in% comnames) %>% 
    ggplot(aes(y = {{coef_sym}}, x = {{x_col_sym}})) +
    geom_line(group = 1, linetype = 3) +
    geom_point(size = 0.6) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    facet_wrap(~comname, ncol = 1, scales = "free") +
    labs(x = "Date Period", y = coef_label)
  
}



#  Test the plots
plot_vonbert_coef(comnames = c("atlantic cod", "silver hake", "red hake"), 
                  x_col = "five_yr_group", 
                  coef_name = "K")

plot_vonbert_coef(comnames = c("atlantic cod", "silver hake", "red hake"), 
                  x_col = "five_yr_group", 
                  coef_name = "n_aged")

plot_vonbert_coef(comnames = c("atlantic cod",  "haddock", "red hake"), 
                  x_col = "five_yr_group", 
                  coef_name = "Linf")

plot_vonbert_coef(comnames = c("atlantic cod", "silver hake"), 
                  x_col = "five_yr_group", 
                  coef_name = "t0")



####  Plotting Age Cohorts  ####

# Plotting age distribution
plot_age_bubbleplot <- function(species){
  species_data[[species]] %>% 
    count(comname, est_year, age) %>% 
    ggplot(aes(y = est_year, x = age, size = n)) +
      geom_point(shape = 21, color = "white", fill = gmri_cols("gmri blue")) +
      scale_y_reverse() +
      facet_wrap(~comname) +
      labs(y = "Year", x = "Age", size = "Number Measured") +
      guides(size = guide_legend(title.position = "top", title.hjust = 0.5))
  }


# # Plot them
bubble_plots <- map(.x = vonbert_species, .f = plot_age_bubbleplot)
bubble_plots$`atlantic cod`
bubble_plots$`silver hake`
bubble_plots$`acadian redfish`
bubble_plots$haddock


# Plotting age distribution as ggridges
plot_age_ridgeplot <- function(species){
  species_data[[species]] %>% 
    mutate(yr = factor(est_year),
           yr = fct_rev(yr)) %>% 
    ggplot(aes(x = age, y = yr)) +
      geom_density_ridges(fill = gmri_cols("gmri blue"), color = "white") +
      facet_wrap(~comname) +
      labs(x = "Age", y = "Year") +
      guides(size = guide_legend(title.position = "top", title.hjust = 0.5))
  
}

# # Plot
# plot_age_ridgeplot("atlantic cod")
# plot_age_ridgeplot("haddock")
ridge_plots <- map(vonbert_species, plot_age_ridgeplot)
ridge_plots$`atlantic cod`
ridge_plots$`silver hake`
ridge_plots$`acadian redfish`
ridge_plots$haddock

