remotes::install_github("NOAA-EDAB/ecodata")

## Loads EPU as sf objects
epu_sf <- ecodata::epu_sf
sf::plot_sf(epu_sf)



## Run VAST on a NE EPU, here, Gulf of Maine
epu <- "GOM"
settings <- FishStatsUtils::make_settings(n_x = 50,
                                          Region = "northwest_atlantic",
                                          purpose = "index2",
                                          strata.limits = "EPU")
## Add the appropriate "EPU" to the settings
settings$epu_to_use <- switch(epu,
                              "MAB" = "Mid_Atlantic_Bight",
                              "GOM" = "Gulf_of_Maine",
                              "GB" = "Georges_Bank",
                              "SS" = "Scotian_Shelf",
                              "All" = "All")

# Run model
fit <- FishStatsUtils::fit_model(settings = settings,
                                 ## add data here
                                 epu_to_use = settings$epu_to_use) # Must be included here, or else not exported to extrapolation_args_*


####  Accessing Ecodata  ####

agg_bio <- ecodata::aggregate_biomass
agg_bio %>% 
  filter(Source == "NEFSC bottom trawl survey (survdat)",
         Var == "Apex Predator Fall Biomass Index") %>% 
  ggplot(aes(Time, Value)) +
    geom_line(show.legend = FALSE) +
    facet_wrap(~EPU) 



####  Aggregate Biomass plots

# Source Code here: https://github.com/NOAA-EDAB/ecodata/blob/master/chunk-scripts/macrofauna_NE.Rmd-agg-bio.R

