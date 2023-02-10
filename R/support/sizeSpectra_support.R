####_____________________####
####  sizeSpectra Support Functions  ####


#  Load Libraries
library(ggpmisc)
library(scales)
library(patchwork)
library(sizeSpectra)
library(tidyverse)



####_____________________####
####  Prep NEFSC Data  ####




# Build dataOrig for sizeSpectra
#' @title Prep Data for sizeSpectra - convert to grams
#' 
#' @description Renames fields and prepares things in grams to use with sizeSpectra function. I hate 
#' how this is written, so confusing with so many columns that are too similar.
#' 
#' 
#' 
#' Re-work 07/21/2021
#' 
#' Changes from kg to grams for wmin_g and wmax
#'
#' @param lw_trawl_data data that contains length weight relationships and numbers caught
#'
#' @return
#' @export
#'
#' @examples
prep_sizeSpectra_data <- function(lw_trawl_data){
  
  # Filter out empty weights
  ss_dat <- lw_trawl_data  %>% 
    filter(is.na(ind_weight_kg) == FALSE)
  
  
  # Format columns for consistency downstream
  ss_dat <- ss_dat %>%
    rename(
      Year         = est_year,                  # year
      # SpecCode   = comname,                   # common name
      # LngtClass  = length_cm,                 # length bin
      # Number     = numlen_adj,                # adjusted number at length
      # LWa        = a,                         # length/weight parameter
      # LWb        = b                          # length/weight parameter
    )
  
  # Get bodymass in g, and the biomass from the catch frequencies
  ss_dat <- ss_dat %>%                     
    mutate(
      ind_weight_g   = ind_weight_kg * 1000,              # weight in grams of individual
      sum_weight_g   = numlen_adj * ind_weight_g,                 # number at size times the mass
      lw_group       = str_c(comname, season, catchsex) # the source of l-w coefficients
    ) 
  
  
  # Get Upper Length Bin Weights
  # Assumes that for all species the lengths only jump in 1 cm increments
  # Captures min/max predicted weight across that length increment
  ss_dat <- ss_dat %>%
    mutate(
      lngt_max  = length_cm + 1,                # maximum length for a fish of wmin
      ln_wmax  = (ln_a + b * log(lngt_max)),  # max weight for fish of that wmin
      wmin_g    = ind_weight_g,                     # minimum weight of individual in grams
      wmax_g    = exp(ln_wmax) * 1000           # max weight of individual in grams
    )                
  
  
  # If available do stratified abundance as well
  if("strat_total_abund_s" %in% names(ss_dat)){
    ss_dat <- ss_dat %>% 
      mutate(
        stratified_sum_bodymass = strat_total_abund_s * ind_weight_g,       
        wmin_area_strat         = strat_total_abund_s * wmin_g,           
        wmax_area_strat         = strat_total_abund_s * wmax_g           
      )
  }
  
  # Arrange
  ss_dat <- ss_dat %>% 
    arrange(Year, season, comname, length_cm)
  
  # Filter zero and NA weights
  ss_dat <- filter(ss_dat,
                   wmin_g != 0,
                   is.na(wmin_g) == FALSE,
                   wmax_g!= 0,
                   is.na(wmax_g) == FALSE)
  
  
  
  
  # Rename the functional groups
  ss_dat <- ss_dat %>% 
    mutate(
      spec_class = case_when(
        spec_class == "gf"   ~ "Groundfish",
        spec_class == "pel"  ~ "Pelagic",
        spec_class == "dem"  ~ "Demersal",
        spec_class == "el"   ~ "Elasmobranch",
        spec_class == "<NA>" ~ "NA"))
  
  # Make regions go N->S
  ss_dat <- ss_dat %>% 
    mutate(
      survey_area = factor(survey_area, levels = c(
        "GoM", "GB", "SNE", "MAB")),
      season = factor(season, levels = c(
        "Spring", "Fall")))
  
  
  # change texxt for fishery status
  ss_dat <- ss_dat %>% 
    mutate(fishery = case_when(
      fishery == "com" ~ "Commercially Targeted",
      fishery == "nc"  ~ "Not Commercially Targeted",
      TRUE             ~ "Not Labelled"))
  
  # return data prepped for sizeSpectra
  return(ss_dat)
  
}



#### Assign Functional Groups:



#' @title Assign Remaining Functional Groups
#'
#' @description Supplement the GMRI cleanup step's species name and functional group assignments to 
#' pick up remaining stragglers.
#'
#' @param species_dat Dataframe containing common name "comname", species class "hare_group",
#' and functional group "hare_group" with which to supplement.
#'
#' @return
#' @export
#'
#' @examples
fill_func_groups <- function(species_dat){
  
  # Assign the Rest of the Unclassified Values
  spec_class_swap <- species_dat %>% 
    distinct(comname, spec_class, hare_group) %>% 
    mutate(
      spec_class = case_when(
        str_detect(comname, "ray|shark|sand tiger|skate")      ~ "Elasmobranch",
        str_detect(comname, "kingfish|sturgeon|weakfish|spot") ~ "Demersal",
        str_detect(comname, "kingfish|sturgeon|weakfish|spot") ~ "Demersal",
        str_detect(comname, "mackerel|herring|sardine")        ~ "Pelagic",
        str_detect(comname, "amberjack|spadefish")             ~ "Reef",
        TRUE ~ spec_class),
      hare_group = case_when(
        str_detect(comname, "flounder")       ~ "groundfish",
        str_detect(comname, "scup")           ~ "coastal",
        str_detect(comname, "thread herring") ~ "coastal",
        str_detect(comname, "dory")           ~ "pelagic",
        str_detect(comname, "sturgeon")       ~ "diadromous",
        TRUE ~ hare_group
      )
    )
  
  # Join the new ones in:
  species_complete <- species_dat %>% 
    select(-spec_class, -hare_group) %>% 
    left_join(spec_class_swap, by = "comname")
  
  return(species_complete)
  
  # # For reference, Kathy's group's gaps
  # spec_class_gaps <- species_dat %>% 
  #   distinct(comname, spec_class) %>% 
  #   filter(is.na(spec_class))
  # 
  # # For Reference:
  # # Gaps in hare paper's functional groups
  # hare_group_gaps <- species_dat %>% 
  #   distinct(comname, hare_group)  %>% 
  #   filter(is.na(hare_group))
  
}




####  Filter Min Length

# data for testing:
# library(targets); tar_load("catch_stratified")


#' @title Set minimum-size cutoff for size spectrum analyses
#'
#' @param trawl_lens The input dataset to filter, should contain length in cm as "length"
#' @param cutoff_cm Length in cm to use as filter
#'
#' @return
#' @export
#'
#' @examples
min_length_cutoff <- function(trawl_lens, cutoff_cm = 1){
  dplyr::filter(trawl_lens,
                is.na(length_cm) == FALSE, 
                length_cm >= cutoff_cm) }


#' @title Set minimum-weight cutoff for size spectrum analyses
#'
#' @description Minimum weight for size spectrum analysis should
#' reflect the minimum size to be adequately caught by the sampling
#' gear.
#' 
#' @param catch_lw The input dataset to filter, should contain length in cm as "length"
#' @param min_weight_g Weight in grams to use as filter
#'
#' @return
#' @export
#'
#' @examples
min_weight_cutoff <- function(catch_lw, min_weight_g = 1){
  dplyr::filter(catch_lw,
                is.na(length_cm) == FALSE, 
                ind_weight_kg >= (min_weight_g/1000))  }




#' @title Set maximum-weight cutoff for size spectrum analyses
#'
#' @description Maximum weight for size spectrum analysis should
#' reflect the maximum size to be adequately caught by the sampling
#' gear.
#'
#' @param catch_lw The input dataset to filter, should contain length in cm as "length"
#' @param min_weight_g Weight in grams to use as filter
#'
#' @return
#' @export
#'
#' @examples
max_weight_cutoff <- function(catch_lw, max_weight_g = 10^5){
  dplyr::filter(catch_lw,
                ind_weight_kg < (max_weight_g/1000))  }






#' @title Add Size Group Formatting of Plotting
#'
#' @description Sets up the text labels and formatting for grouping and 
#' plotting groups by size, weight, or season. Size groups do not follow log10 increments of 
#' size spectrum stuff
#'
#' @param catch_1g 
#'
#' @return
#' @export
#'
#' @examples
size_bin_formatting <- function(catch_1g){
  
  # Cut up discrete length and weight bins
  catch_size_bins <- catch_1g %>% 
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
  catch_size_bins <- catch_size_bins %>% 
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
  
  
  
  
return(catch_size_bins)
  
  
  }




# # testing
# prep_sizeSpectra_data(lw_trawl_data = catch_stratified)







####______________________####
####  Spectrum Group Configuration  ####



#' @title Add Missing Group Columns that drop with group_by or Split
#' 
#' @description Adds factor columns back in to the table that were not expressly
#' included as a grouping column. These columns are added in with values of "all" to
#' communicate that all data across that factor are included.
#'
#' @param group_dataframe Table of distinct factor combinations used to provide summary details.
#' This table is checked to see what columns are missing, which then will be added here.
#'
#' @return
#' @export
#'
#' @examples
add_missing_groups <- function(group_dataframe){
  
  # 4. Make flags for group levels not specified
  # for groups not listed, all levels are used
  # ex. if year not used, then all years were used together
  yr_flag     <- "Year" %in% names(group_dataframe)
  szn_flag    <- "season" %in% names(group_dataframe)
  area_flag   <- "survey_area" %in% names(group_dataframe)
  decade_flag <- "decade" %in% names(group_dataframe)
  class_flag  <- "hare_group" %in% names(group_dataframe)
  
  # flag for which ones to build
  flag_checks <- c("Year"        = yr_flag, 
                   "season"      = szn_flag, 
                   "survey_area" = area_flag, 
                   "decade"      = decade_flag,
                   "hare_group"  = class_flag)
  
  # Columns to add as "all"
  cols_2_add <- flag_checks[which(flag_checks == FALSE)] 
  
  # Add the columns that don't exist to grouping table
  # the data value is "all" as all group levels are included
  for (col_flag in seq_along(cols_2_add)) {
    col_name <- names(cols_2_add[col_flag])
    group_dataframe[, col_name] <- "all"}
  
  # Return the table with the new columns
  return(group_dataframe)
  
}




####_____________________####
#### MLE ISD Spectra Analysis ####




# Takes minimum weight bin data split by a grouping variable .$group_var
# returns the power law coefficients for that group

#' @title Process maximum likelihood size-spectra slope for group of catch data
#'
#' @param dataBinForLike Catch data prepared for mle calculation, use prep_wmin_wmax
#' @param group_var Name of the grouping variable, useful for maintaining source grouping
#' @param abundance_vals Option to use "observed" catch or "stratified" catch 
#' @param vecDiff Optional tuning parameter for mle estimation, sometimes needed when n is small
#'
#' @return
#' @export
#'
#' @examples
group_isd_calc <- function(dataBinForLike, 
                           group_var, 
                           abundance_vals = "stratified",
                           isd_xmin = NULL,
                           isd_xmax = NULL,
                           vecDiff = 0.5){
  
  # 1. toggle which abundance calue to use with switch
  abundance_col <- switch(abundance_vals,
                          "observed"   = sym("numlen_adj"),
                          "strat_mean" = sym("strat_mean_abund_s"),
                          "stratified" = sym("strat_total_abund_s"))
  
  
  # 2. Select only the columns we need:
  # Rename them to match the code for sizeSpectra
  dataBinForLike <- dataBinForLike %>% 
    dplyr::select(
      SpecCode = comname,
      wmin = wmin_g,
      wmax = wmax_g,
      Number = !!abundance_col)
  
  
  
  # 3. Set the power-law limits:
  # n, xmin, xmax
  
  # Total individuals
  n    <- sum(dataBinForLike$Number)
  
  # left bound, grams
  if(is.null(isd_xmin)){
    isd_xmin <- min(dataBinForLike$wmin, na.rm = T)
  }
  
  # right bound, grams
  if(is.null(isd_xmax)){
    isd_xmax <- max(dataBinForLike$wmax, na.rm = T)
  }
  
  
  
  # 4. Estimate the Maximum likelihood calculation for the Individual size distribution
  # previously named sizeSpectra::MLEbins.nSeaFung.new
  mle_group_bins <- calcLike(
    negLL.fn          = negLL.PLB.bins.species,
    p                 = -1.9,
    vecDiff           = vecDiff,
    suppress.warnings = TRUE,
    dataBinForLike    = dataBinForLike,
    n                 = n,
    xmin              = isd_xmin,
    xmax              = isd_xmax)
  
  
  # 5. Store outputs in a dataframe
  mle_group_results <- data.frame(
    group_var = group_var,
    xmin      = isd_xmin,
    xmax      = isd_xmax,
    n         = n,
    b         = mle_group_bins$MLE,
    confMin   = mle_group_bins$conf[1],
    confMax   = mle_group_bins$conf[2]) 
  
  # Process C and standard error
  mle_group_results <- mle_group_results %>% 
    mutate(stdErr = (abs(confMin - b) + abs(confMax - b)) / (2 * 1.96),
           C = (b != -1 ) * (b + 1) / ( xmax^(b + 1) - xmin^(b + 1) ) + (b == -1) * 1 / ( log(xmax) - log(xmin)))
  
  
  
  # Spit out the final output
  return(mle_group_results)
}






#' @title  Plot Individual Size Distribution Estimates for a group of estimates
#'
#' @description Plot the Individual Size Distribution fits from group_isd_calc displaying the 
#' exponent of individual size spectra (b) and confidence intervals on its estimates.
#'
#' @param mle_res output from group_isd_calc
#'
#' @return
#' @export
#'
#' @examples
group_isd_plot <- function(mle_res){
  plot <- mle_res %>% 
    ggplot(aes(Year, b, color = area, shape = season)) +
    geom_pointrange(aes(x = Year, y = b, ymin = confMin, ymax = confMax),
                    alpha = 0.6) +
    labs(x = "",
         y = "Exponent of Size Spectra (b)") +
    theme(axis.text.x = element_text(angle = 90, size = 6, vjust = 0.5))
  
  return(plot)
}









#' @title Grouped MLE of of Individual Size Distribution Fits
#' 
#' @description Estimates ISD exponent and confidence errors for specific factor groups
#' and returns any un-used factors as "all" to be explicit about data was used.
#'
#' @param wmin_grams Data prepped for mle estimation, use prep_wmin_wmax
#' @param min_weight_g Minimum weight cutoff to use for filtering input data
#' @param max_weight_g Maximum weight cutoff to use for filtering input data
#' @param isd_xmin minimum weight in grams for ISD power law fitting
#' @param isd_xmax maximum weight in grams for ISD power law fitting
#' @param abundance_vals Flag to use "observed" or "stratified" abundances
#' @param .group_cols Vector of strings to indicate colummns to group on
#'
#' @return
#' @export
#'
#' @examples
group_isd_estimation <- function(wmin_grams, 
                                     min_weight_g = 1, 
                                     max_weight_g = 10^6,
                                     isd_xmin = NULL,
                                     isd_xmax = NULL,
                                     abundance_vals = "stratified",
                                     .group_cols = "Year"){
  
  # 1. Set bodymass lower limit
  # Used to filter for upper & lower end of gear selectivity
  dbin_truncated <- filter(wmin_grams, 
                           wmin_g >= min_weight_g,
                           wmin_g < max_weight_g) 
  
  
  # Set Bounds to Individual Size Distribution 
  if(is.null(isd_xmin)){
    isd_xmin <- min(dbin_truncated$wmin_g, na.rm = T) 
  }
  if(is.null(isd_xmin)){
   isd_xmax <- max(dbin_truncated$wmin_g, na.rm =T) 
  }
  
  
  
  # 3. Build group_level from desired group columns
  dbin_truncated <- dbin_truncated %>% 
    unite(col = "group_var", {{.group_cols}}, sep = "_", remove = FALSE, na.rm = FALSE)
  
  # 4. Make a table of constituent combinations
  col_syms <- syms(.group_cols)
  grouping_table <- dbin_truncated %>% distinct(!!!col_syms, group_var)
  
  # 5. Add in any missing groups as "all data"
  grouping_table <- add_missing_groups(group_dataframe = grouping_table)
  
  # 6. Run individual size distribution fits for each group
  group_results <- dbin_truncated %>% 
    split(.$group_var) %>% 
    imap_dfr(.f = group_isd_calc, 
             isd_xmin = isd_xmin,
             isd_xmax = isd_xmax,
             abundance_vals = abundance_vals)
  
  # Merge in the group details
  group_results <- full_join(grouping_table, group_results, by = "group_var") %>% 
    mutate(across(c(group_var, Year, season, survey_area, decade), as.character))
  
  # Return the results
  return(group_results)
  

  
}


# Run all the groupings we care about to build single table

#' @title Perform Estimates of Individual Size Distribution Fits for All WARMEM Groups
#' 
#' @description Source script: 06_nmfs_ss_group_estimates.R, original code added to phase 
#' out section. That code was made to function more modular with group_isd_estimation
#'
#' @param wmin_grams Data prepped for mle estimation, use prep_wmin_wmax
#' @param min_weight_g Minimum weight cutoff to usefor filtering input data
#' @param max_weight_g Maximum weight cutoff to use for filtering input data
#' @param isd_xmin minimum weight in grams for ISD power law fitting
#' @param isd_xmax maximum weight in grams for ISD power law fitting
#' @param abundance_vals Flag to use "observed" or "stratified" abundances
#'
#' @return
#' @export
#'
#' @examples
warmem_isd_estimates <- function(wmin_grams, 
                                 min_weight_g = 1, 
                                 max_weight_g = 10^6,
                                 isd_xmin = NULL,
                                 isd_xmax = NULL,
                                 abundance_vals = "stratified"){
  
  
  ####__  Set Bodymass Cutoff and Groups
  
  # 1. Set bodymass lower limit
  # Used to filter for lower end of gear selectivity
  dbin_truncated <- filter(wmin_grams, 
                           wmin_g >= min_weight_g,
                           wmin_g <= max_weight_g) %>% 
    mutate(decade = floor_decade(Year))
  
  
  # 2. Set up the factor groupings we want to compare : 
  
  
  ####_ 5. Every year, entire survey 
  message("Calculating ISD exponent for each year")
  g5_res <- dbin_truncated  %>% 
    group_isd_estimation(min_weight_g = min_weight_g, 
                         max_weight_g = max_weight_g,
                         isd_xmin = isd_xmin,
                         isd_xmax = isd_xmax,
                         abundance_vals = abundance_vals,
                         .group_cols = c("Year")) 
  
  ####_ 6. every year, every region 
  message("Calculating ISD exponent for each year in each region")
  g6_res <- dbin_truncated  %>% 
    group_isd_estimation(min_weight_g = min_weight_g, 
                         max_weight_g = max_weight_g,
                         isd_xmin = isd_xmin,
                         isd_xmax = isd_xmax,
                         abundance_vals = abundance_vals,
                         .group_cols = c("Year", "survey_area")) 
  
  ####_ 7. every year, only seasons 
  message("Calculating ISD exponent for each year in each season")
  g7_res <- dbin_truncated %>% 
    group_isd_estimation(min_weight_g = min_weight_g, 
                         max_weight_g = max_weight_g,
                         isd_xmin = isd_xmin,
                         isd_xmax = isd_xmax,
                         abundance_vals = abundance_vals,
                         .group_cols = c("Year", "season"))
  
  
  ####_ 8. every year, region * season 
  message("Calculating ISD exponent for each year in each region, for every season")
  g8_res <- dbin_truncated %>% 
    group_isd_estimation(min_weight_g = min_weight_g, 
                         max_weight_g = max_weight_g,
                         isd_xmin = isd_xmin,
                         isd_xmax = isd_xmax,
                         abundance_vals = abundance_vals,
                         .group_cols = c("Year", "season", "survey_area")) 
    
  
  # Put the reults in one table with an ID for how they groups are set up
  table_complete <- bind_rows(
    list(
      "single years"                   = g5_res,
      "single years * region"          = g6_res,
      "single years * season "         = g7_res,
      "single years * season * region" = g8_res), 
    .id = "group ID")
  
  # Return the summary table
  return(table_complete)
  
}






####____________________________####
#####  Plotting ISD  ####


# Function to prepare data for Individual Size Distribution Plots
# Takes the wmin_g and wmax_g data split up by the grouping variable
# Gets the number of fishes within each of the size class bins
# Returns the table, wmin_g, wmax_g, abundance,
# and totals for number of fish that fill in, out, or between size bins


#' @title Individual Size Distribution Plot Prep
#' 
#' @description Prepares wmin_grams data to be plotted with the MLE
#' size spectrum slope fit.
#'
#' @param biomass_data Data prepped for mle estimation, use prep_wmin_wmax
#' @param stratified_abundance Flag to use observed or stratified abundances
#' @param min_weight_g Minimum weight cutoff to use to match slope estimation
#'
#' @return
#' @export
#'
#' @examples
isd_plot_prep <- function(biomass_data = databin_split, 
                          stratified_abundance = FALSE,
                          min_weight_g = 1){
  
  # Toggle abundance to stratified abundances
  if(stratified_abundance == TRUE){
    biomass_data <- biomass_data %>% 
      rename(survey_abund = numlen_adj) %>% 
      rename(numlen_adj = strat_total_abund_s) }
  
  # arrange by lower weight 
  biomass_data <- dplyr::arrange(biomass_data, desc(wmin_g)) %>% 
    select(wmin_g, wmax_g, numlen_adj)
  
  # truncate the data to match the results we have
  biomass_data <- biomass_data %>% 
    filter(wmin_g >= min_weight_g,
           wmin_g != 0,
           is.na(wmin_g) == FALSE,
           wmax_g != 0,
           is.na(wmax_g) == FALSE)
  
  # Number of fish for the year
  total_abundance <- sum(biomass_data$numlen_adj)
  
  # Have to do not with dplyr:
  wmin.vec <- biomass_data$wmin_g    # Vector of lower weight bin limits
  wmax.vec <- biomass_data$wmax_g    # Vector of upper weight bin limits
  num.vec  <- biomass_data$numlen_adj  # Vector of corresponding counts for those bins
  
  # doing it with purr and pmap
  biomass_bins <- select(biomass_data, wmin_g, wmax_g) %>% 
    distinct(wmin_g, wmax_g)
  
  # Go row-by-row and get how many of any species falls into each
  # discrete size bin
  out_data <- pmap_dfr(biomass_bins, function(wmin_g, wmax_g){
    
    # 1. Count times wmin.vec >= individual wmin_g, multiply by number of fish for year
    countGTEwmin <- sum( (wmin.vec >= wmin_g) * num.vec)
    
    # 2. Count times wmin.vec >= individual wmax_g, multiply by number of fish for year
    lowCount <- sum( (wmin.vec >= wmax_g) * num.vec)
    
    # 3. Count times wmax.vec > individual wmin_g, multiply by number of fish for year
    highCount <- sum( (wmax.vec >  wmin_g) * num.vec)
    
    # 4. build table
    out_table <- data.frame(
      "wmin_g"       = wmin_g,
      "wmax_g"       = wmax_g,
      "countGTEwmin" = countGTEwmin,
      "lowCount"     = lowCount,
      "highCount"    = highCount
    )
  })
  
  
  
  # return the purr version
  return(out_data)
  
  
  
  
}






# Function for ISD plots 
# should work on whatever the grouping variable is as long as the lists match
#' @title Plot Individual Size Distribution Curves
#'
#' @param isd_data_prepped Data from isd_plot_prep
#' @param mle_results Corresponding group results from group_isd_calc
#' @param abundance_used Label to add for context of what abundance source was used
#' @param plot_rects Flag to plot bodymass rectangles
#' @param show_pl_fit Flag to include powerlaw fit
#' @param xlim_global global xlimits so that other plots can be held side-by-side
#' @param group_name name of the group for the data to use as label
#'
#' @return
#' @export
#'
#' @examples
ggplot_isd <- function(isd_data_prepped, 
                       mle_results, 
                       abundance_vals = "observed",
                       plot_rects = TRUE,
                       show_pl_fit = TRUE,
                       xlim_global = NULL,
                       group_name = NULL){
  
  # Columns and labels
  abundance_lab <- c("observed" = "Survey Abundance", 
                     "stratified" = "Stratified Abundance")
  abundance_label <- abundance_lab[[abundance_vals]]
  
  # Power law parameters and summary details for the group of data:
  b.MLE           <- mle_results$b
  total_abundance <- mle_results$n
  b.confMin       <- mle_results$confMin
  b.confMax       <- mle_results$confMax
  
  
  # Create range of x values from the group to get power law predictions
  # PLB = bounded power-law
  # min and max weights for predictions
  xmin  <- mle_results$xmin
  xmax  <- mle_results$xmax
  
  # Create x values (individual bodymass) to predict across
  # break up the Xlim into pieces between min and max
  x.PLB <- seq(from = xmin, 
               to   = xmax, 
               length.out = 2000)   
  
  # get the length of that vector
  x.PLB.length <- length(x.PLB)  
  
  # remove last entry, add an entry .9999 of the way there, and cap it with the last entry wtf
  x.PLB <- c(x.PLB[-x.PLB.length], 0.99999 * x.PLB[x.PLB.length], x.PLB[x.PLB.length])
  
  
  # Y values for plot limits/bounds/predictions from bounded power law pdf
  y.PLB         <- (1 - sizeSpectra::pPLB(x = x.PLB, b = b.MLE, xmin = min(x.PLB), xmax = max(x.PLB))) * total_abundance
  y.PLB.confMin <- (1 - sizeSpectra::pPLB(x = x.PLB, b = b.confMin, xmin = min(x.PLB), xmax = max(x.PLB))) * total_abundance
  y.PLB.confMax <- (1 - sizeSpectra::pPLB(x = x.PLB, b = b.confMax, xmin = min(x.PLB), xmax = max(x.PLB))) * total_abundance
  
  
  # Put it in a df to make it easier
  PLB_df <- data.frame(
    x.PLB   = x.PLB,
    y.PLB   = y.PLB,
    confMin = y.PLB.confMin,
    confMax = y.PLB.confMax)
  
  
  # Coefficient Labels
  group_name <- str_replace_all(group_name, "_", " ")
  the_slope  <- as.character(round(mle_results$b, 3))
  
  
  ####  First Plot  
  p1 <- ggplot()
  
  # Toggle Rectangles
  if(plot_rects == TRUE) {
    p1 <- p1 + geom_rect(data = isd_data_prepped, 
                         aes(xmin = wmin_g, 
                             xmax = wmax_g, 
                             ymin = lowCount, 
                             ymax = highCount),
                         color = "gray70",
                         fill = "transparent")}
  
  # Add segments for bin widths
  p1 <- p1 + geom_segment(data = isd_data_prepped,
                          aes(x = wmin_g, 
                              xend = wmax_g, 
                              y = countGTEwmin, 
                              yend = countGTEwmin),
                          color = "blue")
  
  
  # Set limits to global limits if desired
  if(is.null(xlim_global) == FALSE){
    p1 <- p1 + scale_x_log10(limits = xlim_global, 
                             labels = trans_format("log10", math_format(10^.x)))
  } else{
    p1 <- p1 + scale_x_log10(labels = trans_format("log10", math_format(10^.x)))
  }
  
  
  # Toggle to turn on the fit lines or not
  if(show_pl_fit == TRUE) {
    p1 <- p1 +
      geom_line(data = PLB_df, aes(x.PLB, y.PLB), color = "darkred") +
      geom_line(data = PLB_df, aes(x.PLB, confMin), color = "darkred", linetype = 2) +
      geom_line(data = PLB_df, aes(x.PLB, confMax), color = "darkred", linetype = 2) }
  
  # Finish off labels
  p1 <- p1 + labs(x = "Individual Bodymass (g)",
                  y =  "Fishes with BodyMass \u2265 x",
                  title = group_name, 
                  tag = str_c("b = ", the_slope),
                  caption = str_c(abundance_label, " used for size spectrum slope estimation.")) +
    theme(plot.tag.position = "topright")
  
  
  
  ####  Second Plot, log transformed y axis
  
  # Start plot
  p2 <- p1
  
  # Change the y-axis, log transformation
  p2 <- p2 + scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    labs(x = "Individual Bodymass (g)",
         y =  "Fishes with BodyMass \u2265 x",
         title = group_name, 
         tag = str_c("b = ", the_slope),
         caption = str_c(abundance_label, " used for size spectrum slope estimation.")) +
    theme(plot.tag.position = "topright")
  
  
  # Build the stacked plot
  p3 <- p1 / p2
  
  # Put figures in list
  plot_list <- list(obs_y   = p1,
                    log10_y = p2,
                    stacked = p3)
  
  return(plot_list)
}




# Streamlined plotting function that ignores the data and will draw the fit
# built for speed
isd_lite <- function(wmin_grams, 
                     mle_results, 
                     abundance_vals = "stratified",
                     xlim_global = NULL,
                     group_name = NULL){
  
  
  # Columns and labels
  abundance_cols <- c("observed" = "numlen_adj", 
                      "stratified" = "strat_total_abund_s")
  abundance_col <- abundance_cols[[abundance_vals]]
  abundance_lab <- c("observed" = "Survey Abundance", "stratified" = "Stratified Abundance")
  abundance_label <- abundance_lab[[abundance_vals]]
  
  # Power law parameters and summary details for the group of data:
  
  # Slope
  b.MLE     <- mle_results$b
  
  #  Total Abundance (observed or stratified)
  total_abundance <- sum(wmin_grams[, abundance_col]) 
  
  # Confidence intervals
  b.confMin <- mle_results$confMin
  b.confMax <- mle_results$confMax
  
  
  # Create range of x values from the group to get power law predictions
  # PLB = bounded power-law
  # min and max weights for predictions
  xmin  <- mle_results$xmin
  xmax  <- mle_results$xmax
  
  # break up the Xlim into pieces to predict across
  x.PLB <- seq(from = xmin, 
               to = xmax, 
               length.out = 1000)   
  
  # get the length of that
  x.PLB.length <- length(x.PLB)  
  
  # remove last entry, add an entry .9999 of the way there, and cap it with the last entry wtf
  x.PLB <- c(x.PLB[-x.PLB.length], 0.99999 * x.PLB[x.PLB.length], x.PLB[x.PLB.length])
  
  
  # Y values for plot limits/bounds/predictions from bounded power law pdf
  y.PLB         <- (1 - sizeSpectra::pPLB(x = x.PLB, b = b.MLE, xmin = min(x.PLB), xmax = max(x.PLB))) * total_abundance
  y.PLB.confMin <- (1 - sizeSpectra::pPLB(x = x.PLB, b = b.confMin, xmin = min(x.PLB), xmax = max(x.PLB))) * total_abundance
  y.PLB.confMax <- (1 - sizeSpectra::pPLB(x = x.PLB, b = b.confMax, xmin = min(x.PLB), xmax = max(x.PLB))) * total_abundance
  
  
  # Put it in a df to make it easier
  PLB_df <- data.frame(
    x.PLB   = x.PLB,
    y.PLB   = y.PLB,
    confMin = y.PLB.confMin,
    confMax = y.PLB.confMax)
  
  
  # Coefficient Labels
  group_label <- str_replace_all(group_name, "_", " ")
  the_slope  <- as.character(round(mle_results$b, 3))
  
  
  ####  First Plot  
  p1 <- ggplot() +
    geom_line(data = PLB_df, aes(x.PLB, y.PLB), color = "darkred") +
    geom_line(data = PLB_df, aes(x.PLB, confMin), color = "darkred", linetype = 2) +
    geom_line(data = PLB_df, aes(x.PLB, confMax), color = "darkred", linetype = 2) 
  
  # Set limits to global limits if desired
  if(is.null(xlim_global) == FALSE){
    p1 <- p1 + scale_x_log10(limits = xlim_global, 
                             labels = trans_format("log10", math_format(10^.x)))
  } else{
    p1 <- p1 + scale_x_log10(labels = trans_format("log10", math_format(10^.x)))
  }
  
  # Finish off labels
  p1 <- p1 + labs(x = "Individual Bodymass (g)",
                  y =  "Fishes with BodyMass \u2265 x",
                  title = group_label, 
                  tag = str_c("b = ", the_slope),
                  caption = str_c(abundance_label, " used for size spectrum slope estimation.")) +
    theme(plot.tag.position = "topright")
  
  
  
  ####  Second Plot, log transformed y axis
  # start plot
  p2 <- p1
  
  # Change the y-axis, log transformation
  p2 <- p2 + 
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    labs(x = "Individual Bodymass (g)",
         y =  "Fishes with BodyMass \u2265 x",
         title = group_name, 
         tag = str_c("b = ", the_slope),
         caption = str_c(abundance_label, " used for size spectrum slope estimation.")) +
    theme(plot.tag.position = "topright")
  
  
  # Build the stacked plot
  p3 <- p1 / p2
  
  # Put figures in list
  plot_list <- list(obs_y   = p1,
                    log10_y = p2,
                    stacked = p3)
  
  return(plot_list)
}








####______________________####
#### Binned Spectra Analysis  ####
####______________________####




####__ Base 2 Binning__  ####

# Follows methodology of LBNbiom from Edwards et al. 2016
# Original Method: Blanchard 2005


#' @title Build Log 2 Bin Structure Dataframe
#' 
#' @description Used to build a dataframe containing equally spaced log2 bins for
#' size spectra analysis. Contains details on the left and right limits, midpoint, bin width, 
#' and a text label for the bins. log2bin number ascends with increasing size for eeasy plotting.
#'
#' @param log2_min 
#' @param log2_max 
#' @param log2_increment 
#'
#' @return
#' @export
#'
#' @examples
define_log2_bins <- function(log2_min = 0, log2_max = 13, log2_increment = 1){
  
  # How many bins
  n_bins  <- length(seq(log2_max, log2_min + log2_increment, by = -log2_increment))
  
  # Build Equally spaced log2 bin df
  log2_bin_structure <- data.frame(
    "log2_bins" = as.character(seq(n_bins, 1, by = -1)),
    "left_lim"  = seq(log2_max - log2_increment, log2_min, by = -log2_increment),
    "right_lim" = seq(log2_max, log2_min + log2_increment, by = -log2_increment)) %>% 
    mutate(
      bin_label    = str_c(round(2^left_lim, 3), " - ", round(2^right_lim, 3), "g"),
      bin_width    = 2^right_lim - 2^left_lim,
      bin_midpoint = (2^right_lim + 2^left_lim) / 2) %>% 
    arrange(left_lim)
  
  return(log2_bin_structure)
}




#' @title Assign Manual log2 Bodymass Bins
#'
#' @description Manually assign log2 bins based on individual length-weight bodymass 
#' in increments of 1 on the log2 scale. Returns data with bins assigned based on individual
#' length-weight biomass
#' 
#' Uses maximum weight, and its corresponding bin as the limit.
#'
#' @param wmin_grams Catch data prepared for mle calculation, use prep_wmin_wmax
#' @param log2_increment Equally spaced increments to use for log 2 bin sizes. Default = 0.5.
#'
#' @return
#' @export
#'
#' @examples
assign_log2_bins <- function(wmin_grams, log2_increment = 1){
  
  
  #### 1. Set up bodymass bins
  
  # filter missing weights
  size_data <- wmin_grams %>% 
    filter(wmin_g > 0,
           is.na(wmin_g) == FALSE,
           wmax_g > 0,
           is.na(wmax_g) == FALSE)
  
  # Get bodymass on log2() scale
  size_data$log2_weight <- log2(size_data$ind_weight_g)
  
  # Set up the bins - Pull min and max weights from data available
  #min_bin <- floor(min(size_data$log2_weight))
  min_bin <- 0
  max_bin <- ceiling(max(size_data$log2_weight))
  
  
  # Build a bin key, could be used to clean up the incremental assignment or for apply style functions
  log2_bin_structure <- define_log2_bins(
    log2_min = min_bin, 
    log2_max = max_bin, 
    log2_increment = log2_increment)
  
  
  
  # Loop through bins to assign the bin details to the data
  log2_assigned <- log2_bin_structure %>%
    split(.$log2_bins) %>%
    map_dfr(function(log2_bin){
      
      # limits and labels
      l_lim   <- log2_bin$left_lim
      r_lim   <- log2_bin$right_lim
      bin_num <- as.character(log2_bin$log2_bin)
      
      # assign the label to the appropriate bodymasses
      size_data %>% mutate(
        log2_bins = ifelse( between(log2_weight, l_lim, r_lim), bin_num, NA),
        log2_bins = as.character(log2_bins)) %>%
        drop_na(log2_bins)
      
    })
  
  # Join in the size bins
  log2_assigned <- left_join(log2_assigned, log2_bin_structure, by = "log2_bins")
  
  # return the data with the bins assigned
  return(log2_assigned)
  
}




#' @title Calculate Normalized and De-Normalized Abundances
#'
#' @description For binned size spectra estimation we use the stratified abundance divided by the
#' bin widths (normalized size spectra). Another way to present the data is to de-normalize, which 
#' takes those values and multiplies them by the mid-point of the log-scale bins.
#' 
#' min/max & bin_increments are used to enforce the presence of a size bin in the event that 
#' there is no abundance. This is done for comparing across different groups/areas that should 
#' conceivably have the same size range sampled.
#'
#' @param log2_assigned size data containing the bin assignments to use
#' @param min_log2_bin Minimum 2^x value for the size spectra being measured (>=)
#' @param max_log2_bin Maximum 2^x value for the size spectra being measured (<)
#' @param bin_increment The bin-width on log scale that separates each bin
#' @param ... Additional grouping factors with which to aggregate on besides the size bins themselves
#'
#' @return
#' @export
#'
#' @examples
aggregate_log2_bins <- function(
    log2_assigned, 
    min_log2_bin = 0, 
    max_log2_bin = 13, 
    bin_increment = 1,
    ...){
  
  # Full Possible Bin Structure
  # Fills in any gaps
  log2_bin_structure <- define_log2_bins(
    log2_min = min_log2_bin, 
    log2_max = max_log2_bin, 
    log2_increment = bin_increment)
  
  
  # Capture all the group levels with a cheeky expand()
  if(missing(...) == FALSE){
    log2_bin_structure <- log2_bin_structure %>% 
      expand(left_lim, distinct(log2_assigned, ...)) %>% 
      full_join(log2_bin_structure)
  }
  
  
  
  # Get bin breaks
  log2_breaks <- sort(unique(c(log2_bin_structure$left_lim, log2_bin_structure$right_lim)))
  
  
  # Get Totals for bodymass and abundances
  log2_aggregates <- log2_assigned %>% 
    group_by(left_lim, ...) %>% 
    summarise(observed_abundance   = sum(numlen_adj, na.rm = T),
              observed_weight_g    = sum(wmin_g, na.rm = T),
              stratified_abundance = sum(strat_total_abund_s, na.rm = T),
              stratified_weight_g  = sum(wmin_area_strat, na.rm = T),
              .groups = "drop")
  
  
  # Join back in what the limits and labels are
  # The defined bins and their labels enforce the size limits
  log2_prepped <- left_join(x = log2_bin_structure, 
                           y = log2_aggregates)
  
  
  #### Fill Gaps with Zero's?? 
  # This ensures that any size bin that is intended to be in use is actually used
  log2_prepped <- log2_prepped %>% 
    mutate(observed_abundance = ifelse(is.na(observed_abundance), 1, observed_abundance),
           stratified_abundance = ifelse(is.na(stratified_abundance), 1, stratified_abundance),
           observed_weight_g = ifelse(is.na(observed_weight_g), 1, observed_weight_g),
           stratified_weight_g = ifelse(is.na(stratified_weight_g), 1, stratified_weight_g))
  
  
  #### normalize abundances using the bin widths
  log2_prepped <- log2_prepped %>% 
    mutate(
      normalized_abund = observed_abundance / bin_width,
      norm_strat_abund = stratified_abundance / bin_width,
      # # Remove Bins Where Normalized Biomass < 0? No!
      # normalized_abund = ifelse(normalized_abund < 2^0, NA, normalized_abund),
      # norm_strat_abund = ifelse(norm_strat_abund < 2^0, NA, norm_strat_abund)
    )
  
  # Add de-normalized abundances (abundance * bin midpoint)
  log2_prepped <- log2_prepped %>% 
    mutate(
      denorm_abund = normalized_abund * bin_midpoint,
      denorm_strat_abund = norm_strat_abund * bin_midpoint)
  
  # Return the aggregations
  return(log2_prepped)
  
}






# Return the slope and intercept for manual estimation of size spectrum slope
#' @title Process Size Spectrum Slope and Intercept using log2 bins for a group of data
#'
#' @param wmin_grams Catch data prepped using prep_wmin_wmax
#' @param .group_cols Vector of string(s) to use as grouping factors that correspond to column names
#' @param min_weight_g Minimum weight cutoff in grams
#' @param min_log2_bin Minimum 2^x value for the size spectra being measured (>=). Left limit of smallest bin.
#' @param max_log2_bin Maximum 2^x value for the size spectra being measured (<). Left limit of largest bin.
#' @param bin_increment The bin-width on log scale that separates each bin
#'
#' @return
#' @export
#'
#' @examples
group_log2_spectra <- function(wmin_grams,
                              .group_cols = "Year",
                              min_weight_g = 1,
                              min_log2_bin = 0,
                              max_log2_bin = 13,
                              bin_increment = 1){
  
  
  # 1. Set bodymass lower limit, and assign bin labels
  log2_assigned <- filter(wmin_grams, wmin_g >= min_weight_g )
  
  
  # 2. Build group_var from desired grouping columns
  log2_assigned <- log2_assigned %>% 
    unite(col = "group_var", 
          {{.group_cols}}, 
          sep = "_", 
          remove = FALSE, 
          na.rm = FALSE)
  
  # 3. Make a table of constituent combinations
  col_syms <- syms(.group_cols)
  grouping_table <- log2_assigned %>% 
    distinct(!!!col_syms, group_var)
  
  
  # 4. Add missing group factors as "all data"
  grouping_table <- add_missing_groups(group_dataframe = grouping_table)
  
  
  
  #### 5.  Running log2 slopes/intercepts
  
  # Run size spectra slopes for each group
  group_results <- log2_assigned %>% 
    split(.$group_var) %>% 
    imap_dfr(function(log2_assigned, group_label){
      
      # Total the abundances for each bin, normalize them
      log2_prepped <- aggregate_log2_bins(
        log2_assigned = log2_assigned, 
        min_log2_bin = min_log2_bin, 
        max_log2_bin = max_log2_bin, 
        bin_increment = bin_increment)
      
      
      ####  Run linear models for slopes
      # Formula for size spectrum slope: log2(abundance) ~ log2 bin, 
      # aka the minimum bodymass for bin
      
      # Abundance from area-stratification
      # Will throw an error if norm_strat_abund = 0
      
      # If the response is abundance then its a normalized abundance spectra
      # If the response were biomass that would be normalized biomass spectra
      lm_abund_strat <- lm(log2(norm_strat_abund) ~ log2(bin_midpoint), data = log2_prepped)
      
      
      #  Pull Slope
      # See edwards et al 2016 methodology Table 2.
      # Really not sure whether to add or subtract values:
      # Depends on what you are comparing to, 
      # the normalized biomass spectra slope is
      # proportional to the abundance spectra (lambda) 
      # as lambda + 1, so we can subtract one to get a comparable value
      
      
      # We fit the normalized abundance spectra, so no adjustment
      lm_b_strat <- lm_abund_strat$coeff[2]
      
      # Pull intercept  
      lm_int_strat <- lm_abund_strat$coeff[1] 
      
      # R squared
      lm_rsqr_strat <- summary(lm_abund_strat)$adj.r.squared
      
      # sig
      lm_sig_strat <- broom::tidy(lm_abund_strat)$p.value[2]
      

      # Put Spectra Details in Table
      log2_results <- data.frame(
        group_var       = group_label,
        log2_slope_strat = lm_b_strat,
        log2_int_strat   = lm_int_strat,
        log2_rsqr_strat  = lm_rsqr_strat,
        log2_sig_strat   = lm_sig_strat
      )
      
      # Return the group results
      return(log2_results)
      
    })
  
  
  # Merge in the group details
  group_results <- full_join(grouping_table, group_results, by = "group_var") %>% 
    mutate(across(c(group_var, Year, season, survey_area, decade), as.character))
  
  # Return the results
  return(group_results)
  
}




#' @title Process log2 binned size spectra for the collections of WARMEM Groups
#' 
#' 
#' @description Takes a dataframe of all data (weight in grams) and performs the
#' group_log2_spectra()
#'
#' @param wmin_grams Catch data prepped using assign_log2_bins
#' @param min_weight_g Minimum weight cutoff in grams
#' @param max_weight_g Maximum weight cutoff in grams
#' @param min_log2_bin Minimum 2^x value for the size spectra being measured (>=)
#' @param max_log2_bin Maximum 2^x value for the size spectra being measured (<)
#' @param bin_increment The bin-width on log scale that separates each bin
#'
#' @return
#' @export
#'
#' @examples
warmem_log2_estimates <- function(wmin_grams,
                                 min_weight_g,
                                 max_weight_g,
                                 min_log2_bin = 0,
                                 max_log2_bin = 13,
                                 bin_increment = 1){
  
  
  ####__  Set Bodymass Cutoff and Groups
  
  # 1. Filter for lower/upper ends of gear selectivity
  log2_assigned_trunc <- filter(
    wmin_grams, 
    wmin_g >= min_weight_g,
    wmin_g < max_weight_g) %>% 
    mutate(decade = floor_decade(Year))
  
  
  ####_ 2. Every year, entire survey 
  message("Calculating log(2) size spectrum slope each year")
  g5_res <- log2_assigned_trunc  %>% 
    group_log2_spectra(min_weight_g = min_weight_g, 
                      .group_cols = c("Year")) 
  
  ####_ 6. every year, every region 
  message("Calculating log(2) size spectrum slope each year in each region")
  g6_res <- log2_assigned_trunc  %>% 
    group_log2_spectra(min_weight_g = min_weight_g, 
                      .group_cols = c("Year", "survey_area")) 
  
  ####_ 7. every year, only seasons 
  message("Calculating log(2) size spectrum slope each year in each season")
  g7_res <- log2_assigned_trunc %>% 
    group_log2_spectra(min_weight_g = min_weight_g, 
                      .group_cols = c("Year", "season"))
  
  
  ####_ 8. every year, region * season 
  message("Calculating log(2) size spectrum slope each year in each region, for every season")
  g8_res <- log2_assigned_trunc %>%
    group_log2_spectra(min_weight_g = min_weight_g,
                      .group_cols = c("Year", "season", "survey_area"))
  
  
  # Put the reults in one table with an ID for how they groups are set up
  table_complete <- bind_rows(
    list(
      "single years"                   = g5_res,
      "single years * region"          = g6_res,
      "single years * season"          = g7_res,
      "single years * season * region" = g8_res), 
    .id = "group ID")
  
  # Return the summary table
  return(table_complete)
  
}












####__ Base 10 Binning__  ####


#' @title Build Log 10 Bin Structure Dataframe
#' 
#' @description Used to build a dataframe containing equally spaced log10 bins for
#' size spectra analysis. Contains details on the left and right limits, midpoint, bin width, 
#' and a text label for the bins. l10bin number ascends with increasing size for eeasy plotting.
#'
#' @param l10_min 
#' @param l10_max 
#' @param l10_increment 
#'
#' @return
#' @export
#'
#' @examples
define_l10_bins <- function(l10_min = 0, l10_max = 4.5, l10_increment = 0.5){
  
  # How many bins
  n_bins  <- length(seq(l10_max, l10_min + l10_increment, by = -l10_increment))
  
  # Build Equally spaced log10 bin df
  l10_bin_structure <- data.frame(
    "log10_bins" = as.character(seq(n_bins, 1, by = -1)),
    "left_lim"  = seq(l10_max - l10_increment, l10_min, by = -l10_increment),
    "right_lim" = seq(l10_max, l10_min + l10_increment, by = -l10_increment)) %>% 
    mutate(
      bin_label    = str_c(round(10^left_lim, 3), " - ", round(10^right_lim, 3), "g"),
      bin_width    = 10^right_lim - 10^left_lim,
      bin_midpoint = (10^right_lim + 10^left_lim) / 2) %>% 
    arrange(log10_bins)
  
  return(l10_bin_structure)
}





#' @title Assign Manual log10 Bodymass Bins
#'
#' @description Manually assign log10 bins based on individual length-weight bodymass 
#' in increments of 0.5 on the log10 scale. Returns data with bins assigned based on individual
#' length-weight biomass
#' 
#' Uses maximum weight, and its corresponding bin as the limit.
#'
#' @param wmin_grams Catch data prepared for mle calculation, use prep_wmin_wmax
#' @param l10_increment Equally spaced increments to use for log 10 bin sizes. Default = 0.5.
#'
#' @return
#' @export
#'
#' @examples
assign_log10_bins <- function(wmin_grams, l10_increment = 1){
  
  
  #### 1. Set up bodymass bins
  
  # filter missing weights
  size_data <- wmin_grams %>% 
    filter(wmin_g > 0,
           is.na(wmin_g) == FALSE,
           wmax_g > 0,
           is.na(wmax_g) == FALSE)

  # Get bodymass on log10() scale
  size_data$log10_weight <- log10(size_data$ind_weight_g)
  
  # Set up the bins - Pull min and max weights from data available
  min_bin <- floor(min(size_data$log10_weight))
  max_bin <- ceiling(max(size_data$log10_weight))
  
  
  # Build a bin key, could be used to clean up the incremental assignment or for apply style functions
  l10_bin_structure <- define_l10_bins(
    l10_min = min_bin, 
    l10_max = max_bin, 
    l10_increment = l10_increment)
  
  
  
  # Loop through bins to assign the bin details to the data
  l10_assigned <- l10_bin_structure %>%
    split(.$log10_bins) %>%
    map_dfr(function(l10_bin){
      
      # limits and labels
      l_lim   <- l10_bin$left_lim
      r_lim   <- l10_bin$right_lim
      bin_num <- as.character(l10_bin$log10_bin)
      
      # assign the label to the appropriate bodymasses
      size_data %>% mutate(
        log10_bins = ifelse( between(log10_weight, l_lim, r_lim), bin_num, NA),
        log10_bins = as.character(log10_bins)) %>%
        drop_na(log10_bins)
      
    })
  
  # Join in the size bins
  l10_assigned <- left_join(l10_assigned, l10_bin_structure, by = "log10_bins")
  
  # return the data with the bins assigned
  return(l10_assigned)
  
}







#' @title Calculate Normalized and De-Normalized Abundances
#'
#' @description For binned size spectra estimation we use the stratified abundance divided by the
#' bin widths (normalized size spectra). Another way to present the data is to de-normalize, which 
#' takes those values and multiplies them by the mid-point of the log-scale bins.
#' 
#' min/max & bin_increments are used to enforce the presence of a size bin in the event that 
#' there is no abundance. This is done for comparing across different groups/areas that should 
#' conceivably have the same size range sampled.
#'
#' @param l10_assigned size data containing the bin assignments to use
#' @param min_l10_bin Minimum 10^x value for the size spectra being measured (>=)
#' @param max_l10_bin Maximum 10^x value for the size spectra being measured (<)
#' @param bin_increment The bin-width on log scale that separates each bin
#' @param ... Additional grouping factors with which to aggregate on besides the size bins themselves
#'
#' @return
#' @export
#'
#' @examples
aggregate_l10_bins <- function(
    l10_assigned, 
    min_l10_bin = 0, 
    max_l10_bin = 4.5, 
    bin_increment = 0.5,
    ...){
  
  # Full Possible Bin Structure
  # Fills in any gaps
  l10_bin_structure <- define_l10_bins(
    l10_min = min_l10_bin, 
    l10_max = max_l10_bin, 
    l10_increment = bin_increment)
  
  
  # Capture all the group levels with a cheeky expand()
  if(missing(...) == FALSE){
    l10_bin_structure <- l10_bin_structure %>% 
      expand(left_lim, distinct(l10_assigned, ...)) %>% 
      full_join(l10_bin_structure)
  }
  
  
  
  # Get bin breaks
  l10_breaks <- sort(unique(c(l10_bin_structure$left_lim, l10_bin_structure$right_lim)))
  
  
  # Get Totals for bodymass and abundances
  l10_aggregates <- l10_assigned %>% 
    # group_by(log10_bins, ...) %>% 
    group_by(left_lim, ...) %>% 
    summarise(observed_abundance   = sum(numlen_adj, na.rm = T),
              observed_weight_g    = sum(wmin_g, na.rm = T),
              stratified_abundance = sum(strat_total_abund_s, na.rm = T),
              stratified_weight_g  = sum(wmin_area_strat, na.rm = T),
              .groups = "drop")
  
  
  # Join back in what the limits and labels are
  # The defined bins and their labels enforce the size limits
  l10_prepped <- left_join(x = l10_bin_structure, 
                           y = l10_aggregates)
  
  
  #### Fill Gaps with Zero's??
  # This ensures that any size bin that is intended to be in use is actually used
  l10_prepped <- l10_prepped %>% 
    mutate(observed_abundance = ifelse(is.na(observed_abundance), 1, observed_abundance),
           stratified_abundance = ifelse(is.na(stratified_abundance), 1, stratified_abundance),
           observed_weight_g = ifelse(is.na(observed_weight_g), 1, observed_weight_g),
           stratified_weight_g = ifelse(is.na(stratified_weight_g), 1, stratified_weight_g))
  
  
  #### normalize abundances using the bin widths
  l10_prepped <- l10_prepped %>% 
    mutate(
      normalized_abund = observed_abundance / bin_width,
      norm_strat_abund = stratified_abundance / bin_width,
      # # Remove Bins Where Normalized Biomass < 0? No!
      # normalized_abund = ifelse(normalized_abund < 10^0, NA, normalized_abund),
      # norm_strat_abund = ifelse(norm_strat_abund < 10^0, NA, norm_strat_abund)
    )
  
  # Add de-normalized abundances (abundance * bin midpoint)
  l10_prepped <- l10_prepped %>% 
    mutate(
      denorm_abund = normalized_abund * bin_midpoint,
      denorm_strat_abund = norm_strat_abund * bin_midpoint)
  
  # Return the aggregations
  return(l10_prepped)
  
}


# Return the slope and intercept for manual estimation of size spectrum slope
#' @title Process Size Spectrum Slope and Intercept using log10 bins for a group of data
#'
#' @param wmin_grams Catch data prepped using prep_wmin_wmax
#' @param .group_cols Vector of string(s) to use as grouping factors that correspond to column names
#' @param min_weight_g Minimum weight cutoff in grams
#' @param min_l10_bin Minimum 10^x value for the size spectra being measured (>=). Left limit of smallest bin.
#' @param max_l10_bin Maximum 10^x value for the size spectra being measured (<). Left limit of largest bin.
#' @param bin_increment The bin-width on log scale that separates each bin
#'
#' @return
#' @export
#'
#' @examples
group_l10_spectra <- function(wmin_grams,
                              .group_cols = "Year",
                              min_weight_g = 1,
                              min_l10_bin = 0,
                              max_l10_bin = 4.5,
                              bin_increment = 0.5){
  
  
  # 1. Set bodymass lower limit, and assign bin labels
  l10_assigned <- filter(wmin_grams, wmin_g >= min_weight_g )
  
  
  # 2. Build group_var from desired grouping columns
  l10_assigned <- l10_assigned %>% 
    unite(col = "group_var", 
          {{.group_cols}}, 
          sep = "_", 
          remove = FALSE, 
          na.rm = FALSE)
  
  # 3. Make a table of constituent combinations
  col_syms <- syms(.group_cols)
  grouping_table <- l10_assigned %>% 
    distinct(!!!col_syms, group_var)
  
  
  # 4. Add missing group factors as "all data"
  grouping_table <- add_missing_groups(group_dataframe = grouping_table)
  
  
  
  #### 5.  Running log10 slopes/intercepts
  
  # Run size spectra slopes for each group
  group_results <- l10_assigned %>% 
    split(.$group_var) %>% 
    imap_dfr(function(l10_assigned, group_label){
      
      # Total the abundances for each bin, normalize them
      l10_prepped <- aggregate_l10_bins(
        l10_assigned = l10_assigned, 
        min_l10_bin = min_l10_bin, 
        max_l10_bin = max_l10_bin, 
        bin_increment = bin_increment)
      
      
      ####  Run linear models for slopes
      # Formula for size spectrum slope: log10(abundance) ~ log10 bin, 
      # aka the minimum bodymass for bin
      
      # Abundance from area-stratification
      # Will throw an error if norm_strat_abund = 0
      # lm_abund_strat <- lm(log10(norm_strat_abund) ~ left_lim, data = l10_prepped)
      lm_abund_strat <- lm(log10(norm_strat_abund) ~ log10(bin_midpoint), data = l10_prepped)
  
      #  Pull Slope
      # See edwards et al 2016 methodology Table 2.
      # Really not sure whether to add or subtract values:
      # Depends on what you are comparing to, 
      # the normalized biomass spectra slope is
      # proportional to the abundance spectra (lambda) 
      # as lambda + 1, so we can subtract one to get a comparable value
      
      
      # We fit the normalized abundance spectra, so no adjustment
      lm_b_strat <- lm_abund_strat$coeff[2]
      
      # Pull intercept  
      lm_int_strat <- lm_abund_strat$coeff[1] 
      
      # R squared
      lm_rsqr_strat <- summary(lm_abund_strat)$adj.r.squared
      
      # sig
      lm_sig_strat <- broom::tidy(lm_abund_strat)$p.value[2]
      
      # Run for abundance from the survey catch
      #lm_abund       <- lm(log10(normalized_abund) ~ left_lim, data = l10_prepped)
      # pull out slope coefficient
      #lm_b <- lm_abund$coeff[2] #- 1
      #lm_int <- lm_abund$coeff[1] 
      #lm_rsqr       <- summary(lm_abund)$adj.r.squared
      #lm_sig       <- broom::tidy(lm_abund)$p.value[2]
      
      
      
      # Put Spectra Details in Table
      l10_results <- data.frame(
        group_var       = group_label,
        l10_slope_strat = lm_b_strat,
        l10_int_strat   = lm_int_strat,
        l10_rsqr_strat  = lm_rsqr_strat,
        l10_sig_strat   = lm_sig_strat#,
        #l10_slope       = lm_b,
        #l10_rsqr        = lm_rsqr,
        #l10_int         = lm_int,
        #l10_sig         = lm_sig,
      )
      
      # Return the group results
      return(l10_results)
  
  })
  
  
  # Merge in the group details
  group_results <- full_join(grouping_table, group_results, by = "group_var") %>% 
    mutate(across(c(group_var, Year, season, survey_area, decade), as.character))
  
  # Return the results
  return(group_results)
  
}



#' @title Process log10 binned size spectra for the collections of WARMEM Groups
#' 
#' 
#' @description Takes a dataframe of all data (weight in grams) and performs the
#' group_log10_spectra()
#'
#' @param wmin_grams Catch data prepped using assign_log10_bins
#' @param min_weight_g Minimum weight cutoff in grams
#' @param max_weight_g Maximum weight cutoff in grams
#' @param min_l10_bin Minimum 10^x value for the size spectra being measured (>=)
#' @param max_l10_bin Maximum 10^x value for the size spectra being measured (<)
#' @param bin_increment The bin-width on log scale that separates each bin
#'
#' @return
#' @export
#'
#' @examples
warmem_l10_estimates <- function(wmin_grams,
                                min_weight_g,
                                max_weight_g,
                                min_l10_bin = 0,
                                max_l10_bin = 4.5,
                                bin_increment = 0.5){
  
  
  ####__  Set Bodymass Cutoff and Groups
  
  # 1. Filter for lower/upper ends of gear selectivity
  l10_assigned_trunc <- filter(
    wmin_grams, 
    wmin_g >= min_weight_g,
    wmin_g < max_weight_g) %>% 
    mutate(decade = floor_decade(Year))
  
  
  # 2. Set up the factor groupings we want to compare : 
  
  ####_ 5. Every year, entire survey 
  message("Calculating log(10) size spectrum slope each year")
  g5_res <- l10_assigned_trunc  %>% 
    group_l10_spectra(min_weight_g = min_weight_g, 
                     .group_cols = c("Year")) 
  
  ####_ 6. every year, every region 
  message("Calculating log(10) size spectrum slope each year in each region")
  g6_res <- l10_assigned_trunc  %>% 
    group_l10_spectra(min_weight_g = min_weight_g, 
                     .group_cols = c("Year", "survey_area")) 
  
  ####_ 7. every year, only seasons 
  message("Calculating log(10) size spectrum slope each year in each season")
  g7_res <- l10_assigned_trunc %>% 
    group_l10_spectra(min_weight_g = min_weight_g, 
                     .group_cols = c("Year", "season"))
  
  
  ####_ 8. every year, region * season 
  message("Calculating log(10) size spectrum slope each year in each region, for every season")
  g8_res <- l10_assigned_trunc %>%
    group_l10_spectra(min_weight_g = min_weight_g,
                     .group_cols = c("Year", "season", "survey_area"))
  
  
  # Put the reults in one table with an ID for how they groups are set up
  table_complete <- bind_rows(
    list(
      "single years"                   = g5_res,
      "single years * region"          = g6_res,
      "single years * season"          = g7_res,
       "single years * season * region" = g8_res), 
    .id = "group ID")
  
  # Return the summary table
  return(table_complete)
  
}




####______________________####
####  Plotting l10 Spectra  ####

# Plot the bins and the size spectrum slope
#' @title Plot log-10 bin size spectrum with slope fits
#'
#' @param l10_assigned Catch data prepped using assign_log10_bins
#'
#' @return
#' @export
#'
#' @examples
plot_log10_ss <- function(l10_assigned){
  
  # Bin aggregation moved to aggregate_log10_bins
  
  # Get totals for each bin:
  l10_prepped <- aggregate_l10_bins(l10_assigned)
  
  
  #### Plots Correcting for the bin widths
  norm_abund_plot <- l10_prepped %>% 
    ggplot(aes(left_lim, normalized_abund)) +
    geom_col(fill = gmri_cols("green")) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) + 
    labs(x = "Log10(bodymass)", y = "Normalized Abundances",
         subtitle = "Abundances Divided by Bin Widths")
  
  norm_strat_abund_plot <- l10_prepped  %>% 
    ggplot(aes(left_lim, norm_strat_abund)) +
    geom_col(fill = gmri_cols("green")) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) + 
    labs(x = "Log10(bodymass)", y = "Normalized Stratified Abundances")
  
  # Abundances from Survey
  p1 <- norm_abund_plot +
    geom_smooth(formula = y ~ x,
                method = "lm",
                color = gmri_cols("orange")) +
    stat_poly_eq(formula = y ~ x,
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 label.y = 1.1, 
                 parse = TRUE) +
    labs(caption = "Survey Abundance", subtitle = "", title = "Manual log10 Binning")
  
  # Stratified Abundances
  p2 <- norm_strat_abund_plot +
    geom_smooth(formula = y ~ x,
                method = "lm",
                color = gmri_cols("orange")) +
    stat_poly_eq(formula = y ~ x,
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 label.y = 1.1, 
                 parse = TRUE) +
    labs(caption = "Stratified Abundance")
  
  # assemble with patchwork
  plot_out <- (p1 | p2)
  return(plot_out)
  
  
}


#' @title Plot normalized abundance of Size Spectra
#' 
#' @description Single panel plot of the distribution of abundance on log size bins.
#'
#' @param l10_assigned Catch data prepped using assign_log10_bins
#' @param stratified Stratified or survey abundances
#'
#' @return
#' @export
#'
#' @examples
plot_normalized_ss <- function(l10_assigned, stratified = TRUE){
  
  # Get totals for each bin:
  l10_prepped <- aggregate_l10_bins(l10_assigned)
  
  #### Plots Correcting for the bin widths
  norm_strat_abund_plot <- l10_prepped  %>% 
    ggplot(aes(left_lim, norm_strat_abund)) +
    geom_col(fill = gmri_cols("green")) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) + 
    labs(x = "Log10(bodymass)", y = "Normalized Stratified Abundances")
  
  
  # Stratified Abundances
  p1 <- norm_strat_abund_plot +
    geom_smooth(formula = y ~ x,
                method = "lm",
                color = gmri_cols("orange")) +
    stat_poly_eq(formula = y ~ x,
                 aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 label.y = 1.1, 
                 parse = TRUE) +
    labs(caption = "Abundance normalized by bin-width.")
  
  # Return plot
  return(p1)
}


#' @title Plot De-normalized abundances of size spectra
#'
#' @param l10_assigned Catch data prepped using assign_log10_bins
#' @param stratified Stratified or survey abundances
#'
#' @return
#' @export
#'
#' @examples
plot_denormalized_ss <- function(l10_assigned, stratified = TRUE){
  
  # Get totals for each bin:
  l10_prepped <- aggregate_l10_bins(l10_assigned)
  
  #### Plots Correcting for the bin widths
  p1 <- l10_prepped  %>% 
    ggplot(aes(left_lim, denorm_strat_abund)) +
    geom_point(color = gmri_cols("green")) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) + 
    labs(x = "Log10(bodymass)", 
         y = "De-Normalized Stratified Abundances",
         caption = "Normalized abundance, multiplied by bin mid-point.")
  
  # Return plot
  return(p1)
  
  
}



####______________________####
####  Size Composition ####




# Get weighted mean lengths and weights using 
# tow-level and total stratified abundances
# Note: uses numlen because individuals were ID'd and measured
group_size_metrics <- function(
    size_data, 
    .group_cols = "Year", 
    abund_vals = "numlen_adj"){
  
  
  
  # 1. Build group_level from desired group columns
  group_size_data <- size_data %>% 
    filter(is.na(ind_weight_kg) == FALSE) %>% 
    mutate(decade = floor_decade(Year)) %>% 
    unite(col = "group_var", 
          {{.group_cols}}, 
          sep = "_", 
          remove = FALSE, 
          na.rm = FALSE)
  
  # 2. Make a table of constituent combinations
  col_syms <- syms(.group_cols)
  grouping_table <- group_size_data %>% 
    distinct(!!!col_syms, group_var)
  
  # 3. Add missing groups as "all data"
  grouping_table <- add_missing_groups(group_dataframe = grouping_table)
  
  # Set what abundance column to use
  #we can use stratified rates instead of strat abundance because the proportions are the same
  weighting_col <- switch(abund_vals,
    "numlen_adj" = "numlen_adj",
    "stratified" = "strat_mean_abund_s")
  
  
  # Run Min/Max/Avg. Size for the group
  group_results <- group_size_data %>% 
    split(.$group_var) %>% 
    imap_dfr(function(group_data, group_name){
      
      # Length (measured for all)
      mean_len    <- weighted.mean(group_data[, "length_cm"], group_data[, weighting_col], na.rm = T)
      min_len     <- min(group_data[, "length_cm"], na.rm = T)
      max_len     <- max(group_data[, "length_cm"], na.rm = T)
      
      # Weights (derived from length)
      mean_weight <- weighted.mean(group_data[,"ind_weight_kg"], group_data[, weighting_col], na.rm = T)
      min_weight  <- min(group_data[, "ind_weight_kg"], na.rm = T)
      max_weight  <- max(group_data[, "ind_weight_kg"], na.rm = T)
      
      # Abundance
      total_abund <- sum(group_data[, "numlen_adj"], na.rm = T)
      
      # number of species
      num_species <- group_data %>% 
        filter(is.na(comname) == FALSE) %>% 
        distinct(comname) %>% 
        length()
     
      # Put in table
      table_out <- data.frame(
          "group_var"    = group_name,
          "n_species"    = num_species,  
          "survey_abund" = total_abund,
          "mean_len_cm"  = mean_len,
          "min_len_cm"   = min_len,
          "max_len_cm"   = max_len,
          "mean_wt_kg"   = mean_weight,
          "min_wt_kg"    = min_weight,
          "max_wt_kg"    = max_weight) %>% 
        
        # replace Inf with NA
        mutate(across(.cols = survey_abund:max_wt_kg, .fns = ~ifelse(is.infinite(abs(.x)), NA, .x)))
      
      # return the table
      return(table_out)
      })
  
  # Merge in the group details 
  group_results <- full_join(grouping_table, group_results, by = "group_var") %>% 
    mutate(across(c(group_var, 
                    Year, 
                    decade, 
                    season, 
                    survey_area,
                    hare_group), as.character))
  
  # Return the results
  return(group_results)
  
  
}





# Run all the groups, preserve the groups not stated for "overall" levels
# Direct match to the groups in warmem_isd_estimates
mean_sizes_all_groups <- function(size_data, 
                                  min_weight_g = 0, 
                                  abund_vals = "numlen_adj"){
  
  
  ####__  Set Bodymass Cutoff and Groups
  
  # 1. Set bodymass lower limit
  # Used to filter for lower end of gear selectivity
  group_size_data <- filter(size_data, 
                           ind_weight_kg >= min_weight_g * 1000) %>% 
    mutate(decade = floor_decade(Year))
  
  
  # 2. Set up the factor groupings we want to compare : 

  
  #####__ 1. Every year, entire survey
  message("Processing body size change for: Each year")
  g1_res <- group_size_data  %>% 
    group_size_metrics(.group_cols = c("Year"),
                       abund_vals = abund_vals) 
  
  #####__ 2. every year, every region
  message("Processing body size change for: Each year in each region")
  g2_res <- group_size_data  %>% 
    group_size_metrics(.group_cols = c("Year", "survey_area"),
                       abund_vals = abund_vals) 
  
  ####__3. Year * Region * Functional Group
  message("Processing body size change for: Each Year in each area, for each functional group")
  g3_res <- group_size_data %>% 
    group_size_metrics(.group_cols = c("Year", "survey_area", "hare_group"),
                       abund_vals = abund_vals) 
  
  
  # Put the reults in one table with an ID for how they groups are set up
  table_complete <- bind_rows(
    list(
      "single years"                             = g1_res,
      "single years * region"                    = g2_res,
      "single years * region * functional group" = g3_res), 
    .id = "group ID")
  
  # Return the summary table
  return(table_complete)
  
}


####______________________####
####  Growth Characteristics ####

# pick species
select_vonbert_species <- function(survdat_biological, rank_cutoff = 17){
  
  # Status
  "work in progress"
  
  # Rank species by how many measurements there are
  species_abunds <- survdat_biological %>% 
    count(comname) %>% 
    arrange(desc(n)) # ordered by number measured
  
  # Reorder alphabetically
  vonbert_species <- sort(species_abunds$comname[1:rank_cutoff])
  
  # Name the list so it carries through
  names(vonbert_species) <- vonbert_species
  
  # Put each species into its own list by common name
  vonbert_species_data <- map(vonbert_species, function(species_name){
    
    # Drop NA ages, set increment labels
    survdat_biological %>%
      filter(comname == species_name) %>% 
      mutate(yearclass = est_year - age) 
  })
  
  return(vonbert_species_data)
}




# set time increments
set_vonbert_groups <- function(cutoff_rank){
  "work in progress"
}



# process growth coef
estimate_vonbert_coef <- function(species_size_age_data){
  "work in progress"
}





####______________________####
####  Species Omission Sensitivity Testing ####





#' @title Perform Leave-One-Out Size Spectra Estimation
#'
#' @description Iterates through the input data for the size spectrum analysis,
#' performing binned size spectrum slope estimation with a single species removed. 
#'
#'
#' @param start_dat Input dataframe that could be passed to `group_l10_spectra()`
#' @param .group_cols Vector of string(s) to use as grouping factors that correspond to column names
#' @param min_weight_g Minimum weight cutoff in grams
#' @param min_l10_bin Minimum 10^x value for the size spectra being measured (>=). Left limit of smallest bin.
#' @param max_l10_bin Maximum 10^x value for the size spectra being measured (<). Left limit of largest bin.
#' @param bin_increment The bin-width on log scale that separates each bin
#'
#' @return
#' @export
#'
#' @examples
species_omit_spectra <- function(
    start_dat = catch_1g_binned,
    .group_cols = c("Year", "survey_area"),
    min_weight_g = 1,
    min_l10_bin = 0,
    max_l10_bin = 4.5,
    bin_increment = 0.5){
  
  
  # 1. Get a list of all the species to iterate through
  spec_list <- start_dat  %>% 
    distinct(comname) %>% 
    pull(comname)
  
  # Set names for that vector
  spec_list <- setNames(spec_list, spec_list)
  
  # 2. For each one, filter comname != species
  # label group as the species omitted
  # Perform group_l10_spectra
  omit_dat <- map_dfr(spec_list, function(species_id){
    
    # Filter the species out
    filtered_dat <- start_dat %>% 
      filter(comname != species_id)
    
    
    # Run the bare essential spectra information for memory sake
    spectra_res <-  group_l10_spectra(
      .group_cols   = .group_cols,
      wmin_grams    = filtered_dat,
      min_weight_g  = min_weight_g,
      min_l10_bin   = min_l10_bin,
      max_l10_bin   = max_l10_bin,
      bin_increment = bin_increment)
    
    # Rename the columns so we can join them in to the data 
    # That had all species more easily
    spectra_res <- spectra_res %>% 
      rename(omit_slope_strat = l10_slope_strat,
             omit_int_strat   = l10_int_strat,
             omit_rsqr_strat  = l10_rsqr_strat,
             omit_sig_strat   = l10_sig_strat)
    
    # Return the results
    return(spectra_res)
    
    # Label what species was left out
  }, .id = "spec_omit")
  
  # Return the list
  return(omit_dat)
  
  
}

# testing
# tar_load(catch_1g_binned)
# Get those sensitivity results
# spec_omit_results <- species_omit_spectra(catch_1g_binned)







species_omit_changes <- function(spec_omit_results,
                                 all_spec_results){
  
  
  # join thew species ommissions to the all species results
  changes_d <- left_join(spec_omit_results, all_spec_results,
                         by = c("group_var", "Year", "season", "survey_area", "hare_group", "decade"))
  
  # Get the variance in the original data for each region
  region_var_df <- all_spec_results %>% 
    select(survey_area, l10_slope_strat, l10_int_strat) %>% 
    group_by(survey_area) %>% 
    summarise(slope_mu = mean(l10_slope_strat, na.rm = T),
              slope_sd = sd(l10_slope_strat, na.rm = T),
              int_mu   = mean(l10_int_strat, na.rm = T),
              int_sd   = sd(l10_int_strat, na.rm = T), 
              .groups = "drop")
  
  # Calculate the difference between
  changes_d <- left_join(changes_d, region_var_df, by = "survey_area") %>% 
    mutate(
      # Subtract the change that omitting a species makes
      slope_shift = l10_slope_strat - omit_slope_strat,
      int_shift   = l10_int_strat - omit_int_strat,
      # standardize the shift with the variance
      slope_shift_z = slope_shift / slope_sd,
      int_shift_z   = int_shift / int_sd
    )
  
  # Return the changes
  return(changes_d)
  
}


# Testing
# Load the slope results for all species, 
# tar_load(nmfs_log10_slopes)
# filter out the groupings well run sensistivity for:
# l10_reg_year <- filter(nmfs_log10_slopes, `group ID` == "single years * region")

# Get Results
# sens_results <- species_omit_changes(spec_omit_results, l10_reg_year)


# # plot a species
# sens_results %>% 
#   mutate(Year = as.numeric(Year)) %>% 
#   filter(spec_omit == "haddock",
#          survey_area == "GoM") %>% 
#   ggplot(aes(Year, l10_slope_strat)) + 
#   geom_point(aes(color = "All Species")) +
#   geom_line(aes(color = "All Species"), group = 1) +
#   geom_point(aes(y = omit_slope_strat, color = "Haddock Removed")) +
#   geom_line(aes(y = omit_slope_strat, color = "Haddock Removed"), group = 1) +
#   facet_wrap(~survey_area, ncol = 1) +
#   theme(legend.position = "bottom") +
#   labs(y = "Biomass Size Spectrum Slope")
# 
# 
# # plot heatmap
# sens_results %>% 
#   mutate(Year = as.numeric(Year)) %>% 
#   filter(survey_area == "GoM") %>% 
#   ggplot(aes(Year, fct_rev(spec_omit), fill = slope_shift_z)) + 
#   geom_tile() +
#   scale_fill_distiller(palette = "RdBu", 
#                        limits = c(-1, 1), 
#                        oob = scales::oob_squish) +
#   facet_wrap(~survey_area, ncol = 1) +
#   theme(legend.position = "bottom") +
#   labs(y = "Species omitted", fill = "Shift in Slope (z-score)")
# 






####______________________####
####  Community Composition ####







####______________________####
####  Phased Out Functions ####



# # Original function to run all the ss groups, manual group creation
# and manual group deconstruction phased out
# all_group_ss_slopes <- function(wmin_grams, min_weight_g, abundance_vals = "stratified"){
#   
#   # Source script
#   # 06_nmfs_ss_group_estimates.R
#   
#   
#   ####__  Set Bodymass Cutoff and Groups
#   
#   # Set bodymass lower limit
#   # Used to filter for lower end of gear selectivity
#   dbin_truncated <- filter(wmin_grams, wmin_g>= min_weight_g) 
#   
#   # Set which function to operate on: stratified or not stratified
#   mle_function <- switch(abundance_vals,
#                          "observed" = group_isd_calc,
#                          "stratified"  = strat_abund_mle_calc)
#   
#   
#   # Set up the factor groupings we want to compare : 
#   
#   #####__ 1.  All years, every region 
#   g1 <- dbin_truncated %>% 
#     mutate(group_level = "all_data") %>% 
#     split(.$group_level) 
#   
#   
#   # get SS results
#   g1_res <- g1 %>% 
#     imap_dfr(mle_function) %>% 
#     mutate(Year   = "all",
#            season = "all",
#            area   = "all")
#   
#   
#   
#   #####__ 2. All Years, each season 
#   
#   # get SS results
#   g2_res <- dbin_truncated  %>% 
#     mutate(group_level = season) %>% 
#     split(.$group_level) %>% 
#     imap_dfr(mle_function) %>% 
#     mutate(Year = "all",
#            season = group_var,
#            area = "all")
#   
#   
#   
#   
#   #####__ 3. All Years, regions  
#   
#   # get SS results
#   g3_res <- dbin_truncated  %>% 
#     mutate(group_level = survey_area) %>% 
#     split(.$group_level) %>% 
#     imap_dfr(mle_function, vecDiff = 2) %>% 
#     mutate(Year = "all",
#            season = "all",
#            area = group_var)
#   
#   
#   
#   
#   #####__ 3. All Years, seasons * regions
#   
#   # get SS results
#   g4_res <- dbin_truncated  %>% 
#     mutate(group_level = str_c(season, survey_area)) %>% 
#     split(.$group_level) %>% 
#     imap_dfr(.f = mle_function) %>% 
#     mutate(Year = "all",
#            season = case_when(
#              str_detect(group_var, "Fall") ~ "Fall",
#              str_detect(group_var, "Spring") ~ "Spring"),
#            area = case_when(
#              str_detect(group_var, "GoM") ~ "GoM",
#              str_detect(group_var, "SNE") ~ "SNE",
#              str_detect(group_var, "MAB") ~ "MAB",
#              str_detect(group_var, "GB") ~ "GB"))
#   
#   
#   #####__ 4. Every year, entire survey
#   
#   # get SS results
#   g5_res <- dbin_truncated  %>% 
#     mutate(group_level = Year) %>% 
#     split(.$group_level) %>% 
#     imap_dfr(.f = mle_function) %>% 
#     mutate(Year = group_var,
#            season = "all",
#            area = "all")
#   
#   
#   #####__ 5. every year, every region
#   
#   # get SS results
#   g6_res <- dbin_truncated  %>% 
#     mutate(group_level = str_c(Year, survey_area)) %>% 
#     split(.$group_level) %>% 
#     imap_dfr(.f = mle_function) %>% 
#     mutate(Year = str_sub(group_var, 1, 4),
#            season = "all",
#            area = str_sub(group_var, 5, -1))
#   
#   #####__ 6. every year, only seasons
#   
#   # get SS results
#   g7_res <- dbin_truncated  %>% 
#     mutate(group_level = str_c(Year, season)) %>% 
#     split(.$group_level) %>% 
#     imap_dfr(.f = mle_function, vecDiff = 2) %>% 
#     mutate(Year = str_sub(group_var, 1, 4),
#            season = str_sub(group_var, 5, -1),
#            area = "all")
#   
#   
#   #####__ 7. every year, region * season
#   
#   
#   # get SS results
#   g8_res <- dbin_truncated  %>% 
#     mutate(group_level = str_c(Year, season, survey_area)) %>% 
#     split(.$group_level) %>% 
#     imap_dfr(.f = mle_function, vecDiff = 2) %>% 
#     mutate(Year = str_sub(group_var, 1, 4),
#            season = case_when(
#              str_detect(group_var, "Fall") ~ "Fall",
#              str_detect(group_var, "Spring") ~ "Spring"),
#            area = case_when(
#              str_detect(group_var, "GoM") ~ "GoM",
#              str_detect(group_var, "SNE") ~ "SNE",
#              str_detect(group_var, "MAB") ~ "MAB",
#              str_detect(group_var, "GB") ~ "GB"))
#   
#   # Put the reults in one table with an ID for how they groups are set up
#   table_complete <- bind_rows(
#     list(
#       "Overall"                        = g1_res,
#       "only seasons"                   = g2_res,
#       "only regions"                   = g3_res,
#       "region * seasons"               = g4_res,
#       "single years"                   = g5_res,
#       "single years * region"          = g6_res,
#       "single years * season"         = g7_res,
#       "single years * season * region" = g8_res), 
#     .id = "group ID")
#   
#   
#   
# }





# # Takes a list of each species, pulls distinct length bins
# # returns a key of what the wmin_gand wmax_g is for each
# # corrects wmin_gand wmax_g to kg for trouble species
# make_bin_key <- function(species_df){
#   
#   #pull the distinct length bins
#   species_df <- species_df %>% 
#     distinct(length_cm, .keep_all = T) %>% 
#     arrange(length_cm)
#   
#   
#   
#   # Add the max length for the bin, and its weight
#   binned_df <- species_df %>% 
#     mutate(
#       lngt_max = length_cm + 1, 
#       wmax_g    = exp(log(LWa) + LWb * log(lngt_max)) * 1000,
#       wmax_g    = ifelse(comname %in% gram_species, wmax_g / 1000, wmax))  %>%
#     select(comname, LWa, LWb, LngtMin = length_cm, wmin_g = ind_weight_g, lngt_max, wmax_g,
#            -c(numlen, Biomass))
#   
#   # return the clean data
#   return(binned_df)}



# # PHASED OUT - Stratified abundance preparation
# # Use to prep data for Individual Size Distribution Plots
# strat_isd_prep <- function(x = dataBin){
#   
#   #tester: x <- strat_year_list[[1]]
#   
#   # arrange by bin lower weight
#   biomass_data <- dplyr::arrange(x, desc(wmin))
#   
#   # Number of fish for the year
#   total_abundance <- sum(biomass_data$strat_total_abund_s)
#   
#   # Have to do not with dplyr:
#   wmin.vec <- biomass_data$wmin_g    # Vector of lower weight bin limits
#   wmax.vec <- biomass_data$wmax_g    # Vector of upper weight bin limits
#   num.vec  <- biomass_data$strat_total_abund_s  # Vector of corresponding counts for those bins
#   
#   # to do a manual count, start with NA's for everything
#   countGTEwmin <- rep(NA, length(num.vec)) 
#   lowCount     <- countGTEwmin
#   highCount    <- countGTEwmin
#   
#   # Loop through
#   for(iii in 1:length(countGTEwmin)) {
#     
#     # Count times wmin.vec >= individual wmin_g, multiply by number of fish for year
#     countGTEwmin[iii]  <- sum( (wmin.vec >= wmin.vec[iii]) * num.vec)
#     
#     # Count times wmin.vec >= individual wmax_g, multiply by number of fish for year
#     lowCount[iii]      <- sum( (wmin.vec >= wmax.vec[iii]) * num.vec)
#     
#     #Count times wmax.vec > individual wmin_g, multiply by number of fish for year
#     highCount[iii]     <- sum( (wmax.vec >  wmin.vec[iii]) * num.vec)
#   }
#   
#   #combine vectors from loop
#   biomass_data = cbind(biomass_data,
#                        "countGTEwmin" = countGTEwmin,
#                        "lowCount"     = lowCount,
#                        "highCount"    = highCount)
#   
#   # format as a tibble
#   biomass_data <- tibble::as_tibble(biomass_data)
#   return(biomass_data)
#   
#   
#   
# }


# # Stratified Abundance Plot
# ggplot_strat_isd <- function(isd_data_prepped, 
#                              mle_results,
#                              plot_rects = TRUE, 
#                              show_pl_fit = TRUE){
#   
#   theme_set(theme_minimal())
#   
#   # Power law parameters and summary details for that year
#   b.MLE     <- mle_results$b
#   total_abundance <- sum(isd_data_prepped$strat_total_abund_s) 
#   b.confMin <- mle_results$confMin
#   b.confMax <- mle_results$confMax
#   
#   # Create range of x values from that year's data
#   # 1. Drop last entry
#   # 2. add it in again as .99999 * what it would be
#   # 3. sandwich the actual last entry on the end...
#   xmin  <- mle_results$xmin
#   xmax  <- mle_results$xmax
#   x.PLB <- seq(xmin, xmax, length = 10000)   # break up the Xlim into 10000 pieces
#   x.PLB.length <- length(x.PLB)              # get the length of that
#   x.PLB <- c(x.PLB[-x.PLB.length], 0.99999 * x.PLB[x.PLB.length], x.PLB[x.PLB.length])
#   
#   
#   # Y values for plot limits/bounds/predictions from bounded power law pdf
#   y.PLB         <- (1 - pPLB(x = x.PLB, b = b.MLE, xmin = min(x.PLB), xmax = max(x.PLB))) * total_abundance
#   y.PLB.confMin <- (1 - pPLB(x = x.PLB, b = b.confMin, xmin = min(x.PLB), xmax = max(x.PLB))) * total_abundance
#   y.PLB.confMax <- (1 - pPLB(x = x.PLB, b = b.confMax, xmin = min(x.PLB), xmax = max(x.PLB))) * total_abundance
#   
#   
#   # Put it in a df to make it easier
#   PLB_df <- data.frame(
#     x.PLB   = x.PLB,
#     y.PLB   = y.PLB,
#     confMin = y.PLB.confMin,
#     confMax = y.PLB.confMax)
#   
#   
#   # Coefficient Labels
#   group_name <- as.character(isd_data_prepped$group_var[1])
#   group_name <- str_replace_all(group_name, "_", " ")
#   the_slope  <- as.character(round(mle_results$b, 3))
#   
#   
#   ####  First Plot  
#   p1 <- ggplot()
#   
#   # Toggle Rectangles
#   if(plot_rects == TRUE) {
#     p1 <- p1 +
#       geom_rect(data = isd_data_prepped, 
#                 aes(xmin = wmin_g, 
#                     xmax = wmax_g, 
#                     ymin = lowCount, 
#                     ymax = highCount),
#                 color = "gray70",
#                 fill = "transparent")}
#   # Add line segments
#   p1 <- p1 + 
#     geom_segment(data = isd_data_prepped,
#                  aes(x = wmin_g, 
#                      xend = wmax_g, 
#                      y = countGTEwmin, 
#                      yend = countGTEwmin),
#                  color = "blue") +
#     scale_x_log10(limits = xlim.global, labels = scales::label_number()) +
#     labs(x = "Body mass (x), g",
#          y =  "Number of values \u2265 x",
#          tag = "(a)",
#          title = group_name, 
#          caption = str_c("b = ", the_slope)) +
#     theme(plot.tag.position = "topright")
#   
#   # Toggle to turn on the fit lines or not
#   if(show_pl_fit == TRUE) {
#     p1 <- p1 +
#       geom_line(data = PLB_df, aes(x.PLB, y.PLB), color = "darkred") +
#       geom_line(data = PLB_df, aes(x.PLB, confMin), color = "darkred", linetype = 2) +
#       geom_line(data = PLB_df, aes(x.PLB, confMax), color = "darkred", linetype = 2)
#     
#   }
#   
#   
#   
#   ####  Second Plot, lm  
#   
#   # Second plot, log10 transformed y scale
#   # lm for size spectra
#   p2 <- ggplot()
#   
#   # Toggle Rectangles
#   if(plot_rects == TRUE) {
#     p2 <- p2 +
#       geom_rect(data = isd_data_prepped, 
#                 aes(xmin = wmin_g, 
#                     xmax = wmax_g, 
#                     ymin = lowCount, 
#                     ymax = highCount),
#                 color = "gray70",
#                 fill = "transparent")}
#   
#   # Add the line segments for biomass bins
#   p2 <- p2 + 
#     geom_segment(data = isd_data_prepped,
#                  aes(x = wmin_g, 
#                      xend = wmax_g, 
#                      y = countGTEwmin, 
#                      yend = countGTEwmin),
#                  color = "blue") +
#     scale_x_log10(limits = xlim.global, labels = scales::label_number()) +
#     labs(x = "Body mass (x), g",
#          y =  "Number of values \u2265 x",
#          tag = "(a)") +
#     theme(plot.tag.position = "topright")
#   
#   
#   
#   # Change the y-axis, log transformation
#   p2 <- p2 +
#     # Add log tranformation on y axis
#     scale_y_log10(labels = scales::label_number()) +
#     labs(y = "Number of values \u2265 x \n(log10 transformed)",
#          tag = "(b)")
#   
#   # Toggle to turn on the fit lines or not
#   if(show_pl_fit == TRUE) {
#     p2 <- p2 +
#       geom_line(data = PLB_df, aes(x.PLB, y.PLB), color = "darkred") +
#       geom_line(data = PLB_df, aes(x.PLB, confMin), color = "darkred", linetype = 2) +
#       geom_line(data = PLB_df, aes(x.PLB, confMax), color = "darkred", linetype = 2)
#   }
#   
#   
#   
#   
#   # Some alterations for the stacked plot
#   p1_stack <- p1 
#   p2_stack <- (p2 + labs(caption = "lm on log tranformation (mizer)"))
#   p3 <- p1_stack / p2_stack
#   
#   the_year  <- as.character(mle_results$Year)
#   plot_list <- list(obs_y   = p1,
#                     log10_y = p2,
#                     stacked = p3)
#   
#   return(plot_list)
#   
# }
# 
# 

