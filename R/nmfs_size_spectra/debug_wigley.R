##




####

analysis_options <- list(
  # 1. Controls the data that reaches the analysis, and is used for context
  # Also Sets min/max for ISD exponent estimates
  min_input_weight_g = 2^0,
  max_input_weight_g = 10^4, # To pick a reasonable sounding limit
  # max_input_weight_g = 2^13, # To match log2 bins
  
  # 2. Set/enforce the bin structure used for spectra analysis
  # These enforce what bins go to binned size spectra analysis
  # min and max set the left limits of the bin range: 
  # i.e. max_l10_bin of 3 = 14^3 up to 10^4
  l10_bin_width = 1,
  min_l10_bin   = 0,
  max_l10_bin   = 3,
  
  # 3. log2 bin limits - LBNbiom method
  log2_bin_width = 1,
  min_log2_bin   = 0,
  max_log2_bin   = 12
)



# run all possible functions to maximize debugging
box_location <- "cloudstorage"
catch_data <- gmri_survdat_prep(
  survdat = NULL, 
  survdat_source = "most recent", 
  box_location = box_location) %>% 
  add_lw_info(
    cutoff = T, 
    box_location = box_location) %>% 
  add_area_stratification(
    include_epu = F, 
    box_location = box_location) %>% 
  fill_func_groups(species_dat = .) 


# Just read it in, no* species filtering
catch_basic <- gmri_survdat_prep(
  survdat = NULL, 
  survdat_source = "most recent", 
  box_location = box_location)


# Totls
catch_basic %>% 
  distinct(cruise6, station, tow, est_towdate, season, comname, biomass_kg) %>% 
  pull(biomass_kg) %>% sum()

# just the lw species added
catch_nostrat <- catch_basic %>% 
  add_lw_info(
    cutoff = F, 
    box_location = box_location)

# Weight of species with wigley coef
catch_nostrat %>% 
  distinct(cruise6, station, tow, est_towdate, season, comname, biomass_kg) %>% 
  pull(biomass_kg) %>% sum()





















# ----= RUN MLE

#' @title {MLE Size Spectra Estimation}
#'
#' @param ss_input Dataframe containing a column of abundance and a column of sizes
#' @param grouping_vars string identifiers of columns to group_by prior to analysis
#' @param abundance_vals string indicating column of abundances
#' @param weight_vals string indicating column with individual weight data
#' @param isd_xmin lower limit for size distribution fitting
#' @param isd_xmax upper limit for size distribution fitting
#'
#' @return
#' @export
#'
#' @examples
group_mle_estimates <-  function(
    ss_input, 
    grouping_vars, 
    abundance_vals = "numlen_adj",
    weight_vals = "ind_weight_g",
    isd_xmin = NULL,
    isd_xmax = NULL,
    global_min = TRUE,
    global_max = TRUE){
  
  
  # 1. toggle which abundance/weight columns to use with switch
  .abund_col  <- sym(abundance_vals)
  .weight_col <- sym(weight_vals)
  .group_cols <- grouping_vars
  .agg_cols   <- c(.group_cols, "comname", weight_vals)
  
  
  # 2. Select only the columns we need:
  # Now aggregate by body-weight within the groups we are measuring spectra for:
  ss_input_summs <- ss_input %>% 
    group_by(!!!syms(.agg_cols)) %>% 
    summarise(
      Number = sum(!!.abund_col),
      LWa = unique(a), 
      LWb = unique(b),
      .groups = "drop")  %>% 
    # Create a group column that we can loop through
    unite(col = "group_var", {{.group_cols}}, sep = "-", remove = FALSE, na.rm = FALSE)
  
  
  
  # Drop columns we don't need, rename to match edwards code
  # edwards calls this df databinforlike
  mle_input <- ss_input_summs  %>% 
    select(
      group_var, 
      Number, 
      bodyMass = !!.weight_col)
  
  
  
  #------ Set all-group constants
  # set these one time
  
  # 3. Set Global power-law limits:
  # set left bound & right bounds, grams
  if(global_min == TRUE){
    if(is.null(isd_xmin)){ isd_xmin <- min(mle_input$bodyMass, na.rm = T)}
  }
  if(global_max == TRUE){
    if(is.null(isd_xmax)){ isd_xmax <- max(mle_input$bodyMass, na.rm = T)}
  }
  
  
  
  # Set global controls if we're doing that:
  if(global_min == TRUE & global_max == TRUE){
    
    # Filter the range of sizes we want to include
    mle_input <- mle_input %>% 
      filter(bodyMass >= isd_xmin,
             bodyMass <= isd_xmax)
    
    # # Get many total individuals
    n <- sum(ceiling(mle_input$Number) )
    
    # Vector of body weights, repeat bodymass x abundance
    x <- rep(mle_input$bodyMass, ceiling(mle_input$Number))
    
    # Sum( log(bodymass) )
    sum_log_x <- sum( log(x) )
    
    # Analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
    PL.bMLE_global  <- 1/( log(min(x)) - sum_log_x / length(x)) - 1
  }
  
  
  #------ Loop through Groups 
  group_results_df <- mle_input %>% 
    split(.$group_var) %>% 
    map_df(
      function(ss_input_i){
        
        
        # 3. Tune subgroup the power-law limits:
        # set left bound & right bounds, grams
        min_i <- min(ss_input_i$bodyMass, na.rm = T)
        max_i <- max(ss_input_i$bodyMass, na.rm = T)
        if(global_min == FALSE){
          if(is.null(isd_xmin)){ isd_xmin <- min_i}
        }
        if(global_max == FALSE){
          if(is.null(isd_xmax)){ isd_xmax <- max_i}
        }
        
        # Filter the range of sizes we want to include
        ss_input_i <- ss_input_i %>% 
          filter(bodyMass >= isd_xmin,
                 bodyMass <= isd_xmax)
        
        
        # Get the subgroup inputs:
        # Vector of body weights, repeat bodymass x abundance
        x_i <- rep(ss_input_i$bodyMass, ceiling(ss_input_i$Number))
        
        # Get many total individuals
        n_i <- sum(ceiling(ss_input_i$Number) )
        
        # Sum( log(bodymass) )
        sum_log_xi <- sum( log(x_i)) 
        
        # Analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
        PL.bMLE_i  <- 1/( log(min(x_i)) - sum_log_xi / length(x_i)) - 1
        
        # Toggle which starting value to use
        PL.bMLE <- ifelse(global_max == T, PL.bMLE_global, PL.bMLE_i)
        
        
        
        # Do the exponent estimation for group i
        group_est <- calcLike(
          negLL.fn = negLL.PLB,  # MLE function, takes a vector of weights
          p        = PL.bMLE, 
          x        = x_i,
          xmin     = isd_xmin, 
          xmax     = isd_xmax,
          n        = n_i,
          sumlogx  = sum_log_xi)
        
        
        
        # Put it into a dataframe to rejoin neatly
        mle_group_results <- data.frame(
          xmin_fit = isd_xmin,
          xmax_fit = isd_xmax,
          xmin_actual = min_i,
          xmax_actual = max_i,
          n = n_i,
          b = group_est$MLE,
          confMin = group_est$conf[1],
          confMax = group_est$conf[2])
        
        
        # Process C and standard error
        mle_group_results <- mle_group_results %>% 
          mutate(
            stdErr = (abs(confMin - b) + abs(confMax - b)) / (2 * 1.96),
            C = (b != -1 ) * (b + 1) / ( xmax_fit^(b + 1) - xmin_fit^(b + 1) ) + (b == -1) * 1 / ( log(xmax_fit) - log(xmin_fit)))
        return(mle_group_results)
      }, 
      .id = "group_var") %>% 
    separate(
      group_var, 
      sep = "-", 
      into = grouping_vars, 
      remove = F)
  
  
  
  #---
  #  Make a table of constituent combinations
  #  Adding in any missing groups as "all data"
  grouping_table <- group_results_df %>% 
    distinct(group_var, !!!syms(.group_cols))
  grouping_table <- add_missing_groups(
    group_dataframe = grouping_table)
  
  
  # Merge in the group details
  group_results_df <- full_join(
    grouping_table, group_results_df, join_by(group_var, !!!syms(.group_cols))) %>% 
    mutate(across(c(group_var, .group_cols), as.character))
  
  #---
  
  
  # Spit it out
  return(group_results_df)
  
}


season_res <- catch_complete  %>%
  mutate(ind_weight_g = ind_weight_kg * 1000) %>% 
  group_mle_estimates(
    ss_input = .,
    isd_xmin = 1,
    isd_xmax = 10000,
    abundance_vals = "numlen_adj",
    weight_vals = "ind_weight_g",
    grouping_vars =  c("Year", "season", "survey_area"),
    global_min = TRUE,
    global_max = TRUE)