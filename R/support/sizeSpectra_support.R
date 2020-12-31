####  sizeSpectra Support Functions  ####

####_____________________####
####  Functions  ####


library(sizeSpectra)
library(tidyverse)



# Takes a list of each species, pulls distinct length bins
# returns a key of what the wmin and wmax is for each
# corrects wmin and wmax to kg for trouble species
make_bin_key <- function(species_df){
  
  #pull the distinct length bins
  species_df <- species_df %>% 
    distinct(LngtClass, .keep_all = T) %>% 
    arrange(LngtClass)
  
  
  
  # Add the max length for the bin, and its weight
  binned_df <- species_df %>% 
    mutate(
      LngtMax = LngtClass + 1, 
      wmax    = exp(log(LWa) + LWb * log(LngtMax)) * 1000,
      # Correction for the species that were in grams
      wmax    = ifelse(SpecCode %in% gram_species, wmax / 1000, wmax))  %>%
    select(SpecCode, LWa, LWb, LngtMin = LngtClass, wmin = bodyMass, LngtMax, wmax,
           -c(Number, Biomass))
  
  # return the clean data
  return(binned_df)}





####_____________________####
#### Survey Abundance MLE  ####

# Takes databin split by a grouping variable .$group_var
# returns the power law coefficients for that group
group_mle_calc <- function(dataBinForLike, group_var, vecDiff = 0.5){
  
  # Select the right columns
  dataBinForLike = dplyr::select(dataBinForLike,
                                 SpecCode,
                                 wmin,
                                 wmax,
                                 Number)
  
  # Set n, xmin, xmax
  n    = sum(dataBinForLike$Number)
  xmin = min(dataBinForLike$wmin)
  xmax = max(dataBinForLike$wmax)
  
  
  
  
  
  
  
  
  # Get the likelihood calculation for the bins
  # previously named MLEbins.nSeaFung.new
  mle_group_bins = calcLike(negLL.fn = negLL.PLB.bins.species,
                            p = -1.9,
                            vecDiff = vecDiff,
                            suppress.warnings = TRUE,
                            dataBinForLike = dataBinForLike,
                            n = n,
                            xmin = xmin,
                            xmax = xmax)
  
  
  # Store outputs in a dataframe
  mle_group_bins = data.frame(group_var = group_var,
                              xmin = xmin,
                              xmax = xmax,
                              n = n,
                              b = mle_group_bins$MLE,
                              confMin = mle_group_bins$conf[1],
                              confMax = mle_group_bins$conf[2]) %>% 
    mutate(stdErr = (abs(confMin - b) + abs(confMax - b)) / (2 * 1.96),
           C = (b != -1 ) * (b + 1) / ( xmax^(b + 1) - xmin^(b + 1) ) + (b == -1) * 1 / ( log(xmax) - log(xmin)))
  
  
  
  # Spit out the final output
  return(mle_group_bins)
}






#plotting group maximum likelihood size specttrum results
group_mle_plot <- function(mle_res){
  plot <- mle_res %>% 
    ggplot(aes(Year, b, color = area, shape = season)) +
    geom_pointrange(aes(x = Year, y = b, ymin = confMin, ymax = confMax),
                    alpha = 0.6) +
    labs(x = "",
         y = "Size Spectrum Slope (b)") +
    theme(axis.text.x = element_text(angle = 90, size = 6, vjust = 0.5))
  
  return(plot)
}



# Function to prepare data for Individual Slope Plots
# Takes a list of the data split up by the grouping variable
# Gets the number of fishes within each of the size class bins
# Returns the table, with columns for GTE the wmin, a low count, and a high count for each bin
isd_plot_prep <- function(x = dataBin){
  
  # arrange by bin lower weight
  data.year <- dplyr::arrange(x, desc(wmin))
  
  # Number of fish for the year
  sumNumber <- sum(data.year$Number)
  
  # Have to do not with dplyr:
  wmin.vec <- data.year$wmin    # Vector of lower weight bin limits
  wmax.vec <- data.year$wmax    # Vector of upper weight bin limits
  num.vec  <- data.year$Number  # Vector of corresponding counts for those bins
  
  # to do a manual count, start with NA's for everything
  countGTEwmin <- rep(NA, length(num.vec)) 
  lowCount     <- countGTEwmin
  highCount    <- countGTEwmin
  
  # Loop through
  for(iii in 1:length(countGTEwmin)) {
    # Count times wmin.vec >= individual wmin, multiply by number of fish for year
    countGTEwmin[iii]  <- sum( (wmin.vec >= wmin.vec[iii]) * num.vec)
    
    # Count times wmin.vec >= individual wmax, multiply by number of fish for year
    lowCount[iii]      <- sum( (wmin.vec >= wmax.vec[iii]) * num.vec)
    
    #Count times wmax.vec > individual wmin, multiply by number of fish for year
    highCount[iii]     <- sum( (wmax.vec >  wmin.vec[iii]) * num.vec)
  }
  
  #combine vectors from loop
  data.year = cbind(data.year,
                    "countGTEwmin" = countGTEwmin,
                    "lowCount"     = lowCount,
                    "highCount"    = highCount)
  
  # format as a tibble
  data.year <- tibble::as_tibble(data.year)
  return(data.year)
  
  
  
}





# Function for ISD plots 
# should work on whatever the grouping variable is as long as the lists match
ggplot_isd <- function(data_group, mle.bins.res, plot_rects = TRUE, show_pl_fit = TRUE){
  
  theme_set(theme_minimal())
  
  # Power law parameters and summary details for that year
  b.MLE     <- mle.bins.res$b
  sumNumber <- sum(data_group$Number) 
  b.confMin <- mle.bins.res$confMin
  b.confMax <- mle.bins.res$confMax
  
  # Create range of x values from that year's data
  # 1. Drop last entry
  # 2. add it in again as .99999 * what it would be
  # 3. sandwich the actual last entry on the end...
  xmin  <- mle.bins.res$xmin
  xmax  <- mle.bins.res$xmax
  x.PLB <- seq(xmin, xmax, length = 10000)   # break up the Xlim into 10000 pieces
  x.PLB.length <- length(x.PLB)              # get the length of that
  x.PLB <- c(x.PLB[-x.PLB.length], 0.99999 * x.PLB[x.PLB.length], x.PLB[x.PLB.length])

  
  # Y values for plot limits/bounds/predictions from bounded power law pdf
  y.PLB         <- (1 - pPLB(x = x.PLB, b = b.MLE, xmin = min(x.PLB), xmax = max(x.PLB))) * sumNumber
  y.PLB.confMin <- (1 - pPLB(x = x.PLB, b = b.confMin, xmin = min(x.PLB), xmax = max(x.PLB))) * sumNumber
  y.PLB.confMax <- (1 - pPLB(x = x.PLB, b = b.confMax, xmin = min(x.PLB), xmax = max(x.PLB))) * sumNumber
  
  
  # Put it in a df to make it easier
  PLB_df <- data.frame(
    x.PLB = x.PLB,
    y.PLB = y.PLB,
    confMin = y.PLB.confMin,
    confMax = y.PLB.confMax)
  
  
  # Coefficient Labels
  group_name <- as.character(data_group$group_var[1])
  group_name <- str_replace_all(group_name, "_", " ")
  the_slope <- as.character(round(mle.bins.res$b, 3))
  
  
  ####  First Plot  
  p1 <- ggplot()
  
  # Toggle Rectangles
  if(plot_rects == TRUE) {
    p1 <- p1 +
    geom_rect(data = data_group, 
              aes(xmin = wmin, xmax = wmax, ymin = lowCount, ymax = highCount),
              color = "gray70",
              fill = "transparent")}
  p1 <- p1 + 
    geom_segment(data = data_group,
                 aes(x = wmin, xend = wmax, y = countGTEwmin, yend = countGTEwmin),
                 color = "blue") +
    scale_x_log10(limits = xlim.global, labels = scales::label_number()) +
    labs(x = "Body mass (x), g",
         y =  "Number of values \u2265 x",
         tag = "(a)",
         title = group_name, 
         caption = str_c("b = ", the_slope)) +
    theme(plot.tag.position = "topright")
    
  # Toggle to turn on the fit lines or not
  if(show_pl_fit == TRUE) {
    p1 <- p1 +
      geom_line(data = PLB_df, aes(x.PLB, y.PLB), color = "darkred") +
      geom_line(data = PLB_df, aes(x.PLB, confMin), color = "darkred", linetype = 2) +
      geom_line(data = PLB_df, aes(x.PLB, confMax), color = "darkred", linetype = 2)
      
  }
    
  
  
  ####  Second Plot, lm  
 
   # Second plot, log10 transformed y scale
   # lm for size spectra
  p2 <- ggplot()
  
  # Toggle Rectangles
  if(plot_rects == TRUE) {
    p2 <- p2 +
      geom_rect(data = data_group, 
                aes(xmin = wmin, xmax = wmax, ymin = lowCount, ymax = highCount),
                color = "gray70",
                fill = "transparent")}
  
  # Add the line segments for biomass bins
  p2 <- p2 + 
    geom_segment(data = data_group,
                 aes(x = wmin, xend = wmax, y = countGTEwmin, yend = countGTEwmin),
                 color = "blue") +
    scale_x_log10(limits = xlim.global, labels = scales::label_number()) +
    labs(x = "Body mass (x), g",
         y =  "Number of values \u2265 x",
         tag = "(a)") +
    theme(plot.tag.position = "topright")
  
  
  
  # Change the y-axis, log transformation
  p2 <- p2 +
    # Add log tranformation on y axis
    scale_y_log10(labels = scales::label_number()) +
    labs(y = "Number of values \u2265 x \n(log10 transformed)",
         tag = "(b)")
  
  # Toggle to turn on the fit lines or not
  if(show_pl_fit == TRUE) {
    p2 <- p2 +
      # geom_smooth(data = data_group,
      #             aes(x = wmin, 
      #                 y = lowCount),
      #             formula = y ~ x,
      #             method = "lm", 
      #             se = T,
      #             color = "darkred")
      # Adding the power law fit here also:
      geom_line(data = PLB_df, aes(x.PLB, y.PLB), color = "darkred") +
      geom_line(data = PLB_df, aes(x.PLB, confMin), color = "darkred", linetype = 2) +
      geom_line(data = PLB_df, aes(x.PLB, confMax), color = "darkred", linetype = 2)
    
  }
  
  
  # Some alterations for the stacked plot
  p1_stack <- p1 
  p2_stack <- (p2 + labs(caption = "lm on log tranformation (mizer)"))
  p3 <- p1_stack / p2_stack
  
  the_year  <- as.character(mle.bins.res$Year)
  plot_list <- list(obs_y   = p1,
                    log10_y = p2,
                    stacked = p3)
  
  return(plot_list)
}







####_____________________####
# Stratified abundance MLE  ####



# Stratified abundance MLE calculation
strat_abund_mle_calc <- function(dataBinForLike, group_var, vecDiff = 0.5){
  
  # Select the right columns
  dataBinForLike = dplyr::select(dataBinForLike,
                                 SpecCode,
                                 wmin,
                                 wmax,
                                 Number = expanded_abund_s)
  
  # Set n, xmin, xmax
  n    = sum(dataBinForLike$Number)
  xmin = min(dataBinForLike$wmin)
  xmax = max(dataBinForLike$wmax)
  
  
  
  # Get the likelihood calculation for the bins
  # previously named MLEbins.nSeaFung.new
  mle_group_bins = calcLike(negLL.fn = negLL.PLB.bins.species,
                            p = -1.9,
                            vecDiff = vecDiff,
                            suppress.warnings = TRUE,
                            dataBinForLike = dataBinForLike,
                            n = n,
                            xmin = xmin,
                            xmax = xmax)
  
  
  # Store outputs in a dataframe
  mle_group_bins = data.frame(group_var = group_var,
                              xmin = xmin,
                              xmax = xmax,
                              n = n,
                              b = mle_group_bins$MLE,
                              confMin = mle_group_bins$conf[1],
                              confMax = mle_group_bins$conf[2]) %>% 
    mutate(stdErr = (abs(confMin - b) + abs(confMax - b)) / (2 * 1.96),
           C = (b != -1 ) * (b + 1) / ( xmax^(b + 1) - xmin^(b + 1) ) + (b == -1) * 1 / ( log(xmax) - log(xmin)))
  
  
  
  # Spit out the final output
  return(mle_group_bins)
}





# Stratified abundance preparation
# Use to prep data for Individual Size Distribution Plots
strat_isd_prep <- function(x = dataBin){
  
  #tester: x <- strat_year_list[[1]]
  
  # arrange by bin lower weight
  data.year <- dplyr::arrange(x, desc(wmin))
  
  # Number of fish for the year
  sumNumber <- sum(data.year$expanded_abund_s)
  
  # Have to do not with dplyr:
  wmin.vec <- data.year$wmin    # Vector of lower weight bin limits
  wmax.vec <- data.year$wmax    # Vector of upper weight bin limits
  num.vec  <- data.year$expanded_abund_s  # Vector of corresponding counts for those bins
  
  # to do a manual count, start with NA's for everything
  countGTEwmin <- rep(NA, length(num.vec)) 
  lowCount     <- countGTEwmin
  highCount    <- countGTEwmin
  
  # Loop through
  for(iii in 1:length(countGTEwmin)) {
    # Count times wmin.vec >= individual wmin, multiply by number of fish for year
    countGTEwmin[iii]  <- sum( (wmin.vec >= wmin.vec[iii]) * num.vec)
    
    # Count times wmin.vec >= individual wmax, multiply by number of fish for year
    lowCount[iii]      <- sum( (wmin.vec >= wmax.vec[iii]) * num.vec)
    
    #Count times wmax.vec > individual wmin, multiply by number of fish for year
    highCount[iii]     <- sum( (wmax.vec >  wmin.vec[iii]) * num.vec)
  }
  
  #combine vectors from loop
  data.year = cbind(data.year,
                    "countGTEwmin" = countGTEwmin,
                    "lowCount"     = lowCount,
                    "highCount"    = highCount)
  
  # format as a tibble
  data.year <- tibble::as_tibble(data.year)
  return(data.year)
  
  
  
}




# Stratified Abundance Plot
ggplot_strat_isd <- function(data_group, group_ss_res, plot_rects = TRUE, show_pl_fit = TRUE){
  
    theme_set(theme_minimal())
    
    # Power law parameters and summary details for that year
    b.MLE     <- group_ss_res$b
    sumNumber <- sum(data_group$expanded_abund_s) 
    b.confMin <- group_ss_res$confMin
    b.confMax <- group_ss_res$confMax
    
    # Create range of x values from that year's data
    # 1. Drop last entry
    # 2. add it in again as .99999 * what it would be
    # 3. sandwich the actual last entry on the end...
    xmin  <- group_ss_res$xmin
    xmax  <- group_ss_res$xmax
    x.PLB <- seq(xmin, xmax, length = 10000)   # break up the Xlim into 10000 pieces
    x.PLB.length <- length(x.PLB)              # get the length of that
    x.PLB <- c(x.PLB[-x.PLB.length], 0.99999 * x.PLB[x.PLB.length], x.PLB[x.PLB.length])
    
    
    # Y values for plot limits/bounds/predictions from bounded power law pdf
    y.PLB         <- (1 - pPLB(x = x.PLB, b = b.MLE, xmin = min(x.PLB), xmax = max(x.PLB))) * sumNumber
    y.PLB.confMin <- (1 - pPLB(x = x.PLB, b = b.confMin, xmin = min(x.PLB), xmax = max(x.PLB))) * sumNumber
    y.PLB.confMax <- (1 - pPLB(x = x.PLB, b = b.confMax, xmin = min(x.PLB), xmax = max(x.PLB))) * sumNumber
    
    
    # Put it in a df to make it easier
    PLB_df <- data.frame(
      x.PLB   = x.PLB,
      y.PLB   = y.PLB,
      confMin = y.PLB.confMin,
      confMax = y.PLB.confMax)
    
    
    # Coefficient Labels
    group_name <- as.character(data_group$group_var[1])
    group_name <- str_replace_all(group_name, "_", " ")
    the_slope  <- as.character(round(group_ss_res$b, 3))
    
    
    ####  First Plot  
    p1 <- ggplot()
    
    # Toggle Rectangles
    if(plot_rects == TRUE) {
      p1 <- p1 +
        geom_rect(data = data_group, 
                  aes(xmin = wmin, 
                      xmax = wmax, 
                      ymin = lowCount, 
                      ymax = highCount),
                  color = "gray70",
                  fill = "transparent")}
    # Add line segments
    p1 <- p1 + 
      geom_segment(data = data_group,
                   aes(x = wmin, 
                       xend = wmax, 
                       y = countGTEwmin, 
                       yend = countGTEwmin),
                   color = "blue") +
      scale_x_log10(limits = xlim.global, labels = scales::label_number()) +
      labs(x = "Body mass (x), g",
           y =  "Number of values \u2265 x",
           tag = "(a)",
           title = group_name, 
           caption = str_c("b = ", the_slope)) +
      theme(plot.tag.position = "topright")
    
    # Toggle to turn on the fit lines or not
    if(show_pl_fit == TRUE) {
      p1 <- p1 +
        geom_line(data = PLB_df, aes(x.PLB, y.PLB), color = "darkred") +
        geom_line(data = PLB_df, aes(x.PLB, confMin), color = "darkred", linetype = 2) +
        geom_line(data = PLB_df, aes(x.PLB, confMax), color = "darkred", linetype = 2)
      
    }
    
    
    
    ####  Second Plot, lm  
    
    # Second plot, log10 transformed y scale
    # lm for size spectra
    p2 <- ggplot()
    
    # Toggle Rectangles
    if(plot_rects == TRUE) {
      p2 <- p2 +
        geom_rect(data = data_group, 
                  aes(xmin = wmin, 
                      xmax = wmax, 
                      ymin = lowCount, 
                      ymax = highCount),
                  color = "gray70",
                  fill = "transparent")}
    
    # Add the line segments for biomass bins
    p2 <- p2 + 
      geom_segment(data = data_group,
                   aes(x = wmin, 
                       xend = wmax, 
                       y = countGTEwmin, 
                       yend = countGTEwmin),
                   color = "blue") +
      scale_x_log10(limits = xlim.global, labels = scales::label_number()) +
      labs(x = "Body mass (x), g",
           y =  "Number of values \u2265 x",
           tag = "(a)") +
      theme(plot.tag.position = "topright")
    
    
    
    # Change the y-axis, log transformation
    p2 <- p2 +
      # Add log tranformation on y axis
      scale_y_log10(labels = scales::label_number()) +
      labs(y = "Number of values \u2265 x \n(log10 transformed)",
           tag = "(b)")
    
    # Toggle to turn on the fit lines or not
    if(show_pl_fit == TRUE) {
      p2 <- p2 +
        # geom_smooth(data = data_group,
        #             aes(x = wmin, 
        #                 y = lowCount),
        #             formula = y ~ x,
        #             method = "lm", 
        #             se = T,
        #             color = "darkred")
        # Adding the power law fit here also:
        geom_line(data = PLB_df, aes(x.PLB, y.PLB), color = "darkred") +
        geom_line(data = PLB_df, aes(x.PLB, confMin), color = "darkred", linetype = 2) +
        geom_line(data = PLB_df, aes(x.PLB, confMax), color = "darkred", linetype = 2)
    }
    
    
    
    
    # Some alterations for the stacked plot
    p1_stack <- p1 
    p2_stack <- (p2 + labs(caption = "lm on log tranformation (mizer)"))
    p3 <- p1_stack / p2_stack
    
    the_year  <- as.character(group_ss_res$Year)
    plot_list <- list(obs_y   = p1,
                      log10_y = p2,
                      stacked = p3)
    
    return(plot_list)
  
}
