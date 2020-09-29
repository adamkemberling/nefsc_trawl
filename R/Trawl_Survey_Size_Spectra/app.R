
#################  Packages  #################
library(shiny)
library(shinymaterial)
library(here)
library(gmRi)
library(janitor)
library(sizeSpectra)
library(patchwork)
library(tidyverse)

####  Support Functions  ####
source(here("R/support/sizeSpectra_support.R"))


#################  Data  #################

# Load the Catch of each station, with the length and weight bins for each species
# there has been no biomass cutoff at this stage
nefsc_master <- read_csv(here("data/NEFSC/nefsc_ss_bins.csv"), col_types = cols()) #%>% filter(Year >= 2000)
menh_master <- read_csv(here("data/MENH/menh_ss_bins.csv"), col_types = cols())


####  Attach CSS, Headers, Footers  ####
gmRi::use_gmri_stylesheets(stylesheet = "gmri rmarkdown", 
                           header = "gmri logo right", 
                           footer = "akemberling")





#################  Plot Pre-Processing  #################


####__ Static NEFSC Plots  ####


# Set the mass cutoff in grams for the static plots
mass_cutoff <- 400

# Filter biomass by cutoff
nefsc_dbin <- nefsc_master %>% 
    filter(wmin >= mass_cutoff) %>% 
    mutate(group_var = str_c(Year, season, sep = "_"))



# # Map through instead of looping
# nefsc_400g_ss <- nefsc_dbin %>% 
#     split(.$group_var) %>% 
#     imap_dfr(group_mle_calc) %>% 
#     mutate(
#         stdErr = (abs(confMin - b) + abs(confMax - b)) / (2 * 1.96),
#         Year = str_sub(group_var, 1, 4),
#         Year = as.numeric(as.character(Year)),
#         season = str_sub(group_var, 6, -1),
#         season = factor(season, levels = c("SPRING", "SUMMER", "FALL", "WINTER")),
#         C = (b != -1 ) * (b + 1) / ( xmax^(b + 1) - xmin^(b + 1) ) + (b == -1) * 1 / ( log(xmax) - log(xmin)))


# Save the ss results and just load them here:
# write_csv(nefsc_400g_ss, here("data/NEFSC/nefsc_400gss_results.csv"))


# Load the pre-processed 400g data
nefsc_400g_ss <- read_csv(here("data/NEFSC/nefsc_400gss_results.csv"),
                          col_types = cols())


# Data inputs
dbin_group_list <- nefsc_dbin %>% 
    filter(wmin >= 400) %>% # filter 400g to be consistent with the data prep
    split(.$group_var) %>% 
    map(isd_plot_prep)


# global limits for plot
xlim.global <- c( min(nefsc_dbin$wmin), max(nefsc_dbin$wmax) )

# Vector of years, for naming
group_names <- names(dbin_group_list)
names(group_names) <- names(dbin_group_list)

# Size Spectra Coefficients for each group
ss_results_list  <- nefsc_400g_ss %>% split(.$group_var)      



# Map through the groups for each plot
seasonal_isd <- map2(dbin_group_list, ss_results_list, .f = ggplot_isd) %>% setNames(group_names)


# Put them into a list by year
nefsc_years <- as.character(min(nefsc_dbin$Year):max(nefsc_dbin$Year))
nefsc_static_plots <- map(nefsc_years, function(the_year){
    
    # Index which list entries match the year
    list_key <- which(str_detect(names(seasonal_isd), the_year)) 
    
    # Take length for assembly of patchwork
    num_plots <- length(list_key) 
    
    if(num_plots == 1){
        p1 <- seasonal_isd[[list_key[1]]][["obs_y"]]
        plot_patchwork <- p1
    } else if(num_plots == 2){
        p1 <- seasonal_isd[[list_key[1]]][["obs_y"]]
        p2 <- seasonal_isd[[list_key[2]]][["obs_y"]]
        plot_patchwork <- p1 | p2
    } else if(num_plots == 3) {
        p1 <- seasonal_isd[[list_key[1]]][["obs_y"]]
        p2 <- seasonal_isd[[list_key[2]]][["obs_y"]]
        p3 <- seasonal_isd[[list_key[3]]][["obs_y"]]
        plot_patchwork <- wrap_plots(p1, p2, p3, ncol =2)
    } else if(num_plots == 4) {
        p1 <- seasonal_isd[[list_key[1]]][["obs_y"]]
        p2 <- seasonal_isd[[list_key[2]]][["obs_y"]]
        p3 <- seasonal_isd[[list_key[3]]][["obs_y"]]
        p4 <- seasonal_isd[[list_key[4]]][["obs_y"]]
        plot_patchwork <- wrap_plots(p1, p2, p3, p4, ncol = 2)
    }
    
    return(plot_patchwork)
}) %>% setNames(nefsc_years)







####__ Static MENH Plots  ####

# Filter biomass by cutoff
menh_dbin <- menh_master %>% 
    filter(wmin >= mass_cutoff) %>% 
    mutate(group_var = str_c(Year, season, sep = "_"))
# 
# 
# # Map through instead of looping
# menh_400g_ss <- menh_dbin %>% 
#     split(.$group_var) %>% 
#     imap_dfr(group_mle_calc) %>% 
#     mutate(
#         stdErr = (abs(confMin - b) + abs(confMax - b)) / (2 * 1.96),
#         Year = str_sub(group_var, 1, 4),
#         Year = as.numeric(as.character(Year)),
#         season = str_sub(group_var, 6, -1),
#         season = factor(season, levels = c("Spring", "Fall")),
#         C = (b != -1 ) * (b + 1) / ( xmax^(b + 1) - xmin^(b + 1) ) + (b == -1) * 1 / ( log(xmax) - log(xmin)))
# 
# 
# # # Save that out so we don't need to wait on the pre-processing
# write_csv(menh_400g_ss, here("data/MENH/menh_400gss_results.csv"))

# Load 400g pre-processed results
menh_400g_ss <- read_csv(here("data/MENH/menh_400gss_results.csv"),
                         col_types = cols())

# Data inputs
dbin_group_list <- menh_dbin %>% 
    filter(wmin >= 400) %>% # filter 400g to be consistent with the data prep
    split(.$group_var) %>% 
    map(isd_plot_prep)


# global limits for plot
xlim.global <- c( min(menh_dbin$wmin), max(menh_dbin$wmax) )

# Vector of years, for naming
group_names <- names(dbin_group_list)
names(group_names) <- names(dbin_group_list)

# Size Spectra Coefficients for each group
ss_results_list  <- menh_400g_ss %>% split(.$group_var)      


# Map through the groups for each plot
seasonal_isd <- map2(dbin_group_list, ss_results_list, .f = ggplot_isd) %>% setNames(group_names)


# Put them into a list by year
menh_years <- as.character(min(menh_dbin$Year):max(menh_dbin$Year))
menh_static_plots <- map(menh_years, function(the_year){
    
    # Index which list entries match the year
    list_key <- which(str_detect(names(seasonal_isd), the_year)) 
    
    # Take length for assembly of patchwork
    num_plots <- length(list_key) 
    
    if(num_plots == 1){
        p1 <- seasonal_isd[[list_key[1]]][["obs_y"]]
        plot_patchwork <- p1
    } else if(num_plots == 2){
        p1 <- seasonal_isd[[list_key[1]]][["obs_y"]]
        p2 <- seasonal_isd[[list_key[2]]][["obs_y"]]
        plot_patchwork <- p1 | p2
    } else if(num_plots == 3) {
        p1 <- seasonal_isd[[list_key[1]]][["obs_y"]]
        p2 <- seasonal_isd[[list_key[2]]][["obs_y"]]
        p3 <- seasonal_isd[[list_key[3]]][["obs_y"]]
        plot_patchwork <- wrap_plots(p1, p2, p3, ncol =2)
    } else if(num_plots == 4) {
        p1 <- seasonal_isd[[list_key[1]]][["obs_y"]]
        p2 <- seasonal_isd[[list_key[2]]][["obs_y"]]
        p3 <- seasonal_isd[[list_key[3]]][["obs_y"]]
        p4 <- seasonal_isd[[list_key[4]]][["obs_y"]]
        plot_patchwork <- wrap_plots(p1, p2, p3, p4, ncol = 2)
    }
    
    return(plot_patchwork)
}) %>% setNames(menh_years)




#################  User Interface  #####################
ui <- material_page(
    nav_bar_fixed = TRUE,
    primary_theme_color = "#00695c", 
    secondary_theme_color = "#00796b",
    
    # Application title
    title = " Comparing Trawl Survey Size Spectra", 

    
    ####  Define tabs  ####
    material_tabs(
        tabs = c(
            "NEFSC Overall" = "first_tab",
            "NEFSC Annual"  = "second_tab",
            "MENH Overall"  = "third_tab",
            "MENH Annual"   = "fourth_tab"
        )
    ),


    ####__ Tab 1 Content - NEFSC  ####
    material_tab_content(
        tab_id = "first_tab",
        material_card(
            # Card Title
            title = "Northeast Groundfish Survey | Overall",
            
            # Selection Controls:
            tags$h6("Step 1. - Choose a Minimum Biomass"),
            material_number_box(input_id = "biomass_cutoff",
                           label = "Minimum Fish Size (g) for Inclusion:", 
                           min_value = 10,
                           max_value = 800, 
                           step_size = 10, 
                           initial_value = 400),
            tags$h6("Step 2. - Choose the Reference Period for the Calculations:"),
            material_slider(input_id = "nefsc_start_year", 
                            label = "Start Year for Estimating Summary Metrics:", 
                            min_value = 1963,
                            max_value = 2018, 
                            step_size = 1, 
                            initial_value = 2000),
            material_slider(input_id = "nefsc_end_year", 
                            label = "End Year for Estimating Summary Metrics:", 
                            min_value = 1963,
                            max_value = 2018, 
                            step_size = 1, 
                            initial_value = 2018),
            
            material_button(
                input_id = "nefsc_button",
                label = "Update Year Range"
            ),
            
            # Comparison Plots
            #plotOutput("comparison_timeline")
            plotOutput("nefsc_plot")
        )
    ),
    
    
    
    
    
    
    ####__ Tab 2 Content - NEFSC Twin Plots  ####
    material_tab_content(
        tab_id = "second_tab",
        material_card(
            # Card Title
            title = "Northeast Groundfish Survey | Single Years Comparisons",
            
            # Selection Controls:
            tags$h6("Choose a Year to Plot:"),
            material_number_box(
                input_id = "nefsc_year", 
                label = "Start Year for Estimating Summary Metrics:", 
                min_value =  1963, 
                max_value = 2018, 
                step_size = 1,
                initial_value = 2000),
            material_button(
                input_id = "nefsc_button_2",
                label = "Change Selected Year"
            ),
            
            # Comparison Plots
            #plotOutput("comparison_timeline")
            plotOutput("nefsc_plot_2")
        )
    ),
    
    
    
    
    
    ####__ Tab 3 Content - MENH  ####
    material_tab_content(
        tab_id = "third_tab",
        material_card(
            # Card Title
            title = "Maine & New Hampshire Inshore Survey | Overall",
            
            # Selection Controls:
            tags$h6("Step 1. - Choose a Minimum Biomass"),
            material_number_box(input_id = "biomass_cutoff",
                           label = "Minimum Fish Size (g) for Inclusion:", 
                           min_value = 10,
                           max_value = 800, 
                           step_size = 10, 
                           initial_value = 400),
            tags$h6("Step 2. - Choose the Reference Period for the Calculations:"),
            material_slider(input_id = "menh_start_year", 
                            label = "Start Year for Estimating Summary Metrics:", 
                            min_value = 2000,
                            max_value = 2019, 
                            step_size = 1, 
                            initial_value = 2000),
            material_slider(input_id = "menh_end_year", 
                            label = "End Year for Estimating Summary Metrics:", 
                            min_value = 2000,
                            max_value = 2019, 
                            step_size = 1, 
                            initial_value = 2019),
            
            material_button(
                input_id = "menh_button",
                label = "Update Year Range"
            ),
            
            # Comparison Plots
            #plotOutput("comparison_timeline")
            plotOutput("menh_plot")
        )
    ),
    
    
    
    
    ####__ Tab 4 Content - MENH Twin Plots  ####
    material_tab_content(
        tab_id = "fourth_tab",
        material_card(
            # Card Title
            title = "Maine & New Hampshire Inshore Survey | Single Years Comparisons",
            
            # Selection Controls:
            tags$h6("Choose a Year to Plot:"),
            material_number_box(
                input_id = "menh_year", 
                label = "Start Year for Estimating Summary Metrics:", 
                min_value =  2000, 
                max_value = 2019, 
                step_size = 1,
                initial_value = 2000),
            material_button(
                input_id = "menh_button_2",
                label = "Change Selected Year"
            ),
            
            # Comparison Plots
            #plotOutput("comparison_timeline")
            plotOutput("menh_plot_2")
        )
    )


)
#################  Server  #######################
server <- function(input, output) {

    
    ####  Reactive Data for NEFSC Overall Plot  ####
    nefsc_window_plot <- eventReactive(eventExpr = input$nefsc_button, {

        # name for group
        group_tag <- str_c(input$nefsc_start_year, " - ", input$nefsc_end_year)
        
        # Filter the Data
        nefsc_filtered <- nefsc_master %>%
            filter(Year >= as.numeric(input$nefsc_start_year),
                   Year <= as.numeric(input$nefsc_end_year),
                   wmin >= as.numeric(input$biomass_cutoff)) %>% 
            mutate(group_var = group_tag)

        # Run Size Spectra Calculation
        nefsc_results <- nefsc_filtered %>% 
                split(.$group_var) %>%
                imap_dfr(group_mle_calc) %>%
                mutate(
                    stdErr = (abs(confMin - b) + abs(confMax - b)) / (2 * 1.96),
                    Year = str_sub(group_var, 1, 4),
                    Year = as.numeric(as.character(Year)),
                    season = str_sub(group_var, 6, -1),
                    season = factor(season, levels = c("SPRING", "SUMMER", "FALL", "WINTER")),
                    C = (b != -1 ) * (b + 1) / ( xmax^(b + 1) - xmin^(b + 1) ) + (b == -1) * 1 / ( log(xmax) - log(xmin)))

        
        # Data inputs
        dbin_group_list <- nefsc_filtered %>% 
            split(.$group_var) %>% 
            map(isd_plot_prep)
        
        
        # global limits for plot
        xlim.global <- c( min(nefsc_filtered$wmin), max(nefsc_filtered$wmax) )
        
        # Vector of years, for naming
        group_names <- names(dbin_group_list)
        names(group_names) <- names(dbin_group_list)
        
        # Size Spectra Coefficients for each group
        ss_results_list  <- nefsc_results %>% split(.$group_var)      
        
        
        
        # Create the list of plots for the single group
        seasonal_isd <- map2(dbin_group_list, ss_results_list, .f = ggplot_isd) %>% setNames(group_names)
        
        
        # Add
        plot_out <- seasonal_isd[[group_tag]][["obs_y"]] + (seasonal_isd[[group_tag]][["log10_y"]] + labs(title = str_c("Biomass Cutoff set to ", input$biomass_cutoff, "g")))
        
        # Return that list
        return(plot_out)
        
        
        
        
        
        
        
        
        })
        
        
   ####  Reactive Data for MENH Overall Plot  ####

    menh_window_plot <- eventReactive(eventExpr = input$menh_button, {
        
        # name for group
        group_tag <- str_c(input$menh_start_year, " - ", input$menh_end_year)
        
        # Filter the Data
        menh_filtered <- menh_master %>%
            filter(Year >= as.numeric(input$menh_start_year),
                   Year <= as.numeric(input$menh_end_year),
                   wmin >= as.numeric(input$biomass_cutoff)) %>% 
            mutate(group_var = group_tag)
        
        # Run Size Spectra Calculation
        menh_results <- menh_filtered %>% 
            split(.$group_var) %>%
            imap_dfr(group_mle_calc) %>%
            mutate(
                stdErr = (abs(confMin - b) + abs(confMax - b)) / (2 * 1.96),
                Year = str_sub(group_var, 1, 4),
                Year = as.numeric(as.character(Year)),
                season = str_sub(group_var, 6, -1),
                season = factor(season, levels = c("SPRING", "SUMMER", "FALL", "WINTER")),
                C = (b != -1 ) * (b + 1) / ( xmax^(b + 1) - xmin^(b + 1) ) + (b == -1) * 1 / ( log(xmax) - log(xmin)))
        
        
        # Data inputs
        dbin_group_list <- menh_filtered %>% 
            split(.$group_var) %>% 
            map(isd_plot_prep)
        
        
        # global limits for plot
        xlim.global <- c( min(menh_filtered$wmin), max(menh_filtered$wmax) )
        
        # Vector of years, for naming
        group_names <- names(dbin_group_list)
        names(group_names) <- names(dbin_group_list)
        
        # Size Spectra Coefficients for each group
        ss_results_list  <- menh_results %>% split(.$group_var)      
        
        
        
        # Create the list of plots for the single group
        seasonal_isd <- map2(dbin_group_list, ss_results_list, .f = ggplot_isd) %>% setNames(group_names)
        
        
        # Add
        plot_out <- seasonal_isd[[group_tag]][["obs_y"]] + (seasonal_isd[[group_tag]][["log10_y"]] + labs(title = str_c("Biomass Cutoff set to ", input$biomass_cutoff, "g")))
        
        # Return that list
        return(plot_out)
        
    })
    
    
    
    
    ####  Reactive Plot of NEFSC Overall  ####
    output$nefsc_plot <- renderPlot({
        
        
        # #--- Show the spinner ---#
        # material_spinner_show(session, "comparison_plot_reactive")
        # 
        
        # Just plot the stecked version
        p <- nefsc_window_plot()
        
        # #--- Hide the spinner ---#
        # material_spinner_hide(session, "comparison_plot_reactive")
        
        p
        
    })
    
    
    ####  Reactive Plot of MENH Overall  ####
    output$menh_plot <- renderPlot({
        
        # #--- Show the spinner ---#
        # material_spinner_show(session, "comparison_plot_reactive")
        
        # Just plot the stecked version
        p <- menh_window_plot()
        
        # #--- Hide the spinner ---#
        # material_spinner_hide(session, "comparison_plot_reactive")
        
        # Plot it
        p
        
    })
    
    
    ####  Static Plot of Single NEFSC Year  ####
    
    nefsc_year <- eventReactive(eventExpr = input$nefsc_button_2, {
        nefsc_year <- as.character(input$nefsc_year)
    })
    
    
    output$nefsc_plot_2 <- renderPlot({
        
        # Pull from list
        nefsc_static_plots[[nefsc_year()]]
        
        
    })
    
    
    ####  Static Plot of Single NEFSC Year  ####
    
    menh_year <- eventReactive(eventExpr = input$menh_button_2, {
        menh_year <- as.character(input$menh_year)
    })
    
    output$menh_plot_2 <- renderPlot({
        
        # Pull from list
        menh_static_plots[[menh_year()]]
        #ggplot()
        
    })
    
    
    
    
    
   
}

# Run the application 
shinyApp(ui = ui, server = server)
