# Date 3/21/2023
# Subject: Northeast US SST, decadal variability paper:
# Premise: Provide figures on SST trends for the area for warming context

# Packages
library(lubridate)
library(scales)
library(gmRi)
library(targets)
library(tidyverse)
library(patchwork)
library(heatwaveR)


# Set theme
theme_set(
  theme_classic() + 
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 16)  
    )
)

# Degree symbol
deg_c <- "\u00b0C"


####  Load Data  ####


# Load SST
shelf_sst <- oisst_access_timeseries("trawl", "inuse strata", box_location = "cloudstorage")

# Annual averages
annual_sst  <- shelf_sst %>% 
  group_by(year = lubridate::year(time)) %>% 
  summarise(across(-c(time, modified_ordinal_day), mean)) %>% 
  filter(year %in% c(1982:2021)) %>% 
  mutate(period = ifelse(year < 2010, "early", "late"))


# Averages for comparison periods
early_avg <- annual_sst %>% filter(year %in% c(1982:2009)) %>% summarise(across(-c(year,period), mean))
late_avg <- annual_sst %>% filter(year %in% c(2010:2021)) %>% summarise(across(-c(year,period), mean))



# Load global temps, get avgs
# # Global Temperature Anomalies
oisst_path <- cs_path("res", "OISST/oisst_mainstays")
global_anoms <- read_csv(
  paste0(oisst_path, "global_timeseries/global_anoms_1982to2011.csv"), 
  guess_max = 1e6,
  col_types = cols()) %>% 
  mutate(year = lubridate::year(time)) %>% 
  filter(between(year, 1982, 2021))

# summarize by year again
global_annual <- global_anoms %>% 
  group_by(year) %>% 
  summarise(across(-c(time, MOD), mean)) %>% 
  filter(year %in% c(1982:2021)) 




# Figure 2. 

"
Sea surface temperature on the Northeast Shelf. (A) Annual average SST anomalies
(dots and thin lines) and means (solid thick lines) for the time periods 1982-2009 (blue) and
2010-2021 (orange). The 1982-2021 warming trend for the region (0.38 °C per decade, solid
black line) is compared to the global SST trend (0.15 °C/decade, dashed black line).
"





# Build plot showcasing the Decadal Temperature Jump
shelf_warming_plot <- ggplot(annual_sst, aes(year, area_wtd_anom)) +
  # geom_ribbon(
  #   color = "#F4A582", linetype = 2,
  #   aes(
  #     ymin = -Inf, 
  #     ymax = predict(lm(area_wtd_anom ~ year))),
  #   alpha = 0.8, fill = "#FDDBC7") +
  # geom_ribbon(
  #   data = global_annual,
  #   color = "#92C5DE", linetype = 2,
  #   aes(
  #     ymin = -Inf, 
  #     ymax = predict(lm(area_wtd_anom ~ year))),
  #   alpha = 0.8, fill = "#D1E5F0") +
  geom_line(linewidth = 0.5, linetype = 3) +
  geom_point(aes(color = period), size = 2.5) +
  geom_segment(data = early_avg, aes(x = 1982, xend = 2009, y = area_wtd_anom, yend = area_wtd_anom, color = "early"), linewidth = 1.5) +
  geom_segment(data = late_avg, aes(x = 2010, xend = 2019, y = area_wtd_anom, yend = area_wtd_anom, color = "late"), linewidth = 1.5) +
  # annotate(x = 2015, y = -0.5, geom = "text", label = "Global Oceans\nBaseline Warming", fontface = "bold") +
  # annotate(x = 2015, y = 0.5, geom = "text", label = "Northeast US\nShelf Warming", fontface = "bold") +
  guides(color = "none") +
  scale_y_continuous(labels = number_format(suffix = deg_c)) +
  scale_color_gmri() +
  labs(y = "SST Anomaly", x = "Year")




# Shelf Warming Plot
shelf_warming_plot


####  Figure 2b.  Heatmap  ####

"
(B) A
heatmap of the daily temperature anomalies from Jan. 1, 1982 to Dec. 31, 2021, with days that
qualify as marine heatwaves indicated with the black horizontal lines.

"


##### Function to Identify MHW  ####
pull_heatwave_events <- function(temperature_timeseries, 
                                 clim_ref_period = c("1982-01-01", "2011-12-31"),
                                 date_col = "time",
                                 temp_col = "sst",
                                 threshold = 90,
                                 detrend = FALSE) {
  
  # temperature_timeseries <- gom_sst
  
  
  # Pull the two column dataframe for mhw estimation
  test_ts <- data.frame(t = as.Date(temperature_timeseries[[date_col]]), 
                        temp = temperature_timeseries[[temp_col]])
  
  
  # Calculate seasonally varying climatology with threshold w/ smoothing window
  ts  <- heatwaveR::ts2clm(data = test_ts, 
                           climatologyPeriod = clim_ref_period, 
                           pctile = threshold) %>% 
    mutate(sst_anom = temp - seas,
           yr = lubridate::year(t))
  
  
  
  # Perform linear detrending on anomalies
  if(detrend){
    
    # Detrend day of year temperature trends:
    ts <- ts %>% 
      split(.$doy) %>% 
      map_dfr(detrend_sst, vals = "sst_anom", yr_col = "yr") %>% 
      mutate(detrend_temp = seas + detrend_vals) %>% 
      arrange(t)
    
  }
  
  
  # Perform Heatwave Detection
  mhw <- ifelse(detrend,
                heatwaveR::detect_event(ts, x = t, y = detrend_temp),
                heatwaveR::detect_event(ts, x = t, y = temp))
  
  
  
  # Select and rename critical heatwave data
  mhw_out <- mhw[[1]] %>% 
    #mutate(sst_anom = temp - seas) %>% 
    rename(time = t,
           sst = temp,
           mhw_thresh = thresh,
           mhw_threshCriterion = threshCriterion,
           mhw_durationCriterion = durationCriterion,
           mhw_event = event,
           mhw_event_no = event_no)
  
  
  # Repeat for cold spells
  # 2. Detect cold spells
  # coldSpells = TRUE flips boolean to < thresh
  ts <- ts2clm(data = test_ts, 
               climatologyPeriod = clim_ref_period, 
               pctile = (100 - threshold)) %>% 
    mutate(sst_anom = temp - seas,
           yr = lubridate::year(t))
  
  
  # Perform linear detrending on anomalies
  if(detrend){
    
    # Detrend day of year temperature trends:
    ts <- ts %>%
      split(.$doy) %>%
      map_dfr(detrend_sst, vals = "sst_anom", yr_col = "yr") %>%
      mutate(detrend_temp = seas + detrend_vals) %>% 
      arrange(t)
    
  }
  
  
  
  # Perform Cold Spell Detection
  mcs <- ifelse(detrend,
                heatwaveR::detect_event(ts, x = t, y = detrend_temp, coldSpells = T),
                heatwaveR::detect_event(ts, x = t, y = temp, coldSpells = T))
  
  
  
  # Prepare cold spell data to join
  # Remove columns that are shared with heatwaves
  mcs_out <- mcs[[1]] %>%
    dplyr::select(time = t,
                  mcs_thresh = thresh,
                  mcs_threshCriterion = threshCriterion,
                  mcs_durationCriterion = durationCriterion,
                  mcs_event = event,
                  mcs_event_no = event_no)
  
  
  # join heatwave detection results to coldspell results
  hot_and_cold <- left_join(mhw_out, mcs_out, by = "time")
  
  
  # 3. Data formatting for plotting, 
  # adds columns to plot hw and cs seperately
  events_out <- hot_and_cold %>% 
    mutate(
      # Set up status to combine labelling for heatwaves and cold spells:
      status   = ifelse(mhw_event == TRUE, "Marine Heatwave Event", "Sea Surface Temperature"),
      status   = ifelse(mcs_event == TRUE, "Marine Cold Spell Event", status),
      event_type = ifelse(detrend, "Jacox Method", "Hobday Method"),
      # Corrective measures for where event flagging is off:
      # status   = ifelse(sst > mhw_thresh, "Marine Heatwave Event", status),
      # status   = ifelse(sst < mcs_thresh, "Marine Cold Spell Event", status),
      # Heatwave event temperature values:
      hwe      = ifelse(mhw_event == TRUE, sst, NA),
      cse      = ifelse(mcs_event == TRUE, sst, NA),
      nonevent = ifelse(mhw_event == FALSE & mcs_event == FALSE, sst, NA)) 
  
  # Close the gaps between a mhw event and sst (might not need if full line for temp exists)
  events_out <- events_out %>% 
    mutate(hwe = ifelse( (is.na(hwe) & is.na(lag(hwe, n = 1))) == FALSE, sst, hwe),
           cse = ifelse( (is.na(cse) & is.na(lag(cse, n = 1))) == FALSE, sst, cse)) %>% 
    distinct(time, .keep_all = T)
  
  
  return(events_out)
}



####  Prepare MHW Data  ####

base_date <- as.Date("2000-01-01")
shelf_hw <- pull_heatwave_events(shelf_sst, clim_ref_period = c("1982-01-01", "2011-12-31")) %>% 
  mutate(yday = yday(time),
         flat_date = as.Date(yday-1, origin = base_date),
         year = year(time))




####  Heatmap  ####


# Color limit for palettes
temp_limits <- c(-3, 3)
temp_suff <- deg_c
temp_breaks <- c(temp_limits[[1]], temp_limits[1]/2,  0, temp_limits[2]/2, temp_limits[2])
temp_labels <- str_c(c(str_c("< ", temp_limits[1]), temp_limits[1]/2, 0, temp_limits[2]/2, str_c("> ", temp_limits[2])), temp_suff)




# Assemble heatmap plot
hw_dat <- shelf_hw

# Make plot
shelf_heatmap <- ggplot(hw_dat, aes(x = flat_date, y = year)) +
  
  # background box fill for missing dates
  geom_rect(xmin = base_date, 
            xmax = base_date + 365, 
            ymin = min(hw_dat$year) - .5, 
            ymax = max(hw_dat$year) + .5, 
            fill = "gray75", color = "transparent") +
  
  # tile for sst colors
  geom_tile(aes(fill = sst_anom)) +
  # points for heatwave events
  geom_point(data = filter(hw_dat, mhw_event == TRUE),
             aes(x = flat_date, y = year), size = .25)  +
  # # Points if over threshold only
  scale_x_date(date_labels = "%b", date_breaks = "1 month", expand = expansion(add = c(0,-1))) +
  scale_y_continuous(limits = c(1982 + .5, 2021 + .5), expand = expansion(add = c(0,0))) +
  scale_fill_distiller(palette = "RdBu", 
                       na.value = "transparent", 
                       limit = temp_limits, 
                       oob = scales::squish,
                       breaks = temp_breaks, 
                       labels = temp_labels) +
  
  #5 inches is default rmarkdown height for barheight
  theme(legend.title = element_text(angle = 90)) +
  guides(fill = guide_colorbar(
    title.position = "right", title.hjust = 0.5,
    barheight = unit(10, "cm"), 
    frame.colour = "black", 
    ticks.colour = "black")) +
  labs(x = "Month", 
       y = "Year",
       fill = "Sea Surface Temperature Anomaly")





# Show it
shelf_heatmap





####  Warming Map  ####
library(sf)
library(terra)
library(tidyterra)
library(rcartocolor)

# Shapefile
countries <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
                                         
# Warming rates
list.files(str_c(oisst_path, "warming_rates/"))
sst_rates <- rast(str_c(oisst_path, "warming_rates/annual_warming_rates1982to2021.nc"))

# Set the extent
ext(sst_rates) <- ext(0, 360, -90, 90)

# Flip and rotata
sst_rates <- terra::flip(sst_rates)
sst_rates <- rotate(sst_rates)


# Rotate
#sst_rates %>% plot()

# Build the map
rates_map <- ggplot() +
  geom_spatraster(data=sst_rates$annual_warming_rate) +
  geom_sf(data = countries) +
  scale_y_continuous(breaks = c(36,40,44)) + 
  scale_x_continuous(breaks = c(-78,-72,-66)) +
  #scale_fill_carto_c(palette = "Geyser", na.value = "transparent") +
  #scale_fill_distiller(palette = "Spectral", na.value = "transparent", direction = -1, limits = c(0, 0.1)) +
  scale_fill_distiller(palette = "RdYlBu", na.value = "transparent", direction = -1, limits = c(0, 0.1)) +
  theme_bw() +
  map_theme(legend.position = "bottom", legend.title = element_text(angle = 0)) +
  guides(fill = guide_colorbar(
    ticks = T, ticks.colour = "black", frame.colour = "black", barwidth = unit(8, "cm"), 
    title.position = "top", title.hjust = 0.5))+
  #guides(fill = guide_colorsteps(show.limits = F, ticks = T, tick.color = "black", frame.color = "black"))+
  coord_sf(
    xlim=c(-76, -66), 
    ylim=c(37.5,47)) +
  labs(x = NULL, y = NULL, 
       fill = expression("Annual Warming Rate \u00b0C Year"^-1))
  

rates_map


####  Saving Figures  #####

# Path to decadal folder:
decadal_folder <- cs_path(
  box_group = "mills", 
  subfolder = "Projects/Decadal Variability/Revisions/SST/")


ggsave(filename = str_c(decadal_folder, "ne_shelf_warming_rate.png"), 
       plot = shelf_warming_plot, height = 4, width = 6, units = "in", dpi = "retina")
ggsave(filename = str_c(decadal_folder, "ne_shelf_heatwaves.png"),
       plot = shelf_heatmap, height = 8, width = 10, units = "in", dpi = "retina")
ggsave(filename = str_c(decadal_folder, "shelf_warming_rate_map.png"), 
       plot = rates_map, height = 6, width = 5, units = "in", dpi = "retina")



# Make a mosaic

heatmap_smol <- ggplot(hw_dat, aes(x = flat_date, y = year)) +
  
  # background box fill for missing dates
  geom_rect(xmin = base_date, 
            xmax = base_date + 365, 
            ymin = min(hw_dat$year) - .5, 
            ymax = max(hw_dat$year) + .5, 
            fill = "gray75", color = "transparent") +
  
  # tile for sst colors
  geom_tile(aes(fill = sst_anom)) +
  # points for heatwave events
  geom_point(data = filter(hw_dat, mhw_event == TRUE),
             aes(x = flat_date, y = year), size = .25)  +
  # # Points if over threshold only
  scale_x_date(date_labels = "%b", date_breaks = "1 month", expand = expansion(add = c(0,-1))) +
  scale_y_continuous(limits = c(1982 + .5, 2021 + .5), expand = expansion(add = c(0,0))) +
  scale_fill_distiller(palette = "RdBu", 
                       na.value = "transparent", 
                       limit = temp_limits, 
                       oob = scales::squish,
                       breaks = temp_breaks, 
                       labels = temp_labels) +
  
  #5 inches is default rmarkdown height for barheight
  theme(legend.title = element_text(angle = 0),
        legend.position = "bottom") +
  guides("fill" = guide_colorbar(
    title = "Sea Surface Temperature Anomaly", 
    title.position = "top", 
    title.hjust = 0.5,
    barwidth = unit(8, "cm"), 
    frame.colour = "black", 
    ticks.colour = "black")) +  
  labs(x = "Month", 
       y = "Year",
       fill = "Sea Surface Temperature Anomaly")




sst_mosaic <- (shelf_warming_plot/ heatmap_smol) | rates_map
sst_mosaic

# # Design tweaking
# design <- "
#   113
#   223
#   223
#   22#
# "
# shelf_warming_plot + heatmap_smol + rates_map + plot_layout(design = design)

# Save that jabroni
ggsave(filename = str_c(decadal_folder, "shelf_warming_mosaic.png"), 
       plot = sst_mosaic, height = 10, width = 10, units = "in", dpi = "retina")





####
# Climate Velocity Map
library(climetrics)

# Annual Averages
sst <- rast(str_c(oisst_path, "yearly_averages/annual_sst_1982to2021.nc"))

# Set the extent
ext(sst) <- ext(0, 360, -90, 90)

# Flip and rotate
sst <- terra::flip(sst)
sst <- rotate(sst)




