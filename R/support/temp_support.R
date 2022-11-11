
# Temperature data:
library(tidyverse)

# make daily sst tables yearly
make_yearly <- function(daily_sst){
  daily_sst %>% 
    mutate(
      time = as.Date(time),
      yr = lubridate::year(time)) %>% 
    filter(between(yr, 1982, 2021)) %>% 
    group_by(yr) %>% 
    summarise(across(c(sst, sst_clim, sst_anom), mean, na.rm = T)) %>% 
    ungroup() %>% 
    mutate(decade = floor_decade(yr))

}
  

