
# Temperature data:
library(tidyverse)

# make daily sst tables yearly
make_yearly <- function(daily_sst){
  daily_sst %>% mutate(yr = lubridate::year(time)) %>% 
    group_by(yr) %>% 
    summarise(across(c(sst, sst, sst_clim, sst_anom), mean)) %>% 
    ungroup() %>% 
    mutate(decade = floor_decade(yr))

}
  

