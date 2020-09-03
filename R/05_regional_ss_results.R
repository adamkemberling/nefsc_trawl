####

# Goal: attribute ocean region to nefsc trawl size spectra
# Do areas of inherently different temperature reflect differences in slope/shape
# Is this holding up, does GOM resemble mid atlantic now
# Want to end up with a chronological record of what the slope is for every year/season/region and how that region did for the entire year as well



####  Packages  ####
library(here)
library(gmRi)
library(sf)
library(janitor)
library(sizeSpectra)
library(patchwork)
library(tidyverse)

####  Support Functions  ####
source(here("R/support/sizeSpectra_support.R"))



#File Paths
mills_path <- shared.path(os.use = "unix", group = "Mills Lab", folder = NULL)
res_path <- shared.path(os.use = "unix", group = "RES Data", folder = NULL)

####  Data  ####




#####__ a. NEFSC bottom trawl  ####
nefsc_master <- read_csv(here("data/NEFSC/2020survdat_nye.csv"), 
                        col_types = cols(), 
                        guess_max = 1e6) %>% 
  clean_names() 



# Cleanup Code from Kathy
nefsc_master <- nefsc_master%>% 
  mutate(comname = tolower(comname),
         id = format(id, scientific = FALSE)) %>%
  # select(id, est_year, season, stratum, decdeg_beglat, decdeg_beglon,
  #        svspp, comname, catchsex, biomass, avgdepth, abundance, length, numlen) %>%
  filter(stratum >= 01010,
         stratum <= 01760,
         stratum != 1310,
         stratum != 1320,
         stratum != 1330,
         stratum != 1350,
         stratum != 1410,
         stratum != 1420,
         stratum != 1490,
         season %in% c("SPRING", "FALL"),
         est_year >= 1970,
         est_year < 2020) %>%
  mutate(biomass = ifelse(is.na(biomass) == TRUE & abundance > 0, 0.01, biomass),
         biomass = ifelse(biomass == 0 & abundance > 0, 0.01, biomass),
         abundance = ifelse(is.na(abundance) == TRUE & biomass > 0, 1, abundance),
         abundance = ifelse(abundance == 0 & biomass > 0, 1, abundance),
         tstrat = str_sub(stratum, 2, 3)) %>%  
  filter(!is.na(biomass),
         !is.na(abundance))


####__ b. L/W coefficients  ####

# Load the lw key with coefficients from fishbase
nefsc_lw <- read_csv(here("data/NEFSC/nefsc_lw_key_filled.csv"),
                     guess_max = 1e3,
                     col_types = cols())

# Fill in NA growth coefficients using related species
nefsc_lw <- nefsc_lw %>% 
  mutate(a = as.numeric(a),
         b = as.numeric(b),
         ln_a = log(a),
         related_ln_a = log(related_a),
         ln_a = ifelse(is.na(ln_a), related_ln_a, ln_a))



# Wrigley Paper, length weight coefficients
load("/Users/akemberling/Box/RES Data/NMFS_trawl/lwreg.Rdata")
lwreg <- lwreg %>% clean_names()  
lwreg <- lwreg %>%  
  mutate(comname = str_to_lower(common_name), .before = scientific_name,
         common_name = NULL,
         scientific_name = str_to_lower(scientific_name),
         svspp = str_pad(svspp, width = 3, side = "left", pad = "0"))



#CHECKING PORTION OF FISH REPRESENTED IN WIGLEY ET AL 2003 BASED ON TOW-LEVEL DATA
spp_classes <- read_csv(here("data/kmills/sppclass.csv"),
                        col_types = cols())
spp_classes <- spp_classes %>% 
  clean_names() %>% 
  mutate(common_name = str_to_lower(common_name),
         scientific_name = str_to_lower(scientific_name)) %>% 
  pivot_longer(gf:inv, names_to = "spec_class", values_to = "flag") %>% 
  pivot_longer(com:nc, names_to = "fishery", values_to = "flag2") %>% 
  filter(is.na(flag) == FALSE,
         is.na(flag2) == FALSE) %>% 
  select(svspp, comname = common_name, everything(), -flag, -flag2)
 

# Merge the growth coefficients with the species classes
bn <- left_join(lwreg, spp_classes, by = c("svspp", "comname", "scientific_name"))



####  Next Steps:   ####

# 1. use these growth coefficients to get biomasses
# 2. Potentially merge in my fishbase ones so that there is just one master
# 3. From there we can do the rest of the size spectra steps





#CHECKING PORTION OF FISH REPRESENTED IN WIGLEY ET AL 2003 BASED ON EXPNUMLEN

sum(bn$NUMLEN,na.rm=TRUE)

bn.wigley<-bn[bn$wigley=="X",]

sum(bn.wigley$NUMLEN,na.rm=TRUE)

sum(bn.wigley$NUMLEN,na.rm=TRUE)/sum(bn$NUMLEN,na.rm=TRUE)

#82%



head(bn)

fish<-bn[bn$INV==0,]

sum(fish$NUMLEN,na.rm=TRUE)

fish.wigley<-fish[fish$wigley=="X",]

sum(fish.wigley$NUMLEN,na.rm=TRUE)

sum(fish.wigley$NUMLEN,na.rm=TRUE)/sum(fish$NUMLEN,na.rm=TRUE)

#94.7% of fish have L-W coefficients in Wigley



########################



#GETTING WEIGHTS FROM L-W COEFFICIENTS



fish.cut<-fish.wigley[fish.wigley$LENGTH!=0,]

fish.cut<-fish.cut[!(is.na(fish.cut$LENGTH)),]

#fish.cut[fish.cut$LENGTH==0,]



fish.cut$LLEN<-log(fish.cut$LENGTH)

class<-merge(fish.cut,sppclass2,by.x="SVSPP",by.y="SVSPP")

class.wigley<-class[class$wigley=="X",]

len.reg<-merge(class.wigley,lwreg,all.x=T)



len.reg$LWEIGHT<-len.reg$lna + (len.reg$b*len.reg$LLEN)

len.reg$WT<-exp(len.reg$LWEIGHT)

len.reg$TOTWT<-len.reg$WT * len.reg$NUMLEN




####_____________________________________________________####
####  Data  Prep  ####


#####__ 1.  Stratum-Area Key  ####


# If we know which strata go where we can use this:
# Stratum Key for filtering specific areas
strata_key <- list(
  "Georges Bank"          = as.character(13:23),
  "Gulf of Maine"         = as.character(24:40),
  "Southern New England"  = as.character(1:10),
  "Mid-Atlantic Bight"    = as.character(61:76))

# Pull the strata we want
strata_select <- c(
  strata_key$`Georges Bank`, 
  strata_key$`Gulf of Maine`,
  strata_key$`Southern New England`,
  strata_key$`Mid-Atlantic Bight`)





####__ 2. Temporal Filtering  ####
#Filtering
nefsc <- nefsc_master %>% 
  mutate(area =  case_when(
           tstrat %in% strata_key$`Georges Bank`         ~ "GB", 
           tstrat %in% strata_key$`Gulf of Maine`        ~ "GoM",
           tstrat %in% strata_key$`Southern New England` ~ "SNE",
           tstrat %in% strata_key$`Mid-Atlantic Bight`   ~ "MAB",
           TRUE                                          ~ "not found")) %>% 
  filter(est_year %in% c(1982:2019),
         season   %in% c("SPRING", "FALL"),
         tstrat   %in% strata_select)


# Check what is going on with numlen and abundance, make conversion factor
conv_factor <- nefsc %>% 
  group_by(id, comname, catchsex, numlen, abundance) %>% 
  summarise(abundance_raw = sum(numlen),
            abundance = mean(abundance), 
            convers = abundance/abundance_raw) 


# Merge back and convert the numlen field
nefsc <- nefsc %>% 
  left_join(conv_factor, by = c("id", "comname", "catchsex", "numlen", "abundance")) %>% 
  mutate(numlen_adj = numlen * convers,
         numlen_adj = round(numlen_adj, 2)) 

# Formatting individual Lengths

nefsc <- nefsc %>% 
  select(c(id,
           est_year:svvessel, 
           decdeg_beglat, 
           decdeg_beglon, 
           catchsex,
           comname:area,
           catchsex,station:numlen_adj))


glimpse(nefsc)


####__ 3. Get Individual Lengths  ####

# Record of unique station catches# row for each species*sex*length
nefsc_lens <- nefsc %>% 
  filter(is.na(numlen) == FALSE,
         numlen > 0) %>% 
  distinct(id, comname, catchsex, length, numlen_adj) 


# recombine with the distinct station info
nefsc_spectra <- nefsc %>% 
  # drop the columns with individual catch info, we want all distinct records merging in ok
  select(-c(est_julian_day, catchsex, comname, numlen_adj, biomass, abundance, length, numlen, abundance_raw, convers)) %>% 
  distinct() %>% left_join(nefsc_lens, by = "id")


glimpse(nefsc_spectra)


# Combine with the length weight coefficients
####__ 4.  Calculate Weights  ####
nefsc_weights <- nefsc_spectra %>% 
  left_join(nefsc_lw[,1:5], by = "comname")  %>% 
  mutate(b = as.numeric(b),
         a = as.numeric(a),
         ln_a = log(a),
         ln_weight = (ln_a + b * log(length)),
         weight_g = exp(ln_weight),
         freq_weight = weight_g * numlen_adj) %>% 
  drop_na(weight_g)


# Check weights
nefsc_weights %>% 
  group_by(comname) %>% 
  summarise(mean_weight =  mean(weight_g)) %>% 
  arrange(desc(mean_weight)) 


 nefsc_weights %>% 
  ggplot( aes( y = fct_reorder(comname, weight_g, .fun = mean, .desc = TRUE), x = weight_g)) + 
  geom_boxplot() + 
  labs(x = "Mean Weight (g)", y = "Common Name")






####__ 5. Estimate Biomass  ####

 
 
 
 ####_____________________________________________________####
 
 ####  Size Spectra Prep  ####
 
 
 # Starting point for sizeSpectra steps
 dataOrig <- nefsc_weights %>% filter(is.na(weight_g) == FALSE)
 
 
 # Keep desired columns, name them for the vignette
 data <- dataOrig %>%  
   mutate(season = factor(season, levels = c("SPRING", "SUMMER", "FALL", "WINTER"))) %>%
   select(
     Year = est_year,                     # year
     samp_month = est_month,
     samp_dat = est_day,
     samp_time = est_time,
     vessel = svvessel,
     season,
     SpecCode = comname,          # common name
     hare_group,
     LngtClass = length,              # length bin
     Number = numlen_adj,              # CPUE in n/effort
     LWa = a,                         # length/weight param
     LWb = b,                         # length/weight param
     bodyMass = weight_g,             # weight of an individual from that size class
     CPUE_bio_per_hour = freq_weight  # CPUE in biomass
   ) %>% arrange(Year, season, SpecCode, LngtClass) %>% 
   mutate(Biomass = Number * bodyMass)
 
 ####__  SS Species Size Bins  ####
 
 # Split for every species
 species_splits <- data %>% 
   split(.$SpecCode)
 
 
 # Get a key for each species and size class
 # assumes 1cm bins for all species as written
 data_bin_key <- species_splits %>% 
   map_dfr(function(species_df){
     
     #pull the distinct length bins
     species_df <- species_df %>% 
       distinct(LngtClass, .keep_all = T) %>% 
       arrange(LngtClass)
     
     
     # Add the max length for the bin, and its weight
     binned_df <- species_df %>% 
       mutate(
         LngtMax = LngtClass + 1, 
         wmax    = exp(log(LWa) + LWb * log(LngtMax)))  %>%
       select(SpecCode, LWa, LWb, LngtMin = LngtClass, wmin = bodyMass, LngtMax, wmax,
              -c(Number, Biomass))
     
     # return the clean data
     return(binned_df)})
 
 
 # Add the bins back into the original and clean up
 dataBin <- data %>% 
   select(Year, season, SpecCode, LngtMin = LngtClass, Number, Biomass) %>% 
   left_join(data_bin_key, by = c("SpecCode", "LngtMin")) 

 
 ####__  Set Bodymass Cutoff and Groups  ####
 
 # Set bodymass lower limit
 # Filter for lower end of gear selectivity
 mass_cutoff <- 400 #grams
 dataBin <- filter(dataBin, wmin >= mass_cutoff)
 

 
 
 ####_____________________________________________________#### 
 ####  Size Spectra Group Calculations  ####
 
 #####__ 1.  All years, every region  ####
 g1 <- data %>% 
   mutate(group_level = "all_data") %>% 
   group_by(group_level)
   summarise(
     #Divide Number / grouping_variable (numAreas),
     Number = sum(Number), 
     LWa = unique(LWa),
     LWb = unique(LWb),
     bodyMass = unique(bodyMass)) %>% 
   ungroup()
 
 #####__ 2. All Years, each season  ####
 
 
 #####__ 3. All Years, region * Season  ####
 
 
 #####__ 4. Every year, entire survey  ####
 
 
 #####__ 5. every year, every region  ####
 
 
 #####__ 6. every year, only seasons  ####
 
 
 #####__ 7. every year, region * season  ####




