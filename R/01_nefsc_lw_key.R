####  NEFSC Trawl Data Access  ####

# COMMENTS:
# This script contains the processing code for the length weight coefficient key
# it combines the Wigley06 coefficients with the values found on fishbase

####  Packages  ####
library(janitor)
library(here)
library(gmRi)
library(tidyverse)
library(patchwork)
library(rfishbase)


####  Data  ####

####  Paths  
box_paths  <- research_access_paths(os.use = "unix")
mills_path <- box_paths$mills
nsf_path   <- box_paths$okn
res_path   <- box_paths$res

# 2019 data From Janet Nye
#load(str_c(mills_path, "Data/Survdat_Nye_allseason.RData"))

# Gonna save it locally for convenience/laziness
#write_csv(survdat, here("data/NEFSC/2019survdat_nye.csv"))

# 2020 edition
load(paste0(res_path, "NMFS_trawl/Survdat_Nye_Aug 2020.RData"))

# Gonna save it locally for convenience/laziness
# write_csv(survdat, here("data/NEFSC/2020survdat_nye.csv"))



####____________________####


####  Combined L-W Key  ####
# Sources:
# Wigley Paper
# Fishbase




# Fishbase Growth Coefficients - Manually found
fishbase_lw <- read_csv(here("data/NEFSC/nefsc_lw_key_filled.csv"),
                     guess_max = 1e3,
                     col_types = cols())


# wigley Paper, Load length weight coefficients
load(paste0(res_path, "/NMFS_trawl/lwreg.Rdata"))
wigley_lw <- lwreg; rm(lwreg)

# Save as csv to double check coefficients against paper
write_csv(wigley_lw, here("data/kmills/wigley_06_lwreg.csv"))


# Species class/groups based on life history
spp_classes <- read_csv(here("data/kmills/sppclass.csv"),
                        col_types = cols())



####  Growth Coefficients Prep  ####

####__ 1.  Fishbase Coefficients  ####
# Fill in NA growth coefficients using related species (typically same genera)
fishbase_lw <- fishbase_lw %>% 
  mutate(a = as.numeric(a),
         b = as.numeric(b),
         ln_a = log(a),
         related_ln_a = log(related_a),
         ln_a = ifelse(is.na(ln_a), related_ln_a, ln_a),
         a = ifelse(is.na(a), exp(ln_a), a))

# # validation plots
# (nefsc_a <- ggplot(fishbase_lw, aes(x = a)) + geom_histogram() + labs(caption = "source: fishbase"))
# (nefsc_lna <- ggplot(fishbase_lw, aes(x = ln_a)) + geom_histogram() + labs(caption = "source: fishbase"))
# (nefsc_b <- ggplot(fishbase_lw, aes(x = b)) + geom_histogram() + labs(caption = "source: fishbase"))
# 
# # Check that ln_a is consistent with what it would be if we calculated them right now from a
# fishbase_lw %>% 
#   mutate(lna_check = log(a)) %>% 
#   ggplot(aes(ln_a, lna_check)) +
#   geom_point()







####__ 2. Susan Wigley Paper Coefficients  ####
wigley_lw <- wigley_lw %>% clean_names()  
wigley_lw <- wigley_lw %>%  
  mutate(comname = str_to_lower(common_name), .before = scientific_name,
         common_name = NULL,
         scientific_name = str_to_lower(scientific_name),
         svspp = str_pad(svspp, width = 3, side = "left", pad = "0")) %>% 
  mutate(ln_a = lna, 
         a = exp(ln_a),
         lna = NULL) %>% 
  select(svspp, comname, scientific_name, season, catchsex, a, ln_a, everything())




# # Validation Plots
# (wigley_a   <- ggplot(wigley_lw, aes(x = a)) + geom_histogram() + labs(caption = "source: wigley"))
# (wigley_lna <- ggplot(wigley_lw, aes(x = ln_a)) + geom_histogram() + labs(caption = "source: wigley"))
# (wigley_b   <- ggplot(wigley_lw, aes(x = b)) + geom_histogram() + labs(caption = "source: wigley"))
# wigley_lw %>% 
#   mutate(lna_check = log(a)) %>% 
#   ggplot(aes(ln_a, lna_check)) +
#   geom_point()







####__ 3. Classes & Economic Groups  ####

# Clean up Species Life History & Fisheries Relevance
spp_classes <- spp_classes %>% 
  clean_names() %>% 
  mutate(common_name = str_to_lower(common_name),
         scientific_name = str_to_lower(scientific_name)) %>% 
  pivot_longer(gf:inv, names_to = "spec_class", values_to = "flag") %>% 
  pivot_longer(com:nc, names_to = "fishery", values_to = "flag2") %>% 
  filter(is.na(flag) == FALSE,
         is.na(flag2) == FALSE) %>% 
  select(svspp, comname = common_name, everything(), -flag, -flag2)




# Merge the growth coefficient tables with the species classes  

# Fishbase
fishbase_lw <- left_join(fishbase_lw, spp_classes, by = "comname")

# wigley
wigley_lw <- left_join(wigley_lw, spp_classes, by = c("svspp", "comname", "scientific_name"))




####____  test: compare value ranges   ####

nefsc_a | wigley_a
nefsc_lna | wigley_lna
nefsc_b | wigley_b

# weird outlier for a
wigley_lw[which(wigley_lw$a == max(wigley_lw$a)) ,]

# same for ln_a?
wigley_lw[which(wigley_lw$ln_a == max(wigley_lw$ln_a)) ,]






####__ 4.  Combining Growth Coefficient Sources  ####


# 1. use these growth coefficients to get biomasses
# 2. Potentially merge in my fishbase ones so that there is just one master

wigley_lw <- wigley_lw %>% mutate(source = "wigley") %>% select(source, everything())
fishbase_lw  <- fishbase_lw %>% mutate(source = "fishbase") %>% select(source, everything())
lw_combined <- full_join(fishbase_lw, wigley_lw, by = c("source", "comname", "b", "a", "ln_a", "svspp", "scientific_name", "spec_class", "fishery")) 
#lw_combined <- full_join(fishbase_lw, wigley_lw, by = c("source", "comname", "b", "a", "ln_a")) 






#### Gap Fill for NA's  ####
# Fill in gaps for things that should be consistent
# Lot of matching by common name and pulling the first value where it matches to fill the NA values

# Create named vectors to use as lookup keys
hare_lookup    <- setNames(fishbase_lw$hare_group, fishbase_lw$comname )
svspp_lookup   <- setNames(wigley_lw$svspp, wigley_lw$comname )
sci_lookup     <- setNames(wigley_lw$scientific_name, wigley_lw$comname)
fishery_lookup <- setNames(wigley_lw$fishery, wigley_lw$comname)
class_lookup   <- setNames(wigley_lw$spec_class, wigley_lw$comname)

# This is/was broken when not split first
lw_combined <- lw_combined %>% 
  arrange(comname) %>%
  split(.$comname) %>% 
  map_dfr(function(x){
    
    # Grab the common name, first one in case of repeats
    species_selector <- x$comname[1]
  
    # Now fill in the NA values by using the species as a key
    x <- x %>% 
      mutate(
        # if hare group missing, match the common names, pull hare group, else leave blank
        hare_group = ifelse(is.na(hare_group) & species_selector %in% names(hare_lookup),
                            hare_lookup[[species_selector]][1],
                            hare_group),
    
        # if svspp missing, match the common names, pull svspp value, else leave blank
        svspp = ifelse(is.na(svspp)& species_selector %in% names(svspp_lookup),
                       svspp_lookup[[species_selector]],
                       svspp),
    
        # if scientific name missing, match the common names, pull sci name, else leave blank
        scientific_name = ifelse(is.na(scientific_name) & species_selector %in% names(sci_lookup),
                                 sci_lookup[[species_selector]],
                                 scientific_name),
    
        # if fishery missing, match the common names, pull fishery, else leave blank
        fishery = ifelse(is.na(fishery) & species_selector %in% names(fishery_lookup),
                         fishery_lookup[[species_selector]],
                         fishery),
    
        # if fishery missing, match the common names, pull fishery, else leave blank
        spec_class = ifelse(is.na(spec_class) & species_selector %in% names(class_lookup),
                            class_lookup[[species_selector]],
                            spec_class))

  })
  
  
  
  
# Grab and order columns we need
lw_combined <- lw_combined %>% 
  select(source, svspp, comname, scientific_name, spec_class, hare_group, fishery, season, catchsex,
         b, a, ln_a, units, related_species, related_a, related_b, related_ln_a,
         wigleymin, wigleymax, wigleyvalue)


####____  test: compare value ranges   ####

# Double check that the values are in the same ballpark:
(combined_a   <- ggplot(lw_combined, aes(x = a)) + geom_histogram() + labs(caption = "source: combined") + facet_wrap(~source, scales = "free"))
(combined_lna <- ggplot(lw_combined, aes(x = ln_a)) + geom_histogram() + labs(caption = "source: combined") + facet_wrap(~source))
(combined_b   <- ggplot(lw_combined, aes(x = b)) + geom_histogram() + labs(caption = "source: combined") + facet_wrap(~source))





# Remove the lw building components
rm(wigley_lw, fishbase_lw, spp_classes, wigley_lw)


# Export the key for the size spectra build
#write_csv(lw_combined, here::here("data/biomass_key_combined.csv"))






####  Loading the Combined Key for Investigation  ####
lw_combined <- read_csv(here::here("data/biomass_key_combined.csv"))

# Species with strange fits

lw_combined %>% 
  filter(comname == "spiny dogfish") %>% 
  View("Spiny dogfish")

lw_combined %>% 
  filter(scientific_name == "scophthalmus aquosus" | comname == "windowpane") %>% 
  ggplot(aes(comname, a)) +
  geom_col()






####  Double checking Fishbase coefficients  ####


fb_manual <- filter(lw_combined, source == "fishbase")
fb_species <- data.frame("ComName" = unique(fb_manual$comname))

#### Step 1: get scientific names and SpecCodes  ####
fb_sci <- rfishbase::common_to_sci(fb_species$ComName)


#### Step 2: Get length-weight coefficients  ####
fb_coefficients <- rfishbase::length_weight(species_list = fb_sci$Species)
fb_coefficients

ggplot() +
  geom_density(data = fb_coefficients, aes(a, color = "rfishbase")) + 
  geom_density(data = fb_manual, aes(a, color = "Manual Lookup"))


ggplot() +
  geom_density(data = fb_coefficients, aes(b, color = "rfishbase")) + 
  geom_density(data = fb_manual, aes(b, color = "Manual Lookup"))
