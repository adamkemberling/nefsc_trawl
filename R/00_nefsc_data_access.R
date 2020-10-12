####  NEFSC Trawl Data Access  ####



####  Packages  ####
library(janitor)
library(here)
library(gmRi)
library(tidyverse)
library(patchwork)


####  Data  ####

# 2019 data From Janet Nye
mills_path <- shared.path(group = "Mills Lab", folder = "")
#load(str_c(mills_path, "Data/Survdat_Nye_allseason.RData"))

# Gonna save it locally for convenience/laziness
#write_csv(survdat, here("data/NEFSC/2019survdat_nye.csv"))

# 2020 edition
load(here("data/NEFSC/Survdat_Nye_Aug 2020.RData"))

# Gonna save it locally for convenience/laziness
# write_csv(survdat, here("data/NEFSC/2020survdat_nye.csv"))



####____________________####


####  Combined L-W Key  ####
# Sources:
# Wigley Paper
# Fishbase

res_path <- shared.path(os.use = "unix", group = "RES Data", folder = NULL)


# Fishbase Growth Coefficients
nefsc_lw <- read_csv(here("data/NEFSC/nefsc_lw_key_filled.csv"),
                     guess_max = 1e3,
                     col_types = cols())


# Wrigley Paper, Load length weight coefficients
load(paste0(res_path, "/NMFS_trawl/lwreg.Rdata"))


# Species class/groups based on life history
spp_classes <- read_csv(here("data/kmills/sppclass.csv"),
                        col_types = cols())



####  Growth Coefficients Prep  ####

####__ 1.  Fishbase Coefficients  ####
# Fill in NA growth coefficients using related species (typically same genera)
nefsc_lw <- nefsc_lw %>% 
  mutate(a = as.numeric(a),
         b = as.numeric(b),
         ln_a = log(a),
         related_ln_a = log(related_a),
         ln_a = ifelse(is.na(ln_a), related_ln_a, ln_a),
         a = ifelse(is.na(a), exp(ln_a), a))

(nefsc_a <- ggplot(nefsc_lw, aes(x = a)) + geom_histogram() + labs(caption = "source: fishbase"))
(nefsc_lna <- ggplot(nefsc_lw, aes(x = ln_a)) + geom_histogram() + labs(caption = "source: fishbase"))
(nefsc_b <- ggplot(nefsc_lw, aes(x = b)) + geom_histogram() + labs(caption = "source: fishbase"))



nefsc_lw %>% 
  mutate(lna_check = log(a)) %>% 
  ggplot(aes(ln_a, lna_check)) +
  geom_point()


####__ 2. Wrigley Paper Coefficients  ####
lwreg <- lwreg %>% clean_names()  
lwreg <- lwreg %>%  
  mutate(comname = str_to_lower(common_name), .before = scientific_name,
         common_name = NULL,
         scientific_name = str_to_lower(scientific_name),
         svspp = str_pad(svspp, width = 3, side = "left", pad = "0")) %>% 
  mutate(ln_a = lna, 
         a = exp(ln_a),
         lna = NULL) %>% 
  select(svspp, comname, scientific_name, season, catchsex, a, ln_a, everything())



####__ 3. Classes & Economic Groups  ####
# Species Life History & Fisheries Relevance
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
wigley_lw <- left_join(lwreg, spp_classes, by = c("svspp", "comname", "scientific_name"))



(wigley_a   <- ggplot(wigley_lw, aes(x = a)) + geom_histogram() + labs(caption = "source: wigley"))
(wigley_lna <- ggplot(wigley_lw, aes(x = ln_a)) + geom_histogram() + labs(caption = "source: wigley"))
(wigley_b   <- ggplot(wigley_lw, aes(x = b)) + geom_histogram() + labs(caption = "source: wigley"))

wigley_lw %>% 
  mutate(lna_check = log(a)) %>% 
  ggplot(aes(ln_a, lna_check)) +
  geom_point()



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
nefsc_lw  <- nefsc_lw %>% mutate(source = "fishbase") %>% select(source, everything())
lw_combined <- full_join(nefsc_lw, wigley_lw, by = c("source", "comname", "b", "a", "ln_a")) 

# Fill in gaps for things that should be consistent
# Lot of matching by common name and pulling the first value where it matches to fill the NA values
lw_combined <- lw_combined %>% 
  arrange(comname) %>% 
  mutate(
    hare_group = ifelse(is.na(hare_group),
                        nefsc_lw[ match(comname,  nefsc_lw$comname)[1], "hare_group"][[1]],
                        hare_group),
    svspp = ifelse(is.na(svspp),
                   wigley_lw[ match(comname, wigley_lw$comname)[1], "svspp"][[1]],
                   svspp),
    scientific_name = ifelse(is.na(scientific_name),
                             wigley_lw[ match(comname, wigley_lw$comname)[1], "scientific_name"][[1]],
                             scientific_name),
    fishery = ifelse(is.na(fishery),
                     wigley_lw[ match(comname, wigley_lw$comname)[1], "fishery"][[1]],
                     fishery),
    spec_class = ifelse(is.na(spec_class),
                        wigley_lw[ match(comname, wigley_lw$comname)[1], "spec_class"][[1]],
                        spec_class))

# Grab and order columns we need
lw_combined <- lw_combined %>% 
  select(source, svspp, comname, scientific_name, spec_class, hare_group, fishery, season, catchsex,
         b, a, ln_a, units, related_species, related_a, related_b, related_ln_a,
         wigleymin, wigleymax, wigleyvalue)


####____  test: compare value ranges   ####

# Double check that the values are in the same ballpark:
(combined_a   <- ggplot(lw_combined, aes(x = a)) + geom_histogram() + labs(caption = "source: combined") + facet_wrap(~source))
(combined_lna <- ggplot(lw_combined, aes(x = ln_a)) + geom_histogram() + labs(caption = "source: combined") + facet_wrap(~source))
(combined_b   <- ggplot(lw_combined, aes(x = b)) + geom_histogram() + labs(caption = "source: combined") + facet_wrap(~source))





# Remove the lw building components
rm(wigley_lw, nefsc_lw, spp_classes, lwreg)


# Export the key for the size spectra build
write_csv(lw_combined, here::here("data/biomass_key_combined.csv"))

