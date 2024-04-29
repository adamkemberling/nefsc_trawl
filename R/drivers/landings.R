####  Landings Data EDA  ####
# 6/7/2022



####  Packages  ####
{
  library(readxl)
  library(here)
  library(tidyverse)
  library(gmRi)
  library(janitor)
  library(gt)
  library(sf)
  library(targets)
  library(scales)
}


####  Data  ####

# Organized by port
by_port <- read_xlsx(
  path = here("data/landings_data/KMills_landings by area 1964-2021 - FINFISH ONLY_MAY 2022.xlsx"), 
  sheet = 1, 
  skip = 0) %>% 
  clean_names()

# Organized by stat zone
by_szone <- read_xlsx(
  path = here("data/landings_data/KMills_landings by area 1964-2021 - FINFISH ONLY_MAY 2022.xlsx"), 
  sheet = 2, 
  skip = 0) %>% 
  clean_names()

# List of species
spec_list <- read_xlsx(
  path = here("data/landings_data/KMills_landings by area 1964-2021 - FINFISH ONLY_MAY 2022.xlsx"), 
  sheet = 3, 
  skip = 0) %>% 
  clean_names()


# Landings of finfish* sheet 5
res_path <- cs_path("res")
landings <- read_xlsx(
  path = str_c(res_path, "GARFO_landings/KMills_landings by area 1964-2021_JUN 2022.xlsx"), sheet = 5) %>% 
  rename_all(tolower)

####  Data Exploration  ####


# 1. Aggregate Statzones by the region definitions of the project

# # Shapefiles for the fisheries stat zones
# stat_zones <- read_sf(str_c(res_path, "Shapefiles/Statistical_Areas/Statistical_Areas_2010_withNames.shp"))


# Make a list of zones to roughly match the survey areas:
fish_zones <- list(
  "Gulf of Maine"        = c(511:515, 464, 465),
  "Georges Bank"         = c(521, 522, 525, 561, 562),
  "Southern New England" = c(611, 612, 613, 616, 526, 537, 538, 539),
  "Mid-Atlantic Bight"   = c(614:615, 621, 622, 625, 626, 631, 632))


# Join to landings
# Add the labels into the landings data and remove what we don't need there:
area_levels_long <- c("Northeast Shelf", "Gulf of Maine", "Georges Bank", "Southern New England", "Mid-Atlantic Bight")
landings <- landings %>% 
  mutate(
    survey_area = case_when(
      `stat area` %in% fish_zones$"Gulf of Maine" ~ "Gulf of Maine",
      `stat area` %in% fish_zones$"Georges Bank" ~ "Georges Bank",
      `stat area` %in% fish_zones$"Southern New England" ~ "Southern New England",
      `stat area` %in% fish_zones$"Mid-Atlantic Bight" ~ "Mid-Atlantic Bight")) %>% 
  filter(survey_area %in% c("Georges Bank", "Gulf of Maine", "Southern New England", "Mid-Atlantic Bight")) %>% 
  mutate(survey_area = factor(survey_area, area_levels_long)) %>% 
  filter(between(year, 1960, 2019),
         !str_detect(sppname, "CONFIDENTIAL")) %>% 
  mutate(decade = floor_decade(year),
         sppname = str_to_title(sppname))
  



####  Landings Breakdowns  ####
# Goal:
# Separate landings into groups that might better explain in a mechanistic sense
# how they would impact a trawl-sampled community



# Add the labels into the landings data and remove what we don't need there:
landings %>% 
  group_by(survey_area, decade,  sppname) %>% 
  summarise(avg_landings_lb = mean(`landed lbs`, na.rm = T),
            across(.cols = c(landings_lb = `landed lbs`, 
                             value = value), 
                   .fns = sum, 
                   .names = "total_{.col}"),
            .groups = "drop") %>% 
  group_by(survey_area, decade) %>% 
  slice_max(order_by = total_landings_lb, n = 3) %>% 
  gt() %>% 
  gt::tab_header(title = md("**Top Commercial Fisheries Landings of Northeastern US (by weight)**")) %>% 
  tab_stubhead(label = md("*Harvest Region*")) %>% 
  fmt_number(columns = c(avg_landings_lb, total_landings_lb, total_value),
             use_seps = T, 
             sep_mark = ",",
             suffixing = T) %>% 
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_row_groups()) %>% 
  cols_label(
    decade            = md("*Decade*"),
    sppname           = md(""),
    avg_landings_lb   = md("*Avg. Annual Landings (lb.)*"),
    total_landings_lb = md("*Total Landings (lb.)*"),
    total_value       = md("*Total Value ($)*")) %>%
  gt::tab_source_note(source_note = md("*Landings data obtained from the Greater Atlantic Regional Fishing Office (GARFO)*"))



# Maybe a plot is better:

# Total landings by region
landings_annual <- landings %>% 
  group_by(survey_area, year) %>% 
  summarise(
    avg_landings_lb = mean(`landed lbs`, na.rm = T),
    across(.cols = c(landings_lb = `landed lbs`, 
                     value = value), 
           .fns = sum, 
           .names = "total_{.col}"),
    .groups = "drop")

landings_annual %>% 
  ggplot() +
  geom_line(aes(year, total_landings_lb), show.legend = F) +
  facet_wrap(~survey_area) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(x = "Year", y = "All Species Yearly Total Landings (lb.)")


# Does Georges Bank have an NA this time?
landings %>% filter(survey_area == "Georges Bank", is.na(`landed lbs`)) 


# Assign groundfish label (Jonathan Labaree, market based)
gf_species <- tolower(
  c("cod, atlantic", 
    "Flounder, American Plaice", 
    "Flounder, Winter", 
    "Flounder, witch",
    "flounder, yellowtail", 
    "haddock", 
    "hake, white", 
    "hake, red", 
    "hake, silver", 
    "halibut, atlantic", 
    "pollock", 
    "redfish, acadian"))


# What does that give us if the groundfish are labeled?
# Split out species and label groundfish
landings_species <-  landings %>% 
  expand(year, survey_area, sppname) %>% 
  left_join(landings) %>% 
  filter(survey_area!= "Northeast Shelf") %>% 
  group_by(survey_area, year,  sppname) %>% 
  summarise(avg_landings_lb = mean(`landed lbs`, na.rm = T),
            across(.cols = c(landings_lb = `landed lbs`, 
                             value = value), 
                   .fns = sum, 
                   .names = "total_{.col}"),
            .groups = "drop") %>% 
  mutate(across(4:6, ~ifelse(is.na(.x), 0, .x)),
         is_gf = ifelse(tolower(sppname) %in% gf_species, "groundfish", "other"))



# Plots
landings_species %>% 
  ggplot() +
  geom_line(aes(year, total_landings_lb, group = sppname, color = is_gf), 
            show.legend = F) +
  facet_grid(is_gf~survey_area, scales = "free_y") +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(x = "Year", y = "Single Species Yearly Total Landings (lb.)")





####  Supplementing Landings information  ####




##### 1.  Limiting Landings to species in the spectrum  ####

# The following table has all the species names from the GARFO landings
# in an attempt to add context, I have been trying to match these names 
# to the names used in the hare 2016 paper about functional groupd
# another way to add context could be with fishbase

# This are quickly matched by hand so we can move forward

# landings types:
landings_match <- tribble(
  ~"sppname",                      ~"comname",                 ~"fishery_type" ,
  # landings name                  # trawl survey name  # is commercial?
  "Alewife"                        ,"alewife"                  ,NA,
  # "Barracudas"                     ,                         ,NA,
  "Bass, Striped"                  ,"striped bass"             ,NA,
  "Bay Anchovy"                    ,"anchovies"                ,NA,
  "Bluefish"                       ,"bluefish"                 ,NA,
  "Bonito, Atlantic"               ,"atlantic bonito"          ,NA,
  #"Bullhead"                       ,                          ,NA,
  "Butterfish"                     ,"butterfish"               ,NA,
  "Butterfish, Gulf"               ,"butterfish"               ,NA,
  #"Carp, Common"                   ,                          ,NA,
  #"Catfish, Blue"                  ,                          ,NA,
  #"Catfish, Freshwater"            ,                          ,NA,
  #"Catfish, Sea"                   ,                          ,NA,
  "Cobia"                          ,"cobia"                    ,NA,
  "Cod, Atlantic"                  ,"atlantic cod"             ,NA,
  #"Crappie"                        ,                          ,NA,
  "Croaker, Atlantic"              ,"Atlantic Croaker"  ,      NA,
  "Cunner"                         ,"cunner"                  ,NA,
  "Cusk"                           ,"cusk"                    ,NA,
  "Dolphinfish"                    ,"dolphinfish"             ,NA,
  # "Dory, American John"            ,""                      ,NA,
  "Drum, Black"                    ,"black drum"              ,NA,
  #"Drum, Nk"                       ,                         ,NA,
  "Drum, Red"                      ,"Red Drum"                ,NA,
  "Eel, American"                  ,"american eel"            ,NA,
  "Eel, Conger"                    ,"conger eel"              ,NA,
  #"Eel, Nk"                        ,                         ,NA,
  # "Eel, Sand/Launce"               ,                        ,NA,
  # "Escolar"                        ,                        ,NA,
  "Flounder, American Plaice"      ,"american plaice"         ,NA,
  "Flounder, Fourspot"             ,"fourspot flounder"       ,NA,
  # "Flounder, Nk"                   ,                        ,NA,
  "Flounder, Sand-Dab/Windowpane"  , "windowpane flounder"    ,NA,
  "Flounder, Southern"             , "southern flounder"      ,NA,
  "Flounder, Summer"               ,"summer flounder"         ,NA,
  "Flounder, Winter"               ,"winter flounder"         ,NA,
  "Flounder, Witch"                ,"witch flounder"          ,NA,
  "Flounder, Yellowtail"           ,"yellowtail flounder"     ,NA,
  # "Garfish"                        ,                        ,NA,
  # "Grenadier/Rat Tail"             ,                        ,NA,
  # "Grouper, Nk"                    ,                        ,NA,
  # "Grunt, Haemulidae"              ,                        ,NA,
  "Haddock"                        , "haddock"                ,NA,
  "Hagfish"                        , "atlantic hagfish"       ,NA,
  "Hake, Offshore"                 , "Offshore Hake"          ,NA,
  "Hake, Red"                      , "red hake"               ,NA,
  "Hake, Red & White Mix"          , "red hake"               ,NA,
  "Hake, Silver"                   , "silver hake"            ,NA,
  "Hake, White"                    , "white hake"             ,NA,
  "Halibut, Atlantic"              , "atlantic halibut"       ,NA,
  # "Halibut, Greenland"             ,                        ,NA,
  # "Harvestfish"                    ,                        ,NA,
  "Herring, Atlantic"              ,"atlantic herring"        ,NA,
  "Herring, Blue Back"             ,"blueback herring"        ,NA,
  # "Herring, Nk"                    ,                        ,NA,
  # "Hogfish"                        ,                        ,NA,
  "Jack, Amberjack, Nk"            , "greater amberjack"      ,NA,
  # "Jack, Crevalle"                 ,                        ,NA,
  "Kingfish, Northern"             ,"Northern kingfish"       ,NA,
  # "Ladyfish"                       ,                        ,NA,
  # "Lumpfish"                       ,                        ,NA,
  "Mackerel, Atlantic"             ,"Atlantic mackerel"       ,NA,
  # "Mackerel, Chub"                 ,                        ,NA,
  "Mackerel, Frigate"              ,"frigate mackerel"        ,NA,
  "Mackerel, King"                 ,"king mackerel"           ,NA,
  "Mackerel, Spanish"              ,"Spanish mackerel"        ,NA,
  # "Marlin, Blue"                   ,                        ,NA,
  # "Marlin, Nk"                     ,                        ,NA,
  # "Marlin, White"                  ,                        ,NA,
  "Menhaden, Atlantic"             ,"Atlantic Menhaden"       ,NA,
  "Monkfish/Angler/Goosefish"      ,"goosefish"               ,NA,
  "Mullet, Striped"                , "striped mullet"         ,NA,
  "Mullet, White"                  , "white mullet"           ,NA,
  "Mummichog"                      , "mummichog"              ,NA,
  "Needlefish, Atlantic"           , "atlantic needlefish"    ,NA,
  # "Opah"                           ,                        ,NA,
  # "Other Fish, Bony"               ,                        ,NA,
  # "Other Fish, Groundfish"         ,                        ,NA,
  # "Other Fish, Pelagics"           ,                        ,NA,
  "Perch, Sand"                    , "sand perch"             ,NA,
  "Perch, White"                   , "white perch"            ,NA,
  "Perch, Yellow"                  , "yellow perch"           ,NA,
  "Pigfish"                        , "pigfish"                ,NA,
  "Pollock"                        , "pollock"                ,NA,
  # "Pompano, Common"                ,                        ,NA,
  # "Pompano, Killifish"             ,                        ,NA,
  "Porgy, Red"                     , "red porgy"              ,NA,
  "Pout, Ocean"                    , "Ocean Pout"             ,NA,
  # "Puffer, Northern"               ,                        ,NA,
  "Ray, Cownose"                   , "cownose ray"            ,NA,
  "Redfish, Acadian"               , "acadian redfish"        ,NA,
  # "Ribbonfish/Dealfish"            ,                        ,NA,
  "Rosefish, Black Bellied"        , "blackbelly rosefish"    ,NA,
  "Runner, Blue"                   , "blue runner"            ,NA,
  "Salmon, Atlantic"               , "atlantic salmon"        ,NA,
  "Sculpin"                        , "sculpin"                ,NA,
  "Scup"                           , "scup"                   ,NA,
  "Sea Bass, Black"                , "Black Sea Bass"         ,NA,
  "Sea Raven"                      , "sea raven"              ,NA,
  # "Searobin, Nk"                   ,                        ,NA,
  "Shad, American"                 , "american shad"          ,NA,
  "Shad, Gizzard"                  , "gizzard shad"           ,NA,
  "Shad, Hickory"                  , "hickory shad"           ,NA,
  "Shark, Atl Sharpnose"           , "atlantic sharpnose shark",NA,
  "Shark, Bignose"                 , "bignose shark"          ,NA,
  "Shark, Black Tip"               , "blacktip shark"         ,NA,
  "Shark, Blue"                    , "blue shark"             ,NA,
  "Shark, Bonnethead"              , "bonnethead shark"       ,NA,
  "Shark, Bull"                    , "bull shark"             ,NA,
  # "Shark, Dogfish, Nk"             ,                        ,NA,
  "Shark, Dogfish, Smooth"         , "smooth dogfish"         ,NA,
  "Shark, Dogfish, Spiny"          , "spiny dogfish"          ,NA,
  "Shark, Dusky"                   , "dusky shark"            ,NA,
  # "Shark, Hammerhead, Nk"          ,                        ,NA,
  # "Shark, Large Coastal"           ,                        ,NA,
  "Shark, Lemon"                   , "lemon shark"            ,NA,
  "Shark, Mako, Longfin"           , "longfin mako"           ,NA,
  #"Shark, Mako, Nk"                ,                 ,       ,NA,
  "Shark, Mako, Shortfin"          , "shortfin mako"          ,NA,
  "Shark, Night"                   , "night shark"            ,NA,
  # "Shark, Nk"                      ,                        ,NA,
  "Shark, Nurse"                   , "nurse shark"            ,NA,
  # "Shark, Other Pelagic"           ,                        ,NA,
  "Shark, Porbeagle"               , "porbeagle"              ,NA,
  "Shark, Sand Tiger"              , "sand tiger"             ,NA,
  "Shark, Sandbar"                 , "sandbar shark"          ,NA,
  "Shark, Silky"                   , "silky shark"            ,NA,
  "Shark, Spinner"                 , "spinner shark"          ,NA,
  "Shark, Thresher, Bigeye"        , "bigeye thresher shark"  ,NA,
  # "Shark, Thresher, Nk"            ,                        ,NA,
  "Shark, Tiger"                   , "tiger shark"            ,NA,
  "Shark, White"                   , "white shark"            ,NA,
  "Sheepshead, Atlantic"           , "sheepshead"             ,NA,
  "Silverside, Atlantic"           , "atlantic silverside"    ,NA,
  # "Silverside, Nk"                 ,                        ,NA,
  "Skate, Little"                  , "little skate"           ,NA,
  # "Skate, Nk"                      ,                        ,NA,
  "Skate, Winter"                  , "winter skate"           ,NA,
  "Smelt"                          , "rainbow smelt"          ,NA,
  # "Snakehead, Northern"            , "northern snakehead"   ,NA,
  # "Snapper, Nk"                    ,                        ,NA,
  "Snapper, Red"                   , "red snapper"            ,NA,
  "Spadefish"                      , "atlantic spadefish"     ,NA,
  "Spot"                           , "spot"                   ,NA,
  # "Starfish"                       ,                        ,NA,
  "Sturgeon, Atlantic"             , "atlantic sturgeon"      ,NA,
  # "Sturgeon, Nk"                   ,                        ,NA,
  # "Suckerfish, Nk"                 ,                        ,NA,
  "Sunfish"                        , "sunfish"                ,NA,
  "Swordfish"                      , "swordfish"              ,NA,
  "Tarpon"                         , "tarpon"                 ,NA,
  "Tautog"                         , "tautog"                 ,NA,
  "Tilefish, Blueline"             , "tilefish"               ,NA,
  "Tilefish, Golden"               , "tilefish"               ,NA,
  # "Tilefish, Nk"                   ,                        ,NA,
  "Tom Cod"                        , "tomcod"                 ,NA,
  "Triggerfish"                    , "triggerfish"            ,NA,
  "Triggerfish, Gray"              , "gray triggerfish"       ,NA,
  "Tripletail"                     , "tripletail"             ,NA,
  "Tuna, Albacore"                 , "albacore tuna"          ,NA,
  "Tuna, Bigeye"                   , "bigeye tuna"            ,NA,
  "Tuna, Bluefin"                  , "bluefin tuna"           ,NA,
  "Tuna, Little Tunny"             , "little tunny"           ,NA,
  # "Tuna, Nk"                       ,                        ,NA,
  "Tuna, Skipjack"                 , "skipjack tuna"          ,NA,
  "Tuna, Yellowfin"                , "yellowfin tuna"         ,NA,
  # "Wahoo"                          ,                        ,NA,
  "Weakfish/Sea Trout, Spotted"    , "spotted seatrout"       ,NA,
  "Weakfish/Sea Trout, Squeteague" , "weakfish"               ,NA,
  # "Whiting, King"                  ,                        ,NA,
  "Wolffish, Atlantic"             , "atlantic wolffish"      ,NA,
  # "Wreckfish"                                       
  )  %>% mutate(comname = tolower(comname))
  



# Using it
# 1. Biological data used as input
tar_load(catch_log2_labelled)  
gom_landings <- landings_yrly %>% filter(survey_area == "Gulf of Maine") 
gom_trawl_only <- gom_landings %>% 
  inner_join(landings_match) %>% 
  filter(comname %in% tolower(unique(catch_log2_labelled$comname)))

gom_all_summ <- group_by(gom_landings, year, survey_area) %>% summarise(total_lb = sum(total_landings_lb))
nohaden_summ <- gom_landings %>% filter(sppname != "Menhaden, Atlantic") %>% group_by(year, survey_area) %>% summarise(total_lb = sum(total_landings_lb))
gom_trawl_summ <- group_by(gom_trawl_only, year, survey_area) %>% summarise(total_lb = sum(total_landings_lb))


ggplot() +
  geom_line(data = gom_all_summ, aes(year, total_lb, color = "All Species")) +
  geom_line(data = nohaden_summ, aes(year, total_lb, color = "Menhaden Removed")) +
  geom_line(data = gom_trawl_summ, aes(year, total_lb, color = "Trawl Sampled Only"))  +
  scale_y_continuous(labels = label_comma()) +
  facet_wrap(~survey_area) +
  labs(y = "Finfish Landings")



##### 2.  Adding Functional Group Information  ####



# From Hare 2016, functional groups - ty data pasta
hare_functional_groups <- tibble::tribble(
                   ~hare_group,         ~hare_common_name,                          ~scientific_name,
                "Coastal Fish",        "Atlantic Croaker",                 "Micropogonias undulates",
                "Coastal Fish",       "Atlantic Menhaden",                     "Brevoortia tyrannus",
                "Coastal Fish",          "Black Sea Bass",                   "Centropristis striata",
                "Coastal Fish",       "Northern Kingfish",                  "Menticirrhus saxatilis",
                "Coastal Fish",                "Red Drum",                     "Sciaenops ocellatus",
                "Coastal Fish",                    "Scup",                     "Stenotomus chrysops",
                "Coastal Fish",        "Spanish Mackerel",                 "Scomberomorus maculatus",
                "Coastal Fish",                    "Spot",                    "Leiostomus xanthurus",
                "Coastal Fish",        "Spotted Seatrout",                     "Cynoscion nebulosus",
                "Coastal Fish",            "Striped Bass",                        "Morone saxatilis",
                "Coastal Fish",         "Summer Flounder",                   "Paralichthys dentatus",
                "Coastal Fish",                  "Tautog",                          "Tautoga onitis",
                "Coastal Fish",                "Weakfish",                       "Cynoscion regalis",
                "Coastal Fish",         "Winter Flounder",           "Pseudopleuronectes americanus",
             "Diadromous Fish",                 "Alewife",                    "Alosa pseudoharengus",
             "Diadromous Fish",              "Conger Eel",                       "Anguilla oceanica",
             "Diadromous Fish",            "American Eel",                       "Anguilla rostrata",
             "Diadromous Fish",           "American Shad",                       "Alosa sapidissima",
             "Diadromous Fish",         "Atlantic Salmon",                             "Salmo salar",
             "Diadromous Fish",       "Atlantic Sturgeon",                   "Acipenser oxyrhynchus",
             "Diadromous Fish",        "Blueback Herring",                        "Alosa aestivalis",
             "Diadromous Fish",            "Hickory Shad",                         "Alosa mediocris",
             "Diadromous Fish",           "Rainbow Smelt",                          "Osmerus mordax",
             "Diadromous Fish",      "Shortnose Sturgeon",                  "Acipenser brevirostrum",
               "Elasmobranchs",          "Barndoor Skate",                         "Dipturus laevis",
               "Elasmobranchs",         "Clearnose Skate",                         "Raja eglanteria",
               "Elasmobranchs",             "Dusky Shark",                   "Carcharhinus obscurus",
               "Elasmobranchs",            "Little Skate",                      "Leucoraja erinacea",
               "Elasmobranchs",               "Porbeagle",                             "Lamna nasus",
               "Elasmobranchs",           "Rosette Skate",                       "Leucoraja garmani",
               "Elasmobranchs",              "Sand Tiger",                       "Carcharias taurus",
               "Elasmobranchs",          "Smooth Dogfish",                          "Mustelus canis",
               "Elasmobranchs",            "Smooth Skate",                        "Malacoraja senta",
               "Elasmobranchs",           "Spiny Dogfish",                       "Squalus acanthias",
               "Elasmobranchs",            "Thorny Skate",                       "Amblyraja radiata",
               "Elasmobranchs",            "Winter Skate",                      "Leucoraja ocellata",
                  "Groundfish",         "Acadian Redfish",                      "Sebastes fasciatus",
                  "Groundfish",         "American Plaice",            "Hippoglossoides platessoides",
                  "Groundfish",            "Atlantic Cod",                            "Gadus morhua",
                  "Groundfish",        "Atlantic Hagfish",                        "Myxine glutinosa",
                  "Groundfish",        "Atlantic Halibut",               "Hippoglossus hippoglossus",
                  "Groundfish",       "Atlantic Wolffish",                        "Anarhichas lupus",
                  "Groundfish",                    "Cusk",                           "Brosme brosme",
                  "Groundfish",                 "Haddock",                "Melanogrammus aeglefinus",
                  "Groundfish",    "Monkfish (Goosefish)",                      "Lophius americanus",
                  "Groundfish",              "Ocean Pout",                      "Zoarces americanus",
                  "Groundfish",           "Offshore Hake",                      "Merluccius albidus",
                  "Groundfish",                 "Pollock",                       "Pollachius virens",
                  "Groundfish",                "Red Hake",                         "Urophycis chuss",
                  "Groundfish",             "Silver Hake",                   "Merluccius bilinearis",
                  "Groundfish",                "Tilefish",           "Lopholatilus chamaeleonticeps",
                  "Groundfish",              "White Hake",                        "Urophycis tenuis",
                  "Groundfish",              "Windowpane",                    "Scophthalmus aquosus",
                  "Groundfish",          "Witch Flounder",              "Glyptocephalus cynoglossus",
                  "Groundfish",     "Yellowtail Flounder",                      "Limanda ferruginea",
"Pelagic Fish and Cephalopods",               "Anchovies",      "Anchoa hepsetus / Anchoa mitchilli",
"Pelagic Fish and Cephalopods",        "Atlantic Herring",                         "Clupea harengus",
"Pelagic Fish and Cephalopods",       "Atlantic Mackerel",                        "Scomber scombrus",
"Pelagic Fish and Cephalopods",          "Atlantic Saury",                      "Scomberesox saurus",
"Pelagic Fish and Cephalopods",                "Bluefish",                     "Pomatomus saltatrix",
"Pelagic Fish and Cephalopods",              "Butterfish",                    "Peprilus triacanthus",
"Pelagic Fish and Cephalopods",   "Longfin Inshore Squid",                     "Doryteuthis pealeii",
"Pelagic Fish and Cephalopods",             "Sand Lances", "Ammodytes americanus & Ammodytes dubius",
"Pelagic Fish and Cephalopods", "Northern Shortfin Squid",                      "Illex illecebrosus",
       "Benthic Invertebrates",        "American Lobster",                      "Homarus americanus",
       "Benthic Invertebrates",    "Atlantic Sea Scallop",                "Placopecten magellanicus",
       "Benthic Invertebrates",       "Atlantic Surfclam",                     "Spisula solidissima",
       "Benthic Invertebrates",             "Bay Scallop",                    "Argopecten irradians",
       "Benthic Invertebrates",               "Bloodworm",                    "Glycera dibranchiata",
       "Benthic Invertebrates",               "Blue Crab",                     "Callinectes sapidus",
       "Benthic Invertebrates",             "Blue Mussel",                          "Mytilus edulis",
       "Benthic Invertebrates",            "Cancer Crabs",      "Cancer borealis / Cancer irroratus",
       "Benthic Invertebrates",         "Channeled Whelk",               "Busycotypus canaliculatus",
       "Benthic Invertebrates",       "Deep-sea Red Crab",                     "Chaceon quinquedens",
       "Benthic Invertebrates",          "Eastern Oyster",                   "Crassostrea virginica",
       "Benthic Invertebrates",        "Green Sea Urchin",       "Strongylocentrotus droebachiensis",
       "Benthic Invertebrates",          "Horseshoe Crab",                      "Limulus polyphemus",
       "Benthic Invertebrates",           "Knobbed Whelk",                          "Busycon carica",
       "Benthic Invertebrates",         "Northern Shrimp",                       "Pandalus borealis",
       "Benthic Invertebrates",            "Ocean Quahog",                       "Arctica islandica",
       "Benthic Invertebrates",         "Northern Quahog",                   "Mercenaria mercenaria",
       "Benthic Invertebrates",          "Softshell Clam",                            "Mya arenaria"
  )



#### Modeling Prep by Area: ####
# Aggregate the landings by survey regions:


# Get Summaries
landings_summ <- landings %>% 
  drop_na() %>% 
  rename(
    "weight_lb" = `landed lbs`,
    "live_lb" = `live lbs`) %>% 
  group_by(year, survey_area) %>% 
  summarise( 
    across(
      .cols = c(value, weight_lb, live_lb), 
      .fns = list(mean = ~mean(.x , na.rm = T), 
                  total = ~sum(.x , na.rm = T)), 
      .names = "{.fn}_{.col}"), 
    .groups = "drop")  


# Scale the landings to create an index by area
landings_summ <- landings_summ %>% 
  group_by(survey_area) %>% 
  mutate(total_wt_z = as.numeric(base::scale(total_weight_lb))) %>% 
  ungroup()


# Prepare timeseries for use in modeling
landings_idx <- landings_summ %>% 
  group_by(year, area = survey_area) %>% 
  summarise(
    value = mean(total_weight_lb, na.rm = T),
    #value = mean(total_wt_z, na.rm = T),
    .groups = "drop") %>% 
  mutate(var = "landings")



