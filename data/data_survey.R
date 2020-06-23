###########################
## Extract Survye Data from master datasets
###########################

library(dplyr)

############################
## Setup directories
############################
setwd("C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/data/")

############################
## Import CPUE survey data
############################
data_wd = 'C:/Users/zack.oyafuso/Desktop/'

data = read.csv(paste0(data_wd, "data-raw/cpue_GOA_selected_spp.csv"), 
                stringsAsFactors = FALSE) # CPUE is (num or kg / km^2)

#############################
## join haul data to get coordinates, depth, bottom and surface temperature
#############################
haul <- read.csv(paste0(data_wd, "data-raw/haul.csv"), 
                 stringsAsFactors = FALSE)
haul <- cbind(haul, 
              # get haul midpoints
              geosphere::midPoint(cbind(haul$START_LONGITUDE, 
                                        haul$START_LATITUDE), 
                                  cbind(haul$END_LONGITUDE, 
                                        haul$END_LATITUDE))) 
haul$DATE <- as.Date(haul$START_TIME, format = "%d-%b-%y")
haul$MONTH <- lubridate::month(haul$DATE)
haul$DAY <- lubridate::day(haul$DATE)
haul <- haul %>% select(HAULJOIN, GEAR_DEPTH, SURFACE_TEMPERATURE, 
                        GEAR_TEMPERATURE, LATITUDE = lat, LONGITUDE = lon, 
                        DATE, DAY, MONTH)
data <- inner_join(data, haul)

###############################
## join species names
################################
species_codes =  read.csv(paste0(data_wd, "data-raw/species.csv"), 
                          stringsAsFactors = FALSE)
species_codes = select(species_codes, -YEAR_ADDED)
data <- inner_join(data, species_codes)

###############################
## select and rename columns, dropping rows with mising depths
###############################
data <- data %>% select(YEAR, SURVEY, BOTTOM_DEPTH = GEAR_DEPTH, 
                        SURFACE_TEMPERATURE, GEAR_TEMPERATURE, 
                        # CPUE = WGTCPUE,
                        EFFORT, WEIGHT,
                        LATITUDE, LONGITUDE, DATE, DAY, MONTH, 
                        SPECIES_NAME, COMMON_NAME) %>%
  tidyr::drop_na(BOTTOM_DEPTH,LATITUDE,LONGITUDE) 

###############################
# filter to GOA survey, remove tows with 0 bottom depth, and drop 2001 year when
# the survey was incomplete and years before 1996 when a different net/soak time
# was used
###############################
data <- data %>% filter(SURVEY == "GOA", 
                        BOTTOM_DEPTH > 0, YEAR != 2001 & YEAR >= 1996)

###############################
## sum catches of northern and southern rock sole with rock sole unid. 
## (not distinguished until 1996)
###############################
rock_soles <- data %>% dplyr::filter(COMMON_NAME %in% c("rock sole unid.", 
                                                        "southern rock sole", 
                                                        "northern rock sole")) %>%
  group_by_at(vars(-WEIGHT, -COMMON_NAME, -SPECIES_NAME)) %>%
  summarise(WEIGHT = sum(WEIGHT)) %>%
  ungroup() %>%
  mutate(SPECIES_NAME = "Lepidopsetta spp.", COMMON_NAME = "rock soles")

data <- as.data.frame(rbind(data, rock_soles))

################################
## sum catches of blackspooted and rougheye rocks with rougheye and 
## blackspotted rockfish unid. 
################################
B_R_rockfishes <- data %>% dplyr::filter(
  COMMON_NAME %in% c("blackspotted rockfish", 
                     "rougheye rockfish", 
                     "rougheye and blackspotted rockfish unid.")) %>%
  group_by_at(vars(-WEIGHT, -COMMON_NAME, -SPECIES_NAME)) %>%
  summarise(WEIGHT = sum(WEIGHT)) %>%
  ungroup() %>%
  mutate(SPECIES_NAME = "Sebastes B_R", COMMON_NAME = "B_R_rockfishes")
data <- as.data.frame(rbind(data, B_R_rockfishes))

##########################################
## scale observed bottom depth, calculate depth^2
##########################################
data$DEPTH = scale(x = log(data$BOTTOM_DEPTH))
data$DEPTH2 = data$DEPTH^2

# Filter species to make it easier to import later
# 1. Arrowtooth Flounder (Atherestes stomias, code 10110)
# 2. Pacific Cod (Gadus macrocephalus, code 21720)
# 3. Pacific Ocean Perch (Sebastes alutus, code 30060)
# 4. Sablefish (Anoplopoma fimbria, code 20510)
# 5. Walleye pollock (Gadus chalcogrammus, code 21740)
# 6. Dover sole (Solea solea, code 10180)
# 7. Pacific halibut (Hippoglossus stenolepis, code 10120)
# 8. Flathead sole (Hippoglossoides elassodon, code 10130)
# 9. Rex sole (Glyptocephalus zachirus, code 10200)
# 10. Dusky rockfish (Sebastes variabilis, code 30152)
# 11. Northern rockfish (Sebastes polyspinis, code 30420)
# 12. Silvergray Rockfish (Sebastes brevispinis, code 30100)
# 13. Yellowfin sole (Limanda aspera, code 10210)
# 14. Shortspine thornyhead (Sebastolobus alascanus, code 30020

# 15. Rougheye and blackspotted rockfishes (Sebastes aleutianus and Sebastes 
# melanostictus, respectively, codes 30050,30051,30052)

# 16/17 Northern and Southern rock sole (Lepidopsetta polyxystra and L.
# bilineata, respectivity, codes 10260,10261,10262)

data = subset(data,
              COMMON_NAME %in% c('Pacific ocean perch', 
                                 'arrowtooth flounder', 
                                 'Pacific cod', 
                                 'walleye pollock', 
                                 'Pacific halibut', 
                                 'rex sole', 
                                 'Dover sole',
                                 'flathead sole', 
                                 'sablefish', 
                                 'dusky rockfish',
                                 'northern rockfish',
                                 "northern rock sole", 
                                 "southern rock sole",
                                 'B_R_rockfishes',
                                 'shortspine thornyhead',
                                 'yellowfin sole',
                                 'silvergray rockfish'))

#Save
data = data[order(data$YEAR, data$SPECIES_NAME),]
write.csv(x = data, file = "data/GOA_multspp.csv", row.names = F)
