###############################################################################
## Project:         GOA Groundfish CPUE Data Synthesis 
## Author:          Zack Oyafuso (zack.oyafuso@noaa.gov)
## Contributors:    Lewis Barnett (lewis.barnett@noaa.gov)
## Description:     Create CPUE dataset used for VAST for species of interest
###############################################################################

##################################################
#### Require Packages
##################################################
library(dplyr)

##################################################
####   Set up directories
##################################################
rm(list = ls())
# setwd("C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/data/")
data_wd = "C:/Users/Zack Oyafuso/Desktop/AK_BTS/"
github_wd = "C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/data/"

##################################################
#### Import CPUE survey data
#### Import EFH bathhymetry raster
##################################################
data = read.csv(paste0(data_wd, "data-raw/cpue_GOA_selected_spp.csv"), 
                stringsAsFactors = FALSE) # CPUE is (num or kg / km^2)

bathy <- raster::raster(
  paste0(github_wd, "EFH_bathymetry/aigoa_bathp1c/dblbnd.adf"))

##################################################
#### Join haul data to get coordinates, depth, bottom and surface temperature
##################################################
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

##################################################
####   Join species names
##################################################
species_codes =  read.csv(paste0(data_wd, "data-raw/species.csv"), 
                          stringsAsFactors = FALSE)
species_codes = select(species_codes, -YEAR_ADDED)
data <- inner_join(data, species_codes)

##################################################
####  Select and rename columns, dropping rows with mising depths
##################################################
data <- data %>% select(YEAR, SURVEY, 
                        BOTTOM_DEPTH = GEAR_DEPTH,
                        SURFACE_TEMPERATURE, GEAR_TEMPERATURE, 
                        # CPUE = WGTCPUE,
                        EFFORT, WEIGHT,
                        LATITUDE, LONGITUDE, DATE, DAY, MONTH, 
                        SPECIES_NAME, COMMON_NAME) %>%
  tidyr::drop_na(BOTTOM_DEPTH,LATITUDE,LONGITUDE) 

##################################################
####  Filter to GOA survey, remove tows with 0 bottom depth, and drop 2001,
####  the year when the survey was incomplete and years before 1996 when a 
####  different net/soak time was used
##################################################
data <- data %>% filter(SURVEY == "GOA", 
                        BOTTOM_DEPTH > 0, 
                        YEAR != 2001 & YEAR >= 1996)

##################################################
####  Sum catches of northern and southern rock sole with rock sole unid.
#### (not distinguished until 1996), rename species complex
##################################################
rock_soles <- data %>% 
  dplyr::filter(COMMON_NAME %in% c("rock sole unid.", 
                                   "southern rock sole", 
                                   "northern rock sole")) %>%
  group_by_at(vars(-WEIGHT, -COMMON_NAME, -SPECIES_NAME)) %>%
  summarise(WEIGHT = sum(WEIGHT)) %>%
  ungroup() %>%
  mutate(SPECIES_NAME = "Lepidopsetta spp.", COMMON_NAME = "rock soles")

data <- as.data.frame(rbind(data, rock_soles))

##################################################
####  Sum catches of blackspooted and rougheye rocks with rougheye and 
####  blackspotted rockfish unid.,  rename species complex
##################################################
B_R_rockfishes <- data %>% dplyr::filter(
  COMMON_NAME %in% c("blackspotted rockfish", 
                     "rougheye rockfish", 
                     "rougheye and blackspotted rockfish unid.")) %>%
  group_by_at(vars(-WEIGHT, -COMMON_NAME, -SPECIES_NAME)) %>%
  summarise(WEIGHT = sum(WEIGHT)) %>%
  ungroup() %>%
  mutate(SPECIES_NAME = "Sebastes B_R", COMMON_NAME = "B_R_rockfishes")
data <- as.data.frame(rbind(data, B_R_rockfishes))

##################################################
#### Filter species to make it easier to import later
#### 1. Arrowtooth Flounder (Atherestes stomias, code 10110)
#### 2. Pacific Cod (Gadus macrocephalus, code 21720)
#### 3. Pacific Ocean Perch (Sebastes alutus, code 30060)
#### 4. Walleye pollock (Gadus chalcogrammus, code 21740)
#### 5. Dover sole (Solea solea, code 10180)
#### 6. Pacific halibut (Hippoglossus stenolepis, code 10120)
#### 7. Flathead sole (Hippoglossoides elassodon, code 10130)
#### 8. Rex sole (Glyptocephalus zachirus, code 10200)
#### 0. Dusky rockfish (Sebastes variabilis, code 30152)
#### 10. Northern rockfish (Sebastes polyspinis, code 30420)
#### 11. Silvergray Rockfish (Sebastes brevispinis, code 30100)
#### 12. Shortspine thornyhead (Sebastolobus alascanus, code 30020
#### 13. Rougheye and blackspotted rockfishes (Sebastes aleutianus and  
####     Sebastes melanostictus, respectively, codes 30050,30051,30052)

#### 14/15 Northern and Southern rock sole (Lepidopsetta polyxystra and L.
#### bilineata, respectivity, codes 10260,10261,10262)
##################################################
data = subset(data,
              COMMON_NAME %in% c("Pacific ocean perch", 
                                 "arrowtooth flounder", 
                                 "Pacific cod", 
                                 
                                 "walleye pollock", 
                                 "Pacific halibut", 
                                 "rex sole", 
                                 
                                 "Dover sole",
                                 "flathead sole", 
                                 "dusky rockfish",
                                 
                                 "northern rockfish",
                                 "northern rock sole", 
                                 "southern rock sole",
                                 
                                 "B_R_rockfishes",
                                 "shortspine thornyhead",
                                 "silvergray rockfish"))

##################################################
####   Assign station depths from EFH layer
##################################################
cpue_shape = sp::SpatialPointsDataFrame(
  coords = data[, c("LONGITUDE", "LATITUDE")],
  data = data,
  proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

cpue_shape_aea <- sp::spTransform(x = cpue_shape,
                                  CRSobj = crs(bathy))
cpue_shape_aea@data$depth =  raster::extract(x = bathy,
                                             y = cpue_shape_aea,
                                             method = "simple")

##################################################
####   Plot bathymetry and station locations along with stations without 
####   assigned depths
####
##################################################
plot(bathy)
plot(cpue_shape_aea, 
     add = T,
     pch = ".")
plot(cpue_shape_aea[is.na(cpue_shape_aea@data$depth),], 
     add = T,
     pch = 16,
     col = 'red')

##################################################
####   Plot correlation between EFH depths and reported depth from BTS
##################################################
plot(depth ~ BOTTOM_DEPTH,
     data = cpue_shape_aea@data,
     subset = COMMON_NAME == "arrowtooth flounder",
     xlab = "Depth recorded by the BTS",
     ylab = "Depth extracted from EFH layer")
with(subset(cpue_shape_aea@data,
            subset = COMMON_NAME == "arrowtooth flounder"), 
     cor(depth, BOTTOM_DEPTH, use = "complete.obs"))

##################################################
####   There are 645 observations (43 stations) without depths, 
####   near the edges of the data layer. For these stations, we extract raster
####   values averaged within an iteratively increasing buffer radius. The 
####   farthest radius was 7.6 km. 
##################################################
how_many_NAs = sum(is.na(cpue_shape_aea@data$depth)) #643 stations

distance = 200 #Initial distance
while (how_many_NAs) {
  cpue_shape_aea@data$depth[is.na(cpue_shape_aea@data$depth)] <- 
    raster::extract(x = bathy,
                    y = cpue_shape_aea[is.na(cpue_shape_aea@data$depth),],
                    buffer = distance,
                    na.rm = T,
                    fun = mean)
  
  how_many_NAs = sum(is.na(cpue_shape_aea@data$depth))
  print(distance)
  print(how_many_NAs)
  distance = distance + 200 #Increase buffer by 200 m in next iteration
}

##################################################
####   Attach depths to dataset, scaled
##################################################
data$DEPTH_EFH = cpue_shape_aea@data$depth
data$LOG_DEPTH_EFH = log(cpue_shape_aea@data$depth)
data$LOG_DEPTH_EFH_CEN = scale(data$LOG_DEPTH_EFH)
data$LOG_DEPTH_EFH_CEN_SQ = data$LOG_DEPTH_EFH_CEN^2

##################################################
#### Save
##################################################
data = data[order(data$YEAR, data$SPECIES_NAME),]
write.csv(x = data, 
          file = paste0(github_wd, "GOA_multspp.csv"), 
          row.names = F)
