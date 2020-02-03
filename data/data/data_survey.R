## Get AK groundfish bottom trawl survey data for 3 primary surveys
# Result is cleaned cpue (kg/km^2) by haul, with zeros included

# Note: EBS includes all species, whereas GOA and AI are only a subset
# here, we will filter the data to only include species represented in all regions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/data/")

library(dplyr)
library(geosphere)
library(marmap)
library(RANN)
library(rgdal)

# or from flat files exported from AFSC database
GOA = read.csv("C:/Users/zack.oyafuso/Desktop/data-raw/cpue_GOA_selected_spp.csv", 
               stringsAsFactors = FALSE) # CPUE is (num or kg / km^2)

# Filter Species: 
# Arrowtooth Flounder (Atherestes stomias, code 10110)
# Pacific Cod (Gadus macrocephalus, code 21720)
# Pacific Ocean Perch (Sebastes alutus, code 30060)
# Sablefish (Anoplopoma fimbria, code 20510)
# Walleye pollock (Gadus chalcogrammus, code 21740)
# Dover sole (Solea solea, code 10180)
# Pacific halibut (Hippoglossus stenolepis, code 10120)
# Flathead sole (Hippoglossoides elassodon, code 10130)
# Rex sole (Glyptocephalus zachirus, code 10200)
# Dusky rockfish (Sebastes variabilis, code 30152)
# Northern rockfish (Sebastes polyspinis, code 30420)

# Rougheye and blackspotted rockfishes (Sebastes aleutianus and Sebastes melanostictus, respectively, codes 30050,30051,30052)
# Northern and Southern rock sole (Lepidopsetta polyxystra and Lepidopseta bilineata, respectivity, codes 10260,10261,10262)

data <- filter(GOA, SPECIES_CODE %in% c(10110, 21720, 30060, 20510, 21740, 
                                        10180, 10120, 10130, 10200, 30152,
                                        30420, 30050,30051,30052, ))

#Add species names
data$SCI = NA
for(i in 1:nrow(data)){
  data$SCI[i] = switch(paste(data$SPECIES_CODE[i]), 
                       '10110' = 'Atherestes stomias',  
                       '21720' = 'Gadus macrocephalus', 
                       '30060' = 'Sebastes alutus', 
                       '20510' = 'Anoplopoma fimbria', 
                       '21740' = 'Gadus chalcogrammus', 
                       '10180' = 'Solea solea', 
                       '10120' = 'Hippoglossus stenolepis', 
                       '10130' = 'Hippoglossoides elassodon', 
                       '10200' = 'Glyptocephalus zachirus', 
                       '30152' = 'Sebastes variabilis'
  )
}


# join haul data to get coordinates
haul <- read.csv("data-raw/haul.csv", stringsAsFactors = FALSE)
haul <- cbind(haul, 
              geosphere::midPoint(cbind(haul$START_LONGITUDE, haul$START_LATITUDE), 
                                  cbind(haul$END_LONGITUDE, haul$END_LATITUDE))) # get haul midpoints
haul <- haul %>% select(HAULJOIN, LATITUDE = lat, LONGITUDE = lon)
data <- inner_join(data, haul)

# get biomass from cpue and effort, remove extraneous variables, remove NAs
data <- data %>% mutate(BIOMASS = WGTCPUE * EFFORT) %>%
  select(HAULJOIN, YEAR, BIOMASS, EFFORT, LATITUDE, LONGITUDE, SCI)%>%
  tidyr::drop_na()

data = data[order(data$YEAR), ]

#####################################
## Bathymetry
#######################################
xmin <- 180
xmax <- 240
ymin <- 50
ymax <- 62
bathy <- getNOAA.bathy(lon1 = xmin-360, lon2 = xmax-360, 
                       lat1 = ymin, lat2 = ymax, resolution = 1)

bathy =  fortify.bathy(bathy)

#####################################
## Attach a bathmetry value to the survey points
#####################################
cord.dec = SpatialPoints(cbind(data$LONGITUDE, data$LATITUDE), 
                         proj4string = CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+proj=utm +zone=5 +datum=WGS84"))
data[,c('X', 'Y')] = cord.UTM@coords

cord.dec = SpatialPoints(cbind(bathy$x, bathy$y), 
                         proj4string = CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+proj=utm +zone=5 +datum=WGS84"))
bathy[,c('X', 'Y')] = cord.UTM@coords

nearest = nn2(query = data[,c('X', 'Y')], 
              data = bathy[,c('X', 'Y')],
              k = 150)

data$DEPTH = bathy$z[nearest$nn.idx[,1]]

#If the closest bathymetry point to a survey station is land (depth > 0), 
#find the next nearest neighbor and assign that bathymetry point. 
#Continue until all survey points are associated with a negative bathy value. 

inn = 2
while(sum(data$DEPTH > 0) !=0 ){
  for(irow in which(data$DEPTH > 0)){
    data$DEPTH[irow] = bathy$z[nearest$nn.idx[irow,inn]]
  }
  inn = inn + 1
}

#Convert bathymetry to positive meters, center, and provide the square of depth
data$DEPTH = -data$DEPTH 
data$DEPTH = scale(x = data$DEPTH)
data$DEPTH2 = data$DEPTH^2

#Save Data Objects
write.csv(x = bathy, file = 'bathymetry/bathy.csv', row.names = F)
write.csv(x = data, file = "data/GOA_multspp.csv", row.names = F)
