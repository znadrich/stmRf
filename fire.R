library(sf)
library(ggplot2)
library(lubridate)
library(dplyr)
library(gganimate)
library(sf)
library(maps)
library(maptools)

latlong2state <- function(pointsDF) {
  # Prepare SpatialPolygons object with one SpatialPolygon
  # per state (plus DC, minus HI & AK)
  states <- map('state', fill=TRUE, col="transparent", plot=FALSE)
  IDs <- sapply(strsplit(states$names, ":"), function(x) x[1])
  states_sp <- map2SpatialPolygons(states, IDs=IDs,
                                   proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  # Convert pointsDF to a SpatialPoints object 
  pointsSP <- SpatialPoints(pointsDF, 
                            proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  # Use 'over' to get _indices_ of the Polygons object containing each point 
  indices <- over(pointsSP, states_sp)
  
  # Return the state names of the Polygons object containing each point
  stateNames <- sapply(states_sp@polygons, function(x) x@ID)
  stateNames[indices]
}

#https://fsapps.nwcg.gov/afm/gisdata.php
setwd('C:/Users/Zack/Downloads/viirs_iband_fire_2018_365_conus_shapefile/')
dbf <- st_read('viirs_iband_fire_2018_365_conus.dbf')
df <- as.data.frame(dbf)
df$latitude <- df$LAT
df$longitude <- df$LONG


# Test the function using points in Wisconsin and Oregon.
testPoints <- data.frame(x = df$LONG, y = df$LAT)

df$state <- latlong2state(testPoints)
df <- df[!is.na(df$state), ]
df <- df[df$state == 'california', ]

df[df$DATE == '2018-09-20', ] %>%
  ggplot(aes(x=LONG, y=LAT)) +
  geom_point()

df %>%
  mutate(mnth = month(DATE)) %>%
  group_by(mnth) %>%
  summarize(n=n())

df[month(df$DATE) == 9, ] %>%
  mutate(mnth = month(DATE),
         begin = day(DATE),
         length = as.integer(1),
         exit = as.integer(1)) %>%
  ggplot(aes(x=LONG, y=LAT)) +
  geom_point() +
  transition_events(start = begin,
                    end = begin + length,
                    enter_length = as.integer(0),
                    exit_length = as.integer(0)) +
  ggtitle('{frame_time}')

df %>%
  group_by(JULIAN) %>%
  summarize(n=n()) %>%
  ggplot(aes(x=JULIAN, y=n)) +
  geom_line()
