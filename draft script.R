library(tidyverse)
library(tidycensus)
library(sf)
library(kableExtra)
library(patchwork)
library(ggplot2)
library(viridis)
library(gridExtra)
library(knitr)
library(dplyr)
library(osmdata)
library(tigris)

st_c <- st_coordinates


options(scipen=999)
options(tigris_class = "sf")
options(tigris_use_cache = TRUE)
source("functions.r")

p5 = inferno(5)
p2 = p5[c(4,2)]
crs = 7131
v = mapview::mapview

mt = read.csv("data/meter_trans_2021-09.csv")
mt%>%nrow()

meters = st_read("data/Parking Meters.geojson")

meters = meters%>%
  select(post_id, ms_space_num, on_offstreet_type, active_meter_flag, meter_type,
         street_id, street_seg_ctrln_id, longitude, latitude, geometry)%>%
  st_transform(crs)%>%
  mutate(
    street_seg_ctrln_id = street_seg_ctrln_id%>%as.numeric()%>%as.character(),
    ms_space_num = as.numeric(ms_space_num),
    ms_space_num = if_else(ms_space_num==0,1,ms_space_num)
  )

## SF Open Data

parking_census <- read.csv("data/On-Street_Parking_Census.csv")%>%
  select(CNN, PRKG_SPLY)

parking_census$CNN <- as.character(parking_census$CNN)

parking_seg_census <- meters%>%
  select(street_seg_ctrln_id, geometry)%>%
  distinct(street_seg_ctrln_id, .keep_all = T)

parking_seg_census <- parking_seg_census%>%
  left_join(parking_census, by = c("street_seg_ctrln_id"="CNN"))


transit <- st_read("data/Muni Stops.geojson")%>%
  st_transform(crs)

park_public <- st_read("data/Recreation and Parks Properties.geojson")%>%
  st_transform(crs)%>%
  select(propertytype, geometry)%>%
  rename(type=propertytype)

park_private <- st_read("data/Privately Owned Public Open Spaces.geojson")%>%
  st_transform(crs)%>%
  select(type, geometry)

cultural_district <- st_read("data/Cultural Districts.geojson")%>%
  st_transform(crs)

parking_seg_census <-
  parking_seg_census %>%
  mutate(
    transit.nn =
      nn_function(st_c(parking_seg_census), st_c(transit),1),
    park_public.nn =
      nn_function(st_c(parking_seg_census), st_c(st_centroid(park_public)),1),
    cultural_district.nn =
      nn_function(st_c(parking_seg_census), st_c(st_centroid(cultural_district)),3))

## OSM Data

sf_boundary <- st_read("data/Bay Area Counties.geojson")%>%
  filter(county == "San Francisco")%>%
  st_transform(crs)%>%
  select(geometry)

q0 <- opq(bbox = c(-122.3505,37.7025,-122.5171,37.8364)) 

restaurant <- add_osm_feature(opq = q0, key = 'amenity', value = "restaurant") %>%
  osmdata_sf(.)

restaurant.sf <- st_geometry(restaurant$osm_points) %>%
  st_transform(4326) %>%
  st_sf() %>%
  cbind(., restaurant$osm_points$amenity) %>%
  rename(NAME = restaurant.osm_points.amenity)%>%
  st_transform(crs)%>%
  st_intersection(sf_boundary)%>%
  dplyr::select(geometry)%>%
  distinct()

