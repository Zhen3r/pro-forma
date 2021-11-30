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

v = mapview::mapview

setwd("/Users/inordia/Desktop/UPenn搞起来/592/pro-forma")

mt = read.csv("/Users/inordia/Desktop/UPenn搞起来/592/pro-forma/data/meter_trans_2021-09.csv")
mt%>%nrow()

meters = st_read("/Users/inordia/Desktop/UPenn搞起来/592/pro-forma/data/Parking Meters.geojson")

meters <-
  meters%>%
  select(post_id, ms_space_num, on_offstreet_type, active_meter_flag, meter_type,
         street_id, street_seg_ctrln_id, longitude, latitude, geometry)%>%
  st_transform('EPSG:7131')

meters$street_seg_ctrln_id<-as.numeric(meters$street_seg_ctrln_id)
round(meters$street_seg_ctrln_id)
meters$street_seg_ctrln_id<-as.character(meters$street_seg_ctrln_id)
meters$ms_space_num<-as.numeric(meters$ms_space_num)
round(meters$ms_space_num)
meters$ms_space_num[meters$ms_space_num==0]<-1

## SF Open Data

parking_census <- read.csv("/Users/inordia/Desktop/UPenn搞起来/592/pro-forma/data/On-Street_Parking_Census.csv")%>%
  select(CNN, PRKG_SPLY)

parking_census$CNN <- as.character(parking_census$CNN)

parking_seg_census <- meters%>%
  select(street_seg_ctrln_id, geometry)%>%
  distinct(street_seg_ctrln_id, .keep_all = T)

parking_seg_census <- parking_seg_census%>%
  left_join(parking_census, by = c("street_seg_ctrln_id"="CNN"))


transit <- st_read("/Users/inordia/Desktop/UPenn搞起来/592/pro-forma/data/Muni Stops.geojson")%>%
  st_transform('EPSG:7131')

park_public <- st_read("/Users/inordia/Desktop/UPenn搞起来/592/pro-forma/data/Recreation and Parks Properties.geojson")%>%
  st_transform('EPSG:7131')%>%
  select(propertytype, geometry)%>%
  rename(type=propertytype)

park_private <- st_read("/Users/inordia/Desktop/UPenn搞起来/592/pro-forma/data/Privately Owned Public Open Spaces.geojson")%>%
  st_transform('EPSG:7131')%>%
  select(type, geometry)

cultural_district <- st_read("/Users/inordia/Desktop/UPenn搞起来/592/pro-forma/data/Cultural Districts.geojson")%>%
  st_transform('EPSG:7131')

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

sf_boundary <- st_read("/Users/inordia/Desktop/UPenn搞起来/592/pro-forma/data/Bay Area Counties.geojson")%>%
  filter(county == "San Francisco")%>%
  st_transform('EPSG:7131')%>%
  select(geometry)

q0 <- opq(bbox = c(-122.3505,37.7025,-122.5171,37.8364)) 
restaurant <- add_osm_feature(opq = q0, key = 'amenity', value = "restaurant") %>%
  osmdata_sf(.)

restaurant.sf <- st_geometry(restaurant$osm_points) %>%
  st_transform(4326) %>%
  st_sf() %>%
  cbind(., restaurant$osm_points$amenity) %>%
  rename(NAME = restaurant.osm_points.amenity)%>%
  st_transform('EPSG:7131')%>%
  st_intersection(sf_boundary)%>%
  dplyr::select(geometry)%>%
  distinct()

