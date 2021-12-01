library(tidyverse)
library(tidycensus)
library(sf)
library(lubridate)
library(kableExtra)
library(patchwork)
library(ggplot2)
library(viridis)
library(gridExtra)
library(knitr)
library(dplyr)
library(osmdata)
library(tigris)
library(riem)

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
  left_join(parking_census, by = c("street_seg_ctrln_id"="CNN"))%>%
  drop_na(PRKG_SPLY)

off_parking_census <- read.csv("/Users/inordia/Desktop/UPenn搞起来/592/pro-forma/data/SFMTA_Managed_Off-street_Parking.csv")%>%
  select(STREET_SEG_CTRLN_ID, CAPACITY)


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
      nn_function(st_c(parking_seg_census), st_c(st_centroid(cultural_district)),1))

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

theatre <- add_osm_feature(opq = q0, key = 'amenity', value = "theatre") %>%
  osmdata_sf(.)
theatre.sf <- st_geometry(theatre$osm_points) %>%
  st_transform(4326) %>%
  st_sf() %>%
  cbind(., theatre$osm_points$amenity) %>%
  rename(NAME = theatre.osm_points.amenity)%>%
  st_transform('EPSG:7131')%>%
  st_intersection(sf_boundary)%>%
  dplyr::select(geometry)%>%
  distinct()

cafe <- add_osm_feature(opq = q0, key = 'amenity', value = "cafe") %>%
  osmdata_sf(.)
cafe.sf <- st_geometry(cafe$osm_points) %>%
  st_transform(4326) %>%
  st_sf() %>%
  cbind(., cafe$osm_points$amenity) %>%
  rename(NAME = cafe.osm_points.amenity)%>%
  st_transform('EPSG:7131')%>%
  st_intersection(sf_boundary)%>%
  dplyr::select(geometry)%>%
  distinct()

clothes <- add_osm_feature(opq = q0, key = 'shop', value = "clothes") %>%
  osmdata_sf(.)
clothes.sf <- st_geometry(clothes$osm_points) %>%
  st_transform(4326) %>%
  st_sf() %>%
  cbind(., clothes$osm_points$shop) %>%
  rename(NAME = clothes.osm_points.shop)%>%
  st_transform('EPSG:7131')%>%
  st_intersection(sf_boundary)%>%
  dplyr::select(geometry)%>%
  distinct()

parking_seg_census <-
  parking_seg_census %>%
  mutate(
    restaurant.nn =
      nn_function(st_c(parking_seg_census), st_c(restaurant.sf),3),
    cafe.nn =
      nn_function(st_c(parking_seg_census), st_c(cafe.sf),3),
    theatre.nn =
      nn_function(st_c(parking_seg_census), st_c(theatre.sf),3),
    clothes.nn =
      nn_function(st_c(parking_seg_census), st_c(clothes.sf),3))

## Weather

weather.SF <- 
  riem_measures(station = "SFO", date_start = "2021-09-01", date_end = "2021-10-01")

weather.Panel <-  
  weather.SF %>%
  mutate_if(is.character, list(~replace(as.character(.), is.na(.), "0"))) %>% 
  replace(is.na(.), 0) %>%
  mutate(interval60 = ymd_h(substr(valid, 1, 13))) %>%
  mutate(week = week(interval60),
         dotw = wday(interval60, label=TRUE)) %>%
  group_by(interval60) %>%
  summarize(Temperature = max(tmpf),
            Percipitation = sum(p01i),
            Wind_Speed = max(sknt)) %>%
  mutate(Temperature = ifelse(Temperature == 0, 42, Temperature))

grid.arrange(top = "Weather Data - San Francisco - Sptember, 2021",
             ggplot(weather.Panel, aes(interval60,Percipitation)) + geom_line() + 
               labs(title="Percipitation", x="Hour", y="Percipitation") + plotTheme(),
             ggplot(weather.Panel, aes(interval60,Wind_Speed)) + geom_line() + 
               labs(title="Wind Speed", x="Hour", y="Wind Speed") + plotTheme(),
             ggplot(weather.Panel, aes(interval60,Temperature)) + geom_line() + 
               labs(title="Temperature", x="Hour", y="Temperature") + plotTheme())







