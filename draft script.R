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
  left_join(parking_census, by = c("street_seg_ctrln_id"="CNN"))%>%
  drop_na(PRKG_SPLY)

off_parking_census <- read.csv("data/SFMTA_Managed_Off-street_Parking.csv")%>%
  select(STREET_SEG_CTRLN_ID, CAPACITY)


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
      nn_function(st_c(parking_seg_census), st_c(st_centroid(cultural_district)),1))

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





# add up the time of parking in each hour
for (row in 1:nrow(meter_trans)) {
  start.interval = meter_trans[row, "interval60"]
  start.time = meter_trans[row, "session_start_dt"]
  duration  = meter_trans[row, "parking_time"]
  post_id_row = meter_trans[row, "post_id"]
  row_number = parking.time.panel%>%
    filter(post_id==post_id_row,
           interval60==start.interval)%>%
    pull(index)
  time.first.hour = as.numeric(ymd_hms(start.interval)+3600-ymd_hms(start.time),
                               units="secs")
  duration=duration - time.first.hour
  parking.time.panel[row_number,"parking_time_in_60m"] = parking.time.panel[row_number,"parking_time_in_60m"] + time.first.hour
  
  while(duration>0){
    row_number=row_number+1
    time.this.hour = min(duration, 3600)
    parking.time.panel[row_number,"parking_time_in_60m"] = parking.time.panel[row_number,"parking_time_in_60m"] + time.this.hour
    duration = duration - time.this.hour
  }
}


meters.bySegId = meters%>%st_drop_geometry()%>%
  group_by(street_seg_ctrln_id)%>%
  summarise(n = n())%>%
  # filter(on_offstreet_type=="ON")%>%
  left_join(parking_seg_census,by="street_seg_ctrln_id")%>%
  mutate(sub = n-PRKG_SPLY)

meters.bySegId%>%
  filter(abs(sub)<=2)

meters.bySegId%>%
  ggplot()+
  geom_histogram(aes(sub),bins=50)

meters.bySegId%>%
  ggplot()+
  geom_point(aes(ms_space_num,PRKG_SPLY))+
  geom_abline(intercept = 0, slope = 1,color="red")

study.panel =
  expand.grid(interval60 = unique(weather.Panel$interval60), 
              start.station.id = unique(parking_seg_census$street_seg_ctrln_id)) 

plotData.lag <-
  as.data.frame(park.engineer)%>%
  dplyr::select(starts_with("lag"), parking_time_in_60m) %>%
  gather(Variable, Value, -parking_time_in_60m) %>%
  mutate(Variable = fct_relevel(Variable, "lagHour","lag2Hours","lag3Hours",
                                "lag4Hours","lag12Hours","lag1day"))
correlation.lag <-
  group_by(plotData.lag, Variable) %>%
  summarize(correlation = round(cor(Value, parking_time_in_60m, use = "complete.obs"), 2)) %>%
  kable(caption = "Parking Time in One Hour") %>%
  kable_styling("striped", full_width = F)

correlation.lag

park.engineer <-
  park.engineer%>%
  mutate(week = week(interval60),
         dotw = wday(interval60, label=TRUE))
sundays <- 
  mutate(park.engineer,
         sunday = ifelse(dotw == "Sun" & hour(interval60) == 1,
                         interval60, 0)) %>%
  filter(sunday != 0) 

park.engineer.test <-
  park.engineer%>%
  st_as_sf()

ggplot()+
  geom_sf(data =sf_boundary)+
  geom_sf(data = st_as_sf(park.engineer) %>%
               group_by(street.id, week)%>%
               tally(),
             aes(color = n), 
             fill = "transparent", alpha = 0.8, size = 1)+
  scale_colour_viridis(direction = -1,
                       discrete = FALSE, option = "D")+
  ylim(min(dat_census$start_lat), max(dat_census$start_lat))+
  xlim(min(dat_census$start_lon), max(dat_census$start_lon))+
  facet_grid(~week)+
  labs(title="Sum of rideshare trips by station and week",
       subtitle = "Philadelphia, October 8th - November 11st, 2019")+
  mapTheme()



library(hms)

int.ampeak <- interval(as_hms("07:00:00"),as_hms("10:00:00"))
int.pmpeak <- interval(as_hms("15:00:00"),as_hms("18:00:00"))
int.midday <- interval(as_hms("10:00:00"),as_hms("15:00:00"))
int.night <- interval(as_hms("18:00:00"),as_hms("23:00:00"))
int.overnight <- interval(as_hms("00:00:00"),as_hms("07:00:00"))

park.engineer <-park.engineer%>%
  mutate(period = case_when(ymd_hms(paste0('1970-01-01',str_sub(interval60, 12)))%within% int.ampeak ~ "AM",
                            ymd_hms(paste0('1970-01-01',str_sub(interval60, 12)))%within% int.pmpeak ~ "PM",
                            ymd_hms(paste0('1970-01-01',str_sub(interval60, 12)))%within% int.midday ~ "Mid_Day",
                            ymd_hms(paste0('1970-01-01',str_sub(interval60, 12)))%within% int.night ~ "Night",
                            ymd_hms(paste0('1970-01-01',str_sub(interval60, 12)))%within% int.overnight ~ "Overnight"))


parking.time.panel.sum %>%
  group_by(interval60, street_seg_ctrln_id, period) %>%
  summarize(mean_time = mean(parking_time_in_60m))%>%
  ggplot()+
  geom_histogram(aes(mean_time), binwidth = 10000)+
  labs(title="Mean Number of Hourly Trips Per Station",
       subtitle="Philadelphia, October 8th - November 11st, 2019",
       x="Mean Parking Time", 
       y="Frequency")+
  facet_wrap(~period)+
  plotTheme()

weather.Panel$interval60<- as.character(weather.Panel$interval60)

park.engineer <- park.engineer%>%
  left_join(weather.Panel, by = "interval60")

sf_neighborhood <-st_read("/Users/inordia/Desktop/UPenn搞起来/592/pro-forma/data/SF Find Neighborhoods.geojson")%>%
  st_transform('EPSG:7131')
  

park.engineer %>%
  group_by(interval60, street.id, period) %>%
  summarize(mean_time = mean(parking_time_in_60m))%>%
  ungroup()%>%
  left_join(park.engineer)%>%
  st_as_sf()%>%
ggplot()+
  geom_sf(aes(color=mean_time),
          fill = "transparent", alpha=0.8, size=1)+
  geom_sf(data = sf_neighborhood)+
  facet_grid(~period)

  geom_point(data = dat_census %>%
               group_by(start_station, start_lat, start_lon, week)%>%
               tally(),
             aes(x=start_lon, y = start_lat, color = n), 
             fill = "transparent", alpha = 0.8, size = 1)+
  scale_colour_viridis(direction = -1,
                       discrete = FALSE, option = "D")+
  ylim(min(dat_census$start_lat), max(dat_census$start_lat))+
  xlim(min(dat_census$start_lon), max(dat_census$start_lon))+
  facet_grid(~week)+
  labs(title="Sum of rideshare trips by station and week",
       subtitle = "Philadelphia, October 8th - November 11st, 2019")+
  mapTheme()

parking.time.panel.sum<-
  parking.time.panel.sum%>%
mutate(week = week(interval60),
       dotw = wday(interval60, label=TRUE))



park.engineer %>%
  group_by(interval60, street.id, period) %>%
  summarize(mean_time = mean(parking_time_in_60m))%>%
  ungroup()%>%
  left_join(park.engineer)%>%
  st_as_sf()%>%
  ggplot()+
  geom_sf(aes(color=mean_time),
          fill = "transparent", alpha=0.8, size=0.2)+
  geom_sf(data = sf_neighborhood)+
  facet_grid(~period)  
