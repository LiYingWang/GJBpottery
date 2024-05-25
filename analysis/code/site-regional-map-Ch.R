library(ggplot2)
theme_set(theme_bw(base_size = 6))
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
#devtools::install_github('3wen/legendMap')
library(legendMap)
library(tmaptools)
library(shadowtext)
library(tidyverse)

# points for site location
site_location <-
  data.frame(location = c("Guijuabao", "Gaoshan"), # Chengdu(104.0587, 30.5899)
             lon = c(101.3617, 103.3446),
             lat = c(27.2657, 30.2709))

# get topographic map from the Stadia Maps
library(ggmap)
# register_stadiamaps("API goes here") # API
# from https://stackoverflow.com/questions/77432892/problem-visualizing-maps-with-staten-maps
# we don't want to download every time, so let's save the map locally
China_map <- get_stadiamap(rbind(as.numeric(c(80, 10, #ggmap
                                              125, 45))), zoom = 5, maptype = "stamen_terrain") #stamen_terrain_background
#saveRDS(China_map, here("analysis", "data", "derived_data", "China_map.rds"))
GJB_map <- get_stadiamap(rbind(as.numeric(c(99, 25, #ggmap
                                            106, 32))), zoom = 10, maptype = "stamen_terrain_background")
#saveRDS(GJB_map, here("analysis", "data", "derived_data", "GJB_map.rds"))
China_map <- readRDS(here::here("analysis", "data", "derived_data", "china_map.rds"))
GJB_map <- readRDS(here::here("analysis", "data", "derived_data", "GJB_map.rds"))

# site map
China_map_with_site <-
  GJB_map +
  geom_point(data = site_location,
             aes(x = lon, y = lat), size = 1, color = "red") +
  geom_shadowtext(data = site_location,  #geom_shadowtext
                  aes(x = lon, y = lat, label = location),
                  color='black', bg.colour='white', size = 2.8,
                  position = position_nudge(y = - 0.2),
                  check.overlap = TRUE) +
  annotate(geom = "text", x = 104, y = 30.7, label = "Chengdu\nPlain",
           fontface = "italic", color = "grey1", size = 2.5) +
  annotate(geom = "text", x = 104.8, y = 30, label = "Sichuan\nBasin",
           fontface = "italic", color = "grey1", size = 2.5) +
  annotate(geom = "text", x = 102.7, y = 27.5, label = "Yanyuan\nBasin",
           fontface = "italic", color = "grey1", size = 2.5) +
  annotate(geom = "text", x = 101.5, y = 29.3, label = "Hengduan\nMountains",
           fontface = "italic", color = "grey1", size = 2.5) +
  annotate("segment", x = 102.2, xend = 101.8, y = 27.5, yend = 27.5,
           color = "grey5", arrow = arrow(length = unit(.15,"cm"))) +
  coord_sf(xlim = c(99.5, 105.5), ylim = c(25.5, 31.5),
           expand = FALSE) +
  #scale_x_continuous(breaks = c(121.0, 121.5, 122.0),
  #limits = c(120.9, 122.7)) +
  legendMap::scale_bar(
    # edit these numbers to select a suitable location, not cover important details on the map
    lon = 104,
    lat = 25.7,
    legend_size = 1.7, # size of scale legend
    distance_lon = 50, # distance of one section of scale bar, in km
    distance_lat = 1.8, # height of the scale bar, in km
    distance_legend = 17, # distance between scale bar and units, in km
    dist_unit = "km", # units of scale bar
    orientation = TRUE, # add the north arrow
    arrow_length = 10, # length of N arrow, in km
    arrow_distance = 35, # distance between scale bar & base of N arrow, in km
    arrow_north_size = 2.5) + # size of letter 'N' on N arrow, in km
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

# leaflet
library(leaflet)

# Satellite: World Imagery
satellite <- leaflet() %>%
  addTiles() %>%
  setView(lng = 99.5, lat = 30.5, zoom = 5) %>%
  addProviderTiles("Esri.WorldImagery") #save this as satellite_sw_china.png

library(magick)
satellite <- image_read(here::here("analysis","figures", "satellite_sw_china.png"))
satellite_crop <- image_crop(satellite, geometry = "1000x748+84-20") # 1084,768

satellite_map <-
  ggdraw() +
  draw_image(satellite_crop)

library(cowplot)
plt1 <-
  plot_grid(satellite_map,
            China_map_with_site,
            ncol = 2)

ggsave(here::here("analysis","figures", "Sichuan-sites-map.png"),
       width = 6,
       height = 3,
       dpi = 300,
       units = "in")
