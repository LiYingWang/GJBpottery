library(ggplot2)
theme_set(theme_bw(base_size = 6))
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(here)
library(ggspatial)
#devtools::install_github('3wen/legendMap')
library(legendMap)
library(tmaptools)
library(shadowtext)
library(tidyverse)

world <- ne_countries(scale = "medium", returnclass = "sf")
# class(world)
world_points <- sf::st_point_on_surface(world)
world_points <- cbind(world, st_coordinates(st_point_on_surface(world$geometry))) %>%
  filter(brk_name %in% c("Thailand", "China", "Myanmar", "Laos",
                         "Vietnam", "Cambodia", "Bangladesh", "Bhutan")) %>%
  mutate(Y = ifelse(brk_name %in% c("Bhutan", "Bangladesh"), Y+2, Y)) %>%
  mutate(Y = ifelse(brk_name %in% c("Thailand", "Myanmar"), Y+3.5, Y))

# add site location
site_location <-
  data.frame(location = c("Guijuabao", "Gaoshan"), # Chengdu(104.0587, 30.5899)
             lon = c(101.3617, 103.3446),
             lat = c(27.2657, 30.2709))

China_SE_Asia <-
  ggplot(data = world) +
  geom_sf(fill = "antiquewhite") + #antiquewhite
  geom_rect(xmin = 99.5, xmax = 105.5, ymin = 25.5, ymax = 31.5,
            fill = NA, colour = "red", size = 0.5) +
  geom_shadowtext(data= world_points, #geom_shadowtext
                  aes(x = X, y = Y,
                      label = brk_name),
                  color='black',
                  bg.colour='white',
                  size = 2,
                  position = position_nudge(y = - 1.7, x = 0.5)) +
  #annotate(geom = "text", x = 102.3, y = 9.4, label = "Gulf of\nTailand",
           #fontface = "italic", color = "grey22", size = 2) +
  coord_sf(xlim = c(85, 120), ylim = c(18, 45), expand = FALSE) + #add datum = NA to remove
  scale_x_continuous(breaks = seq(80, 125, by = 10)) +
  scale_y_continuous(breaks = seq(5, 45, by = 10)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

# leaflet
library(leaflet)

# Background 2: World Imagery
satellite <- leaflet() %>%
  addTiles() %>%
  setView( lng = 101, lat = 33, zoom = 5) %>%
  addProviderTiles("Esri.WorldImagery")

library(png)
satellite <- readPNG(here::here("analysis","figures", "satellite_sw_china.png"),
                     native = FALSE, info = FALSE)

region_satellite_map <-
  China_map %>%
  ggmap() +
  geom_rect(xmin = 99, xmax = 105.5, ymin = 25, ymax = 31.5,
            fill = NA, colour = "red", size = 0.5) +
  coord_sf(xlim = c(88, 116.5), ylim = c(13.5, 40), expand = FALSE) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

# Topographic map, get map from the Stadia Maps
library(ggmap)
register_stadiamaps("707a2ccb-7cfa-4aca-b2e0-5c517c4b9d66") # API
# from https://stackoverflow.com/questions/77432892/problem-visualizing-maps-with-staten-maps
# we don't want to download every time, so let's save the map locally
China_map <- get_stadiamap(rbind(as.numeric(c(80, 10, #ggmap
                                              125, 45))), zoom = 5, maptype = "stamen_terrain") #stamen_terrain_background
#saveRDS(China_map, here("analysis", "data", "derived_data", "China_map.rds"))
#GJB_map <- get_stadiamap(rbind(as.numeric(c(99, 25, #ggmap
                                            #106, 32))), zoom = 10, maptype = "stamen_terrain_background")
#saveRDS(GJB_map, here("analysis", "data", "derived_data", "SW_china_map.rds"))
China_map <- readRDS(here::here("analysis", "data", "derived_data", "china_map.rds"))
SW_china_map <- readRDS(here::here("analysis", "data", "derived_data", "SW_china_map.rds"))

region_map <-
  China_map %>%
  ggmap() +
  geom_rect(xmin = 99, xmax = 105.5, ymin = 25, ymax = 31.5,
            fill = NA, colour = "red", size = 0.5) +
  coord_sf(xlim = c(88, 116.5), ylim = c(13.5, 40), expand = FALSE) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

China_map_with_site <-
  SW_china_map +
  geom_point(data = site_location,
             aes(x = lon,
                 y = lat),
             size = 1,
             color = "red") +
  geom_shadowtext(data = site_location,  #geom_shadowtext
                  aes(x = lon,
                      y = lat,
                      label = location),
                  color='black',
                  bg.colour='white',
                  size = 2.8,
                  position = position_nudge(y = - 0.2),
                  check.overlap = TRUE) +
  annotate(geom = "text", x = 104, y = 30.7, label = "Chengdu\nPlain",
           fontface = "italic", color = "grey1", size = 2.5) +
  annotate(geom = "text", x = 104.8, y = 30, label = "Sichuan\nBasin",
           fontface = "italic", color = "grey1", size = 2.5) +
  annotate(geom = "text", x = 102.7, y = 27.5, label = "Yanyuan\nBasin",
           fontface = "italic", color = "grey1", size = 2.5) +
  annotate(geom = "text", x = 101.5, y = 29.3, label = "Hengduan\nMountain",
           fontface = "italic", color = "grey1", size = 2.5) +
  annotate("segment", x = 102.2, xend = 101.8, y = 27.5, yend = 27.5,
           color = "grey5", arrow = arrow(length = unit(.15,"cm"))) +
  coord_sf(xlim = c(99.5, 105.5),
           ylim = c(25.5, 31.5),
           expand = FALSE) +
  #scale_x_continuous(breaks = c(121.0, 121.5, 122.0),
  #limits = c(120.9, 122.7)) +
  legendMap::scale_bar(
    # edit these numbers to select a suitable location
    # for the scale bar where it does not cover
    # important details on the map
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
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

ggplot() +
  scale_x_continuous(limits = c(0, 2.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  coord_equal() +
  annotation_custom(ggplotGrob(region_map), #China_SE_Asia, #region_map
                    xmin = 0,
                    xmax = 1.5,
                    ymin = 0,
                    ymax = 1) +
  annotation_custom(ggplotGrob(China_map_with_site),
                    xmin = 1.2,
                    xmax = 2.5,
                    ymin = 0,
                    ymax = 1) +
  theme_void() +
  theme(plot.background = element_rect(color = "white", fill = "white"))


ggsave(here::here("analysis","figures", "Sichuan-sites-map2.png"),
       width = 6,
       height = 3,
       dpi = 300,
       units = "in")

