library(ggplot2)
theme_set(theme_bw(base_size = 6))
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(here)
library(ggspatial)
# devtools::install_github('3wen/legendMap')
library(legendMap)
library(tmaptools)
library(shadowtext)
library(tidyverse)

world <- ne_countries(scale = "medium", returnclass = "sf")
# class(world)
world_points <- sf::st_point_on_surface(world)
world_points <- cbind(world, st_coordinates(st_point_on_surface(world$geometry))) %>%
  filter(brk_name %in% c("Thailand", "China", "Myanmar",
                         "Laos", "Vietnam", "Cambodia")) %>%
  mutate(Y = Y+2, X = X-1)

# add site location
site_location <-
  data.frame(location = c("Guijuabao", "Gaoshan"), # Chengdu(104.0587, 30.5899)
             lon = c(101.3617, 103.3446),
             lat = c(27.2657, 30.2709))

China_SE_Asia <-
  ggplot(data = world) +
  geom_sf( fill= "antiquewhite") +
  geom_rect(xmin = 100, xmax = 106, ymin = 26, ymax = 32,
            fill = NA, colour = "red", size = 0.5) +
  geom_shadowtext(data= world_points,
                  aes(x = X, y = Y,
                      label = brk_name),
                  color='black',
                  bg.colour='white',
                  size = 2,
                  position = position_nudge(y = - 1.7, x = 0.5)) +
  #annotate(geom = "text", x = 102.3, y = 9.4, label = "Gulf of\nTailand",
           #fontface = "italic", color = "grey22", size = 2) +
  coord_sf(xlim = c(80, 123), ylim = c(9, 45), expand = FALSE) + #add datum = NA to remove
  scale_x_continuous(breaks = seq(80, 125, by = 10)) +
  scale_y_continuous(breaks = seq(5, 45, by = 10)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

# Topographic map
library(ggmap)
# we don't want to download every time, so let's save the map locally
# from https://stackoverflow.com/a/52710855/1036500
GJB_map <- ggmap(get_stamenmap(rbind(as.numeric(c(100, 26,
                                                  106, 32))), zoom = 10))
# saveRDS(tw_map, here("analysis", "data", "raw_data", "tw_map.rds"))
# tw_map <- readRDS(here("analysis", "data", "raw_data", "tw_map.rds"))
pg <- ggplot_build(GJB_map)

China_map_with_site <-
 GJB_map +
  geom_point(data = site_location,
             aes(x = lon,
                 y = lat),
             size = 1,
             color = "red") +
  geom_shadowtext(data = site_location,
                  aes(x = lon,
                      y = lat,
                      label = location),
                  color='black',
                  bg.colour='white',
                  size = 2.6,
                  position = position_nudge(y = - 0.2),
                  check.overlap = TRUE) +
  annotate(geom = "text", x = 104, y = 30.7, label = "Chengdu\nPlain",
           fontface = "italic", color = "grey15", size = 2) +
  coord_sf(xlim = c(100, 106),
           ylim = c(26, 32),
           expand = FALSE) +
  #scale_x_continuous(breaks = c(121.0, 121.5, 122.0),
  #limits = c(120.9, 122.7)) +
  legendMap::scale_bar(
    # edit these numbers to select a suitable location
    # for the scale bar where it does not cover
    # important details on the map
    lon = 104.6,
    lat = 26.2,
    legend_size = 1.6, # size of scale legend
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
  annotation_custom(ggplotGrob(China_SE_Asia),
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

ggsave(here::here("analysis","figures", "Sichuan-sites-map.png"),
       width = 6,
       height = 3,
       dpi = 300,
       units = "in")

