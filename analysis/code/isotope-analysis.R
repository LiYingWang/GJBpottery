library(tidyverse)
isotope_data <-
  readxl::read_excel(here::here("analysis","data", "raw_data","12_21_potshard_13C.xlsx"),
                     sheet = "Summary", col_names =  FALSE)

# filter the blank corrected data
delta_blank_cor <-
  isotope_data %>%
  select(1, 15, 17) %>%
  slice(2:8) %>%
  dplyr::mutate(...1 = ifelse(is.na(...1), paste0("C13"), ...1))

names(delta_blank_cor) <- delta_blank_cor[1,]

delta_blank_cor <-
  delta_blank_cor %>%
  slice(2:7) %>%
  dplyr::mutate(across(starts_with("13C"), as.numeric))

# filter the Meth corrected data
delta_meth_cor <-
  isotope_data %>%
  select(1, 7, 9) %>%
  slice(2:8) %>%
  dplyr::mutate(...1 = ifelse(is.na(...1), paste0("C13"), ...1))

names(delta_meth_cor) <- delta_meth_cor[1,]

delta_meth_cor <-
  delta_meth_cor %>%
  slice(2:7) %>%
  dplyr::mutate(across(starts_with("13C"), as.numeric))

# import a SVG for later ggplot overlaying with ellipses
tem <- rsvg::rsvg(here::here("analysis", "figures","carbon-isotope-ellipses-template.svg"))
tem2 <- rsvg::rsvg(here::here("analysis", "figures","C16-18-ellipse.svg"))
png::writePNG(tem2, "tem2.png", dpi =300)
browseURL("tem2.png") # take a look

# plot Carbon isotopes of 16 and C18
ggplot(delta_blank_cor, # delta_meth_cor or delta_blank_cor
       aes(`13C C16:0`,`13C C18:0`)) +
  geom_point(size = 1, alpha = 0.9) +
  ggrepel::geom_text_repel(aes(label = C13)) +
  theme_minimal(base_size = 14) +
  labs(x = bquote(delta*{}^13*"C 16:0 \u2030"),
       y = bquote(delta*{}^13*"C 18:0 \u2030")) +
  xlim(-40,-10) +
  ylim(-40,-10) +
  coord_fixed(ratio=1) +
  #scale_colour_viridis_d(direction = -1) +
  annotation_raster(tem, ymin = -42.8, ymax= -8 ,
                    xmin = -44.5 , xmax = -5.65)

# plot Carbon isotopes of 16 and C18 using a different reference figure
delta_blank_cor %>% # delta_meth_cor or delta_blank_cor
  filter(!C13 == "Pot 6") %>% # remove pot from CDG
  ggplot(aes(`13C C16:0`,`13C C18:0`)) +
  geom_point(size = 2, alpha = 0.9, color = "red") +
  ggrepel::geom_text_repel(aes(label = C13), size = 5) +
  theme_minimal(base_size = 14) +
  labs(x = bquote(delta*{}^13*"C 16:0 \u2030"),
       y = bquote(delta*{}^13*"C 18:0 \u2030")) +
  xlim(-40,-20) +
  ylim(-40,-20) +
  coord_fixed(ratio=1) +
  #scale_colour_viridis_d(direction = -1) +
  annotation_raster(tem2, ymin = -45.1, ymax= -15.1 ,
                    xmin = -45.55, xmax = -15.55)

ggsave(here::here("delta_C16_C18_remove_pot6.png"),
       width = 8,
       height = 8,
       dpi = 300,
       units = "in")
