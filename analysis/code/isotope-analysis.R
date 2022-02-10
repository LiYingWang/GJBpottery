library(tidyverse)
isotope_data <-
  readxl::read_excel(here::here("Data", "12_21_potshard_13C.xlsx"),
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
tem <- rsvg::rsvg(here::here("carbon-isotope-ellipses-template.svg"))
png::writePNG(tem, "tem.png", dpi =300)
browseURL("tem.png") # take a look

# plot C16 and C18
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

ggsave(here::here("delta_C16_C18.png"),
       width = 8,
       height = 8,
       dpi = 300,
       units = "in")

