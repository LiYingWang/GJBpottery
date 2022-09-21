library(tidyverse)

soil <-
  readxl::read_excel(here::here("analysis","data", "raw_data","All_test_data_soil_0728_FAs.xlsx"))

# tidy data
ixd <- which(str_detect(soil$`Peak #`, "D"))
soils <- split(soil, cumsum(1:nrow(soil) %in% ixd))
names(soils) <- map(soils, ~.x$`Peak #`[2])
soils <- map_df(soils, ~.x[-(1:2), ] %>%
                  mutate_all(as.numeric), .id = 'sample')

soil_test <-
  soils %>%
  mutate(retention = round(`Ret Time`, 1)) %>%
  group_by(sample, retention) %>%
  distinct(retention, .keep_all = TRUE) %>%
  mutate(FA = case_when(
    retention %in% c(25.1, 25.2) ~ "C16:0",
    retention %in% c(29.6, 29.7) ~ "C18:0",
    retention == 33.9 ~ "C20:0",
    retention %in% c(36.0, 36.1) ~ "C21:0",
    retention == 38.0 ~ "C22:0",
    retention == 41.8 ~ "C24:0"))

# plot all
soil_test %>%
  filter(!sample == "C4_C24_stds") %>%
  filter(Area > 0) %>%
  filter(between(`Ret Time`, 5, 43)) %>%
  ggplot(aes(`Ret Time`, Area)) +
  geom_segment(aes(xend = `Ret Time`, yend = 0), size = 0.8, lineend = "butt") +
  geom_text(aes(label = FA), #ggrepel::geom_text_repel
            size = 2.5,
            nudge_x = 0,
            nudge_y = 500000,
            show.legend = FALSE) +
  facet_wrap(~sample, ncol = 3) +
  #scale_y_continuous(labels = scales::comma_format(),
                     #limits = c(0, 7500000),
                     #breaks = seq(0, 7500000, 1000000)) +
  scale_x_continuous(limits = c(20, 40),
                     expand = c(0, 0.5)) + # don't log, many peaks and distort real counts
  labs(x = "Retention Time", y = "Intensity")

ggsave(here::here("analysis","figures", "GC_MS_soil_test.png"),
       width = 6,
       height = 6,
       dpi = 360,
       units = "in")
