library(tidyverse)

# GC-FID results after methyl
meth <- read_csv(here::here("analysis/data/raw_data/LY_21_08_20_pot_chrom.CSV"))

# tidy meth
ixd <- which(str_detect(meth$Path, "C"))
meths <- split(meth, cumsum(1:nrow(meth) %in% ixd))
meth_data <- meths[c(FALSE, TRUE)]
meth_name <- meths[c(TRUE, FALSE)]
names(meth_data) <- map(meth_data, ~.x$Path[1])
names(meth_name) <-map(meth_name, ~.x$Sample[1])
meths <- map_df(meth_data, ~.x[-1, ] %>%
                   mutate_all(as.numeric), .id = 'sample')

# join two datasets and get sample number
meth_all <-
  meths %>%
  mutate(sample = paste0(parse_number(sample),
                         str_extract(sample, pattern = "_00[0-9]"))) %>%
  mutate(name = case_when(
    sample %in% c("210729_001") ~ "CDG-062",
    sample %in% c("210729_002") ~ "SYG-3",
    sample %in% c("210729_003") ~ "SYG-TN13-E22-2#1",
    sample %in% c("210729_004") ~ "SYG-TN13-E22-1",
    sample %in% c("210729_005") ~ "SYG-TN13-E22-2#2",
    sample %in% c("210729_006") ~ "SYG-TN13-E23-3")) %>%
  mutate(name = factor(name, levels = c("SYG-TN13-E22-1",
                                        "SYG-TN13-E22-2#1",
                                        "SYG-TN13-E22-2#2",
                                        "SYG-TN13-E23-3",
                                        "SYG-3",
                                        "CDG-062")))

# highlight the compounds
highlight_GJB_four <-
  meth_all %>%
  mutate(retention = round(Path, 1)) %>%
  group_by(sample, retention) %>%
  arrange(desc(File)) %>%
  #mutate(retention = ifelse(Sample %in% c("45 cm", "95 cm", "145 cm")& retention == 47.2,
   #                         47.3, retention)) %>%
  distinct(retention, .keep_all = TRUE) %>%
  mutate(FA = case_when(
    retention == 21.6 ~ "C14:0",
    retention == 23.6 ~ "C15:0",
    retention == 25.5 ~ "C16:0",
    retention == 27.3 ~ "C17:0",
    retention == 28.4 ~ "C18:1",
    retention == 29.1 ~ "C18:0",
    retention == 32.4 ~ "C20:0",
    retention == 35.4 ~ "C22:0")) %>%
  filter(name %in% c("SYG-3", "SYG-TN13-E22-2#1", "SYG-TN13-E22-1", "SYG-TN13-E23-3")) %>%
  filter(File > 100000)

highlight_IS <-
  meth_all %>%
  mutate(retention = round(Path, 1)) %>%
  group_by(sample, retention) %>%
  arrange(desc(File)) %>%
  distinct(retention, .keep_all = TRUE) %>%
  mutate(FA = case_when(
    retention == 33.9 ~ "IS",
    retention == 41.3 ~ "IS")) %>%
  filter(name %in% c("SYG-3", "SYG-TN13-E22-2#1", "SYG-TN13-E22-1", "SYG-TN13-E23-3"))

# highlight the compounds in the potsherd with high concentration
highlight_SYG2_2 <-
  meth_all %>%
  mutate(retention = round(Path, 1)) %>%
  group_by(sample, retention) %>%
  arrange(desc(File)) %>%
  distinct(retention, .keep_all = TRUE) %>%
  mutate(FA = case_when(
    retention == 21.6 ~ "C14:0",
    retention == 25.5 ~ "C16:0",
    retention == 27.3 ~ "C17:0",
    retention == 28.4 ~ "C18:1",
    retention == 29.1 ~ "C18:0",
    retention == 31.0 ~ "C19:0",
    retention == 32.4 ~ "C20:0",
    retention == 35.4 ~ "C22:0",
    retention == 38.2 ~ "C24:0",
    retention == 40.9 ~ "C26:0",
    retention == 33.9 ~ "IS",
    retention == 41.3 ~ "IS",
    retention == 47.8 ~ "C:32",)) %>%
  filter(name == "SYG-TN13-E22-2#2")

# plot the first four
meth_all %>%
  filter(name %in% c("SYG-3", "SYG-TN13-E22-2#1", "SYG-TN13-E22-1", "SYG-TN13-E23-3")) %>%
  mutate(File = ifelse(File > 550000, 550000, File)) %>%
  ggplot(aes(Path, File)) +
  geom_line(size = 0.3) +
  geom_text(data = highlight_GJB_four, #ggrepel::geom_text_repel
            aes(label = FA, angle = 90),
            size = 3,
            nudge_y = 87000,
            show.legend = FALSE) +
  geom_text(data = highlight_IS,
            aes(label = FA),
            size = 3,
            nudge_y = 40000,
            show.legend = FALSE) +
  geom_segment(aes(x= 41.1, xend= 41.5, y= 550000, yend= 550000), size =0.3) +
  geom_segment(aes(x= 41.1, xend= 41.5, y= 560000, yend= 560000), size =0.3) +
  annotate("text", x = 41.3, y = 600000, label = "IS", size =3) +
  facet_wrap(~name, ncol = 1) +
  scale_y_continuous(labels = scales::comma_format(),
                     limits = c(0, 600000),
                     breaks = seq(0, 600000, 200000)) +
  scale_x_continuous(limits = c(20, 45),
                     expand = c(0, 0.5)) + # don't log, many peaks and distort real counts
  labs(x = "retention time", y = "Relative Intensity") +
  theme_minimal()

ggsave(here::here("analysis","figures", "chromatograms_four.png"),
       width = 8,
       height = 8,
       dpi = 360,
       units = "in")

# plot the one with more compounds
meth_all %>%
  filter(name == "SYG-TN13-E22-2#2") %>%
  ggplot(aes(Path, File)) +
  geom_line(size = 0.3) +
  geom_text(data = highlight_SYG2_2, #ggrepel::geom_text_repel
            aes(label = FA),
            size = 3,
            nudge_x = 0,
            nudge_y = 250000,
            show.legend = FALSE) +
  scale_y_continuous(labels = scales::comma_format(),
                     limits = c(0, 7500000),
                     breaks = seq(0, 7500000, 1000000)) +
  scale_x_continuous(limits = c(15, 50),
                     expand = c(0, 0.5)) + # don't log, many peaks and distort real counts
  theme_minimal()

ggsave(here::here("analysis","figures", "chromatograms_specific.png"),
       width = 8,
       height = 6,
       dpi = 360,
       units = "in")

