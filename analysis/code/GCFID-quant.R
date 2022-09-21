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
    sample %in% c("210729_001") ~ "Pot 1",
    sample %in% c("210729_002") ~ "Pot 2",
    sample %in% c("210729_003") ~ "Pot 3",
    sample %in% c("210729_004") ~ "Pot 4",
    sample %in% c("210729_005") ~ "Pot 5",
    sample %in% c("210729_006") ~ "Pot 6",
    sample %in% c("210729_007") ~ "Pot 7",
    sample %in% c("210729_008") ~ "blank")) %>%
  mutate(name = factor(name, levels = c("Pot 1",
                                        "Pot 2",
                                        "Pot 3",
                                        "Pot 4",
                                        "Pot 5",
                                        "Pot 6",
                                        "Pot 7",
                                        "blank")))

  #mutate(sample_label = paste0(period, " (", Sample, ")")) %>%
  #mutate(sample_label = factor(sample_label, levels = c("post-Angkor (0 cm)", "post-Angkor (10 cm)",
                                                    #    "Angkor: after temple construction (25 cm)",
                                                   #     "Angkor: after temple construction (45 cm)",
                                                  #      "Angkor: before temple construction (75 cm)",
                                                  #      "Angkor: before temple construction (95 cm)",
                                                  #      "pre-Angkor (135 cm)", "pre-Angkor (145 cm)")))

# highlight the compounds
highlight_meth_four <-
  meth_all %>%
  mutate(retention = round(Path, 1)) %>%
  group_by(sample, retention) %>%
  arrange(desc(File)) %>%
  #mutate(retention = ifelse(Sample %in% c("45 cm", "95 cm", "145 cm")& retention == 47.2,
   #                         47.3, retention)) %>%
  distinct(retention, .keep_all = TRUE) %>%
  mutate(FA = case_when(
    retention == 25.5 ~ "C16:0",
    retention == 29.1 ~ "C18:0",
    retention == 33.9 ~ "IS",
    retention == 41.3 ~ "IS")) %>%
  filter(name %in% c("Pot 1", "Pot 2", "Pot 3", "Pot 4", "blank"))

highlight_S5 <-
  meth_all %>%
  mutate(retention = round(Path, 1)) %>%
  group_by(sample, retention) %>%
  arrange(desc(File)) %>%
  #mutate(retention = ifelse(Sample %in% c("45 cm", "95 cm", "145 cm")& retention == 47.2,
  #                         47.3, retention)) %>%
  distinct(retention, .keep_all = TRUE) %>%
  mutate(FA = case_when(
    retention == 21.6 ~ "C14:0",
    retention == 25.5 ~ "C16:0",
    retention == 29.1 ~ "C18:0",
    retention == 32.4 ~ "C20:0",
    retention == 35.4 ~ "C22:0",
    retention == 38.2 ~ "C24:0",
    retention == 33.9 ~ "IS",
    retention == 41.3 ~ "IS")) %>%
  filter(name == "Pot 5")

# plot the first four
meth_all %>%
  filter(name %in% c("Pot 1", "Pot 2", "Pot 3", "Pot 4", "blank")) %>%
  filter(between(Path, 20, 45)) %>%
  ggplot(aes(Path, File)) +
  geom_line(size = 0.3) +
  geom_text(data = highlight_meth_four, #ggrepel::geom_text_repel
            aes(label = FA),
            size = 3,
            nudge_x = 0,
            nudge_y = 100000,
            show.legend = FALSE) +
  facet_wrap(~name, ncol = 1) +
  scale_y_continuous(labels = scales::comma_format(),
                     limits = c(0, 1250000),
                     breaks = seq(0, 1250000, 250000)) +
  scale_x_continuous(limits = c(20, 45),
                     expand = c(0, 0.5)) + # don't log, many peaks and distort real counts
  labs(x = "retention time", y = "intensity") +
  theme_minimal()

# plot the last one
meth_all %>%
  filter(name == "Pot 5") %>%
  filter(between(Path, 15, 50)) %>%
  ggplot(aes(Path, File)) +
  geom_line(size = 0.3) +
  geom_text(data = highlight_S5, #ggrepel::geom_text_repel
            aes(label = FA),
            size = 3,
            nudge_x = 0,
            nudge_y = 250000,
            show.legend = FALSE) +
  facet_wrap(~name, ncol = 1) +
  scale_y_continuous(labels = scales::comma_format(),
                     limits = c(0, 7500000),
                     breaks = seq(0, 7500000, 1000000)) +
  scale_x_continuous(limits = c(15, 50),
                     expand = c(0, 0.5)) + # don't log, many peaks and distort real counts
  theme_minimal()

ggsave(here::here("analysis","figures", "GC_FID_meth_all.png"),
       width = 12,
       height = 8,
       dpi = 360,
       units = "in")

#--------------------------------------------------------------
# from LiYing_meth_NF2_INT.CSV, after silylation and before methyl, only peaks
meth_int <- read_csv(here::here("analysis/data/raw_data/LiYing_meth_NF2_INT.CSV"),
                     col_names = paste0("column", 1:9))

meth_int <-
  meth_int %>%
  filter(!str_detect(column1, "INT") & !str_detect(column1, "Mon") & !str_detect(column1, "Path"))

names(meth_int) <- meth_int[2,]
ixd <- which(str_detect(meth_int$Peak, paste(c("P", "C"), collapse = "|")))
meth_int_dfs <- split(meth_int, cumsum(1:nrow(meth_int) %in% ixd))
meth_int_data <- meth_int_dfs[c(FALSE, TRUE)]

meth_int_name <- meth_int_dfs[c(TRUE, FALSE)]
names(meth_int_data) <- map(meth_int_name, ~.x$End[1])
names(meth_int_name) <- map(meth_int_name, ~.x$End[1])
meth_int_dfs <- map_df(meth_int_data, ~.x[-1, ] %>%
                         #mutate_all(as.numeric),
                         mutate(across(c(Peak:End, Height:`Pct Total`), as.numeric)),
                       .id = 'sample')

meth_int_dfs <-
  meth_int_dfs %>%
  mutate(sample = paste(parse_number(sample), "cm")) %>%
  mutate(sample = factor(sample, levels = c("0 cm", "10 cm", "25 cm", "45 cm",
                                            "75 cm", "95 cm", "135 cm", "145 cm"))) %>%
  mutate(period = case_when(
    sample %in% c("0 cm", "10 cm") ~ "post-Angkor",
    sample %in% c("25 cm", "45 cm") ~ "Angkor: after temple construction",
    sample %in% c("75 cm", "95 cm") ~ "Angkor: before temple construction",
    sample %in% c("135 cm", "145 cm") ~ "pre-Angkor")) %>%
  mutate(period = factor(period, levels = c("post-Angkor",
                                            "Angkor: after temple construction",
                                            "Angkor: before temple construction",
                                            "pre-Angkor")))

meth_int_dfs %>%
  filter(R.T.>= 10) %>%
  ggplot(aes(R.T., Height)) +
  geom_line() +
  labs(x = "Retention time",
       y = "Intensity") +
  theme_minimal() +
  facet_wrap(~sample + period, ncol = 2)
