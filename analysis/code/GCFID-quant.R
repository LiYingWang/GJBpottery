library(tidyverse)

# load the GC data
meth <- read_csv(here::here("analysis/data/raw_data/LY_21_08_20_pot_chrom.CSV"))

# tidy data
ixd <- which(str_detect(meth$Path, "C"))
meths <- split(meth, cumsum(1:nrow(meth) %in% ixd))
meth_data <- meths[c(FALSE, TRUE)]
meth_name <- meths[c(TRUE, FALSE)]
names(meth_data) <- map(meth_data, ~.x$Path[1])
names(meth_name) <-map(meth_name, ~.x$Sample[1])
meths <- map_df(meth_data, ~.x[-1, ] %>%
                   mutate_all(as.numeric), .id = 'sample')

# assign sample name
meth_all <-
  meths %>%
  mutate(sample = paste0(parse_number(sample),
                         str_extract(sample, pattern = "_00[0-9]"))) %>%
  mutate(name = case_when(sample %in% c("210729_001") ~ "CDG-062",
                          sample %in% c("210729_002") ~ "SYG-3",
                          sample %in% c("210729_003") ~ "SYG-TN13-E22-2#1",
                          sample %in% c("210729_004") ~ "SYG-TN13-E22-1",
                          sample %in% c("210729_005") ~ "SYG-TN13-E22-2#2",
                          sample %in% c("210729_006") ~ "SYG-TN13-E23-3",
                          sample %in% c("210729_008") ~ "blank"))%>%
  mutate(name = factor(name, levels = c("SYG-TN13-E22-1",
                                        "SYG-TN13-E22-2#1",
                                        "SYG-TN13-E22-2#2",
                                        "SYG-TN13-E23-3",
                                        "SYG-3",
                                        "CDG-062",
                                        "blank")))

# round all peaks in order
all_peaks <-
  meth_all %>%
  mutate(retention = round(Path, 1)) %>%
  group_by(sample, retention) %>%
  arrange(desc(File)) %>%
  distinct(retention, .keep_all = TRUE)

# highlight the compounds
highlight_common <-
  all_peaks %>%
  mutate(FA = case_when(
    retention == 21.6 ~ "C[14:0]",
    retention == 23.6 ~ "C[15:0]",
    retention == 25.5 ~ "C[16:0]",
    retention == 27.3 ~ "C[17:0]",
    retention == 28.4 ~ "C[18:1]",
    retention == 29.1 ~ "C[18:0]",
    retention == 32.3 ~ "C[20:0]",
    retention == 35.4 ~ "C[22:0]")) %>%
  filter(name %in% c("SYG-3", "SYG-TN13-E22-2#1", "SYG-TN13-E22-1", "SYG-TN13-E23-3", "CDG-062")) %>%
  filter(File > 100000)

highlight_IS <-
  all_peaks %>%
  mutate(FA = case_when(
    retention == 33.9 ~ "IS",
    retention == 41.3 ~ "IS")) %>%
  filter(name %in% c("SYG-3", "SYG-TN13-E22-2#1", "SYG-TN13-E22-1", "SYG-TN13-E23-3", "CDG-062"))

# highlight the compounds in the potsherd with high concentration
highlight_SYG2_2 <-
  all_peaks %>%
  mutate(FA = case_when(
    retention == 21.6 ~ "C[14:0]",
    retention == 25.5 ~ "C[16:0]",
    retention == 27.3 ~ "C[17:0]",
    retention == 28.4 ~ "C[18:1]",
    retention == 29.1 ~ "C[18:0]",
    retention == 31.0 ~ "C[19:0]",
    retention == 32.4 ~ "C[20:0]",
    retention == 35.4 ~ "C[22:0]",
    retention == 36.8 ~ "C[23:0]",
    retention == 38.2 ~ "C[24:0]",
    retention == 40.9 ~ "C[26:0]",
    retention == 33.9 ~ "IS",
    retention == 41.3 ~ "IS",
    retention == 45.6 ~ "C[30:0]",
    retention == 47.8 ~ "C[32:0]")) %>%
  filter(name == "SYG-TN13-E22-2#2")

insert_minor <- function(major_labs, n_minor) {labs <-
  c( sapply( major_labs, function(x) c(x, rep("", 4) ) ) )
labs[1:(length(labs)-n_minor)]}

# plot the first five
meth_all %>%
  filter(name %in% c("SYG-3","SYG-TN13-E22-2#1","SYG-TN13-E22-1","SYG-TN13-E23-3","CDG-062")) %>%
  mutate(File = ifelse(File > 1000000, 1000000, File)) %>%
  ggplot(aes(Path, File)) +
  geom_line(size = 0.3) +
  ggrepel::geom_text_repel(data = highlight_common, #ggrepel::geom_text_repel
                           aes(label = FA, angle = 0), size = 3, nudge_y = 115000,
                           show.legend = FALSE, parse = TRUE) +
  geom_text(data = highlight_IS,
            aes(label = FA), size = 3, nudge_y = 100000, show.legend = FALSE) +
  geom_segment(aes(x= 41.1, xend= 41.5, y= 1000000, yend= 1000000), size =0.3) +
  geom_segment(aes(x= 41.1, xend= 41.5, y= 1030000, yend= 1030000), size =0.3) +
  annotate("text", x = 41.3, y = 1100000, label = "IS", size =3) +
  facet_wrap(~name, ncol = 1) +
  scale_y_continuous(limits = c(0, 1160000),
                     breaks = seq(0, 1100000, 200000),
                     expand = c(0, 0)) +
  scale_x_continuous(limits = c(10, 47),
                     breaks = seq(10, 50, 1),
                     labels = insert_minor(seq(10, 50, 5), 4),
                     expand = c(0, 0)) +
  labs(x = "Retention time (mins)", y = expression("Relative abundance " %->% "")) +
  theme_classic() +
  theme(axis.line.y = element_blank(), strip.background = element_blank(), # remove title box
        axis.ticks.y=element_blank(), axis.text.y = element_blank())

ggsave(here::here("analysis","figures", "chromatograms_four.png"),
       bg = "white", width = 7, height = 8, dpi = 300, units = "in")

# plot the one with high yield (GC program 1)
chroma_specific <-
  meth_all %>%
  filter(name %in% c("SYG-TN13-E22-2#2")) %>%
  ggplot(aes(Path, File)) +
  geom_line(size = 0.3) +
  ggrepel::geom_text_repel(data = highlight_SYG2_2, #ggrepel::geom_text_repel
                           aes(label = FA),
                           size = 3.5,
                           nudge_y = 250000,
                           show.legend = FALSE,
                           parse = TRUE) +
  scale_y_continuous(labels = scales::comma_format(),
                     limits = c(0, 7500000),
                     breaks = seq(0, 7500000, 1000000),
                     expand = c(0, 0.5)) +
  scale_x_continuous(limits = c(10, 50),
                     breaks = seq(10, 50, 5),
                     minor_breaks = seq(10, 45, 1),
                     expand = c(0, 0.5)) +
  labs(title = "SYG-TN13-E22-2#2", x = "Retention time (mins)", y = "Relative abundance") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave(here::here("analysis","figures", "chromatograms_specific.png"),
       bg = "white", width = 8, height = 5, dpi = 360, units = "in")

# SYG-E22-22-sim and scan mode for miliacin (GC program 2)
meth_SYG22 <- read_csv(here::here("analysis/data/raw_data/LY_pot_SYG_E22_22_chrom.CSV"))

# tidy data
ixd_2 <- which(str_detect(meth_SYG22$Path, "D"))
meths_2 <- split(meth_SYG22, cumsum(1:nrow(meth_SYG22) %in% ixd))
meth_data_2 <- meths_2[c(FALSE, TRUE)]
meth_name_2 <- meths_2[c(TRUE, FALSE)]
names(meth_data_2) <- map(meth_data_2, ~.x$Path[1])
names(meth_name_2) <-map(meth_name_2, ~.x$Sample[1])
meths_2 <- map_df(meth_data_2, ~.x[-1, ] %>%
                  mutate_all(as.numeric), .id = 'sample')

# peaks in SYG-E22-2#2
peaks_SYGE22 <-
  meths_2 %>%
  mutate(retention = round(Path, 1)) %>%
  group_by(sample, retention) %>%
  arrange(desc(File)) %>%
  distinct(retention, .keep_all = TRUE)

highlight_SYG_FA <-
  peaks_SYGE22 %>%
  mutate(FA = case_when(
    retention == 17.7 ~ "C[15:0]",
    retention == 18.7 ~ "C[16:0]",
    retention == 19.7 ~ "C[17:0]",
    retention == 20.4 ~ "C[18:1]",
    retention == 20.6 ~ "C[18:0]",
    retention == 21.5 ~ "C[19:0]",
    retention == 22.4 ~ "C[20:0]",
    retention == 24.0 ~ "C[22:0]",
    retention == 24.8 ~ "C[23:0]",
    retention == 25.5 ~ "C[24:0]",
    retention == 26.3 ~ "C[25:0]",
    retention == 27.0 ~ "C[26:0]",
    retention == 28.3 ~ "C[28:0]"))

highlight_std <-
  peaks_SYGE22 %>%
  mutate(FA = case_when(
    retention == 23.3 ~ "IS",
    retention == 30.5~ "IS"))

highlight_small <-
  peaks_SYGE22 %>%
  mutate(FA = case_when(
    retention == 14.2 ~ "C[12:0]",
    retention == 15.4 ~ "C[13:0]",
    retention == 16.6 ~ "C[14:0]"))

highlight_mili <-
 peaks_SYGE22 %>%
 mutate(FA = case_when(retention == 29.3 ~ "beta",
                       retention == 29.4 ~ "M"))

## plot the one with high yield (GC program 2)
chroma_specific_2 <-
  meths_2 %>%
  ggplot(aes(Path, File)) +
  geom_line(size = 0.3) +
  ggrepel::geom_text_repel(data = highlight_SYG_FA , #ggrepel::geom_text_repel
                           aes(label = FA, angle = 90), size = 3.5, nudge_y = 70000,
                           show.legend = FALSE, parse = TRUE) +
  ggrepel::geom_text_repel(data = highlight_std, #ggrepel::geom_text_repel
                           aes(label = FA), size = 3.5, nudge_y = 70000,
                           show.legend = FALSE, parse = TRUE) +
  ggrepel::geom_text_repel(data = highlight_small, #ggrepel::geom_text_repel
                           aes(label = FA, angle = 90), size = 3.5, segment.alpha = .25,
                           nudge_y = 350000, #min.segment.length = 300000,
                           show.legend = FALSE, parse = TRUE) +
  ggrepel::geom_text_repel(data = highlight_mili, #ggrepel::geom_text_repel
                           aes(label = FA), size = 3.5, segment.alpha = .25,
                           nudge_y = 350000, min.segment.length = 0,
                           show.legend = FALSE, parse = TRUE) +
  scale_y_continuous(labels = scales::comma_format(),
                     limits = c(0, 4600000),
                     breaks = seq(0, 4600000, 1000000),
                     expand = c(0, 0)) +
  scale_x_continuous(limits = c(13.5, 32.5),
                     labels = insert_minor(seq(10, 35, 5), 4),
                     breaks = seq(10, 35, 1),
                     expand = c(0, 0)) +
  labs(x = "Retention time (mins)", y = expression("Relative abundance " %->% "")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), axis.line.y = element_blank(),
        axis.ticks.y=element_blank(), axis.text.y = element_blank())

## miliacin
miliacin <- read_csv(here::here("analysis/data/raw_data/LY_miliacin_std_chrom.CSV"))
ixd_3 <- which(str_detect(miliacin$Path, "D"))
meths_3 <- split(miliacin, cumsum(1:nrow(miliacin) %in% ixd))
meth_data_3 <- meths_3[c(FALSE, TRUE)]
meth_name_3 <- meths_3[c(TRUE, FALSE)]
names(meth_data_3) <- map(meth_data_3, ~.x$Path[1])
names(meth_name_3) <-map(meth_name_3, ~.x$Sample[1])
meths_3 <- map_df(meth_data_3, ~.x[-1, ] %>%
                    mutate_all(as.numeric), .id = 'sample')

highlight_m <-
  meths_3  %>%
  mutate(retention = round(Path, 1)) %>%
  group_by(sample, retention) %>%
  arrange(desc(File)) %>%
  distinct(retention, .keep_all = TRUE) %>%
  mutate(std = case_when(retention == 29.4 ~ "M", retention == 30.5 ~ "IS"))

std_miliacin <-
  meths_3%>%
  ggplot(aes(Path, File)) +
  geom_line(size = 0.3) +
  ggrepel::geom_text_repel(data = highlight_m, #ggrepel::geom_text_repel
                           aes(label = std, angle = 0), size = 3.5, nudge_y = 100000,
                           show.legend = FALSE, parse = TRUE) +
  scale_y_continuous(labels = scales::comma_format(),
                     limits = c(0, 4600000),
                     breaks = seq(0, 4600000, 1000000),
                     expand = c(0,0)) +
  scale_x_continuous(limits = c(13.5, 32.5),
                     labels = insert_minor(seq(10, 35, 5), 4),
                     breaks = seq(10, 35, 1),
                     expand = c(0, 0)) +
  labs(x = "Retention time (mins)", y = expression("Relative abundance " %->% "")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), axis.line.y = element_blank(),
        axis.ticks.y=element_blank(), axis.text.y = element_blank())

library(cowplot)
SYG_with_miliacin <-
  plot_grid(std_miliacin,
            chroma_specific_2 ,
            labels = c('A', 'B'),
            ncol =1,
            align = 'v',
            label_size = 12)
