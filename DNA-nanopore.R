# read data
nano <- read.csv(here::here("analysis/data/raw_data/Combined_Renamed_Seq.csv"),
                      header = F)

# tidy data
ixd <- which(str_detect(nano$V1, "seq"))
nano_split <- split(nano, cumsum(1:nrow(nano) %in% ixd))
names(nano_split) <- map(nano_split, ~.x$V1[1])
nano_all <- map_df(nano_data_split, ~.x[-1,])

nano_long <-
  nano_all %>%
