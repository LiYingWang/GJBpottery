# read data
nano <- read.csv(here::here("analysis/data/raw_data/Combined_Renamed_Seq.csv"),
                      header = F)

# tidy data
ixd <- which(str_detect(nano$V1, "seq"))
nano_split <- split(nano, cumsum(1:nrow(nano) %in% ixd))
names(nano_split) <- map(nano_split, ~.x$V1[1])
nano_all <- map_df(nano_split, ~.x[-1,])

# transfer the data layout
nano_long <-
  nano_all %>%
  pivot_longer(cols = starts_with(">"),
               names_to = "read",
               values_to = "sequence") %>%
  mutate(read = gsub('>', '', read)) %>%
  mutate(dup_w = lapply(str_extract_all(sequence, "(?i)([a-z])\\1+"), nchar)) %>%
  mutate(dup_w_combo = lapply(str_extract_all(sequence, "([ACGT][ACGT])\\1+"), nchar)) %>%
  mutate(length_w = nchar(sequence)) %>%
  mutate(dup_w_max = lapply(dup_w, max)) %>%
  mutate(dup_w_com_max = lapply(dup_w_combo, max)) %>%
  filter(dup_w_max < 10 & dup_w_com_max < 10) %>%
  filter(length_w >= 200)

# load a BLAST database (replace db with the location + name of the BLAST DB)
bl <- blast(db = ".16SMicrobialDB/16S_ribosomal_RNA")
bl
