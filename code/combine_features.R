library(tidyverse)

enname <- function(x) {
  names(x) <- x
  x
}

sim <- read_csv("analysis/features_normalized/mrxcat_simulation/features.csv") %>%
  separate(file, into=c("noise", "seed"), sep="_") %>%
  mutate(patient="sim")

ranks_sim <- read_csv("analysis/tables/mrxcat_simulation/ranks.csv") %>% select(-filter, -feature_class, -feature) %>% mutate(patient="sim")

acdc <- dir("analysis/features_normalized/ACDC/", full.names = TRUE) %>% map_df(read_csv) %>%
  separate(file, into=c("patient", "noise", "seed"), sep="_")

ranks_acdc <- dir("analysis/tables/ACDC", full.names = TRUE) %>%
  enname() %>%
  map_df(read_csv, .id="patient") %>%
  mutate(patient = str_remove(patient, ".*ranks_"), patient=str_remove(patient, ".csv"))

bae <- dir("analysis/features_normalized/subject/", full.names = TRUE) %>% map_df(read_csv) %>%
  separate(file, into=c("patient", "noise", "seed"), sep="_")

ranks_bae <- dir("analysis/tables/subject", full.names = TRUE) %>%
  enname() %>%
  map_df(read_csv, .id="patient") %>%
  mutate(patient = str_remove(patient, ".*ranks_"), patient=str_remove(patient, ".csv"))

combined <- bind_rows(sim, acdc, bae)
ranks_combined <- bind_rows(ranks_sim, select(ranks_acdc, colnames(ranks_sim)), select(ranks_bae, colnames(ranks_sim)))

combined %>% write_tsv("analysis/combined/features.tsv")
ranks_combined %>% write_tsv("analysis/combined/ranks.tsv")