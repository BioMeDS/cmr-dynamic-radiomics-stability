library(tidyverse)

mmae_and_ranks <- function(table){
  data <- read_csv(table, col_types = cols(
    ...1 = col_double(),
    Feature_name = col_character(),
    Dataframe_1 = col_character(),
    Dataframe_2 = col_character(),
    mae = col_double()
  ), name_repair = "unique_quiet")
  total_ranks <- data %>%
      filter(!str_detect(Feature_name, "shape")) %>%
      summarise(mae = mean(mae), .by=Feature_name) %>%
      mutate(rank = dense_rank((mae)))

  within_ranks <- data %>%
      filter(!str_detect(Feature_name, "shape")) %>%
      mutate(
        Dataframe_1 = str_remove(Dataframe_1, "\\d+.csv"),
        Dataframe_2 = str_remove(Dataframe_2, "\\d+.csv")
      ) %>%
      filter(Dataframe_1 == Dataframe_2) %>%
      summarise(within_mae = mean(mae), .by=Feature_name) %>%
      mutate(within_rank = dense_rank((within_mae)))

  left_join(total_ranks, within_ranks, by="Feature_name") %>%
      arrange(rank) %>%
      mutate(file=table, .before=1)
}

data <- dir("analysis/calculated_mae/", pattern="*.csv", recursive = TRUE, full.names = TRUE) %>%
    map(mmae_and_ranks) %>%
    bind_rows() 

summarized <- data %>%
  summarize(mean_mae=mean(mae), sd=sd(mae), rank=median(rank), mean_within_mae=mean(within_mae), within_sd=sd(within_mae), within_rank=median(within_rank), .by=Feature_name)

summarized %>% write_tsv("figures/supp_tab1.tsv")
