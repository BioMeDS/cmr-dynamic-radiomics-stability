library(tidyverse)

save_single_rank_table <- function(table, output_path) {
  read_csv(table) %>%
    separate(Feature_name,
             into = c("filter", "feature_class", "feature"),
             sep = "_",
             remove = FALSE) %>%
    filter(feature_class != "shape") %>%
    group_by(Feature_name, filter, feature_class, feature) %>%
    summarise(mae = mean(mae)) %>%
    mutate(rank = dense_rank((mae))) %>%
    group_by(filter, feature_class, feature, mae, rank, Feature_name) %>%
    distinct(Feature_name, mae, rank) %>%
    write_csv(output_path)
}


save_single_rank_table(snakemake@input[[1]], snakemake@output[[1]])

# save_single_rank_table("analysis/calculated_mae/subject/mae_Proband X11.csv", "analysis/tables/subject/mae_Proband X11.csv")