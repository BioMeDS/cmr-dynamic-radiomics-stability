library(tidyverse)
library(corrr)
theme_set(theme_minimal())

load_csvs_in_folder <- function(paths) {
  df <- paths %>%
    map_df(
      ~read_csv(paste0(.x), show_col_types = FALSE) %>% mutate(slice = (.x %>% file.path %>% dirname %>% basename)), 
      .id = "source")
  return(df)
}
df <- load_csvs_in_folder(snakemake@input)
df <- df %>% group_by(Feature_name, slice) %>% summarise(mae = mean(mae)) %>% group_by(slice) %>% mutate(rank = dense_rank(mae))
plot <- df %>%
  filter(!is.na(rank)) %>%
  select(-rank) %>%
  #mutate(mae = pmin(mae, 30)) %>%
  pivot_wider(names_from=slice, values_from=mae) %>%
  correlate(method = "spearman") %>%
  #shave() %>%
  stretch(na.rm = TRUE) %>%
  rename(`Spearman r` = r) %>%
  ggplot(aes(x,y,fill=`Spearman r`)) +
    geom_tile() +
    theme_minimal() +
    scale_fill_gradient2(low="blue", high="red") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab("") + ylab("")

plot1 <- df %>%
  select(-mae) %>%
  pivot_wider(names_from=slice, values_from=rank) %>%
  group_by(Feature_name) %>%
  ggplot(aes(x = `6`, y = `7`)) +
	 geom_point() +
	 theme_minimal()

plot2 <- df %>%
  select(-mae) %>%
  pivot_wider(names_from=slice, values_from=rank) %>%
  group_by(Feature_name) %>%
  ggplot(aes(x = `7`, y = `8`)) +
	 geom_point() +
	 theme_minimal()

plot3 <- df %>%
  select(-mae) %>%
  pivot_wider(names_from=slice, values_from=rank) %>%
  group_by(Feature_name) %>%
  ggplot(aes(x = `6`, y = `8`)) +
	 geom_point() +
	 theme_minimal()

ggsave("figures/rev/spearman.svg", plot=plot)
ggsave("figures/rev/slice6vs7.svg", plot=plot1)
ggsave("figures/rev/slice7vs8.svg", plot=plot2)
ggsave("figures/rev/slice6vs8.svg", plot=plot3)
