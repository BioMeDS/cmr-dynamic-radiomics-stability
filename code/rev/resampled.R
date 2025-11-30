library(tidyverse)
library(ggridges)
library(munsell)
library(corrr)
theme_set(theme_minimal())

load_csvs_in_folder <- function(paths) {
  df <- paths %>%
    map_df(
      ~read_csv(paste0(.x), show_col_types = FALSE) %>% mutate(path = .x, filename = basename(.x), resampled = grepl(pattern="resampled", .x)))
  return(df)
}
df <- load_csvs_in_folder(snakemake@input) %>% mutate(resampled = ifelse(resampled, "resampled", "non-resampled"))

df1 <- df %>% separate(filename, into=c("_", "patient"), sep="_") %>%
           group_by(Feature_name, resampled, patient) %>%
           summarise(mae = mean(mae)) %>% group_by(resampled) %>%
	   mutate(rank = dense_rank(mae)) %>% filter(!is.na(rank)) %>% select(-rank) %>%
	   group_by(resampled) %>% group_split()
plot <- correlate(df1[[1]] %>% pivot_wider(names_from=c(patient), values_from=mae), df1[[2]] %>% pivot_wider(names_from=c(patient), values_from=mae), method = "spearman") %>%
        stretch(na.rm = TRUE) %>%
        rename(`Spearman r` = r) %>% write_csv("test.csv")
        ggplot(aes(x,`Spearman r`)) +
        geom_boxplot() +
        theme_minimal() +
          #scale_fill_gradient2(low="blue", high="red") +
          #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        xlab("") + ylab("")

plot1 <- df %>% group_by(Feature_name, resampled) %>% summarise(mae = mean(mae)) %>% group_by(resampled) %>% mutate(rank = dense_rank(mae)) %>%
  select(-mae) %>%
  pivot_wider(names_from=resampled, values_from=rank) %>%
  group_by(Feature_name) %>%
  ggplot(aes(x = `resampled`, y = `non-resampled`)) +
	 geom_point() +
	 theme_minimal()

ggsave("figures/rev/resampled_vs_unresampled.svg", plot=plot)
ggsave("figures/rev/resampled_unresampled_rank_comparison.svg", plot=plot1)
