library(tidyverse)
theme_set(theme_minimal())

ranks_combined <- read_tsv("analysis/combined/ranks.tsv")

supp_fig1 <- ranks_combined %>%
  filter(!is.na(rank)) %>%
  mutate(mae = pmin(mae, 30)) %>%
  mutate(Feature_name = fct_reorder(Feature_name, rank, median)) %>%
  ggplot(aes(x=Feature_name, y=rank)) +
    geom_point(size=.4, alpha=.2) +
    stat_summary(geom="pointrange", color="red", size=.2, fun.data=ggpubr::median_mad) +
    # remove x-axis tick labels
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    xlab("Feature") +
    ylab("Rank")

ggsave("figures/supp_fig1.svg", plot=supp_fig1, width=10, height=5)
