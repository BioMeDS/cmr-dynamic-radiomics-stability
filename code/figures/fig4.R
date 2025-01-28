library(tidyverse)
library(corrr)
theme_set(theme_minimal())

ranks_combined <- read_tsv("analysis/combined/ranks.tsv")

fig4 <- ranks_combined %>%
  filter(!is.na(rank)) %>%
  select(-rank) %>%
  #mutate(mae = pmin(mae, 30)) %>%
  pivot_wider(names_from=patient, values_from=mae) %>%
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

ggsave("figures/fig4.svg", plot=fig4, width=10, height=8)
