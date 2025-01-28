library(tidyverse)
library(ggridges)
theme_set(theme_light())

ranks_combined <- read_tsv("analysis/combined/ranks.tsv")

ranks_summary <- ranks_combined %>%
  filter(!is.na(rank)) %>%
  mutate(mae = pmin(mae, 30)) %>%
  summarize(mae=mean(mae), rank=median(rank), .by=Feature_name) %>%
  arrange(rank)

fig5 <- ranks_summary %>%
  separate(Feature_name, into=c("wavelet", "feature_class", "feature"), sep="_") %>%
    ggplot(aes(x = rank, y = feature_class, fill = feature_class)) +
    geom_density_ridges(
      jittered_points = TRUE,
      position = position_points_jitter(width = 0.05, height = 0),
      point_shape = "|", point_size = 3, point_alpha = 1, alpha = 0.7,
      scale = 0.7
    ) +
    xlab("Rank of Features") +
    ylab("Feature Class") +
    labs(fill = "Feature Class") +
    theme(text = element_text(size = 14), legend.position = "none") +
    scale_fill_brewer(palette = "Set2")

ggsave("figures/fig5.svg", plot=fig5, width=10, height=8)
