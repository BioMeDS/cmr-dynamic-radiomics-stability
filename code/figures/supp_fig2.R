library(tidyverse)
theme_set(theme_light())

summarized <- read_csv("analysis/calculated_mae/summary.csv")

supp_fig2 <- summarized |>
  ggplot(aes(x=mean_mae, y=mean_within_mae)) +
    geom_point() +
    scale_x_continuous(limits=c(0,1)) +
    scale_y_continuous(limits=c(0,1)) +
    geom_hline(yintercept=.4) +
    geom_vline(xintercept = .4) +
    xlab("Mean average error (total)") +
    ylab("Mean average error (within noise class)")

ggsave("figures/supp_fig2.svg", plot=supp_fig2, width=10, height=8)
