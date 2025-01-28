library(tidyverse)
theme_set(theme_minimal())

combined <- read_tsv("analysis/combined/features.tsv")
ranks_combined <- read_tsv("analysis/combined/ranks.tsv")

make_comparison_plot <- function(feature_name, ylabel="{feature_name}"){
  combined %>%
    left_join(filter(ranks_combined, Feature_name==.env$feature_name), by="patient") %>%
    mutate(patient = glue::glue("{patient} ({rank})")) %>%
    ggplot(aes(x=ID, y=.data[[feature_name]], group=interaction(noise, seed))) +
      geom_line(alpha=.2) +
      geom_line(data = ~filter(., parse_number(noise)==0 | (noise=="snr50")), color="red") +
      facet_wrap(~patient, scales="free_x") +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
      xlab("frame") +
      ylab(str_glue(ylabel))
}

fig3 <- make_comparison_plot("wavelet-LLL_glcm_Idn", ylabel="normalized feature value ({feature_name})")
ggsave("figures/fig3.svg", plot=fig3, width=10, height=10)

rank_order_by_median <- ranks_combined %>%
  filter(!is.na(rank)) %>%
  group_by(Feature_name) %>%
  summarize(median_rank=median(rank), mad_rank=mad(rank)) %>%
  mutate(Feature_name = fct_reorder(Feature_name, median_rank)) %>%
  pull(Feature_name) %>%
  levels()

pdf("figures/fig3_supplement.pdf", onefile = TRUE)
for (i in rank_order_by_median) {
  print(make_comparison_plot(i))
}
dev.off()