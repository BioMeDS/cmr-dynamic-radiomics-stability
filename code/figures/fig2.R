library(tidyverse)
library(patchwork)

theme_set(theme_light())
data <- read_csv("analysis/calculated_mae/mrxcat_simulation/mae.csv")
curves <- read_csv("analysis/features_normalized/mrxcat_simulation/features.csv")

mae <- data |>
  filter(!str_detect(Feature_name, "_shape_")) |>
  filter(!is.na(mae)) |>
  group_by(Feature_name) |>
  summarize(mae=mean(mae))

features_of_interest <- c(
  "wavelet-LLL_glcm_Idn",
  "wavelet-LHL_glcm_JointAverage",
  "wavelet-HLH_glcm_Imc2",
  "wavelet-HHH_firstorder_10Percentile"
)

mae_of_interest <- mae |>
  filter(Feature_name %in% features_of_interest) |>
  arrange(mae)

p1 <- mae |>
  ggplot(aes(x=mae)) + geom_histogram(bins=100) +
  geom_segment(data=mae_of_interest, aes(x=mae,xend=mae), y=-10, yend=-1, color="black", arrow = arrow(length = unit(0.2, "cm")))

p2 <- curves |>
  mutate(file=str_remove(file,".csv")) |>
  mutate(file=str_remove(file,"snr")) |>
  separate(file, into=c("snr","rep"), sep="_") |>
  mutate(snr = fct_reorder(snr, as.numeric(snr))) |>
  select(snr, rep, ID, all_of(features_of_interest)) %>%
  pivot_longer(4:ncol(.), names_to="Feature_name", values_to="value") |>
  left_join(mae, by="Feature_name") |>
  mutate(fn = str_glue("{Feature_name} (mae={round(mae,2)})")) |>
  mutate(fn=fct_reorder(fn,mae)) |>
  ggplot(aes(x=ID, y=value, color=snr, group=interaction(snr, rep))) +
    geom_line() +
    geom_line(data=~filter(., snr=="50"), size=1) +
    facet_wrap(~fn, scales="free_y") +
    scale_color_brewer(palette = "Set1") +
    ylab("normalized feature value") +
    xlab("frame") +
    # make facet text black
    theme(strip.text = element_text(color="black"))

ggsave("figures/fig2.svg", plot=p1/p2, width=10, height=10)