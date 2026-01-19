library(tidyverse)
theme_set(theme_light())

# define minimum validation accuracy for a feature to be considered informative
min_performance = 0.33

data <- read_tsv("figures/supp_tab2.tsv")

supp_fig5 <- data |>
	filter(val_acc>=min_performance) |>
	ggplot(aes(rank, test_acc)) +
		geom_point() +
		geom_smooth(method="lm") +
		xlab("feature stability (rank)") +
		ylab("accuracy (test set)")

ggsave("figures/supp_fig5.svg", plot = supp_fig5, width = 10, height = 8)

data |>
	filter(val_acc>=min_performance) |>
	summarize(cor = cor(rank, test_acc, method = "spearman")) |>
	write_tsv("analysis/ml/correlation.tsv")

data |>
	filter(val_acc>=min_performance) |>
	arrange(rank) |>
	mutate(
		rank_notie = row_number(),
		range = case_when(
		rank_notie <= 15 ~ "top",
		length(rank_notie) - rank_notie <= 15 ~ "bottom",
		.default = "middle",
	)) |>
	summarize(median(test_acc), .by=range) |>
	write_tsv("analysis/ml/performance.tsv")