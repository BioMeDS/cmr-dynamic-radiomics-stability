library(tidyverse)
library(tidymodels)

set.seed(20260119)

files <- map_chr(snakemake@input, c)

data <- read_csv(files)

# define minimum validation accuracy for a feature to be considered informative
min_performance = 0.33

# split each group into train/val/test as 12/4/4
metadata <- read_tsv("data/ACDC_dataset/group.tsv", col_names = c("patient", "group")) |>
	mutate(set = c(rep(c(rep("train", 12), rep("val", 4), rep("test", 4)),5)))

ranks <- read_tsv("figures/supp_tab1.tsv")

# function to subsample each curve to 12 (equally spaced) timepoints
subsampleTo <- function(n, to=12){
	round(seq(0,n,length.out=to))
}

data_subsampled <- data |>
	select(-c(image_path, segmentation_path, extraction_ID)) |>
	filter(ID %in% subsampleTo(max(ID),12), .by=file) |>
	mutate(ID = str_glue("f{row_number()}"), .by=file) |>
	mutate(
		file=str_remove(file, ".csv"),
	) |>
	separate(file, into=c("patient","noise","seed"), sep="_") |>
	left_join(metadata, by="patient") |>
	# Remove noisy curves from the test set
	filter(set != "test" | noise == "0.000") |>
	mutate(group = as.factor(group))

# remove features containing NA
complete_features <- data_subsampled |>
	pivot_longer(names_to="feature", contains("_")) |>
	count(feature, nac=is.na(value)) |>
	pivot_wider(names_from=nac, values_from=n, values_fill=0) |>
	filter(`TRUE`==0) |>
	pull(feature)

data_nested <- data_subsampled |>
	pivot_longer(names_to="feature", contains("_")) |>
	filter(feature %in% complete_features) |>
	nest_by(feature) |>
    mutate(data = list(pivot_wider(data, names_from=ID, values_from=value))) |>
	mutate(
		train = list(filter(data, set=="train")),
		val = list(filter(data, set=="val")),
		test = list(filter(data, set=="test")),
	) |>
	select(-data)

decision_tree_fit <- workflow() |>
	add_model(decision_tree(mode="classification")) |>
	add_variables(outcomes = group, predictors = f1:f12)

data_acc <- data_nested |>
	mutate(
		my_fit = list(fit(decision_tree_fit, train)),
		val = list(augment(my_fit, val)),
		test = list(augment(my_fit, test)),
		val_acc = accuracy(val, group, .pred_class) |> pull(.estimate),
		test_acc = accuracy(test, group, .pred_class) |> pull(.estimate)
	) |>
	select(feature, val_acc, test_acc) |>
    ungroup() |>
	left_join(ranks |> select(feature=Feature_name, rank), by="feature")

write_tsv(data_acc, "figures/supp_tab2.tsv")