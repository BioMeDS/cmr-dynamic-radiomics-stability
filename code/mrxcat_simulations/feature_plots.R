library(tidyverse)
library(ggridges)
library(munsell)

s1_palette <- c("#CC29C7", "#993D96", "#FF01F6", "#332933")
s2_palette <- c("#B86725", "#855935", "#EB6800", "#332D29")
s3_palette <- c("#26BDBD", "#378A8A", "#00F0EE", "#293333")
s4_palette <- c("#2BBD26", "#3A8A37", "#09F000", "#293329")
s5_palette <- c("#B8B025", "#858135", "#EBE100", "#333229")
s6_palette <- c("#B82A25", "#853835", "#EB0800", "#332929")
s7_palette <- c("#f2f0f7", "#cbc9e2", "#9e9ac8", "#6a51a3")
class_colors = c("#FF01F6","#EB6800","#00F0EE", "#09F000", "#EBE100", "#EB0800")

palettes <- list(s1_palette,
                 s2_palette,
                 s3_palette,
                 s4_palette,
                 s5_palette,
                 s6_palette,
                 s7_palette)

load_csv_and_factorize <- function(path_csv) {
  data <- read_csv(path_csv) %>%
    separate(file, into = c("snr", "ending"),
             sep = "_", remove = FALSE) %>%
    mutate(snr = as.numeric(str_remove(snr, "snr"))) %>%
    mutate(snr = as_factor(snr))
  return(data)
}

plot_top_12_features <- function(data, filename) {
  top_12_candidates <- data %>%
    select(snr, `original_shape_Elongation`:`wavelet-LLL_ngtdm_Strength`) %>%
    pivot_longer(names_to = "feature", values_to = "value", -snr) %>%
    group_by(snr, feature) %>%
    summarize(value = mean(value)) %>%
    pivot_wider(names_from = snr, values_from = value) %>%
    mutate(distance = (`30` - `5`) / `30`) %>%
    arrange(-abs(distance)) %>%
    head(12) %>%
    pull(feature)
  data %>%
    select(ID, snr, file, top_12_candidates) %>%
    pivot_longer(names_to = "feature",
                 values_to = "value",
                 top_12_candidates) %>%
    ggplot(aes(x = ID, y = value, color = snr)) +
    geom_line(aes(group = file)) +
    facet_wrap(~feature, scale = "free_y")
  ggsave(filename = filename %>%
           file.path %>%
           basename,
         path = filename %>%
           file.path %>%
           dirname,
         device = "png",
         width = 20)
}

plot_snr_curves <- function(data, output_path, palettes) {
  data %>%
    select(ID, snr, file,
           `original_shape_Elongation`:`wavelet-LLL_ngtdm_Strength`) %>%
    pivot_longer(names_to = "feature",
                 values_to = "value",
                 `original_shape_Elongation`:`wavelet-LLL_ngtdm_Strength`) %>%
    group_split(feature) %>%
    map(~{
          if (str_detect(first(.x$feature), pattern = "firstorder")) {
            color_palette <- palettes[[1]]
          } else if (str_detect(first(.x$feature), pattern = "glcm")) {
            color_palette <- palettes[[2]]
          } else if (str_detect(first(.x$feature), pattern = "gldm")) {
            color_palette <- palettes[[3]]
          } else if (str_detect(first(.x$feature), pattern = "glrlm")) {
            color_palette <- palettes[[4]]
          } else if (str_detect(first(.x$feature), pattern = "glszm")) {
            color_palette <- palettes[[5]]
          } else if (str_detect(first(.x$feature), pattern = "ngtdm")) {
            color_palette <- palettes[[6]]
          } else if (str_detect(first(.x$feature), pattern = "shape")) {
            color_palette <- palettes[[7]]
          } else {
            color_palette <- palettes[[7]]
          }
          plot <- ggplot(.x) + aes(x = ID, y = value, color = snr) +
            geom_line(aes(group = file)) +
            labs(y = first(.x$feature)) +
            scale_color_manual(values = color_palette) +
            theme_classic()
          ggsave(paste0(first(.x$feature), ".png"),
                 plot,
                 path = output_path,
                 width = 10)})
}

create_snr_only_table <- function(path) {
  data <- read_csv(path)
  data <- data %>%
    separate(Dataframe_1,
             into = c("snr_a", "ending_a"),
             sep = "_",
             remove = FALSE) %>%
    mutate(snr_a = as_factor(as.numeric(str_remove(snr_a, "snr"))))
  data <- data %>%
    separate(Dataframe_2,
             into = c("snr_b", "ending_b"),
             sep = "_",
             remove = FALSE) %>%
    mutate(snr_b = as_factor(as.numeric(str_remove(snr_b, "snr"))))
  # Define the SNR values
  snr_values <- c(5, 10, 20, 30)
  # Create an empty list to store the data frames
  snr_data <- list()
  # Iterate over the SNR values
  for (x in seq(from = 1, to = 4, by = 1)) {
    # Filter the data based on SNR values
    filtered_data <- filter(data, snr_a == snr_values[[x]], snr_b == snr_values[[x]])
    # Group and summarize the filtered data
    summarized_data <- filtered_data %>%
      group_by(Feature_name) %>%
      summarise(mae = mean(mae)) %>%
      group_by(mae)
    # Add the SNR column to the summarized data
    summarized_data$SNR <- as.character(snr_values[[x]])
    # Store the summarized data in the list
    snr_data[[x]] <- summarized_data
  }
  snr_mean <- rbind(snr_data[[1]], snr_data[[2]], snr_data[[3]], snr_data[[4]])
  return(snr_mean)
}

create_table_for_snr_plots <- function(path, snr_mean) {
  mae_data <- read_csv(path) %>%
    group_by(Feature_name) %>%
    summarise(mae = mean(mae)) %>%
    group_by(mae)
  mae_data <- mae_data %>%
    separate(Feature_name,
             into = c("filter", "feature_class", "feature"),
             sep = "_",
             remove = FALSE)
  snr_mean <- snr_mean %>%
    separate(Feature_name,
             into = c("filter", "feature_class", "feature"),
             sep = "_",
             remove = FALSE)
  combined_table <- left_join(snr_mean,
                              mae_data,
                              by = c("Feature_name",
                                     "filter",
                                     "feature_class",
                                     "feature"))
  return(combined_table)
}

save_rank_table <- function(data, output_path) {
  data %>%
    filter(feature_class != "shape") %>%
    mutate(rank = dense_rank((mae.y))) %>%
    group_by(feature_class, mae.y, rank, Feature_name) %>%
    write_csv(output_path)
}

create_plots <- function(data, output1, output2, output3) {
  data %>%
    ggplot(aes(x = mae.y, y = `mae.x`, group = SNR, color = SNR)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    xlab("Mean average error (MAE) per signal to noise (SNR) class") +
    ylab("Mean average error (MAE) total")
  ggsave(filename = output1 %>%
           file.path %>%
           basename,
         path = output1 %>%
           file.path %>%
           dirname,
         device = "png",
         width = 10)
  data %>%
    slice_min(mae.y, prop = .90) %>%
    group_by(Feature_name, mae.y) %>%
    summarize(mae.x = mean(mae.x)) %>%
    ggplot(aes(x = mae.y, y = `mae.x`)) +
    geom_point() + geom_vline(xintercept = 0.25) +
    geom_hline(yintercept = 0.3) +
    xlab("Mean average error (MAE) per signal to noise (SNR) class") +
    ylab("Mean average error (MAE) total")
  ggsave(filename = output2 %>%
           file.path %>%
           basename,
         path = output2 %>%
           file.path %>%
           dirname,
         device = "png",
         width = 10)
  data %>%
    mutate(rank = dense_rank((mae.y))) %>%
    group_by(feature_class, mae.y, rank, Feature_name) %>%
    filter(feature_class != "shape") %>%
    summarize(mae.x = mean(mae.x)) %>%
    ggplot(aes(x = rank, y = feature_class, fill = feature_class)) +
    geom_density_ridges(
      jittered_points = TRUE,
      position = position_points_jitter(width = 0.05, height = 0),
      point_shape = "|", point_size = 3, point_alpha = 1, alpha = 0.7,
      scale = 0.7
    ) +
    xlab("Rank of features") +
    ylab("Feature Class") +
    labs(fill = "Feature Class") +
    theme(text = element_text(size = 14)) +
    scale_fill_manual(values = class_colors)
  ggsave(filename = output3 %>%
           file.path %>%
           basename,
         path = output3 %>%
           file.path %>%
           dirname,
         device = "png",
         width = 10)
}


data_curve_plots <- load_csv_and_factorize(snakemake@input[[1]]) # input: features.csv (normalized)

plot_top_12_features(data_curve_plots,
                     snakemake@output[[1]])

plot_snr_curves(data_curve_plots,
                snakemake@output[[1]] %>%
                  file.path %>%
                  dirname,
                palettes)

data_summary_plots <- create_table_for_snr_plots(snakemake@input[[2]],
                                                 create_snr_only_table(snakemake@input[[2]])) # input mae.csv

save_rank_table(data_summary_plots, snakemake@output[[2]])

create_plots(data_summary_plots,
             snakemake@output[[2]],
             snakemake@output[[3]],
             snakemake@output[[4]])

# test <-create_snr_only_table("analysis/calculated_mae/mrxcat_simulation/mae.csv")

# table <- create_table_for_snr_plots("analysis/calculated_mae/mrxcat_simulation/mae.csv", test)
