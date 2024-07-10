library(tidyverse)
library(ggridges)
library(munsell)

# This section defines color palettes used for plotting in the script.
# Each palette is stored as a vector of hexadecimal color codes.
# The class_colors vector stores additional color codes used for class-specific coloring.

s1_palette <- c("#CC29C7", "#993D96", "#FF01F6", "#332933")
s2_palette <- c("#B86725", "#855935", "#EB6800", "#332D29")
s3_palette <- c("#26BDBD", "#378A8A", "#00F0EE", "#293333")
s4_palette <- c("#2BBD26", "#3A8A37", "#09F000", "#293329")
s5_palette <- c("#B8B025", "#858135", "#EBE100", "#333229")
s6_palette <- c("#B82A25", "#853835", "#EB0800", "#332929")
s7_palette <- c("#f2f0f7", "#cbc9e2", "#9e9ac8", "#6a51a3")

class_colors <- c("#FF01F6","#EB6800","#00F0EE", "#09F000", "#EBE100", "#EB0800")

palettes <- list(s1_palette,
         s2_palette,
         s3_palette,
         s4_palette,
         s5_palette,
         s6_palette,
         s7_palette)

#' Function to prepare data for analysis
#'
#' This function takes a path as input and loads all CSV files in the folder specified by the path.
#' It then separates the file names into three columns: snr and ending.
#' The "ending" column is removed.
#'
#' @param path A string specifying the path to the folder containing the CSV file.
#'
#' @return A data frame with the modified data

prepare_data <- function(path) {
  data <- read_csv(path) %>%
    separate(file, into = c("snr", "ending"),
             sep = "_", remove = FALSE) %>%
    mutate(ending = NULL) %>%
    mutate(snr = as.numeric(str_remove(snr, "snr"))) %>%
    mutate(snr = as_factor(snr))
  return(data)
}

#' Plot Top 12 Features
#'
#' This function plots the top 12 features based on their stability across different noise levels.
#'
#' @param path The path to the CSV file.
#' @param filename The name of the output file.
#'
#' @return None
#'
#' @examples
#' plot_top_12_features("/path/to/data", "output.png")
plot_top_12_features <- function(path, filename) {
  top_12_candidates <- prepare_data(path) %>%
    select(snr, `original_shape_Elongation`:`wavelet-LLL_ngtdm_Strength`) %>%
    pivot_longer(names_to = "feature", values_to = "value", -snr) %>%
    group_by(snr, feature) %>%
    summarize(value = mean(value)) %>%
    pivot_wider(names_from = snr, values_from = value) %>%
    mutate(distance = (`30` - `5`) / `30`) %>%
    arrange(-abs(distance)) %>%
    head(12) %>%
    pull(feature)
  prepare_data(path) %>%
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

plot_snr_curves <- function(path, output_path, palettes) {
  prepare_data(path) %>%
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

#' Create SNR Only Table
#'
#' This function reads a CSV file from the specified path and creates a table
#' that summarizes the mean absolute error (MAE) values for different features
#' at different signal-to-noise ratio (SNR) levels.
#'
#' @param path The path to the CSV file.
#'
#' @return A data frame containing the summarized data.
#'
#' @examples
#' create_snr_only_table("/path/to/data.csv")
create_snr_only_table <- function(path) {
  data <- read_csv(path) %>%
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

#' Create a table for SNR plots
#'
#' This function reads a CSV file from the given path and creates a table
#' for signal-to-noise ratio (SNR) plots. It calculates the mean average error (MAE)
#' for each feature and groups them by feature name. It then joins this data with
#' the SNR mean data obtained from the \code{create_snr_only_table} function.
#' The resulting combined table is returned.
#'
#' @param path The path to the CSV file.
#' @return A data frame representing the combined table for SNR plots.
#' @export
create_table_for_snr_plots <- function(path) {
  mae_data <- read_csv(path) %>%
    group_by(Feature_name) %>%
    summarise(mae = mean(mae)) %>%
    group_by(mae)
  snr_mean <- create_snr_only_table(path)
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

#' Save rank table
#'
#' This function takes a data frame and an output path and saves a rank table
#' based on the given data. It filters out the "shape" feature class, calculates
#' the rank based on the "mae.y" column, and groups the data by feature class,
#' "mae.y", rank, and Feature_name. The resulting table is saved as a CSV file
#' at the specified output path.
#'
#' @param data The data frame to create the rank table from.
#' @param output_path The path to save the rank table CSV file.
#' @export
save_rank_table <- function(data, output_path) {
  data %>%
    filter(feature_class != "shape") %>%
    mutate(rank = dense_rank((mae.y))) %>%
    group_by(feature_class, mae.y, rank, Feature_name) %>%
    write_csv(output_path)
}

#' Create plots
#'
#' This function takes a data frame and three output paths and creates plots
#' based on the given data. It uses the ggplot2 package to create scatter plots
#' and density ridges plots. The resulting plots are saved as PNG files at the
#' specified output paths.
#'
#' @param data The data frame to create the plots from.
#' @param output1 The path to save the first plot PNG file.
#' @param output2 The path to save the second plot PNG file.
#' @param output3 The path to save the third plot PNG file.
#' @export
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


plot_top_12_features(snakemake@input[[1]],   # input: features.csv (normalized)
                     snakemake@output[[1]])

plot_snr_curves(snakemake@input[[1]],
                snakemake@output[[1]] %>%
                  file.path %>%
                  dirname,
                palettes)

data_summary_plots <- create_table_for_snr_plots(snakemake@input[[2]]) # input mae.csv

save_rank_table(data_summary_plots, snakemake@output[[2]])

create_plots(data_summary_plots,
             snakemake@output[[2]],
             snakemake@output[[3]],
             snakemake@output[[4]])
