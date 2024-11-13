library(tidyverse)
library(ggridges)
library(munsell)

#' Function to prepare data for analysis
#'
#' This function takes a path as input and loads all CSV files in the folder specified by the path.
#' It then separates the file names into three columns: patient, noise, and seed.
#' The "seed" column is modified to remove the ".csv" extension.
#'
#' @param path A string specifying the path to the folder containing the CSV files.
#'
#' @return A data frame with the modified data
prepare_data <- function(path) {
  data <- read_csv(path) %>%
    separate(file, into = c("patient", "noise", "seed"),
             sep = "_", remove = FALSE) %>%
    mutate(seed = str_remove(seed, ".csv"))
  return(data)
}

#' Plot Top 12 Features
#'
#' This function plots the top 12 features based on their stability across different noise levels.
#'
#' @param path The path to the folder containing the CSV files
#' @param filename The name of the output file.
#'
#' @return None
#'
#' @examples
#' plot_top_12_features("/path/to/data", "output.png")
plot_top_12_features <- function(path, filename) {
  top_12_candidates <- prepare_data(path) %>%
    select(noise, `original_shape_Elongation`:`wavelet-LLL_ngtdm_Strength`) %>%
    pivot_longer(names_to = "feature", values_to = "value", -noise) %>%
    group_by(noise, feature) %>%
    summarize(value = mean(value)) %>%
    pivot_wider(names_from = noise, values_from = value) %>%
    mutate(distance = (`0.000` - `0.040`) / `0.000`) %>%
    arrange(-abs(distance)) %>%
    head(12) %>%
    pull(feature)
  prepare_data(path) %>%
    select(ID, noise, file, top_12_candidates) %>%
    pivot_longer(names_to = "feature",
                 values_to = "value",
                 top_12_candidates) %>%
    ggplot(aes(x = ID, y = value, color = noise)) +
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

#' Plot Noise Curves
#'
#' This function plots noise curves for different features using ggplot2.
#'
#' @param path The path to the folder containing the CSV files
#' @param output_path The path to save the output plots.
#'
#' @return None
#'
#' @examples
#' plot_noise_curves("data", "output/plots")
plot_noise_curves <- function(path, output_path) {
  prepare_data(path) %>%
    select(ID, noise, file,
           `original_shape_Elongation`:`wavelet-LLL_ngtdm_Strength`) %>%
    pivot_longer(names_to = "feature",
                 values_to = "value",
                 `original_shape_Elongation`:`wavelet-LLL_ngtdm_Strength`) %>%
    group_split(feature) %>%
    map(~{
          plot <- ggplot(.x) + aes(x = ID, y = value, color = noise) +
            geom_line(aes(group = file)) +
            labs(y = first(.x$feature)) +
            scale_color_brewer(palette="Set1") +
            theme_classic()
          ggsave(paste0(first(.x$feature), ".png"),
                 plot,
                 path = output_path,
                 width = 10)})
}




plot_top_12_features(snakemake@input[[1]],
                     snakemake@output[[1]])  # input: features.csv (normalized)

plot_noise_curves(snakemake@input[[1]],
                  snakemake@output[[1]] %>%
                    file.path %>%
                    dirname)