library(tidyverse)
library(ggridges)
library(munsell)

# This section defines color palettes used for plotting in the script.
# Each palette is stored as a vector of hexadecimal color codes.
# The class_colors vector stores additional color codes used for class-specific coloring.


s1_palette <- c("#818589", "#CC29C7", "#993D96", "#FF01F6", "#332933")
s2_palette <- c("#818589", "#B86725", "#855935", "#EB6800", "#332D29")
s3_palette <- c("#818589", "#26BDBD", "#378A8A", "#00F0EE", "#293333")
s4_palette <- c("#818589", "#2BBD26", "#3A8A37", "#09F000", "#293329")
s5_palette <- c("#818589", "#B8B025", "#858135", "#EBE100", "#333229")
s6_palette <- c("#818589", "#B82A25", "#853835", "#EB0800", "#332929")
s7_palette <- c("#818589", "#f2f0f7", "#cbc9e2", "#9e9ac8", "#6a51a3")
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
#' @param palettes A list of color palettes for different feature groups.
#'
#' @return None
#'
#' @examples
#' plot_noise_curves("data", "output/plots", palettes)
plot_noise_curves <- function(path, output_path, palettes) {
  prepare_data(path) %>%
    select(ID, noise, file,
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
          plot <- ggplot(.x) + aes(x = ID, y = value, color = noise) +
            geom_line(aes(group = file)) +
            labs(y = first(.x$feature)) +
            scale_color_manual(values = color_palette) +
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
                    dirname,
                  palettes)