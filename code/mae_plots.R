library(tidyverse)
library(ggridges)
library(munsell)

#' Function to load CSV files in a folder and combine them into a single data frame
#'
#' @param path The path to the folder containing the CSV files
#'
#' @return A data frame containing the combined data from all the CSV files
load_csvs_in_folder <- function(path) {
  files <- dir(path)
  names(files) <- files
  data <- files %>%
    map_df(
      ~read_csv(paste0(str_glue("{path}/"), .x), show_col_types = FALSE),
      .id = "source"
    )
  return(data)
}


#' Create noise Only Table
#'
#' This function reads a CSV file from the specified path and creates a table
#' that summarizes the mean absolute error (MAE) values for different features
#' at different noise levels.
#'
#' @param path The path to the CSV file.
#'
#' @return A data frame containing the summarized data.
#'
#' @examples
#' create_noise_only_table("/path/to/data.csv")
#'
#' @param path The path to the folder containing the CSV files.
#' @return A data frame representing the noise mean table.


create_noise_only_table <- function(path) {
  data <- load_csvs_in_folder(path) %>%
    separate(Dataframe_1,
             into = c("patient_a", "noise_a", "seed_a"),
             sep = "_", remove = FALSE) %>%
    mutate(seed_a = str_remove(seed_a, ".csv")) #%>%
    # mutate(patient_a = str_remove(patient_a, "patient"))
  data <- data %>%
    separate(Dataframe_2,
             into = c("patient_b", "noise_b", "seed_b"),
             sep = "_",
             remove = FALSE) %>%
    mutate(seed_b = str_remove(seed_b, ".csv")) #%>%
    # mutate(patient_b = str_remove(patient_b, "patient"))
  # Define the noise values
  noise_values <- c("0.000", "0.010", "0.020", "0.030", "0.040")
  # Create an empty list to store the data frames
  noise_data <- list()
  # Iterate over the noise values
  for (x in seq(from = 1, to = 5, by = 1)) {
    # Filter the data based on noise values
    filtered_data <- filter(data, noise_a == noise_values[[x]], noise_b == noise_values[[x]])
    # Group and summarize the filtered data
    summarized_data <- filtered_data %>%
      group_by(Feature_name) %>%
      summarise(mae = mean(mae)) %>%
      group_by(mae)
    # Add the noise column to the summarized data
    summarized_data$noise <- as.character(noise_values[[x]])
    # Store the summarized data in the list
    noise_data[[x]] <- summarized_data
  }
   noise_mean <- bind_rows(noise_data)

  return(noise_mean)
}

#' Creates a table for noise plots
#'
#' This function takes a path as input and loads CSV files in the specified folder.
#' It groups the data by the feature name and calculates the mean absolute error (MAE) for each group.
#' The resulting table is then further processed by separating the feature name into filter, feature class, and feature columns.
#' The function also creates a separate table for noise data by calling the create_noise_only_table function and modifies the
#' table the same way as the previous table. Finally, it joins the noise table and the MAE table based on the feature name,
#' filter, feature class, and feature columns.
#'
#' @param path The path to the folder containing the CSV files
#' @return A combined table with noise specific MAE data and general MAE data
#' @export
create_table_for_noise_plots <- function(path) {
  mae_data <- load_csvs_in_folder(path) %>%
    group_by(Feature_name) %>%
    summarise(mae = mean(mae)) %>%
    group_by(mae)
  noise_mean <- create_noise_only_table(path)
  mae_data <- mae_data %>%
    separate(Feature_name,
             into = c("filter", "feature_class", "feature"),
             sep = "_",
             remove = FALSE)
  noise_mean <- noise_mean %>%
    separate(Feature_name,
             into = c("filter", "feature_class", "feature"),
             sep = "_",
             remove = FALSE)
  combined_table <- left_join(noise_mean,
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
    group_by(filter, feature_class, mae.y, rank, Feature_name) %>%
    distinct(Feature_name, mae.y, rank) %>%
    write_csv(output_path)
}

create_plots <- function(data, output1, output2, output3) {
  data %>%
    slice_min(mae.y, prop = .90) %>%
    ggplot(aes(x = mae.y, y = `mae.x`, group = noise, color = noise)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    xlab("Mean average error (MAE) per noise class") +
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
    xlab("Mean average error (MAE) per noise class") +
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


data_summary_plots <- create_table_for_noise_plots(snakemake@input[[1]] %>%
                                                     file.path %>%
                                                     dirname) # input mae.csv

save_rank_table(data_summary_plots, snakemake@output[[1]])

create_plots(data_summary_plots,
             snakemake@output[[2]],
             snakemake@output[[3]],
             snakemake@output[[4]])
