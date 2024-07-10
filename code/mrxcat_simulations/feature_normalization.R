library(tidyverse)

#' Normalize feature dataframe
#'
#' This function normalizes a feature dataframe by subtracting the mean and dividing by the standard deviation.
#'
#' @param df_location The location of the feature dataframes.
#' @param save_location The location to save the normalized dataframe.
#' @param standard_df The name of the standard dataframe to use for normalization.
#' @param reg Regular expression to filter the files in the df_location directory.
#'
#' @return None
#'
#' @examples
#' normalize_feature_df("/path/to/feature_data", "/path/to/save_location", "standard_df", "")
#'
normalize_feature_df <- function(df_location, save_location, standard_df, reg) {
  files <- dir(df_location)
  names(files) <- files
  if (reg != "") {
    files <- grep(reg, files, value = TRUE)
  }
  data <- files %>%
    map_df(
      ~read_csv(paste0(str_glue("{df_location}/"), .x), show_col_types = FALSE),
      .id = "file"
    ) %>%
    pivot_longer(
      names_to = "feature",
      values_to = "value",
      `original_shape_Elongation`:`wavelet-LLL_ngtdm_Strength`
    )
  data_left_join <- data %>%
    filter(file == standard_df) %>%
    group_by(feature) %>%
    summarize(mean_value = mean(value), sd_value = sd(value))
  data_long <- left_join(data, data_left_join, by = c("feature")) %>%
    mutate(value = (value - mean_value) / sd_value, .keep = "unused") %>%
    pivot_wider(names_from = feature, values_from = value)
  data_long %>%
    write_csv(file.path(save_location))
}

normalize_feature_df(
  snakemake@input[[1]] %>%
    file.path %>%
    dirname,
  snakemake@output[[1]],
  snakemake@params[[1]],
  snakemake@params[[2]]
)
