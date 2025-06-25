library(dplyr)
library(ggplot2)

setwd("C:/Users/david9456/Desktop/env_vars_Global")
data_dir <- getwd()
file_list <- list.files(path = data_dir, pattern = "\\.txt$", full.names = TRUE)

interp <- function(file) {
  df <- read.table(file, header = TRUE, stringsAsFactors = FALSE)
  age <- colnames(df)[1]
  value_col <- colnames(df)[2]
  
  time_intervals <- diff(df[[age]])
  avg_time_interval <- mean(time_intervals, na.rm = TRUE)
  
  if (avg_time_interval <= 1) {
    cat("Using linear interpolation for", basename(file), ", average time interval <= 1.\n")
    linear_interp <- approx(df[[age]], df[[value_col]], xout = seq(ceiling(min(df[[age]])), floor(max(df[[age]])), by = 1))
    df_new <- data.frame(age = linear_interp$x)
    df_new[[value_col]] <- linear_interp$y
  } else {
    cat("Using spline interpolation for", basename(file), ", average time interval > 1.\n")
    new_time <- seq(ceiling(min(df[[age]])), floor(max(df[[age]])), by = 1)
    spline_interp <- spline(df[[age]], df[[value_col]], xout = new_time)
    df_new <- data.frame(age = spline_interp$x)
    df_new[[value_col]] <- spline_interp$y
  }
  
  new_file_name <- paste0(tools::file_path_sans_ext(file), "_interpolated.txt")
  write.table(df_new, file = new_file_name, row.names = FALSE, quote = FALSE, sep = "\t")
  
  p <- ggplot() +
    geom_point(data = df, aes(x = .data[[age]], y = .data[[value_col]]), color = "red", size = 2, alpha = 0.7) +
    geom_line(data = df_new, aes(x = age, y = .data[[value_col]]), color = "blue", linewidth = 1) +
    labs(title = paste(basename(file)),
         x = "Age", y = value_col) +
    theme_minimal()
  
  print(p)
}

lapply(file_list, interp)
