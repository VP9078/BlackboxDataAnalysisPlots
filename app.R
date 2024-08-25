library(shiny)
library(shinyjs)
library(patchwork)
library(knitr)
library(colourpicker)
library(randomcoloR)
library(ggplot2)
library(rhdf5)
library(stringr)
library(dplyr)
library(shinyWidgets)

##########################
#### USEFUL FUNCTIONS ####
##########################

# Increase the maximum upload size to 500 MB
options(shiny.maxRequestSize = 500 * 1024^2)
# Create a new environment to store plots
plot_storage <- new.env()
# Conversion factor from pixels to cm
pixel_to_cm <- 0.1

remove_periods_and_underscores <- function(text) {
  text <- gsub("_", " ", text)  # Replace underscores with spaces
  text <- gsub("\\.{2,}", ".", text)  # Replace multiple periods with a single period
  text <- gsub("\\.", " ", text)  # Replace periods with spaces
  text <- gsub(" +", " ", text)  # Replace multiple spaces with a single space
  text <- trimws(text)  # Trim leading and trailing whitespace
  
  if (grepl("average", text)) {
    avgidx <- grep("average", text) + nchar("average") + 1
    text <- str_sub(text, start = (avgidx))
  }
  
  return(text)
}

cleanup <- function(x) {
  str_to_title(remove_periods_and_underscores(x))
}

rid0s <- function(df, vec) {
  # Get the column data as a vector
  clmn <- df[[vec]]
  
  # Get the last 3 elements of the column
  if (length(clmn) >= 3) {
    last_3 <- tail(clmn, 3)
  } else {
    last_3 <- clmn
  }
  
  # Check if the last 3 elements are all 0s
  if (all(last_3 == 0)) {
    # Remove rows with 0s at the end
    while(length(clmn) > 0 && clmn[length(clmn)] == 0) {
      df <- df[-nrow(df),]
      clmn <- df[[vec]]
    }
  }
  
  return(df)
}

get_mouse_name <- function(file) {
  if(grepl("features.h5", file))
    return(str_extract(file, "(?<=_)[^_]+(?=-features\\.h5)"))
  if(grepl("tracking.h5", file))
    return(str_extract(file, "(?<=_)[^_]+(?=-tracking\\.h5)"))
}

#######################
#### TRACKING CODE ####
#######################

tracking_load <- function(file_path, tag) {
  
  suppressWarnings(assign(paste0(tag, "_tracking_data"), H5Fopen(file_path), envir = .GlobalEnv))
  tracking_data <- get(paste0(tag, "_tracking_data"))
  
  suppressWarnings(assign(paste0(tag, "_tracking_coordinates"), as.data.frame(h5read(tracking_data, name = "/df_with_missing/table", compoundAsDataFrame = FALSE)$values_block_0), envir = .GlobalEnv))
  tracking_coordinates <- get(paste0(tag, "_tracking_coordinates"))
  
  tracking_names <- h5readAttributes(tracking_data, name = "/df_with_missing/table")$values_block_0_kind
  pattern <- "(?<=\\nV)(.*?)(?=\\np)"
  all_matches <- str_extract_all(tracking_names, pattern)[[1]]
  exclude_indices <- c(2, 3, 4)
  body_parts <<- all_matches[-exclude_indices]
  
  assign(paste0(tag, "_row_names"),
         unlist(lapply(body_parts, function(part) {
           c(paste(part, "x", sep = " "),
             paste(part, "y", sep = " "),
             paste(part, "likelihood", sep = " "))
         }), use.names = FALSE), envir = .GlobalEnv)
  
  rownames(tracking_coordinates) <- get(paste0(tag, "_row_names"))
  assign(paste0(tag, "_tracking_coordinates"), tracking_coordinates, envir = .GlobalEnv)
  
  rm(tracking_coordinates, tracking_data)
}

plot_speed <- function(tag, body_part, start_sec = NULL, end_sec = NULL) {
  tracking_coordinates <- get(paste0(tag, "_tracking_coordinates"))
  fps <- 45  # Frames per second
  
  start_idx <- ifelse(is.null(start_sec), 1, start_sec * fps)
  end_idx <- ifelse(is.null(end_sec) || end_sec == 0, ncol(tracking_coordinates), end_sec * fps)
  
  subset_size <- start_idx:end_idx
  
  distances <- numeric()
  speeds <- numeric()
  
  time_between_frames <- 1 / fps
  
  x_row_index <- which(grepl("x", get(paste0(tag, "_row_names"))) & grepl(body_part, get(paste0(tag, "_row_names"))))
  x_row <- tracking_coordinates[x_row_index, ]
  x_row <- unlist(x_row)
  
  y_row_index <- which(grepl("y", get(paste0(tag, "_row_names"))) & grepl(body_part, get(paste0(tag, "_row_names"))))
  y_row <- tracking_coordinates[y_row_index, ]
  y_row <- unlist(y_row)
  
  for (i in 2:length(subset_size)) {
    distances[i - 1] <- sqrt((x_row[subset_size[i]] - x_row[subset_size[i - 1]])^2 + (y_row[subset_size[i]] - y_row[subset_size[i - 1]])^2)
  }
  distances <- distances * pixel_to_cm
  
  for (i in seq_along(distances)) {
    speeds[i] <- distances[i] / time_between_frames
  }
  assign(paste0(tag, "_", body_part, "_speed_values"), speeds, envir = .GlobalEnv)
  
  time_vector <- seq(from = start_sec + time_between_frames, to = (start_sec + (length(subset_size) - 1) * time_between_frames), by = time_between_frames)
  
  # Ensure vectors are of the same length
  if (length(time_vector) != length(speeds)) {
    min_length <- min(length(time_vector), length(speeds))
    time_vector <- time_vector[1:min_length]
    speeds <- speeds[1:min_length]
  }
  
  speed_df <- data.frame(time_vector, speeds)
  
  # Get rid of 0s at the end
  speed_df <- rid0s(df = speed_df, vec = "speeds")
  
  # Get rid of NAs
  speed_df <- na.omit(speed_df)
  
  average_speed <<- mean(speed_df$speeds)
  
  assign(paste(tag, tracking_mouse_name, body_part, "Speed_df", sep = "_"), speed_df, envir = .GlobalEnv)
  
   return(
     ggplot(speed_df, aes(x = time_vector, y = speeds)) +
     geom_line(color = "blue") +
     labs(
       title = str_to_title(paste0(tag, " ", tracking_mouse_name, " ", remove_periods_and_underscores(body_part), " Speed")),
       x = "Time (s)",
       y = "Speed"
     ) +
     theme_minimal() +
     theme(
       plot.title = element_text(face = "bold", hjust = 0.5, size = 18)
     ))
}

plot_likelihood <- function(tag, body_part, start_sec = NULL, end_sec = NULL) {
  tracking_coordinates <- get(paste0(tag, "_tracking_coordinates"))
  fps <- 45  # Frames per second
  
  start_idx <- ifelse(is.null(start_sec), 1, start_sec * fps)
  end_idx <- ifelse(is.null(end_sec) || end_sec == 0, ncol(tracking_coordinates), end_sec * fps)
  
  subset_size <- start_idx:end_idx
  
  time_between_frames <- 1 / fps
  
  likelihood_row_index <- which(grepl("likelihood", get(paste0(tag, "_row_names"))) & grepl(body_part, get(paste0(tag, "_row_names"))))
  likelihood_row <- tracking_coordinates[likelihood_row_index, ]
  likelihood_row <- unlist(likelihood_row)
  
  time_vector <- seq(from = start_sec + time_between_frames, to = (start_sec + (length(subset_size) - 1) * time_between_frames), by = time_between_frames)
  
  # Ensure vectors are of the same length
  if (length(time_vector) != length(likelihood_row)) {
    min_length <- min(length(time_vector), length(likelihood_row))
    time_vector <- time_vector[1:min_length]
    likelihood_row <- likelihood_row[1:min_length]
  }
  
  likelihood_df <- data.frame(time_vector, likelihood_row)
  
  # Get rid of 0s at the end
  likelihood_df <- rid0s(df = likelihood_df, vec = "likelihood_row")
  
  average_likelihood <<- mean(likelihood_df$likelihood_row)
  
  # Get rid of NAs
  likelihood_df <- na.omit(likelihood_df)
  
  assign(paste(tag, tracking_mouse_name,body_part, "Likelihood_df", sep = "_"), likelihood_df, envir = .GlobalEnv)
  
  return(ggplot(likelihood_df, aes(x = time_vector, y = likelihood_row)) +
           geom_line(color = "green") +
           labs(
             title = str_to_title(paste0(tag, " ", tracking_mouse_name, " ", remove_periods_and_underscores(body_part), " Likelihood")),
             x = "Time (s)",
             y = "Likelihood"
           ) +
           theme_minimal() +
           theme(
             plot.title = element_text(face = "bold", hjust = 0.5, size = 18)
           ) +
           coord_cartesian(ylim = c(0, max(likelihood_row))))
}

plot_trajectory <- function(tag, body_part, start_sec = NULL, end_sec = NULL, grid_size = "full") {
  tracking_coordinates <- get(paste0(tag, "_tracking_coordinates"))
  fps <- 45  # Frames per second
  
  start_idx <- ifelse(is.null(start_sec), 1, start_sec * fps)
  end_idx <- ifelse(is.null(end_sec) || end_sec == 0, ncol(tracking_coordinates), end_sec * fps)
  
  subset_size <- start_idx:end_idx
  
  x_row_index <- which(grepl("x", get(paste0(tag, "_row_names"))) & grepl(body_part, get(paste0(tag, "_row_names"))))
  x_row <- tracking_coordinates[x_row_index, ]
  x_row <- unlist(x_row) * pixel_to_cm
  
  y_row_index <- which(grepl("y", get(paste0(tag, "_row_names"))) & grepl(body_part, get(paste0(tag, "_row_names"))))
  y_row <- tracking_coordinates[y_row_index, ]
  y_row <- unlist(y_row) * pixel_to_cm
  
  trajectory_df <- data.frame(x_row, y_row)[subset_size, ]
  trajectory_df$index <- 1:nrow(trajectory_df)
  
  if (grid_size == "full")
    grid_size <- max(x_row)
  
  # Get rid of 0s at the end
  trajectory_df <- rid0s(df = trajectory_df, vec = "x_row")
  
  assign(paste(tag, tracking_mouse_name, body_part, "Trajectory_df", sep = "_"), trajectory_df, envir = .GlobalEnv)

  return(trajectory_plot <- ggplot(trajectory_df, aes(x = x_row, y = y_row)) +
    geom_path(aes(color = index)) +
    scale_color_gradientn(colors = c("blue", "green", "yellow", "red")) +
    labs(
      title = str_to_title(paste0(tag, " ", tracking_mouse_name, " ", remove_periods_and_underscores(body_part), " Trajectory")),
      x = "Box Length (cm)",
      y = "Box Width (cm)"
    ) +
    theme_minimal() +
    coord_cartesian(xlim = c(0, grid_size), ylim = c(0, grid_size)) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18)
    ))
}

plot_heatmap <- function(tag, body_part, start_sec = NULL, end_sec = NULL, grid_size = "full") {
  tracking_coordinates <- get(paste0(tag, "_tracking_coordinates"))
  fps <- 45  # Frames per second
  
  start_idx <- ifelse(is.null(start_sec), 1, start_sec * fps)
  end_idx <- ifelse(is.null(end_sec) || end_sec == 0, ncol(tracking_coordinates), end_sec * fps)
  
  subset_size <- start_idx:end_idx
  
  x_row_index <- which(grepl("x", get(paste0(tag, "_row_names"))) & grepl(body_part, get(paste0(tag, "_row_names"))))
  x_row <- tracking_coordinates[x_row_index, ]
  x_row <- unlist(x_row) * pixel_to_cm
  
  y_row_index <- which(grepl("y", get(paste0(tag, "_row_names"))) & grepl(body_part, get(paste0(tag, "_row_names"))))
  y_row <- tracking_coordinates[y_row_index, ]
  y_row <- unlist(y_row) * pixel_to_cm
  
  if (grid_size == "full")
    grid_size <- 50
  
  trajectory_df <- data.frame(x_row, y_row)[subset_size, ]
  
  # Get rid of 0s at the end
  trajectory_df <- rid0s(df = trajectory_df, vec = "x_row")
  
  assign(paste(tag, tracking_mouse_name, body_part, "Trajectory_df", sep = "_"), trajectory_df, envir = .GlobalEnv)

  trajectory_df <- as.matrix(trajectory_df)
  
  return(heatmap <- ggplot(trajectory_df, aes(x = x_row, y = y_row)) +
    geom_bin_2d(bins = 50) +
    scale_fill_gradient(low = "blue", high = "red", name = "Density") +
    theme_minimal() +
    coord_cartesian(xlim = c(0, grid_size), ylim = c(0, grid_size)) +
    labs(title = paste0(tag, " ", tracking_mouse_name, " Cage Occupancy Heatmap"),
         x = "Box Length (cm)",
         y = "Box Width (cm)") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18)
    ))
}

### Batch extract speed averages ###
load_multiple_tracking_files <- function(file_paths, file_names) {
  process_tracking_file <- function(file_path, file_name) {
    
    tracking_data <- H5Fopen(file_path)
    tracking_coordinates <- as.data.frame(h5read(tracking_data, name = "/df_with_missing/table", compoundAsDataFrame = FALSE)$values_block_0)
    tracking_names <- h5readAttributes(tracking_data, name = "/df_with_missing/table")$values_block_0_kind
    pattern <- "(?<=\\nV)(.*?)(?=\\np)"
    all_matches <- str_extract_all(tracking_names, pattern)[[1]]
    exclude_indices <- c(2, 3, 4)
    batch_body_parts <<- all_matches[-exclude_indices]
    assign(paste0(file_name, "_row_names"),
           unlist(lapply(batch_body_parts, function(part) {
             c(paste(part, "x", sep = " "),
               paste(part, "y", sep = " "),
               paste(part, "likelihood", sep = " "))
           }), use.names = FALSE), envir = .GlobalEnv)
    
    rownames(tracking_coordinates) <- get(paste0(file_name, "_row_names"))
    assign(paste0(file_name, "_tracking_coordinates"), tracking_coordinates, envir = .GlobalEnv)
    rm(tracking_coordinates, tracking_data)
  }
  for (i in seq_along(file_paths)) {
    process_tracking_file(file_paths[i], file_names[i])
  }
}

compute_speed_avgs <- function(file_names, body_part, start_sec = NULL, end_sec = NULL) {
  get_avg <- function(file_name) {
    tracking_coordinates <- get(paste0(file_name, "_tracking_coordinates"))
    fps <- 45
    start_idx <- ifelse(is.null(start_sec), 1, start_sec * fps)
    end_idx <- ifelse(is.null(end_sec) || end_sec == 0, ncol(tracking_coordinates), end_sec * fps)
    subset_size <- start_idx:end_idx
    distances <- numeric()
    speeds <- numeric()
    time_between_frames <- 1 / fps
    x_row_index <<- which(grepl("x", get(paste0(file_name, "_row_names"))) & grepl(body_part, get(paste0(file_name, "_row_names"))))
    x_row <- tracking_coordinates[x_row_index, ]
    x_row <- unlist(x_row)
    y_row_index <<- which(grepl("y", get(paste0(file_name, "_row_names"))) & grepl(body_part, get(paste0(file_name, "_row_names"))))
    y_row <- tracking_coordinates[y_row_index, ]
    y_row <- unlist(y_row)
    for (i in 2:length(subset_size)) {
      distances[i - 1] <- sqrt((x_row[subset_size[i]] - x_row[subset_size[i - 1]])^2 + (y_row[subset_size[i]] - y_row[subset_size[i - 1]])^2)
    }
    distances <- distances * pixel_to_cm
    for (i in seq_along(distances)) {
      speeds[i] <- distances[i] / time_between_frames
    }
    time_vector <- seq(from = start_sec + time_between_frames, to = (start_sec + (length(subset_size) - 1) * time_between_frames), by = time_between_frames)
    if (length(time_vector) != length(speeds)) {
      min_length <- min(length(time_vector), length(speeds))
      time_vector <- time_vector[1:min_length]
      speeds <- speeds[1:min_length]
    }
    speed_df <- data.frame(time_vector, speeds)
    speed_df <- rid0s(df = speed_df, vec = "speeds")
    speed_df <- na.omit(speed_df)
    average_speed <- mean(speed_df$speeds)
    return(average_speed)
  }
  avgs <- list()
  for (file in file_names) {
    avgs[[file]] <- get_avg(file)
  }
  avgs_df <- data.frame(Averages <- unlist(avgs))
  
  colnames(avgs_df) <- c("average_speeds")
  
  return(avgs_df)
}

#######################
#### FEATURES CODE ####
#######################

features_load_and_plot <- function(bboxdataPath, bboxdsetnm, metric, start_sec = NULL, end_sec = NULL) {
  suppressWarnings(bboxdata <<- H5Fopen(bboxdataPath))
  groupls <- h5ls(bboxdata)
  bboxgroupname <- unique(groupls$group)[2]
  
  behavioral_metrics <<- names(h5read(bboxdataPath, bboxgroupname))
  filtered_metrics <- behavioral_metrics[!grepl("fps|frame_count|recording_time", behavioral_metrics)]
  
  objectpath <- paste(bboxgroupname, metric, sep="/")
  assign(paste0(metric,"_", bboxdsetnm), h5read(bboxdataPath, objectpath))
  
  obj <- get(paste0(metric,"_", bboxdsetnm))
  
  fps <- 45  # Frames per second
  start_idx <- ifelse(is.null(start_sec), 1, start_sec * fps)
  end_idx <- ifelse(is.null(end_sec) || end_sec == 0, length(obj), end_sec * fps)
  
  obj <- obj[start_idx:end_idx]
  objframe_count <- length(obj)
  vid_length <- objframe_count / fps
  time_increments <<- 1 / fps
  timelist <- seq(from = start_sec + time_increments, to = (start_sec + (length(obj) - 1) * time_increments), by = time_increments)
  
  obj_df <- data.frame(values = obj, times = timelist[1:objframe_count])
  
  if (any(is.na(obj_df))){
    obj_df <- na.omit(obj_df)
  }
  
  # Ensure vectors are of the same length
  if (length(obj_df$times) != length(obj_df$values)) {
    min_length <- min(length(obj_df$times), length(obj_df$values))
    obj_df <- obj_df[1:min_length, ]
  }
  
  df_name <- paste(metric, features_mouse_name, bboxdsetnm, "features_df", sep = "_")
  assign(df_name, obj_df, envir = .GlobalEnv)

  if (is.factor(obj) || is.logical(obj)) {

    features_true_tot_pct <<- (length(which(obj_df$values == T)) / length(obj)) * 100
    
    obj_df$values <- as.numeric(obj_df$values)
    return(ggplot(obj_df, aes(x = times, y = values)) +
      geom_tile(aes(width = time_increments, height = Inf, fill = factor(values))) +
      scale_fill_manual(values = c("1" = "transparent", "0" = "red"), na.value = "purple") +
      labs(x = "Time (s)", y = str_to_title(remove_periods_and_underscores(metric)), fill = "State") +
      ggtitle(str_to_title(paste0(remove_periods_and_underscores(metric), " ", remove_periods_and_underscores(bboxdsetnm), " ", features_mouse_name, " Time Series Plot"))) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18)))
  }
  else if (is.array(obj)) {
    
    # Get rid of 0s at the end
    obj_df <- rid0s(df = obj_df, vec = "values")
    
    features_average <<- mean(obj_df$values)
    
    return(ggplot(obj_df, aes(x = times, y = values)) +
      geom_line(color = "red") +
      labs(x = "Time (s)",
           y = str_to_title(remove_periods_and_underscores(metric))) +
      ggtitle(str_to_title(paste0(remove_periods_and_underscores(metric), " ", remove_periods_and_underscores(bboxdsetnm), " ", features_mouse_name, " Time Series Plot"))) +
      theme_minimal() +
        theme(
          
        plot.title = element_text(face = "bold", hjust = 0.5, size = 18)
      ))
  }
}

# Batch extract features averages
compute_features_avgs <- function(file_paths, file_names, features_metric, start_sec = NULL, end_sec = NULL) {
  
  get_avg <- function(file_path, file_name) {
    suppressWarnings(bboxdata <<- H5Fopen(file_path))
    groupls <- h5ls(bboxdata)
    bboxgroupname <- unique(groupls$group)[2]
    
    objectpath <- paste(bboxgroupname, features_metric, sep="/")
    assign(paste0(features_metric,"_", file_name), h5read(file_path, objectpath))
    
    obj <- get(paste0(features_metric,"_", file_name))
    
    fps <- 45  # Frames per second
    start_idx <- ifelse(is.null(start_sec), 1, start_sec * fps)
    end_idx <- ifelse(is.null(end_sec) || end_sec == 0, length(obj), end_sec * fps)
    
    obj <- obj[start_idx:end_idx]
    objframe_count <- length(obj)
    time_increments <- 1 / fps
    timelist <- seq(from = start_sec + time_increments, to = (start_sec + (length(obj) - 1) * time_increments), by = time_increments)
    
    obj_df <- data.frame(values = obj, times = timelist[1:objframe_count])
    
    if (any(is.na(obj_df))) {
      obj_df <- na.omit(obj_df)
    }
    
    # Ensure vectors are of the same length
    if (length(obj_df$times) != length(obj_df$values)) {
      min_length <- min(length(obj_df$times), length(obj_df$values))
      obj_df <- obj_df[1:min_length, ]
    }
    
    if (is.factor(obj_df$values) || is.logical(obj_df$values)) {
      obj_factor <- TRUE
      result <- (length(which(obj_df$values == TRUE)) / length(obj)) * 100
    } else if (is.array(obj) || is.numeric(obj)) {
      obj_factor <- FALSE
      result <- mean(obj_df$values)
    }
    
    return(list(result = result, obj_factor = obj_factor))
  }
  
  avgs <- numeric(length(file_names))
  obj_factors <- logical(length(file_names))
  
  for (i in seq_along(file_names)) {
    result <- get_avg(file_path = file_paths[i], file_name = file_names[i])
    avgs[i] <- result$result
    obj_factors[i] <- result$obj_factor
  }
  
  avgs_df <- data.frame(avgs)
  rownames(avgs_df) <- file_names
  
  if (any(obj_factors)) {
    colnames(avgs_df) <- paste0("percentages_true_", features_metric)
  } else {
    colnames(avgs_df) <- paste0("average_", features_metric)
  }
  
  return(avgs_df)
}

#######################
#### SUMMARY CODE #####
#######################

summary_load_and_plot <- function(tag, file_path, mouse_name, metric, y_limit = NULL) {
  
  suppressWarnings({
    assign(paste0(tag, "_summary_data"), read.csv(file_path, row.names = 1, header = TRUE, fileEncoding = "UTF-8"), envir = .GlobalEnv)
  })
  
  metrics <<- colnames(get(paste0(tag, "_summary_data")))
  mouse_names <<- rownames(get(paste0(tag, "_summary_data")))
  
  summary_data <- get(paste0(tag, "_summary_data"))
  
  # Validate if the selected mouse name and metric exist in the data
  validate(
    need(mouse_name %in% rownames(summary_data), "Selected mouse name does not exist in the data."),
    need(metric %in% colnames(summary_data), "Selected metric does not exist in the data.")
  )
  
  row_index <- which(rownames(summary_data) == mouse_name)
  col_index <- which(colnames(summary_data) == metric)
  
  datapt <- data.frame(Value = summary_data[row_index, col_index])
  datapt$Treat_group <- rownames(summary_data)[row_index]
  datapt$Color <- ifelse(datapt$Value >= 0, "positive", "negative")
  colors <- c("positive" = "#00BFC4", "negative" = "#F8766D")
  
  plot <- ggplot(data = datapt, aes(x = Treat_group, y = Value, fill = Color)) +
    geom_bar(stat = "identity", width = 0.7, color = "black", size = 0.3) +
    scale_fill_manual(values = colors) +
    geom_text(aes(label = round(Value, 2)), vjust = ifelse(datapt$Value >= 0, -0.5, 1.5), color = "black", size = 5) +  # Adjusted text size
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    labs(
      title = str_to_title(paste0("Average ", tag, " ", mouse_name, " ", remove_periods_and_underscores(metric), " Bar Plot")),
      y = str_to_title(remove_periods_and_underscores(metric)),
      x = ""
    ) +
    theme_minimal(base_size = 15) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.title.y = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank()
    ) +
    coord_cartesian(ylim = y_limit)
  
  return(plot)
}

#######################
####### PCA CODE ######
#######################

PCA_2D <- function(data, features, colors) {
  # Perform PCA
  pca <- prcomp(data[, features], scale = TRUE)
  pca_df <- as.data.frame(pca$x)
  pca_df$Sheet <- data$Sheet  # Add the sheet information for color coding
  
  # Create the 2D PCA plot with color coding
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Sheet)) +
    geom_point(size = 3) +
    scale_color_manual(values = colors) +
    labs(title = "PCA",
         x = "Principal Component 1",
         y = "Principal Component 2") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  
  list(plot = p, pca = pca)  # Return PCA plot and object for further analysis
}

# Function to create a scree plot
scree_plot <- function(pca) {
  pca.var <- pca$sdev^2
  pca.var.per <- round(pca.var / sum(pca.var) * 100, 1)
  
  # Create a data frame for plotting
  pca.var.per.df <- data.frame(
    Number = seq_along(pca.var.per),
    Variance = pca.var.per
  )
  
  # Create the scree plot
  ggplot(data = pca.var.per.df, aes(x = Number, y = Variance)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Variance), position = position_dodge(width = 0.9), vjust = -0.25) +
    labs(x = "Principal Component", y = "Percentage of Variance Explained") +
    ggtitle("Principal Component Scree Plot") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
}

# Function to get the top 5 contributing variables to PC1
get_top_5 <- function(pca) {
  loading_scores <- pca$rotation[,1]
  gene_scores <- abs(loading_scores)
  gene_score_ranked <- sort(gene_scores, decreasing = TRUE)
  top_5_variables <- names(gene_score_ranked[1:5])
  print("Top 5 contributing variables to PC1:")
  print(top_5_variables)
  print("Scores:")
  print(pca$rotation[top_5_variables,1])
}

############
#### UI ####
############

ui <- navbarPage(
  title = "Blackbox Data Analysis Plots",
  tabPanel("Documentation",
           fluidPage(
             titlePanel("Documentation"),
             div(style = "text-align: center; margin-bottom: 20px;",
                 img(src = 'Blackbox_data_analysis_app_pipeline.png', 
                     style = "max-width: 75%; height: auto; margin-bottom: 10px;")),
             p("Image created using BioRender.com.",
               style = "font-size: 14px; font-style: italic; text-align: center;"),
             
             h3("Overview", style = "color: #2C3E50;"),
             p("This Shiny application allows for the analysis and visualization of Blackbox output data from tracking, features, and summary datasets.",
               style = "margin-bottom: 10px;"),
             p("Also allows for simple usage of Principal Component Analysis (PCA) for dimensionality reduction.",
               style = "margin-bottom: 20px;"),
             p("Make sure to input the correct file type.",
               style = "font-weight: bold; color: #E74C3C;"),
             
             h3("Troubleshooting Tips", style = "color: #2C3E50;"),
             tags$ul(
               tags$li("Reload the website if it crashes or the screen becomes unresponsive."),
               tags$li(HTML("Raw scripts available at: <a href='https://github.com/VP9078/Blackbox-Plotting-in-R.git' target='_blank'>GitHub Repository</a>.")),
               tags$li(HTML("Contact <a href='mailto:vihaanpande08@gmail.com'>vihaanpande08@gmail.com</a> for any queries."))
             ),
             
             tags$hr(),
             
             h3("Tracking Data Analysis", style = "color: #2C3E50;"),
             tags$ol(
               tags$li("Upload tracking HDF5 (tracking.h5) file."),
               tags$li("(Optional) Enter a tag (identifier) for the dataset."),
               tags$li("Select the body part to analyze and the type of plot (Speed, Likelihood, Trajectory, Heatmap)."),
               tags$li("(Optional) Specify start and end indices for data subsetting."),
               tags$li("Press 'Plot Data' to generate the plot."),
               tags$li("Press 'Download Plot' or 'Download Data as CSV'.")
             ),
             h4("Batch Extract Speed Averages", style = "color: #2C3E50;"),
             tags$ol(
               tags$li("Upload as many tracking HDF5 (tracking.h5) files as desired."),
               tags$li("Select the body part to analyze."),
               tags$li("(Optional) Specify start and end indices for data subsetting."),
               tags$li("Press 'Extract Averages' to extract the speed averages for the specified data.")
             ),
             
             h3("Features Data Analysis", style = "color: #2C3E50;"),
             tags$ol(
               tags$li("Upload features HDF5 (features.h5) file."),
               tags$li("(Optional) Enter a tag (identifier) for the dataset."),
               tags$li("Select the metric to analyze."),
               tags$li("(Optional) Specify start and end indices for data subsetting."),
               tags$li("Press 'Plot Data' to generate the plot."),
               tags$li("Press 'Download Plot' or 'Download Data as CSV'.")
             ),
             p("NOTE: Most plots will contain numeric data; however some (ethogram metrics) will be in the form of boolean values, and thus plotted as raster plots."),
             
             h4("Batch Extract Feature Averages", style = "color: #2C3E50;"),
             tags$ol(
               tags$li("Upload as many features HDF5 (features.h5) files as desired."),
               tags$li("Select the metric to analyze."),
               tags$li("(Optional) Specify start and end indices for data subsetting."),
               tags$li("Press 'Extract Averages' to extract the feature averages for the specified data.")
             ),
             p("NOTE: As previously mentioned, some data will be in the form of boolean (TRUE or FALSE) values. In this case, the extract averages feature will extract the percentage of values equal to TRUE for that metric."),
             
             h3("Summary Data Analysis", style = "color: #2C3E50;"),
             tags$ol(
               tags$li("Upload CSV file containing summary data (Try to avoid editing the CSV as much as possible; can edit the mouse names if desired)."),
               tags$li("(Optional) Enter a tag (identifier) for the dataset."),
               tags$li("Select the metric and mouse name to analyze."),
               tags$li("Press 'Plot Data' to generate the plot."),
               tags$li("Press 'Download Plot' or 'Download Data as CSV'.")
             ),
             
             h3("Principal Component Analysis (PCA)", style = "color: #2C3E50;"),
             tags$ol(
               tags$li("Upload CSV file containing summary data."),
               tags$li("(REQUIRED) Enter a UNIQUE identifier for the dataset."),
               tags$li("Press 'Load Sheet'."),
               tags$li("Repeat steps 1 & 2 until all desired files have been uploaded (must be >=2 files)."),
               tags$li("Change colors for PCA plots if desired (can highlight color hex code and replace with color name (eg blue))"),
               tags$li("Select >=2 features from the 'Select Features for PCA' dropdown menu."),
               tags$li("Press 'Plot PCA' to generate the PCA plot."),
               tags$li("Press 'Plot Scree' to generate the Scree plot.")
             ),
             p("NOTE: If having trouble deciding on which features to select, can 1. slim down the options by only considering those which have relevance to your case / condition and/or 2. run many (or all) of the metrics, and view the Top 5 contributing variables to PC1 box (top 5 most significant features)."),
             p("NOTE: Ignore 'Error: subscript out of bounds' if using fewer than 5 features."),
             p(HTML("For more information on PCA, visit: <a href='https://youtu.be/HMOI_lkzW08?si=pvO3uilHs-c1G2fj' target='_blank'>PCA Video Tutorial</a>.")),
             
             h3("Side by Side Plots", style = "color: #2C3E50;"),
             tags$ol(
               tags$li("Select multiple saved plots from the 'Saved Plots' dropdown. ('Saved Plots' dropdown will only retain plots of the TRACKING, FEATURES, or SUMMARY analysis type."),
               tags$li("Press 'Side by Side Plots' to compare them side by side (plots of the same type will scale with each other).")
             ),
             
             h3("Overlay Plots", style = "color: #2C3E50;"),
             tags$ol(
               tags$li("Select multiple saved plots from the 'Saved Plots' dropdown. ('Saved Plots' dropdown will only retain plots of the TRACKING, FEATURES, or SUMMARY analysis type."),
               tags$li("Press 'Overlay Plots' to compare them overlayed on top of eachother, using a checkbox system to select or deselect which data to view after initial plotting (press boxes with data names)."),
             ),
             p("NOTE: Only SPEED, LIKELIHOOD (both speed and likelihood coming from tracking datasets), and FEATURES plots can be overlayed."),
             p("If interested in manually editing any plots, download the plot using the download plots feature, and load the plot into Adobe Illustrator, where it will be available as a vectorized image.")
           )
  ),
  tabPanel("Analysis",
           sidebarLayout(
             sidebarPanel(
               selectInput("analysis_type", "Select Analysis Type", choices = c("Tracking Data", "Features Data", "Summary Data", "PCA")),
               uiOutput("dynamic_ui"),
               uiOutput("loaded_files"),  # Display loaded CSV files
               selectInput("saved_plots", "Saved Plots", choices = NULL, multiple = TRUE),
               actionButton("side_by_side_plots", "Side by Side Plots"),
               actionButton("overlay_plots", "Overlay Plots"),
             ),
             mainPanel(
               plotOutput("plot"),
               conditionalPanel(
                 condition = "output.tracking_plot_generated",
                 verbatimTextOutput("tracking_average")
               ),
               conditionalPanel(
                 condition = "output.features_plot_generated",
                 verbatimTextOutput("features_average")
               ),
               conditionalPanel(
                 condition = "output.pca_plot_generated",
                 verbatimTextOutput("top_5_output")
               ),
               uiOutput("dynamic_plot_ui"),
               downloadButton("downloadPlot", "Download Plot"),
               downloadButton("downloadData", "Download Data as CSV")
             )
           )
  )
)

################
#### SERVER ####
################

server <- function(input, output, session) {
  combined_plot <- NULL
  sheets <- reactiveValues(data = list())  # Use reactiveValues to store sheets
  
###################
#### SERVER UI ####
###################
  
  observeEvent(input$analysis_type, {
    if (input$analysis_type == "Tracking Data") {
      output$dynamic_ui <- renderUI({
        tagList(
          fileInput("tracking_file", "Choose HDF5 File", accept = ".h5"),
          textInput("tracking_tag", "Enter Tag", value = ""),
          selectInput("tracking_body_part", "Select Body Part", choices = NULL),
          selectInput("tracking_plot_type", "Select Plot Type", choices = c("Speed", "Likelihood", "Trajectory", "Heatmap")),
          numericInput("start_sec", "Starting Time (seconds) (Leave as 0 for full values)", value = 0),
          numericInput("end_sec", "Ending Time (seconds) (Leave as 0 for full values)", value = 0),
          actionButton("tracking_plot_data", "Plot Data"),
          h4("Batch Extract Speed Averages"),
          fileInput("batch_speed_avg_files", "Choose HDF5 Files", multiple = T, accept = ".h5"),
          selectInput("batch_speed_avg_body_part", "Select Body Part", choices = NULL),
          numericInput("batch_speed_avg_start_sec", "Starting Time (seconds) (Leave as 0 for full values)", value = 0),
          numericInput("batch_speed_avg_end_sec", "Ending Time (seconds) (Leave as 0 for full values)", value = 0),
          downloadButton("downloadAvgSpeed", "Extract Averages")
        )
      })
    } else if (input$analysis_type == "Features Data") {
      output$dynamic_ui <- renderUI({
        tagList(
          fileInput("features_file", "Choose HDF5 File", accept = ".h5"),
          textInput("features_dataset_name", "Enter Tag", value = ""),
          selectInput("features_metric", "Select Metric", choices = NULL),
          numericInput("start_sec", "Starting Time (seconds) (Leave as 0 for full values)", value = 0),
          numericInput("end_sec", "Ending Time (seconds) (Leave as 0 for full values)", value = 0),
          actionButton("features_load_data", "Plot Data"),
          h4("Batch Extract Feature Averages"),
          fileInput("batch_features_avg_files", "Choose HDF5 Files", multiple = T, accept = ".h5"),
          selectInput("batch_features_avg_metric", "Select Metric", choices = NULL),
          numericInput("batch_features_avg_start_sec", "Starting Time (seconds) (Leave as 0 for full values)", value = 0),
          numericInput("batch_features_avg_end_sec", "Ending Time (seconds) (Leave as 0 for full values)", value = 0),
          downloadButton("downloadAvgFeature", "Extract Averages")
        )
      })
    } else if (input$analysis_type == "Summary Data") {
      output$dynamic_ui <- renderUI({
        tagList(
          fileInput("summary_file", "Choose CSV File", accept = ".csv"),
          textInput("summary_tag", "Enter Tag", value = "summary"),
          selectInput("summary_metric", "Select Metric", choices = NULL),
          selectInput("summary_mouse_name", "Select Mouse Name", choices = NULL),
          actionButton("summary_plot_data", "Plot Data")
        )
      })
    } else if (input$analysis_type == "PCA") {
      output$dynamic_ui <- renderUI({
        tagList(
          fileInput("sheet_file", "Choose CSV File", accept = ".csv"),
          textInput("sheet_tag", "Enter Tag", value = ""),
          uiOutput("colorPickers"),  # Add color pickers output
          actionButton("load_sheet", "Load Sheet"),
          h4("Loaded Sheets"),
          uiOutput("loaded_sheets"),
          selectInput("selected_features", "Select Features for PCA", choices = NULL, multiple = TRUE),
          actionButton("plot_pca", "Plot PCA"),
          actionButton("plot_scree", "Plot Scree")
        )
      })
    }
  })
  
  output$tracking_plot_generated <- reactive({
    tracking_plot_generated()
  })
  
  output$features_plot_generated <- reactive({
    features_plot_generated()
  })
  
  output$pca_plot_generated <- reactive({
    pca_plot_generated()
  })
  
  output$overlay_plot_generated <- reactive({
    overlay_plot_generated()
  })
  
  outputOptions(output, "tracking_plot_generated", suspendWhenHidden = FALSE)
  outputOptions(output, "features_plot_generated", suspendWhenHidden = FALSE)
  outputOptions(output, "pca_plot_generated", suspendWhenHidden = FALSE)
  
  
  observeEvent(input$analysis_type, {
    tracking_plot_generated(FALSE)
    features_plot_generated(FALSE)
    pca_plot_generated(FALSE)
  })
  
####################
##### PCA LOGIC ####
####################
  observeEvent(input$load_sheet, {
    req(input$sheet_file, input$sheet_tag)
    
    validate(
      need(grepl("\\.csv$", input$sheet_file$name), "Please upload a valid CSV file.")
    )
    
    # Load the selected sheet
    filepath <- input$sheet_file$datapath
    filename <- input$sheet_file$name
    sheet_tag <- input$sheet_tag
    
    if (sheet_tag == "" ) {
      sheet_tag <- filename  # Use filename as tag if no custom tag is provided
    }
    if (grepl(" ", sheet_tag)) {
      sheet_tag <- gsub(" ", "_", sheet_tag)
    }
    
    # Save the sheet in reactiveValues
    sheets$data[[sheet_tag]] <- read.csv(filepath)
    
    # Update the feature selection based on the loaded sheets
    combined_columns <- unique(unlist(lapply(sheets$data, colnames)))
    
    # Remove "X", "X.1", and "recording_time..min." if they exist
    combined_columns <- combined_columns[!combined_columns %in% c("X", "X.1", "recording_time..min.", "total.recording_time..min.", "PV..bin.start.end..min.", "bin.duration..min.")]
    
    # Apply cleanup
    cleaned_columns <- sapply(combined_columns, cleanup, USE.NAMES = FALSE)
    
    # Make cleaned & combined column df
    all_PCA_columns <<- data.frame(combined_columns, cleaned_columns)
    
    # Update the feature selection input
    updateSelectInput(session, "selected_features", choices = all_PCA_columns[,2])
    
    # Automatically select all loaded sheets
    selected_sheets <- names(sheets$data)
    
    output$colorPickers <- renderUI({
      lapply(selected_sheets, function(sheet) {
        colourInput(paste0("color_", sheet), paste("Choose color for", sheet), value = randomColor())
      })
    })
    
    # Update the list of loaded sheets with a selectizeInput
    output$loaded_sheets <- renderUI({
      selectizeInput(
        inputId = "loaded_sheets_select",
        label = "Loaded Sheets",
        choices = selected_sheets,
        selected = selected_sheets,  # Automatically select all loaded sheets
        multiple = TRUE,
        options = list(placeholder = 'Select loaded sheets...')
      )
    })
  })
  
  # Reactive value to track whether the plot has been generated
  pca_plot_generated <- reactiveVal(F)
  
  # Observe when the plot type changes to reset the plot_generated flag
  observeEvent(input$side_by_side_plots, {
    pca_plot_generated(F)
  })
  
  observeEvent(input$overlay_plots, {
    pca_plot_generated(F)
  })
  
  observeEvent(input$plot_pca, {
    req(input$selected_features, sheets$data)
    selected_features <- match(input$selected_features, all_PCA_columns[,2])
    selected_features <- all_PCA_columns[selected_features, 1]
    
    # Only use the sheets that are currently selected in the loaded_sheets_select input
    selected_sheets <- input$loaded_sheets_select
    combined_df <- bind_rows(sheets$data[selected_sheets], .id = "Sheet")
    
    combined_df$Sheet <- as.factor(combined_df$Sheet)
    
    selected_data <<- combined_df[, selected_features]
    
    # Get the colors for the selected sheets
    colors <- sapply(selected_sheets, function(sheet) {
      input[[paste0("color_", sheet)]]
    })
    names(colors) <- selected_sheets
    
    pca_result <<- PCA_2D(combined_df, selected_features, colors)
    
    output$plot <- renderPlot({
      pca_result$plot
    })
    
    pca_plot_generated(T)
  })
  
  output$top_5_output <- renderPrint({
    if (!pca_plot_generated()) {
      return("")
    }
    get_top_5(pca_result$pca)
  })
  
  observeEvent(input$plot_scree, {
    req(input$selected_features, sheets$data)
    selected_features <- match(input$selected_features, all_PCA_columns[,2])
    selected_features <- all_PCA_columns[selected_features, 1]
    
    # Only use the sheets that are currently selected in the loaded_sheets_select input
    selected_sheets <- input$loaded_sheets_select
    combined_df <- bind_rows(sheets$data[selected_sheets], .id = "Sheet")
    
    combined_df$Sheet <- as.factor(combined_df$Sheet)
    pca_result <- PCA_2D(combined_df, selected_features)
    
    selected_data <<- combined_df[, selected_features]
    
    output$plot <- renderPlot({
      scree_plot(pca_result$pca)
    })
  })

########################
#### TRACKING LOGIC ####
########################
  
  observeEvent(input$tracking_file, {
    req(input$tracking_file)
    
    validate(
      need(grepl("tracking\\.h5$", input$tracking_file$name), "Please upload a valid tracking HDF5 (.h5) file.")
    )
    
    tracking_load(input$tracking_file$datapath, input$tracking_tag)
    
    # Apply cleanup
    cleaned_body_parts <- sapply(body_parts, cleanup, USE.NAMES = FALSE)
    
    # Make cleaned & combined column df
    all_body_parts <<- data.frame(body_parts, cleaned_body_parts)
    
    updateSelectInput(session, "tracking_body_part", choices = all_body_parts[,2])
  })
  
  observeEvent(input$tracking_tag, {
    req(input$tracking_file)
    
    tracking_load(input$tracking_file$datapath, input$tracking_tag)
    
    # Apply cleanup
    cleaned_body_parts <- sapply(body_parts, cleanup, USE.NAMES = FALSE)
    
    # Make cleaned & combined column df
    all_body_parts <<- data.frame(body_parts, cleaned_body_parts)
    
    updateSelectInput(session, "tracking_body_part", choices = all_body_parts[,2])
  })
  
  # Reactive value to track whether the plot has been generated
  tracking_plot_generated <- reactiveVal(F)
  
  # Observe when the plot type changes to reset the plot_generated flag
  observeEvent(input$tracking_plot_type, {
    tracking_plot_generated(F)
  })
  
  # Observe when the body part changes to reset the plot_generated flag
  observeEvent(input$tracking_body_part, {
    tracking_plot_generated(F)
  })
  
  observeEvent(input$tracking_file, {
    tracking_plot_generated(F)
  })
  
  observeEvent(input$side_by_side_plots, {
    tracking_plot_generated(F)
  })
  
  observeEvent(input$overlay_plots, {
    tracking_plot_generated(F)
  })
  
  # Observe the plot data button to generate the plot
  observeEvent(input$tracking_plot_data, {
    req(input$tracking_body_part, input$tracking_plot_type)
    
    tracking_body_part <<- match(input$tracking_body_part, all_body_parts[,2])
    tracking_body_part <<- all_body_parts[tracking_body_part, 1]
    
    tracking_mouse_name <<- get_mouse_name(input$tracking_file$name)
    
    plot <- switch(input$tracking_plot_type,
                   "Speed" = plot_speed(input$tracking_tag, tracking_body_part, input$start_sec, input$end_sec),
                   "Likelihood" = plot_likelihood(input$tracking_tag, tracking_body_part, input$start_sec, input$end_sec),
                   "Trajectory" = plot_trajectory(input$tracking_tag, tracking_body_part, input$start_sec, input$end_sec),
                   "Heatmap" = plot_heatmap(input$tracking_tag, tracking_body_part, input$start_sec, input$end_sec))
    
    plot_name <- paste(input$tracking_tag, tracking_mouse_name, tracking_body_part, input$tracking_plot_type, "plot", sep = "_")
    assign(plot_name, plot, envir = plot_storage)
    
    updateSelectInput(session, "saved_plots", choices = ls(envir = plot_storage))
    
    # Render the plot
    output$plot <- renderPlot({
      plot
    })
    
    # Set plot_generated to TRUE after the plot is rendered
    tracking_plot_generated(T)
  })
  
  # Render the average output only if the plot has been generated
  output$tracking_average <- renderPrint({
    if (!tracking_plot_generated()) {
      return("")
    }
    if (input$tracking_plot_type == "Speed") {
      paste0("Average Speed: ", average_speed)
    }
    else if (input$tracking_plot_type == "Likelihood") {
      paste0("Average Likelihood (Coordinate Accuracy): ", average_likelihood)
    }
  })
  
  
  ### Batch Extract Speed Averages ###
  observeEvent(input$batch_speed_avg_files, {
    req(input$batch_speed_avg_files)
    
    batch_body_parts <<- NULL
    
    load_multiple_tracking_files(input$batch_speed_avg_files$datapath, input$batch_speed_avg_files$name)
    
    # Apply cleanup
    cleaned_batch_body_parts <- sapply(batch_body_parts, cleanup, USE.NAMES = FALSE)
    
    # Make cleaned & combined column df
    all_batch_body_parts <<- data.frame(batch_body_parts, cleaned_batch_body_parts)
    
    updateSelectInput(session, "batch_speed_avg_body_part", choices = all_batch_body_parts[,2])
  })
  
  output$downloadAvgSpeed <- downloadHandler(
    filename = function() { tolower(paste0("average_", input$batch_speed_avg_body_part,"_speeds.csv")) },
    content = function(file) {
      req(input$batch_speed_avg_files, input$batch_speed_avg_body_part)
      
      batch_body_part <- match(input$batch_speed_avg_body_part, all_batch_body_parts[,2])
      batch_body_part <- all_batch_body_parts[batch_body_part, 1]
      
      avgs_df <- compute_speed_avgs(
        file_names = input$batch_speed_avg_files$name,
        body_part = batch_body_part,
        start_sec = input$batch_speed_avg_start_sec,
        end_sec = input$batch_speed_avg_end_sec
      )
      
      if (!is.null(avgs_df)) {
        write.csv(avgs_df, file, row.names = T)
      }
    }
  )
  
########################
#### FEATURES LOGIC ####
########################
  
  observeEvent(input$features_file, {
    req(input$features_file)
    
    validate(
      need(grepl("features\\.h5$", input$features_file$name), "Please upload a valid features HDF5 (.h5) file.")
    )
    
    bboxdataPath <- input$features_file$datapath
    suppressWarnings(bboxdata <<- H5Fopen(bboxdataPath))
    
    groupls <- h5ls(bboxdata)
    bboxgroupname <- unique(groupls$group)[2]
    
    behavioral_metrics <<- names(h5read(bboxdataPath, bboxgroupname))
    filtered_metrics <- behavioral_metrics[!grepl("fps|frame_count|recording_time", behavioral_metrics)]
    
    # Apply cleanup
    cleaned_metrics <- sapply(filtered_metrics, cleanup, USE.NAMES = FALSE)
    
    # Make cleaned & combined column df
    all_metrics <<- data.frame(filtered_metrics, cleaned_metrics)
    
    updateSelectInput(session, "features_metric", choices = all_metrics[,2])
  })
  
  # Reactive value to track whether the plot has been generated
  features_plot_generated <- reactiveVal(F)
  
  # Observe when the features metric changes to reset the plot_generated flag
  observeEvent(input$features_metric, {
    features_plot_generated(F)
  })
  
  observeEvent(input$side_by_side_plots, {
    features_plot_generated(F)
  })
  
  observeEvent(input$overlay_plots, {
    features_plot_generated(F)
  })
  
  observeEvent(input$features_file, {
    features_plot_generated(F)
  })
  
  observeEvent(input$features_load_data, {
    req(input$features_file, input$features_metric)
    
    metric <<- match(input$features_metric, all_metrics[,2])
    metric <<- all_metrics[metric, 1]
    
    bboxdataPath <- input$features_file$datapath
    bboxdsetnm <- input$features_dataset_name
    
    features_mouse_name <<- get_mouse_name(input$features_file$name)
    
    plot <- features_load_and_plot(bboxdataPath, bboxdsetnm, metric, input$start_sec, input$end_sec)
    
    plot_name <- paste(metric, features_mouse_name, bboxdsetnm, "features_plot", sep = "_")
    
    assign(plot_name, plot, envir = plot_storage)
    
    updateSelectInput(session, "saved_plots", choices = ls(envir = plot_storage))
    
    output$plot <- renderPlot({
      get(plot_name, envir = plot_storage)
    })
    
    # Set plot_generated to TRUE after the plot is rendered
    features_plot_generated(T)
  })

  # Render the average output only if the plot has been generated
  output$features_average <- renderPrint({
    if (!features_plot_generated()) {
      return("")
    }
    if (is.factor(get(paste(metric, features_mouse_name, input$features_dataset_name, "features_df", sep = "_"), envir = .GlobalEnv)$values)) {
      paste0("Percentage TRUE: ", features_true_tot_pct)
    }
    else if (is.numeric(get(paste(metric, features_mouse_name, input$features_dataset_name, "features_df", sep = "_"), envir = .GlobalEnv)$values)) {
      paste0("Average: ", features_average)
    }
  })
  
  ### Batch Extract Features Averages ###
  observeEvent(input$batch_features_avg_files, {
    req(input$batch_features_avg_files)
    
    bboxdataPath <- input$batch_features_avg_files$datapath[1]
    
    suppressWarnings(bboxdata <<- H5Fopen(bboxdataPath))
    
    groupls <- h5ls(bboxdata)
    bboxgroupname <- unique(groupls$group)[2]
    
    behavioral_metrics <<- NULL
    
    behavioral_metrics <- names(h5read(bboxdataPath, bboxgroupname))
    filtered_metrics <- behavioral_metrics[!grepl("fps|frame_count", behavioral_metrics)]
    
    # Apply cleanup
    cleaned_metrics <- sapply(filtered_metrics, cleanup, USE.NAMES = FALSE)
    
    # Make cleaned & combined column df
    all_batch_metrics <<- data.frame(filtered_metrics, cleaned_metrics)
    
    updateSelectInput(session, "batch_features_avg_metric", choices = all_batch_metrics[,2])
  })
  
  output$downloadAvgFeature <- downloadHandler(
    filename = function() {tolower(paste0("batch_extraction_", input$batch_features_avg_metric,".csv"))},
    content = function(file) {
      req(input$batch_features_avg_files, input$batch_features_avg_metric)
      
      batch_metric <- match(input$batch_features_avg_metric, all_batch_metrics[,2])
      batch_metric <- all_batch_metrics[batch_metric, 1]
      
      features_avgs_df <- compute_features_avgs(
        file_paths = input$batch_features_avg_files$datapath,
        file_names = input$batch_features_avg_files$name,
        features_metric = batch_metric,
        start_sec = input$batch_features_avg_start_sec,
        end_sec = input$batch_features_avg_end_sec
      )
      
      if (!is.null(features_avgs_df)) {
        write.csv(features_avgs_df, file, row.names = T)
      }
    }
  )
  
########################
#### SUMMARY LOGIC #####
########################
  
  observeEvent(input$summary_file, {
    req(input$summary_file)
    
    validate(
      need(grepl("\\.csv$", input$summary_file$name), "Please upload a valid CSV file.")
    )
    
    suppressWarnings({
      summary_data <- read.csv(input$summary_file$datapath, row.names = 1, header = TRUE, fileEncoding = "UTF-8")
    })
    
    metrics <<- colnames(summary_data)
    mouse_names <<- rownames(summary_data)
    
    # Apply cleanup
    filtered_sum_metrics <- metrics[!grepl("PV..bin.start.end..min.|bin.duration..min.", metrics)]
    cleaned_sum_metrics <- sapply(filtered_sum_metrics, cleanup, USE.NAMES = FALSE)
    
    # Make cleaned & combined column df
    all_sum_metrics <<- data.frame(filtered_sum_metrics, cleaned_sum_metrics)
    
    updateSelectInput(session, "summary_metric", choices = all_sum_metrics[,2])
    updateSelectInput(session, "summary_mouse_name", choices = mouse_names)
  })
  
  observeEvent(input$summary_plot_data, {
    req(input$summary_file, input$summary_metric, input$summary_mouse_name)
    
    sum_metric <<- match(input$summary_metric, all_sum_metrics[,2])
    sum_metric <<- all_sum_metrics[sum_metric, 1]
    
    plot <- summary_load_and_plot(input$summary_tag, input$summary_file$datapath, input$summary_mouse_name, sum_metric)
    
    plot_name <- paste0(input$summary_tag, "_", input$summary_mouse_name, "_", sum_metric, "_summary_plot")
    assign(plot_name, plot, envir = plot_storage)
    
    updateSelectInput(session, "saved_plots", choices = ls(envir = plot_storage))
    
    output$plot <- renderPlot({
      plot
    })
  })
  
######################
#### SIDE BY SIDE ####
######################
  
  # Reactive value to track whether the the user is plotting an side_by_side plot
  side_by_side_plot_generated <- reactiveVal(F)
  
  observeEvent(input$side_by_side_plots, {
    req(input$saved_plots)
    selected_plots <- input$saved_plots
    
    # Initialize lists to group plots by type
    summary_plots <- list()
    tracking_plots <- list()
    features_plots <- list()
    other_plots <- list()
    
    for (plot_name in selected_plots) {
      plot <- get(plot_name, envir = plot_storage)
      
      # Ensure the plot is a ggplot object
      if (inherits(plot, "gg")) {
        # Determine the type of plot based on the name
        if (grepl("summary", plot_name, ignore.case = TRUE)) {
          summary_plots <- c(summary_plots, list(plot))
        } else if (grepl("speed", plot_name, ignore.case = TRUE)) {
          tracking_plots <- c(tracking_plots, list(plot))
        } else if (grepl("features", plot_name, ignore.case = TRUE)) {
          features_plots <- c(features_plots, list(plot))
        } else {
          other_plots <- c(other_plots, list(plot))
        }
      } else {
        warning(paste("The selected item", plot_name, "is not a ggplot object and will be skipped."))
      }
    }
    
    # Function to calculate max y limit for a list of plots
    calculate_max_y_limit <- function(plots) {
      max(sapply(plots, function(plot) {
        tryCatch({
          y_range <- ggplot_build(plot)$layout$panel_params[[1]]$y.range
          if (is.numeric(y_range)) {
            max(y_range, na.rm = TRUE)
          } else {
            NA
          }
        }, error = function(e) NA)
      }, simplify = TRUE, USE.NAMES = FALSE), na.rm = TRUE)
    }
    
    calculate_min_y_limit <- function(plots) {
      min(sapply(plots, function(plot) {
        tryCatch({
          y_range <- ggplot_build(plot)$layout$panel_params[[1]]$y.range
          if (is.numeric(y_range)) {
            min(y_range, na.rm = TRUE)
          } else {
            NA
          }
        }, error = function(e) NA)
      }, simplify = TRUE, USE.NAMES = FALSE), na.rm = TRUE)
    }
    
    calculate_max_x_limit <- function(plots) {
      max(sapply(plots, function(plot) {
        tryCatch({
          x_range <- ggplot_build(plot)$layout$panel_params[[1]]$x.range
          if (is.numeric(x_range)) {
            max(x_range, na.rm = TRUE)
          } else {
            NA
          }
        }, error = function(e) NA)
      }, simplify = TRUE, USE.NAMES = FALSE), na.rm = TRUE)
    }
    
    calculate_min_x_limit <- function(plots) {
      min(sapply(plots, function(plot) {
        tryCatch({
          x_range <- ggplot_build(plot)$layout$panel_params[[1]]$x.range
          if (is.numeric(x_range)) {
            min(x_range, na.rm = TRUE)
          } else {
            NA
          }
        }, error = function(e) NA)
      }, simplify = TRUE, USE.NAMES = FALSE), na.rm = TRUE)
    }
    
    # Initialize combined plots
    combined_summary_plot <- NULL
    combined_tracking_plot <- NULL
    combined_features_plot <- NULL
    
    # Apply scaling to summary plots
    if (length(summary_plots) > 0) {
      max_y_limit_summary <- calculate_max_y_limit(summary_plots)
      min_y_limit_summary <- calculate_min_y_limit(summary_plots)
      min_x_limit_summary <- calculate_min_x_limit(summary_plots)
      max_x_limit_summary <- calculate_max_x_limit(summary_plots)
      
      # Apply the max y limit to all summary plots
      scaled_summary_plots <- lapply(summary_plots, function(plot) {
        plot + coord_cartesian(xlim = c(min_x_limit_summary, max_x_limit_summary),
                               ylim = c(min_y_limit_summary, max_y_limit_summary))
      })
      
      # Combine the summary plots
      combined_summary_plot <- patchwork::wrap_plots(scaled_summary_plots)
    }
    
    # Apply scaling to tracking plots
    if (length(tracking_plots) > 0) {
      max_y_limit_tracking <- calculate_max_y_limit(tracking_plots)
      min_y_limit_tracking <- calculate_min_y_limit(tracking_plots)
      max_x_limit_tracking <- calculate_max_x_limit(tracking_plots)
      min_x_limit_tracking <- calculate_min_x_limit(tracking_plots)
      
      # Apply the max y limit to all tracking plots
      scaled_tracking_plots <- lapply(tracking_plots, function(plot) {
        plot + coord_cartesian(xlim = c(min_x_limit_tracking, max_x_limit_tracking),
                               ylim = c(min_y_limit_tracking, max_y_limit_tracking))
      })
      
      # Combine the tracking plots
      combined_tracking_plot <- patchwork::wrap_plots(scaled_tracking_plots)
    }
    
    # Apply scaling to features plots
    if (length(features_plots) > 0) {
      max_y_limit_features <- calculate_max_y_limit(features_plots)
      min_y_limit_features <- calculate_min_y_limit(features_plots)
      max_x_limit_features <- calculate_max_x_limit(features_plots)
      min_x_limit_features <- calculate_min_x_limit(features_plots)
      
      # Apply the max y limit to all features plots
      scaled_features_plots <- lapply(features_plots, function(plot) {
        plot + coord_cartesian(xlim = c(min_x_limit_features, max_x_limit_features),
                               ylim = c(min_y_limit_features, max_y_limit_features))
      })
      
      # Combine the features plots
      combined_features_plot <- patchwork::wrap_plots(scaled_features_plots)
    }
    
    # Combine all the combined plots, ensuring non-null plots are included
    combined_plot <<- patchwork::wrap_plots(
      c(
        if (!is.null(combined_summary_plot)) list(combined_summary_plot) else NULL,
        if (!is.null(combined_tracking_plot)) list(combined_tracking_plot) else NULL,
        if (!is.null(combined_features_plot)) list(combined_features_plot) else NULL,
        if (length(other_plots) > 0) other_plots else NULL
      ), ncol = 1
    )
    
    side_by_side_plot_generated(T)

    # Render the combined plot
    output$plot <- renderPlot({
      if (!is.null(combined_plot)) {
        combined_plot
      } else {
        plot(NULL)
      }
    })
  })
  
  # If any other plot is generated, make side_by_side_plot_generated FALSE
  observeEvent(input$features_load_data, {
    side_by_side_plot_generated(F)
  })
  observeEvent(input$tracking_plot_data, {
    side_by_side_plot_generated(F)
  })
  observeEvent(input$summary_plot_data, {
    side_by_side_plot_generated(F)
  })
  observeEvent(input$overlay_plots, {
    side_by_side_plot_generated(F)
  })
  
#################
#### OVERLAY ####
#################
  
  overlay_plot_generated <- reactiveVal(F)
  
  observeEvent(input$overlay_plots, {
    req(input$saved_plots)
    selected_plots <- input$saved_plots
    
    # Initialize lists to group plots by type
    speed_dfs <- list()
    likelihood_dfs <- list()
    features_num_dfs <- list()
    features_factor_dfs <- list()
    other_dfs <- list()
    
    for (plot_name in selected_plots) {
      
      df_name <- str_replace(plot_name, "plot", "df")
      
      if (grepl("summary", df_name, ignore.case = TRUE)) {
        next
      }
      
      df <- get(df_name, envir = .GlobalEnv)
      
      # Split up the plots into their respective lists
      if (grepl("speed", df_name, ignore.case = TRUE)) {
        speed_dfs[[df_name]] <- df
      } else if (grepl("likelihood", df_name, ignore.case = TRUE)) {
        likelihood_dfs[[df_name]] <- df
      } else if (grepl("features", df_name, ignore.case = TRUE) && is.numeric(df[,1])) {
        features_num_dfs[[df_name]] <- df
      } else if (grepl("features", df_name, ignore.case = TRUE) && is.factor(df[,1])) {
        features_factor_dfs[[df_name]] <- df
      } else {
        other_dfs[[df_name]] <- df
      }
    }
    
    all_df_names <- na.omit(c(names(speed_dfs), names(likelihood_dfs), names(features_num_dfs), names(features_factor_dfs)))
    
    if (length(all_df_names) == 0) {
      showNotification("Please select SPEED, LIKELIHOOD, or FEATURES plots to overlay.", type = "warning")
      return()
    }
    
    color_palette <- as.list(rainbow(length(all_df_names)))
        
    names(color_palette) <- all_df_names
    
    output$dynamic_plot_ui <- renderUI({
      fluidRow(
        column(12,
               uiOutput("checkbox_ui")  # Checkboxes take full width below the plot
        )
      )
    })
    
  output$checkbox_ui <- renderUI({
    if (overlay_plot_generated()){
      checkboxGroupButtons(
      inputId = "selected_dfs",
      label = NULL,
      choices = all_df_names,
      selected = all_df_names,
      direction = "vertical",
      individual = TRUE,
      size = "sm"
    )}
  })
    
  # Render the combined plot
  output$plot <- renderPlot({
    color_list <- unlist(color_palette)
    overlay_plot <<- ggplot() +
      labs(x = "Time", y = "Values") +
      ggtitle("Combined Overlay Plot") +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 18)
      ) +
      scale_colour_manual(name = "", values = color_list)
    
    # Speed Dataframes
    if (length(speed_dfs) > 0) {
      for (i in seq_along(speed_dfs)) {
        speed_dfs[[i]]$df_name <- names(speed_dfs)[i]
        if (speed_dfs[[i]]$df_name[1] %in% input$selected_dfs) {
          overlay_plot <<- overlay_plot + 
            geom_line(data = speed_dfs[[i]], aes(x = time_vector, y = speeds, color = df_name))
        }
      }
    }
    
    # Likelihood Dataframes
    if (length(likelihood_dfs) > 0) {
      for (i in seq_along(likelihood_dfs)) {
        likelihood_dfs[[i]]$df_name <- names(likelihood_dfs)[i]
        if (likelihood_dfs[[i]]$df_name[1] %in% input$selected_dfs) {
          overlay_plot <<- overlay_plot + 
            geom_line(data = likelihood_dfs[[i]], aes(x = time_vector, y = likelihood_row, color = df_name))
        }
      }
    }
    
    # Features Factor Dataframes
    if (length(features_factor_dfs) > 0) {
      combined_data <- data.frame()
      
      for (i in seq_along(features_factor_dfs)) {
        features_factor_dfs[[i]]$df_name <- names(features_factor_dfs)[i]
        if (features_factor_dfs[[i]]$df_name[1] %in% input$selected_dfs) {
          features_factor_dfs[[i]]$values <- as.numeric(features_factor_dfs[[i]]$values)
          combined_data <- rbind(combined_data, features_factor_dfs[[i]])
        }
      }
      
      if (nrow(combined_data) > 0) {
        
        combined_data$values <- combined_data$values - 1
        # Define color mapping, assigning colors to 1 and a fully transparent color to 0
        unique_dfs <- unique(combined_data$df_name)
        
        # Ensure color_palette has names matching unique_dfs
        color_mapping <- setNames(color_palette[unique_dfs], unique_dfs)
        
        # Assign colors based on values
        combined_data$fill_key <- ifelse(combined_data$values == 1, 
                                         combined_data$df_name,
                                         NA
                                         )
        
        # Add the tiles to the base plot with proper transparency handling
        overlay_plot <<- overlay_plot +
          geom_tile(data = combined_data, 
                    aes(x = times, y = values, width = time_increments, height = Inf, fill = fill_key)) +  # Set alpha for the geom_tile
          scale_fill_manual(values = c(color_mapping, "NA" = "transparent"), 
                            name = "",
                            breaks = unique_dfs,
                            labels = names(color_palette[unique_dfs]),
                            na.value = "transparent")
      }
    }
    
    # Features Number Dataframes
    if (length(features_num_dfs) > 0) {
      for (i in seq_along(features_num_dfs)) {
        features_num_dfs[[i]]$df_name <- names(features_num_dfs)[i]
        if (features_num_dfs[[i]]$df_name[1] %in% input$selected_dfs) {
          overlay_plot <<- overlay_plot + 
            geom_line(data = features_num_dfs[[i]], aes(x = times, y = values, color = df_name))
        }
      }
    }
    
    overlay_plot_generated(T)
    
    print(overlay_plot)
  })
})
  
  # If any other plot is generated, make side_by_side_plot_generated FALSE
  observeEvent(input$features_load_data, {
    overlay_plot_generated(F)
  })
  observeEvent(input$tracking_plot_data, {
    overlay_plot_generated(F)
  })
  observeEvent(input$summary_plot_data, {
    overlay_plot_generated(F)
  })
  observeEvent(input$pca_plot_generated, {
    overlay_plot_generated(F)
  })
  observeEvent(input$side_by_side_plots, {
    overlay_plot_generated(F)
  })
  
##################
#### DOWNLOAD ####
##################
  #########################
  ### Plot PDF Download ###
  #########################
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      if (side_by_side_plot_generated()) {
        "side_by_side_plot.pdf"
      } else if (overlay_plot_generated()) {
        "overlay_plot.pdf"
      } else if (input$analysis_type == "PCA") {
        "PCA_plot.pdf"
      } else if (input$analysis_type == "Tracking Data") {
        paste(input$tracking_tag, tracking_mouse_name, input$tracking_body_part, input$tracking_plot_type, "plot.pdf", sep = "_")
      } else if (input$analysis_type == "Features Data") {
        paste(input$features_metric, features_mouse_name, input$features_dataset_name, "plot.pdf", sep = "_")
      } else if (input$analysis_type == "Summary Data") {
        paste(input$summary_tag, input$summary_mouse_name, input$summary_metric, "plot.pdf", sep = "_")
      }
    },
    content = function(file) {
      plot <-
      if (side_by_side_plot_generated()) {
        combined_plot
      } else if (overlay_plot_generated()) {
        overlay_plot
      } else if (input$analysis_type == "PCA") {
        pca_result$plot
      } else if (input$analysis_type == "Tracking Data") {
      plot_name <- paste(input$tracking_tag, tracking_mouse_name, tracking_body_part, input$tracking_plot_type, "plot", sep = "_")
      get(plot_name, envir = plot_storage)
      } else if (input$analysis_type == "Features Data") {
        plot_name <- paste(metric, features_mouse_name, input$features_dataset_name, "features_plot", sep = "_")
        get(plot_name, envir = plot_storage)
      } else if (input$analysis_type == "Summary Data") {
        plot_name <- paste0(input$summary_tag, "_", input$summary_mouse_name, "_", sum_metric, "_summary_plot")
        get(plot_name, envir = plot_storage)
      }
      ggsave(file, plot = plot, device = "pdf", width = 12, height = 8, units = "in")
    }
  )
  
  #####################
  ### Data Download ###
  #####################
  output$downloadData <- downloadHandler(
    filename = function() {
      if (input$analysis_type == "Tracking Data") {
        paste(input$tracking_tag, tracking_mouse_name, tracking_body_part, input$tracking_plot_type, "data.csv", sep = "_")
      } else if (input$analysis_type == "Features Data") {
        paste(input$features_dataset_name, features_mouse_name, input$features_metric, "data.csv", sep = "_")
      } else if (input$analysis_type == "Summary Data") {
        paste(input$summary_tag, input$summary_mouse_name, input$summary_metric, "data.csv", sep = "_")
      } else if (input$analysis_type == "PCA") {
        "PCA_data.csv"
      } else {
        "data.csv"
      }
    },
    content = function(file) {
      data <- NULL
      
      if (input$analysis_type == "Tracking Data") {
        if (input$tracking_plot_type == "Speed") {
          data <- get(paste(input$tracking_tag, tracking_mouse_name, tracking_body_part, "Speed_df", sep = "_"), envir = .GlobalEnv)
        } else if (input$tracking_plot_type == "Likelihood") {
          data <- get(paste(input$tracking_tag, tracking_mouse_name, tracking_body_part, "Likelihood_df", sep = "_"), envir = .GlobalEnv)
        } else if (input$tracking_plot_type == "Trajectory" || input$tracking_plot_type == "Heatmap") {
          data <- get(paste(input$tracking_tag, tracking_mouse_name, tracking_body_part, "Trajectory_df", sep = "_"), envir = .GlobalEnv)
        }
      } else if (input$analysis_type == "Features Data") {
        data <- get(paste(metric, features_mouse_name, input$features_dataset_name, "features_df", sep = "_"), envir = .GlobalEnv)
      } else if (input$analysis_type == "PCA") {
        data <- selected_data
      }
      
      if (!is.null(data)) {
        write.csv(data, file)
      }
    }
  )
}

shinyApp(ui = ui, server = server)
