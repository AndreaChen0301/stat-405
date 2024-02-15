##Yahan Chen; chen2254@wisc.edu

require("FITSio")

rm(list=ls())

# Standardize the target spectrum
target <- readFrameFromFITS("cB58_Lyman_break.fit")$FLUX
standardize <- function(v) {
  (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
}
target_standardized <- standardize(target)

# Function to find the best slice
find_best_slice <- function(data, target_standardized, target_length) {
  min_distance <- Inf
  best_start <- 1
  distances <- length(data$flux) - target_length + 1
  if (length(data$flux) < target_length) {
    min_distance = sqrt(sum((standardize(data$flux) - target_standardized)^2, na.rm = TRUE))
    best_start = 0
    list(best_distance = min_distance, best_start = best_start)
  }
  else{
    for (start_index in 1:(length(data$flux) - target_length + 1)) {
      slice <- data$flux[start_index:(start_index + target_length - 1)]
      slice_standardized <- standardize(slice)
      distance <- sqrt(sum((slice_standardized - target_standardized)^2, na.rm = TRUE))
      distances[start_index] <- distance
      
      if (distance < min_distance) {
        min_distance <- distance
        best_start <- start_index
      }
    }
  }
  list(best_distance = min_distance, best_start = best_start)
}

# Apply the function to each of the 100 spectrum
spectra_files <- list.files(path = "data", pattern = "\\.fits$", full.names = TRUE)

results <- lapply(spectra_files, function(file_path) {
  data <- readFrameFromFITS(file_path)
  good_data <- data[data$and_mask == 0,]
  best_slice_info <- find_best_slice(good_data, target_standardized, length(target_standardized))
  
  list(distance = best_slice_info$best_distance, spectrumID = basename(file_path), i = best_slice_info$best_start)
})

# Convert results to a data frame
results_df <- do.call(rbind, lapply(results, data.frame, stringsAsFactors = FALSE))
results_df <- results_df[order(results_df$distance),]

write.csv(results_df, "hw2.csv", row.names = FALSE)


