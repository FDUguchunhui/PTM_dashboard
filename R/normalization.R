batch_effect_norm <- function(data, metadata){

  # Step 1: Compute Total Quantity Per Well
  summarized_data <- data.frame(
    plate = metadata$plate,
    total_intensity = colSums(data[, metadata$Run, drop = FALSE], na.rm = TRUE))

  # Step 2: Calculate Median of Each Well
  plate_medians <- summarized_data %>%
    group_by(plate) %>%
    summarise(median_sum = median(total_intensity), .groups = 'drop')

  # Step 3: Compute Correction Factor for Each Well
  max_median <- max(plate_medians$median_sum)
  plate_medians$cor_factor <- max_median / plate_medians$median_sum

  # Step 4: Normalize Each plate in Wide Format
  normalized_data <- data
  for (i in seq_along(metadata$Run)) {
    sample_name <- metadata$Run[i]
    plate_name <- metadata$plate[i]
    correction_factor <- plate_medians$cor_factor[plate_medians$plate == plate_name]
    normalized_data[, sample_name] <- data[, sample_name] * correction_factor
  }

  # only keep 2 digits after the decimal point
  normalized_data <- round(normalized_data, digits = 2)

  return(normalized_data)

}

batch_effect_norm_by_assay <- function(data, metadata) {
  # Ensure metadata contains 'assay' column
  if (!"assay" %in% colnames(metadata)) {
    stop("Error: 'assay' column is missing in metadata.")
  }

  # Initialize normalized data storage
  normalized_data <- NULL

  # Process each assay separately
  assay_types <- unique(metadata$assay)

  for (assay in assay_types) {
    # Subset data and metadata for the current assay type
    assay_metadata <- metadata[metadata$assay == assay, ]
    assay_data <- data[, assay_metadata$Run, drop = FALSE]

    # Perform batch effect normalization
    normalized_assay_data <- batch_effect_norm(assay_data, assay_metadata)

    # if normalized_data is NULL, assign normalized_assay_data to normalized_data
    if (is.null(normalized_data)) {
      normalized_data <- normalized_assay_data
    }
    else {
      # Merge normalized data back
      normalized_data <- cbind(normalized_data, normalized_assay_data)
    }
  }

  return(normalized_data)
}



median_normalization <- function(mat) {
  sample_medians <- apply(mat, 2, median, na.rm = TRUE)  # Median of each column
  global_median <- median(sample_medians, na.rm = TRUE)  # Global median

  # Normalize: Scale each column by median ratio
  norm_mat <- sweep(mat, 2, sample_medians, "/") * global_median
  return(norm_mat)
}

# Example Matrix (Columns = Samples, Rows


