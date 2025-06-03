# Function to perform statistical analysis (t-test or Wilcoxon rank-sum test)
# This function performs statistical comparison between two groups of cancer types
# Returns a data frame with log2 fold change, p-values, significance indicators, and AUC
# Parameters:
#   - annotated_data: AnnotatedData object containing the experimental data
#   - group1_types: vector of cancer types for group 1
#   - group2_types: vector of cancer types for group 2
#   - assay: assay type to analyze ('FT' or 'IgB')
#   - normalization: normalization method ('none', 'median', 'plate_median')
#   - test_type: statistical test to perform ('t_test' or 'wilcoxon')
#   - pvalue_threshold: threshold for significance testing (default 0.05)
#   - missing_fill: value to fill missing values with (default NA)
perform_statistical_analysis <- function(annotated_data, group1_types, group2_types, assay, normalization, test_type = "t_test", pvalue_threshold = 0.05, missing_fill = NA) {
  # Load required library
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("pROC package is required but not installed. Please install it with: install.packages('pROC')")
  }

  # Validate inputs
  if (is.null(annotated_data)) {
    stop("annotated_data cannot be NULL")
  }

  if (length(group1_types) == 0 || length(group2_types) == 0) {
    stop("Both group1_types and group2_types must contain at least one cancer type")
  }

  if (!test_type %in% c("t_test", "wilcoxon")) {
    stop("test_type must be either 't_test' or 'wilcoxon'")
  }

  if (pvalue_threshold <= 0 || pvalue_threshold >= 1) {
    stop("pvalue_threshold must be between 0 and 1")
  }

  # Check if annotated_data has required methods
  if (!exists("get_annotated_data", where = annotated_data) || !exists("get_row_metadata", where = annotated_data)) {
    stop("annotated_data must have get_data() and get_row_metadata() methods")
  }

  # Get data for each group
  tryCatch({
    group1_data <- annotated_data$create_filtered_view(
      cancer_type = group1_types,
      assay = assay,
      normalization = normalization,
      missing_fill = missing_fill
    )$get_annotated_data()

    group2_data <- annotated_data$create_filtered_view(
      cancer_type = group2_types,
      assay = assay,
      normalization = normalization,
      missing_fill = missing_fill
    )$get_annotated_data()
  }, error = function(e) {
    stop(paste("Error retrieving data:", e$message))
  })

  # Validate that data was retrieved
  if (is.null(group1_data) || is.null(group2_data)) {
    stop("Failed to retrieve data for one or both groups")
  }

  if (nrow(group1_data) == 0 || nrow(group2_data) == 0) {
    stop("No data found for one or both groups")
  }

  if (nrow(group1_data) != nrow(group2_data)) {
    stop("Group 1 and Group 2 data have different numbers of features")
  }

  # Get row metadata columns to identify the intensity columns
  row_metadata <- annotated_data$get_row_metadata()
  row_meta_cols <- colnames(row_metadata)

  # Get intensity column names (exclude row metadata columns)
  group1_intensity_cols <- setdiff(colnames(group1_data), row_meta_cols)
  group2_intensity_cols <- setdiff(colnames(group2_data), row_meta_cols)

  # Check if there are any intensity columns
  if (length(group1_intensity_cols) == 0 || length(group2_intensity_cols) == 0) {
    stop("No intensity columns found after removing metadata columns")
  }

  # Extract intensity matrices for vectorized operations
  group1_matrix <- as.matrix(group1_data[, group1_intensity_cols, drop = FALSE])
  group2_matrix <- as.matrix(group2_data[, group2_intensity_cols, drop = FALSE])

  # Convert to numeric and replace 0 with NA for more robust analysis
  group1_matrix <- apply(group1_matrix, 2, as.numeric)
  group2_matrix <- apply(group2_matrix, 2, as.numeric)
  group1_matrix[group1_matrix <= 0] <- NA
  group2_matrix[group2_matrix <= 0] <- NA

  n_rows <- nrow(group1_matrix)

  # Pre-allocate result vectors for better performance
  log2_fc <- numeric(n_rows)
  p_value <- numeric(n_rows)
  mean_group1 <- numeric(n_rows)
  mean_group2 <- numeric(n_rows)
  n_group1 <- integer(n_rows)
  n_group2 <- integer(n_rows)
  num_non_NA_1 <- integer(n_rows)
  num_non_NA_2 <- integer(n_rows)
  auc <- numeric(n_rows)

  # Vectorized calculations
  # Count non-NA values (before filtering zeros)
  group1_orig <- as.matrix(group1_data[, group1_intensity_cols, drop = FALSE])
  group2_orig <- as.matrix(group2_data[, group2_intensity_cols, drop = FALSE])
  group1_orig <- apply(group1_orig, 2, as.numeric)
  group2_orig <- apply(group2_orig, 2, as.numeric)

  num_non_NA_1 <- rowSums(!is.na(group1_orig) & group1_orig != 0)
  num_non_NA_2 <- rowSums(!is.na(group2_orig) & group2_orig != 0)

  # Count valid values (after removing zeros and NAs)
  n_group1 <- ncol(group1_matrix)
  n_group2 <- ncol(group2_matrix)

  # Calculate means using vectorized operations
  mean_group1 <- rowMeans(group1_matrix, na.rm = TRUE)
  mean_group2 <- rowMeans(group2_matrix, na.rm = TRUE)

  # Calculate log2 fold change vectorized
  valid_means <- (mean_group1 > 0) & (mean_group2 > 0) & !is.na(mean_group1) & !is.na(mean_group2)
  log2_fc[valid_means] <- log2(mean_group1[valid_means] / mean_group2[valid_means])
  log2_fc[!valid_means] <- NA

  # Optimized AUC calculation function
  calculate_auc_fast <- function(g1_vals, g2_vals) {
    # Remove NAs
    g1_vals <- g1_vals[!is.na(g1_vals)]
    g2_vals <- g2_vals[!is.na(g2_vals)]

    if (length(g1_vals) == 0 || length(g2_vals) == 0) {
      return(NA)
    }

    # Check if all values are the same
    all_vals <- c(g1_vals, g2_vals)
    if (length(unique(all_vals)) <= 1) {
      return(0.5)
    }

    # Fast Mann-Whitney U calculation
    n1 <- length(g1_vals)
    n2 <- length(g2_vals)

    # Use vectorized comparison
    comparisons <- sum(outer(g1_vals, g2_vals, ">"))
    ties <- sum(outer(g1_vals, g2_vals, "=="))

    u_stat <- comparisons + 0.5 * ties
    return(u_stat / (n1 * n2))
  }

  # Calculate AUC for each row using apply (faster than for loop)
  auc <- apply(cbind(group1_matrix, group2_matrix), 1, function(row) {
    g1_vals <- row[1:ncol(group1_matrix)]
    g2_vals <- row[(ncol(group1_matrix) + 1):length(row)]
    calculate_auc_fast(g1_vals, g2_vals)
  })

  # Optimized statistical testing
  if (test_type == "t_test") {
    # Vectorized t-test calculation
    p_value <- apply(cbind(group1_matrix, group2_matrix), 1, function(row) {
      g1_vals <- row[1:ncol(group1_matrix)]
      g2_vals <- row[(ncol(group1_matrix) + 1):length(row)]

      # Remove NAs
      g1_vals <- g1_vals[!is.na(g1_vals)]
      g2_vals <- g2_vals[!is.na(g2_vals)]

      if (length(g1_vals) < 2 || length(g2_vals) < 2) {
        return(NA)
      }

      tryCatch({
        t.test(g1_vals, g2_vals, var.equal = FALSE)$p.value
      }, error = function(e) NA)
    })
  } else if (test_type == "wilcoxon") {
    # Vectorized Wilcoxon test calculation
    p_value <- apply(cbind(group1_matrix, group2_matrix), 1, function(row) {
      g1_vals <- row[1:ncol(group1_matrix)]
      g2_vals <- row[(ncol(group1_matrix) + 1):length(row)]

      # Remove NAs
      g1_vals <- g1_vals[!is.na(g1_vals)]
      g2_vals <- g2_vals[!is.na(g2_vals)]

      if (length(g1_vals) < 1 || length(g2_vals) < 1) {
        return(NA)
      }

      tryCatch({
        wilcox.test(g1_vals, g2_vals, exact = FALSE)$p.value
      }, error = function(e) NA)
    })
  }

  # Create results dataframe efficiently
  results <- data.frame(
    log2_fc = log2_fc,
    p_value = p_value,
    mean_group1 = mean_group1,
    mean_group2 = mean_group2,
    n_group1 = n_group1,
    n_group2 = n_group2,
    num_non_NA_1 = num_non_NA_1,
    num_non_NA_2 = num_non_NA_2,
    auc = auc,
    stringsAsFactors = FALSE
  )

  # Add row metadata to results
  results <- cbind(row_metadata, results)

  # Add significance indicator
  results$significant <- results$p_value < pvalue_threshold & !is.na(results$p_value)

  # Sort by p-value
  results <- results[order(results$p_value, na.last = TRUE), ]

  return(results)
}
