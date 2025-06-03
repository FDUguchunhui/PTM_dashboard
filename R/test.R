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
perform_statistical_analysis <- function(annotated_data, group1_types, group2_types, assay, normalization, test_type = "t_test", pvalue_threshold = 0.05) {
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
      normalization = normalization
    )$get_annotated_data()

    group2_data <- annotated_data$create_filtered_view(
      cancer_type = group2_types,
      assay = assay,
      normalization = normalization
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

  # Initialize results dataframe
  n_rows <- nrow(group1_data)
  results <- data.frame(
    feature_id = 1:n_rows,
    log2_fold_change = numeric(n_rows),
    p_value = numeric(n_rows),
    mean_group1 = numeric(n_rows),
    mean_group2 = numeric(n_rows),
    n_group1 = integer(n_rows),
    n_group2 = integer(n_rows),
    num_non_NA_1 = integer(n_rows),
    num_non_NA_2 = integer(n_rows),
    auc = numeric(n_rows),
    stringsAsFactors = FALSE
  )

  # Add row metadata to results
  results <- cbind(row_metadata, results[, c("log2_fold_change", "p_value", "mean_group1", "mean_group2", "n_group1", "n_group2", "num_non_NA_1", "num_non_NA_2", "auc")])

  # Function to calculate AUC using pROC package
  calculate_auc <- function(group1_vals, group2_vals) {
    if (length(group1_vals) == 0 || length(group2_vals) == 0) {
      return(NA)
    }

    # Check if all values are the same (would cause pROC to fail)
    all_values <- c(group1_vals, group2_vals)
    if (length(unique(all_values)) <= 1) {
      return(0.5)  # AUC = 0.5 when all values are identical
    }

    # Combine values and create labels
    labels <- c(rep(1, length(group1_vals)), rep(0, length(group2_vals)))

    # Use pROC to calculate AUC
    tryCatch({
      roc_obj <- pROC::roc(response = labels, predictor = all_values, quiet = TRUE)
      return(as.numeric(pROC::auc(roc_obj)))
    }, error = function(e) {
      return(NA)
    })
  }

  # Perform statistical test for each row (peptide/modified peptide)
  for (i in 1:n_rows) {
    group1_values <- as.numeric(group1_data[i, group1_intensity_cols])
    group2_values <- as.numeric(group2_data[i, group2_intensity_cols])

    # Count non-missing values before filtering
    results$num_non_NA_1[i] <- sum(!is.na(group1_values))
    results$num_non_NA_2[i] <- sum(!is.na(group2_values))

    # Remove NAs and zeros for more robust analysis
    group1_values <- group1_values[!is.na(group1_values) & group1_values > 0]
    group2_values <- group2_values[!is.na(group2_values) & group2_values > 0]

    # Record number of non-missing observations
    results$n_group1[i] <- length(group1_values)
    results$n_group2[i] <- length(group2_values)

    # Handle edge cases where one or both groups are empty
    if (length(group1_values) == 0 || length(group2_values) == 0) {
      results$mean_group1[i] <- ifelse(length(group1_values) == 0, NA, mean(group1_values))
      results$mean_group2[i] <- ifelse(length(group2_values) == 0, NA, mean(group2_values))
      results$log2_fold_change[i] <- NA
      results$p_value[i] <- NA
      results$auc[i] <- NA
      next
    }

    # Calculate means
    mean_group1 <- mean(group1_values)
    mean_group2 <- mean(group2_values)

    results$mean_group1[i] <- mean_group1
    results$mean_group2[i] <- mean_group2

    # Calculate log2 fold change
    if (mean_group1 > 0 & mean_group2 > 0) {
      results$log2_fold_change[i] <- log2(mean_group1 / mean_group2)
    } else {
      results$log2_fold_change[i] <- NA
    }

    # Calculate AUC using pROC
    results$auc[i] <- calculate_auc(group1_values, group2_values)

    # Perform statistical test if we have enough data points
    min_samples_required <- ifelse(test_type == "wilcoxon", 1, 2)

    if (length(group1_values) >= min_samples_required & length(group2_values) >= min_samples_required) {
      tryCatch({
        if (test_type == "t_test") {
          test_result <- t.test(group1_values, group2_values, var.equal = FALSE)
          results$p_value[i] <- test_result$p.value
        } else if (test_type == "wilcoxon") {
          # Wilcoxon rank-sum test (Mann-Whitney U test)
          test_result <- wilcox.test(group1_values, group2_values, exact = FALSE)
          results$p_value[i] <- test_result$p.value
        }
      }, error = function(e) {
        results$p_value[i] <- NA
      })
    } else {
      results$p_value[i] <- NA
    }
  }

  # Add significance indicator
  results$significant <- results$p_value < pvalue_threshold & !is.na(results$p_value)

  # Sort by p-value
  results <- results[order(results$p_value, na.last = TRUE), ]

  return(results)
}
