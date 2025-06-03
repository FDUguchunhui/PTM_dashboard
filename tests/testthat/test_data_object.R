library(testthat)
library(R6)

# Create the test data once using the shared helper function
test_objects <- create_test_data()

# Test suite for AnnotatedData class
test_that("AnnotatedData initialization works correctly", {
  # Use shared test data
  annotated_data <- AnnotatedData$new(
    test_objects$data,
    test_objects$col_metadata,
    test_objects$row_metadata
  )

  expect_true(inherits(annotated_data, "AnnotatedData"))
  expect_equal(dim(annotated_data$data), c(10, 20))
  expect_equal(nrow(annotated_data$col_metadata), 20)
  expect_equal(nrow(annotated_data$row_metadata), 10)
})

test_that("AnnotatedData validation catches input errors", {
  # Test non-matrix data
  expect_error(
    AnnotatedData$new(data.frame(test_objects$data), test_objects$col_metadata, test_objects$row_metadata),
    "Data must be a matrix"
  )

  # Test non-data.frame col_metadata
  expect_error(
    AnnotatedData$new(test_objects$data, as.matrix(test_objects$col_metadata), test_objects$row_metadata),
    "col_metadata must be a data.frame"
  )

  # Test non-data.frame row_metadata
  expect_error(
    AnnotatedData$new(test_objects$data, test_objects$col_metadata, as.matrix(test_objects$row_metadata)),
    "row_metadata must be a data.frame"
  )

  # Test mismatched dimensions
  wrong_col_metadata <- test_objects$col_metadata[1:15, ]
  expect_error(
    AnnotatedData$new(test_objects$data, wrong_col_metadata, test_objects$row_metadata)
  )

  wrong_row_metadata <- test_objects$row_metadata[1:5, ]
  expect_error(
    AnnotatedData$new(test_objects$data, test_objects$col_metadata, wrong_row_metadata)
  )

  # Test missing required columns
  incomplete_col_metadata <- test_objects$col_metadata[, c("Run", "assay")]  # Missing "Cancer Type"
  expect_error(
    AnnotatedData$new(test_objects$data, incomplete_col_metadata, test_objects$row_metadata),
    "col_metadata is missing required columns: Cancer Type"
  )
})

test_that("AnnotatedData accessor methods work correctly", {
  annotated_data <- AnnotatedData$new(
    test_objects$data,
    test_objects$col_metadata,
    test_objects$row_metadata
  )

  # Test get_col_metadata
  expect_equal(annotated_data$get_col_metadata(), test_objects$col_metadata)

  # Test get_row_metadata
  expect_equal(annotated_data$get_row_metadata(), test_objects$row_metadata)

  # Test get_original_data
  original <- annotated_data$get_original_data()
  expect_equal(original$data, test_objects$data)
  expect_equal(original$col_metadata, test_objects$col_metadata)
  expect_equal(original$row_metadata, test_objects$row_metadata)

  # Test get_info
  info <- annotated_data$get_info()
  expect_equal(info$n_features, 10)
  expect_equal(info$n_samples, 20)
  expect_equal(sort(info$cancer_types), c("Breast", "Lung"))
  expect_equal(sort(info$assays), c("Proteomics", "RNA-seq"))
  expect_equal(info$sample_ids, paste0("Run_", 1:20))
})

test_that("AnnotatedData sample checking methods work correctly", {
  annotated_data <- AnnotatedData$new(
    test_objects$data,
    test_objects$col_metadata,
    test_objects$row_metadata
  )

  # Test has_samples
  expect_true(annotated_data$has_samples("Breast", "RNA-seq"))
  expect_true(annotated_data$has_samples("Lung", "Proteomics"))
  expect_false(annotated_data$has_samples("Pancreatic", "RNA-seq"))
  expect_false(annotated_data$has_samples("Breast", "Metabolomics"))

  # Test get_sample_counts
  counts <- annotated_data$get_sample_counts()
  expect_equal(nrow(counts), 4)  # 2 cancer types x 2 assays
  expect_true(all(counts$count == 5))  # Each combination should have 5 samples
})

test_that("AnnotatedData get_annotated_data works correctly", {
  annotated_data <- AnnotatedData$new(
    test_objects$data,
    test_objects$col_metadata,
    test_objects$row_metadata
  )

  # Test get_annotated_data
  combined_data <- annotated_data$get_annotated_data()

  # Should combine row_metadata with data
  expect_equal(nrow(combined_data), 10)
  expect_equal(ncol(combined_data), ncol(test_objects$row_metadata) + ncol(test_objects$data))  # 2 metadata cols + 20 data cols
  expect_true(all(colnames(test_objects$row_metadata) %in% colnames(combined_data)))
  expect_true(all(colnames(test_objects$data) %in% colnames(combined_data)))
})

test_that("AnnotatedData create_filtered_view works correctly", {
  annotated_data <- AnnotatedData$new(
    test_objects$data,
    test_objects$col_metadata,
    test_objects$row_metadata
  )

  # Test filtering for Breast cancer RNA-seq samples
  filtered_view <- annotated_data$create_filtered_view("Breast", "RNA-seq")

  expect_true(inherits(filtered_view, "AnnotatedData"))
  expect_equal(ncol(filtered_view$data), 5)  # Should have 5 Breast RNA-seq samples
  expect_equal(nrow(filtered_view$col_metadata), 5)
  expect_true(all(filtered_view$col_metadata$`Cancer Type` == "Breast"))
  expect_true(all(filtered_view$col_metadata$assay == "RNA-seq"))

  # Test with multiple cancer types
  multi_filtered <- annotated_data$create_filtered_view(c("Breast", "Lung"), "RNA-seq")
  expect_equal(ncol(multi_filtered$data), 10)  # Should have 10 RNA-seq samples
})

test_that("AnnotatedData filtering handles edge cases", {
  annotated_data <- AnnotatedData$new(
    test_objects$data,
    test_objects$col_metadata,
    test_objects$row_metadata
  )

  # Test filtering with non-existent cancer type
  expect_error(
    annotated_data$create_filtered_view("Pancreatic", "RNA-seq"),
    "No matching samples found for the given cancer_type and assay"
  )

  # Test filtering with non-existent assay
  expect_error(
    annotated_data$create_filtered_view("Breast", "Metabolomics"),
    "No matching samples found for the given cancer_type and assay"
  )
})

# Mock normalization functions for testing
mock_median_normalization <- function(data) {
  return(data * 2)  # Simple mock transformation
}

# Test normalization parameter passing
test_that("AnnotatedData normalization parameter works", {
  annotated_data <- AnnotatedData$new(
    test_objects$data,
    test_objects$col_metadata,
    test_objects$row_metadata
  )

  # Test with "none" normalization (should not error)
  filtered_view <- annotated_data$create_filtered_view("Breast", "RNA-seq", normalization = "none")
  expect_true(inherits(filtered_view, "AnnotatedData"))

  # Test with invalid normalization method
  expect_error(
    annotated_data$create_filtered_view("Breast", "RNA-seq", normalization = "invalid_method"),
    "Unknown normalization method: invalid_method"
  )
})

# Run all tests
cat("Running AnnotatedData test suite...\n")
