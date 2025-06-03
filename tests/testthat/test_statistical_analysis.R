# Test file for perform_statistical_analysis function
# This file contains unit tests and a mock AnnotatedData class

# Load required libraries
library(testthat)

# Create the test data once using the shared helper function
test_objects <- create_test_data()

# create AnnotationData
annotated_data <- AnnotatedData$new(
  test_objects$data,
  test_objects$col_metadata,
  test_objects$row_metadata
)

# Run tests
test_that("Basic functionality works", {

  result <- perform_statistical_analysis(
    annotated_data = annotated_data,
    group1_types = "Breast",
    group2_types = "Lung",
    assay = "FT",
    normalization = "none",
    test_type = "t_test"
  )

  expect_s3_class(result, "data.frame")
  expect_true("log2_fold_change" %in% colnames(result))
  expect_true("p_value" %in% colnames(result))
  expect_true("auc" %in% colnames(result))
  expect_true("significant" %in% colnames(result))
  expect_equal(nrow(result), 10)
})

test_that("Wilcoxon test works", {

  result <- perform_statistical_analysis(
    annotated_data = annotated_data,
    group1_types = "Breast",
    group2_types = "Lung",
    assay = "FT",
    normalization = "none",
    test_type = "wilcoxon"
  )

  expect_s3_class(result, "data.frame")
  expect_true(any(!is.na(result$p_value)))
})

test_that("Input validation works", {

  # Test NULL annotated_data
  expect_error(
    perform_statistical_analysis(NULL, "Breast", "Lung", "FT", "none"),
    "annotated_data cannot be NULL"
  )

  # Test empty group types
  expect_error(
    perform_statistical_analysis(annotated_data, character(0), "Lung", "FT", "none"),
    "Both group1_types and group2_types must contain at least one cancer type"
  )

  expect_error(
    perform_statistical_analysis(annotated_data, "Breast", character(0), "FT", "none"),
    "Both group1_types and group2_types must contain at least one cancer type"
  )

  # Test invalid test_type
  expect_error(
    perform_statistical_analysis(annotated_data, "Lung", "Breast", "FT", "none", test_type = "invalid"),
    "test_type must be either 't_test' or 'wilcoxon'"
  )

  # Test invalid p-value threshold
  expect_error(
    perform_statistical_analysis(annotated_data, "Breast", "Lung", "FT", "none", pvalue_threshold = 0),
    "pvalue_threshold must be between 0 and 1"
  )

  expect_error(
    perform_statistical_analysis(annotated_data, "Breast", "Lung", "FT", "none", pvalue_threshold = 1.5),
    "pvalue_threshold must be between 0 and 1"
  )
})


test_that("Handles missing methods in annotated_data", {
  # Create object without required methods
  bad_mock <- list(some_other_method = function() {})

  expect_error(
    perform_statistical_analysis(bad_mock, "Breast", "Lung", "FT", "none"),
    "annotated_data must have get_data\\(\\) and get_row_metadata\\(\\) methods"
  )
})


test_that("Significance calculation works", {

  result <- perform_statistical_analysis(
    annotated_data = annotated_data,
    group1_types = "Breast",
    group2_types = "Lung",
    assay = "FT",
    normalization = "none",
    pvalue_threshold = 0.1
  )

  # Check that significance is calculated correctly
  expected_significant <- result$p_value < 0.1 & !is.na(result$p_value)
  expect_equal(result$significant, expected_significant)
})

test_that("Results are sorted by p-value", {

  result <- perform_statistical_analysis(
    annotated_data = annotated_data,
    group1_types = "Breast",
    group2_types = "Lung",
    assay = "FT",
    normalization = "none"
  )

  # Remove NAs for sorting check
  non_na_pvals <- result$p_value[!is.na(result$p_value)]
  expect_true(all(diff(non_na_pvals) >= 0))
})

