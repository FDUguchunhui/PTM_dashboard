# Helper file for shared test data creation
# This file is automatically loaded before running tests
# Create shared test data for data object tests
create_test_data <- function() {
  test_data <- matrix(rnorm(200), nrow = 10, ncol = 20)
  colnames(test_data) <- paste0("Run_", 1:20)
  rownames(test_data) <- paste0("Gene_", 1:10)

  col_metadata <- data.frame(
    Run = paste0("Run_", 1:20),
    `Cancer Type` = rep(c("Breast", "Lung"), each = 10),
    assay = rep(c("RNA-seq", "Proteomics"), times = 10),
    plate = rep(c("P1", "P2"), each = 10),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  row_metadata <- data.frame(
    Protein = paste0("Protein_", 1:10),
    Gene = paste0("Gene_", 1:10),
    Type = rep(c("Oncogene", "Tumor Suppressor"), times = 5),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  return(list(
    data = test_data,
    col_metadata = col_metadata,
    row_metadata = row_metadata
  ))
}

# Generic function that allows custom assay values
create_test_data <- function(assays = c("FT", "IgB")) {
  test_data <- matrix(rnorm(200), nrow = 10, ncol = 20)
  colnames(test_data) <- paste0("Run_", 1:20)
  rownames(test_data) <- paste0("Gene_", 1:10)

  col_metadata <- data.frame(
    Run = paste0("Run_", 1:20),
    `Cancer Type` = rep(c("Breast", "Lung"), each = 10),
    assay = rep(assays, times = 10),
    plate = rep(c("P1", "P2"), each = 10),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  row_metadata <- data.frame(
    Protein = paste0("Protein_", 1:10),
    Gene = paste0("Gene_", 1:10),
    Type = rep(c("Oncogene", "Tumor Suppressor"), times = 5),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  return(list(
    data = test_data,
    col_metadata = col_metadata,
    row_metadata = row_metadata
  ))
} 