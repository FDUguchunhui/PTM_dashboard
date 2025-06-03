library(R6)
library(dplyr)

#' AnnotatedData: A Class for Managing Biological Data with Metadata
#'
#' @title AnnotatedData R6 Class
#' @description
#' An R6 class designed for handling biological/medical data matrices with associated metadata.
#' Provides functionality for data validation, filtering, normalization, and creating processed 
#' views of the data while preserving original data integrity.
#'
#' @details
#' The AnnotatedData class follows an immutable design pattern where the original data is 
#' preserved privately, and operations create new views or instances. This ensures data 
#' integrity and allows for multiple analysis workflows from the same source data.
#'
#' The class supports:
#' \itemize{
#'   \item Data validation and integrity checks
#'   \item Sample filtering by cancer type and assay
#'   \item Multiple normalization methods
#'   \item Metadata management for both samples and features
#'   \item Creation of filtered views without modifying original data
#' }
#'
#' @section Required Metadata Columns:
#' The col_metadata parameter must contain the following columns:
#' \itemize{
#'   \item \code{Run}: Sample identifier (character)
#'   \item \code{Cancer Type}: Type of cancer (character)
#'   \item \code{assay}: Assay type (character)
#' }
#'
#' @section Normalization Methods:
#' \describe{
#'   \item{\code{"none"}}{No normalization applied}
#'   \item{\code{"median"}}{Median normalization using median_normalization() function}
#'   \item{\code{"plate_median"}}{Batch effect normalization by assay using batch_effect_norm_by_assay() function}
#' }
#'
#' @examples
#' \dontrun{
#' # Create sample data
#' data_matrix <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' colnames(data_matrix) <- paste0("Sample_", 1:10)
#' rownames(data_matrix) <- paste0("Feature_", 1:100)
#'
#' # Create metadata
#' col_metadata <- data.frame(
#'   Run = paste0("Sample_", 1:10),
#'   `Cancer Type` = rep(c("Breast", "Lung"), each = 5),
#'   assay = rep(c("RNA-seq", "Proteomics"), 5),
#'   stringsAsFactors = FALSE
#' )
#'
#' row_metadata <- data.frame(
#'   Feature_ID = paste0("Feature_", 1:100),
#'   Gene_Name = paste0("Gene_", 1:100),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Create AnnotatedData object
#' annotated_data <- AnnotatedData$new(
#'   data = data_matrix,
#'   col_metadata = col_metadata,
#'   row_metadata = row_metadata
#' )
#'
#' # Explore the data
#' info <- annotated_data$get_info()
#' counts <- annotated_data$get_sample_counts()
#'
#' # Create filtered view
#' breast_rna <- annotated_data$create_filtered_view(
#'   cancer_type = "Breast",
#'   assay = "RNA-seq",
#'   normalization = "median"
#' )
#' }
#'
#' @importFrom R6 R6Class
#' @importFrom dplyr group_by summarise
#' @export
AnnotatedData <- R6Class(
  "AnnotatedData",

  private = list(
    # Private fields to store original data
    .original_data = NULL,
    .original_col_metadata = NULL,
    .original_row_metadata = NULL,
    
    #' @description Validate constructor inputs
    #' @param data A numeric matrix
    #' @param col_metadata A data.frame with sample metadata
    #' @param row_metadata A data.frame with feature metadata
    #' @return TRUE if validation passes, throws error otherwise
    validate_inputs = function(data, col_metadata, row_metadata) {
      if (!is.matrix(data)) stop("Data must be a matrix")
      if (!is.data.frame(col_metadata)) stop("col_metadata must be a data.frame")
      if (!is.data.frame(row_metadata)) stop("row_metadata must be a data.frame")
      
      # Check for required columns in col_metadata
      required_col_cols <- c("Run", "Cancer Type", "assay")
      missing_cols <- setdiff(required_col_cols, colnames(col_metadata))
      if (length(missing_cols) > 0) {
        stop(paste("col_metadata is missing required columns:", paste(missing_cols, collapse = ", ")))
      }
      
      if (ncol(data) != nrow(col_metadata)) {
        stop("Number of columns in data must match rows in col_metadata")
      }
      if (nrow(data) != nrow(row_metadata)) {
        stop("Number of rows in data must match rows in row_metadata")
      }
      
      return(TRUE)
    },
    
    #' @description Filter samples based on cancer type and assay
    #' @param cancer_type Vector of cancer types to include
    #' @param assay Vector of assay types to include
    #' @return Vector of selected run IDs
    filter_samples = function(cancer_type, assay) {
      # Get runs that match the filtering criteria
      matching_rows <- self$col_metadata[
        self$col_metadata$`Cancer Type` %in% cancer_type & 
        self$col_metadata$assay %in% assay, 
      ]
      
      if (nrow(matching_rows) == 0) {
        stop("No matching samples found for the given cancer_type and assay.")
      }
      
      selected_runs <- matching_rows$Run
      return(selected_runs)
    },
    
    #' @description Subset data matrix to selected runs
    #' @param selected_runs Vector of run IDs to include
    #' @return Filtered data matrix
    subset_data = function(selected_runs) {
      # Subset the data matrix using selected runs
      if (!all(selected_runs %in% colnames(self$data))) {
        missing_runs <- setdiff(selected_runs, colnames(self$data))
        stop(paste("The following runs are not found in data columns:", 
                  paste(missing_runs, collapse = ", ")))
      }
      
      return(self$data[, selected_runs, drop = FALSE])
    },
    
    #' @description Apply normalization to data subset
    #' @param data_subset Matrix to normalize
    #' @param normalization Normalization method
    #' @param selected_runs Vector of run IDs for the subset
    #' @return Normalized data matrix
    apply_normalization = function(data_subset, normalization, selected_runs) {
      if (normalization == "none") {
        return(data_subset)
      }
      
      if (normalization == "median") {
        return(median_normalization(data_subset))
      }
      
      if (normalization == "plate_median") {
        # Get metadata for selected samples only
        sample_metadata <- self$col_metadata[
          self$col_metadata$Run %in% selected_runs, 
        ]
        return(batch_effect_norm_by_assay(data_subset, sample_metadata))
      }
      
      stop(paste("Unknown normalization method:", normalization))
    },
    
    #' @description Combine processed data with row metadata
    #' @param processed_data Processed data matrix
    #' @return Combined data.frame
    combine_with_metadata = function(processed_data) {
      return(cbind(self$row_metadata, as.data.frame(processed_data)))
    }
  ),

  public = list(
    #' @field data Current data matrix (may be filtered/processed)
    data = NULL,
    
    #' @field col_metadata Current column metadata (may be filtered)
    col_metadata = NULL,
    
    #' @field row_metadata Current row metadata
    row_metadata = NULL,

    #' @description Initialize a new AnnotatedData object
    #' @param data A numeric matrix where rows are features and columns are samples
    #' @param col_metadata A data.frame with sample metadata containing required columns:
    #'   Run, Cancer Type, and assay
    #' @param row_metadata A data.frame with feature metadata
    #' @return A new AnnotatedData object
    #' @examples
    #' \dontrun{
    #' data_matrix <- matrix(rnorm(100), nrow = 10, ncol = 10)
    #' col_metadata <- data.frame(
    #'   Run = paste0("Sample_", 1:10),
    #'   `Cancer Type` = rep(c("Breast", "Lung"), each = 5),
    #'   assay = rep("RNA-seq", 10)
    #' )
    #' row_metadata <- data.frame(Gene = paste0("Gene_", 1:10))
    #' annotated_data <- AnnotatedData$new(data_matrix, col_metadata, row_metadata)
    #' }
    initialize = function(data, col_metadata, row_metadata) {
      # Validate inputs
      private$validate_inputs(data, col_metadata, row_metadata)
      
      # Store original data privately
      private$.original_data <- data
      private$.original_col_metadata <- col_metadata
      private$.original_row_metadata <- row_metadata
      
      # Set public fields (these represent the current state)
      self$data <- data
      self$col_metadata <- col_metadata
      self$row_metadata <- row_metadata
    },

    #' @description Get annotated data combining row metadata with current data
    #' @return A data.frame combining row metadata with the current data matrix
    #' @examples
    #' \dontrun{
    #' annotated_data <- AnnotatedData$new(data, col_metadata, row_metadata)
    #' combined_data <- annotated_data$get_annotated_data()
    #' }
    get_annotated_data = function() {
      # Process data through the pipeline
      final_data <- private$combine_with_metadata(self$data)
      return(final_data)
    },

    #' @description Get current column metadata
    #' @return A data.frame of column metadata
    get_col_metadata = function() {
      return(self$col_metadata)
    },

    #' @description Get current row metadata
    #' @return A data.frame of row metadata
    get_row_metadata = function() {
      return(self$row_metadata)
    },
    
    #' @description Get original, unprocessed data and metadata
    #' @return A list containing original data, col_metadata, and row_metadata
    #' @examples
    #' \dontrun{
    #' original <- annotated_data$get_original_data()
    #' original_matrix <- original$data
    #' }
    get_original_data = function() {
      return(list(
        data = private$.original_data,
        col_metadata = private$.original_col_metadata,
        row_metadata = private$.original_row_metadata
      ))
    },
    
    #' @description Get summary information about the current data
    #' @return A list containing data dimensions and available categories
    #' @examples
    #' \dontrun{
    #' info <- annotated_data$get_info()
    #' print(paste("Features:", info$n_features))
    #' print(paste("Samples:", info$n_samples))
    #' }
    get_info = function() {
      return(list(
        n_features = nrow(self$data),
        n_samples = ncol(self$data),
        cancer_types = unique(self$col_metadata$`Cancer Type`),
        assays = unique(self$col_metadata$assay),
        sample_ids = colnames(self$data)
      ))
    },
    
    #' @description Check if samples matching criteria exist
    #' @param cancer_type Vector of cancer type(s) to check
    #' @param assay Vector of assay type(s) to check
    #' @return Logical indicating if matching samples exist
    #' @examples
    #' \dontrun{
    #' has_breast <- annotated_data$has_samples("Breast", "RNA-seq")
    #' if (has_breast) {
    #'   # Proceed with analysis
    #' }
    #' }
    has_samples = function(cancer_type, assay) {
      matching_samples <- self$col_metadata[
        self$col_metadata$`Cancer Type` %in% cancer_type & 
        self$col_metadata$assay %in% assay, 
      ]
      return(nrow(matching_samples) > 0)
    },
    
    #' @description Get sample counts by cancer type and assay
    #' @return A data.frame with columns Cancer Type, assay, and count
    #' @examples
    #' \dontrun{
    #' counts <- annotated_data$get_sample_counts()
    #' print(counts)
    #' }
    get_sample_counts = function() {
      return(
        self$col_metadata %>%
          group_by(`Cancer Type`, assay) %>%
          summarise(count = n(), .groups = "drop") %>%
          as.data.frame()
      )
    },
    
    #' @description Create a filtered view with optional normalization
    #' @param cancer_type Vector of cancer type(s) to include
    #' @param assay Vector of assay type(s) to include
    #' @param normalization Normalization method: "none", "median", or "plate_median"
    #' @param missing_fill Value to fill missing data with. Default is NA (no filling)
    #' @return A new AnnotatedData instance with filtered and optionally normalized data
    #' @examples
    #' \dontrun{
    #' # Create filtered view for Breast cancer RNA-seq samples
    #' breast_rna <- annotated_data$create_filtered_view(
    #'   cancer_type = "Breast",
    #'   assay = "RNA-seq",
    #'   normalization = "median"
    #' )
    #' 
    #' # Multiple criteria with missing value filling
    #' multi_filtered <- annotated_data$create_filtered_view(
    #'   cancer_type = c("Breast", "Lung"),
    #'   assay = c("RNA-seq", "Proteomics"),
    #'   normalization = "plate_median",
    #'   missing_fill = 0
    #' )
    #' }
    create_filtered_view = function(cancer_type, assay, normalization = "none", missing_fill = NA) {
      selected_runs <- private$filter_samples(cancer_type, assay)
      data_subset <- private$subset_data(selected_runs)
      normalized_data <- private$apply_normalization(data_subset, normalization, selected_runs)
      
      # Fill missing values if missing_fill is not NA
      if (!is.na(missing_fill)) {
        normalized_data[is.na(normalized_data)] <- missing_fill
      }
      
      # Filter col_metadata for selected runs
      filtered_col_metadata <- self$col_metadata[
        self$col_metadata$Run %in% selected_runs, 
      ]
      
      # Create new instance with corrected variable names
      AnnotatedData$new(
        data = normalized_data,
        col_metadata = filtered_col_metadata,
        row_metadata = self$row_metadata
      )
    }
  )
) 