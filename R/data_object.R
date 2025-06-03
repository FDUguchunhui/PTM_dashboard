library(R6)
library(dplyr)

AnnotatedData <- R6Class(
  "AnnotatedData",

  private = list(
    # Private fields to store original data
    .original_data = NULL,
    .original_col_metadata = NULL,
    .original_row_metadata = NULL,
    
    # Private methods for data processing
    validate_inputs = function(data, col_metadata, row_metadata) {
      if (!is.matrix(data)) stop("Data must be a matrix")
      if (!is.data.frame(col_metadata)) stop("col_metadata must be a data.frame")
      if (!is.data.frame(row_metadata)) stop("row_metadata must be a data.frame")
      
      if (ncol(data) != nrow(col_metadata)) {
        stop("Number of columns in data must match rows in col_metadata")
      }
      if (nrow(data) != nrow(row_metadata)) {
        stop("Number of rows in data must match rows in row_metadata")
      }
      
      return(TRUE)
    },
    
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
    
    subset_data = function(selected_runs) {
      # Subset the data matrix using selected runs
      if (!all(selected_runs %in% colnames(self$data))) {
        missing_runs <- setdiff(selected_runs, colnames(self$data))
        stop(paste("The following runs are not found in data columns:", 
                  paste(missing_runs, collapse = ", ")))
      }
      
      return(self$data[, selected_runs, drop = FALSE])
    },
    
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
    
    combine_with_metadata = function(processed_data) {
      return(cbind(self$row_metadata, as.data.frame(processed_data)))
    }
  ),

  public = list(
    # Public fields (read-only access to processed data)
    data = NULL,
    col_metadata = NULL,
    row_metadata = NULL,

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

    # Main data retrieval method
    get_data = function(cancer_type = NULL, assay = NULL, normalization = "none") {
      # Validate required parameters
      if (is.null(cancer_type) || is.null(assay)) {
        stop("Error: Both 'cancer_type' and 'assay' must be provided.")
      }
      
      # Validate normalization parameter
      valid_normalizations <- c("none", "median", "plate_median")
      if (!normalization %in% valid_normalizations) {
        stop(paste("Normalization must be one of:", 
                  paste(valid_normalizations, collapse = ", ")))
      }
      
      # Process data through the pipeline
      selected_runs <- private$filter_samples(cancer_type, assay)
      data_subset <- private$subset_data(selected_runs)
      normalized_data <- private$apply_normalization(data_subset, normalization, selected_runs)
      final_data <- private$combine_with_metadata(normalized_data)
      
      return(final_data)
    },

    # Accessor methods for metadata
    get_col_metadata = function() {
      return(self$col_metadata)
    },

    get_row_metadata = function() {
      return(self$row_metadata)
    },
    
    # Method to get original data (useful for debugging or reprocessing)
    get_original_data = function() {
      return(list(
        data = private$.original_data,
        col_metadata = private$.original_col_metadata,
        row_metadata = private$.original_row_metadata
      ))
    },
    
    # Method to get data dimensions and summary info
    get_info = function() {
      return(list(
        n_features = nrow(self$data),
        n_samples = ncol(self$data),
        cancer_types = unique(self$col_metadata$`Cancer Type`),
        assays = unique(self$col_metadata$assay),
        sample_ids = colnames(self$data)
      ))
    },
    
    # Method to check if specific samples exist
    has_samples = function(cancer_type, assay) {
      matching_samples <- self$col_metadata[
        self$col_metadata$`Cancer Type` %in% cancer_type & 
        self$col_metadata$assay %in% assay, 
      ]
      return(nrow(matching_samples) > 0)
    },
    
    # Method to get available sample counts
    get_sample_counts = function() {
      return(
        self$col_metadata %>%
          group_by(`Cancer Type`, assay) %>%
          summarise(count = n(), .groups = "drop") %>%
          as.data.frame()
      )
    }
  )
)
