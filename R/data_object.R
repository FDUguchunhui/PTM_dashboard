library(R6)

AnnotatedData <- R6Class(
  "AnnotatedData",

  public = list(
    data = NULL,
    col_metadata = NULL,
    row_metadata = NULL,

    initialize = function(data, col_metadata, row_metadata) {
      if (!is.matrix(data)) stop("Data must be a matrix")
      if (!is.data.frame(col_metadata)) stop("col_metadata must be a data.frame")
      if (!is.data.frame(row_metadata)) stop("row_metadata must be a data.frame")

      if (ncol(data) != nrow(col_metadata)) stop("Number of columns in data must match rows in col_metadata")
      if (nrow(data) != nrow(row_metadata)) stop("Number of rows in data must match rows in row_metadata")

      self$data <- data
      self$col_metadata <- col_metadata
      self$row_metadata <- row_metadata
    },

    # Accessor Methods
    get_data = function(cancer_type = NULL, assay = NULL, normalization = "none") {

      if (is.null(cancer_type) || is.null(assay)) {
        stop("Error: Both 'cancer_type' and 'assay' must be provided.")
      }

      # Filter column metadata for selected cancer type and assay
      selected_columns <- self$col_metadata[
        self$col_metadata$`Cancer Type` %in% cancer_type & self$col_metadata$assay %in% assay,
        "Run"] %>% pull()

      # Check if selected_columns is empty
      if (length(selected_columns) == 0) {
        stop("No matching samples found for the given cancer_type and assay.")
      }

      # Subset the data matrix using selected columns
      selected_data <- self$data[, selected_columns, drop = FALSE]

      # Apply normalization if specified
      if (normalization == "median") {
        # browser()
        normalized_data <- median_normalization(selected_data)  # Normalize only selected data
      }

      else if (normalization == 'plate_median') {
        # browser()
        normalized_data <- batch_effect_norm_by_assay(selected_data, self$col_metadata[self$col_metadata$Run %in% selected_columns, ])
      }

      else {
        normalized_data <- selected_data
      }

      # Combine with row metadata before returning
      return(self$combine_with_row_metadata(normalized_data))
    },

    get_col_metadata = function() {
      return(self$col_metadata)
    },

    get_row_metadata = function() {
      return(self$row_metadata)
    },

    combine_with_row_metadata = function(data) {
      return(cbind(self$row_metadata, as.data.frame(data)))
    }
  )
)
