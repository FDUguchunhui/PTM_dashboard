# Test script for the redesigned AnnotatedData class
source('R/data_object.R')
source('R/normalization.R')

# Create sample data for testing
set.seed(42)
n_features <- 100
n_samples <- 20

# Create sample data matrix
sample_data <- matrix(
  rnorm(n_features * n_samples, mean = 1000, sd = 200),
  nrow = n_features,
  ncol = n_samples
)
colnames(sample_data) <- paste0("Sample_", 1:n_samples)

# Create sample row metadata (features)
row_metadata <- data.frame(
  Feature_ID = paste0("Feature_", 1:n_features),
  Gene = paste0("Gene_", sample(1:50, n_features, replace = TRUE)),
  Type = sample(c("Modified", "Unmodified"), n_features, replace = TRUE)
)

# Create sample column metadata (samples)
col_metadata <- data.frame(
  Run = colnames(sample_data),
  `Cancer Type` = rep(c("Breast", "Lung", "Control"), length.out = n_samples),
  assay = rep(c("FT", "IgB"), length.out = n_samples),
  plate = rep(1:4, length.out = n_samples),
  check.names = FALSE
)

print("=== Sample Data Overview ===")
print("Data matrix dimensions:")
print(dim(sample_data))
print("Row metadata:")
print(head(row_metadata))
print("Column metadata:")
print(head(col_metadata))

print("\n=== Pattern 1: Original backward-compatible usage ===")
# Pattern 1: Create instance without filtering (original usage)
annotated_data_full <- AnnotatedData$new(
  data = sample_data,
  col_metadata = col_metadata,
  row_metadata = row_metadata
)

print("Full data instance created successfully")
print(paste("Stored data dimensions:", paste(dim(annotated_data_full$get_data_matrix()), collapse=" x ")))
print(paste("Stored col_metadata rows:", nrow(annotated_data_full$get_col_metadata())))

# Use get_data method to filter (original pattern)
filtered_data_traditional <- annotated_data_full$get_data(
  cancer_type = c("Breast", "Lung"),
  assay = "FT",
  normalization = "none"
)

print("Filtered data using get_data method:")
print(paste("Dimensions:", paste(dim(filtered_data_traditional), collapse=" x ")))
print("Column names (samples):")
print(colnames(filtered_data_traditional)[-(1:ncol(row_metadata))])

print("\n=== Pattern 2: New filtered constructor usage ===")
# Pattern 2: Create instance with filtering during construction (new usage)
annotated_data_filtered <- AnnotatedData$new(
  data = sample_data,
  col_metadata = col_metadata,
  row_metadata = row_metadata,
  cancer_type = c("Breast", "Lung"),
  assay = "FT",
  normalization = "none"
)

print("Filtered data instance created successfully")
print(paste("Stored data dimensions:", paste(dim(annotated_data_filtered$get_data_matrix()), collapse=" x ")))
print(paste("Stored col_metadata rows:", nrow(annotated_data_filtered$get_col_metadata())))

# Get full data from filtered instance
filtered_data_new <- annotated_data_filtered$get_full_data()

print("Full data from filtered instance:")
print(paste("Dimensions:", paste(dim(filtered_data_new), collapse=" x ")))
print("Column names (samples):")
print(colnames(filtered_data_new)[-(1:ncol(row_metadata))])

print("\n=== Verification: Both patterns should give same results ===")
# Compare results from both patterns
traditional_data_cols <- colnames(filtered_data_traditional)[-(1:ncol(row_metadata))]
new_data_cols <- colnames(filtered_data_new)[-(1:ncol(row_metadata))]

print(paste("Traditional pattern sample count:", length(traditional_data_cols)))
print(paste("New pattern sample count:", length(new_data_cols)))
print(paste("Same samples selected:", identical(sort(traditional_data_cols), sort(new_data_cols))))

print("\n=== Dimensional consistency check ===")
print("Filtered instance metadata consistency:")
filtered_instance <- annotated_data_filtered
print(paste("Data matrix columns:", ncol(filtered_instance$get_data_matrix())))
print(paste("Column metadata rows:", nrow(filtered_instance$get_col_metadata())))
print(paste("Row metadata rows:", nrow(filtered_instance$get_row_metadata())))
print(paste("Data matrix rows:", nrow(filtered_instance$get_data_matrix())))

# Check that column names match
data_cols <- colnames(filtered_instance$get_data_matrix())
meta_runs <- filtered_instance$get_col_metadata()$Run
print(paste("Column names match metadata Run:", identical(data_cols, meta_runs)))

print("\n=== Test with different assays ===")
# Test filtering for different assay
igb_data <- AnnotatedData$new(
  data = sample_data,
  col_metadata = col_metadata,
  row_metadata = row_metadata,
  cancer_type = "Control",
  assay = "IgB",
  normalization = "none"
)

print("IgB Control data:")
print(paste("Dimensions:", paste(dim(igb_data$get_full_data()), collapse=" x ")))
print("Selected samples:")
print(igb_data$get_col_metadata()$Run)

print("\n=== Test completed successfully ===") 