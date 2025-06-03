# AnnotatedData Class Documentation

## Overview

The `AnnotatedData` class is an R6 class designed for handling biological/medical data matrices with associated metadata. It provides functionality for data validation, filtering, normalization, and creating processed views of the data.

## Dependencies

```r
library(R6)
library(dplyr)
source("R/normalization.R")
```

## Class Structure

### Constructor

```r
AnnotatedData$new(data, col_metadata, row_metadata)
```

#### Parameters

- **data**: A numeric matrix where rows represent features (e.g., genes, proteins) and columns represent samples
- **col_metadata**: A data.frame containing metadata for each sample (column) with required columns:
  - `Run`: Sample identifier
  - `Cancer Type`: Type of cancer
  - `assay`: Assay type
- **row_metadata**: A data.frame containing metadata for each feature (row)

#### Validation

The constructor validates that:
- `data` is a matrix
- `col_metadata` and `row_metadata` are data.frames
- Required columns exist in `col_metadata`
- Dimensions match between data and metadata

## Public Fields

- **data**: Current data matrix (may be filtered/processed)
- **col_metadata**: Current column metadata (may be filtered)
- **row_metadata**: Current row metadata

## Public Methods

### Core Data Access

#### `get_annotated_data()`
Returns the processed data combined with row metadata.

**Returns**: A data.frame combining row metadata with the current data matrix

#### `get_col_metadata()`
Returns the current column metadata.

**Returns**: A data.frame of column metadata

#### `get_row_metadata()`
Returns the current row metadata.

**Returns**: A data.frame of row metadata

#### `get_original_data()`
Returns the original, unprocessed data and metadata.

**Returns**: A list containing:
- `data`: Original data matrix
- `col_metadata`: Original column metadata
- `row_metadata`: Original row metadata

### Data Information

#### `get_info()`
Returns summary information about the current data.

**Returns**: A list containing:
- `n_features`: Number of features (rows)
- `n_samples`: Number of samples (columns)
- `cancer_types`: Unique cancer types in the data
- `assays`: Unique assay types in the data
- `sample_ids`: Column names (sample IDs)

#### `get_sample_counts()`
Returns counts of samples by cancer type and assay.

**Returns**: A data.frame with columns `Cancer Type`, `assay`, and `count`

#### `has_samples(cancer_type, assay)`
Checks if samples matching the specified criteria exist.

**Parameters**:
- **cancer_type**: Vector of cancer type(s) to check
- **assay**: Vector of assay type(s) to check

**Returns**: Boolean indicating if matching samples exist

### Data Processing

#### `create_filtered_view(cancer_type, assay, normalization = "none")`
Creates a new `AnnotatedData` instance with filtered and optionally normalized data.

**Parameters**:
- **cancer_type**: Vector of cancer type(s) to include
- **assay**: Vector of assay type(s) to include
- **normalization**: Normalization method (default: "none")
  - `"none"`: No normalization
  - `"median"`: Median normalization
  - `"plate_median"`: Batch effect normalization by assay

**Returns**: A new `AnnotatedData` instance with filtered data

## Private Methods

### `validate_inputs(data, col_metadata, row_metadata)`
Validates constructor inputs and ensures data integrity.

### `filter_samples(cancer_type, assay)`
Filters samples based on cancer type and assay criteria.

**Returns**: Vector of selected run IDs

### `subset_data(selected_runs)`
Subsets the data matrix to include only specified runs.

**Returns**: Filtered data matrix

### `apply_normalization(data_subset, normalization, selected_runs)`
Applies the specified normalization method to the data.

**Returns**: Normalized data matrix

### `combine_with_metadata(processed_data)`
Combines processed data with row metadata.

**Returns**: Combined data.frame

## Usage Examples

### Basic Usage

```r
# Load required libraries
library(R6)
library(dplyr)
source("R/normalization.R")
source("R/data_object.R")

# Create sample data
data_matrix <- matrix(rnorm(1000), nrow = 100, ncol = 10)
colnames(data_matrix) <- paste0("Sample_", 1:10)
rownames(data_matrix) <- paste0("Feature_", 1:100)

# Create metadata
col_metadata <- data.frame(
  Run = paste0("Sample_", 1:10),
  `Cancer Type` = rep(c("Breast", "Lung"), each = 5),
  assay = rep(c("RNA-seq", "Proteomics"), 5),
  stringsAsFactors = FALSE
)

row_metadata <- data.frame(
  Feature_ID = paste0("Feature_", 1:100),
  Gene_Name = paste0("Gene_", 1:100),
  stringsAsFactors = FALSE
)

# Create AnnotatedData object
annotated_data <- AnnotatedData$new(
  data = data_matrix,
  col_metadata = col_metadata,
  row_metadata = row_metadata
)
```

### Data Exploration

```r
# Get basic information
info <- annotated_data$get_info()
print(info)

# Check sample counts
counts <- annotated_data$get_sample_counts()
print(counts)

# Check if specific samples exist
has_breast_rna <- annotated_data$has_samples("Breast", "RNA-seq")
print(has_breast_rna)
```

### Data Filtering and Normalization

```r
# Create filtered view for Breast cancer RNA-seq samples with median normalization
breast_rna <- annotated_data$create_filtered_view(
  cancer_type = "Breast",
  assay = "RNA-seq",
  normalization = "median"
)

# Get the processed data
processed_data <- breast_rna$get_annotated_data()

# Get information about the filtered data
filtered_info <- breast_rna$get_info()
print(filtered_info)
```

### Multiple Filter Criteria

```r
# Filter for multiple cancer types and assays
multi_filtered <- annotated_data$create_filtered_view(
  cancer_type = c("Breast", "Lung"),
  assay = c("RNA-seq", "Proteomics"),
  normalization = "plate_median"
)
```

## Error Handling

The class includes comprehensive error checking:

- **Data validation**: Ensures input data types are correct
- **Dimension matching**: Verifies data and metadata dimensions align
- **Required columns**: Checks for required metadata columns
- **Sample existence**: Validates that filtered samples exist
- **Normalization methods**: Ensures valid normalization methods are used

## Normalization Methods

The class supports three normalization methods:

1. **"none"**: No normalization applied
2. **"median"**: Median normalization using `median_normalization()` function
3. **"plate_median"**: Batch effect normalization by assay using `batch_effect_norm_by_assay()` function

Note: The normalization functions must be available from `R/normalization.R`.

## Design Patterns

- **Immutability**: Original data is preserved privately; operations create new views
- **Validation**: Comprehensive input validation at construction
- **Encapsulation**: Private methods handle internal data processing
- **Fluent Interface**: Methods return objects that can be chained or used independently

## Best Practices

1. Always validate your input data before creating an `AnnotatedData` object
2. Use `get_info()` and `get_sample_counts()` to explore your data
3. Check sample existence with `has_samples()` before filtering
4. Create filtered views for specific analyses rather than modifying the original object
5. Use appropriate normalization methods based on your data type and analysis needs 