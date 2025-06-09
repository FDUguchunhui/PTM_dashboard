# Statistical Analysis Function - Fixes and Tests

## Bugs Fixed in `test.R`

### 1. **Incorrect Function Name**
- **Issue**: Used `wilcoxon.test()` instead of `wilcox.test()`
- **Fix**: Changed to `wilcox.test()` which is the correct R function name

### 2. **Missing Input Validation**
- **Issue**: No validation for null inputs, empty groups, or invalid parameters
- **Fix**: Added comprehensive input validation:
  - Check for NULL `annotated_data`
  - Validate non-empty group types
  - Validate `test_type` parameter
  - Validate `pvalue_threshold` range (0 < p < 1)
  - Check that `annotated_data` has required methods

### 3. **No Error Handling for Data Retrieval**
- **Issue**: No error handling if `get_data()` fails
- **Fix**: Added `tryCatch()` around data retrieval with informative error messages

### 4. **Missing Data Validation**
- **Issue**: No checks for retrieved data validity
- **Fix**: Added validation for:
  - Non-null data retrieval
  - Non-empty datasets
  - Matching number of features between groups
  - Presence of intensity columns after removing metadata

### 5. **AUC Calculation Edge Case**
- **Issue**: pROC fails when all values are identical
- **Fix**: Added check for identical values and return AUC = 0.5 in such cases

### 6. **Inefficient Metadata Handling**
- **Issue**: `get_row_metadata()` called multiple times
- **Fix**: Call once and reuse the result

## Test Coverage

The test file `test_statistical_analysis.R` includes:

### Mock Classes
- `MockAnnotatedData`: Simulates the real AnnotatedData object
- `ErrorMockAnnotatedData`: Tests error conditions
- Helper functions for creating test data

### Test Scenarios
1. **Basic Functionality**
   - Normal t-test operation
   - Normal Wilcoxon test operation
   - Correct output structure

2. **Input Validation**
   - NULL inputs
   - Empty group types
   - Invalid test types
   - Invalid p-value thresholds
   - Missing methods in annotated_data

3. **Error Handling**
   - Data retrieval failures
   - Method availability checks

4. **Edge Cases**
   - Data with NAs and zeros
   - Identical values (AUC edge case)
   - Empty groups after filtering

5. **Result Validation**
   - Significance calculation
   - P-value sorting
   - AUC calculation

## Running the Tests

### Option 1: Run the test script directly
```r
source("run_tests.R")
```

### Option 2: Run tests manually
```r
# Install required packages if needed
install.packages(c("testthat", "pROC"))

# Load libraries
library(testthat)
library(pROC)

# Run tests
source("test_statistical_analysis.R")
```

### Option 3: Using testthat directly
```r
library(testthat)
test_file("test_statistical_analysis.R")
```

## Dependencies

- `testthat`: For unit testing framework
- `pROC`: For AUC calculation (already required by main function)

## Function Usage

```r
# Example usage of the fixed function
result <- perform_statistical_analysis(
  annotated_data = your_annotated_data_object,
  group1_types = c("Cancer_Type_A"),
  group2_types = c("Cancer_Type_B"), 
  assay = "FT",
  normalization = "none",
  test_type = "t_test",
  pvalue_threshold = 0.05
)
```

## Output Structure

The function returns a data.frame with the following columns:
- Row metadata columns (from `get_row_metadata()`)
- `log2_fold_change`: Log2 ratio of group1 mean / group2 mean
- `p_value`: Statistical test p-value
- `mean_group1`, `mean_group2`: Group means
- `n_group1`, `n_group2`: Number of valid samples per group
- `num_non_NA_1`, `num_non_NA_2`: Number of non-NA values before filtering
- `auc`: Area under the ROC curve
- `significant`: Boolean indicator based on p-value threshold 