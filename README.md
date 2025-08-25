# DIANN Analysis Platform

A comprehensive platform for biological data management and analysis, featuring both R package functionality for data processing and interactive dashboards for visualization. This project combines automated protein searching pipelines with robust data management tools for biological and medical research.

## 🚀 Overview

This platform provides:

- **R Package (`AnnotatedDataPackage`)**: Advanced biological data management with metadata handling
- **Interactive Dashboards**: Web-based visualization tools for data exploration
- **Containerized Environment**: Docker support for consistent deployment
- **DIA-NN Integration**: Automated protein search pipeline integration

## 📋 Table of Contents

- [Features](#features)
- [Installation](#installation)
- [R Package Usage](#r-package-usage)
- [Dashboard](#dashboard)
- [Docker Deployment](#docker-deployment)
- [Project Structure](#project-structure)
- [Documentation](#documentation)
- [Testing](#testing)
- [Contributing](#contributing)
- [License](#license)

## ✨ Features

### R Package Features
- **Data Validation**: Robust input validation for biological data matrices
- **Metadata Management**: Comprehensive metadata handling for samples and features
- **Data Filtering**: Advanced filtering capabilities by cancer type, assay, and other criteria
- **Normalization**: Multiple normalization methods including median and batch effect correction
- **Data Views**: Create processed views while preserving original data integrity
- **Type Safety**: R6 class-based architecture for reliable data operations

### Dashboard Features
- **Interactive Visualization**: Web-based dashboards for data exploration
- **Real-time Analysis**: Dynamic data processing and visualization
- **Multi-platform Support**: Cross-platform compatibility

### Infrastructure Features
- **Containerization**: Docker support for easy deployment
- **DIA-NN Integration**: Automated protein search pipeline
- **Cross-language Support**: R and Python components working together

## 🛠 Installation

### Prerequisites

- **R** (>= 3.5.0)
- **Docker** (for containerized deployment)
- **Git** (for version control)

### R Package Installation

```r
# Install dependencies
install.packages(c("R6", "dplyr", "testthat", "knitr", "rmarkdown"))

# Clone the repository
# git clone <repository-url>

# Install the package
devtools::install()
# or simulate installation
devtools::load_all()
```

### Docker Installation

```bash
# Build the Docker image
docker build -t diann-analysis .

# Run the container
docker run -it diann-analysis
```

### Development Setup

```bash
# Clone the repository
git clone <repository-url>
cd diann

# Set up R environment (if using renv)
R -e "renv::restore()"
```

## 📖 R Package Usage

### Basic Setup

```r
library(AnnotatedDataPackage)
library(R6)
library(dplyr)

# Load your data
data_matrix <- your_data_matrix
col_metadata <- your_column_metadata  # Must include Run, Cancer Type, assay columns
row_metadata <- your_row_metadata

# Create AnnotatedData object
annotated_data <- AnnotatedData$new(
  data = data_matrix,
  col_metadata = col_metadata,
  row_metadata = row_metadata
)
```

### Data Exploration

```r
# Get basic information about your data
info <- annotated_data$get_info()
print(info)
# Returns: n_features, n_samples, cancer_types, assays, sample_ids

# Check sample distribution
sample_counts <- annotated_data$get_sample_counts()
print(sample_counts)

# Verify specific sample types exist
has_samples <- annotated_data$has_samples(
  cancer_type = "Breast", 
  assay = "RNA-seq"
)
```

### Data Processing

```r
# Create filtered and normalized views
breast_rna <- annotated_data$create_filtered_view(
  cancer_type = "Breast",
  assay = "RNA-seq",
  normalization = "median"
)

# Multiple cancer types and assays
multi_view <- annotated_data$create_filtered_view(
  cancer_type = c("Breast", "Lung"),
  assay = c("RNA-seq", "Proteomics"),
  normalization = "plate_median"
)

# Access processed data
processed_data <- breast_rna$get_annotated_data()
metadata <- breast_rna$get_col_metadata()
```

### Available Normalization Methods

- **`"none"`**: No normalization applied
- **`"median"`**: Median normalization across samples
- **`"plate_median"`**: Batch effect normalization by assay type

## 🎛 Dashboard

The project includes interactive dashboards built with Shiny for R-based visualization:

### Running the Dashboard

```r
# Navigate to dashboard directory
setwd("diann_dashboard")

# Run the Shiny application
shiny::runApp()
```

### Dashboard Features

- **Data Upload**: Upload and validate biological datasets
- **Interactive Filtering**: Dynamic filtering by cancer type and assay
- **Normalization Options**: Real-time normalization with visual feedback
- **Export Capabilities**: Download processed data and visualizations

### License Compliance

⚠️ **Important**: Users are responsible for ensuring compliance with DIA-NN license terms when using or distributing the Docker image.

## 📁 Project Structure

```
diann/
├── R/                          # R package source code
│   ├── data_object.R          # Core AnnotatedData class
│   ├── normalization.R        # Normalization functions
│   ├── visualization.R        # Plotting utilities
│   ├── helper_functions.R     # Utility functions
│   └── PTM_definition.R       # Post-translational modification definitions
├── diann_dashboard/           # Shiny dashboard application
│   ├── app.R                  # Main application entry point
│   ├── ui.R                   # User interface definition
│   └── server.R               # Server logic
├── docs/                      # Documentation
│   └── AnnotatedData_Documentation.md
├── tests/                     # Test suite
│   ├── testthat.R
│   └── testthat/
├── data/                      # Sample/test data
├── output/                    # Generated outputs
├── rmarkdown/                 # R Markdown documents
├── Dockerfile                 # Container configuration
├── DESCRIPTION               # R package metadata
├── NAMESPACE                 # R package namespace
├── pyproject.toml           # Python project configuration
└── README.md                # This file
```

## 📚 Documentation

### Available Documentation

- **[Complete API Documentation](docs/AnnotatedData_Documentation.md)**: Comprehensive guide to the AnnotatedData class
- **R Package Help**: Use `?AnnotatedData` in R for built-in help
- **Function Documentation**: Individual function documentation via `?function_name`

### Generating Documentation

```r
# Generate package documentation
devtools::document()

# Build package manual
devtools::build_manual()
```

## 🧪 Testing

The project includes a comprehensive test suite using `testthat`:

### Running Tests

```r
# Run all tests
devtools::test()

# Run specific test files
testthat::test_file("tests/testthat/test-data-object.R")

# Check package
devtools::check()
```

### Test Coverage

Tests cover:
- Data validation and error handling
- Filtering and subsetting operations
- Normalization methods
- Metadata management
- Edge cases and boundary conditions

## 🤝 Contributing

We welcome contributions! Please follow these guidelines:

1. **Fork the repository**
2. **Create a feature branch**: `git checkout -b feature/your-feature`
3. **Write tests** for new functionality
4. **Ensure documentation** is updated
5. **Submit a pull request**

### Development Guidelines

- Follow R package development best practices
- Use roxygen2 for documentation
- Maintain test coverage above 80%
- Follow consistent coding style

## 📄 License

This project is licensed under the MIT License. See the `DESCRIPTION` file for details.

### Third-party Licenses

- **DIA-NN**: Users must comply with DIA-NN licensing terms
- **R Dependencies**: See individual package licenses

## 🆘 Support

For questions, issues, or contributions:

1. **Check existing issues** in the GitHub repository
2. **Create a new issue** with detailed information
3. **Include reproducible examples** when reporting bugs

## 🔄 Version History

- **v1.0.0**: Initial release with core AnnotatedData functionality
- **v0.1.0**: Development version with dashboard integration

---

**Note**: This project combines R and Python components for comprehensive biological data analysis. The focus is on providing robust, well-tested tools for researchers working with biological and medical datasets.
