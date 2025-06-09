#!/usr/bin/env Rscript

# Script to generate R documentation from roxygen2 comments
# This script will create .Rd files in the man/ directory

# Install roxygen2 if not already installed
if (!require(roxygen2, quietly = TRUE)) {
  install.packages("roxygen2")
  library(roxygen2)
}

# Create man directory if it doesn't exist
if (!dir.exists("man")) {
  dir.create("man")
}

# Generate documentation
cat("Generating roxygen2 documentation...\n")

# Document the package
roxygen2::roxygenise()

cat("Documentation generated successfully!\n")
cat("Check the man/ directory for .Rd files\n")

# Optionally, you can also generate a NAMESPACE file
cat("NAMESPACE file updated\n") 