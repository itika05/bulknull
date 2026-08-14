# Generate NUMBERS_TABLE.tsv from config/dilution_parameters.csv
# This ensures all numbers in the manuscript have documented provenance

config_file <- "config/dilution_parameters.csv"

if (!file.exists(config_file)) {
  stop("ERROR: config/dilution_parameters.csv not found. Run from project root.")
}

# Read parameters (keep as character strings to preserve precision/notation)
params <- read.csv(config_file, stringsAsFactors = FALSE, colClasses = "character")

# Extract unique values as character strings
unique_values <- unique(params$value)

# Sort numerically but preserve string format
numeric_sort_order <- order(suppressWarnings(as.numeric(unique_values)))
sorted_values <- unique_values[numeric_sort_order]

# Create output table
numbers_table <- data.frame(value = sorted_values, stringsAsFactors = FALSE)
rownames(numbers_table) <- NULL

# Write to project root for use in check_orphan_numbers.R
output_file <- "NUMBERS_TABLE.tsv"
# Use write.csv with explicit settings to avoid scientific notation
write.csv(numbers_table, file = output_file, quote = FALSE, row.names = FALSE)

cat(sprintf("Generated %s with %d unique values\n", output_file, nrow(numbers_table)))
cat("Values sourced from config/dilution_parameters.csv:\n")
cat("  - Derived: qnorm critical values\n")
cat("  - Thresholds: DDS/DI interpretation cutoffs, applicability gate\n")
cat("  - Parameters: synthetic validation design, portfolio analysis ranges\n")
cat("  - Observations: DPP9/M8 case study measurements\n")
cat("  - Results: synthetic validation outcomes\n")
cat("  - Package: test tolerances, assertion counts\n")
