#!/usr/bin/env Rscript

# Load the data.table package
library(data.table)

# Main execution block
args <- commandArgs(trailingOnly = TRUE)

# Check if arguments are provided
if (length(args) != 3) {
  stop("Usage: Rscript script.R <input_file> <output_file> <colnames>")
}

input_file <- args[1]
output_file <- args[2]

# read in data as data table
data <- fread(input_file, sep = "\t")

# Split the first column based on '_'
split_columns <- strsplit(as.character(data[[1]]), "_")

# Create a new data frame with the split columns
split_df <- do.call(rbind, lapply(split_columns, function(x) {
  length(x) <- 4  # Ensure there are always four columns, filling with NA if necessary
  return(x)
}))

# Convert the split_df to a data table
split_df <- as.data.table(split_df)

# Rename the columns using the user input
setnames(split_df, unlist(strsplit(args[3], "_")))

# Combine the new split columns with the original data (dropping the first column)
data <- data[, -1, with = FALSE]  # Remove the original first column
data <- cbind(split_df, data)  # Combine the split columns with the remaining data

# Save the modified data table (optional)
fwrite(data, output_file, sep = "\t", row.names = FALSE) 