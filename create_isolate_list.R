#!/usr/bin/env Rscript
# Script to extract isolates from metadata and create file paths for pipeline
# This script can be scaled from 10 isolates to all 2000

library(readr)
library(dplyr)

# Read the metadata CSV file
cat("Reading metadata file...\n")
metadata <- read_csv("../final_merged_data.csv")

# Get unique isolates (remove duplicates from paired-end files)
unique_isolates <- metadata %>%
  select(Main_isolate_name, Path_on_legserv) %>%
  distinct(Main_isolate_name, .keep_all = TRUE) %>%
  arrange(Main_isolate_name)

cat("Total isolates found:", nrow(unique_isolates), "\n")

# Extract first 10 isolates for testing
n_isolates <- 10  # Change this number to scale up
test_isolates <- unique_isolates[1:n_isolates, ]

cat("Selected", n_isolates, "isolates for testing:\n")
print(test_isolates$Main_isolate_name)

# For each isolate, find the paired-end files
isolate_files <- metadata %>%
  filter(Main_isolate_name %in% test_isolates$Main_isolate_name) %>%
  select(Main_isolate_name, File_name, Path_on_legserv) %>%
  arrange(Main_isolate_name, File_name)

# Group by isolate to create file pairs
isolate_pairs <- isolate_files %>%
  group_by(Main_isolate_name) %>%
  summarise(
    file_paths = paste(Path_on_legserv, collapse = " "),
    file_count = n(),
    .groups = 'drop'
  )

cat("\nFile pairs found:\n")
print(isolate_pairs)

# Create the file path pattern for nextflow
# We need to create a pattern that matches the paired-end files
base_dirs <- unique(dirname(isolate_files$Path_on_legserv))
cat("\nBase directories found:\n")
print(base_dirs)

# Save results for pipeline use
write_csv(test_isolates, "selected_isolates.csv")
write_csv(isolate_files, "isolate_file_paths.csv")

cat("\nFiles created:\n")
cat("- selected_isolates.csv: List of selected isolates\n")
cat("- isolate_file_paths.csv: Detailed file paths\n")

cat("\nNext step: Use these files to run the pipeline\n")