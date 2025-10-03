#!/usr/bin/env Rscript

# ============================================================================
# Merge BLAST Hit Tables with Sample Metadata from Excel
# ============================================================================
# 
# This script merges three TSV BLAST hit summary files (wagA, wagB, wagC) 
# with sample metadata from an Excel file using Sample_ID as the key.
#
# Usage:
#   Rscript merge_blast_with_metadata.R
#
# Input files (auto-detected):
#   - nf_output/wagA/wagA_all_blast_hits_summary.tsv
#   - nf_output/wagB/wagB_all_blast_hits_summary.tsv
#   - nf_output/wagC/wagC_all_blast_hits_summary.tsv
#   - 41467_2023_36674_MOESM4_ESM.xlsx (Excel metadata file)
#
# Output:
#   - nf_output/merged_blast_results.tsv
#   - nf_output/merged_blast_results.xlsx
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
})

# Check if writexl is available
has_writexl <- requireNamespace("writexl", quietly = TRUE)

# ============================================================================
# Configuration
# ============================================================================

# Use absolute paths
project_root <- "/data/daniel/genomedenovo_nf"
base_dir <- file.path(project_root, "nf_output")
excel_file <- file.path(project_root, "41467_2023_36674_MOESM4_ESM.xlsx")

# Input files
wagA_file <- file.path(base_dir, "wagA", "wagA_all_blast_hits_summary.tsv")
wagB_file <- file.path(base_dir, "wagB", "wagB_all_blast_hits_summary.tsv")
wagC_file <- file.path(base_dir, "wagC", "wagC_all_blast_hits_summary.tsv")

# Output files
output_tsv <- file.path(project_root, "merged_blast_results.tsv")
output_xlsx <- file.path(project_root, "merged_blast_results.xlsx")

# ============================================================================
# Helper Functions
# ============================================================================

read_blast_summary <- function(filepath, query_name) {
  cat(sprintf("Reading %s...\n", basename(filepath)))
  
  if (!file.exists(filepath)) {
    warning(sprintf("File not found: %s", filepath))
    return(NULL)
  }
  
  df <- read.delim(filepath, stringsAsFactors = FALSE)
  
  # Add query name as prefix to relevant columns
  colnames(df) <- gsub("Query_ID", paste0(query_name, "_Query"), colnames(df))
  colnames(df) <- gsub("Subject_ID", paste0(query_name, "_Subject"), colnames(df))
  colnames(df) <- gsub("Percent_Identity", paste0(query_name, "_Identity"), colnames(df))
  colnames(df) <- gsub("Alignment_Length", paste0(query_name, "_AlignLen"), colnames(df))
  colnames(df) <- gsub("Start", paste0(query_name, "_Start"), colnames(df))
  colnames(df) <- gsub("End", paste0(query_name, "_End"), colnames(df))
  colnames(df) <- gsub("E_value", paste0(query_name, "_Evalue"), colnames(df))
  colnames(df) <- gsub("Bit_Score", paste0(query_name, "_BitScore"), colnames(df))
  
  cat(sprintf("  - Found %d samples with hits\n", nrow(df)))
  
  return(df)
}

# ============================================================================
# Main Processing
# ============================================================================

cat("\n========================================\n")
cat("BLAST Hit Table Merger\n")
cat("========================================\n\n")

# Read BLAST summary files
cat("Step 1: Reading BLAST summary files\n")
cat("------------------------------------\n")

wagA_data <- read_blast_summary(wagA_file, "wagA")
wagB_data <- read_blast_summary(wagB_file, "wagB")
wagC_data <- read_blast_summary(wagC_file, "wagC")

# Read Excel metadata
cat("\nStep 2: Reading Excel metadata file\n")
cat("------------------------------------\n")

if (!file.exists(excel_file)) {
  stop(sprintf("Excel file not found: %s", excel_file))
}

# Read first sheet, skipping the first 2 rows
# Row 1 is the title/description, Row 2 is all NA, Row 3 contains the actual column names
metadata <- read_excel(excel_file, sheet = 1, skip = 2)
cat(sprintf("Read metadata for %d samples\n", nrow(metadata)))
cat(sprintf("Metadata columns: %s\n", paste(head(colnames(metadata), 10), collapse = ", ")))

# Identify the sample ID column (first column or column named Isolate_name/Sample_ID)
sample_col <- colnames(metadata)[1]
cat(sprintf("Using '%s' as sample ID column\n", sample_col))

# Rename first column to Sample_ID if not already
if (sample_col != "Sample_ID") {
  colnames(metadata)[1] <- "Sample_ID"
  cat(sprintf("Renamed '%s' to 'Sample_ID'\n", sample_col))
}

# ============================================================================
# Merge Data
# ============================================================================

cat("\nStep 3: Merging data\n")
cat("------------------------------------\n")

# Start with metadata
merged_data <- metadata

# Merge wagA data
if (!is.null(wagA_data)) {
  merged_data <- left_join(merged_data, wagA_data, by = "Sample_ID")
  cat(sprintf("Merged wagA data: %d/%d samples have wagA hits\n", 
              sum(!is.na(merged_data$wagA_Identity)), nrow(merged_data)))
}

# Merge wagB data
if (!is.null(wagB_data)) {
  merged_data <- left_join(merged_data, wagB_data, by = "Sample_ID")
  cat(sprintf("Merged wagB data: %d/%d samples have wagB hits\n", 
              sum(!is.na(merged_data$wagB_Identity)), nrow(merged_data)))
}

# Merge wagC data
if (!is.null(wagC_data)) {
  merged_data <- left_join(merged_data, wagC_data, by = "Sample_ID")
  cat(sprintf("Merged wagC data: %d/%d samples have wagC hits\n", 
              sum(!is.na(merged_data$wagC_Identity)), nrow(merged_data)))
}

# ============================================================================
# Create Summary Columns
# ============================================================================

cat("\nStep 4: Creating summary columns\n")
cat("------------------------------------\n")

# Add presence/absence columns
merged_data <- merged_data %>%
  mutate(
    wagA_Present = !is.na(wagA_Identity),
    wagB_Present = !is.na(wagB_Identity),
    wagC_Present = !is.na(wagC_Identity),
    Total_Genes_Present = wagA_Present + wagB_Present + wagC_Present,
    Gene_Pattern = paste(
      ifelse(wagA_Present, "A", "-"),
      ifelse(wagB_Present, "B", "-"),
      ifelse(wagC_Present, "C", "-"),
      sep = ""
    )
  )

# Summary statistics
cat(sprintf("\nSummary Statistics:\n"))
cat(sprintf("  Total samples: %d\n", nrow(merged_data)))
cat(sprintf("  Samples with wagA: %d (%.1f%%)\n", 
            sum(merged_data$wagA_Present), 
            100 * sum(merged_data$wagA_Present) / nrow(merged_data)))
cat(sprintf("  Samples with wagB: %d (%.1f%%)\n", 
            sum(merged_data$wagB_Present), 
            100 * sum(merged_data$wagB_Present) / nrow(merged_data)))
cat(sprintf("  Samples with wagC: %d (%.1f%%)\n", 
            sum(merged_data$wagC_Present), 
            100 * sum(merged_data$wagC_Present) / nrow(merged_data)))
cat(sprintf("  Samples with all three genes: %d (%.1f%%)\n", 
            sum(merged_data$Total_Genes_Present == 3),
            100 * sum(merged_data$Total_Genes_Present == 3) / nrow(merged_data)))

# Count gene patterns
cat("\nGene Patterns:\n")
pattern_counts <- merged_data %>%
  count(Gene_Pattern) %>%
  arrange(desc(n))
print(pattern_counts)

# ============================================================================
# Save Output
# ============================================================================

cat("\nStep 5: Saving output files\n")
cat("------------------------------------\n")

# Save as TSV
write.table(merged_data, output_tsv, 
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("✓ Saved TSV: %s\n", output_tsv))

# Save as Excel (if writexl is available)
if (has_writexl) {
  writexl::write_xlsx(merged_data, output_xlsx)
  cat(sprintf("✓ Saved Excel: %s\n", output_xlsx))
} else {
  cat(sprintf("⚠ writexl package not available, skipping Excel output\n"))
  cat(sprintf("  Install with: install.packages('writexl')\n"))
}

cat("\n========================================\n")
cat("Merge Complete!\n")
cat("========================================\n\n")

cat("Output files:\n")
cat(sprintf("  - %s\n", output_tsv))
cat(sprintf("  - %s\n", output_xlsx))
cat("\nColumn structure:\n")
cat("  - Original metadata columns\n")
cat("  - wagA columns: wagA_Query, wagA_Subject, wagA_Identity, wagA_AlignLen, wagA_Start, wagA_End, wagA_Evalue, wagA_BitScore\n")
cat("  - wagB columns: (same structure as wagA)\n")
cat("  - wagC columns: (same structure as wagA)\n")
cat("  - Summary columns: wagA_Present, wagB_Present, wagC_Present, Total_Genes_Present, Gene_Pattern\n")
cat("\n")
