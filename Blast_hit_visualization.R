# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)

# Set the output directory path
output_dir <- "./nf_output"

# Read the Excel file with metadata
# Skip first row (description) and use second row as headers
xlsx_file <- "41467_2023_36674_MOESM4_ESM.xlsx"
cat("Reading metadata from:", xlsx_file, "\n")
metadata <- read_excel(xlsx_file, skip = 1)

# The first row after skip becomes the header row
# Extract actual column names from the second row
metadata_headers <- as.character(metadata[1, ])
colnames(metadata) <- metadata_headers
metadata <- metadata[-1, ]  # Remove the header row from data

# Display metadata column names to verify
cat("Metadata columns:", paste(colnames(metadata)[1:5], collapse = ", "), "...\n")
cat("Number of metadata rows:", nrow(metadata), "\n")

# Get the name of the first column (should be "Isolate_name" or similar)
metadata_id_col <- colnames(metadata)[1]
cat("Metadata ID column:", metadata_id_col, "\n\n")

# Find all query-specific summary files
# Pattern matches: queryname_all_blast_hits_summary.tsv
blast_files <- list.files(
  path = output_dir, 
  pattern = "_all_blast_hits_summary\\.tsv$", 
  full.names = TRUE,
  recursive = TRUE
)

# Print found files for verification
cat("Found", length(blast_files), "summary files:\n")
print(blast_files)

# Load and combine all summary files
blast_data_list <- lapply(blast_files, function(file) {
  # Read the file
  df <- read.delim(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  # Extract query name from file path (e.g., wagA from wagA_all_blast_hits_summary.tsv)
  query_folder <- basename(dirname(file))
  
  # Add query folder as a column for tracking
  df$Query_Folder <- query_folder
  
  return(df)
})

# Combine all data frames
blast_data <- do.call(rbind, blast_data_list)

# Display summary before merge
cat("\nTotal samples across all queries:", nrow(blast_data), "\n")
cat("Queries found:", unique(blast_data$Query_Folder), "\n")
cat("\nSamples with hits:", sum(blast_data$Query_ID != "no_hits_found"), "\n")
cat("Samples without hits:", sum(blast_data$Query_ID == "no_hits_found"), "\n")

# Merge with metadata using the first column
# Rename metadata ID column to match for easier merging
metadata_renamed <- metadata
colnames(metadata_renamed)[1] <- "Sample_ID"

# Perform left join to keep all BLAST data and add metadata where available
blast_data_merged <- blast_data %>%
  left_join(metadata_renamed, by = "Sample_ID")

cat("\nAfter merge with metadata:\n")
cat("Total rows:", nrow(blast_data_merged), "\n")
cat("Samples with metadata matched:", sum(!is.na(blast_data_merged[[colnames(metadata)[2]]])), "\n")

# Use the merged data for downstream analysis
blast_data <- blast_data_merged

# Filter out samples without hits for visualization
blast_hits <- blast_data %>%
  filter(Query_ID != "no_hits_found") %>%
  mutate(
    Percent_Identity = as.numeric(Percent_Identity),
    Alignment_Length = as.numeric(Alignment_Length),
    Bit_Score = as.numeric(Bit_Score)
  )

# ===========================
# Visualization 1: Hit Rate by Query
# ===========================
hit_rate <- blast_data %>%
  group_by(Query_Folder) %>%
  summarise(
    Total_Samples = n(),
    Samples_with_Hits = sum(Query_ID != "no_hits_found"),
    Hit_Rate = Samples_with_Hits / Total_Samples * 100
  )

ggplot(hit_rate, aes(x = Query_Folder, y = Hit_Rate, fill = Query_Folder)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)", Hit_Rate, Samples_with_Hits, Total_Samples)),
            vjust = -0.5, size = 3) +
  theme_minimal() +
  labs(
    title = "BLAST Hit Rate by Query",
    x = "Query",
    y = "Hit Rate (%)",
    fill = "Query"
  ) +
  ylim(0, 100) +
  theme(legend.position = "none")

ggsave("nf_output/hit_rate_by_query.png", width = 8, height = 6)

# ===========================
# Visualization 2: Percent Identity Distribution
# ===========================
if (nrow(blast_hits) > 0) {
  ggplot(blast_hits, aes(x = Percent_Identity, fill = Query_Folder)) +
    geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
    facet_wrap(~Query_Folder, ncol = 1, scales = "free_y") +
    theme_minimal() +
    labs(
      title = "Distribution of Percent Identity by Query",
      x = "Percent Identity (%)",
      y = "Count",
      fill = "Query"
    )
  
  ggsave("nf_output/percent_identity_distribution.png", width = 10, height = 8)
}

# ===========================
# Visualization 3: Alignment Length vs Percent Identity
# ===========================
if (nrow(blast_hits) > 0) {
  ggplot(blast_hits, aes(x = Alignment_Length, y = Percent_Identity, color = Query_Folder)) +
    geom_point(alpha = 0.6, size = 3) +
    theme_minimal() +
    labs(
      title = "Alignment Length vs Percent Identity",
      x = "Alignment Length (bp)",
      y = "Percent Identity (%)",
      color = "Query"
    )
  
  ggsave("nf_output/alignment_length_vs_identity.png", width = 10, height = 6)
}

# ===========================
# Visualization 4: Hit Frequency by Continent and Query
# ===========================
# Create hit/no-hit status column
blast_data_with_status <- blast_data %>%
  filter(!is.na(Continent) & Continent != "") %>%  # Remove samples without continent data
  mutate(
    Hit_Status = ifelse(Query_ID == "no_hits_found", "No Hit", "Hit")
  )

# Count hits and no-hits by continent and query
continent_summary <- blast_data_with_status %>%
  group_by(Query_Folder, Continent, Hit_Status) %>%
  summarise(Count = n(), .groups = "drop")

cat("\n=== Hit frequency by continent ===\n")
print(continent_summary)

# Create stacked barplot
ggplot(continent_summary, aes(x = Continent, y = Count, fill = Hit_Status)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Query_Folder, ncol = 1) +
  scale_fill_manual(values = c("Hit" = "#2E7D32", "No Hit" = "#C62828")) +
  theme_minimal() +
  labs(
    title = "BLAST Hit Frequency by Continent and Query",
    x = "Continent",
    y = "Number of Samples",
    fill = "Status"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

ggsave("nf_output/hit_frequency_by_continent.png", width = 10, height = 10)

# Also create a grouped barplot showing hit rates
continent_hit_rate <- blast_data_with_status %>%
  group_by(Query_Folder, Continent) %>%
  summarise(
    Total = n(),
    Hits = sum(Hit_Status == "Hit"),
    Hit_Rate = (Hits / Total) * 100,
    .groups = "drop"
  )

ggplot(continent_hit_rate, aes(x = Continent, y = Hit_Rate, fill = Query_Folder)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  labs(
    title = "BLAST Hit Rate by Continent and Query",
    x = "Continent",
    y = "Hit Rate (%)",
    fill = "Query"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  ylim(0, 100)

ggsave("nf_output/hit_rate_by_continent_query.png", width = 12, height = 6)

# ===========================
# Visualization 5: Summary Statistics Table
# ===========================
if (nrow(blast_hits) > 0) {
  summary_stats <- blast_hits %>%
    group_by(Query_Folder) %>%
    summarise(
      N_Hits = n(),
      Mean_Identity = mean(Percent_Identity, na.rm = TRUE),
      SD_Identity = sd(Percent_Identity, na.rm = TRUE),
      Min_Identity = min(Percent_Identity, na.rm = TRUE),
      Max_Identity = max(Percent_Identity, na.rm = TRUE),
      Mean_Length = mean(Alignment_Length, na.rm = TRUE),
      Mean_BitScore = mean(Bit_Score, na.rm = TRUE)
    ) %>%
    mutate(across(where(is.numeric), ~round(., 2)))
  
  cat("\n=== Summary Statistics by Query ===\n")
  print(summary_stats)
  
  # Save to file
  write.table(summary_stats, 
              "nf_output/blast_summary_statistics.tsv", 
              sep = "\t", 
              row.names = FALSE, 
              quote = FALSE)
}

# ===========================
# Final Output: Comprehensive Sample Table
# ===========================
cat("\n=== Creating comprehensive sample table ===\n")

# Create a comprehensive table with one row per sample
# First, get all unique samples from metadata to ensure we don't miss any
all_samples <- metadata_renamed %>%
  select(Sample_ID, everything()) %>%
  distinct(Sample_ID, .keep_all = TRUE)

cat("Total samples in metadata:", nrow(all_samples), "\n")

# For each sample, summarize BLAST hits across all queries
blast_summary_per_sample <- blast_data %>%
  group_by(Sample_ID) %>%
  summarise(
    # Count hits per query
    Total_Queries_Tested = n_distinct(Query_Folder),
    Queries_With_Hits = n_distinct(Query_Folder[Query_ID != "no_hits_found"]),
    Total_Hits = sum(Query_ID != "no_hits_found"),
    
    # Individual query hit percentages
    wagA_Percent = ifelse(any(Query_Folder == "wagA" & Query_ID != "no_hits_found"),
                         round(as.numeric(Percent_Identity[Query_Folder == "wagA" & Query_ID != "no_hits_found"][1]), 2),
                         ifelse(any(Query_Folder == "wagA"), "No_Hit", "Not_Tested")),
    
    wagB_Percent = ifelse(any(Query_Folder == "wagB" & Query_ID != "no_hits_found"),
                         round(as.numeric(Percent_Identity[Query_Folder == "wagB" & Query_ID != "no_hits_found"][1]), 2),
                         ifelse(any(Query_Folder == "wagB"), "No_Hit", "Not_Tested")),
    
    wagC_Percent = ifelse(any(Query_Folder == "wagC" & Query_ID != "no_hits_found"),
                         round(as.numeric(Percent_Identity[Query_Folder == "wagC" & Query_ID != "no_hits_found"][1]), 2),
                         ifelse(any(Query_Folder == "wagC"), "No_Hit", "Not_Tested")),
    
    # Individual query hit lengths
    wagA_Length = ifelse(any(Query_Folder == "wagA" & Query_ID != "no_hits_found"),
                        as.numeric(Alignment_Length[Query_Folder == "wagA" & Query_ID != "no_hits_found"][1]),
                        ifelse(any(Query_Folder == "wagA"), "No_Hit", "Not_Tested")),
    
    wagB_Length = ifelse(any(Query_Folder == "wagB" & Query_ID != "no_hits_found"),
                        as.numeric(Alignment_Length[Query_Folder == "wagB" & Query_ID != "no_hits_found"][1]),
                        ifelse(any(Query_Folder == "wagB"), "No_Hit", "Not_Tested")),
    
    wagC_Length = ifelse(any(Query_Folder == "wagC" & Query_ID != "no_hits_found"),
                        as.numeric(Alignment_Length[Query_Folder == "wagC" & Query_ID != "no_hits_found"][1]),
                        ifelse(any(Query_Folder == "wagC"), "No_Hit", "Not_Tested")),
    
    # BLAST hit details (for samples with hits)
    Best_Hit_Identity = ifelse(any(Query_ID != "no_hits_found"), 
                              max(as.numeric(Percent_Identity[Query_ID != "no_hits_found"]), na.rm = TRUE), 
                              NA),
    Best_Hit_BitScore = ifelse(any(Query_ID != "no_hits_found"), 
                              max(as.numeric(Bit_Score[Query_ID != "no_hits_found"]), na.rm = TRUE), 
                              NA),
    Average_Hit_Identity = ifelse(any(Query_ID != "no_hits_found"), 
                                 mean(as.numeric(Percent_Identity[Query_ID != "no_hits_found"]), na.rm = TRUE), 
                                 NA),
    
    .groups = "drop"
  ) %>%
  mutate(
    # Round numeric columns
    Best_Hit_Identity = round(Best_Hit_Identity, 2),
    Best_Hit_BitScore = round(Best_Hit_BitScore, 2),
    Average_Hit_Identity = round(Average_Hit_Identity, 2)
  )

# Merge BLAST summary with all metadata
comprehensive_table <- all_samples %>%
  left_join(blast_summary_per_sample, by = "Sample_ID") %>%
  # Fill NA values for samples not tested
  mutate(
    Total_Queries_Tested = ifelse(is.na(Total_Queries_Tested), 0, Total_Queries_Tested),
    Queries_With_Hits = ifelse(is.na(Queries_With_Hits), 0, Queries_With_Hits),
    Total_Hits = ifelse(is.na(Total_Hits), 0, Total_Hits),
    wagA_Percent = ifelse(is.na(wagA_Percent), "Not_Tested", wagA_Percent),
    wagB_Percent = ifelse(is.na(wagB_Percent), "Not_Tested", wagB_Percent),
    wagC_Percent = ifelse(is.na(wagC_Percent), "Not_Tested", wagC_Percent),
    wagA_Length = ifelse(is.na(wagA_Length), "Not_Tested", wagA_Length),
    wagB_Length = ifelse(is.na(wagB_Length), "Not_Tested", wagB_Length),
    wagC_Length = ifelse(is.na(wagC_Length), "Not_Tested", wagC_Length)
  )

# Reorder columns to put BLAST summary columns first after Sample_ID
blast_cols <- c("wagA_Percent", "wagB_Percent", "wagC_Percent", 
                "wagA_Length", "wagB_Length", "wagC_Length", "Total_Hits", "Queries_With_Hits", "Total_Queries_Tested",
                "Best_Hit_Identity", "Best_Hit_BitScore", "Average_Hit_Identity")

metadata_cols <- setdiff(colnames(comprehensive_table), c("Sample_ID", blast_cols))

comprehensive_table_final <- comprehensive_table %>%
  select(Sample_ID, all_of(blast_cols), all_of(metadata_cols))

# Display summary
cat("Comprehensive table created with", nrow(comprehensive_table_final), "samples\n")
cat("Columns included:", ncol(comprehensive_table_final), "\n")
cat("BLAST summary columns:", paste(blast_cols, collapse = ", "), "\n")

# Show summary of each query results
cat("\nwagA results distribution:\n")
print(table(comprehensive_table_final$wagA_Percent, useNA = "ifany"))
cat("\nwagB results distribution:\n")
print(table(comprehensive_table_final$wagB_Percent, useNA = "ifany"))
cat("\nwagC results distribution:\n")
print(table(comprehensive_table_final$wagC_Percent, useNA = "ifany"))

# Show summary statistics for hit lengths
cat("\nHit length summaries:\n")
numeric_wagA_length <- as.numeric(comprehensive_table_final$wagA_Length[comprehensive_table_final$wagA_Length != "No_Hit" & comprehensive_table_final$wagA_Length != "Not_Tested"])
numeric_wagB_length <- as.numeric(comprehensive_table_final$wagB_Length[comprehensive_table_final$wagB_Length != "No_Hit" & comprehensive_table_final$wagB_Length != "Not_Tested"])
numeric_wagC_length <- as.numeric(comprehensive_table_final$wagC_Length[comprehensive_table_final$wagC_Length != "No_Hit" & comprehensive_table_final$wagC_Length != "Not_Tested"])

if(length(numeric_wagA_length) > 0) cat("wagA hit lengths - Mean:", round(mean(numeric_wagA_length), 1), "Range:", range(numeric_wagA_length), "\n")
if(length(numeric_wagB_length) > 0) cat("wagB hit lengths - Mean:", round(mean(numeric_wagB_length), 1), "Range:", range(numeric_wagB_length), "\n")
if(length(numeric_wagC_length) > 0) cat("wagC hit lengths - Mean:", round(mean(numeric_wagC_length), 1), "Range:", range(numeric_wagC_length), "\n")

# Save the comprehensive table
write.table(comprehensive_table_final, 
            "nf_output/comprehensive_sample_blast_summary.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

cat("\n✓ Comprehensive sample table saved to: nf_output/comprehensive_sample_blast_summary.tsv\n")

cat("\n✓ Analysis complete! Check the nf_output/ folder for plots and summary statistics.\n")
