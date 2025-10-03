# Feature: Merge BLAST Summaries Process

**Date**: October 3, 2025  
**Module**: `modules/summarize_blast.nf`  
**Process**: `MergeBlastSummaries`

## Overview

Added a new Nextflow process to automatically merge all per-query BLAST summary files into a single wide-format TSV table. This replaces the need for manual R script execution and integrates summary table generation directly into the pipeline.

## Implementation

### New Module: `modules/summarize_blast.nf`

Created new process `MergeBlastSummaries` that:
1. Takes all `*_all_blast_hits_summary.tsv` files as input
2. Processes each query's results independently
3. Creates query-prefixed columns (e.g., `wagA_Present`, `wagA_Identity`, `wagB_Present`, etc.)
4. Merges all queries into a single wide-format table
5. Handles samples without hits (fills with NA/FALSE values)

### Column Structure

For each query (e.g., `wagA`, `wagB`, `wagC`), the following 5 essential columns are created:
- `<query>_Present` - TRUE/FALSE flag indicating if hit was found
- `<query>_Identity` - Percent identity of alignment
- `<query>_AlignLen` - Alignment length in bp
- `<query>_Evalue` - E-value
- `<query>_BitScore` - Bit score

**Note**: To reduce table width and improve readability, Query ID, Subject ID, Start, and End coordinates are not included in the merged table. These details are available in the per-query summary files (`{query_id}/{query_id}_all_blast_hits_summary.tsv`).

### Workflow Integration

Modified `workflows/denovo_query.nf`:
```groovy
// Collect all summary TSV files and merge into a single wide-format table
all_summaries = CombineBlastHits.out.summary
    .map { query_id, summary_file -> summary_file }
    .collect()

MergeBlastSummaries(all_summaries)
```

### Output

- **File**: `nf_output/blast_summary_merged.tsv`
- **Format**: Wide-format TSV with Sample_ID in first column
- **Location**: Published to main output directory

## Key Features

1. **Automatic Integration**: Runs as part of the main pipeline workflow
2. **Multi-Query Support**: Handles any number of query sequences
3. **Missing Data Handling**: Samples without hits receive NA/FALSE values
4. **Query-Specific Columns**: Each query gets its own set of prefixed columns
5. **Wide Format**: Easy to analyze in spreadsheet software or R/Python

## Technical Details

### Bash Script Logic

The process uses a multi-step approach:

1. **Process Each Query File**:
   - Extract query ID from filename
   - Parse TSV and create query-prefixed columns
   - Detect "no_hits_found" entries and convert to FALSE/NA

2. **Collect Metadata**:
   - Get all unique Sample IDs across all queries
   - Get all unique column names (sorted by query)

3. **Build Merged Table**:
   - Create header with Sample_ID + all column names
   - For each sample, collect values for all columns
   - Use grep to find values from processed files

4. **Add Summary Statistics** (optional):
   - Count total samples
   - Count samples with hits for each query

### Error Handling

- Validates that summary files exist
- Handles empty BLAST results
- Uses `set -euo pipefail` for bash error detection
- Fills missing values with "NA"

## Comparison with R Script

### Similarities
- Same column naming convention (query_Present, query_Identity, etc.)
- Same wide-format output structure
- Handles missing data appropriately

### Differences
- **Integration**: Built into pipeline vs. standalone R script
- **Input**: Works directly with pipeline outputs vs. requires manual file paths
- **Metadata**: R script can merge with Excel metadata; this process focuses on BLAST data only
- **Statistics**: R script adds Gene_Pattern column; this process focuses on raw BLAST results

## Usage

The process runs automatically when BLAST queries are provided:

```bash
nextflow run main.nf \
  -config nextflow.config \
  -profile local_highCPU \
  --scaffolds '/path/to/assemblies' \
  --blast_query './queries.fasta' \
  --outdir './nf_output'
```

Output will be in: `nf_output/blast_summary_merged.tsv`

## Future Enhancements

Potential improvements:
1. Add Gene_Pattern column similar to R script
2. Add summary statistics as separate output file
3. Support optional metadata merging from external file
4. Add filtering options (e.g., minimum identity threshold)
5. Generate visualization-ready formats

## Testing

To verify the process works correctly:

1. Run pipeline with multiple queries
2. Check that `blast_summary_merged.tsv` exists in output directory
3. Verify column count = 1 + (9 × number_of_queries)
4. Verify all samples are present
5. Verify Present/Absent flags match individual query summaries

## Related Files

- Module: `modules/summarize_blast.nf`
- Workflow: `workflows/denovo_query.nf`
- Input: Output from `CombineBlastHits` process
- Related: `merge_blast_with_metadata.R` (for adding metadata)
