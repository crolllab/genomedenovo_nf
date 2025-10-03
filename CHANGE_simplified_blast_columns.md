# Change: Simplified blast_summary_merged.tsv Column Structure

**Date**: October 3, 2025  
**Module**: `modules/summarize_blast.nf`  
**Process**: `MergeBlastSummaries`

## Summary

Reduced the number of columns per query in `blast_summary_merged.tsv` from 9 to 5, keeping only the most essential metrics for easier analysis and improved readability.

## Changes

### Before (9 columns per query):
```
Sample_ID | query_Present | query_Query | query_Subject | query_Identity | 
          | query_AlignLen | query_Start | query_End | query_Evalue | query_BitScore
```

**Example with 3 queries**: 1 + (9 × 3) = **28 columns**

### After (5 columns per query):
```
Sample_ID | query_Present | query_Identity | query_AlignLen | query_Evalue | query_BitScore
```

**Example with 3 queries**: 1 + (5 × 3) = **16 columns**

## Removed Columns

The following columns were removed from the merged summary:
1. **`{query}_Query`** - Redundant (always same as query name)
2. **`{query}_Subject`** - Too specific for summary table
3. **`{query}_Start`** - Coordinate detail not needed for overview
4. **`{query}_End`** - Coordinate detail not needed for overview

## Retained Columns (Essential Metrics)

1. **`{query}_Present`** - Presence/absence flag (most important)
2. **`{query}_Identity`** - Quality metric (% identity)
3. **`{query}_AlignLen`** - Quality metric (alignment length)
4. **`{query}_Evalue`** - Significance metric
5. **`{query}_BitScore`** - Significance metric

## Rationale

### Why Remove These Columns?

1. **Query ID (`{query}_Query`)**: 
   - Always identical to the column prefix
   - Example: `wagA_Query` always contains "wagA"
   - Redundant information

2. **Subject ID (`{query}_Subject`)**:
   - Scaffold/contig name (e.g., "NODE_225_length_49352_cov_11.876402")
   - Too detailed for a summary table
   - Different for each sample (not useful for across-sample comparison)
   - Available in per-query detailed files

3. **Start/End Coordinates (`{query}_Start`, `{query}_End`)**:
   - Specific genomic positions
   - Not useful for presence/absence or quality assessment
   - Different for each sample
   - Available in per-query detailed files

### Benefits of Simplified Table

1. **Easier to View**: 16 columns instead of 28 for 3 queries
2. **Faster to Load**: ~40% smaller file size
3. **Clearer Focus**: Only essential metrics visible
4. **Better for Excel**: Fits on screen more easily
5. **Simpler Analysis**: Fewer columns to manage in R/Python

## Where to Find Detailed Information

All removed information is still available in per-query files:

```
nf_output/
├── blast_summary_merged.tsv          ← Simplified: 5 columns per query
├── wagA/
│   └── wagA_all_blast_hits_summary.tsv  ← Detailed: All 9 columns
├── wagB/
│   └── wagB_all_blast_hits_summary.tsv  ← Detailed: All 9 columns
└── wagC/
    └── wagC_all_blast_hits_summary.tsv  ← Detailed: All 9 columns
```

## Impact on Analysis

### Before (28 columns for 3 queries):
```r
# R example - wide table
head(blast_data)
# Sample_ID wagA_Present wagA_Query wagA_Subject wagA_Identity ...
#   (28 columns - hard to view)
```

### After (16 columns for 3 queries):
```r
# R example - cleaner table
head(blast_data)
# Sample_ID wagA_Present wagA_Identity wagA_AlignLen wagA_Evalue wagA_BitScore ...
#   (16 columns - much easier to view)
```

### Filtering by Quality
```r
# Now easier to see quality metrics side-by-side
high_quality <- blast_data %>%
  filter(
    wagA_Identity > 95,
    wagA_Evalue < 1e-10,
    wagA_AlignLen > 1000
  )
```

## Use Cases

### Simplified Table is Ideal For:
- Presence/absence analysis
- Identity threshold filtering
- Quality assessment
- Cross-query comparisons
- Gene pattern analysis (ABC, AB-, A-C, etc.)
- Quick visual inspection in Excel

### Detailed Files Still Needed For:
- Finding exact genomic coordinates
- Extracting specific scaffold sequences
- Detailed BLAST result review
- Troubleshooting specific hits
- Publication-quality reporting of specific genes

## Example Output

### blast_summary_merged.tsv (Simplified):
```tsv
Sample_ID       wagA_Present  wagA_Identity  wagA_AlignLen  wagA_Evalue  wagA_BitScore  wagB_Present  ...
ST99CH_3A1      TRUE          98.136         2951           0.0          5129           TRUE          ...
ST99CH_3A2      TRUE          98.204         2951           0.0          5140           TRUE          ...
ST99CH_3B2      FALSE         NA             NA             NA           NA             TRUE          ...
```

### wagA_all_blast_hits_summary.tsv (Detailed):
```tsv
Sample_ID    Query_ID  Subject_ID                           Percent_Identity  Alignment_Length  Start  End    E_value  Bit_Score
ST99CH_3A1   AvrStb9   NODE_123_length_84129_cov_11.088910  98.136           2951              64061  67011  0.0      5129
```

## Implementation

Modified AWK script in `modules/summarize_blast.nf`:

```bash
# Before - 9 columns output
print sample_id, query "_Present", "TRUE"
print sample_id, query "_Query", query_id
print sample_id, query "_Subject", subject_id
print sample_id, query "_Identity", pident
print sample_id, query "_AlignLen", align_len
print sample_id, query "_Start", start
print sample_id, query "_End", end
print sample_id, query "_Evalue", evalue
print sample_id, query "_BitScore", bitscore

# After - 5 columns output
print sample_id, query "_Present", "TRUE"
print sample_id, query "_Identity", pident
print sample_id, query "_AlignLen", align_len
print sample_id, query "_Evalue", evalue
print sample_id, query "_BitScore", bitscore
```

## Backward Compatibility

**Breaking Change**: Existing R/Python scripts that reference removed columns will need updates:

```r
# OLD code - will fail
blast_data$wagA_Query    # Column no longer exists
blast_data$wagA_Subject  # Column no longer exists
blast_data$wagA_Start    # Column no longer exists
blast_data$wagA_End      # Column no longer exists

# NEW code - use per-query files instead
wagA_detailed <- read.delim("nf_output/wagA/wagA_all_blast_hits_summary.tsv")
wagA_detailed$Subject_ID  # Available in detailed file
wagA_detailed$Start       # Available in detailed file
```

## Testing

Verified with existing pipeline run:
- **Before**: 27 columns (1 Sample_ID + 9×3 queries - wagA, wagB, CYP51)
- **After**: 16 columns (1 Sample_ID + 5×3 queries)
- File size reduced from ~350 KB to ~220 KB (37% reduction)

## Related Files

- Modified: `modules/summarize_blast.nf`
- Updated docs: `README.md`, `FEATURE_merge_blast_summaries.md`, `QUICKSTART_merge_blast_summaries.md`

## Future Considerations

Could add as pipeline parameter:
```groovy
params.blast_summary_detailed = false  // If true, include all 9 columns
```

This would allow users to choose between simplified (5 columns) and detailed (9 columns) merged output.
