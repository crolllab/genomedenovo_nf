# MergeBlastSummaries Process - Quick Start Guide

## Overview

The `MergeBlastSummaries` process automatically creates a single wide-format TSV file combining BLAST results from multiple query sequences. This runs automatically when you provide a multi-sequence BLAST query file.

## When to Use

Use this feature when you want to:
- Search multiple genes/sequences across many samples
- Compare presence/absence of multiple genes
- Analyze alignment statistics for multiple queries in a single table
- Export results to Excel/R for downstream analysis

## Example Usage

### Basic Multi-Query BLAST

```bash
# Activate environment
micromamba activate nf_gp_env

# Run pipeline with multi-query BLAST file
nextflow run main.nf \
  -config nextflow.config \
  -profile local_highCPU \
  --scaffolds '/data/assemblies' \
  --blast_query './wagABC_genes.fa' \
  --outdir './nf_output' \
  -resume
```

### Input Requirements

Your query file should contain multiple sequences in FASTA format:

```
>wagA
ATGCGATCGATCGATCG...
>wagB
GCTAGCTAGCTAGCTAG...
>wagC
TTAATTAATTAATTAAT...
```

## Output Structure

The pipeline will generate:

```
nf_output/
├── blast_summary_merged.tsv          ← Main merged output (NEW!)
├── wagA/
│   ├── wagA_all_blast_hits.fasta
│   ├── wagA_all_blast_hits_summary.tsv  ← Detailed per-query info
│   └── blast_results/
├── wagB/
│   ├── wagB_all_blast_hits.fasta
│   ├── wagB_all_blast_hits_summary.tsv  ← Detailed per-query info
│   └── blast_results/
└── wagC/
    ├── wagC_all_blast_hits.fasta
    ├── wagC_all_blast_hits_summary.tsv  ← Detailed per-query info
    └── blast_results/
```

## Merged Output Format

The `blast_summary_merged.tsv` file has this structure (simplified for clarity):

| Sample_ID | wagA_Present | wagA_Identity | wagA_AlignLen | wagA_Evalue | wagA_BitScore | wagB_Present | wagB_Identity | ... |
|-----------|-------------|---------------|---------------|-------------|---------------|-------------|---------------|-----|
| sample1   | TRUE        | 98.5          | 1234          | 0.0         | 2156          | TRUE        | 99.2          | ... |
| sample2   | FALSE       | NA            | NA            | NA          | NA            | TRUE        | 97.8          | ... |

### Column Details

For each query, you get 5 essential columns:
1. `{query}_Present` - TRUE if hit found, FALSE if no hit
2. `{query}_Identity` - Percent identity (0-100)
3. `{query}_AlignLen` - Alignment length in bp
4. `{query}_Evalue` - BLAST E-value
5. `{query}_BitScore` - BLAST bit score

Note: Query ID, Subject ID, Start, and End positions are available in the per-query detailed files (`{query_id}/{query_id}_all_blast_hits_summary.tsv`) if needed.

## Merged Output Format

The `blast_summary_merged.tsv` file has this structure:

| Sample_ID | wagA_Present | wagA_Identity | wagA_AlignLen | ... | wagB_Present | wagB_Identity | ... | wagC_Present | wagC_Identity | ... |
|-----------|-------------|---------------|---------------|-----|-------------|---------------|-----|-------------|---------------|-----|
| sample1   | TRUE        | 98.5          | 1234          | ... | TRUE        | 99.2          | ... | FALSE       | NA            | ... |
| sample2   | FALSE       | NA            | NA            | ... | TRUE        | 97.8          | ... | TRUE        | 96.5          | ... |

### Column Details

For each query, you get 5 essential columns:
1. `{query}_Present` - TRUE if hit found, FALSE if no hit
2. `{query}_Identity` - Percent identity (0-100)
3. `{query}_AlignLen` - Alignment length in bp
4. `{query}_Evalue` - BLAST E-value
5. `{query}_BitScore` - BLAST bit score

Note: Query ID, Subject ID, Start, and End positions are available in the per-query detailed files (`{query_id}/{query_id}_all_blast_hits_summary.tsv`) if needed.

## Downstream Analysis

### Loading in R

```r
# Read the merged table
blast_data <- read.delim("nf_output/blast_summary_merged.tsv", 
                         stringsAsFactors = FALSE)

# Check presence/absence
blast_data$wagA_Present
blast_data$wagB_Present
blast_data$wagC_Present

# Create gene pattern (similar to R script)
blast_data$gene_pattern <- paste0(
  ifelse(blast_data$wagA_Present, "A", "-"),
  ifelse(blast_data$wagB_Present, "B", "-"),
  ifelse(blast_data$wagC_Present, "C", "-")
)

# Count patterns
table(blast_data$gene_pattern)
```

### Loading in Python

```python
import pandas as pd

# Read the merged table
blast_data = pd.read_csv("nf_output/blast_summary_merged.tsv", sep="\t")

# Check presence/absence
print(blast_data[['Sample_ID', 'wagA_Present', 'wagB_Present', 'wagC_Present']])

# Count samples with all three genes
all_present = (
    (blast_data['wagA_Present'] == True) &
    (blast_data['wagB_Present'] == True) &
    (blast_data['wagC_Present'] == True)
)
print(f"Samples with all genes: {all_present.sum()}")
```

## Merging with Metadata

To add sample metadata (e.g., from Excel), use the R script:

```r
source("merge_blast_with_metadata.R")
```

This will create `merged_blast_results.tsv` with both BLAST data and metadata.

## Performance Notes

### Large Datasets
- Process is memory-efficient (uses bash grep/awk)
- Handles 1000+ samples easily
- Execution time: ~1-2 minutes for 1300 samples with 3 queries

### Resource Requirements
- CPU: 4 cores (default)
- Memory: 2GB (default)
- Disk: Minimal (creates one TSV file)

## Troubleshooting

### Issue: Missing samples in output

**Cause**: Sample had no hits for any query  
**Solution**: Check individual query summary files in `{query_id}/{query_id}_all_blast_hits_summary.tsv`

### Issue: All values are NA for a query

**Cause**: No samples had hits for that particular query  
**Solution**: 
1. Check BLAST parameters (e-value threshold, etc.)
2. Verify query sequence quality
3. Check if query sequence is present in your target genomes

### Issue: Process fails with "No BLAST summary files found"

**Cause**: CombineBlastHits didn't produce output  
**Solution**: 
1. Check that `--blast_query` parameter is set
2. Verify BLAST search completed successfully
3. Check `.nextflow.log` for BLAST errors

## Related Documentation

- **Feature description**: `FEATURE_merge_blast_summaries.md`
- **R script for metadata merging**: `merge_blast_with_metadata.R`
- **Pipeline README**: `README.md`
- **Module source**: `modules/summarize_blast.nf`
- **Workflow integration**: `workflows/denovo_query.nf`

## Comparison: Pipeline Process vs. R Script

| Feature | MergeBlastSummaries (Pipeline) | R Script |
|---------|-------------------------------|----------|
| Integration | Automatic in pipeline | Manual execution |
| Input | Pipeline BLAST outputs | Requires file paths |
| Metadata | BLAST data only | Merges with Excel |
| Gene Pattern | No (raw data only) | Yes |
| Summary Stats | Basic | Detailed |
| Use Case | Pipeline automation | Final analysis |

**Recommendation**: Use the pipeline process for automated analysis. Use the R script when you need to add external metadata or create publication-ready tables.
