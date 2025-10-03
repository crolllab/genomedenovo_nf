# Feature: Include Query Sequence in Combined BLAST Hits Output

## Overview
Modified the `CombineBlastHits` process to prepend the original query sequence to each combined FASTA file containing BLAST hits. This makes the output files self-contained and easier to use for downstream analysis like multiple sequence alignment.

## Changes Made

### 1. Updated `CombineBlastHits` Process (`modules/extract_blast_hits.nf`)

**Input Change:**
```groovy
// Before:
tuple val(query_id), path(hit_sequences), path(hit_info_files)

// After:
tuple val(query_id), path(query_fasta), path(hit_sequences), path(hit_info_files)
```

**Script Change:**
```bash
# First, add the query sequence to the combined FASTA
if [ -f "${query_fasta}" ]; then
    cat "${query_fasta}" > ${query_id}_all_blast_hits.fasta
    echo "" >> ${query_id}_all_blast_hits.fasta
else
    touch ${query_id}_all_blast_hits.fasta
fi

# Then combine all BLAST hit FASTA files (skip empty files)
for fasta in *_top_hits.fasta; do
    if [ "$fasta" ] && [ -s "$fasta" ]; then
        cat "$fasta" >> ${query_id}_all_blast_hits.fasta
    fi
done
```

### 2. Updated Workflow (`workflows/denovo_query.nf`)

**Channel Preparation:**
```groovy
// Join with query sequences to include in output
queries_for_combine = individual_queries
    .map { query_id, query_file -> [query_id, query_file] }

// Combine sequences and hit info for each query
combined_by_query = queries_for_combine
    .join(sequences_by_query)
    .join(hit_info_by_query)

// Combine all extracted sequences for each query (with query sequence prepended)
CombineBlastHits(combined_by_query)
```

## Output Format

Each `<query_id>_all_blast_hits.fasta` file now contains:

1. **Query sequence** (first entry)
   - Original FASTA header from query file
   - Complete query sequence

2. **BLAST hits** (subsequent entries)
   - All extracted sequences from assemblies
   - With padding if `--padding > 0`
   - Sorted by sample

### Example Output Structure

```
>wagA
ATGCGATCGATCG...

>sample001_wagA_hit_scaffold123:1500-1600
ATGCGATCGATCG...

>sample002_wagA_hit_NODE_45:2300-2450
GCTAGCTAGCTAG...

>sample003_wagA_hit_contig_78:5100-5200
TAGCTAGCTAGCT...
```

## Benefits

1. **Self-contained files**: Each output contains both query and hits
2. **Ready for alignment**: Can directly pass to MAFFT or other aligners
3. **Better documentation**: Query sequence preserved with results
4. **Easier manual inspection**: Can see what was searched for
5. **Reproducibility**: Query sequence archived with results

## Usage

No command-line changes needed. This feature is automatically applied to all BLAST search workflows:

```bash
nextflow run main.nf \
    --scaffolds /path/to/assemblies \
    --blast_query query.fa \
    --padding 100
```

## Output Files

For each query in `query.fa`, you'll get:
- `<query_id>_all_blast_hits.fasta` - Query + all hits (with query prepended)
- `<query_id>_all_blast_hits_summary.tsv` - TSV with hit coordinates and metrics

## Notes

- Query sequence is added **first** in the FASTA file
- Empty line separator added after query sequence
- No changes to TSV summary files (still contain only hit information)
- Works with both nucleotide and protein queries
- Compatible with `--padding` parameter

## Testing

Validate syntax:
```bash
nextflow run main.nf -config nextflow.config --help
```

Run with example data:
```bash
./run.sh  # Uses existing configuration
```

Check output structure:
```bash
head -20 nf_output/wagA/wagA_all_blast_hits.fasta
# Should show query sequence first, then BLAST hits
```
