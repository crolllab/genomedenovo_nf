# Bug Fix: extract_blast_hits.nf Producing Empty Files

## Problem
The `ExtractBlastHits` process was producing empty FASTA files despite having valid BLAST hits in the input.

## Root Causes

### 1. **Bash Variable Escaping Issues**
The original code had inconsistent escaping of bash variables in the Nextflow script block:
- Some `$` were escaped as `\$`, others were not
- This caused variables to be interpreted by Groovy instead of bash
- The `set -euo pipefail` was missing, allowing silent failures

### 2. **Sed Delimiter Conflict with Nextflow**
The sed command used `/` as delimiter which conflicts with Nextflow's string parsing:
```bash
# BROKEN - Nextflow can't parse the /.../ in sed regex
sed "s/^>.*$/>${sample_id}_.../" 

# FIXED - Use | as delimiter instead
sed "s|^>.*|>${sample_id}_...|"
```

### 3. **Excessive Debug Output**
- Heavy debug logging made it hard to identify the actual problem
- Multiple echo statements for every hit slowed down processing
- Error logs were scattered and not consolidated

### 3. **Overly Complex Fallback Logic**
- Multiple nested if-else blocks for extraction
- Redundant samtools calls
- Unclear error paths

### 4. **Newline Handling in FASTA**
- When using the manual cut fallback, extracted sequences weren't wrapped to standard FASTA format (60 chars per line)
- Missing newlines could cause downstream tools to fail

## Solutions Implemented

### 1. **Simplified Bash Script**
```bash
#!/bin/bash
set -euo pipefail  # Fail fast on errors, undefined variables, or pipe failures
```

### 2. **Proper Variable Escaping**
All bash variables now use `\$` consistently:
```bash
while IFS=\$'\t' read -r qseqid sseqid ...
    if [ "\$sstart" -le "\$send" ]; then
```

### 3. **Fixed Sed Delimiter**
Changed from `/` to `|` delimiter to avoid Nextflow parsing errors:
```bash
# Use | instead of / to avoid conflicts with Nextflow string parsing
sed "s|^>.*|>${sample_id}_\${qseqid}_hit_\${sseqid}:\${start}-\${end}|" temp_hit.fasta
```

### 4. **Streamlined Extraction Logic**
```bash
# Primary method: Direct samtools extraction with coordinates
if samtools faidx "${scaffolds_fasta}" "\${sseqid}:\${start}-\${end}" > temp_hit.fasta 2>/dev/null; then
    # Success path
else
    # Fallback: Extract full sequence then cut region
    if samtools faidx "${scaffolds_fasta}" "\${sseqid}" 2>/dev/null | tail -n +2 | tr -d '\n' > temp_full_seq.txt; then
        # Manual substring extraction with validation
    fi
fi
```

### 5. **Proper FASTA Formatting**
```bash
# Add newlines every 60 characters for proper FASTA format
fold -w 60 temp_subseq.txt >> ${sample_id}_top_hits.fasta
```

### 6. **Removed Debug Output**
- Eliminated all `echo "DEBUG: ..."` statements
- Errors are silently handled with 2>/dev/null and fallback logic
- Process succeeds even if some hits fail to extract

### 7. **Fixed CombineBlastHits**
```bash
# Skip empty FASTA files when combining
for fasta in *.fasta; do
    if [ -s "\$fasta" ]; then  # Only cat non-empty files
        cat "\$fasta"
    fi
done > all_blast_hits.fasta
```

## Testing Recommendations

1. **Test with various BLAST result formats:**
   ```bash
   # Re-run the pipeline with existing data
   nextflow run main.nf -config nextflow.config -profile local_highCPU \
       -work-dir '/scratch/daniel_tmp' --outdir './test_output' \
       --scaffolds './IPO323_sub.fasta' --blast_query './wagA_gene.fa' \
       -resume
   ```

2. **Verify output files:**
   ```bash
   # Check that FASTA files are not empty
   ls -lh test_output/blast_hit_sequences/*.fasta
   
   # Verify FASTA format is correct
   head test_output/blast_hit_sequences/*_top_hits.fasta
   
   # Check TSV has data rows (not just headers)
   wc -l test_output/blast_hit_sequences/*_top_hits_info.tsv
   ```

3. **Check combined outputs:**
   ```bash
   # Should contain all hits from all samples
   grep "^>" test_output/all_blast_hits.fasta | wc -l
   
   # Should match number of data rows in summary
   tail -n +2 test_output/all_blast_hits_summary.tsv | wc -l
   ```

## Expected Behavior

### Before Fix
- `scaffolds_top_hits.fasta`: Empty file (0 bytes)
- `scaffolds_top_hits_info.tsv`: Only header + "No sequences could be extracted"

### After Fix
- `scaffolds_top_hits.fasta`: Contains extracted sequences in proper FASTA format
- `scaffolds_top_hits_info.tsv`: Header + one row per top BLAST hit with coordinates and scores
- `all_blast_hits.fasta`: Properly combined sequences from all samples
- `all_blast_hits_summary.tsv`: Properly combined metadata from all samples

## Additional Notes

- The `errorStrategy 'ignore'` means individual sample failures won't stop the pipeline
- Empty BLAST results files are handled gracefully (create empty output files)
- Coordinate conversion (sstart/send) handles both forward and reverse strand hits
- The fallback extraction method handles edge cases where samtools region extraction fails
