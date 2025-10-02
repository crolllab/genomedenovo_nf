# Multi-Sample Scaffold Input Feature

## Summary

Modified the pipeline to accept a directory containing multiple sample folders, each with a `scaffolds.fasta` file. The pipeline now processes each sample independently through BLAST analysis and combines all results into a single output file.

## Changes Made

### 1. Modified `workflows/denovo_query.nf`

**Before**: Only supported direct file patterns (e.g., `*.fasta`)

**After**: Now detects if `--scaffolds` parameter is a directory and searches for `scaffolds.fasta` in subdirectories:

```groovy
if (file(scaffolds_path).isDirectory()) {
    // Pattern: parent_dir/*/scaffolds.fasta
    scaffolds_pattern = "${scaffolds_path}/*/scaffolds.{fasta,fa,fna}"
    
    local_scaffolds_ch = Channel
        .fromPath(scaffolds_pattern, checkIfExists: true)
        .map { file ->
            // Use parent directory name as sample_id
            def sample_id = file.parent.name
            [sample_id, file]
        }
}
```

**Sample ID extraction**: Uses the subdirectory name (e.g., `05STB001_1`) as the sample identifier.

### 2. Updated Documentation

#### README.md
- Added explanation of two input modes for `--scaffolds`
- Provided example directory structure
- Added example run command for multi-sample BLAST analysis

#### .github/copilot-instructions.md
- Documented scaffold input modes
- Updated data flow description
- Added sample ID extraction patterns

### 3. Created Test Script

**`test_multisample_scaffolds.sh`**: Comprehensive test script that:
- Counts samples in input directory
- Runs pipeline on all samples
- Validates output files
- Displays summary statistics
- Shows samples with BLAST hits

## Usage

### Input Directory Structure
```
assemblies/
├── sample1/
│   └── scaffolds.fasta
├── sample2/
│   └── scaffolds.fasta
└── sample3/
    └── scaffolds.fasta
```

### Command
```bash
nextflow run main.nf -config nextflow.config \
    -profile local_highCPU \
    -work-dir '/scratch/daniel_tmp' \
    --outdir './nf_output' \
    --scaffolds '/path/to/assemblies' \
    --blast_query './wagA_gene.fa'
```

### Output Files

**Per-sample outputs** (in `nf_output/blast_hit_sequences/`):
- `{sample_id}_top_hits.fasta` - Extracted hit sequences for each sample
- `{sample_id}_top_hits_info.tsv` - Metadata for each hit

**Combined outputs** (in `nf_output/`):
- `all_blast_hits.fasta` - All hit sequences from all samples concatenated
- `all_blast_hits_summary.tsv` - Combined metadata table with Sample_ID column

**QC outputs**:
- `assembly_summary_metrics.tsv` - QUAST metrics for all samples
- `quast_reports/{sample_id}_quast/` - Individual QUAST reports

## Example Workflow

1. **Prepare assemblies directory**:
   ```bash
   ls /data/daniel/wagABC_blast/assemblies/
   # Output: 05STB001_1/ 05STB003_1/ 05STB011_1/ ...
   ```

2. **Run pipeline**:
   ```bash
   ./test_multisample_scaffolds.sh
   ```

3. **Check combined results**:
   ```bash
   # See all samples with hits
   grep "^>" nf_output/all_blast_hits.fasta
   
   # View summary table
   column -t -s $'\t' nf_output/all_blast_hits_summary.tsv | less
   ```

## Benefits

1. **Scalable**: Process hundreds of samples in parallel
2. **Organized**: Sample IDs automatically extracted from directory structure
3. **Comprehensive**: Single combined FASTA with all hits across all samples
4. **Traceable**: Each sequence header includes sample ID and coordinates
5. **Flexible**: Still supports direct file patterns for backwards compatibility

## Testing

Run the test script:
```bash
./test_multisample_scaffolds.sh
```

This will:
- Process all samples in `/data/daniel/wagABC_blast/assemblies`
- Extract BLAST hits for each sample
- Generate combined output files
- Display summary statistics

## Backwards Compatibility

✅ Old usage still works:
```bash
--scaffolds '/path/to/specific_file.fasta'
--scaffolds '/path/to/files/*.fasta'
```

✅ New usage for multi-sample directories:
```bash
--scaffolds '/path/to/assemblies'  # Automatically finds */scaffolds.fasta
```
