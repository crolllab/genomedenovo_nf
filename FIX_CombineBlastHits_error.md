# Fixed: CombineBlastHits Output File Conflict

## Error
```
cat: all_blast_hits.fasta: input file is output file
```

Process `CombineBlastHits` failed with exit status 1 at the end of the pipeline run.

## Root Cause

The bash script used a glob pattern `*.fasta` to combine all FASTA files:

```bash
for fasta in *.fasta; do
    cat "$fasta"
done > all_blast_hits.fasta
```

**Problem**: When using `-resume` or if files are staged from cache, the glob `*.fasta` can match:
- `sample1_top_hits.fasta` ✓ (input)
- `sample2_top_hits.fasta` ✓ (input)
- `all_blast_hits.fasta` ✗ (output file!)

When bash tries to `cat all_blast_hits.fasta > all_blast_hits.fasta`, it fails because you can't read from and write to the same file simultaneously.

## Solution

Changed the glob pattern to be more specific:

```bash
# BEFORE (too broad):
for fasta in *.fasta; do

# AFTER (specific pattern):
for fasta in *_top_hits.fasta; do
```

This pattern only matches the individual sample output files (`{sample}_top_hits.fasta`) and excludes the combined output file (`all_blast_hits.fasta`).

## Additional Improvements

Added file existence check:
```bash
if [ -f "$fasta" ] && [ -s "$fasta" ]; then
```

This ensures:
- File exists (`-f`)
- File is not empty (`-s`)

## Re-running the Pipeline

Since most of the work completed successfully (3773 tasks succeeded), just re-run with `-resume`:

```bash
cd /data/daniel/genomedenovo_nf
micromamba activate nf_gp_env

SCAFFS_DIR='/data/daniel/wagABC_blast/assemblies'
QUERY='./wagA_gene.fa'
outdir='./nf_output'

nextflow run main.nf -config nextflow.config \
    -profile local_highCPU \
    -work-dir '/scratch/daniel_tmp' \
    --outdir $outdir \
    --scaffolds $SCAFFS_DIR \
    --blast_query $QUERY \
    -resume
```

The `-resume` flag will:
- ✓ Skip the 3773 already completed tasks
- ✓ Skip the 2 ignored tasks  
- ✓ Retry the 1 failed CombineBlastHits task
- ✓ Complete the 104 pending tasks that were killed

## Status from Last Run

From the log output:
- ✅ **3773 tasks succeeded** (BLAST, QUAST, ExtractBlastHits for most samples)
- ⚠️ **3 tasks failed** (likely memory errors, will retry with maxForks limits)
- ℹ️ **2 tasks ignored** (errorStrategy 'ignore')
- 🔄 **104 tasks pending** (were queued when pipeline aborted)
- ❌ **15 tasks aborted** (killed when pipeline terminated)
- ❌ **1 task failed** (CombineBlastHits - now fixed!)

## Expected Result

After re-running with `-resume`, you should get:
- `nf_output/all_blast_hits.fasta` - Combined sequences from all samples
- `nf_output/all_blast_hits_summary.tsv` - Combined metadata table
- `nf_output/assembly_summary_metrics.tsv` - QUAST metrics for all samples

## File Modified

- `modules/extract_blast_hits.nf` - Changed glob pattern in `CombineBlastHits` process
