# Memory Exhaustion Issue with Multi-Sample Processing

## Problem

When running the pipeline on a large number of samples (`--scaffolds` pointing to directory with 1300+ samples), BLASTn processes failed with memory allocation errors:

```
blastn: error while loading shared libraries: cannot allocate symbol search list: Cannot allocate memory
```

## Root Cause

**Excessive parallelization**: With the `local_highCPU` profile, each task gets 24 CPUs. When processing hundreds of samples:
- Multiple BLASTn processes run simultaneously (one per sample)
- Each BLAST process allocates memory for:
  - Database indexing
  - Query processing  
  - Thread management (24 threads × memory per thread)
- System ran out of available memory

**Example**: With 50+ concurrent BLAST processes, each using ~4GB, total memory demand exceeded system capacity.

## Solution

Added `maxForks` limits to resource-intensive processes in `nextflow.config`:

```groovy
withName:BLASTn {
    container = './singularity/blast%3A2.17.0--h66d330f_0'
    memory = '4 GB'
    maxForks = 10  // Limit concurrent BLAST processes
}

withName:QUAST {
    container = './singularity/quast:5.3.0--py39pl5321heaaa4ec_0'
    memory = '4 GB'
    maxForks = 15  // Limit concurrent QUAST processes
}

withName:ExtractBlastHits {
    container = './singularity/samtools%3A1.22.1--h96c455f_0'
    memory = '2 GB'
}
```

## Explanation of maxForks

`maxForks = 10` means:
- **Maximum 10 BLASTn processes** run simultaneously
- Additional samples queue until a slot opens
- Prevents memory exhaustion
- Still allows good parallelization

**Memory calculation**:
- 10 BLAST processes × 4 GB = 40 GB maximum memory usage for BLAST
- Plus system overhead and other processes
- Safe for systems with 64-256 GB RAM

## Alternative Solutions

### 1. Use SLURM Profile (Recommended for Large Datasets)

```bash
nextflow run main.nf -config nextflow.config \
    -profile slurm \
    -work-dir '/scratch/daniel_tmp' \
    --outdir './nf_output' \
    --scaffolds '/data/daniel/wagABC_blast/assemblies' \
    --blast_query './wagA_gene.fa'
```

**Benefits**:
- Distributes jobs across multiple compute nodes
- Each job gets dedicated resources
- Better for 100+ samples

### 2. Reduce CPUs per Task

Modify `nextflow.config` to use fewer CPUs per BLAST process:

```groovy
withName:BLASTn {
    cpus = 4  // Instead of 24
    maxForks = 20  // Can run more in parallel with less CPUs each
}
```

### 3. Process in Batches

Split samples into batches and process separately:

```bash
# Batch 1: First 100 samples
--scaffolds '/path/to/assemblies/batch1'

# Batch 2: Next 100 samples  
--scaffolds '/path/to/assemblies/batch2'
```

## Monitoring

To check how many processes are running:

```bash
# Count current BLAST processes
ps aux | grep blastn | wc -l

# Monitor memory usage
watch -n 5 'free -h'

# Check Nextflow task status
nextflow log <run_name> -f "name,status,cpus,memory,duration"
```

## Recommendations

### For Small Datasets (< 50 samples)
```bash
-profile local_highCPU  # Fast, no queue delays
```

### For Medium Datasets (50-200 samples)
```bash
-profile local_highCPU
# With maxForks limits (already added)
```

### For Large Datasets (200+ samples)
```bash
-profile slurm  # Best option
```

## Testing

After applying the fix, test with your full dataset:

```bash
cd /data/daniel/genomedenovo_nf

micromamba activate nf_gp_env

nextflow run main.nf -config nextflow.config \
    -profile local_highCPU \
    -work-dir '/scratch/daniel_tmp' \
    --outdir './nf_output_all' \
    --scaffolds '/data/daniel/wagABC_blast/assemblies' \
    --blast_query './wagA_gene.fa' \
    -resume
```

The `-resume` flag will skip already completed samples and retry the failed ones.

## Expected Behavior

With `maxForks = 10`:
- Pipeline will queue samples
- Maximum 10 BLASTn processes run at once
- Progress may be slower but won't crash
- Estimated time: ~1-2 minutes per sample for BLAST + extraction

For 1300 samples with maxForks=10:
- Theoretical time: 1300 samples ÷ 10 parallel = 130 batches × 2 min ≈ 4-5 hours total

## Verification

After completion, verify all samples were processed:

```bash
# Count output files
ls nf_output_all/blast_hit_sequences/*_top_hits.fasta | wc -l

# Check for samples with hits (non-empty FASTAs)
find nf_output_all/blast_hit_sequences -name "*_top_hits.fasta" -size +0 | wc -l

# View combined results
grep "^>" nf_output_all/all_blast_hits.fasta | wc -l
```
