# Copilot Instructions for genomedenovo_nf

## Project Overview
This is a **Nextflow DSL2 pipeline** for de novo genome assembly and BLAST-based querying of fungal genomes (primarily *Zymoseptoria tritici*). The pipeline handles three input modes: local FASTQ files, NCBI SRA accessions, or pre-assembled scaffolds.

## Architecture & Data Flow

### Pipeline Structure
- **Entry point**: `main.nf` → includes `workflows/denovo_query.nf`
- **Modular design**: All processes defined in `modules/*.nf` files
- **Single workflow**: `denovo_query` orchestrates all steps via channel mixing and conditional execution

### Data Flow Pattern
1. **Input routing** (3 parallel paths):
   - Local FASTQ pairs → `read_pairs_local_ch`
   - SRA accessions → `SRAresolve` → split PE/SE → `SRAdownloadPE/SE`
   - Pre-assembled scaffolds → `local_scaffolds_ch`

2. **Channel mixing**: All paths converge into `all_scaffolds_ch` using `.mix()` operations
3. **Parallel processing**: QUAST (QC) and BLASTn run independently on all scaffolds
4. **Aggregation**: Results collected with `.collect()` for summary generation

### Critical Architectural Patterns

**Channel mixing for multi-source input**:
```groovy
combined_pe_ch = sra_pe_formatted.mix(local_pe_formatted)
all_scaffolds_ch = spades_scaffolds.mix(local_scaffolds_ch)
```

**Tuple structure standardization**: All channels emit `[sample_id, file(s)]` tuples for consistency across modules.

**Conditional workflow execution**: Input validation at workflow start, not in processes:
```groovy
if (!params.reads && !params.SRA_index && !params.scaffolds) {
    exit 1, "ERROR: No input provided..."
}
```

## Singularity Container Management

**Critical**: All bioinformatics tools run in **local Singularity containers** stored in `./singularity/` directory.

- Container paths are **relative** (e.g., `'./singularity/spades%3A4.2.0--haf24da9_0'`)
- Containers are pre-built from Galaxy Project depot (see README)
- Each process in `nextflow.config` maps to a specific container using `withName:` directive
- URL-encoded characters in filenames (e.g., `%3A` for `:`) must match exactly

**Container assignment pattern** in `nextflow.config`:
```groovy
withName:SPAdes {
    container = './singularity/spades%3A4.2.0--haf24da9_0'
}
```

## Critical Developer Workflows

### Running the Pipeline

**Always use `micromamba` environment**:
```bash
micromamba activate nf_gp_env
```

**Profile selection** (`-profile` in `nextflow.config`):
- `local`: 4 CPUs per task, no SLURM
- `local_highCPU`: 24 CPUs per task (for development/testing)
- `slurm`: Submit jobs to SLURM cluster (production)

**Essential parameters**:
- `--reads`: Glob pattern for local FASTQ pairs (must use single quotes!)
- `--SRA_index`: File with SRA accessions (one per line)
- `--scaffolds`: Pre-assembled genomes (skip assembly step)
- `--blast_query`: FASTA file for BLAST search
- `--NCBI_API_key`: Required for SRA operations
- `-work-dir`: Nextflow temp directory (use `/scratch` for large runs)
- `--outdir`: Final output location (default: `./nf_output`)

**Example invocation** (see `run.sh`):
```bash
nextflow run main.nf -config nextflow.config \
    -profile local_highCPU \
    -work-dir '/scratch/daniel_tmp' \
    --outdir './nf_output' \
    --reads '/legserv/NGS_data/**_{1,2}.fq.gz' \
    --SRA_index './SRA_accessions.txt' \
    --blast_query './wagA_gene.fa' \
    --NCBI_API_key $NCBI_API_KEY
```

### Debugging & Development

**Nextflow resume**: Use `-resume` to skip completed processes (requires intact `work-dir`)

**Error handling**: Most processes use `errorStrategy 'ignore'` to continue pipeline on individual sample failures

**Logging**: Nextflow logs rotate (`.nextflow.log.1`, `.nextflow.log.2`, etc.)

**Testing workflow changes**: Use `local_highCPU` profile with small datasets before submitting to SLURM

## Project-Specific Conventions

### FASTQ File Naming
The pipeline expects **paired-end** reads with flexible naming:
- `*_1.fq.gz` / `*_2.fq.gz`
- `*_R1.fq.gz` / `*_R2.fq.gz`
- `*_R1_001.fq.gz` / `*_R2_001.fq.gz`

Pattern: `--reads '/path/**_{,R}{1,2}{,_001,_001_*}.{fq,fastq}.gz'`

### SRA Resolution Pattern
`SRAresolve` process splits accessions into PE/SE based on NCBI metadata using `esearch` + `efetch`:
- Outputs: `NCBI_SRR_PE_accessions.txt` and `NCBI_SRR_SE_accessions.txt`
- Accepts any SRA accession type (SRR, SRX, SRP, PRJNA)
- Uses `maxForks = 50` in config to limit concurrent downloads

### Resource Management
**SPAdes-specific cleanup**: Module deletes intermediate directories (`K*/`, `corrected/`, `tmp/`) and input FASTQ files after assembly to save disk space. This is critical for large-scale runs.

**BLAST workflow**: Creates temporary database per sample, searches, then deletes database files immediately.

### Output Organization
```
nf_output/
├── spades_assemblies/         # Contigs and scaffolds per sample
├── spades_stats/              # Assembly logs
├── quast_reports/             # Per-sample QUAST reports
├── assembly_summary_metrics.tsv  # Aggregated QC metrics
├── blast_results/             # Raw BLAST TSV outputs
├── blast_hit_sequences/       # Extracted hit regions (FASTA + info TSV)
├── all_blast_hits.fasta       # Combined hits from all samples
└── all_blast_hits_summary.tsv # Combined hit metadata
```

## Integration Points

### External Services
- **NCBI SRA**: Accessed via `entrez-direct` and `sra-tools` in containers
- **LEGserv**: Lab file server for genomic data (paths like `/legserv/NGS_data/`)

### Tool Versions (in containers)
- SPAdes: 4.2.0
- QUAST: 5.3.0
- BLAST+: 2.17.0
- entrez-direct: 24.0
- sra-tools: 3.2.1
- samtools: 1.22.1 (for FASTA indexing in BLAST hit extraction)

### Key Files Referenced Across Modules
- `params.outdir`: All publishDir paths use this variable
- `params.NCBI_API_key`: Implicitly used by entrez-direct (environment variable injection expected)
- Process tag: Always format as `"ProcessName: ${sample_id}"` for logging

## Current Development Status

**Status**: EARLY DRAFT - core assembly + BLAST pipeline functional

**Working modules**:
- SRA download and resolution (PE/SE splitting)
- SPAdes assembly with cleanup
- QUAST QC with metric summarization
- BLAST search and hit extraction

**Known limitations** (per README):
- Single-end SRA reads are downloaded but not yet fully integrated into SPAdes assembly path
- No variant calling or reference-based mapping implemented
- BLAST hit extraction debugging output still present in `extract_blast_hits.nf`

## When Editing This Codebase

1. **Adding new processes**: Create in `modules/`, include in `workflows/denovo_query.nf`, add container mapping in `nextflow.config`
2. **Channel modifications**: Maintain `[sample_id, file(s)]` tuple structure
3. **New parameters**: Add to `main.nf` params list, document in README example runs
4. **Container updates**: Update both the `singularity/` directory filename and `nextflow.config` container path
5. **Output files**: Use `publishDir` with `mode: 'copy'` (not symlinks) for final outputs
