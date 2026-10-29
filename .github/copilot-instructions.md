
# Copilot Instructions for genomedenovo_nf

## Project Overview
This project is a **Nextflow DSL2 pipeline** for de novo genome assembly and BLAST-based querying of fungal genomes (mainly *Zymoseptoria tritici*). It supports three input modes: local FASTQ files, NCBI SRA accessions, or pre-assembled scaffolds. The pipeline is modular, with all processes in `modules/`, a single workflow in `workflows/denovo_query.nf`, and the entry point in `main.nf`.

## Architecture & Data Flow
- **Input routing**: Three parallel paths (local FASTQ, SRA, or scaffolds) are routed and normalized to `[sample_id, file(s)]` tuples. SRA accessions are resolved and split into PE/SE using `resolve_SRA.nf` and `download_SRA.nf`.
- **Channel mixing**: All input types are merged into a single channel for downstream processing using `.mix()` (see `main.nf`).
- **Parallel processing**: Each sample is processed independently for assembly (SPAdes), QC (QUAST), and BLAST search. Results are aggregated with `.collect()` for summary outputs.
- **Conditional execution**: Input validation is performed at the workflow level (not in processes). Example:
  ```groovy
  if (!params.reads && !params.SRA_index && !params.scaffolds) {
      exit 1, "ERROR: No input provided..."
  }
  ```

## Developer Workflows
- **Nextflow availability**: Always use `module load Nextflow`.
- **Running**: Use `nextflow run main.nf -config nextflow.config ...` with appropriate `-profile` (`local`, `local_highCPU`, `slurm`). See `run.sh` for examples.
- **Testing**: Use `local_highCPU` for development/testing with small datasets. Use `-resume` to skip completed steps. Run `./validate_syntax.sh` for quick syntax checks.
- **Debugging**: Most processes use `errorStrategy 'ignore'` to allow the pipeline to continue on sample failure. Logs rotate as `.nextflow.log.*`.
- **Memory management**: Critical for large sample sets - use `maxForks` limits in `nextflow.config` to prevent memory exhaustion on concurrent processes.

## Project-Specific Conventions
- **FASTQ naming**: Paired-end reads must match patterns like `*_1.fq.gz`, `*_R1.fq.gz`, or `*_R1_001.fq.gz`. See `README.md` for glob patterns.
- **Scaffold input**: `--scaffolds` can be a directory of subfolders (each with `scaffolds.fasta`) or a glob of files. Sample IDs are extracted from subdirectory or filename.
- **Channel tuple structure**: All channels between modules use `[sample_id, file(s)]` for consistency.
- **Containerization**: All tools run in local Singularity containers from `./singularity/`. Container paths are relative and must match the URL-encoded filenames in `nextflow.config` (e.g., `%3A` for `:`).
- **Resource cleanup**: SPAdes modules delete intermediate files after assembly. BLAST modules remove temporary databases after search.
- **Bash in Nextflow**: Use proper variable escaping (`\$` for bash vars) and avoid `/` delimiters in sed commands (use `|` instead). Always include `set -euo pipefail` for error handling.

## Output Organization
All outputs are written to `nf_output/` (configurable via `--outdir`). Key subfolders:
- `spades_assemblies/`, `spades_stats/`, `quast_reports/`, `blast_results/`, `blast_hit_sequences/`
- Aggregated files: `assembly_summary_metrics.tsv`, `all_blast_hits.fasta`, `all_blast_hits_summary.tsv`, `all_sequences_aligned.fasta`, `mafft_alignment_log.txt`

## Integration & External Dependencies
- **NCBI SRA**: Accessed via `entrez-direct` and `sra-tools` (in containers)
- **LEGserv**: Lab file server for input data (e.g., `/legserv/NGS_data/`)
- **Tool versions**: See `singularity/` for exact versions (SPAdes 4.2.0, QUAST 5.3.0, BLAST+ 2.17.0, MAFFT 7.525, etc.)

## When Modifying This Codebase
1. Add new processes in `modules/`, include in `workflows/denovo_query.nf`, and map containers in `nextflow.config`.
2. Maintain `[sample_id, file(s)]` tuple structure for all channels.
3. Add new parameters to `main.nf` and document them in `docs/configuration.md` (the GitHub Pages site).
4. Update both the `singularity/` directory and `nextflow.config` for new/updated containers.
5. Use `publishDir` with `mode: 'copy'` for all final outputs.

## Known Limitations
- Single-end SRA reads are downloaded but not yet fully integrated into assembly.
- No variant calling or reference-based mapping.
- Some debug output remains in `extract_blast_hits.nf`.

## Common Pitfalls & Solutions
- **Memory exhaustion**: With 100+ samples, set `maxForks` limits in `nextflow.config` for BLAST/QUAST processes to prevent system overload.
- **Empty BLAST outputs**: Ensure proper bash variable escaping (`\$var`) and use `|` instead of `/` as sed delimiters in process scripts.
- **Container path issues**: Singularity container names must be URL-encoded in `nextflow.config` (`:` becomes `%3A`).
- **Channel mixing failures**: All input channels must emit `[sample_id, file(s)]` tuples before `.mix()` operations.
