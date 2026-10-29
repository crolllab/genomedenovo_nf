# HPC usage & utilities

## Notes on HPC execution

### Storage requirements

Each task runs in a node-local scratch folder (`process.scratch = true`), and only its outputs are copied back to the work directory. SPAdes intermediates are removed after assembly and BUSCO keeps only its summaries and single-copy sequences unless `--keep_full_busco` is set. Large panels can still require substantial scratch space. Always point `-work-dir` to a fast, high-capacity scratch filesystem:

```bash
nextflow run main.nf ... -work-dir '/path/to/scratch/genomedenovo_work'
```

The work directory is **not** deleted at the end of a run (`cleanup = false`), so `-resume` remains available, for example to add a query or switch on `--search_tool exonerate` without re-assembling. Delete it yourself once you have checked the outputs.

Keep the Nextflow process alive in a `tmux` or `screen` session — even when using `-profile slurm`, the Nextflow process must remain running until all jobs complete.

### SLURM time limits

Assembly (SPAdes) and SRA download tasks request a **7-day** time limit, most other tasks 1 hour to 2 days. All requests are capped by `--max_time` (default `7.d`). If your cluster enforces shorter queue limits, lower `--max_time` or edit the `time` directives in `nextflow.config`.

### Task resources

Resources are set per process in `nextflow.config`, matched to what each tool used in SLURM accounting of large runs. Memory grows on retry for the steps that can run out of memory. Every request is capped by `--max_cpus`, `--max_memory` and `--max_time`.

| Process | CPUs | Memory (first attempt) | Time |
|---------|------|------------------------|------|
| SPAdes | profile | 16 GB, × attempt | 7 d |
| fastp trimming | 8 | 2 GB, × attempt | 1 d |
| QUAST | 1 | 1 GB | 1 d |
| BUSCO | 8 | 8 GB, × attempt | 2 d |
| BLAST (query, paralogs, reference locus) | 1 | 4 GB | 1 d |
| Exonerate | 1 | 4 GB | 1 d |
| MAFFT (BLAST hits) | profile | 8 GB, doubled per attempt | 1 d |
| IQ-TREE (per-query trees) | 4 | 8 GB, × attempt | 1 d |
| RAxML-NG / IQ-TREE (BUSCO tree) | profile | 8 / 16 GB, × attempt | 2 d |
| RepeatMasker | 4 | 4 GB, × attempt | 2 d |
| bwa-mem2 mapping | profile | 8 GB, × attempt | 2 d |
| SAMtools sort | 4 | 4 GB, × attempt | 4 h |

"profile" means the CPU count of the execution profile (`local` 4, `local_highCPU` and `slurm` 16).

### Singularity / Apptainer

**Image cache** — by default, images are pulled to `$HOME/.singularity/cache`. Worker nodes must have read access to this path. If `/home` is not shared across nodes, set `cacheDir` in the `singularity {}` block of `nextflow.config` to a path on a shared filesystem:

```groovy
singularity {
    cacheDir = "/scratch/$USER/.singularity/cache"
}
```

**Bind mounts** — on clusters running Apptainer (the successor to Singularity), automatic bind-mounting can be disabled by the sysadmin via `/etc/apptainer/apptainer.conf` (`mount hostfs = no`). If tasks fail with "file not found" errors inside the container, add explicit bind paths in `nextflow.config`:

```groovy
singularity.runOptions = "--bind /scratch,/data"
```

**BUSCO lineages** — without `--busco_db_path`, the lineage is downloaded once per run into `<outdir>/busco_db_cache/`. On clusters whose compute nodes have no internet access, download it beforehand with `busco --download <lineage> --download_path <folder>` and pass `--busco_db_path <folder>`.

---

### Concurrency and rate limits

Three settings in `nextflow.config` govern how much work the pipeline puts on the scheduler:

| Setting | Default | Applies to |
|---------|---------|------------|
| `executor.queueSize` | 500 | All processes combined |
| `executor.submitRateLimit` | `'240/1min'` | All processes combined |
| `maxForks` | per process, see below | A single process |

`maxForks` additionally protects resources that are not the CPU:

| Process | Default `maxForks` | Notes |
|---------|-------------------|-------|
| SRA download | 10 | NCBI will throttle or block connections with too many parallel downloads (PE and SE separately) |
| fastp trimming | 20 | Very I/O-intensive; reduce if storage I/O is a bottleneck (PE and SE separately) |

---

## Utilities

### Download SRA files manually

If `sra-tools` encounters connection resets or other SRA-side errors, download files manually with `fastq-dump` and use `--reads` to pass the downloaded files instead:

```bash
# Install sra-tools via micromamba
micromamba install -c bioconda sra-tools

# Create a download script (one accession per fastq-dump call)
cat > SRA_download.sh << 'EOF'
#!/bin/bash
fastq-dump --split-files --gzip SRR24910574
fastq-dump --split-files --gzip SRR24910575
fastq-dump --split-files --gzip SRR25074049
EOF

# Parallelise the download (limit to 10 concurrent)
parallel -j 10 < SRA_download.sh
```

!!! tip
    Do not exceed ~10 parallel downloads; NCBI will stall the connection.

Then run the pipeline pointing to the downloaded files, using `--SRR_sample_map` to name the samples:

```bash
nextflow run main.nf ... --reads '/path/to/downloads/*_{1,2}.fastq.gz'
```

---

### Re-use assemblies from a previous run

The assemblies of a run are written as `2_assemblies/{sample_id}/scaffolds.fasta`, the layout `--scaffolds` expects. To search them for new queries without re-assembling (and without keeping the work directory), pass the folder:

```bash
nextflow run main.nf -profile slurm --slurm_queue long \
  --scaffolds 'previous_run/2_assemblies' \
  --blast_query new_queries.fa --search_tool exonerate
```

---

### Planned features

Planned work is listed in [FUTURE_FEATURES.md](https://github.com/crolllab/genomedenovo_nf/blob/main/FUTURE_FEATURES.md). Feature requests are welcome on the [GitHub Issues](https://github.com/crolllab/genomedenovo_nf/issues) page.
