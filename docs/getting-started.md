# Getting started

## Requirements

- **Nextflow** ≥ 25.0
- **Singularity** (or Apptainer) — all pipeline software runs in containers
- A POSIX-compatible file system accessible from all compute nodes (for SLURM execution)
- An [NCBI API key](https://account.ncbi.nlm.nih.gov/) (recommended when using `--SRA_index`)

---

## Step 1: Download the latest release

```bash
curl -sL https://github.com/crolllab/genomedenovo_nf/archive/refs/tags/latest.tar.gz | tar xz
mv genomedenovo_nf-latest genomedenovo_nf
cd genomedenovo_nf
```

The `latest` tag is automatically updated with every new release.

!!! tip "Alternative"
    You can also browse all releases at [github.com/crolllab/genomedenovo_nf/releases](https://github.com/crolllab/genomedenovo_nf/releases) and download a specific version, or clone the repository with `git clone https://github.com/crolllab/genomedenovo_nf.git`.

---

## Step 2: Set up Nextflow

=== "Module (recommended on HPC)"

    ```bash
    module load Nextflow
    ```

    !!! note "Croll lab"
        On LEGcompute you can skip Steps 1–2 entirely by loading the module:
        ```bash
        module load genomedenovo_nf
        ```

=== "micromamba / conda"

    ```bash
    micromamba create -n nf_gd_env
    micromamba activate nf_gd_env
    micromamba install -c conda-forge -c bioconda nextflow
    ```

---

## Step 3: Singularity images

All container images are pulled automatically by Nextflow on the first run. No manual download is required. Images are fetched from [quay.io/biocontainers](https://quay.io/organization/biocontainers) (and `rocker/tidyverse` from Docker Hub) and cached locally so subsequent runs reuse them without re-downloading.

By default the cache is stored in `$HOME/.singularity/cache`. On HPC clusters where `/home` is not shared across worker nodes, point the cache to a shared filesystem in `nextflow.config`:

Example `nextflow.config` snippet if only `/scratch` is shared:

```groovy
singularity {
    cacheDir = "/scratch/$USER/.singularity/cache"
}
```

See [HPC usage & utilities](resources.md#singularity-apptainer) for further Singularity/Apptainer configuration notes.

---

## Step 4: Try the example dataset

The repository ships with a small *E. coli* dataset from the Long-Term Evolution Experiment (LTEE, BioProject [PRJNA295606](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA295606/)) that lets you verify your setup end-to-end before working with your own data. The input files live in `example/`:

| File | Content |
|------|---------|
| `fastq/SRR*_{1,2}.fastq.gz` | Three LTEE clones (Ara−3, 5,000 to 50,000 generations), truncated to 300,000 read pairs (~19× coverage) |
| `sra_accessions.txt` | Three further LTEE runs, downloaded at runtime with `--SRA_index` |
| `lacZ_WP_000177906.fasta` | *lacZ* beta-galactosidase protein of *E. coli* REL606, the query |
| `ecoli_K12_MG1655.fasta` | *E. coli* K-12 MG1655 reference genome, for `--reference_genome` |
| `ecoli_IS_elements.fasta` | *E. coli* IS element library, for `--repeat_library` |

Run from the repository root. The smallest run assembles the three local samples and searches them for *lacZ*:

```bash
nextflow run main.nf -profile local \
    --reads 'example/fastq/SRR*_{1,2}.fastq.gz' \
    --blast_query example/lacZ_WP_000177906.fasta \
    --outdir example/output
```

Every option can be tried on the same data. This run switches all of them on:

```bash
nextflow run main.nf -profile local \
    --reads 'example/fastq/SRR*_{1,2}.fastq.gz' \
    --SRA_index example/sra_accessions.txt \
    --blast_query example/lacZ_WP_000177906.fasta \
    --search_tool exonerate \
    --busco_lineage bacteria_odb10 \
    --phylo_tree \
    --reference_genome example/ecoli_K12_MG1655.fasta \
    --locus_coverage \
    --repeat_library example/ecoli_IS_elements.fasta \
    --outdir example/output
```

Output is written to `example/output/`. The example [README](https://github.com/crolllab/genomedenovo_nf/blob/main/example/README.md) has one command per option and the sources of all files.

---

## Step 5: Prepare your inputs

At least one of `--reads`, `--SRA_index` or `--scaffolds` is required; they can be combined in one run. A query FASTA (`--blast_query`) is needed for everything beyond assembly and assembly QC.

### Local FASTQ files

The `--reads` glob pattern must be bracketed by **single quotes** and must capture paired files.
Without them the shell expands the pattern first and only the first file reaches the pipeline —
see [File paths and glob patterns](configuration.md#file-paths-and-glob-patterns). Common patterns:

```bash
# All files ending with _1.fq.gz and _2.fq.gz
--reads '/path/to/reads/*_{1,2}.fq.gz'

# Including sub-directories
--reads '/path/to/reads/**_{1,2}.fq.gz'

# Files ending with _R1.fastq.gz and _R2.fastq.gz (sample name gets a trailing _R, see below)
--reads '/path/to/reads/*_R{1,2}.fastq.gz'
```

To read from more than one location, keep everything inside one pair of quotes and separate
the patterns with a **semicolon**:

```bash
--reads '/data/runA/*_{1,2}.fq.gz;/data/runB/*_{1,2}.fq.gz'
```

The sample name is the common prefix of the two read files, see [Sample names](configuration.md#sample-names).

### SRA/ENA accessions

Create a plain-text file with one accession per line. The pipeline accepts `SRR`, `SRX`, `SRP`, `PRJNA`, `ERR` and similar accessions; every run of an accession is downloaded, paired-end or single-end:

```
SRR2591045
SRR2584867
SRR2591036
```

[SRA Explorer](https://sra-explorer.info) helps to collect accessions. The repository includes an example file `example/sra_accessions.txt`.

### Existing assemblies

`--scaffolds` takes a folder with one subfolder per sample, each holding `scaffolds.fasta` (the SPAdes output layout), or a glob of FASTA files:

```bash
--scaffolds '/path/to/assemblies'            # /path/to/assemblies/<sample>/scaffolds.fasta
--scaffolds '/path/to/assemblies/*.fasta'    # sample name = file name
```

See [Scaffold input](configuration.md#scaffold-input) for how sample names are derived.

### Sample name mapping (optional)

Use a CSV file (no header) to rename samples, e.g. SRR runs to isolate names:

```
SRR1234567,Isolate_A
E1,Iran_2024_A3ghp1i
```

Pass it with `--SRR_sample_map sample_map.csv`. Samples not listed keep their name.

### Query sequences

`--blast_query` is a FASTA with one or more query genes, each searched separately. Queries can be nucleotide (searched with `blastn`) or protein (`tblastn`); the type is detected per sequence. With `--search_tool exonerate`, nucleotide queries must be a **mature CDS** without introns.

Optional companion files:

- `--query_paralog`: known paralogs of the queries, named `<query>_paralog1`, `<query>_paralog2`, … (see [Paralog check](configuration.md#query-search)).
- `--reference_genome`: a reference genome FASTA, to map the reads and measure coverage at the query locus.
- `--repeat_library`: a TE consensus library (FASTA, e.g. from RepeatModeler or [Dfam](https://www.dfam.org)).

---

## Step 6: Run the pipeline

Start the pipeline inside a `tmux` or `screen` session — the Nextflow process must stay alive until completion, even with `-profile slurm`.

```bash
# Run with SLURM, mixed local fastq + SRA input, exact CDS and isoforms
export NCBI_API_KEY=your_ncbi_api_key_here

nextflow run main.nf -config nextflow.config -profile slurm --slurm_queue long \
  -work-dir /scratch/$USER/genomedenovo_work --outdir my_run \
  --NCBI_API_key $NCBI_API_KEY \
  --reads '/path/to/reads/*_{1,2}.fq.gz' \
  --SRA_index './SRA_accessions.txt' \
  --blast_query query.fa --search_tool exonerate \
  --busco_lineage fungi_odb10 --phylo_tree
```

To resume a previous run (if the `work-dir` is still intact), simply add `-resume` to the command. Nextflow will skip completed steps and only execute missing ones.

The pipeline checks its parameters and input files before starting, and stops with an error naming the problem, for example:

- no input given, or a query, reference or repeat library file that does not exist;
- an option that needs another one (`--phylo_tree` without `--busco_lineage`, `--locus_coverage` without `--reference_genome`, `--blast_query` or read input);
- a `--query_paralog` file whose sequences are not named after a query, or a nucleotide query with in-frame stop codons under `--search_tool exonerate`;
- `-profile slurm` without `--slurm_queue`;
- two inputs that end up with the same sample name (checked before assembly, see [Sample names](configuration.md#sample-names)).

!!! warning "Storage"
    De novo assembly of many genomes can require large amounts of temporary storage. Point `-work-dir` to a fast scratch filesystem. The work directory is not deleted automatically, so `-resume` remains possible; delete it yourself once the outputs are checked.

See the [Configuration](configuration.md) page for a full list of all available parameters.
