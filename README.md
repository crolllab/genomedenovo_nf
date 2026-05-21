![alt text](https://github.com/crolllab/genomedenovo_nf/blob/main/logo.png?raw=true)

Pipeline to assemble and query draft genomes
========================================

***This is work in progress and not yet fully tested!***

The `genomedenovo_nf` Nextflow pipeline performs de novo assembly of Illumina paired-end reads, assesses assembly quality, and searches for user-provided query sequences using BLAST. The pipeline supports providing assembled scaffolds, in addition local FASTQ files and downloading from NCBI SRA can be used to start SPAdes assemblies.

Implemented steps:
- `entrez-direct`: query NCBI SRA for metadata (optional)
- `sratools`: download SRA files (optional)
- `spades`: de novo assembly of Illumina paired-end and single-end reads
- `quast`: assembly quality assessment
- `ncbi-blast+`: BLAST search of user-provided query sequence(s) against assembled genomes
- `mafft`: multiple sequence alignment of extracted BLAST hits

## Step 1: Repository, singularity containers and nextflow environment

Cloning the repository

```bash
git clone git@github.com:crolllab/genomedenovo_nf.git
cd genomedenovo_nf
```

Create a conda environment for Nextflow

```bash
micromamba create -n nf_gp_env -c conda-forge -c bioconda nextflow
micromamba activate nf_gp_env
```

Alternatively, if your HPC provides Nextflow as an environment module:

```bash
module load Nextflow
```

### Singularity containers

Singularity is required to run the pipeline. All container images are pulled automatically by Nextflow on first use and cached in `~/.singularity/cache`. No manual setup is needed.

## Step 2: Configure `genomedenovo_nf` options  

Available parameters

- `-profile`: Optional. `local` skips the SLURM queuing system and runs tasks with up to 4 threads each. `local_highCPU` runs tasks with up to 24 threads each. `slurm` submits all tasks as SLURM jobs. Default: `slurm`.

- `--outdir`: Optional. Folder to save final output files. Default: `./nf_output`.

- `-work-dir`: Optional. Defines where to store temporary files (often many TB). Consider `/scratch/work` for large datasets. Default: `./work`. 

- `-resume`: Optional. If set, the pipeline will resume from the last completed step, skipping already completed steps. The `work-dir` needs to be intact for this.

- `--NCBI_API_key`: Required for querying and downloading from NCBI SRA. You can get your key by creating an [account on NCBI](https://account.ncbi.nlm.nih.gov/). After registration/login, find on the top right the link to the "Account settings". Click on "Create API key" and copy it.

- `--blast_query`: Optional. Provide a FASTA file containing query sequence(s) for BLAST search against assembled genomes. The pipeline automatically detects whether queries are **nucleotide** or **protein** sequences:
  - **Nucleotide queries** (DNA/RNA): Uses `blastn` for direct nucleotide comparison
  - **Protein queries** (amino acids): Uses `tblastn` to search translated genome sequences
  - Supports multi-sequence FASTA files - each sequence is processed independently
  - Detection: Sequences with ≥95% ATGCUN characters are classified as nucleotide, otherwise protein

- `--search_tool`: Optional. Choose the sequence search tool. Default: `blast`.
  - `blast` — BLAST search only (blastn for nucleotide queries, tblastn for protein queries).
  - `gth` — BLAST search plus **GenomeThreader** spliced alignment. GTH resolves exact exon/intron boundaries and reports CDS nucleotide and protein sequences per sample. **Important:** when providing a nucleotide query with `--search_tool gth`, the sequence must be a **mature CDS** (coding sequence without introns), not genomic DNA. GTH treats the nucleotide query as an EST/cDNA and will produce incorrect results if the query contains introns. Protein queries are always safe to use with GTH.

- `--skip_quast`: Optional flag. When set, skips the QUAST assembly quality assessment step. Useful when running on pre-assembled scaffolds that have already been assessed, or to speed up BLAST/GTH-only runs.

- `--busco_lineage`: Optional. BUSCO lineage dataset for genome completeness assessment using [BUSCO v6](https://busco.ezlab.org/). **If this parameter is not provided, BUSCO is skipped entirely.** The lineage determines which set of near-universal single-copy orthologs (OrthoDB) is used to score the assembly. Select the most specific lineage available for your organism:

  | Lineage | Scope | Recommended for |
  |---|---|---|
  | `bacteria_odb10` | All Bacteria | Broad screening across any bacterial species |
  | `fungi_odb10` | All Fungi | Broad screening across any fungal species |
  | `ascomycota_odb10` | Ascomycota | Most plant pathogens, model yeasts |
  | `dothideomycetes_odb10` | Dothideomycetes | *Zymoseptoria*, *Alternaria*, *Cladosporium*, etc. |

  A full list of available lineages can be found at [https://busco-data.ezlab.org/v5/data/lineages/](https://busco-data.ezlab.org/v5/data/lineages/).

  Example:
  ```bash
  --busco_lineage dothideomycetes_odb10
  ```

- `--busco_db_path`: Optional. Path to a directory containing pre-downloaded BUSCO lineage databases. When provided, BUSCO runs in offline mode and will not attempt to download any data. If omitted, BUSCO automatically downloads the requested lineage on first use and caches it in `~/.busco_db`. Pre-download a lineage with:
  ```bash
  busco --download dothideomycetes_odb10 --download_path /path/to/busco_db
  ```

- `--padding`: Optional. Number of base pairs to extend extracted BLAST hit sequences on both sides. Default: `500`.
  - Example: `--padding 100` extracts hit region plus 100 bp upstream and 100 bp downstream
  - Useful for capturing regulatory regions, promoters, or getting more context around hits
  - Automatically ensures coordinates stay within scaffold boundaries
  - Padded coordinates shown in FASTA headers; original BLAST coordinates saved in TSV files

- `--reads`: Optional. Provide the path to the folder containing the fastq read files. The pipeline will automatically find paired-end read files based on the naming convention. Must be bracketed by single quotes `'`. See below for examples. Note: Local file input currently supports paired-end reads only. Single-end reads are supported via `--SRA_index`.

  Multiple glob patterns can be combined by passing a **semicolon-separated list** (wrap the whole value in single quotes). A semicolon is used as the separator to avoid conflicts with brace expansions like `{1,2}` in glob patterns:
  ```bash
  --reads '/path/to/dir1/**_{1,2}.fq.gz;/path/to/dir2/**_{1,2}.fq.gz'
  ```

- `--SRA_index`: Optional. Instead of local fastq files, you can provide a file listing NCBI SRA accessions (or ENA, etc.) with one accession per line including `SRR...`, `SRP...`, `SRX...`, `PRNJ...`, etc. If the accession includes multiple `SRR...` runs, all included runs are processed. The pipeline automatically handles both paired-end and single-end SRA reads. You can use the [SRA Explorer](https://sra-explorer.info) to collect accession ids. 

- `--scaffolds`: Optional. Instead of assembling reads, you can provide already assembled draft genomes in fasta format. The pipeline will then skip the assembly step and proceed with the quality assessment and BLAST analysis.

  Multiple sources can be combined by passing a **semicolon-separated list** of paths or glob patterns (wrap the whole value in single quotes):
  ```bash
  --scaffolds '/path/to/dir1;/path/to/dir2'
  --scaffolds '/path/to/dir1;/legserv/NGS_data/**/scaffolds.fasta'
  ```

  **Three input modes supported (combinable via semicolon-separated list):**
  1. **Directory with multiple samples**: Provide a directory path where each subdirectory contains a `scaffolds.fasta` file. The subdirectory name will be used as the sample ID.
     ```bash
     --scaffolds '/path/to/assemblies'
     # Expected structure:
     #   assemblies/
     #     sample1/scaffolds.fasta
     #     sample2/scaffolds.fasta
     #     sample3/scaffolds.fasta
     ```

  2. **Recursive glob pattern**: Provide a glob pattern (with `**`) to match `scaffolds.fasta` files across a directory tree. The parent directory name will be used as the sample ID.
     ```bash
     --scaffolds '/legserv/NGS_data/**/scaffolds.fasta'
     ```

  3. **Direct file pattern**: Provide a glob pattern to match specific fasta files. The base filename will be used as the sample ID.
     ```bash
     --scaffolds '/path/to/folder/*.fasta'
     ```

- `--SRR_sample_map`: Optional. A CSV file mapping sample IDs (as derived from FASTQ filenames, SRA accessions, or scaffold directory names) to user-defined names. This is useful for renaming samples to meaningful isolate or strain identifiers before output is written. The CSV format is: `original_id,new_name` (one mapping per line, no header). If a sample ID is not listed in the file, the original ID is kept unchanged. Default: `false` (no mapping).

  Example `sample_map.csv`:
  ```
  E1,Iran_2024_A3ghp1i
  E2,Iran_2024_A3ghp3b
  SRR1234567,MyIsolate_A
  ```

## Step 3: Run `genomedenovo_nf` - example options

The below example is based on the _Zymoseptoria tritici_ Illumina paired-end reads from local file servers and NCBI SRA. Substitute NCBI accessions and local read paths with your own data.

### Select sets of local fastq file pairs

The following examples show different ways to select paired-end fastq files. The `--reads` option must be bracketed by single quotes `'`. Multiple glob patterns can be combined with a semicolon-separated list. The examples assume that the read files are named according to the Illumina convention, i.e. ending with `_1.fq.gz` and `_2.fq.gz` or `_R1.fq.gz` and `_R2.fq.gz`. Adjust the patterns according to your own naming conventions.

```bash
# Option 1 - select all fastq files ending with _1.fq.gz and _2.fq.gz
--reads '/path/to/reads/*{1,2}.fq.gz'

# Option 2 - select all fastq files ending with _1.fq.gz and _2.fq.gz, including all subdirectories (** instead of *)!
--reads '/path/to/reads/**{1,2}.fq.gz'

# Option 3 - select all files ending with _1.fq.gz and _2.fq.gz OR ending with _R1.fq.gz and _R2.fq.gz
--reads '/path/to/reads/*_{,R}{1,2}.fq.gz'

# Option 4 - select all files (including all subdirectories), with optional variation in fq/fastq, 1/R1 (2/R2), optional _001 or _001_ additions
--reads '/path/to/reads/**_{,R}{1,2}{,_001,_001_*}.{fq,fastq}.gz'
````

```bash
# Option 4 will include these files among others:
ST01IR_A48b.cleanData_1.fq.gz and ST01IR_A48b.cleanData_2.fq.gz
ST01IR_A26b.cleanData_R1.fq.gz and ST01IR_A26b.cleanData_R2.fq.gz
ST01IR_A26b.cleanData_R1.fastq.gz and ST01IR_A26b.cleanData_R2.fastq.gz
ST01IR_A26b.cleanData_R1_001.fastq.gz and ST01IR_A26b.cleanData_R2_001.fastq.gz
J9_L2_R1_001_18ku2CAeFgfk.fastq.gz and J9_L2_R2_001_j2kKKZcCX6h0.fastq.gz
```

### Define NCBI SRA accessions

Example file to provide for the `--SRA_index ...` option. 

```bash
ERR13824484
ERR13824571
ERR13824499
```

Save the text file e.g. as `SRA_accessions.txt`. The repository includes an example file.


## Step 4: `genomedenovo_nf` run example

```bash
# read selection example 1 (small set) - selects 9 paired-end reads from LEGserv
READS='/legserv/NGS_data/Zymoseptoria/Illumina_DNAseq/Croll_2013/ST99CH_*_{1,2}.fq.gz'
# read selection example 2 (very large set) - selects nearly all available paired-end DNA-seq datasets on LEGserv
READS='/legserv/NGS_data/Zymoseptoria/Illumina_DNAseq/*_{,R}{1,2}{,_001,_001_*}.{fq,fastq}.gz'
```

Start the pipeline using `slurm` and processing local fastq files and NCBI accessions. Temporary files are written to `/scratch` and the output will be in the `my_nf_run_output` folder. The `SRA_accessions.txt` example file is included in the repository.

```bash
# substitute with your own NCBI API key!
NCBI_API_KEY=abcdef1234567890
# activate conda environment
micromamba activate nf_gp_env
# run nextflow pipeline
nextflow run main.nf -config nextflow.config -profile slurm \
  -work-dir '/scratch/nf_tmp' --outdir './my_nf_run_output' \
  --NCBI_API_key $NCBI_API_KEY \
  --reads $READS \
  --SRA_index './SRA_accessions.txt' \
  --blast_query './query.fa' \
```

### Example: Multi-sample BLAST analysis on pre-assembled genomes

Process multiple pre-assembled genomes and extract BLAST hits:

```bash
# Directory containing multiple sample folders, each with scaffolds.fasta
ASSEMBLIES='/path/to/assemblies'
QUERY='./wagA_gene.fa'

micromamba activate nf_gp_env

nextflow run main.nf -config nextflow.config \
  -profile local_highCPU \
  -work-dir '/scratch/nf_tmp' \
  --outdir './blast_output' \
  --scaffolds $ASSEMBLIES \
  --blast_query $QUERY
```

This will:
1. Process all `scaffolds.fasta` files found in subdirectories of `$ASSEMBLIES`
2. Run BLAST search for each sample
3. Extract hit sequences from each sample
4. Combine all hits into `all_blast_hits.fasta` and `all_blast_hits_summary.tsv`

Notes on the execution:
- Before executing the `nextflow` command, enter e.g. a `tmux` or `screen` session. The session needs to remain active until the end of the pipeline (even if you specify the `slurm` option)
- In `local` and `local_highCPU` modes, the pipeline will run on the local machine. Use the `slurm` option to spread the load to all available nodes (requires SLURM).


## Step 5: Pipeline output

The pipeline will produce the following output files in the specified `--outdir`:

### SRA download outputs (if using --SRA_index)
- `NCBI_SRR_PE_accessions.txt`: List of paired-end SRA accessions
- `NCBI_SRR_SE_accessions.txt`: List of single-end SRA accessions
- `NCBI_download_urls.tsv`: Download URLs and sources (ENA vs NCBI) for each accession

### Assembly outputs
- `assembly_stats.tsv`: Table with assembly statistics for all samples
- `assembly_summary_metrics.tsv`: Table with QUAST metrics for all samples
- `quast_report.html`: Interactive HTML report from QUAST

### BLAST outputs (per query)
When multiple query sequences are provided, outputs are organized by query ID:
- `{query_id}/blast_results/`: Folder containing raw BLAST results for each sample
  - `{sample_id}_blast_results.tsv`: Raw BLAST output in tabular format
- `{query_id}/blast_hit_sequences/`: Folder containing extracted hit sequences
  - `{sample_id}_top_hits.fasta`: FASTA file with extracted hit sequences
  - `{sample_id}_top_hits_info.tsv`: Metadata table for each hit
- `{query_id}/{query_id}_all_blast_hits.fasta`: Combined FASTA file with **query sequence first**, then all hits for this query
- `{query_id}/{query_id}_all_blast_hits_summary.tsv`: Combined metadata table for this query (hits only, not query)

### Merged BLAST summary
- **`blast_summary_merged.tsv`**: Wide-format table merging results from all queries
  - First column: Sample_ID
  - For each query: 5 essential columns with query-prefixed names:
    - `{query}_Present`: TRUE/FALSE flag indicating if hit was found
    - `{query}_Identity`: Percent identity of alignment
    - `{query}_AlignLen`: Alignment length in bp
    - `{query}_Evalue`: E-value
    - `{query}_BitScore`: Bit score
  - Example: With queries wagA, wagB, wagC → 1 + (5 × 3) = 16 columns total
  - Detailed information (Query ID, Subject ID, coordinates) available in per-query files

### Alignment outputs
- `mafft_alignments/`: Folder containing per-query multiple sequence alignments
  - `{query_id}_aligned.fasta`: MAFFT alignment for all BLAST hits for this query
  - `{query_id}_alignment_log.txt`: MAFFT alignment log and statistics


## Description of pipeline steps

- **Hybrid SRA Download** (optional): The pipeline uses a hybrid approach to download FASTQ files from SRA accessions. It first attempts fast direct downloads from ENA (European Nucleotide Archive), then automatically falls back to NCBI if needed. This significantly improves download speed while maintaining reliability. The `--SRA_index` parameter accepts a file with one accession per line. See `FEATURE_hybrid_sra_download.md` for details.

- **Assembly with SPAdes**: The pipeline assembles Illumina paired-end and single-end reads using SPAdes with default parameters optimized for fungal genomes. Each sample is processed independently.

- **Assembly quality assessment with QUAST**: Each assembled genome is evaluated using QUAST, producing both per-sample reports and a combined summary table.

- **BLAST search and hit extraction**: Each assembled genome is searched against a user-provided query sequence using BLAST (automatically selects blastn for nucleotide queries or tblastn for protein queries). The top hits are extracted with optional flanking sequence padding. Combined FASTA files include the original query sequence first, followed by all BLAST hits, making them ready for downstream analysis like multiple sequence alignment.

- **Multiple sequence alignment with MAFFT**: The extracted BLAST hit sequences from each sample are aligned using MAFFT with `--auto` and `--adjustdirection` options. This produces both per-sample alignments and a combined alignment of all sequences from all samples.


## Features to consider / bug fixes


## Utilities

### Download SRA files manually

Use conda/micromamba to install the `sratools` package, which includes the `fastq-dump` command.

Create the following file e.g. `SRA_download.sh` with one accession per line:

```bash
#!/bin/bash
fastq-dump --split-files --gzip SRR24910574
fastq-dump --split-files --gzip SRR24910575
fastq-dump --split-files --gzip SRR25074049
fastq-dump --split-files --gzip SRR25074050
...
```

Parallelize the download using GNU parallel:

```bash
parallel -j 10 < SRA_download.sh
```

NB: Don't try to download too many files at once, as NCBI will stall you.

Proceed with `genomedenovo_nf` as described above, using the `--reads` option to point to the download folder.

