![alt text](https://github.com/nextflow-io/nextflow/raw/master/docs/_static/nextflow-logo-bg-light.png)

Pipeline to assemble and query draft genomes
========================================

# EARLY DRAFT - NOT MUCH FULLY WORKING YET!!


The `genomedenovo_nf` Nextflow pipeline performs ....

Implemented steps:
- `entrez-direct`: query NCBI SRA for metadata (optional)
- `sratools`: download SRA files (optional)
- `spades`: denovo assembly of Illumina paired-end reads
- `quast`: assembly quality assessment


## Step 1: Repository, singularity containers and nextflow environment

Cloning the repository, 
```bash
git clone git@github.com:crolllab/genomedenovo_nf.git
cd genomedenovo_nf
```

Create a conda environment for Nextflow

```bash
micromamba create -n nf_gp_env
micromamba activate nf_gp_env
micromamba install -c bioconda nextflow
```

### Get singularity images  
The `singularity` folder must be in the same directory as the `main.nf` file.

A copy of the compatible images is on LEGserv.

```bash
rsync -va /legserv/Temp/Shared/genomedenovo_nf/singularity .
```

Alternatively, you can pull the images directly from the Galaxy Project depot.

```bash
mkdir -p singularity
cd singularity
# entrez-direct
singularity pull https://depot.galaxyproject.org/singularity/entrez-direct:24.0--he881be0_0
# sratools
singularity pull https://depot.galaxyproject.org/singularity/sra-tools%3A3.2.1--h4304569_1
# SPAdes
singularity pull https://depot.galaxyproject.org/singularity/spades%3A4.2.0--haf24da9_0
# QUAST
singularity pull https://depot.galaxyproject.org/singularity/quast:5.3.0--py39pl5321heaaa4ec_0
# NCBI BLAST+
singularity pull https://depot.galaxyproject.org/singularity/blast%3A2.17.0--h66d330f_0
# MAFFT
singularity pull https://depot.galaxyproject.org/singularity/mafft%3A7.525--h031d066_1
```

## Step 2: Configure `genomedenovo_nf` options  

Available parameters

- `-profile`: Optional. `local` skips the SLURM queuing system and runs tasks with up to 4 threads each. `local_highCPU` runs tasks with up to 24 threads each. `slurm` submits all tasks as SLURM jobs. Default: `slurm`.

- `--outdir`: Optional. Folder to save final output files. Default: `./nf_output`.

- `-work-dir`: Optional. Defines where to store temporary files (often many TB). Consider `/scratch/work` for large datasets. Default: `./work`. 

- `-resume`: Optional. If set, the pipeline will resume from the last completed step, skipping already completed steps. The `work-dir` needs to be intact for this.

- `--NCBI_API_key`: Required for querying and downloading from NCBI SRA. You can get your key by creating an [account on NCBI](https://account.ncbi.nlm.nih.gov/). After registration/login, find on the top right the link to the "Account settings". Click on "Create API key" and copy it.

- `--reads`: Optional. Provide the path to the folder containing the fastq read files. The pipeline will automatically find all paired read files based on the naming convention. Must be bracketed by single quotes `'`. See below for examples. Important: accepts only paired-end reads.

- `--SRA_index`: Optional. Instead of local fastq files, you can provide a file listing NCBI SRA accessions (or ENA, etc.) with one accession per line including `SRR...`, `SRP...`, `SRX...`, `PRNJ...`, etc. If the accession includes multiple `SRR...` runs, all included runs are processed. You can use the [SRA Explorer](https://sra-explorer.info) to collect accession ids. 

- `--blast_query`: Optional. Provide a FASTA file containing query sequence(s) for BLAST search against assembled genomes. The pipeline automatically detects whether queries are **nucleotide** or **protein** sequences:
  - **Nucleotide queries** (DNA/RNA): Uses `blastn` for direct nucleotide comparison
  - **Protein queries** (amino acids): Uses `tblastn` to search translated genome sequences
  - Supports multi-sequence FASTA files - each sequence is processed independently
  - Detection: Sequences with ≥95% ATGCUN characters are classified as nucleotide, otherwise protein
  - See `FEATURE_blast_tool_selection.md` for details on detection logic

- `--padding`: Optional. Number of base pairs to extend extracted BLAST hit sequences on both sides. Default: `0` (extract exact hit region only).
  - Example: `--padding 100` extracts hit region plus 100 bp upstream and 100 bp downstream
  - Useful for capturing regulatory regions, promoters, or getting more context around hits
  - Automatically ensures coordinates stay within scaffold boundaries
  - Padded coordinates shown in FASTA headers; original BLAST coordinates saved in TSV files

- `--scaffolds`: Optional. Instead of assembling reads, you can provide already assembled draft genomes in fasta format. The pipeline will then skip the assembly step and proceed with the quality assessment and BLAST analysis. 
  
  **Two input modes supported:**
  1. **Directory with multiple samples**: Provide a directory path where each subdirectory contains a `scaffolds.fasta` file. The subdirectory name will be used as the sample ID.
     ```bash
     --scaffolds '/path/to/assemblies'
     # Expected structure:
     #   assemblies/
     #     sample1/scaffolds.fasta
     #     sample2/scaffolds.fasta
     #     sample3/scaffolds.fasta
     ```
  
  2. **Direct file pattern**: Provide a glob pattern to match specific fasta files. The base filename will be used as the sample ID.
     ```bash
     --scaffolds '/path/to/folder/*.fasta'
     ```

  Important:
  - Accepts only paired-end reads locally. Both SE and PE are accepted from SRA.
  - If the download from SRA produces errors (connection reset, etc.), an alternative is to download the files separately, and use the `--reads` option to specify the downloaded files. See below for an example of how to download SRA files manually.

## Step 3: Run `genomedenovo_nf` - example options

The below example is based on the _Zymoseptoria tritici_ Illumina paired-end reads from local file servers and NCBI SRA. Substitute NCBI accessions and local read paths with your own data.

### Select sets of local fastq file pairs

The following examples show different ways to select paired-end fastq files. The `--reads` option must be bracketed by single quotes `'`. The examples assume that the read files are named according to the Illumina convention, i.e. ending with `_1.fq.gz` and `_2.fq.gz` or `_R1.fq.gz` and `_R2.fq.gz`. Adjust the patterns according to your own naming conventions.

Note the use of `

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
READS='/legserv/NGS_data/Zymoseptoria/Illumina_DNAseq/_{,R}{1,2}{,_001,_001_*}.{fq,fastq}.gz'
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
- `{query_id}/{query_id}_all_blast_hits.fasta`: Combined FASTA file with all hits for this query
- `{query_id}/{query_id}_all_blast_hits_summary.tsv`: Combined metadata table for this query

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

- **Download from NCBI SRA** (optional): Using `entrez-direct` and `sratools`, the pipeline can download FASTQ files from NCBI SRA based on provided accession IDs. The `--SRA_index` parameter accepts a file with one accession per line.

- **Assembly with SPAdes**: The pipeline assembles Illumina paired-end reads using SPAdes with default parameters optimized for bacterial genomes. Each sample is processed independently.

- **Assembly quality assessment with QUAST**: Each assembled genome is evaluated using QUAST, producing both per-sample reports and a combined summary table.

- **BLAST search and hit extraction**: Each assembled genome is searched against a user-provided query sequence using BLASTn. The top hits are extracted, and both per-sample and combined FASTA files and metadata tables are generated.

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

Proceed with `genomepanel_nf` as described above, using the `--reads` option to point to the download folder.

