![alt text](https://github.com/nextflow-io/nextflow/raw/master/docs/_static/nextflow-logo-bg-light.png)

Pipeline to assemble and query draft genomes
========================================

The `genomedenovo_nf` Nextflow pipeline performs ....

Implemented steps:
- `entrez-direct`: query NCBI SRA for metadata (optional)
- `sratools`: download SRA files (optional)
- `spades`: denovo assembly of Illumina paired-end reads
- `quast`: assembly quality assessment

Current limitations:


## Step 1: Repository, singularity containers and nextflow environment

Cloning the repository, 
```bash
git clone git@github.com:crolllab/genomepanel_nf.git
cd genomepanel_nf
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

- `--scaffolds`: Optional. Instead of assembling reads, you can provide already assembled draft genomes in fasta format. The pipeline will then skip the assembly step and proceed with the quality assessment and BLAST analysis. You can provide multiple fasta files using wildcards, e.g. `/path/to/folder/*.fasta`. The fasta files must have unique names and the file names will be used to name the output files.

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

Notes on the exection:
- Before executing the `nextflow` command, enter e.g. a `tmux` or `screen` session. The session needs to remain active until the end of the pipeline (even if you specify the `slurm` option)
- In `local` and `local_highCPU` modes, the pipeline will run on the local machine. Use the `slurm` option to spread the load to all available nodes (requires SLURM).


## Step 5: Pipeline output

The pipeline will produce the following output files in the specified `--outdir`:


## Description of pipeline steps



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

