![alt text](https://github.com/crolllab/genomedenovo_nf/blob/main/logo.png?raw=true)

Pipeline to assemble and query draft genomes
========================================

***This is work in progress and not yet fully tested!***

The `genomedenovo_nf` Nextflow pipeline performs de novo assembly of Illumina paired-end reads, assesses assembly quality, and searches for user-provided query sequences using BLAST. The pipeline supports providing assembled scaffolds, in addition local FASTQ files and downloading from NCBI SRA can be used to start SPAdes assemblies.

Implemented steps:
- `entrez-direct`: query NCBI SRA for metadata (optional)
- `sratools`: download SRA files (optional)
- `fastp`: adapter trimming and read quality filtering
- `spades`: de novo assembly of Illumina paired-end and single-end reads
- `quast`: assembly quality assessment
- `busco`: genome completeness benchmarking against near-universal single-copy orthologs (optional)
- `ncbi-blast+`: BLAST search of user-provided query sequence(s) against assembled genomes
- `mafft`: multiple sequence alignment of extracted BLAST hits
- `exonerate`: spliced alignment and CDS extraction using protein or CDS queries (optional, `--search_tool exonerate`)
- `bwa-mem2`: reference-guided read mapping for coverage assessment (optional, `--reference_genome`)
- `repeatmasker`: transposable element annotation using a user-provided consensus library (optional, `--repeat_library`)
- `raxml-ng` or `iqtree3`: maximum-likelihood phylogenetic tree from BUSCO single-copy orthologs (optional, `--phylo_tree`; RAxML-NG is the default, switch with `--busco_tree iqtree`)
- `iqtree3`: maximum-likelihood phylogenetic trees per query from the BLAST hit alignment and, with `--search_tool exonerate`, from the CDS/protein alignments and the CDS/protein isoforms
- `ape` (R): PDF rendering of the isoform trees (optional, `--search_tool exonerate`)
- haplotype identification: assigns consecutive isoform IDs to unique CDS and protein sequences per query (optional, `--search_tool exonerate`)

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

- `--slurm_queue`: Required when using `-profile slurm`. Specifies the SLURM partition (queue) name on your cluster (e.g. `long`, `normal`, `compute`). The pipeline validates that this is set before starting when SLURM is the active profile.
  ```bash
  -profile slurm --slurm_queue long
  ```

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
  - `exonerate` — BLAST search plus **Exonerate** spliced alignment (`protein2genome` for protein queries, `coding2genome` for nucleotide CDS queries). Exonerate resolves exon boundaries and reports CDS nucleotide and protein sequences per sample. **Important:** when providing a nucleotide query with `--search_tool exonerate`, the sequence must be a **mature CDS** (coding sequence without introns), not genomic DNA.

- `--exonerate_refine`: Optional. How Exonerate refines the heuristic alignment (`region`, `full` or `none`). Default: `region`, which re-aligns the hit region with full dynamic programming so that exons of divergent alleles are not dropped (which would otherwise look like truncations). `full` refines against the whole extracted window; `none` restores the plain heuristic.

- `--skip_quast`: Optional flag. When set, skips the QUAST assembly quality assessment step. Useful when running on pre-assembled scaffolds that have already been assessed, or to speed up BLAST/Exonerate-only runs.

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

- `--keep_full_busco`: Optional flag. By default, BUSCO output is compacted after each sample finishes to reduce disk usage: only `short_summary.*.txt` (for BUSCO metrics aggregation) and extracted single-copy ortholog protein FASTA files used for BUSCO phylogeny are retained. Set this flag to keep the full BUSCO directory structure (all intermediate and log files).
  ```bash
  --keep_full_busco
  ```

- `--reference_genome`: Optional. Path to a reference genome FASTA file. When provided, trimmed reads are mapped to this reference using **BWA-mem2** and the resulting BAM files are sorted and indexed with SAMtools. Useful for estimating mapping rates and coverage uniformity alongside the de novo assemblies. Only available when reads are supplied via `--reads` or `--SRA_index`.
  ```bash
  --reference_genome /path/to/reference.fasta
  ```

- `--repeat_library`: Optional. Path to a TE consensus library FASTA file. When provided, all assembled (or supplied) scaffolds are annotated with **RepeatMasker** using this library. One GFF file per sample is produced, along with masked FASTA files and a summary table. Combined with `--blast_query`, the TEs on the scaffold carrying each sample's hit are listed with their distance to the hit (see output section 13). A pre-built library can be generated with RepeatModeler or obtained from [Dfam](https://www.dfam.org).
  ```bash
  --repeat_library /path/to/TE_library.fasta
  ```

- `--phylo_tree`: Optional flag. When set, a maximum-likelihood phylogenetic tree is built from a **BUSCO gene supermatrix** (requires `--busco_lineage`): single-copy ortholog **protein** sequences (`.faa`, produced by BUSCO v6 via miniprot) are extracted per sample, aligned gene-by-gene with MAFFT (`--amino`), and concatenated into a supermatrix. By default the tree is inferred with **RAxML-NG** using per-gene LG+G partitions and 200 bootstrap replicates. Switch to IQ-tree3 with `--busco_tree iqtree` (per-partition model selection `-m MFP`, 1000 ultrafast bootstrap). Bootstrap is applied only when ≥4 unique taxa are present. Output written to `5_busco_phylogeny/`.

  Threads are set by the active `-profile` configuration (16 by default for `slurm`/`local_highCPU`).
  ```bash
  --busco_lineage dothideomycetes_odb10 --phylo_tree
  --busco_lineage dothideomycetes_odb10 --phylo_tree --busco_tree iqtree
  ```
  The per-query trees (BLAST hits, Exonerate CDS/protein, isoforms) do not depend on this flag: they are always built when `--blast_query` (and, for the Exonerate trees, `--search_tool exonerate`) is set.

- `--busco_tree`: Optional. Phylogenetic tool used to infer the BUSCO supermatrix tree when `--phylo_tree` is set. Default: `raxml-ng`. Use `iqtree` to switch to IQ-tree3.
  ```bash
  --busco_tree raxml-ng   # default: RAxML-NG with LG+G partitions (protein), 200 bootstraps
  --busco_tree iqtree     # alternative: IQ-tree3 with MFP model selection, 1000 UFboot
  ```

- `--busco_subsample`: Optional. Number of BUSCO single-copy ortholog genes to randomly subsample when building the supermatrix for phylogeny inference (`--phylo_tree`). Default: `100`. Use `all` to include every shared single-copy gene without subsampling. Subsampling is performed after filtering for genes present in ≥2 samples, using a random draw (different genes may be selected each run unless you use `-resume`). Useful for reducing computation time with lineages that have thousands of BUSCOs (e.g. `fungi_odb10` has ~758, `bacteria_odb10` has ~124).
  ```bash
  --busco_subsample 100   # default: use 100 randomly selected genes
  --busco_subsample 50    # use only 50 genes (faster, less signal)
  --busco_subsample all   # use all shared single-copy genes
  ```

- `--padding`: Optional. Number of base pairs to extend extracted BLAST hit sequences on both sides. Default: `500`.
  - Example: `--padding 100` extracts hit region plus 100 bp upstream and 100 bp downstream
  - Useful for capturing regulatory regions, promoters, or getting more context around hits
  - Automatically ensures coordinates stay within scaffold boundaries

- `--locus_coverage`: Optional flag. When set, each query sequence is BLASTed against the reference genome (`--reference_genome`) to identify the best-hit locus. Per-base read depth (mapping quality ≥ 20) is then extracted from the sorted BAM files at the padded locus using SAMtools. In parallel, genome-wide depth is summarised in 1 000 bp windows (also MQ ≥ 20) in a single streaming pass, and the mean and standard deviation of window-mean depths are reported. The output table lists locus depth first, followed by the genome-wide window statistics, enabling a direct comparison of local vs. background coverage. Requires `--blast_query`, `--reference_genome`, and read input (`--reads` or `--SRA_index`).
  ```bash
  --locus_coverage
  ```
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

- `--SRR_sample_map`: Optional. A CSV file mapping sample IDs (as derived from FASTQ filenames, SRA accessions, or scaffold directory names) to user-defined names. This is useful for renaming samples to meaningful isolate or strain identifiers before output is written. The CSV format is: `original_id,new_name` (one mapping per line, no header). If a sample ID is not listed in the file, the original ID is kept unchanged. Default: `false` (no mapping). Final sample names must be unique, also after characters outside `[A-Za-z0-9._-]` are replaced by `_` (e.g. `isolate 1` and `isolate_1` clash). The pipeline checks this for all inputs (reads, SRA runs, scaffolds) before assembly and stops with a list of the clashing inputs.

  Example `sample_map.csv`:
  ```
  E1,Iran_2024_A3ghp1i
  E2,Iran_2024_A3ghp3b
  SRR1234567,MyIsolate_A
  ```

## Step 3: Run `genomedenovo_nf`

For a fully worked example using *E. coli* LTEE data, see [example/README.md](example/README.md). It covers all major input modes (local FASTQ, SRA accessions, Exonerate, BUSCO phylogeny, reference mapping, locus coverage) with ready-to-run commands.

A minimal run using local FASTQ files:

```bash
micromamba activate nf_gp_env
nextflow run main.nf \
  --reads '/path/to/reads/*_{1,2}.fq.gz' \
  --blast_query './query.fa' \
  --outdir './nf_output'
```

For SLURM clusters with SRA downloads:

```bash
NCBI_API_KEY=abcdef1234567890   # substitute your own key
micromamba activate nf_gp_env
nextflow run main.nf -config nextflow.config -profile slurm \
  -work-dir '/scratch/nf_tmp' --outdir './my_nf_run_output' \
  --NCBI_API_key $NCBI_API_KEY \
  --reads '/path/to/reads/*_{1,2}.fq.gz' \
  --SRA_index './SRA_accessions.txt' \
  --blast_query './query.fa'
```

### Multi-sample BLAST analysis on pre-assembled genomes

Process multiple pre-assembled genomes and extract BLAST hits:

```bash
nextflow run main.nf -config nextflow.config \
  -profile local_highCPU \
  --outdir './blast_output' \
  --scaffolds '/path/to/assemblies' \
  --blast_query './query.fa'
```

Notes:
- Enter a `tmux` or `screen` session before running — the session must remain active for the duration of the pipeline
- In `local` and `local_highCPU` modes the pipeline runs on the local machine; use `-profile slurm` to distribute across cluster nodes

## Step 4: Pipeline output

The pipeline will produce the following output files in the specified `--outdir`:

### SRA download outputs (if using --SRA_index)
- `NCBI_SRR_PE_accessions.txt`: List of paired-end SRA accessions
- `NCBI_SRR_SE_accessions.txt`: List of single-end SRA accessions
- `NCBI_download_urls.tsv`: Download URLs and sources (ENA vs NCBI) for each accession

### 1 — Read QC outputs (if using --reads or --SRA_index)
- `1_fastq_stats/fastp_summary.tsv`: Table with read quality metrics before and after trimming for all samples (columns prefixed with `before_` / `after_`)
- `1_fastq_stats/fastp_summary.txt`: Human-readable formatted version of the same metrics

### 2 — Assembly outputs (if using --reads or --SRA_index)
- `2_assemblies/{sample_id}/scaffolds.fasta`: SPAdes scaffold file per sample
- `2_assemblies/{sample_id}/contigs.fasta`: SPAdes contig file per sample
- `2_assemblies/{sample_id}/spades.log`: SPAdes log per sample

### 3 — Assembly QC outputs
- `3_quast/assembly_summary_metrics.tsv`: Combined QUAST metrics table for all samples
- `3_quast/{sample_id}_quast/`: Per-sample QUAST report directory

### 4 — BUSCO outputs (if using --busco_lineage)
- `4_busco/{sample_id}_busco/`: Per-sample BUSCO result directory (compact by default; full when `--keep_full_busco` is set)
- `4_busco/short_summary.*.{sample_id}_busco.txt`: BUSCO short summary per sample
- `4_busco/busco_summary_metrics.tsv`: Combined table with completeness metrics for all samples (Complete %, Single-copy %, Duplicated %, Fragmented %, Missing %)

Default compact mode removes large temporary/intermediate BUSCO outputs (for example `run_*/hmmer_output`, `run_*/miniprot_output`, `run_*/busco_sequences`, and temporary files) after extracting required downstream artifacts.

### 5 — BUSCO gene phylogeny outputs (if using --phylo_tree and --busco_lineage)
- `5_busco_phylogeny/`: Contents depend on `--busco_tree`:

  **Default — RAxML-NG** (`--busco_tree raxml-ng`):
  - `busco_phylogeny.raxml.bestTree`: Best-scoring ML tree in Newick format
  - `busco_phylogeny.raxml.support`: Best tree annotated with bootstrap support values
  - `busco_phylogeny.raxml.bestModel`: Best-fit substitution model per partition
  - `busco_phylogeny.raxml.bootstraps`: All bootstrap trees
  - `busco_phylogeny.raxml.mlTrees`: All ML trees from independent searches
  - `busco_phylogeny.raxml.log`: RAxML-NG log file

  **Alternative — IQ-tree3** (`--busco_tree iqtree`):
  - `busco_phylogeny.treefile`: Best-scoring ML tree in Newick format
  - `busco_phylogeny.contree`: Consensus tree with ultrafast bootstrap support values
  - `busco_phylogeny.iqtree`: Full IQ-tree report (model selection, tree statistics)
  - `busco_phylogeny.log`: IQ-tree log file
  - `busco_phylogeny.best_scheme.nex`: Best partition scheme used

### 6 — BLAST outputs (if using --blast_query)
When multiple query sequences are provided, outputs are organized by query ID:
- `6_blast_hits/{query_id}/blast_results/{sample_id}_blast_results.tsv`: Raw BLAST output in tabular format
- `6_blast_hits/{query_id}/hit_sequences/{sample_id}_top_hits.fasta`: Extracted hit sequences per sample
- `6_blast_hits/{query_id}/hit_sequences/{sample_id}_top_hits_info.tsv`: Metadata table for each hit
- `6_blast_hits/{query_id}/{query_id}_all_blast_hits.fasta`: Combined FASTA with **query sequence first**, then all hits; hits are named by sample ID (renamed IDs when `--SRR_sample_map` is used, whitespace replaced by `_`) with the hit location as description (`>sample_id scaffold:start-end`), so alignments and trees show sample names
- `6_blast_hits/{query_id}/{query_id}_all_blast_hits_summary.tsv`: Combined metadata table (hits only)
- **`6_blast_hits/blast_summary_merged.tsv`**: Wide-format table merging results from all queries
  - First column: Sample_ID; for each query: `{query}_Present`, `{query}_Identity`, `{query}_AlignLen`, `{query}_Evalue`, `{query}_BitScore`

### 7 — Exonerate outputs (if using --search_tool exonerate)
- `7_exonerate_hits/{query_id}/gff3_results/{sample_id}_exonerate.gff3`: Raw Exonerate target GFF per sample
- `7_exonerate_hits/{query_id}/hit_sequences/{sample_id}_exonerate_cds.fasta`: Extracted CDS nucleotide sequences
- `7_exonerate_hits/{query_id}/hit_sequences/{sample_id}_exonerate_protein.fasta`: Extracted protein sequences
- `7_exonerate_hits/{query_id}/hit_sequences/{sample_id}_exonerate_summary.tsv`: Per-sample summary of exon coordinates and lengths
- `7_exonerate_hits/{query_id}/{query_id}_all_exonerate_cds.fasta`: Combined CDS FASTA across samples
- `7_exonerate_hits/{query_id}/{query_id}_all_exonerate_protein.fasta`: Combined protein FASTA across samples
- `7_exonerate_hits/{query_id}/{query_id}_all_exonerate_summary.tsv`: Combined per-query summary table
- `7_exonerate_hits/exonerate_summary_merged.tsv`: Wide-format table merging Exonerate results across all queries and samples

### 8 — Alignment outputs
- `8_alignments/blast/{query_id}/{query_id}_all_blast_hits_aligned_sequences.fasta`: MAFFT alignment of all BLAST hits per query
- `8_alignments/blast/{query_id}/{query_id}_all_blast_hits_mafft_alignment_log.txt`: MAFFT alignment log and statistics
- `8_alignments/exonerate/{query_id}/{query_id}_all_exonerate_cds_aligned_sequences.fasta`: CDS alignment (if `--search_tool exonerate`)
- `8_alignments/exonerate/{query_id}/{query_id}_all_exonerate_protein_aligned_sequences.fasta`: Protein alignment (if `--search_tool exonerate`)

### 9 — Alignment phylogeny outputs (if using --blast_query)
Per-query IQ-tree3 ML trees are built from every alignment with at least 3 unique sequences; model selection (`-m MFP`) is always applied, ultrafast bootstrap (`-B 1000`) when ≥4 unique sequences are present. They do not require `--phylo_tree`.
- `9_phylogenies/blast/{query_id}/{query_id}_blast_phylogeny.*`: Tree of the BLAST hit alignment
- `9_phylogenies/exonerate/{query_id}/` (if using `--search_tool exonerate`): one tree each for the CDS and protein alignments (one tip per sample sequence) and for the CDS and protein isoforms (one tip per isoform):
  - `{query_id}_cds_phylogeny.treefile`: Best ML tree (CDS) in Newick format
  - `{query_id}_cds_phylogeny.contree`: Consensus tree with bootstrap support (if ≥4 unique taxa)
  - `{query_id}_cds_phylogeny.iqtree`: Full IQ-tree report (model selection, tree statistics)
  - `{query_id}_cds_phylogeny.log`: IQ-tree log
  - `{query_id}_protein_phylogeny.treefile`: Best ML tree (protein) in Newick format
  - `{query_id}_protein_phylogeny.contree`: Consensus tree with bootstrap support (if ≥4 unique taxa)
  - `{query_id}_protein_phylogeny.iqtree`: Full IQ-tree report
  - `{query_id}_protein_phylogeny.log`: IQ-tree log
  - `{query_id}_cds_isoform_phylogeny.*` / `{query_id}_protein_isoform_phylogeny.*`: The same IQ-tree files for the isoform trees, built from `13_summary_data/{query_id}_{cds,protein}_isoforms_aligned.fasta`
  - **`{query_id}_cds_isoform_phylogeny.pdf` / `{query_id}_protein_isoform_phylogeny.pdf`**: Isoform tree drawn as a PDF. The tree is midpoint-rooted, tips read `Isoform N (n = samples carrying it)`, isoforms found only in truncated (orange) or stop-lost (blue) sequences are coloured, the query reference is marked with a triangle, and bootstrap support is shown at internal nodes (from the `.contree` when available). Not produced when a query has fewer than 3 isoforms.

### 10 — BWA-mem2 mapping outputs (if using --reference_genome)
- `10_bwa_mapping/{sample_id}_sorted.bam`: Sorted, MAPQ≥10 filtered BAM file
- `10_bwa_mapping/{sample_id}_sorted.bam.bai`: BAM index
- `10_bwa_mapping/{sample_id}_flagstat.json`: SAMtools flagstat report with mapping statistics
- `10_bwa_mapping/bwa_mapping_summary.tsv`: Combined mapping statistics across all samples

### 11 — Locus coverage outputs (if using --locus_coverage)
- `11_locus_coverage/{query_id}/{query_id_safe}_ref_locus.tsv`: Best BLAST hit of the query on the reference genome (scaffold, coordinates)
- `11_locus_coverage/{query_id}/{sample_id}_{query_id_safe}_locus_coverage.tsv`: Per-sample coverage table with columns:
  - `scaffold`, `hit_start`, `hit_end`: BLAST hit coordinates on the reference
  - `locus_start`, `locus_end`: Padded locus boundaries (±`--padding` bp)
  - `locus_mean_depth`, `locus_median_depth`: Read depth at the padded locus (MQ ≥ 20)
  - `genome_mean_depth`, `genome_sd_depth`: Mean and SD of per-base read depth across the whole genome (MQ ≥ 20)
- `11_locus_coverage/{query_id}/{query_id_safe}_all_locus_coverage.tsv`: Combined table across all samples for this query

### 12 — RepeatMasker outputs (if using --repeat_library)
- `12_repeatmasker/{sample_id}_repeatmasker.gff`: TE annotations in GFF format
- `12_repeatmasker/{sample_id}_repeatmasker.masked`: Soft-masked genome FASTA
- `12_repeatmasker/{sample_id}_repeatmasker.tbl`: RepeatMasker summary table
- `12_repeatmasker/{sample_id}_repeatmasker.out`: RepeatMasker annotation table (includes the TE class/family missing from the GFF)

### 13 — Summary data

#### TEs near the hit (if using --repeat_library with --blast_query)
For each sample, the RepeatMasker TEs (simple repeats and low-complexity regions excluded) on the scaffold carrying the query hit. The hit is the Exonerate gene model with `--search_tool exonerate`, otherwise the best BLAST hit (the merged span of all HSPs on the best scaffold, so its ends are not necessarily the start/stop codons).

Distances are edge to edge in bp and follow the **orientation of the hit (5'→3')**: negative = the TE lies upstream (5') of the hit, positive = downstream (3'), 0 = the TE overlaps the hit span (introns included). For a hit on the minus strand, upstream therefore means higher scaffold coordinates.

- **`13_summary_data/{query_id}_scaffold_TEs.tsv`**: One row per TE on the hit scaffold, ordered 5'→3' within each hit. Columns: `sample_id`, `hit_id` (`{sample}_{n}` when a sample has several Exonerate gene models), `query_id`, `scaffold`, `hit_start`, `hit_end`, `hit_strand`, `te_name`, `te_class`, `te_start`, `te_end`, `te_strand`, `te_orientation` (`same` / `opposite` strand as the hit), `te_length`, `position` (`upstream`, `downstream`, `overlap`; `none` when the scaffold carries no TE, `no_hit` for samples without a hit), `distance_bp`, `perc_div` (divergence from the consensus), `sw_score`, `te_id` (RepeatMasker ID).
- `13_summary_data/{query_id}_scaffold_TEs_summary.tsv`: One row per hit: number of TEs on the scaffold, overlapping TEs, and the nearest upstream and downstream TE with their distances.

#### Isoforms / haplotypes (if using --search_tool exonerate)
Isoform numbers are shared by all files of a query and follow **isoform frequency**: isoforms carried by at least one intact sequence come first, ordered by the number of sample sequences carrying them (most frequent = 1; ties keep alignment order; the query reference is not counted). They are followed by the query reference if its sequence is not found in any intact sample, then by isoforms found only in truncated or stop-lost sequences, again ordered by frequency. The intact-only table therefore numbers consecutively, and truncated isoforms continue the series.

- `13_summary_data/{query_id}_isoforms.tsv`: CDS and protein isoform IDs of the query reference and every intact sequence (truncated and stop-lost sequences excluded). Columns: `sample_id`, `cds_isoform`, `protein_isoform`, `cds_sequence`, `protein_sequence` (ungapped).
- **`13_summary_data/{query_id}_isoforms_incl_truncated.tsv`**: Every retrieved sequence (`sample_id` matches the Exonerate FASTA headers, `{sample}_{n}` when a sample has several gene models) with its isoform IDs and full-length/truncated classification. Samples without an Exonerate alignment are listed with status `no_hit`. Key columns:
  - `status`: `reference`, `full_length`, `truncated`, `extended` (stop codon lost only) or `no_hit`
  - `truncation_cause`: `nucleotide_deletion` (part of the query is absent from the contig), `ORF_mutation` (the sequence is there but the reading frame is disrupted), both, or `none`
  - `deletion_types`: `5prime_deletion` / `3prime_deletion` (segment absent, contig continues with unrelated sequence), `contig_end` (contig ends before the segment could fit), `assembly_gap` (N run), `internal_deletion` (in-frame deletion ≥ 10 aa)
  - `orf_mutation_types`: `frameshift`, `premature_stop`, `start_lost`, `stop_lost`, `noncanonical_splice_site`, `segment_outside_ORF` (the unaligned query segment re-aligns to the contig outside the ORF, e.g. an unbridged frameshift or a retained intron)
  - `contig`, `contig_length`, `strand`, `gene_start`, `gene_end`, `exon_count`: location of the gene model
  - `aligned_query_start_aa` / `aligned_query_end_aa`, `query_coverage_pct`, `aa_identity_pct`: extent and similarity of the alignment
  - `observed_protein_length`: length of the protein the contig actually encodes (from the first aligned codon to the first stop, frameshifted tail included); `in_frame_query_aa` / `in_frame_query_pct`: reference residues encoded in the correct frame before that stop
  - `start_codon`, `stop_codon`, `n_frameshifts`, `n_premature_stops`, `n_noncanonical_splice_sites`, `n_inframe_indels`, `first_disruption_aa`
  - `mutation_details`: readable description of every mutation (e.g. `p.Q46* nonsense codon TAG`, `p.C36fs frameshift -1 nt (deletion in contig)`, `C-terminal query aa 69-82 (14 aa, 42 nt) missing: segment absent; ...`)
  - `notes`: e.g. termini that Exonerate left unaligned although the ORF is intact, telomeric repeats or assembly gaps in the flank, low-confidence frameshift calls
  - `cds_sequence`, `protein_sequence`: ungapped sequences (`NA` for `no_hit`)
- **`13_summary_data/{query_id}_isoform_mutations.tsv`**: One row per mutation event with the isoform IDs of the sequence carrying it, the event class (`ORF_mutation`, `nucleotide_deletion`, `in_frame_indel`), type, affected query residues, contig coordinates, reference and observed codon/residue, indel size and description.
- `13_summary_data/{query_id}_cds_isoforms.fasta` / `{query_id}_protein_isoforms.fasta`: One ungapped sequence per unique CDS / protein isoform (truncated ones included), named `cds_isoform_N` / `protein_isoform_N` to match the TSVs; the header also gives the number of samples carrying it (`n_samples=`) and their status (`status=`).
- `13_summary_data/{query_id}_cds_isoforms_aligned.fasta` / `{query_id}_protein_isoforms_aligned.fasta`: One aligned row per isoform (taken from the MAFFT alignment, all-gap columns removed); input for the isoform trees in `9_phylogenies/exonerate/{query_id}/`.


## Description of pipeline steps

- **Hybrid SRA Download** (optional): The pipeline uses a hybrid approach to download FASTQ files from SRA accessions. It first attempts fast direct downloads from ENA (European Nucleotide Archive), then automatically falls back to NCBI if needed. This significantly improves download speed while maintaining reliability. The `--SRA_index` parameter accepts a file with one accession per line. See `FEATURE_hybrid_sra_download.md` for details.

- **Assembly with SPAdes**: The pipeline assembles Illumina paired-end and single-end reads using SPAdes with default parameters optimized for fungal genomes. Each sample is processed independently.

- **Assembly quality assessment with QUAST**: Each assembled genome is evaluated using QUAST, producing both per-sample reports and a combined summary table.

- **BLAST search and hit extraction**: Each assembled genome is searched against a user-provided query sequence using BLAST (automatically selects blastn for nucleotide queries or tblastn for protein queries). The top hits are extracted with optional flanking sequence padding. Combined FASTA files include the original query sequence first, followed by all BLAST hits, making them ready for downstream analysis like multiple sequence alignment.

- **Exonerate spliced alignment and CDS/protein extraction** (optional, `--search_tool exonerate`): BLAST hit regions are used as anchors for Exonerate (`protein2genome` or `coding2genome`) to refine exon boundaries and recover coding sequences per sample. Raw target-side GFF output is parsed into CDS nucleotide FASTA, translated protein FASTA, and per-sample/per-query summary tables.

- **CDS integrity classification** (optional, `--search_tool exonerate`): each Exonerate gene model is walked codon by codon against the query to call frameshifts, premature stop codons, lost start/stop codons, non-canonical splice sites and in-frame indels. Query segments that Exonerate left unaligned are checked on the full contig: an in-frame start/stop codon at the expected distance means a divergent but intact terminus; otherwise the segment is re-aligned to the flanking sequence (local alignment for short segments, Exonerate `est2genome`/`protein2genome` for longer ones). A segment found there is reported as an ORF mutation (the sequence is present but outside the reading frame); a segment not found is a nucleotide-level deletion, with contig ends, assembly gaps and telomeric repeats flagged.

- **Multiple sequence alignment with MAFFT**: The extracted BLAST hit sequences from each sample are aligned using MAFFT with `--auto` and `--adjustdirection` options. This produces both per-sample alignments and a combined alignment of all sequences from all samples.

- **Reference-guided read mapping with BWA-mem2** (optional, `--reference_genome`): Trimmed reads are mapped to the provided reference genome using BWA-mem2. Alignments are filtered at MAPQ≥10, sorted, and indexed with SAMtools. Flagstat reports are generated per sample.

- **Locus coverage analysis** (optional, `--locus_coverage`): Each query sequence is BLASTed against the reference genome to find its best-hit locus. For every sorted BAM file, SAMtools `depth` (MQ ≥ 20, `-a` for zero-coverage positions) is run at the padded locus to report mean and median per-base depth. Genome-wide background depth is assessed by computing the mean and SD of per-base depth across all positions in a single streaming pass (using Welford's computational variance formula). The combined output tables list locus depth alongside genome-wide per-bp depth statistics for each sample and query.

- **Genome completeness with BUSCO** (optional, `--busco_lineage`): Each assembled genome is evaluated against a specified OrthoDB lineage database using BUSCO v6 in genome mode. Results include a per-sample completeness score and a combined `busco_summary_metrics.tsv` table. If `--phylo_tree` is also set, single-copy ortholog nucleotide sequences are additionally extracted from each sample for downstream phylogenomics.

- **Transposable element annotation with RepeatMasker** (optional, `--repeat_library`): All scaffolds (assembled or user-supplied) are annotated with RepeatMasker using the provided TE consensus library. Outputs include GFF annotations, soft-masked FASTA files, and per-sample summary tables. With `--blast_query`, the TEs on each hit scaffold are listed with their position relative to the hit in the hit's 5'→3' orientation.

- **BUSCO gene phylogeny** (optional, `--phylo_tree`): Single-copy BUSCO ortholog sequences are extracted per sample, aligned gene-by-gene with MAFFT, and concatenated into a supermatrix. By default, **RAxML-NG** infers the tree with per-gene LG+G partitions and 200 bootstrap replicates. Switch to IQ-tree3 with `--busco_tree iqtree` (per-partition MFP model selection, 1000 ultrafast bootstrap). Bootstrap is skipped when fewer than 4 unique taxa are present. Output written to `5_busco_phylogeny/`. Requires `--busco_lineage`.

- **Per-query phylogenies with IQ-tree3** (with `--blast_query`; Exonerate trees with `--search_tool exonerate`): ML trees are built per query from the BLAST hit alignment and from the MAFFT-aligned CDS and protein sequences produced by Exonerate, plus isoform trees with one tip per unique CDS / protein isoform, which are also drawn as midpoint-rooted PDFs with R/ape. Model selection (`-m MFP`) and ultrafast bootstrap (`-B 1000`) are applied when ≥4 unique sequences are present. Outputs written to `9_phylogenies/{blast,exonerate}/{query_id}/`.

- **Isoform / haplotype identification** (optional, `--search_tool exonerate`): After alignment, unique CDS and protein sequences are identified and numbered per query by frequency (most frequent intact isoform = 1). The resulting table (`13_summary_data/{query_id}_isoforms.tsv`) allows quick assessment of allelic or strain-level diversity across samples; truncated isoforms are listed with their mutations in `{query_id}_isoforms_incl_truncated.tsv`.


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

