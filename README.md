![alt text](https://github.com/crolllab/genomedenovo_nf/blob/main/logo.png?raw=true)

Pipeline to assemble and query draft genomes
========================================

***This is work in progress and not yet fully tested!***

`genomedenovo_nf` is a Nextflow pipeline that assembles Illumina genomes de novo (or takes existing assemblies), checks assembly quality, and searches every genome for one or more query genes. For each query it extracts the hits, aligns them, builds trees and, optionally, recovers the exact coding sequences, classifies them as intact or truncated, and groups them into isoforms.

## Contents

- [Pipeline steps](#pipeline-steps)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Options](#options)
- [Outputs](#outputs)
- [Methods](#methods)
- [Utilities](#utilities)

## Pipeline steps

| Step | Tool | Runs when |
|---|---|---|
| SRA metadata and download | `entrez-direct`, `sra-tools` (ENA first, NCBI fallback) | `--SRA_index` |
| Read trimming and QC | `fastp` | reads given (`--reads` / `--SRA_index`) |
| De novo assembly | `SPAdes` (`--isolate`) | reads given |
| Assembly QC | `QUAST` | always, unless `--skip_quast` |
| Genome completeness | `BUSCO` v6 | `--busco_lineage` |
| BUSCO species tree | `MAFFT` + `RAxML-NG` (or `IQ-TREE`) | `--phylo_tree` |
| Query search and hit extraction | `BLAST+` (`blastn` / `tblastn`) | `--blast_query` |
| Spliced alignment, CDS and protein | `Exonerate` | `--search_tool exonerate` |
| CDS integrity (full-length / truncated) | custom, `Exonerate` | `--search_tool exonerate` |
| Isoforms and isoform trees | custom, `IQ-TREE`, R `ape` | `--search_tool exonerate` |
| Per-query alignments and trees | `MAFFT`, `IQ-TREE` | `--blast_query` |
| Reference mapping | `BWA-mem2`, `SAMtools` | `--reference_genome` (reads given) |
| Query locus coverage | `BLAST+`, `SAMtools` | `--locus_coverage` |
| TE annotation and TEs near hits | `RepeatMasker` | `--repeat_library` |

## Installation

```bash
git clone git@github.com:crolllab/genomedenovo_nf.git
cd genomedenovo_nf

# Nextflow (>= 25) in a conda environment ...
micromamba create -n nf_gp_env -c conda-forge -c bioconda nextflow
micromamba activate nf_gp_env
# ... or as an HPC module
module load Nextflow
```

Singularity (or Apptainer) must be available. All tools run in containers that Nextflow pulls on first use and caches in `~/.singularity/cache`; nothing else needs to be installed.

## Quick start

A worked example with *E. coli* data, covering every input mode and option, is in [example/README.md](example/README.md).

```bash
# Assemble local reads and search them for a query
nextflow run main.nf -profile local \
  --reads '/path/to/reads/*_{1,2}.fq.gz' \
  --blast_query query.fa

# Search existing assemblies, with exact CDS, isoforms and TEs near the hits
nextflow run main.nf -profile local_highCPU \
  --scaffolds '/path/to/assemblies' \
  --blast_query query.fa --search_tool exonerate \
  --repeat_library TE_library.fasta

# SRA accessions on a SLURM cluster
nextflow run main.nf -profile slurm --slurm_queue long \
  -work-dir /scratch/nf_work --outdir my_run \
  --SRA_index SRA_accessions.txt --NCBI_API_key "$NCBI_API_KEY" \
  --blast_query query.fa --busco_lineage fungi_odb10 --phylo_tree
```

Run the pipeline inside `tmux` or `screen`: the Nextflow process must stay alive until the run ends. Add `-resume` to restart a run from its last completed steps.

## Options

Nextflow options take one dash (`-profile`), pipeline options two (`--reads`). Optional steps are switched on by setting their option; everything else has a default.

### Input

At least one of `--reads`, `--SRA_index` or `--scaffolds` is required; they can be combined.

| Option | Default | Description |
|---|---|---|
| `--reads` | – | Local paired-end FASTQ files as a quoted glob, e.g. `'/data/*_{1,2}.fq.gz'`. Several patterns: separate with `;`. |
| `--SRA_index` | – | Text file with one SRA/ENA accession per line (`SRR`, `SRX`, `SRP`, `PRJNA`, …). Every run of an accession is downloaded; paired-end and single-end runs are supported. [SRA Explorer](https://sra-explorer.info) helps to collect accessions. |
| `--NCBI_API_key` | – | NCBI API key for the SRA lookup. Optional but recommended: it raises NCBI's request limit. Must be valid, or the lookup fails. Create one under *Account settings* of your [NCBI account](https://account.ncbi.nlm.nih.gov/). |
| `--scaffolds` | – | Existing assemblies (FASTA); assembly is skipped for them. See [scaffold input](#scaffold-input). Several sources: separate with `;`. |
| `--SRR_sample_map` | – | CSV `original_id,new_name` (no header) to rename samples. See [sample names](#sample-names). |

#### Sample names

Sample names come from the input:

| Input | Sample name |
|---|---|
| `--reads` | Common prefix of the two read files, without trailing `_`, `.` or `-` (`S1_1.fq.gz` + `S1_2.fq.gz` → `S1`; note that `S1_R1`/`S1_R2` gives `S1_R`) |
| `--SRA_index` | SRR run accession |
| `--scaffolds` | Sample folder or file name (see below) |

`--SRR_sample_map` renames them; samples not listed keep their name:

```
SRR1234567,Isolate_A
E1,Iran_2024_A3ghp1i
```

Final names must be unique, also after characters outside `[A-Za-z0-9._-]` are replaced by `_` (`isolate 1` and `isolate_1` clash). The pipeline checks this before assembly and stops with a list of the clashing inputs.

#### Scaffold input

| `--scaffolds` value | Sample name |
|---|---|
| Folder with one subfolder per sample, each holding `scaffolds.fasta` (or `.fa`/`.fna`): `'/path/to/assemblies'` | Subfolder name |
| Recursive glob: `'/data/**/scaffolds.fasta'` | Parent folder name |
| File glob: `'/path/*.fasta'` | File name without extension (parent folder name if the file is called `scaffolds`, `assembly`, `genome` or `contigs`) |

### Query search

| Option | Default | Description |
|---|---|---|
| `--blast_query` | – | FASTA with one or more query sequences, each searched separately. Nucleotide queries (≥ 95 % `ACGTUN`) use `blastn`, protein queries `tblastn`. |
| `--search_tool` | `blast` | `blast`: BLAST only. `exonerate`: BLAST, then Exonerate spliced alignment for exact CDS and protein, CDS integrity and isoforms. Nucleotide queries must then be a **mature CDS** (no introns); the run stops if the query has in-frame stop codons. |
| `--exonerate_refine` | `region` | Exonerate alignment refinement: `region` (re-align the hit region with full dynamic programming, so exons of divergent alleles are not lost), `full` (refine on the whole extracted window) or `none`. |
| `--padding` | `500` | bp added on both sides of the BLAST hit when extracting it. The padded region is also the window searched by Exonerate and the locus used by `--locus_coverage`. |

### Assembly QC

| Option | Default | Description |
|---|---|---|
| `--skip_quast` | off | Skip QUAST. |
| `--busco_lineage` | – | Run BUSCO with this OrthoDB lineage, e.g. `bacteria_odb10`, `fungi_odb10`, `ascomycota_odb10`, `dothideomycetes_odb10` ([all lineages](https://busco-data.ezlab.org/v5/data/lineages/)). Choose the most specific lineage for your organism. |
| `--busco_db_path` | – | Folder with pre-downloaded lineages (`busco --download <lineage> --download_path <folder>`); BUSCO then runs offline. Without it, the lineage is downloaded once into `<outdir>/busco_db_cache/`. |
| `--keep_full_busco` | off | Keep the complete BUSCO output folders. By default only the short summaries and the single-copy protein sequences are kept. |

### BUSCO species tree

| Option | Default | Description |
|---|---|---|
| `--phylo_tree` | off | Build a species tree from the BUSCO single-copy genes. Requires `--busco_lineage`. |
| `--busco_tree` | `raxml-ng` | `raxml-ng` or `iqtree`. |
| `--busco_subsample` | `100` | Number of BUSCO genes randomly drawn for the supermatrix, or `all`. |

The per-query trees (BLAST hits, Exonerate CDS/protein, isoforms) do not need `--phylo_tree`: they are always built with `--blast_query`.

### Reference mapping and locus coverage

| Option | Default | Description |
|---|---|---|
| `--reference_genome` | – | Map the trimmed reads to this reference with BWA-mem2. Needs read input. |
| `--locus_coverage` | off | Compare read depth at the query locus of the reference with the genome-wide depth. Needs `--reference_genome`, `--blast_query` and read input. |

### Repeats

| Option | Default | Description |
|---|---|---|
| `--repeat_library` | – | TE consensus library (FASTA, e.g. from RepeatModeler or [Dfam](https://www.dfam.org)). All assemblies are annotated with RepeatMasker; with `--blast_query`, the TEs near each query hit are reported. |

### Execution

| Option | Default | Description |
|---|---|---|
| `--outdir` | `./nf_output` | Output folder. |
| `-profile` | none | `local`: local machine, 4 CPUs per task. `local_highCPU`: local machine, 16 CPUs per task. `slurm`: submit each task as a SLURM job (16 CPUs). Without a profile, tasks run locally with their own defaults (mostly 1 CPU). |
| `--slurm_queue` | – | SLURM partition; required with `-profile slurm`. The partition should allow long jobs (up to 7 days for large assemblies). |
| `-work-dir` | `./work` | Folder for intermediate files, which can grow very large; use scratch space for big runs. |
| `-resume` | off | Reuse completed tasks from a previous run (needs the intact work folder). |

## Outputs

Folders are numbered in pipeline order; a folder exists only if its step ran. `{query_id}` is the query FASTA header (characters outside `[A-Za-z0-9._-]` replaced by `_` in file names).

### Run overview (`<outdir>/`)
- `pipeline_execution_stats.tsv` / `.txt`: Run time and resources per step
- `pipeline_trace.txt`: Nextflow trace of every task
- `NCBI_SRR_PE_accessions.txt`, `NCBI_SRR_SE_accessions.txt`, `NCBI_download_urls.tsv`: Resolved SRA runs and their download sources (`--SRA_index`)

### 1 — Read QC (`1_fastq_stats/`)
- `fastp_summary.tsv` / `.txt`: Read metrics before and after trimming for all samples
- `{sample_id}_{PE,SE}_fastp.json`: fastp report per sample

### 2 — Assemblies (`2_assemblies/{sample_id}/`)
- `scaffolds.fasta`, `contigs.fasta`, `spades.log`

### 3 — Assembly QC (`3_quast/`)
- `assembly_summary_metrics.tsv`: QUAST metrics for all samples
- `{sample_id}_quast/`: Full QUAST report per sample

### 4 — BUSCO (`4_busco/`)
- `busco_summary_metrics.tsv`: Complete, single-copy, duplicated, fragmented and missing % for all samples
- `short_summary.*.{sample_id}_busco.txt`, `{sample_id}_busco/`: Per-sample results (compact unless `--keep_full_busco`)

### 5 — BUSCO species tree (`5_busco_phylogeny/`)
- RAxML-NG: `busco_phylogeny.raxml.support` (best tree with bootstrap support), `.bestTree`, `.bestModel`, `.bootstraps`, `.mlTrees`, `.log`
- IQ-TREE: `busco_phylogeny.contree` (tree with bootstrap support), `.treefile`, `.iqtree` (report), `.best_scheme.nex`, `.log`

### 6 — BLAST hits (`6_blast_hits/`)
- **`blast_summary_merged.tsv`**: One row per sample; per query `{query}_Present`, `_Identity`, `_AlignLen`, `_Evalue`, `_BitScore`
- `query_sequences/`: The split queries and their detected type
- `{query_id}/{query_id}_all_blast_hits.fasta`: Query followed by the hit of every sample, named `>sample_id scaffold:start-end`
- `{query_id}/{query_id}_all_blast_hits_summary.tsv`: Hit coordinates, identity, e-value and bit score per sample (`no_hits_found` without a hit)
- `{query_id}/blast_results/{sample_id}_blast_results.tsv`: Raw BLAST tabular output
- `{query_id}/hit_sequences/{sample_id}_top_hits.fasta` / `_top_hits_info.tsv`: Extracted (padded) hit region per sample, named by its scaffold location

### 7 — Exonerate (`7_exonerate_hits/`, with `--search_tool exonerate`)
- `exonerate_summary_merged.tsv`: Exonerate results for all queries and samples
- `{query_id}/{query_id}_all_exonerate_cds.fasta` / `_protein.fasta`: CDS and protein of every sample (and the query), named by sample ID (`{sample}_{n}` when a sample has several gene models)
- `{query_id}/{query_id}_all_exonerate_summary.tsv`: Scaffold, strand, coordinates, exon count and lengths per gene model
- `{query_id}/hit_sequences/`, `{query_id}/gff3_results/`: The same per sample, and the raw Exonerate GFF

### 8 — Alignments (`8_alignments/`)
- `blast/{query_id}/{query_id}_all_blast_hits_aligned_sequences.fasta`: MAFFT alignment of the BLAST hits
- `exonerate/{query_id}/{query_id}_all_exonerate_{cds,protein}_aligned_sequences.fasta`: CDS and protein alignments
- `*_mafft_alignment_log.txt`: MAFFT log, including the sequences MAFFT reverse-complemented

### 9 — Per-query trees (`9_phylogenies/`)
IQ-TREE trees, built when an alignment has at least 3 unique sequences. Each tree comes as `.treefile` (ML tree), `.contree` (bootstrap support, with ≥ 4 unique sequences), `.iqtree` (report) and `.log`.
- `blast/{query_id}/{query_id}_blast_phylogeny.*`: BLAST hits, one tip per sample
- `exonerate/{query_id}/{query_id}_{cds,protein}_phylogeny.*`: Exonerate CDS / protein, one tip per sample
- `exonerate/{query_id}/{query_id}_{cds,protein}_isoform_phylogeny.*`: CDS / protein isoforms, one tip per isoform
- **`exonerate/{query_id}/{query_id}_{cds,protein}_isoform_phylogeny.pdf`**: The isoform trees drawn midpoint-rooted. Tips read `Isoform N (n = samples)`; isoforms found only in truncated (orange) or stop-lost (blue) sequences are coloured, the query reference is marked with a triangle, and bootstrap support is shown at the nodes.

### 10 — Reference mapping (`10_bwa_mapping/`)
- `bwa_mapping_summary.tsv` / `.txt`: Mapping statistics for all samples
- `{sample_id}_sorted.bam` (+ `.bai`): Sorted alignments, MAPQ ≥ 10
- `{sample_id}_flagstat.json`: SAMtools flagstat per sample

### 11 — Locus coverage (`11_locus_coverage/{query_id}/`)
- `{query_id}_all_locus_coverage.tsv`: All samples; per-sample tables `{sample_id}_{query_id}_locus_coverage.tsv`. Columns: `scaffold`, `hit_start`, `hit_end` (query hit on the reference), `locus_start`, `locus_end` (hit ± `--padding`), `locus_mean_depth`, `locus_median_depth`, `genome_mean_depth`, `genome_sd_depth`
- `{query_id}_ref_locus.tsv`: Best hit of the query on the reference

### 12 — RepeatMasker (`12_repeatmasker/`)
- `{sample_id}_repeatmasker.out`: Annotation table with TE class/family
- `{sample_id}_repeatmasker.gff`: Annotations as GFF (no class/family)
- `{sample_id}_repeatmasker.masked`, `.tbl`: Masked assembly and summary

### 13 — Summary tables (`13_summary_data/`)

#### TEs near the hit (`--repeat_library` with `--blast_query`)
The hit is the Exonerate gene model with `--search_tool exonerate`, otherwise the BLAST hit. Distances are edge to edge in bp, in the hit's own **5'→3' orientation**: negative = TE upstream (5') of the hit, positive = downstream (3'), 0 = TE overlaps the hit (introns included). On a minus-strand hit, upstream therefore means higher scaffold coordinates.

- **`{query_id}_scaffold_TEs.tsv`**: One row per TE on the hit scaffold, ordered 5'→3'. Columns: `sample_id`, `hit_id` (`{sample}_{n}` for several gene models), `query_id`, `scaffold`, `hit_start`, `hit_end`, `hit_strand`, `te_name`, `te_class`, `te_start`, `te_end`, `te_strand`, `te_orientation` (`same`/`opposite` strand as the hit), `te_length`, `position` (`upstream`, `downstream`, `overlap`; `none` = no TE on the scaffold, `no_hit` = sample without hit), `distance_bp`, `perc_div` (divergence from the consensus), `sw_score`, `te_id`
- `{query_id}_scaffold_TEs_summary.tsv`: One row per hit: number of TEs on the scaffold, overlapping TEs, nearest upstream and downstream TE with distance

#### Isoforms (`--search_tool exonerate`)
Isoform numbers are shared by all files of a query (see [isoform numbering](#isoforms)).

- `{query_id}_isoforms.tsv`: The query reference and every intact sequence. Columns: `sample_id`, `cds_isoform`, `protein_isoform`, `cds_sequence`, `protein_sequence`
- **`{query_id}_isoforms_incl_truncated.tsv`**: Every sequence, including truncated ones and samples without a hit (`no_hit`). Columns:
  - `sample_id`, `cds_isoform`, `protein_isoform`, `status` (`reference`, `full_length`, `truncated`, `extended` = stop codon lost only, `no_hit`)
  - `truncation_cause`: `nucleotide_deletion` (part of the query is missing from the contig), `ORF_mutation` (the sequence is there but the reading frame is disrupted), both, or `none`
  - `deletion_types`: `5prime_deletion` / `3prime_deletion` (segment missing, contig continues with unrelated sequence), `contig_end`, `assembly_gap` (N run), `internal_deletion` (in-frame deletion ≥ 10 aa)
  - `orf_mutation_types`: `frameshift`, `premature_stop`, `start_lost`, `stop_lost`, `noncanonical_splice_site`, `segment_outside_ORF` (the unaligned query segment is found on the contig outside the reading frame, e.g. an unbridged frameshift or a retained intron)
  - `contig`, `contig_length`, `strand`, `gene_start`, `gene_end`, `exon_count`: gene model location
  - `aligned_query_start_aa`, `aligned_query_end_aa`, `query_coverage_pct`, `aa_identity_pct`: alignment extent and similarity
  - `observed_protein_length`: protein actually encoded (first aligned codon to first stop); `in_frame_query_aa` / `in_frame_query_pct`: query residues encoded in frame before that stop
  - `start_codon`, `stop_codon`, `n_frameshifts`, `n_premature_stops`, `n_noncanonical_splice_sites`, `n_inframe_indels`, `first_disruption_aa`
  - `mutation_details`: readable list of mutations (e.g. `p.Q46* nonsense codon TAG`, `p.C36fs frameshift -1 nt (deletion in contig)`)
  - `notes`: e.g. termini Exonerate left unaligned although the ORF is intact, telomeric repeats or assembly gaps in the flank, low-confidence frameshift calls
  - `cds_sequence`, `protein_sequence`
- **`{query_id}_isoform_mutations.tsv`**: One row per mutation with the isoform IDs of the sequence carrying it, event class (`ORF_mutation`, `nucleotide_deletion`, `in_frame_indel`), type, query residues, contig coordinates, reference and observed codon/residue, indel size and description
- `{query_id}_{cds,protein}_isoforms.fasta`: One sequence per isoform, named `cds_isoform_N` / `protein_isoform_N`, with the number of samples (`n_samples=`) and their status (`status=`) in the header
- `{query_id}_{cds,protein}_isoforms_aligned.fasta`: The isoforms as aligned rows; input of the isoform trees

## Methods

**SRA download.** Accessions are resolved to runs with `entrez-direct`. Each run is downloaded from ENA when available, with `sra-tools` from NCBI as fallback.

**Reads and assembly.** Reads are trimmed with `fastp` (default settings) and assembled per sample with `SPAdes --isolate`.

**Assembly QC.** QUAST reports contiguity metrics. BUSCO v6 (genome mode, miniprot) scores completeness against the chosen lineage.

**BUSCO species tree.** The single-copy BUSCO proteins found in at least 2 samples are subsampled (`--busco_subsample`), aligned per gene with MAFFT and concatenated into a partitioned supermatrix. RAxML-NG uses an LG+G model per gene and 200 bootstraps. IQ-TREE selects a model per partition (`-m MFP`) with 1000 ultrafast bootstraps. Bootstrapping needs at least 4 taxa (4 unique sequences for IQ-TREE).

**BLAST search.** Each query is searched against every assembly (`blastn -task blastn` or `tblastn`). Per sample, the scaffold with the best bit score is kept, and all hits of the query on that scaffold are merged into one span. That way, the separate hits of a multi-exon gene make up one region. The span ± `--padding` is extracted in the orientation of the query.

**Exonerate CDS.** Exonerate aligns the query to the extracted region (`protein2genome` for protein queries, `coding2genome` for CDS queries) to find the exon boundaries. The aligned exons give the CDS, which is translated to protein. Where gene models overlap, the best-scoring one is kept, and antisense models of nucleotide queries are discarded.

**CDS integrity.** Each gene model is compared codon by codon with the query to find frameshifts, premature stops, lost start/stop codons, non-canonical splice sites and in-frame indels. Query segments that Exonerate left unaligned are looked up on the full contig:
- An in-frame start or stop codon at the expected distance means a divergent but intact end.
- Otherwise the segment is re-aligned to the flanking sequence. If it is found, the sequence is present but outside the reading frame: an ORF mutation. If it is not found, it is a nucleotide-level deletion. Contig ends, assembly gaps and telomeric repeats are flagged.

<a id="isoforms"></a>**Isoforms.** Identical aligned CDS (or protein) sequences form one isoform. Numbering follows frequency:
1. Isoforms carried by at least one intact sequence, ordered by the number of sample sequences carrying them; the most frequent is 1. Ties keep alignment order, and the query reference is not counted.
2. The query reference, if no intact sample shares its sequence.
3. Isoforms found only in truncated or stop-lost sequences, again by frequency.

The intact-only table therefore numbers consecutively.

**Alignments and per-query trees.** BLAST hits, Exonerate CDS and Exonerate proteins are aligned per query with MAFFT (`--auto --adjustdirection`). For protein queries, the query is left out of the BLAST hit alignment, since the hits are nucleotide. IQ-TREE builds a tree from every alignment with ≥ 3 unique sequences (`-m MFP`; 1000 ultrafast bootstraps with ≥ 4 unique sequences). Isoform trees use one aligned row per isoform and are drawn as PDFs with R/`ape`.

**Reference mapping and locus coverage.** Trimmed reads are mapped with BWA-mem2, filtered at MAPQ ≥ 10, sorted and indexed. For `--locus_coverage`, the query is BLASTed against the reference to find its locus. Per sample, the mean and median depth over the locus ± `--padding` (MAPQ ≥ 20) are compared with the mean and SD of per-base depth over the whole genome.

**TEs.** RepeatMasker annotates every assembly with the given library. TEs on the scaffold of each query hit are reported with their position relative to the hit (see [section 13](#13--summary-tables-13_summary_data)); simple repeats and low-complexity regions are left out.

## Utilities

### Download SRA files manually

Instead of `--SRA_index`, runs can be downloaded beforehand with `sra-tools` and passed with `--reads`:

```bash
# SRA_download.sh: one command per accession
fastq-dump --split-files --gzip SRR24910574
fastq-dump --split-files --gzip SRR24910575

parallel -j 10 < SRA_download.sh   # NCBI throttles larger numbers of parallel downloads
```

Planned features are listed in [FUTURE_FEATURES.md](FUTURE_FEATURES.md).
