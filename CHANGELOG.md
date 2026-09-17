# Changelog

## v0.9.0 — 2026-09-17

First tagged pre-release. The pipeline is functional but not yet fully tested; see `FUTURE_FEATURES.md` for planned work.

### Features

- **Input options.** Start from local paired-end FASTQ files, NCBI SRA accessions (metadata via `entrez-direct`, download via `sratools`/ENA), or pre-assembled scaffolds.
- **Read QC and assembly.** `fastp` trimming and quality filtering with a per-sample read statistics table; `SPAdes` de novo assembly (`--isolate`).
- **Assembly QC.** `QUAST` metrics and optional `BUSCO` completeness (`--busco_lineage`, `--busco_subsample`).
- **Query search.** `BLAST` (blastn/tblastn, auto-detected from the query FASTA) with hit extraction and `MAFFT` alignments; optional `Exonerate` spliced alignment with CDS/protein extraction (`--search_tool exonerate`).
- **Isoform / haplotype identification.** Consecutive isoform IDs for unique CDS and protein sequences per query.
- **Phylogenetics.** Maximum-likelihood tree from BUSCO single-copy orthologs (`--phylo_tree`; RAxML-NG default, IQ-TREE 3 via `--busco_tree iqtree`) and per-query IQ-TREE trees from Exonerate alignments.
- **Reference mapping.** Optional `bwa-mem2` mapping to a reference genome (`--reference_genome`) with % mapped reads and locus coverage relative to the genomic background (`--locus_coverage`).
- **Repeat annotation.** Optional `RepeatMasker` with a user-provided consensus library (`--repeat_library`).

### Infrastructure

- Numbered output directories following the order of the pipeline.
- Cluster-agnostic SLURM profile (`--slurm_queue`), plus `local` and `local_highCPU` profiles.
- Standardised retry and error strategies with incremental memory on retry.
- Sample IDs sanitised for tools that reject spaces or special characters in paths.
- Nextflow 26 language compatibility.
