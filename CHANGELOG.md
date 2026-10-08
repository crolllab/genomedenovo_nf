# Changelog

## Unreleased

### Added

- **Truncated isoforms.** Every Exonerate CDS is classified as full-length, truncated or stop-lost (`extended`), separating nucleotide-level deletions (query segment absent, contig end, assembly gap) from ORF mutations (frameshift, premature stop, lost start/stop codon, non-canonical splice site). New outputs in `13_summary_data/`: `{query_id}_isoforms_incl_truncated.tsv` (all sequences with isoform IDs, status and mutations; samples without a hit listed as `no_hit`) and `{query_id}_isoform_mutations.tsv` (one row per mutation with query position, contig coordinates and codon change). `{query_id}_isoforms.tsv` now holds intact isoforms only.
- **Isoform sequences.** `13_summary_data/{query_id}_cds_isoforms.fasta` and `{query_id}_protein_isoforms.fasta` hold one sequence per unique isoform, with the number of samples and their status in the header. The isoform tables also have `cds_sequence` and `protein_sequence` columns.
- **Isoform trees.** IQ-TREE trees with one tip per CDS and per protein isoform (`9_phylogenies/exonerate/{query_id}/{query_id}_{cds,protein}_isoform_phylogeny.*`), built from `13_summary_data/{query_id}_{cds,protein}_isoforms_aligned.fasta`. Each tree is also drawn as a midpoint-rooted PDF, with tips labelled by isoform number and sample count, truncated isoforms coloured, the reference marked and bootstrap support shown. Requires at least 3 isoforms. Uses the `bioconductor-ggtree` container (R/ape).
- **TEs near the hit.** With `--repeat_library` and `--blast_query`, `13_summary_data/{query_id}_scaffold_TEs.tsv` lists every TE on the scaffold carrying each sample's hit (Exonerate gene model, or best BLAST hit), with class/family, strand relative to the hit, and edge-to-edge distance in the hit's 5'→3' orientation (negative = upstream, positive = downstream, 0 = overlap). `{query_id}_scaffold_TEs_summary.tsv` gives one row per hit with the overlapping and nearest upstream/downstream TEs. RepeatMasker `.out` tables are now kept in `12_repeatmasker/`.
- **Sample name collision check.** The run stops before assembly when two inputs (reads, SRA runs, scaffolds) end up with the same sample name, including after `--SRR_sample_map` renaming or after characters outside `[A-Za-z0-9._-]` are replaced by `_`. The error lists the clashing inputs. Before, their outputs overwrote each other and IQ-TREE failed on duplicated names.
- **`--exonerate_refine`** (`region` default, `full`, `none`) sets how Exonerate refines its alignment (see Fixed).

### Changed

- **README restructured.** One table of pipeline steps; options grouped by purpose in tables with their defaults and requirements; outputs listed per folder; one methods section instead of descriptions repeated under options and outputs. Corrected: profile CPUs and the default profile, the genome-wide depth statistic of `--locus_coverage`, BUSCO protein (not nucleotide) sequences for the species tree, the BUSCO download location, and outputs that were missing (fastp JSON reports, run statistics, query sequences, BWA summary text).
- **Isoforms are numbered by frequency.** Before, the query reference was isoform 1, and the others were numbered in alphabetical order of the samples. Now intact isoforms come first, ordered by the number of sample sequences carrying them, so the most frequent is isoform 1. They are followed by the query reference if no intact sample shares it, then by isoforms found only in truncated or stop-lost sequences, again by frequency. The intact-only table numbers consecutively. Numbering is shared by all isoform files.
- **Sequence names are sample IDs throughout.** Protein headers no longer carry a `_prot` suffix (they match the CDS: `sample_id`, or `sample_id_N` for several gene models; the translated nucleotide query is named by the first word of its header). In the combined BLAST hit FASTA, each hit is named by its sample ID with its location as description (`>sample_id scaffold:start-end`). MAFFT's `_R_` reverse-complement prefix is removed from all alignments and listed in the MAFFT log. Alignments and trees therefore show sample names.
- **Documentation of the per-query trees.** The README said the Exonerate CDS/protein trees require `--phylo_tree`. They are in fact always built with `--search_tool exonerate`, and `--phylo_tree` requires `--busco_lineage`, so the documented example failed. The README now describes the actual behaviour, and gives the RAxML-NG BUSCO model consistently as LG+G.

### Fixed

- **`--NCBI_API_key` had no effect.** The key was printed but never passed to the SRA lookup. It is now given to `entrez-direct` (`NCBI_API_KEY`), and the startup summary shows only whether a key is set instead of the key itself.
- **BLAST hit tree failed or was unreadable.** BLAST hits were named only by scaffold location, so tree tips did not show the sample, and IQ-TREE failed when two samples shared a scaffold name and location.
- **Exonerate with nucleotide CDS queries returned no hits.** Nucleotide queries now use the `coding2genome` model instead of `est2genome`. `est2genome` reported the aligned region as UTR/exon and emitted no CDS features, so nothing was extracted even when the exon structure was found. Separate alignments on the same target are no longer merged into one gene model, and only the best-scoring model is kept where alignments overlap.
- **Exonerate dropped exons of divergent alleles.** Exonerate now refines its heuristic alignment (`--refine region`). Before, intact but divergent alleles often came out as single-exon fragments and looked truncated.
- **Antisense CDS models.** For nucleotide queries, `coding2genome` alignments of the reverse-complemented query are no longer extracted. At divergent loci they could outscore the sense model and replace it.
- **Nucleotide BLAST missed divergent or short-exon hits.** `blastn` now runs with `-task blastn` instead of the megablast default.

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
