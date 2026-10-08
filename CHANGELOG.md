# Changelog

## Unreleased

### Fixed

- **Exonerate dropped exons of divergent alleles.** Exonerate now refines its heuristic alignment (`--refine region`, configurable with `--exonerate_refine`). Before, intact but divergent alleles often came out as single-exon fragments and looked truncated.
- **Antisense CDS models.** For nucleotide queries, `coding2genome` alignments of the reverse-complemented query are no longer extracted. At divergent loci they could outscore the sense model and replace it.
- **Exonerate with nucleotide CDS queries returned no hits.** Nucleotide queries now use the `coding2genome` model instead of `est2genome`, which reported the aligned region as UTR/exon and emitted no CDS features, so nothing was extracted even when the exon structure was found. Separate alignments on the same target are no longer merged into one gene model, and only the best-scoring model is kept where alignments overlap.
- **Nucleotide BLAST missed divergent or short-exon hits.** `blastn` now runs with `-task blastn` instead of the megablast default.

### Changed

- **Isoforms are numbered by frequency.** The most frequent intact isoform is now isoform 1 (previously the query reference was 1, followed by isoforms in alphabetical order of the samples). Intact isoforms come first by the number of sequences carrying them, then the query reference if no intact sample shares it, then truncated/stop-lost-only isoforms by frequency.
- **No `_prot` suffix in protein headers.** Protein sequences now carry the same name as their CDS (`sample_id`, or `sample_id_N` for several gene models) in all protein FASTA files, alignments and trees; the translated nucleotide query is named by the first word of its header.
- **Documentation of the per-query trees.** The README stated that the Exonerate CDS/protein trees require `--phylo_tree`; they are in fact always built with `--search_tool exonerate` (and `--phylo_tree` requires `--busco_lineage`). The README now describes the actual behaviour, and the RAxML-NG BUSCO model is consistently given as LG+G.

### Added

- **TEs near the hit.** With `--repeat_library` and `--blast_query`, `13_summary_data/{query_id}_scaffold_TEs.tsv` lists every TE on the scaffold carrying each sample's hit (Exonerate gene model, or best BLAST hit), with class/family, strand relative to the hit, and edge-to-edge distance in the hit's 5'→3' orientation (negative = upstream, positive = downstream, 0 = overlap). `{query_id}_scaffold_TEs_summary.tsv` gives one row per hit with overlapping and nearest upstream/downstream TEs. RepeatMasker `.out` tables are now kept in `12_repeatmasker/`.
- **Isoform trees.** IQ-TREE trees with one tip per CDS and per protein isoform (`9_phylogenies/exonerate/{query_id}/{query_id}_{cds,protein}_isoform_phylogeny.*`), drawn as midpoint-rooted PDFs labelled with isoform number and sample count, with truncated isoforms coloured and bootstrap support shown.
- **Sequences in the isoform tables.** `{query_id}_isoforms.tsv` and `{query_id}_isoforms_incl_truncated.tsv` have `cds_sequence` and `protein_sequence` columns.
- **Truncated isoforms.** Every Exonerate CDS is classified as full-length or truncated, separating nucleotide-level deletions (query segment absent, contig end, assembly gap) from ORF mutations (frameshift, premature stop, lost start/stop codon, non-canonical splice site). `13_summary_data/{query_id}_isoforms_incl_truncated.tsv` lists all sequences with their isoform IDs, status and mutations, and `{query_id}_isoform_mutations.tsv` lists each mutation with query position, contig coordinates and codon change. `{query_id}_isoforms.tsv` now holds intact isoforms only; numbering is shared by all files, with truncated isoforms numbered after the intact ones. The query reference is isoform 1 for protein queries too.
- **Isoform sequences.** `13_summary_data/{query_id}_cds_isoforms.fasta` and `{query_id}_protein_isoforms.fasta` hold one sequence per unique isoform, numbered as in `{query_id}_isoforms.tsv`.

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
