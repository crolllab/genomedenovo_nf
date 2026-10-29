# Methods

This page describes what each step of the pipeline does. The options that control the steps are listed on the [Configuration](configuration.md) page, their output files on the [Output](output.md) page.

---

## Reads and assembly

**SRA download.** Accessions are resolved to runs with `entrez-direct`. Each run is downloaded from ENA when available, with `sra-tools` from NCBI as fallback.

**Trimming and assembly.** Reads are trimmed with `fastp` (default settings) and assembled per sample with `SPAdes --isolate`.

**Assembly QC.** QUAST reports contiguity metrics. BUSCO v6 (genome mode, miniprot) scores completeness against the chosen lineage.

**BUSCO species tree.** The single-copy BUSCO proteins found in at least 2 samples are subsampled (`--busco_subsample`), aligned per gene with MAFFT and concatenated into a partitioned supermatrix. RAxML-NG uses an LG+G model per gene and 200 bootstraps. IQ-TREE selects a model per partition (`-m MFP`) with 1000 ultrafast bootstraps. Bootstrapping needs at least 4 taxa (4 unique sequences for IQ-TREE).

---

## Query search

**BLAST search.** Each query is searched against every assembly (`blastn -task blastn` or `tblastn`). Per sample, the scaffold with the best bit score is kept, and all hits of the query on that scaffold are merged into one span. That way, the separate hits of a multi-exon gene make up one region. The span ± `--padding` is extracted in the orientation of the query.

**Exonerate CDS.** Exonerate aligns the query to the extracted region (`protein2genome` for protein queries, `coding2genome` for CDS queries) to find the exon boundaries. The aligned exons give the CDS, which is translated to protein. Where gene models overlap, the best-scoring one is kept, and antisense models of nucleotide queries are discarded.

---

## CDS integrity

Each gene model is compared codon by codon with the query to find frameshifts, premature stops, lost start/stop codons, non-canonical splice sites and in-frame indels. Query segments that Exonerate left unaligned are looked up on the full contig:

- An in-frame start or stop codon at the expected distance means a divergent but intact end.
- If the contig ends before the missing segment (or the stop codon) could fit, the gene runs off the contig: a contig edge. The model is then `partial`, not `truncated`, unless it also carries a deletion or ORF mutation. An assembly gap next to the missing segment is treated alike.
- Otherwise the segment is re-aligned to the flanking sequence. If it is found, the sequence is present but outside the reading frame: an ORF mutation. If it is not found, it is a nucleotide-level deletion. Telomeric repeats in the flank and deletions close to a contig end are noted.

Separately, the model is translated from its first codon to find where the reading frame actually ends (`orf_end`): at the expected stop codon (within 3 codons, more for divergent termini), earlier (cut short), later, or not before the contig end.

---

## Isoforms

Identical aligned CDS (or protein) sequences form one isoform. Numbering follows frequency:

1. Isoforms carried by at least one intact sequence, ordered by the number of sample sequences carrying them; the most frequent is 1. Ties keep alignment order, and the query reference is not counted.
2. The query reference, if no intact sample shares its sequence.
3. Isoforms found only in truncated, stop-lost or partial sequences, again by frequency.

The intact-only table therefore numbers consecutively.

---

## Alignments and per-query trees

BLAST hits, Exonerate CDS and Exonerate proteins are aligned per query with MAFFT (`--auto --adjustdirection`). For protein queries, the query is left out of the BLAST hit alignment, since the hits are nucleotide. IQ-TREE builds a tree from every alignment with ≥ 3 unique sequences (`-m MFP`; 1000 ultrafast bootstraps with ≥ 4 unique sequences). Isoform trees use one aligned row per isoform and are drawn as PDFs with R/`ape`.

### Alignment filter

Before a tree is built, sequences that align only partially are left out: a sequence must have residues in at least `--min_align_prop` (default 0.8) of the core alignment columns. The core columns are those where the query reference has a residue; alignments without the reference (BLAST hits of protein queries, isoform alignments) use the columns where at least half of the sequences have a residue.

Measuring against the core rather than all columns keeps a single long outlier (a TE insertion, or a hit spanning two loci) from pushing every other sequence below the threshold; the hit flanks added by `--padding` do not count either. The query reference is always kept. The decision for every sequence is listed in `.alignment_filter.tsv` next to the tree (`core_basis` says which core was used).

---

## Paralog check

With `--query_paralog`, the paralogs of each query are BLASTed against every assembly with the same program and settings as the query. This runs before the hit extraction, per sample:

1. Every BLAST hit (HSP) of the query and the paralogs is assigned to its reference, unless a hit of another reference overlaps it on the scaffold with a higher bit score; then it belongs to that reference. Ties go to the query.
2. A reference's own hits on one scaffold form a copy; the copy is split where a hit of another reference lies in between (tandem arrays).
3. The query copy to extract is chosen as without paralogs (the one with the best single hit), but among the query's own copies only. It must also win at the locus level: bit scores of non-overlapping hits over the copy's span are summed (greedy by score) per reference, and no paralog may score higher than the query. Otherwise the copy is skipped and the next is tried.
4. Only the hits of the chosen copy go on to the hit extraction. If no query copy remains, the sample is `paralog_only`: nothing is extracted, and it is absent from the hit FASTA, alignments, trees, Exonerate, CDS integrity, isoforms and the TE summary.

Paralog copies on the scaffold of the extracted query copy (e.g. a tandem paralog) are kept out of its flanks: the `--padding` stops short of them, so Exonerate can neither model the paralog nor fuse it to a truncated query copy, and the CDS integrity check masks them when it looks for missing query segments in the flanking contig. A query exon that scores higher against a paralog than against the query is treated as the paralog's, consistent with step 1. As a further guard, Exonerate gene models must overlap the extracted BLAST hit span; a model that lies only in the flank is discarded.

With `--locus_coverage`, the query's locus on the reference is dropped if a paralog scores higher there. Ambiguous calls are kept downstream and flagged in the paralog table. See [the output tables](output.md#paralog-check-query_paralog).

---

## Reference mapping and locus coverage

Trimmed reads are mapped with BWA-mem2, filtered at MAPQ ≥ 10, sorted and indexed. For `--locus_coverage`, the query is BLASTed against the reference to find its locus. Per sample, the mean and median depth over the locus ± `--padding` (MAPQ ≥ 20) are compared with the mean and SD of per-base depth over the whole genome.

Depth and breadth over the hit itself are checked against the sample's assembly: a sample without the query in its assembly (no BLAST hit, or only paralog hits with `--query_paralog`) but with reads covering the locus is flagged in `coverage_conflict`. The `paralog_only` / `no_hit` call then reflects the assembly, not the genome.

---

## Transposable elements

RepeatMasker annotates every assembly with the given library. TEs on the scaffold of each query hit are reported with their position relative to the hit (see [TEs near the hit](output.md#tes-near-the-hit-repeat_library-with-blast_query)); simple repeats and low-complexity regions are left out.
