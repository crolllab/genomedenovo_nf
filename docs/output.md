# Pipeline output

All output files are written to the directory specified by `--outdir` (default: `./nf_output`).

Output folders are **numbered in pipeline order**. A folder is only created when the step that
fills it actually runs, so gaps in the numbering are normal — a `--scaffolds`-only run has no
`1_fastq_stats/` or `2_assemblies/`, one without `--busco_lineage` has no `4_busco/`.

`{query_id}` is the first word of the query FASTA header, with characters outside `[A-Za-z0-9._-]` replaced by `_` in file names.

```
<outdir>/
├── 1_fastq_stats/       fastp reports + summary table              read input
├── 2_assemblies/        SPAdes assemblies per sample               read input
├── 3_quast/             Assembly metrics                           unless --skip_quast
├── 4_busco/             Genome completeness                        --busco_lineage
├── 5_busco_phylogeny/   BUSCO species tree                         --phylo_tree
├── 6_blast_hits/        BLAST results and extracted hits           --blast_query
├── 7_exonerate_hits/    Exonerate CDS and protein                  --search_tool exonerate
├── 8_alignments/        MAFFT alignments per query                 --blast_query
├── 9_phylogenies/       IQ-TREE trees per query and isoform        --blast_query
├── 10_bwa_mapping/      Reads mapped to the reference              --reference_genome
├── 11_locus_coverage/   Read depth at the query locus              --locus_coverage
├── 12_repeatmasker/     TE annotation per assembly                 --repeat_library
└── 13_summary_data/     Isoforms, CDS integrity, paralogs, TEs     depends on options
```

---

## Run overview (`<outdir>/`)

| File | Description |
|------|-------------|
| `pipeline_execution_stats.tsv` / `.txt` | Run time and resources per step. |
| `pipeline_trace.txt` | Nextflow trace of every task, including `%cpu` and `%mem`. |
| `sample_names.tsv` | Every input (`input_id`), its `--SRR_sample_map` name (`renamed_to`) and the sample name used in all outputs (`sample_id`). |
| `NCBI_SRR_PE_accessions.txt`, `NCBI_SRR_SE_accessions.txt` | Resolved paired-end and single-end SRA runs (`--SRA_index`). |
| `NCBI_download_urls.tsv` | Download sources of every resolved run (`--SRA_index`). |
| `busco_db_cache/` | Downloaded BUSCO lineage (`--busco_lineage` without `--busco_db_path`). |

---

## `1_fastq_stats/`

Read-trimming QC. Not generated for `--scaffolds`-only runs.

| File | Description |
|------|-------------|
| `fastp_summary.tsv` / `.txt` | Read metrics before and after trimming for all samples. |
| `{sample_id}_{PE,SE}_fastp.json` | fastp report per sample. |

---

## `2_assemblies/`

| File | Description |
|------|-------------|
| `{sample_id}/scaffolds.fasta` | SPAdes scaffolds; the assembly used by all downstream steps. |
| `{sample_id}/contigs.fasta` | SPAdes contigs. |
| `{sample_id}/spades.log` | SPAdes log. |

---

## `3_quast/`

| File | Description |
|------|-------------|
| `assembly_summary_metrics.tsv` | QUAST metrics for all samples. |
| `{sample_id}_quast/` | Full QUAST report per sample. |

---

## `4_busco/`

| File | Description |
|------|-------------|
| `busco_summary_metrics.tsv` | Complete, single-copy, duplicated, fragmented and missing % for all samples. |
| `short_summary.*.{sample_id}_busco.txt` | BUSCO short summary per sample. |
| `{sample_id}_busco/` | Per-sample results: single-copy protein sequences, or the complete BUSCO folder with `--keep_full_busco`. |

---

## `5_busco_phylogeny/`

| Tool | Files |
|------|-------|
| RAxML-NG (default) | `busco_phylogeny.raxml.support` (best tree with bootstrap support), `.bestTree`, `.bestModel`, `.bootstraps`, `.mlTrees`, `.log` |
| IQ-TREE (`--busco_tree iqtree`) | `busco_phylogeny.contree` (tree with bootstrap support), `.treefile`, `.iqtree` (report), `.best_scheme.nex`, `.log` |

---

## `6_blast_hits/`

| File | Description |
|------|-------------|
| **`blast_summary_merged.tsv`** | One row per sample; per query `{query}_Present`, `_Identity`, `_AlignLen`, `_Evalue`, `_BitScore`. |
| `query_sequences/` | The split queries and their detected type. |
| `{query_id}/{query_id}_all_blast_hits.fasta` | Query followed by the hit of every sample, named `>sample_id scaffold:start-end`. |
| `{query_id}/{query_id}_all_blast_hits_summary.tsv` | Hit coordinates, identity, e-value and bit score per sample (`no_hits_found` without a hit). |
| `{query_id}/blast_results/{sample_id}_blast_results.tsv` | Raw BLAST tabular output. |
| `{query_id}/hit_sequences/{sample_id}_top_hits.fasta` / `_top_hits_info.tsv` | Extracted (padded) hit region per sample, named by its scaffold location (`>scaffold:start-end hit=start-end`, the unpadded hit span; with `--query_paralog` also `paralogs=start-end,...`, the paralog copies on the scaffold that the padding stops short of). |
| `{query_id}/paralog_check/{sample_id}_paralog_blast_results.tsv` | Raw BLAST output of the query's paralogs (`--query_paralog`). |

!!! note
    With `--query_paralog`, samples whose query hits all belong to a paralog have `subject_id`
    `paralog_only` in the hit summaries (BLAST and Exonerate), and no hit sequence.
    `blast_results/` keeps the raw, unfiltered BLAST output.

---

## `7_exonerate_hits/`

Only generated with `--search_tool exonerate`.

| File | Description |
|------|-------------|
| `exonerate_summary_merged.tsv` | Exonerate results for all queries and samples. |
| `{query_id}/{query_id}_all_exonerate_cds.fasta` / `_protein.fasta` | CDS and protein of every sample (and the query), named by sample ID (`{sample}_{n}` when a sample has several gene models). |
| `{query_id}/{query_id}_all_exonerate_summary.tsv` | Scaffold, strand, coordinates, exon count and lengths per gene model. |
| `{query_id}/hit_sequences/` | CDS, protein and summary per sample. |
| `{query_id}/gff3_results/` | Raw Exonerate GFF per sample. |

---

## `8_alignments/`

| File | Description |
|------|-------------|
| `blast/{query_id}/{query_id}_all_blast_hits_aligned_sequences.fasta` | MAFFT alignment of the BLAST hits. |
| `exonerate/{query_id}/{query_id}_all_exonerate_{cds,protein}_aligned_sequences.fasta` | CDS and protein alignments (`--search_tool exonerate`). |
| `*_mafft_alignment_log.txt` | MAFFT log, including the sequences MAFFT reverse-complemented. |

---

## `9_phylogenies/`

IQ-TREE trees, built when an alignment has at least 3 unique sequences after [partially aligned sequences are left out](methods.md#alignment-filter). Each tree comes with:

| File | Description |
|------|-------------|
| `.treefile` | Maximum-likelihood tree. |
| `.contree` | Tree with bootstrap support (≥ 4 unique sequences). |
| `.iqtree` / `.log` | IQ-TREE report and log. |
| `.alignment_filter.tsv` | Per sequence: residues in the core columns, number of core columns and their basis, proportion, kept or not. |
| `.filtered_alignment.fasta` | The tree input. |

| Tree | Files |
|------|-------|
| BLAST hits, one tip per sample | `blast/{query_id}/{query_id}_blast_phylogeny.*` |
| Exonerate CDS / protein, one tip per sample | `exonerate/{query_id}/{query_id}_{cds,protein}_phylogeny.*` |
| CDS / protein isoforms, one tip per isoform | `exonerate/{query_id}/{query_id}_{cds,protein}_isoform_phylogeny.*` |
| **Isoform tree plots** | `exonerate/{query_id}/{query_id}_{cds,protein}_isoform_phylogeny.pdf` |

The isoform tree PDFs are drawn midpoint-rooted. Tips read `Isoform N (n = samples)`; isoforms found only in truncated (orange), stop-lost (blue) or partial (grey, contig edge) sequences are coloured, the query reference is marked with a triangle, and bootstrap support is shown at the nodes.

---

## `10_bwa_mapping/`

Only generated with `--reference_genome` and read input.

| File | Description |
|------|-------------|
| `bwa_mapping_summary.tsv` / `.txt` | Mapping statistics for all samples. |
| `{sample_id}_sorted.bam` (+ `.bai`) | Sorted alignments, MAPQ ≥ 10. |
| `{sample_id}_flagstat.json` | SAMtools flagstat per sample. |

---

## `11_locus_coverage/`

Only generated with `--locus_coverage`. One folder per query.

| File | Description |
|------|-------------|
| `{query_id}/{query_id}_all_locus_coverage.tsv` | All samples (columns below); per-sample tables `{sample_id}_{query_id}_locus_coverage.tsv`. Samples are named as in the other tables (after `--SRR_sample_map`). |
| `{query_id}/{query_id}_ref_locus.tsv` | Best hit of the query on the reference. Empty with `--query_paralog` when a paralog scores higher at that hit: the reference has no copy of the query, and the coverage columns are `NA`. |

### Locus coverage columns

| Column | Description |
|--------|-------------|
| `scaffold`, `hit_start`, `hit_end` | Query hit on the reference. |
| `locus_start`, `locus_end` | Hit ± `--padding`. |
| `locus_mean_depth`, `locus_median_depth` | Read depth over the padded locus. |
| `genome_mean_depth`, `genome_sd_depth` | Per-base depth over the whole genome. |
| `hit_mean_depth`, `hit_breadth_pct` | Depth and breadth (% of bases with ≥ 1 read) over the unpadded hit. |
| `hit_to_genome_depth_ratio` | `hit_mean_depth` / `genome_mean_depth`. |
| `assembly_hit` | The sample's query hit in its assembly: `hit`, `no_hit`, `paralog_only`. |
| `coverage_conflict` | `no_hit_locus_covered` / `paralog_only_locus_covered` when the assembly has no copy of the query but the reads cover its reference locus (thresholds `--locus_conflict_min_depth_ratio`, `--locus_conflict_min_breadth`), e.g. an allele that failed to assemble; otherwise `none`. |

!!! warning "Coverage conflicts"
    Samples with a `coverage_conflict` are also listed as warnings in the run log. Their
    `no_hit` / `paralog_only` call reflects the assembly, not the genome.

---

## `12_repeatmasker/`

Only generated with `--repeat_library`.

| File | Description |
|------|-------------|
| `{sample_id}_repeatmasker.out` | Annotation table with TE class/family. |
| `{sample_id}_repeatmasker.gff` | Annotations as GFF (no class/family). |
| `{sample_id}_repeatmasker.masked`, `.tbl` | Masked assembly and summary. |

---

## `13_summary_data/`

### TEs near the hit (`--repeat_library` with `--blast_query`)

The hit is the Exonerate gene model with `--search_tool exonerate`, otherwise the BLAST hit. Distances are edge to edge in bp, in the hit's own **5'→3' orientation**: negative = TE upstream (5') of the hit, positive = downstream (3'), 0 = TE overlaps the hit (introns included). On a minus-strand hit, upstream therefore means higher scaffold coordinates.

| File | Description |
|------|-------------|
| **`{query_id}_scaffold_TEs.tsv`** | One row per TE on the hit scaffold, ordered 5'→3' (columns below). |
| `{query_id}_scaffold_TEs_summary.tsv` | One row per hit: number of TEs on the scaffold, overlapping TEs, nearest upstream and downstream TE with distance; `hit_status` (`hit`, `no_hit`, `paralog_only`). |

`{query_id}_scaffold_TEs.tsv` columns: `sample_id`, `hit_id` (`{sample}_{n}` for several gene models), `query_id`, `scaffold`, `hit_start`, `hit_end`, `hit_strand`, `te_name`, `te_class`, `te_start`, `te_end`, `te_strand`, `te_orientation` (`same`/`opposite` strand as the hit), `te_length`, `position` (`upstream`, `downstream`, `overlap`; `none` = no TE on the scaffold, `no_hit` = sample without hit, `paralog_only` = the query hit only paralog loci, `--query_paralog`), `distance_bp`, `perc_div` (divergence from the consensus), `sw_score`, `te_id`.

### Paralog check (`--query_paralog`)

Paralog hits never enter the other outputs: a sample's query hit is extracted only where it belongs to the query (see [Methods](methods.md#paralog-check)).

| File | Description |
|------|-------------|
| **`{query_id}_paralog_check.tsv`** | One row per sample (columns below). |
| **`{query_id}_paralog_loci.tsv`** | Every copy of the query and of each paralog in every sample: `sample_id`, `query_id`, `reference`, `reference_label` (`query`, `paralog1`, …), `scaffold`, `start`, `end`, `n_hsps`, `bitscore`, `identity_pct`, `coverage_pct`, `extracted` (`yes` for the query copy used downstream). |

`{query_id}_paralog_check.tsv` columns:

| Column | Description |
|--------|-------------|
| `sample_id`, `query_id` | |
| `call` | `query` (the query's own copy was found and extracted), `ambiguous` (as `query`, but the second-best reference scores ≥ `--paralog_ambiguous_ratio` × the query there), `paralog` (the query hits only paralog loci; nothing is extracted, and the sample is `paralog_only` in all other tables), `no_hit`. |
| `hit_scaffold`, `hit_start`, `hit_end` | The extracted query copy; for `paralog`, the paralog locus that the query hit best. |
| `closest_reference`, `best_bitscore`, `second_reference`, `second_bitscore`, `score_ratio` | Ranking of the query and all paralogs at that locus. |
| `query_` / `paralogN_` `bitscore_at_hit`, `identity_pct_at_hit`, `coverage_pct_at_hit` | Score of every reference at that locus. |
| `skipped_paralog_loci` | Query hit loci that were not extracted because a paralog scores higher there (`scaffold:start-end(paralogN)`), e.g. when the best query hit lies on a paralog but a weaker query copy exists elsewhere. |
| `query_` / `paralogN_` `n_copies`, `best_copy`, `best_copy_bitscore`, `best_copy_identity_pct`, `best_copy_coverage_pct` | Every reference's own copies in the sample (number, and the best one with its location and score). |

### Isoforms and CDS integrity (`--search_tool exonerate`)

Isoform numbers are shared by all files of a query (see [Isoform numbering](methods.md#isoforms)).

| File | Description |
|------|-------------|
| `{query_id}_isoforms.tsv` | The query reference and every intact sequence: `sample_id`, `cds_isoform`, `protein_isoform`, `cds_sequence`, `protein_sequence`. |
| **`{query_id}_isoforms_incl_truncated.tsv`** | Every sequence, including truncated ones and samples without a hit (columns below). |
| **`{query_id}_isoform_mutations.tsv`** | One row per mutation with the isoform IDs of the sequence carrying it, event class (`ORF_mutation`, `nucleotide_deletion`, `contig_edge`, `assembly_gap`, `in_frame_indel`), type, query residues, contig coordinates, reference and observed codon/residue, indel size and description. |
| `{query_id}_{cds,protein}_isoforms.fasta` | One sequence per isoform, named `cds_isoform_N` / `protein_isoform_N`, with the number of samples (`n_samples=`) and their status (`status=`) in the header. |
| `{query_id}_{cds,protein}_isoforms_aligned.fasta` | The isoforms as aligned rows; input of the isoform trees. |

`{query_id}_isoforms_incl_truncated.tsv` columns:

| Column | Description |
|--------|-------------|
| `sample_id`, `cds_isoform`, `protein_isoform` | |
| `status` | `reference`, `full_length`, `truncated`, `extended` (stop codon lost only), `partial` (incomplete only because the contig ends or an assembly gap), `no_hit`, `paralog_only` (the query hit only paralog loci, `--query_paralog`). |
| `truncation_cause` | `contig_edge` (the gene runs off the contig), `assembly_gap` (N run next to the missing part), `nucleotide_deletion` (part of the query is missing although the contig continues), `ORF_mutation` (the sequence is there but the reading frame is disrupted), several joined by `+`, or `none`. The first two are limits of the assembly, the last two evidence of a truncated gene. |
| `contig_edge` | Contig end(s) the gene model or its reading frame runs into: `5prime`, `3prime`, `5prime,3prime` or `none`. |
| `orf_end` | How the encoded reading frame ends: `expected_stop`, `premature_stop` (ORF cut short), `late_stop` (read-through to a later stop), `contig_end` (no stop codon before the contig ends; `observed_protein_length` is then a minimum, `>=N`). |
| `deletion_types` | `5prime_deletion` / `3prime_deletion` (segment missing, contig continues with unrelated sequence), `internal_deletion` (in-frame deletion ≥ 10 aa). |
| `orf_mutation_types` | `frameshift`, `premature_stop`, `start_lost`, `stop_lost`, `noncanonical_splice_site`, `segment_outside_ORF` (the unaligned query segment is found on the contig outside the reading frame, e.g. an unbridged frameshift or a retained intron). |
| `contig`, `contig_length`, `strand`, `gene_start`, `gene_end`, `exon_count` | Gene model location. |
| `contig_flank_5prime_nt`, `contig_flank_3prime_nt` | Contig sequence beyond the first / last coding base, in gene orientation (small values mean the gene sits at a contig end). |
| `aligned_query_start_aa`, `aligned_query_end_aa`, `query_coverage_pct`, `aa_identity_pct` | Alignment extent and similarity. |
| `observed_protein_length` | Protein actually encoded (first aligned codon to first stop). |
| `in_frame_query_aa`, `in_frame_query_pct` | Query residues encoded in frame before that stop. |
| `start_codon`, `stop_codon`, `n_frameshifts`, `n_premature_stops`, `n_noncanonical_splice_sites`, `n_inframe_indels`, `first_disruption_aa` | Codon-level summary. |
| `mutation_details` | Readable list of mutations (e.g. `p.Q46* nonsense codon TAG`, `p.C36fs frameshift -1 nt (deletion in contig)`). |
| `notes` | E.g. termini Exonerate left unaligned although the ORF is intact, telomeric repeats or assembly gaps in the flank, a deletion called close to a contig end (`3prime_contig_end_Nnt_beyond_alignment`), low-confidence frameshift calls. |
| `cds_sequence`, `protein_sequence` | |
