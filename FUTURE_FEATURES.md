# Future Features

## Implemented (testing required)

- **RepeatMasker**: Repeat masking of assemblies using a user-provided consensus FASTA file (e.g., a custom repeat library).
- **Read mapping to reference genome**: To assess % mapped reads, coverage uniformity.
- **Read stats**: Get read stats, used fastp output statistics, creates a tsv file with read quality metrics.
- **IQ-TREE phylogenetics**: Maximum-likelihood phylogenetic tree from BUSCO single-copy ortholog supermatrix (MAFFT + IQ-tree2, activated with `--phylo_tree`).

## Planned

- **Replace or improve on the Genomethreader step**
- **Metadata grouping**: Allow users to specify groups of samples (e.g., by population, treatment) for comparative analyses and visualizations.
- **Visualizations**: Quality stats, assembly metrics, TE landscapes, etc.
- **Provide read mapping statistics**: provide % mapped reads, coverage at focal interval versus genomic background
- **TE analyses**: Consider summarizing content, flipping focal contig to match query orientation, etc.