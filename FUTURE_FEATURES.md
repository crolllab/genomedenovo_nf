# Future Features

- **BUSCO**: Minimize storage use by deleting intermediate files on the fly.

## Implemented (testing required)

- **RepeatMasker**: Repeat masking of assemblies using a user-provided consensus FASTA file (e.g., a custom repeat library).
- **Read mapping to reference genome**: To assess % mapped reads, coverage uniformity.
- **Read stats**: Get read stats, used fastp output statistics, creates a tsv file with read quality metrics.
- **IQ-TREE phylogenetics**: Maximum-likelihood phylogenetic tree from BUSCO single-copy ortholog supermatrix (MAFFT + IQ-tree2, activated with `--phylo_tree`).
- **IQ-TREE for each query sequence alignment**: Placed together with mafft alignments.
- **Replaced Genomethreader with Exonerate**
- **Provide read mapping statistics**: provide % mapped reads, coverage at focal interval (query hit on reference genome) versus genomic background
- **Isoform identification**


## To consider

- **Visualizations**: Quality stats, assembly metrics, TE coverage, phylogenetic tree, etc.
