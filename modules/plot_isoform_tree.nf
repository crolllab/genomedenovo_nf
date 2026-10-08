// Render an IQ-TREE isoform phylogeny (cds_isoform / protein_isoform) as a PDF.
// The tree is midpoint-rooted; tips are labelled with the isoform number and the
// number of samples carrying it (from the isoform FASTA headers written by
// SummarizeIsoforms), truncated/extended isoforms are coloured and the query
// reference is marked. Uses the consensus tree (.contree, with bootstrap
// support) when IQ-TREE produced one, otherwise the ML tree (.treefile).
//
// Input:  [query_id, seq_type, IQ-TREE result files, ungapped isoform FASTA]
// Output: <query_id>_<seq_type>_phylogeny.pdf in 9_phylogenies/exonerate/<query_id>/
process PlotIsoformTree {
    tag "Plot isoform tree"
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    publishDir { "${params.outdir}/9_phylogenies/exonerate/${query_id}" }, mode: 'copy'

    input:
    tuple val(query_id), val(seq_type), path(tree_files, stageAs: 'tree/*'), path(isoform_fasta)

    output:
    tuple val(query_id), path("*_phylogeny.pdf"), emit: pdf

    script:
    """
    #!/bin/bash
    set -euo pipefail

    clean_id="\$(echo '${query_id}' | sed 's|[^a-zA-Z0-9_.-]|_|g')"
    prefix="\${clean_id}_${seq_type}_phylogeny"

    # Prefer the consensus tree with bootstrap support
    tree="tree/\${prefix}.contree"
    [ -s "\${tree}" ] || tree="tree/\${prefix}.treefile"

    plot_isoform_tree.R "\${tree}" '${isoform_fasta}' "\${prefix}.pdf" '${seq_type}'
    """
}
