// Concatenate per-gene MAFFT alignments into a supermatrix FASTA and
// NEXUS partition file for IQ-tree. Runs the Python helper script.
process BuildSupermatrix {
    tag "Build BUSCO supermatrix"
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3

    input:
    path aligned_genes
    path concat_script

    output:
    path "supermatrix.fasta", emit: supermatrix
    path "partitions.nex",    emit: partitions

    script:
    """
    set -euo pipefail
    python3 ${concat_script}
    """
}
