// Build a maximum-likelihood phylogeny from a BUSCO supermatrix with IQ-tree2.
// Input: supermatrix FASTA and NEXUS partition file from BuildSupermatrix.
// Output: IQ-tree result files including the treefile.
process IQtreePhylogeny {
    tag "IQ-tree BUSCO phylogeny"
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    publishDir "${params.outdir}/phylogeny", mode: 'copy'

    input:
    path supermatrix
    path partitions

    output:
    path "busco_phylogeny.*", emit: results

    script:
    """
    set -euo pipefail

    # Only apply ultrafast bootstrap when there are ≥4 taxa (IQ-tree requirement)
    n_taxa=\$(grep -c '^>' ${supermatrix})
    if [ "\$n_taxa" -ge 4 ]; then
        BOOT_ARG="-B 1000"
    else
        BOOT_ARG=""
    fi

    iqtree2 \\
        -s ${supermatrix} \\
        -p ${partitions} \\
        -m MFP \\
        \$BOOT_ARG \\
        -T ${task.cpus} \
        --prefix busco_phylogeny \\
        --redo
    """
}
