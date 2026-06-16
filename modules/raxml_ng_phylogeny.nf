// Build a maximum-likelihood phylogeny from a BUSCO supermatrix with RAxML-NG.
// Input: supermatrix FASTA and NEXUS partition file from BuildSupermatrix.
// Output: RAxML-NG result files including the best tree and support values.
process RaxmlNgPhylogeny {
    tag "RAxML-NG BUSCO phylogeny"
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6
    publishDir "${params.outdir}/5_busco_phylogeny", mode: 'copy'

    input:
    path supermatrix
    path partitions

    output:
    path "busco_phylogeny.*", emit: results

    script:
    """
    set -euo pipefail

    n_taxa=\$(grep -c '^>' ${supermatrix})
    if [ "\$n_taxa" -lt 3 ]; then
        echo "Insufficient taxa (\$n_taxa) for phylogeny inference — need at least 3" > busco_phylogeny.log
        exit 0
    fi

    # Convert NEXUS charset block to RAxML-NG partition format (LG+G per gene — protein model)
    awk '/charset/{gsub(/[;]/,""); split(\$0,a," = "); locus=a[1]; sub(/^ +charset +/,"",locus); print "LG+G, " locus " = " a[2]}' \\
        ${partitions} > partitions.raxml

    if [ "\$n_taxa" -ge 4 ]; then
        BS_ARG="--bs-trees 200"
    else
        BS_ARG=""
    fi

    raxml-ng \\
        --all \\
        --msa ${supermatrix} \\
        --model partitions.raxml \\
        \$BS_ARG \\
        --tree pars{10} \\
        --threads ${task.cpus} \\
        --prefix busco_phylogeny \\
        --redo
    """
}
