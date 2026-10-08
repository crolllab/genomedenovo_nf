// Build a maximum-likelihood phylogeny from a single MAFFT alignment with IQ-tree3.
// Used for BLAST or Exonerate alignments and for the Exonerate isoform alignments
// (one aligned row per isoform, from SummarizeIsoforms) — one tree per query.
// Input:  [query_id, seq_type, alignment_fasta]
//         (seq_type: "blast", "cds", "protein", "cds_isoform" or "protein_isoform")
// Output: IQ-tree result files prefixed <query_id>_<seq_type>_phylogeny, published
//         in 9_phylogenies/<blast|exonerate>/<query_id>/.
process IQtreeAlignmentPhylogeny {
    tag "IQ-tree alignment phylogeny"
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6
    publishDir { "${params.outdir}/9_phylogenies/${seq_type == 'blast' ? 'blast' : 'exonerate'}/${query_id}" }, mode: 'copy'

    input:
    tuple val(query_id), val(seq_type), path(alignment_fasta)

    output:
    tuple val(query_id), path("*_phylogeny.*"), emit: results

    script:
    """
    set -euo pipefail

    clean_id="\$(echo '${query_id}' | sed 's|[^a-zA-Z0-9_.-]|_|g')"
    prefix="\${clean_id}_${seq_type}_phylogeny"

    n_taxa=\$(grep -c '^>' ${alignment_fasta} || echo 0)

    if [ "\$n_taxa" -lt 3 ]; then
        echo "Insufficient taxa (\$n_taxa) for phylogeny inference — need at least 3" > "\${prefix}.log"
        exit 0
    fi

    # Count unique sequences to avoid bootstrap failure when IQ-tree removes duplicates
    n_unique=\$(awk '/^>/{if(seq!="") seen[seq]++; seq=""; next}{seq=seq\$0} END{if(seq!="") seen[seq]++; print length(seen)}' ${alignment_fasta})

    if [ "\$n_unique" -lt 3 ]; then
        echo "Insufficient unique sequences (\$n_unique) for phylogeny inference — need at least 3" > "\${prefix}.log"
        exit 0
    fi

    if [ "\$n_unique" -ge 4 ]; then
        BOOT_ARG="-B 1000"
    else
        BOOT_ARG=""
    fi

    iqtree3 \\
        -s ${alignment_fasta} \\
        -m MFP \\
        \$BOOT_ARG \\
        -T ${task.cpus} \\
        --prefix "\${prefix}" \\
        --redo
    """
}

// Build a maximum-likelihood phylogeny from a BUSCO supermatrix with IQ-tree3.
// Input: supermatrix FASTA and NEXUS partition file from BuildSupermatrix.
// Output: IQ-tree result files including the treefile.
process IQtreePhylogeny {
    tag "IQ-tree BUSCO phylogeny"
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

    # Only apply ultrafast bootstrap when there are ≥4 unique taxa (IQ-tree requirement)
    n_taxa=\$(grep -c '^>' ${supermatrix})
    n_unique=\$(awk '/^>/{if(seq!="") seen[seq]++; seq=""; next}{seq=seq\$0} END{if(seq!="") seen[seq]++; print length(seen)}' ${supermatrix})
    if [ "\$n_unique" -ge 4 ]; then
        BOOT_ARG="-B 1000"
    else
        BOOT_ARG=""
    fi

    iqtree3 \\
        -s ${supermatrix} \\
        -p ${partitions} \\
        -m MFP \\
        \$BOOT_ARG \\
        -T ${task.cpus} \\
        --prefix busco_phylogeny \\
        --redo
    """
}
