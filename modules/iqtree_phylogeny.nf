// Build a maximum-likelihood phylogeny from a single MAFFT alignment with IQ-tree3.
// Used for BLAST or Exonerate alignments and for the Exonerate isoform alignments
// (one aligned row per isoform, from SummarizeIsoforms) — one tree per query.
// Sequences covering less than --min_align_prop of the core alignment are left
// out of the tree; the query reference is always kept. Core columns are the
// columns of the query reference, or, in alignments without it (BLAST hits of
// protein queries, isoforms), the columns with a residue in at least half of
// the sequences.
// Input:  [query_id, seq_type, alignment_fasta]
//         (seq_type: "blast", "cds", "protein", "cds_isoform" or "protein_isoform")
// Output: IQ-tree result files prefixed <query_id>_<seq_type>_phylogeny, published
//         in 9_phylogenies/<blast|exonerate>/<query_id>/, plus
//         <prefix>.alignment_filter.tsv (coverage of every sequence, kept or not)
//         and <prefix>.filtered_alignment.fasta (the tree input).
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
    // The query reference is named by the first word of the query header
    def ref_name = query_id.trim().tokenize(' \t')[0]
    """
    set -euo pipefail

    clean_id="\$(echo '${query_id}' | sed 's|[^a-zA-Z0-9_.-]|_|g')"
    prefix="\${clean_id}_${seq_type}_phylogeny"
    aln="\${prefix}.filtered_alignment.fasta"

    # Leave out partially aligned sequences: the proportion of the core
    # alignment covered by a sequence must reach --min_align_prop. Core columns
    # are those of the query reference or, without it, those where at least
    # half of the sequences have a residue. Core rather than all columns, so
    # that one long outlier (e.g. a TE insertion or two merged loci) cannot
    # push every other sequence below the threshold.
    awk -v minp=${params.min_align_prop} -v ref='${ref_name}' -v report="\${prefix}.alignment_filter.tsv" '
        /^>/ { n++; hdr[n] = \$0; name[n] = substr(\$1, 2); seq[n] = ""; if (name[n] == ref) r = n; next }
        { gsub(/[ \\t\\r]/, ""); seq[n] = seq[n] \$0 }
        END {
            for (i = 1; i <= n; i++) {
                s = seq[i]; len = length(s)
                for (j = 1; j <= len; j++) if (substr(s, j, 1) != "-") occ[j]++
                if (len > L) L = len
            }
            core = 0
            basis = r ? "query_reference" : "half_of_sequences"
            for (j = 1; j <= L; j++) {
                if (r ? (j <= length(seq[r]) && substr(seq[r], j, 1) != "-") : (2 * occ[j] >= n)) { core++; is_core[j] = 1 }
            }
            print "sequence\\tresidues_in_core_columns\\tcore_columns\\tcore_basis\\tproportion\\tkept" > report
            for (i = 1; i <= n; i++) {
                s = seq[i]; len = length(s); c = 0
                for (j = 1; j <= len; j++) if (is_core[j] && substr(s, j, 1) != "-") c++
                prop = core ? c / core : 0
                if (name[i] == ref)     kept = "yes (query reference)"
                else if (prop >= minp)  kept = "yes"
                else                    kept = "no"
                printf "%s\\t%d\\t%d\\t%s\\t%.3f\\t%s\\n", name[i], c, core, basis, prop, kept > report
                if (kept != "no") { print hdr[i]; print s }
            }
        }' ${alignment_fasta} > "\${aln}"

    n_all=\$(grep -c '^>' ${alignment_fasta} || true)
    n_taxa=\$(grep -c '^>' "\${aln}" || true)
    n_dropped=\$(( n_all - n_taxa ))

    if [ "\$n_taxa" -lt 3 ]; then
        echo "Insufficient taxa (\$n_taxa, after leaving out \$n_dropped sequence(s) below --min_align_prop ${params.min_align_prop}) for phylogeny inference — need at least 3" > "\${prefix}.log"
        exit 0
    fi

    # Count unique sequences to avoid bootstrap failure when IQ-tree removes duplicates
    n_unique=\$(awk '/^>/{if(seq!="") seen[seq]++; seq=""; next}{seq=seq\$0} END{if(seq!="") seen[seq]++; print length(seen)}' "\${aln}")

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
        -s "\${aln}" \\
        -m MFP \\
        \$BOOT_ARG \\
        -T AUTO --threads-max ${task.cpus} \\
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
