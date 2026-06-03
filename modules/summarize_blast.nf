// Process to merge all BLAST query summaries into a single long-format table
process MergeBlastSummaries {
    tag "Merging all BLAST summaries"
    publishDir "${params.outdir}/6_blast_hits", mode: 'copy'

    input:
    path summary_files  // All *_all_blast_hits_summary.tsv files

    output:
    path "blast_summary_merged.tsv", emit: merged_summary

    script:
    """
    #!/bin/bash
    set -euo pipefail

    tsv_files=( ${summary_files} )

    if [ \${#tsv_files[@]} -eq 0 ]; then
        echo "ERROR: No BLAST summary files found"
        exit 1
    fi

    echo "Merging \${#tsv_files[@]} BLAST summary files..."

    # Universal lowercase header (one row per sample × query)
    printf 'sample_id\tquery_id\tsubject_id\tstrand\tpercent_identity\talignment_length\tquery_start\tquery_end\tsubject_start\tsubject_end\tevalue\tbitscore\n' \
        > blast_summary_merged.tsv

    # Append data rows from every per-query summary (skip their header lines)
    for tsv_file in "\${tsv_files[@]}"; do
        tail -n +2 "\$tsv_file" >> blast_summary_merged.tsv
    done

    echo "Merge complete: \$(tail -n +2 blast_summary_merged.tsv | wc -l) rows written."
    """
}
