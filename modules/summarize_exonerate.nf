// Process to merge all Exonerate query summaries into a single long-format table
process MergeExonerateSummaries {
    tag "Merging all Exonerate summaries"
    publishDir "${params.outdir}/7_exonerate_hits", mode: 'copy'

    input:
    path summary_files  // All *_all_exonerate_summary.tsv files

    output:
    path "exonerate_summary_merged.tsv", emit: merged_summary

    script:
    """
    #!/bin/bash
    set -euo pipefail

    tsv_files=( ${summary_files} )

    if [ \${#tsv_files[@]} -eq 0 ]; then
        echo "ERROR: No Exonerate summary files found"
        exit 1
    fi

    echo "Merging \${#tsv_files[@]} Exonerate summary files..."

    # Universal lowercase header (one row per sample × query)
    printf 'sample_id\tquery_id\tsubject_id\tstrand\tquery_start\tquery_end\tsubject_start\tsubject_end\texon_count\tcds_length\tprotein_length\texon_coordinates\n' \
        > exonerate_summary_merged.tsv

    # Append data rows from every per-query summary (skip their header lines)
    for tsv_file in "\${tsv_files[@]}"; do
        tail -n +2 "\$tsv_file" >> exonerate_summary_merged.tsv
    done

    echo "Merge complete: \$(tail -n +2 exonerate_summary_merged.tsv | wc -l) rows written."
    """
}
