// Transposable elements on the scaffold carrying the primary hit of a query.
// For every sample, lists the RepeatMasker TEs on the hit scaffold with their
// position and edge-to-edge distance relative to the hit, in the hit's own
// 5'->3' orientation (negative = upstream, positive = downstream, 0 = overlap).
// Hit coordinates come from the Exonerate summary (--search_tool exonerate) or
// the BLAST summary (merged HSP span on the best scaffold).
//
// Input:
//   tuple val(query_id), path(hits_tsv), path(rm_out_files)
//         hits_tsv: *_all_exonerate_summary.tsv or *_all_blast_hits_summary.tsv
//         rm_out_files: RepeatMasker .out tables of all samples
//
// Output (published to 13_summary_data/):
//   path "*_scaffold_TEs.tsv"           one row per TE on the hit scaffold
//   path "*_scaffold_TEs_summary.tsv"   one row per hit (overlapping / nearest TEs)
process SummarizeHitTEs {
    tag "TEs near hits"
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    publishDir "${params.outdir}/13_summary_data", mode: 'copy'

    input:
    tuple val(query_id), path(hits_tsv), path(rm_out_files, stageAs: 'rm/*')

    output:
    tuple val(query_id), path("*_scaffold_TEs.tsv"),         emit: tes
    tuple val(query_id), path("*_scaffold_TEs_summary.tsv"), emit: summary

    script:
    """
    #!/bin/bash
    set -euo pipefail
    shopt -s nullglob
    rm_out=( rm/*_repeatmasker.out )
    summarize_hit_tes.py '${query_id}' '${hits_tsv}' --rm-out "\${rm_out[@]}"
    """
}
