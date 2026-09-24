// Assign consecutive isoform IDs to unique CDS and protein haplotypes
// derived from per-query MAFFT-aligned FASTA files, and flag truncated
// isoforms using the per-sample CDS integrity tables (AssessCdsIntegrity).
//
// Input:
//   tuple val(query_id), path(cds_alignment), path(prot_alignment),
//         path(integrity_tsvs), path(integrity_event_tsvs)
//
// Output (published to 13_summary_data/, isoform numbers shared by all files):
//   path "*_isoforms.tsv"                   intact isoforms (sample_id, cds_isoform, protein_isoform)
//   path "*_isoforms_incl_truncated.tsv"    all sequences with full-length/truncated status
//   path "*_isoform_mutations.tsv"          one row per mutation event
//   path "*_cds_isoforms.fasta"             one ungapped sequence per unique CDS isoform
//   path "*_protein_isoforms.fasta"         one ungapped sequence per unique protein isoform
process SummarizeIsoforms {
    tag "Isoform summary"
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6
    publishDir "${params.outdir}/13_summary_data", mode: 'copy'

    input:
    tuple val(query_id), path(cds_alignment), path(prot_alignment), path(integrity_tsvs, stageAs: 'integrity/*'), path(event_tsvs, stageAs: 'integrity/*')

    output:
    path("*.tsv"), emit: isoforms_tsv
    path("*_isoforms.fasta"), emit: isoforms_fasta

    script:
    """
    #!/bin/bash
    set -euo pipefail
    shopt -s nullglob
    integrity=( integrity/*_cds_integrity.tsv )
    events=( integrity/*_cds_integrity_events.tsv )
    summarize_isoforms.py '${query_id}' '${cds_alignment}' '${prot_alignment}' \
        --integrity "\${integrity[@]}" --events "\${events[@]}"
    """
}
