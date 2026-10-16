// Assign isoform IDs to unique CDS and protein haplotypes derived from
// per-query MAFFT-aligned FASTA files, and flag truncated isoforms using the
// per-sample CDS integrity tables (AssessCdsIntegrity). Isoforms are numbered
// by frequency: intact isoforms by number of samples carrying them, then a
// reference-only isoform, then truncated/extended-only isoforms.
//
// Input:
//   tuple val(query_id), path(cds_alignment), path(prot_alignment),
//         path(integrity_tsvs), path(integrity_event_tsvs)
//
// Output (published to 13_summary_data/, isoform numbers shared by all files):
//   path "*_isoforms.tsv"                   intact isoforms (sample_id, cds_isoform, protein_isoform,
//                                           cds_sequence, protein_sequence)
//   path "*_isoforms_incl_truncated.tsv"    all sequences with full-length/truncated status
//   path "*_isoform_mutations.tsv"          one row per mutation event
//   path "*_cds_isoforms.fasta"             one ungapped sequence per unique CDS isoform
//   path "*_protein_isoforms.fasta"         one ungapped sequence per unique protein isoform
//   path "*_{cds,protein}_isoforms_aligned.fasta"
//                                           one aligned row per isoform (input for the isoform trees)
process SummarizeIsoforms {
    tag "Isoform summary"
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6
    publishDir "${params.outdir}/13_summary_data", mode: 'copy'

    input:
    tuple val(query_id), path(cds_alignment), path(prot_alignment), path(integrity_tsvs, stageAs: 'integrity/*'), path(event_tsvs, stageAs: 'integrity/*'),
          path(blast_summary, stageAs: 'blast/*')

    output:
    path("*.tsv"), emit: isoforms_tsv
    tuple val(query_id), path("*_cds_isoforms.fasta"),             emit: cds_isoforms
    tuple val(query_id), path("*_protein_isoforms.fasta"),         emit: protein_isoforms
    tuple val(query_id), path("*_cds_isoforms_aligned.fasta"),     emit: cds_isoforms_aligned
    tuple val(query_id), path("*_protein_isoforms_aligned.fasta"), emit: protein_isoforms_aligned

    script:
    """
    #!/bin/bash
    set -euo pipefail
    shopt -s nullglob
    integrity=( integrity/*_cds_integrity.tsv )
    events=( integrity/*_cds_integrity_events.tsv )
    summarize_isoforms.py '${query_id}' '${cds_alignment}' '${prot_alignment}' \
        --integrity "\${integrity[@]}" --events "\${events[@]}" --blast-summary ${blast_summary}
    """
}
