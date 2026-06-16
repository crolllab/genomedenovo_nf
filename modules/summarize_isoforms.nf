// Assign consecutive isoform IDs to unique CDS and protein haplotypes
// derived from per-query MAFFT-aligned FASTA files.
//
// Input:
//   tuple val(query_id), path(cds_alignment), path(prot_alignment)
//
// Output:
//   path "*_isoforms.tsv"   published to 12_summary_data/
//
// One output file per query; columns: sample_id, cds_isoform, protein_isoform.
process SummarizeIsoforms {
    tag "Isoform summary"
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6
    publishDir "${params.outdir}/13_summary_data", mode: 'copy'

    input:
    tuple val(query_id), path(cds_alignment), path(prot_alignment)

    output:
    path("*_isoforms.tsv"), emit: isoforms_tsv

    script:
    """
    #!/bin/bash
    set -euo pipefail
    summarize_isoforms.py '${query_id}' '${cds_alignment}' '${prot_alignment}'
    """
}
