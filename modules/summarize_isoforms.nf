// Assign consecutive isoform IDs to unique CDS and protein haplotypes
// derived from per-query MAFFT-aligned FASTA files.
//
// Input:
//   tuple val(query_id), path(cds_alignment), path(prot_alignment)
//
// Output:
//   path "*_isoforms.tsv"            published to 13_summary_data/
//   path "*_cds_isoforms.fasta"      one ungapped sequence per unique CDS isoform
//   path "*_protein_isoforms.fasta"  one ungapped sequence per unique protein isoform
//
// TSV columns: sample_id, cds_isoform, protein_isoform.
process SummarizeIsoforms {
    tag "Isoform summary"
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6
    publishDir "${params.outdir}/13_summary_data", mode: 'copy'

    input:
    tuple val(query_id), path(cds_alignment), path(prot_alignment)

    output:
    path("*_isoforms.tsv"), emit: isoforms_tsv
    path("*_isoforms.fasta"), emit: isoforms_fasta

    script:
    """
    #!/bin/bash
    set -euo pipefail
    summarize_isoforms.py '${query_id}' '${cds_alignment}' '${prot_alignment}'
    """
}
