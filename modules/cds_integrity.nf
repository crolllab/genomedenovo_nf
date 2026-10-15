// CDS integrity: classify every Exonerate gene model as full-length or
// truncated, separating nucleotide-level deletions (query segment absent from
// the contig) from ORF mutations (frameshifts, premature stops, lost start/stop
// codons, non-canonical splice sites, segments present but outside the ORF).
//
// AssessCdsIntegrity -- per sample x query: classify_cds_integrity.py (exonerate container:
//                       python3 + exonerate for re-aligning unaligned query segments).
//                       The tables are not published here; SummarizeIsoforms folds them
//                       into the per-query isoform summaries in 13_summary_data/.

// ---------------------------------------------------------------------------
// AssessCdsIntegrity  (exonerate container)
// ---------------------------------------------------------------------------
process AssessCdsIntegrity {
    tag "CDS integrity ${sample_id}"
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3

    input:
    // The full assembly gives exact contig ends and contig coordinates; it is
    // staged in its own folder so its name cannot clash with the hit FASTA.
    tuple val(query_id), val(sample_id), path(hit_fasta), path(gff3_file), path(query_fasta), val(seq_type), path(assembly, stageAs: 'assembly/*')

    output:
    tuple val(query_id), path("*_cds_integrity.tsv"),        emit: summary
    tuple val(query_id), path("*_cds_integrity_events.tsv"), emit: events

    script:
    def prefix = "${sample_id}.${query_id.replaceAll(/[^A-Za-z0-9_.-]/, '_')}"
    """
    #!/bin/bash
    set -euo pipefail

    # Staged input paths are backslash-escaped by Nextflow and must stay unquoted
    classify_cds_integrity.py \
        --sample "${sample_id}" \
        --query-id "${query_id}" \
        --query ${query_fasta} \
        --seq-type "${seq_type}" \
        --exonerate ${gff3_file} \
        --target ${hit_fasta} \
        --assembly ${assembly} \
        --out-prefix "${prefix}"
    """
}
