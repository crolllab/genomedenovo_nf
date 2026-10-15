// Paralog check (--query_paralog): is the hit reported for a query closer to the
// query or to one of its known paralogs?
//
// BLASTparalogs -- per sample x query: BLAST the paralogs of the query against
//                  the assembly, with the same program and settings as the query
// ParalogCheck  -- per query: at the hit locus of every sample (the merged BLAST
//                  span used for hit extraction), compare the BLAST scores of the
//                  query and of each paralog (paralog_check.py)

// ---------------------------------------------------------------------------
// BLASTparalogs  (blast container)
// ---------------------------------------------------------------------------
process BLASTparalogs {
    tag "BLAST paralogs ${sample_id}"
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6
    memory '4GB'
    cpus 4
    publishDir { "${params.outdir}/6_blast_hits/${query_id}/paralog_check" }, mode: 'copy'

    input:
    tuple val(sample_id), path(scaffolds_fasta), val(query_id), path(paralog_fasta), val(seq_type)

    output:
    tuple val(query_id), val(sample_id), path("${sample_id}_paralog_blast_results.tsv"), emit: results

    script:
    // Same program and settings as BLASTsearch, so that scores are comparable
    def blast_tool = (seq_type == "protein") ? "tblastn" : "blastn"
    def blast_task = (seq_type == "protein") ? "" : "-task blastn"
    """
    set -euo pipefail

    makeblastdb -in ${scaffolds_fasta} -dbtype nucl -out assembly_db

    ${blast_tool} ${blast_task} -query ${paralog_fasta} \\
           -db assembly_db \\
           -out "${sample_id}_paralog_blast_results.tsv" \\
           -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \\
           -num_threads ${task.cpus} \\
           -evalue 1e-5

    rm -f assembly_db.*
    """
}

// ---------------------------------------------------------------------------
// ParalogCheck  (spades container — has python3)
// ---------------------------------------------------------------------------
process ParalogCheck {
    tag "Paralog check"
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    publishDir "${params.outdir}/13_summary_data", mode: 'copy'

    input:
    // Query and paralog BLAST tables are both named after the sample; separate
    // folders keep them apart
    tuple val(query_id), path(query_fasta), path(paralog_fasta), path(hits_summary),
          path(query_blast, stageAs: 'query_blast/*'), path(paralog_blast, stageAs: 'paralog_blast/*')

    output:
    tuple val(query_id), path("*_paralog_check.tsv"), emit: table

    script:
    """
    #!/bin/bash
    set -euo pipefail
    paralog_check.py '${query_id}' ${query_fasta} ${paralog_fasta} ${hits_summary} \\
        --query-blast query_blast/* --paralog-blast paralog_blast/* \\
        --ambiguous-ratio ${params.paralog_ambiguous_ratio}
    """
}
