// Paralog check (--query_paralog): keep the hits of a query that belong to one
// of its known paralogs out of every downstream step.
//
// BLASTparalogs   -- per sample x query: BLAST the paralogs of the query against
//                    the assembly, with the same program and settings as the query
// ResolveParalogs -- per sample x query: assign every BLAST hit to the query or a
//                    paralog and pick the query copy to extract
//                    (paralog_check.py). Its filtered query BLAST table replaces
//                    the raw one as input of the hit extraction, so paralog loci
//                    never reach Exonerate, CDS integrity, isoforms, alignments,
//                    trees or the TE summary. It also lists the paralog copies on
//                    the hit scaffold, which the extraction keeps out of the padding.
// ParalogCheck    -- per query: merge the per-sample tables

// ---------------------------------------------------------------------------
// BLASTparalogs  (blast container)
// ---------------------------------------------------------------------------
process BLASTparalogs {
    tag "BLAST paralogs ${sample_id}"
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6
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
// ResolveParalogs  (spades container — has python3)
// ---------------------------------------------------------------------------
process ResolveParalogs {
    tag "Resolve paralogs ${sample_id}"
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3

    input:
    // Query and paralog BLAST tables are both named after the sample; separate
    // folders keep them apart
    tuple val(query_id), val(sample_id), path(query_blast, stageAs: 'query_blast/*'),
          path(paralog_blast, stageAs: 'paralog_blast/*'), path(query_fasta), path(paralog_fasta)

    output:
    tuple val(query_id), val(sample_id), path("${sample_id}_query_blast_results.tsv"), emit: query_hits
    tuple val(query_id), path("${sample_id}_paralog_check.tsv"),                      emit: check
    tuple val(query_id), path("${sample_id}_paralog_loci.tsv"),                       emit: loci

    script:
    """
    #!/bin/bash
    set -euo pipefail
    paralog_check.py "${sample_id}" '${query_id}' ${query_fasta} ${paralog_fasta} \\
        ${query_blast} ${paralog_blast} --ambiguous-ratio ${params.paralog_ambiguous_ratio}
    """
}

// ---------------------------------------------------------------------------
// ParalogCheck  (spades container)
// ---------------------------------------------------------------------------
process ParalogCheck {
    tag "Paralog check"
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    publishDir "${params.outdir}/13_summary_data", mode: 'copy'

    input:
    tuple val(query_id), path(check_tsvs, stageAs: 'check/*'), path(loci_tsvs, stageAs: 'loci/*')

    output:
    tuple val(query_id), path("*_paralog_check.tsv"), emit: table
    tuple val(query_id), path("*_paralog_loci.tsv"),  emit: loci

    script:
    """
    #!/bin/bash
    set -euo pipefail
    clean_id="\$(echo '${query_id}' | sed 's|[^a-zA-Z0-9_.-]|_|g')"
    # One header, then the sample rows sorted by sample
    for kind in check loci; do
        files=( \$kind/* )
        out="\${clean_id}_paralog_\${kind}.tsv"
        head -n 1 "\${files[0]}" > "\$out"
        for f in "\${files[@]}"; do tail -n +2 "\$f"; done | sort -t \$'\\t' -k1,1 >> "\$out"
    done
    """
}
