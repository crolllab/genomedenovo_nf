process SummarizeFastpStats {
    tag "Summarizing FASTP stats"
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6
    publishDir "${params.outdir}/1_fastq_stats", mode: 'copy'

    input:
    path json_files
    path r_script

    output:
    path "fastp_summary.tsv", emit: tsv
    path "fastp_summary.txt", emit: txt

    script:
    """
    Rscript ${r_script}
    """
}
