process SummarizeFastpStats {
    tag "Summarizing FASTP stats"
    errorStrategy 'retry'
    maxRetries 3
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
