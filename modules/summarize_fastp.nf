process SummarizeFastpStats {
    tag "Summarizing FASTP stats"
    errorStrategy 'retry'
    maxRetries 3
    publishDir "${params.outdir}/assemblies/fastp_reports", mode: 'copy'

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
