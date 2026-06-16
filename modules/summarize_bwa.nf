process SummarizeBwaStats {
    tag "Summarizing BWA mapping stats"
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6
    publishDir "${params.outdir}/10_bwa_mapping", mode: 'copy'

    input:
    path json_files
    path r_script

    output:
    path "bwa_mapping_summary.tsv", emit: tsv
    path "bwa_mapping_summary.txt", emit: txt

    script:
    """
    Rscript ${r_script}
    """
}
