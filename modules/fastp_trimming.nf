process trimSequencesPE {
    tag "FASTP PE trimming"
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6
    publishDir "${params.outdir}/1_fastq_stats", mode: 'copy', pattern: "*.json"
        
    input:
    tuple val(sample_id), path(read1), path(read2), val(source)
    
    output:
    tuple val(sample_id), path("${sample_id}_*_trimmed.fastq.gz"), emit: reads
    path "${sample_id}_PE_fastp.json", emit: report
    
    script:
    // Output names are quoted because sample IDs may contain characters that
    // are syntactically significant to bash, e.g. "Iran_2024_Kch(2)21-2c" or
    // "Iran_2024_RpH2_a_YMDA (Yeast)". Staged input paths ($read1/$read2) are
    // already backslash-escaped by Nextflow and must stay unquoted.
    """
    fastp \
        -w $task.cpus \
        -i $read1 \
        -I $read2 \
        -o "${sample_id}_1_trimmed.fastq.gz" \
        -O "${sample_id}_2_trimmed.fastq.gz" \
        --json "${sample_id}_PE_fastp.json"
    """
}

process trimSequencesSE {
    tag "FASTP SE trimming"
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6
    publishDir "${params.outdir}/1_fastq_stats", mode: 'copy', pattern: "*.json"

    beforeScript """
        mkdir -p "\$PWD/tmp"
        export TMPDIR="\$PWD/tmp"
    """

    input:
    tuple val(sample_id), path(r1), val(source)
    
    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz"), emit: reads
    path "${sample_id}_SE_fastp.json", emit: report
    
    script:
    """
    fastp \
        -w $task.cpus \
        -i $r1 \
        -o "${sample_id}_trimmed.fastq.gz" \
        --json "${sample_id}_SE_fastp.json"
    """
}
