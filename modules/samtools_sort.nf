process samtoolsSort {
    tag "Sort/index BAM: ${sample_id}"
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    publishDir "${params.outdir}/9_bwa_mapping", mode: 'copy'

    input:
    tuple val(sample_id), path(sample_sam)

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai"), emit: bam
    path "${sample_id}_flagstat.json", emit: report

    script:
    """
    set -euo pipefail

    samtools flagstat -@ ${task.cpus} -O json $sample_sam > ${sample_id}_flagstat.json
    samtools view -@ ${task.cpus} -Sb -q 10 $sample_sam | samtools sort -@ ${task.cpus} -o ${sample_id}_sorted.bam -
    samtools index -@ ${task.cpus} ${sample_id}_sorted.bam

    # Delete the SAM file after successful BAM conversion (resolve symlink first)
    sam_target=\$(readlink -f "$sample_sam")
    [ -n "\$sam_target" ] && [ -f "\$sam_target" ] && rm "\$sam_target" || true
    """
}
