process samtoolsSort {
    tag "Sort/index BAM"
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6
    publishDir "${params.outdir}/10_bwa_mapping", mode: 'copy'

    input:
    tuple val(sample_id), path(sample_sam)

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai"), emit: bam
    path "${sample_id}_flagstat.json", emit: report

    script:
    // Output names are quoted so sample IDs containing parentheses or spaces
    // are handled correctly. $sample_sam is escaped by Nextflow and stays bare.
    //
    // The upstream .sam is deliberately NOT deleted here: resolving the staged
    // symlink and removing its target destroys the bwaMap task's cached output,
    // which breaks -resume and any later retry of this task. Reclaim space with
    // `nextflow clean` once the run has completed successfully.
    //
    // samtools -@ counts threads in addition to the main one; view streams
    // uncompressed BAM to sort, so it needs no compression threads.
    def extra_threads = Math.max(0, task.cpus - 1)
    def sort_mem_mb   = Math.max(256, (task.memory.toMega() * 0.6 / task.cpus) as long)
    """
    set -euo pipefail

    samtools flagstat -@ ${extra_threads} -O json $sample_sam > "${sample_id}_flagstat.json"
    samtools view -u -q 10 $sample_sam | samtools sort -@ ${extra_threads} -m ${sort_mem_mb}M -o "${sample_id}_sorted.bam" -
    samtools index -@ ${extra_threads} "${sample_id}_sorted.bam"
    """
}
