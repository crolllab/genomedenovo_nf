process bwaMap {
    tag "BWA-mem2 mapping: ${sample_id}"
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3

    input:
    path reference
    file "${reference}.amb"
    file "${reference}.ann"
    file "${reference}.bwt.2bit.64"
    file "${reference}.pac"
    file "${reference}.0123"
    tuple val(sample_id), path(trimmed_reads)

    output:
    tuple val(sample_id), path("${sample_id}.sam"), emit: sam

    script:
    // Build bwa-mem2 input command: handle both PE (2 files) and SE (1 file)
    def reads_cmd = trimmed_reads.size() == 2 ? "${trimmed_reads[0]} ${trimmed_reads[1]}" : "${trimmed_reads[0]}"
    """
    set -euo pipefail

    bwa-mem2 mem -t $task.cpus $reference $reads_cmd > ${sample_id}.sam
    """
}
