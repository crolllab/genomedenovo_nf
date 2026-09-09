process bwaMap {
    tag "BWA-mem2 mapping"
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6

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
    // Build bwa-mem2 input command: handle both PE (2 files) and SE (1 file).
    // The read paths must stay UNQUOTED: Nextflow backslash-escapes staged path
    // names (parentheses and spaces included), and wrapping that in quotes would
    // turn the escapes into literal characters. Only $sample_id, which is a plain
    // val and is not escaped, needs quoting.
    def reads_cmd = trimmed_reads.size() == 2 ? "${trimmed_reads[0]} ${trimmed_reads[1]}" : "${trimmed_reads[0]}"
    """
    set -euo pipefail

    bwa-mem2 mem -t $task.cpus $reference $reads_cmd > "${sample_id}.sam"
    """
}
