process bwaIndex {
    tag "BWA-mem2 index"
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6

    input:
    path reference

    output:
    file "${reference}.amb"
    file "${reference}.ann"
    file "${reference}.bwt.2bit.64"
    file "${reference}.pac"
    file "${reference}.0123"

    script:
    """
    bwa-mem2 index $reference
    """
}
