process trimSequencesPE {
    tag "FASTP PE trimming"
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    publishDir "${params.outdir}/fastp_stats", mode: 'copy', pattern: "*.json"
        
    input:
    tuple val(sample_id), path(read1), path(read2), val(source)
    
    output:
    tuple val(sample_id), path("${sample_id}_*_trimmed.fastq.gz"), emit: reads
    path "${sample_id}_PE_fastp.json", emit: report
    
    script:
    """
    fastp \
        -w $task.cpus \
        -i $read1 \
        -I $read2 \
        -o ${sample_id}_1_trimmed.fastq.gz \
        -O ${sample_id}_2_trimmed.fastq.gz \
        --json ${sample_id}_PE_fastp.json
    
    # Only delete files if they were downloaded from SRA
    # User-provided files (source='local') are preserved
    if [ "${source}" = "SRA" ]; then
        echo "Deleting SRA-downloaded files after trimming..."
        for file in "$read1" "$read2"; do
            if [ -L "\$file" ]; then
                # Resolve symlink to actual file
                target=\$(readlink -f "\$file" 2>/dev/null)
                [ -n "\$target" ] && [ -f "\$target" ] && rm "\$target" || true
            elif [ -f "\$file" ]; then
                # Direct file deletion
                rm "\$file" || true
            fi
        done
    else
        echo "Preserving user-provided files (source: ${source})"
    fi
    """
}

process trimSequencesSE {
    tag "FASTP SE trimming"
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    publishDir "${params.outdir}/fastp_stats", mode: 'copy', pattern: "*.json"

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
        -o ${sample_id}_trimmed.fastq.gz \
        --json ${sample_id}_SE_fastp.json
    
    # Only delete files if they were downloaded from SRA
    # User-provided files (source='local') are preserved
    if [ "${source}" = "SRA" ]; then
        echo "Deleting SRA-downloaded file after trimming..."
        if [ -L "$r1" ]; then
            # Resolve symlink to actual file
            target=\$(readlink -f "$r1" 2>/dev/null)
            [ -n "\$target" ] && [ -f "\$target" ] && rm "\$target" || true
        elif [ -f "$r1" ]; then
            # Direct file deletion
            rm "$r1" || true
        fi
    else
        echo "Preserving user-provided file (source: ${source})"
    fi
    """
}
