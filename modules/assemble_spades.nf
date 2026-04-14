// Fixed SPAdes process for specific file outputs
process SPAdes {
    tag "SPAdes assembly: ${sample_id}"
    errorStrategy 'retry'
    maxRetries 3

    publishDir "${params.outdir}/assemblies", mode: 'copy', pattern: "${sample_id}/contigs.fasta"
    publishDir "${params.outdir}/assemblies", mode: 'copy', pattern: "${sample_id}/scaffolds.fasta"
    publishDir "${params.outdir}/assemblies", mode: 'copy', pattern: "${sample_id}/spades.log"
    
    input:
    tuple val(sample_id), path(read1), path(read2), val(source)
    
    output:
    tuple val(sample_id), path("${sample_id}/contigs.fasta"), emit: contigs
    tuple val(sample_id), path("${sample_id}/scaffolds.fasta"), emit: scaffolds
    path "${sample_id}/spades.log", emit: log
    
    script:
    """
    spades.py \
        --pe1-1 ${read1} \
        --pe1-2 ${read2} \
        -o ${sample_id} \
        -t ${task.cpus} \
        -m \$((${task.memory.toGiga()})) \
        --isolate
    
    # Clean up intermediate files to save space
    if [ -d "${sample_id}/K21" ]; then
        rm -rf ${sample_id}/K*/
    fi
    if [ -d "${sample_id}/corrected" ]; then
        rm -rf ${sample_id}/corrected/
    fi
    if [ -d "${sample_id}/tmp" ]; then
        rm -rf ${sample_id}/tmp/
    fi
    
    # Only delete input files if they were downloaded from SRA
    # User-provided files (source='local') are preserved
    if [ "${source}" = "SRA" ]; then
        echo "Deleting SRA-downloaded files after assembly..."
        for file in "${read1}" "${read2}"; do
            if [ -L "\$file" ]; then
                target=\$(readlink -f "\$file" 2>/dev/null)
                [ -n "\$target" ] && [ -f "\$target" ] && rm "\$target" || true
            elif [ -f "\$file" ]; then
                rm "\$file" || true
            fi
        done
    else
        echo "Preserving user-provided files (source: ${source})"
    fi
    """
}

// SPAdes process for single-end reads
process SPAdesSE {
    tag "SPAdes SE assembly: ${sample_id}"
    errorStrategy 'retry'
    maxRetries 3

    publishDir "${params.outdir}/assemblies", mode: 'copy', pattern: "${sample_id}/contigs.fasta"
    publishDir "${params.outdir}/assemblies", mode: 'copy', pattern: "${sample_id}/scaffolds.fasta"
    publishDir "${params.outdir}/assemblies", mode: 'copy', pattern: "${sample_id}/spades.log"
    
    input:
    tuple val(sample_id), path(read), val(source)
    
    output:
    tuple val(sample_id), path("${sample_id}/contigs.fasta"), emit: contigs
    tuple val(sample_id), path("${sample_id}/scaffolds.fasta"), emit: scaffolds
    path "${sample_id}/spades.log", emit: log
    
    script:
    """
    spades.py \
        --s1 ${read} \
        -o ${sample_id} \
        -t ${task.cpus} \
        -m \$((${task.memory.toGiga()})) \
        --isolate
    
    # Clean up intermediate files to save space
    if [ -d "${sample_id}/K21" ]; then
        rm -rf ${sample_id}/K*/
    fi
    if [ -d "${sample_id}/corrected" ]; then
        rm -rf ${sample_id}/corrected/
    fi
    if [ -d "${sample_id}/tmp" ]; then
        rm -rf ${sample_id}/tmp/
    fi
    
    # Only delete input file if it was downloaded from SRA
    # User-provided files (source='local') are preserved
    if [ "${source}" = "SRA" ]; then
        echo "Deleting SRA-downloaded file after assembly..."
        if [ -L "${read}" ]; then
            target=\$(readlink -f "${read}" 2>/dev/null)
            [ -n "\$target" ] && [ -f "\$target" ] && rm "\$target" || true
        elif [ -f "${read}" ]; then
            rm "${read}" || true
        fi
    else
        echo "Preserving user-provided file (source: ${source})"
    fi
    """
}

