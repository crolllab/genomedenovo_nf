// Fixed SPAdes process for specific file outputs
process SPAdes {
    tag "SPAdes assembly: ${sample_id}"
    errorStrategy 'ignore'
    memory '64GB'
    cpus 8
    publishDir "${params.outdir}/spades_assemblies", mode: 'copy', pattern: "${sample_id}/contigs.fasta"
    publishDir "${params.outdir}/spades_assemblies", mode: 'copy', pattern: "${sample_id}/scaffolds.fasta"
    publishDir "${params.outdir}/spades_stats", mode: 'copy', pattern: "${sample_id}/spades.log"
    
    input:
    tuple val(sample_id), path(read1), path(read2)
    
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
        --careful
        --fungal
    
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
    
    # Safely delete the input files
    for file in "${read1}" "${read2}"; do
        if [ -L "\$file" ]; then
            # It's a symlink - try to delete the target
            target=\$(readlink -f "\$file" 2>/dev/null)
            if [ -n "\$target" ] && [ -f "\$target" ]; then
                rm "\$target" 2>/dev/null || echo "Warning: Could not delete \$target (write-protected or permission denied)"
            fi
        elif [ -f "\$file" ]; then
            # It's a regular file - try to delete it directly
            rm "\$file" 2>/dev/null || echo "Warning: Could not delete \$file (write-protected or permission denied)"
        fi
    done
    """
}

// SPAdes process for single-end reads
process SPAdesSE {
    tag "SPAdes SE assembly: ${sample_id}"
    errorStrategy 'ignore'
    memory '64GB'
    cpus 8
    publishDir "${params.outdir}/spades_assemblies", mode: 'copy', pattern: "${sample_id}/contigs.fasta"
    publishDir "${params.outdir}/spades_assemblies", mode: 'copy', pattern: "${sample_id}/scaffolds.fasta"
    publishDir "${params.outdir}/spades_stats", mode: 'copy', pattern: "${sample_id}/spades.log"
    
    input:
    tuple val(sample_id), path(read)
    
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
        --careful
    
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
    
    # Safely delete the input file
    if [ -L "${read}" ]; then
        # It's a symlink - try to delete the target
        target=\$(readlink -f "${read}" 2>/dev/null)
        if [ -n "\$target" ] && [ -f "\$target" ]; then
            rm "\$target" 2>/dev/null || echo "Warning: Could not delete \$target (write-protected or permission denied)"
        fi
    elif [ -f "${read}" ]; then
        # It's a regular file - try to delete it directly
        rm "${read}" 2>/dev/null || echo "Warning: Could not delete ${read} (write-protected or permission denied)"
    fi
    """
}

