// Process to perform multiple sequence alignment on combined BLAST hit sequences using MAFFT
process AlignMAFFT {
    tag "MAFFT alignment for query: ${query_id}"
    errorStrategy 'ignore'
    memory '8GB'
    cpus 8
    publishDir "${params.outdir}/${query_id}", mode: 'copy'
    
    input:
    tuple val(query_id), path(sequences_fasta)
    
    output:
    tuple val(query_id), path("${query_id}_aligned_sequences.fasta"), emit: alignment
    tuple val(query_id), path("${query_id}_mafft_alignment_log.txt"), emit: log
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    # Check if sequences file exists and is not empty
    if [ ! -s "${sequences_fasta}" ]; then
        echo "No sequences found for alignment" > ${query_id}_mafft_alignment_log.txt
        touch ${query_id}_aligned_sequences.fasta
        exit 0
    fi
    
    # Count the number of sequences in the FASTA file
    seq_count=\$(grep -c "^>" "${sequences_fasta}" || echo "0")
    
    # Check if we have enough sequences for alignment (at least 2)
    if [ "\$seq_count" -lt 2 ]; then
        echo "Insufficient sequences for alignment (found \$seq_count sequences, need at least 2)" > ${query_id}_mafft_alignment_log.txt
        # Copy the input as output if only one sequence
        cp "${sequences_fasta}" "${query_id}_aligned_sequences.fasta"
        exit 0
    fi
    
    echo "Starting MAFFT alignment with \$seq_count sequences" > ${query_id}_mafft_alignment_log.txt
    echo "Input file: ${sequences_fasta}" >> ${query_id}_mafft_alignment_log.txt
    echo "Command: mafft --auto --adjustdirection ${sequences_fasta}" >> ${query_id}_mafft_alignment_log.txt
    echo "---" >> ${query_id}_mafft_alignment_log.txt
    
    # Run MAFFT with auto mode and adjust directions
    mafft --thread ${task.cpus} --auto --adjustdirection "${sequences_fasta}" > "${query_id}_aligned_sequences.fasta" 2>> ${query_id}_mafft_alignment_log.txt
    
    # Check if alignment was successful
    if [ -s "${query_id}_aligned_sequences.fasta" ]; then
        aligned_count=\$(grep -c "^>" "${query_id}_aligned_sequences.fasta" || echo "0")
        echo "Alignment completed successfully with \$aligned_count sequences" >> ${query_id}_mafft_alignment_log.txt
    else
        echo "WARNING: Alignment output file is empty or missing" >> ${query_id}_mafft_alignment_log.txt
        touch ${query_id}_aligned_sequences.fasta
    fi
    """
}
