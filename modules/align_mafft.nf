// Process to perform multiple sequence alignment on combined BLAST hit sequences using MAFFT
process AlignMAFFT {
    tag "MAFFT alignment: all_blast_hits"
    errorStrategy 'ignore'
    memory '4GB'
    cpus 2
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path sequences_fasta
    
    output:
    path "all_sequences_aligned.fasta", emit: alignment
    path "mafft_alignment_log.txt", emit: log
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    # Check if sequences file exists and is not empty
    if [ ! -s "${sequences_fasta}" ]; then
        echo "No sequences found for alignment" > mafft_alignment_log.txt
        touch all_sequences_aligned.fasta
        exit 0
    fi
    
    # Count the number of sequences in the FASTA file
    seq_count=\$(grep -c "^>" "${sequences_fasta}" || echo "0")
    
    # Check if we have enough sequences for alignment (at least 2)
    if [ "\$seq_count" -lt 2 ]; then
        echo "Insufficient sequences for alignment (found \$seq_count sequences, need at least 2)" > mafft_alignment_log.txt
        # Copy the input as output if only one sequence
        cp "${sequences_fasta}" "all_sequences_aligned.fasta"
        exit 0
    fi
    
    echo "Starting MAFFT alignment with \$seq_count sequences" > mafft_alignment_log.txt
    echo "Input file: ${sequences_fasta}" >> mafft_alignment_log.txt
    echo "Command: mafft --auto --adjustdirection ${sequences_fasta}" >> mafft_alignment_log.txt
    echo "---" >> mafft_alignment_log.txt
    
    # Run MAFFT with auto mode and adjust directions
    mafft --auto --adjustdirection "${sequences_fasta}" > "all_sequences_aligned.fasta" 2>> mafft_alignment_log.txt
    
    # Check if alignment was successful
    if [ -s "all_sequences_aligned.fasta" ]; then
        aligned_count=\$(grep -c "^>" "all_sequences_aligned.fasta" || echo "0")
        echo "Alignment completed successfully with \$aligned_count sequences" >> mafft_alignment_log.txt
    else
        echo "ERROR: MAFFT alignment failed" >> mafft_alignment_log.txt
        touch "all_sequences_aligned.fasta"
    fi
    """
}
