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
    tuple val(query_id), path("*_aligned_sequences.fasta"), emit: alignment
    tuple val(query_id), path("*_mafft_alignment_log.txt"), emit: log
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    # Sanitize query_id for filename use (replace special characters)
    clean_query_id="\$(echo '${query_id}' | sed 's|[^a-zA-Z0-9_.-]|_|g')"
    
    # Check if sequences file exists and is not empty
    if [ ! -s "${sequences_fasta}" ]; then
        echo "No sequences found for alignment" > "\${clean_query_id}_mafft_alignment_log.txt"
        touch "\${clean_query_id}_aligned_sequences.fasta"
        exit 0
    fi
    
    # Count the number of sequences in the FASTA file
    seq_count=\$(grep -c "^>" "${sequences_fasta}" || echo "0")
    
    # Check if we have enough sequences for alignment (at least 2)
    if [ "\$seq_count" -lt 2 ]; then
        echo "Insufficient sequences for alignment (found \$seq_count sequences, need at least 2)" > "\${clean_query_id}_mafft_alignment_log.txt"
        cp "${sequences_fasta}" "\${clean_query_id}_aligned_sequences.fasta"
        exit 0
    fi
    
    echo "Starting MAFFT alignment with \$seq_count sequences" > "\${clean_query_id}_mafft_alignment_log.txt"
    echo "Input file: ${sequences_fasta}" >> "\${clean_query_id}_mafft_alignment_log.txt"

    # Detect whether the first sequence is protein or nucleotide.
    # When the BLAST query is a protein (tblastn), CombineBlastHits prepends
    # the protein query to nucleotide hit sequences.  MAFFT cannot align that
    # mixture in nucleotide mode, so we strip the protein query and align only
    # the nucleotide hit sequences.
    first_seq="\$(awk 'BEGIN{found=0}/^>/{if(found)exit;found=1;next}found{printf "%s",\$0}' "${sequences_fasta}" | tr -d '[:space:]' | tr '[:lower:]' '[:upper:]')"
    first_len=\${#first_seq}
    if [ "\$first_len" -gt 0 ]; then
        nuc_chars=\$(printf '%s' "\${first_seq}" | tr -cd 'ATGCUN' | wc -c)
        nuc_pct=\$((nuc_chars * 100 / first_len))
    else
        nuc_pct=100
    fi

    if [ "\$nuc_pct" -lt 95 ]; then
        # First sequence is a protein query — keep only the nucleotide hit sequences
        echo "Protein query detected as first sequence (nuc% = \${nuc_pct}); aligning nucleotide hit sequences only" >> "\${clean_query_id}_mafft_alignment_log.txt"
        awk '/^>/{count++} count>=2' "${sequences_fasta}" > "\${clean_query_id}_input.fasta"
        align_input="\${clean_query_id}_input.fasta"
        nuc_count=\$(grep -c "^>" "\${clean_query_id}_input.fasta" || echo "0")
        if [ "\${nuc_count}" -lt 2 ]; then
            echo "Insufficient nucleotide sequences after removing protein query (found \${nuc_count})" >> "\${clean_query_id}_mafft_alignment_log.txt"
            touch "\${clean_query_id}_aligned_sequences.fasta"
            exit 0
        fi
    else
        align_input="${sequences_fasta}"
    fi

    echo "Command: mafft --auto --adjustdirection \${align_input}" >> "\${clean_query_id}_mafft_alignment_log.txt"
    echo "---" >> "\${clean_query_id}_mafft_alignment_log.txt"

    # Run MAFFT with auto mode and direction adjustment
    mafft --thread ${task.cpus} --auto --adjustdirection "\${align_input}" > "\${clean_query_id}_aligned_sequences.fasta" 2>> "\${clean_query_id}_mafft_alignment_log.txt"
    
    # Check if alignment was successful
    if [ -s "\${clean_query_id}_aligned_sequences.fasta" ]; then
        aligned_count=\$(grep -c "^>" "\${clean_query_id}_aligned_sequences.fasta" || echo "0")
        echo "Alignment completed successfully with \$aligned_count sequences" >> "\${clean_query_id}_mafft_alignment_log.txt"
    else
        echo "WARNING: Alignment output file is empty or missing" >> "\${clean_query_id}_mafft_alignment_log.txt"
        touch "\${clean_query_id}_aligned_sequences.fasta"
    fi
    """
}
