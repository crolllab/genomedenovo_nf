// Process to perform multiple sequence alignment on combined BLAST hit sequences using MAFFT
process AlignMAFFT {
    tag "MAFFT alignment"
    errorStrategy 'ignore'
    publishDir { "${params.outdir}/alignments/${query_id}" }, mode: 'copy'
    
    input:
    tuple val(query_id), path(sequences_fasta)
    
    output:
    tuple val(query_id), path("*_aligned_sequences.fasta"), emit: alignment
    tuple val(query_id), path("*_mafft_alignment_log.txt"), emit: log
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    # Derive output prefix from the input filename so that multiple MAFFT
    # invocations for the same query (BLAST hits, Exonerate CDS, Exonerate protein) each
    # produce distinct output files and do not overwrite each other.
    file_stem="\$(basename '${sequences_fasta}' .fasta)"
    
    # Check if sequences file exists and is not empty
    if [ ! -s "${sequences_fasta}" ]; then
        echo "No sequences found for alignment" > "\${file_stem}_mafft_alignment_log.txt"
        touch "\${file_stem}_aligned_sequences.fasta"
        exit 0
    fi
    
    # Count the number of sequences in the FASTA file
    seq_count=\$(grep -c "^>" "${sequences_fasta}" || echo "0")
    
    # Check if we have enough sequences for alignment (at least 2)
    if [ "\$seq_count" -lt 2 ]; then
        echo "Insufficient sequences for alignment (found \$seq_count sequences, need at least 2)" > "\${file_stem}_mafft_alignment_log.txt"
        cp "${sequences_fasta}" "\${file_stem}_aligned_sequences.fasta"
        exit 0
    fi
    
    echo "Starting MAFFT alignment with \$seq_count sequences" > "\${file_stem}_mafft_alignment_log.txt"
    echo "Input file: ${sequences_fasta}" >> "\${file_stem}_mafft_alignment_log.txt"

    # Detect whether this is a mixed case (protein query prepended to nucleotide
    # hit sequences, as in tblastn BLAST output).  Only strip the first sequence
    # when the first is protein-type AND the second is nucleotide-type.
    # If both are the same type (all-nucleotide or all-protein), align all.
    first_seq="\$(awk 'BEGIN{found=0}/^>/{if(found)exit;found=1;next}found{printf "%s",\$0}' "${sequences_fasta}" | tr -d '[:space:]' | tr '[:lower:]' '[:upper:]')"
    second_seq="\$(awk 'BEGIN{n=0}/^>/{n++;next} n==2{printf "%s",\$0}' "${sequences_fasta}" | tr -d '[:space:]' | tr '[:lower:]' '[:upper:]')"

    first_len=\${#first_seq}
    second_len=\${#second_seq}

    first_is_nuc=true
    second_is_nuc=true

    if [ "\$first_len" -gt 0 ]; then
        nuc_chars=\$(printf '%s' "\${first_seq}" | tr -cd 'ATGCUN' | wc -c)
        pct=\$((nuc_chars * 100 / first_len))
        [ "\$pct" -lt 95 ] && first_is_nuc=false
    fi
    if [ "\$second_len" -gt 0 ]; then
        nuc_chars=\$(printf '%s' "\${second_seq}" | tr -cd 'ATGCUN' | wc -c)
        pct=\$((nuc_chars * 100 / second_len))
        [ "\$pct" -lt 95 ] && second_is_nuc=false
    fi

    if [ "\${first_is_nuc}" = "false" ] && [ "\${second_is_nuc}" = "true" ]; then
        # Protein query prepended to nucleotide hit sequences — align only hits
        echo "Protein query + nucleotide hits detected; aligning nucleotide hit sequences only" >> "\${file_stem}_mafft_alignment_log.txt"
        awk '/^>/{count++} count>=2' "${sequences_fasta}" > "\${file_stem}_input.fasta"
        align_input="\${file_stem}_input.fasta"
        nuc_count=\$(grep -c "^>" "\${file_stem}_input.fasta" || echo "0")
        if [ "\${nuc_count}" -lt 2 ]; then
            echo "Insufficient nucleotide sequences after removing protein query (found \${nuc_count})" >> "\${file_stem}_mafft_alignment_log.txt"
            touch "\${file_stem}_aligned_sequences.fasta"
            exit 0
        fi
    else
        align_input="${sequences_fasta}"
    fi

    echo "Command: mafft --auto --adjustdirection \${align_input}" >> "\${file_stem}_mafft_alignment_log.txt"
    echo "---" >> "\${file_stem}_mafft_alignment_log.txt"

    # Run MAFFT with auto mode and direction adjustment
    mafft --thread ${task.cpus} --auto --adjustdirection "\${align_input}" > "\${file_stem}_aligned_sequences.fasta" 2>> "\${file_stem}_mafft_alignment_log.txt"
    
    # Check if alignment was successful
    if [ -s "\${file_stem}_aligned_sequences.fasta" ]; then
        aligned_count=\$(grep -c "^>" "\${file_stem}_aligned_sequences.fasta" || echo "0")
        echo "Alignment completed successfully with \$aligned_count sequences" >> "\${file_stem}_mafft_alignment_log.txt"
    else
        echo "WARNING: Alignment output file is empty or missing" >> "\${file_stem}_mafft_alignment_log.txt"
        touch "\${file_stem}_aligned_sequences.fasta"
    fi
    """
}
