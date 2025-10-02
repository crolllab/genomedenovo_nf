// Process to extract top BLAST hit regions from assemblies
process ExtractBlastHits {
    tag "Extract BLAST hits: ${sample_id}"
    errorStrategy 'ignore'
    memory '2GB'
    cpus 1
    publishDir "${params.outdir}/blast_hit_sequences", mode: 'copy'
    
    input:
    tuple val(sample_id), path(blast_results), path(scaffolds_fasta)
    
    output:
    tuple val(sample_id), path("${sample_id}_top_hits.fasta"), emit: sequences
    path "${sample_id}_top_hits_info.tsv", emit: hit_info
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    # Check if BLAST results file exists and is not empty
    if [ ! -s "${blast_results}" ]; then
        echo "No BLAST hits found for ${sample_id}" > ${sample_id}_top_hits_info.tsv
        touch ${sample_id}_top_hits.fasta
        exit 0
    fi
    
    # Create header for hit info file
    echo -e "Sample_ID\tQuery_ID\tSubject_ID\tPercent_Identity\tAlignment_Length\tStart\tEnd\tE_value\tBit_Score" > ${sample_id}_top_hits_info.tsv
    
    # Extract top hit per query (highest bit score)
    awk -F'\t' '
    {
        # BLAST format: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
        query = \$1
        bitscore = \$12
        
        # Keep track of best hit per query (highest bitscore)
        if (!(query in best_bitscore) || bitscore > best_bitscore[query]) {
            best_bitscore[query] = bitscore
            best_hit[query] = \$0
        }
    }
    END {
        for (query in best_hit) {
            print best_hit[query]
        }
    }' "${blast_results}" > top_hits.tsv
    
    # Check if we have any hits
    if [ ! -s top_hits.tsv ]; then
        echo "No BLAST hits passed filtering for ${sample_id}" > ${sample_id}_top_hits_info.tsv
        touch ${sample_id}_top_hits.fasta
        exit 0
    fi
    
    # Create FASTA index
    samtools faidx "${scaffolds_fasta}"
    
    # Process each top hit
    > ${sample_id}_top_hits.fasta
    
    while IFS=\$'\t' read -r qseqid sseqid pident align_length mismatch gapopen qstart qend sstart send evalue bitscore; do
        # Determine start and end coordinates (ensure start < end for samtools)
        if [ "\$sstart" -le "\$send" ]; then
            start=\$sstart
            end=\$send
        else
            start=\$send
            end=\$sstart
        fi
        
        # Extract sequence using samtools faidx with proper coordinate format
        if samtools faidx "${scaffolds_fasta}" "\${sseqid}:\${start}-\${end}" > temp_hit.fasta 2>/dev/null; then
            # Rename the header with informative name
            sed "s|^>.*|>${sample_id}_\${qseqid}_hit_\${sseqid}:\${start}-\${end}|" temp_hit.fasta >> ${sample_id}_top_hits.fasta
            
            # Add hit information to TSV
            echo -e "${sample_id}\t\${qseqid}\t\${sseqid}\t\${pident}\t\${align_length}\t\${start}\t\${end}\t\${evalue}\t\${bitscore}" >> ${sample_id}_top_hits_info.tsv
        else
            # Fallback: Try extracting full sequence and manually cutting
            if samtools faidx "${scaffolds_fasta}" "\${sseqid}" 2>/dev/null | tail -n +2 | tr -d '\n' > temp_full_seq.txt; then
                # Get length of sequence
                seqlen=\$(wc -c < temp_full_seq.txt)
                
                # Validate coordinates
                if [ "\$end" -le "\$seqlen" ] && [ "\$start" -gt 0 ]; then
                    # Extract substring (1-indexed for cut)
                    cat temp_full_seq.txt | cut -c\${start}-\${end} > temp_subseq.txt
                    
                    if [ -s temp_subseq.txt ]; then
                        echo ">${sample_id}_\${qseqid}_hit_\${sseqid}:\${start}-\${end}" >> ${sample_id}_top_hits.fasta
                        # Add newlines every 60 characters for proper FASTA format
                        fold -w 60 temp_subseq.txt >> ${sample_id}_top_hits.fasta
                        
                        echo -e "${sample_id}\t\${qseqid}\t\${sseqid}\t\${pident}\t\${align_length}\t\${start}\t\${end}\t\${evalue}\t\${bitscore}" >> ${sample_id}_top_hits_info.tsv
                    fi
                fi
            fi
        fi
        
        rm -f temp_hit.fasta temp_full_seq.txt temp_subseq.txt
        
    done < top_hits.tsv
    
    # Clean up
    rm -f top_hits.tsv
    
    # If no sequences were extracted, ensure files exist
    if [ ! -s ${sample_id}_top_hits.fasta ]; then
        touch ${sample_id}_top_hits.fasta
    fi
    """
}

// Process to combine all extracted hit sequences
process CombineBlastHits {
    tag "Combining all BLAST hit sequences"
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path hit_sequences
    path hit_info_files
    
    output:
    path "all_blast_hits.fasta"
    path "all_blast_hits_summary.tsv"
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    # Combine all FASTA files (skip empty files)
    for fasta in *.fasta; do
        if [ -s "\$fasta" ]; then
            cat "\$fasta"
        fi
    done > all_blast_hits.fasta
    
    # Combine all TSV files (keep header from first file only)
    first_file=\$(ls *_top_hits_info.tsv 2>/dev/null | head -n 1 || echo "")
    if [ -n "\$first_file" ]; then
        head -n 1 "\$first_file" > all_blast_hits_summary.tsv
        for file in *_top_hits_info.tsv; do
            tail -n +2 "\$file" >> all_blast_hits_summary.tsv
        done
    else
        echo -e "Sample_ID\tQuery_ID\tSubject_ID\tPercent_Identity\tAlignment_Length\tStart\tEnd\tE_value\tBit_Score" > all_blast_hits_summary.tsv
        echo "No BLAST hit files found" >> all_blast_hits_summary.tsv
    fi
    """
}