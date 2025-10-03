// Process to extract top BLAST hit regions from assemblies
process ExtractBlastHits {
    tag "Extract BLAST hits for ${query_id}"
    errorStrategy 'ignore'
    memory '2GB'
    cpus 1
    publishDir "${params.outdir}/${query_id}/blast_hit_sequences", mode: 'copy'
    
    input:
    tuple val(query_id), val(sample_id), path(blast_results), path(scaffolds_fasta)
    
    output:
    tuple val(query_id), val(sample_id), path("${sample_id}_top_hits.fasta"), emit: sequences
    tuple val(query_id), path("${sample_id}_top_hits_info.tsv"), emit: hit_info
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    # Create header for hit info file
    echo -e "Sample_ID\tQuery_ID\tSubject_ID\tPercent_Identity\tAlignment_Length\tStart\tEnd\tE_value\tBit_Score" > ${sample_id}_top_hits_info.tsv
    
    # Check if BLAST results file exists and is not empty
    if [ ! -s "${blast_results}" ]; then
        touch ${sample_id}_top_hits.fasta
        exit 0
    fi
    
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
        touch ${sample_id}_top_hits.fasta
        exit 0
    fi
    
    # Create FASTA index
    samtools faidx "${scaffolds_fasta}"
    
    # Process each top hit
    > ${sample_id}_top_hits.fasta
    
    while IFS=\$'\t' read -r qseqid sseqid pident align_length mismatch gapopen qstart qend sstart send evalue bitscore; do
        # Determine start and end coordinates (ensure start < end for samtools)
        # Also determine if sequence is in reverse orientation
        if [ "\$sstart" -le "\$send" ]; then
            start=\$sstart
            end=\$send
            reverse=false
        else
            start=\$send
            end=\$sstart
            reverse=true
        fi
        
        # Extract sequence using samtools faidx with proper coordinate format
        if samtools faidx "${scaffolds_fasta}" "\${sseqid}:\${start}-\${end}" > temp_hit.fasta 2>/dev/null; then
            # Extract just the sequence (no header)
            tail -n +2 temp_hit.fasta | tr -d '\n' > temp_seq.txt
            
            # Reverse complement if hit is on the reverse strand
            if [ "\$reverse" = true ]; then
                # Reverse complement using tr and rev
                cat temp_seq.txt | tr 'ATGCatgc' 'TACGtacg' | rev > temp_rc.txt
                mv temp_rc.txt temp_seq.txt
            fi
            
            # Write sequence with informative header
            echo ">${sample_id}_\${qseqid}_hit_\${sseqid}:\${start}-\${end}" >> ${sample_id}_top_hits.fasta
            fold -w 60 temp_seq.txt >> ${sample_id}_top_hits.fasta
            echo "" >> ${sample_id}_top_hits.fasta
            
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
                        # Reverse complement if hit is on the reverse strand
                        if [ "\$reverse" = true ]; then
                            cat temp_subseq.txt | tr 'ATGCatgc' 'TACGtacg' | rev > temp_rc.txt
                            mv temp_rc.txt temp_subseq.txt
                        fi
                        
                        echo ">${sample_id}_\${qseqid}_hit_\${sseqid}:\${start}-\${end}" >> ${sample_id}_top_hits.fasta
                        # Add newlines every 60 characters for proper FASTA format
                        fold -w 60 temp_subseq.txt >> ${sample_id}_top_hits.fasta
                        echo "" >> ${sample_id}_top_hits.fasta
                        
                        echo -e "${sample_id}\t\${qseqid}\t\${sseqid}\t\${pident}\t\${align_length}\t\${start}\t\${end}\t\${evalue}\t\${bitscore}" >> ${sample_id}_top_hits_info.tsv
                    fi
                fi
            fi
        fi
        
        rm -f temp_hit.fasta temp_seq.txt temp_rc.txt temp_full_seq.txt temp_subseq.txt
        
    done < top_hits.tsv
    
    # Clean up
    rm -f top_hits.tsv
    
    # If no sequences were extracted, ensure files exist
    if [ ! -s ${sample_id}_top_hits.fasta ]; then
        touch ${sample_id}_top_hits.fasta
    fi
    """
}

// Process to combine all extracted hit sequences for a specific query
process CombineBlastHits {
    tag "Combining BLAST hits for ${query_id}"
    publishDir "${params.outdir}/${query_id}", mode: 'copy'
    
    input:
    tuple val(query_id), path(hit_sequences), path(hit_info_files)
    
    output:
    tuple val(query_id), path("${query_id}_all_blast_hits.fasta"), emit: combined_sequences
    tuple val(query_id), path("${query_id}_all_blast_hits_summary.tsv"), emit: summary
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    # Combine all FASTA files (skip empty files)
    # Use specific pattern to avoid reading the output file
    for fasta in *_top_hits.fasta; do
        if [ -f "\$fasta" ] && [ -s "\$fasta" ]; then
            cat "\$fasta"
        fi
    done > ${query_id}_all_blast_hits.fasta
    
    # Combine all TSV files and add NO HIT entries for samples without hits
    first_file=\$(ls *_top_hits_info.tsv 2>/dev/null | head -n 1 || echo "")
    if [ -n "\$first_file" ]; then
        # Create header
        head -n 1 "\$first_file" > ${query_id}_all_blast_hits_summary.tsv
        
        # Process each TSV file
        for file in *_top_hits_info.tsv; do
            # Extract sample ID from filename
            sample_id=\$(basename "\$file" "_top_hits_info.tsv")
            
            # Check if file has data beyond the header
            line_count=\$(wc -l < "\$file")
            if [ "\$line_count" -gt 1 ]; then
                # File has hits, add them
                tail -n +2 "\$file" >> ${query_id}_all_blast_hits_summary.tsv
            else
                # File has only header (no hits), add no_hits_found entry
                echo -e "\${sample_id}\tno_hits_found\tno_hits_found\tNA\tNA\tNA\tNA\tNA\tNA" >> ${query_id}_all_blast_hits_summary.tsv
            fi
        done
    else
        echo -e "Sample_ID\tQuery_ID\tSubject_ID\tPercent_Identity\tAlignment_Length\tStart\tEnd\tE_value\tBit_Score" > ${query_id}_all_blast_hits_summary.tsv
        echo "No BLAST hit files found" >> ${query_id}_all_blast_hits_summary.tsv
    fi
    """
}