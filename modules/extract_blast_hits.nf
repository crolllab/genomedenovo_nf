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
    
    # Padding value (bp to extend on both sides)
    PADDING=${params.padding}
    
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
        
        # Store original coordinates for TSV output
        orig_start=\$start
        orig_end=\$end
        
        # Apply padding - extend coordinates by PADDING bp on both sides
        # Get scaffold length from FASTA index
        scaffold_len=\$(awk -v seq="\$sseqid" '\$1==seq {print \$2}' "${scaffolds_fasta}.fai")
        
        # Sanity check: Ensure scaffold length was found
        if [ -z "\$scaffold_len" ] || [ "\$scaffold_len" -lt 1 ]; then
            echo "WARNING: Could not determine length for scaffold \$sseqid, skipping..." >&2
            continue
        fi
        
        # Sanity check: Ensure original coordinates are valid
        if [ "\$start" -lt 1 ] || [ "\$end" -gt "\$scaffold_len" ] || [ "\$start" -gt "\$end" ]; then
            echo "WARNING: Invalid coordinates for \$sseqid (\$start-\$end), scaffold length is \$scaffold_len, skipping..." >&2
            continue
        fi
        
        # Calculate padded coordinates
        padded_start=\$(( start - PADDING ))
        padded_end=\$(( end + PADDING ))
        
        # Ensure padded coordinates are within scaffold boundaries
        if [ "\$padded_start" -lt 1 ]; then
            padded_start=1
        fi
        
        if [ "\$padded_end" -gt "\$scaffold_len" ]; then
            padded_end=\$scaffold_len
        fi
        
        # Final sanity check: Ensure padded coordinates are valid
        if [ "\$padded_start" -lt 1 ] || [ "\$padded_end" -gt "\$scaffold_len" ] || [ "\$padded_start" -gt "\$padded_end" ]; then
            echo "WARNING: Invalid padded coordinates for \$sseqid (\$padded_start-\$padded_end), scaffold length is \$scaffold_len, skipping..." >&2
            continue
        fi
        
        # Use padded coordinates for extraction
        start=\$padded_start
        end=\$padded_end
        
        # Log extraction details (useful for debugging)
        if [ \$PADDING -gt 0 ]; then
            echo "Extracting \$sseqid:\$start-\$end (original: \$orig_start-\$orig_end, padding: \$PADDING bp, scaffold_len: \$scaffold_len)" >&2
        fi
        
        # Extract sequence using samtools faidx with proper coordinate format
        if samtools faidx "${scaffolds_fasta}" "\${sseqid}:\${start}-\${end}" > temp_hit.fasta 2>temp_samtools_error.log; then
            # Extract just the sequence (no header)
            tail -n +2 temp_hit.fasta | tr -d '\n' > temp_seq.txt
            
            # Reverse complement if hit is on the reverse strand
            if [ "\$reverse" = true ]; then
                # Reverse complement using tr and rev
                cat temp_seq.txt | tr 'ATGCatgc' 'TACGtacg' | rev > temp_rc.txt
                mv temp_rc.txt temp_seq.txt
            fi
            
            # Write sequence with informative header (show padded coordinates)
            echo ">${sample_id}_\${qseqid}_hit_\${sseqid}:\${start}-\${end}" >> ${sample_id}_top_hits.fasta
            fold -w 60 temp_seq.txt >> ${sample_id}_top_hits.fasta
            echo "" >> ${sample_id}_top_hits.fasta
            
            # Add hit information to TSV (use original BLAST hit coordinates, not padded)
            echo -e "${sample_id}\t\${qseqid}\t\${sseqid}\t\${pident}\t\${align_length}\t\${orig_start}\t\${orig_end}\t\${evalue}\t\${bitscore}" >> ${sample_id}_top_hits_info.tsv
            rm temp_hit.fasta temp_seq.txt temp_samtools_error.log 2>/dev/null || true
        else
            # Log samtools error details
            echo "WARNING: samtools faidx failed for \${sseqid}:\${start}-\${end} (scaffold_len: \${scaffold_len}), trying fallback method..." >&2
            if [ -f temp_samtools_error.log ] && [ -s temp_samtools_error.log ]; then
                echo "  samtools error: \$(cat temp_samtools_error.log)" >&2
            fi
            rm temp_samtools_error.log 2>/dev/null || true
            
            # Fallback: Try extracting full sequence and manually cutting
            if samtools faidx "${scaffolds_fasta}" "\${sseqid}" 2>/dev/null | tail -n +2 | tr -d '\n' > temp_full_seq.txt; then
                # Get length of sequence
                seqlen=\$(wc -c < temp_full_seq.txt)
                
                # Validate padded coordinates
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
                        
                        # Use original coordinates in TSV
                        echo -e "${sample_id}\t\${qseqid}\t\${sseqid}\t\${pident}\t\${align_length}\t\${orig_start}\t\${orig_end}\t\${evalue}\t\${bitscore}" >> ${sample_id}_top_hits_info.tsv
                    else
                        echo "ERROR: Fallback method failed - extracted sequence is empty for \${sseqid}:\${start}-\${end}" >&2
                    fi
                else
                    echo "ERROR: Fallback method failed - coordinates out of bounds for \${sseqid} (start:\$start end:\$end seqlen:\$seqlen)" >&2
                fi
            else
                echo "ERROR: Fallback method failed - could not extract full sequence for \${sseqid}" >&2
            fi
        fi
        
        # Clean up all temporary files
        rm -f temp_hit.fasta temp_seq.txt temp_rc.txt temp_full_seq.txt temp_subseq.txt temp_samtools_error.log
        
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
    tuple val(query_id), path(query_fasta), path(hit_sequences), path(hit_info_files)
    
    output:
    tuple val(query_id), path("${query_id}_all_blast_hits.fasta"), emit: combined_sequences
    tuple val(query_id), path("${query_id}_all_blast_hits_summary.tsv"), emit: summary
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    # First, add the query sequence to the combined FASTA
    if [ -f "${query_fasta}" ]; then
        cat "${query_fasta}" > ${query_id}_all_blast_hits.fasta
        echo "" >> ${query_id}_all_blast_hits.fasta
    else
        touch ${query_id}_all_blast_hits.fasta
    fi
    
    # Then combine all BLAST hit FASTA files (skip empty files)
    # Use specific pattern to avoid reading the output file
    for fasta in *_top_hits.fasta; do
        if [ -f "\$fasta" ] && [ -s "\$fasta" ]; then
            cat "\$fasta" >> ${query_id}_all_blast_hits.fasta
        fi
    done
    
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