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
    
    echo "DEBUG: Processing ${sample_id}"
    echo "DEBUG: BLAST results file: ${blast_results}"
    echo "DEBUG: Scaffolds file: ${scaffolds_fasta}"
    
    # Check if BLAST results file exists and is not empty
    if [ ! -s "${blast_results}" ]; then
        echo "No BLAST hits found for ${sample_id}" > ${sample_id}_top_hits_info.tsv
        touch ${sample_id}_top_hits.fasta
        exit 0
    fi
    
    echo "DEBUG: BLAST file size: \$(wc -l < ${blast_results}) lines"
    echo "DEBUG: First few lines of BLAST results:"
    head -5 ${blast_results}
    
    # Create header for hit info file
    echo -e "Sample_ID\\tQuery_ID\\tSubject_ID\\tPercent_Identity\\tAlignment_Length\\tStart\\tEnd\\tE_value\\tBit_Score" > ${sample_id}_top_hits_info.tsv
    
    # Extract top hit per query (highest bit score)
    awk -F'\\t' '
    {
        # BLAST format: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
        query = \$1
        subject = \$2
        pident = \$3
        align_length = \$4
        sstart = \$9
        send = \$10
        evalue = \$11
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
    }' ${blast_results} > top_hits.tsv
    
    echo "DEBUG: Top hits found: \$(wc -l < top_hits.tsv) hits"
    echo "DEBUG: Top hits content:"
    cat top_hits.tsv
    
    echo "DEBUG: Available sequences in FASTA:"
    grep "^>" ${scaffolds_fasta} | head -10
    
    # Create FASTA index if it doesn't exist
    echo "DEBUG: Creating FASTA index..."
    if [ ! -f "${scaffolds_fasta}.fai" ]; then
        samtools faidx ${scaffolds_fasta}
    fi
    
    # Process each top hit
    > ${sample_id}_top_hits.fasta
    
    while IFS=\$'\\t' read -r qseqid sseqid pident align_length mismatch gapopen qstart qend sstart send evalue bitscore; do
        echo "DEBUG: Processing hit - Query: \$qseqid, Subject: \$sseqid, Coords: \$sstart-\$send"
        
        # Determine start and end coordinates (ensure start < end)
        if [ "\$sstart" -le "\$send" ]; then
            start=\$sstart
            end=\$send
        else
            start=\$send
            end=\$sstart
        fi
        
        echo "DEBUG: Adjusted coordinates: \$start-\$end"
        
        # Check if subject sequence exists in FASTA
        if grep -q ">\$sseqid" ${scaffolds_fasta}; then
            echo "DEBUG: Found subject sequence \$sseqid in FASTA"
        else
            echo "DEBUG: Subject sequence \$sseqid NOT found in FASTA"
            echo "DEBUG: Trying alternative matching..."
            grep "^>" ${scaffolds_fasta} | grep "\$sseqid" || echo "No matches found"
        fi
        
        # Extract sequence using samtools faidx
        echo "DEBUG: Running samtools faidx command..."
        samtools faidx ${scaffolds_fasta} "\$sseqid:\${start}-\${end}" > temp_hit.fasta 2>samtools_error.log
        
        if [ -s temp_hit.fasta ]; then
            echo "DEBUG: Successfully extracted sequence"
            # Rename the extracted sequence with informative header
            sed "s/>.*/>${sample_id}_\${qseqid}_hit_\${sseqid}:\${start}-\${end}/" temp_hit.fasta >> ${sample_id}_top_hits.fasta
            
            # Add hit information to TSV
            echo -e "${sample_id}\\t\${qseqid}\\t\${sseqid}\\t\${pident}\\t\${align_length}\\t\${start}\\t\${end}\\t\${evalue}\\t\${bitscore}" >> ${sample_id}_top_hits_info.tsv
        else
            echo "DEBUG: Failed to extract sequence with coordinates"
            echo "DEBUG: samtools error log:"
            cat samtools_error.log || echo "No error log"
            
            # Try extracting the full sequence and then cutting it
            echo "DEBUG: Trying full sequence extraction..."
            samtools faidx ${scaffolds_fasta} "\$sseqid" > temp_full.fasta 2>samtools_error2.log
            
            if [ -s temp_full.fasta ]; then
                echo "DEBUG: Full sequence extracted, now extracting region..."
                # Extract just the sequence (no header), then get the substring
                tail -n +2 temp_full.fasta | tr -d '\\n' | cut -c\${start}-\${end} > temp_seq.txt
                
                if [ -s temp_seq.txt ]; then
                    echo "DEBUG: Region extraction succeeded"
                    echo ">${sample_id}_\${qseqid}_hit_\${sseqid}:\${start}-\${end}" > temp_hit.fasta
                    cat temp_seq.txt >> temp_hit.fasta
                    cat temp_hit.fasta >> ${sample_id}_top_hits.fasta
                    echo -e "${sample_id}\\t\${qseqid}\\t\${sseqid}\\t\${pident}\\t\${align_length}\\t\${start}\\t\${end}\\t\${evalue}\\t\${bitscore}" >> ${sample_id}_top_hits_info.tsv
                else
                    echo "DEBUG: Region extraction failed"
                fi
            else
                echo "DEBUG: Full sequence extraction also failed"
                cat samtools_error2.log || echo "No error log"
            fi
        fi
        
        rm -f temp_hit.fasta temp_full.fasta temp_seq.txt samtools_error.log samtools_error2.log
        
    done < top_hits.tsv
    
    # Clean up
    rm -f top_hits.tsv
    
    # If no sequences were extracted, create empty files
    if [ ! -s ${sample_id}_top_hits.fasta ]; then
        echo "No sequences could be extracted for ${sample_id}" >> ${sample_id}_top_hits_info.tsv
        touch ${sample_id}_top_hits.fasta
    fi
    
    echo "DEBUG: Final output file size: \$(wc -c < ${sample_id}_top_hits.fasta) bytes"
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
    # Combine all FASTA files
    cat *.fasta > all_blast_hits.fasta
    
    # Combine all TSV files (keep header from first file only)
    first_file=\$(ls *_top_hits_info.tsv | head -n 1)
    if [ -n "\$first_file" ]; then
        head -n 1 "\$first_file" > all_blast_hits_summary.tsv
        for file in *_top_hits_info.tsv; do
            tail -n +2 "\$file" >> all_blast_hits_summary.tsv
        done
    else
        echo -e "Sample_ID\\tQuery_ID\\tSubject_ID\\tPercent_Identity\\tAlignment_Length\\tStart\\tEnd\\tE_value\\tBit_Score" > all_blast_hits_summary.tsv
        echo "No BLAST hit files found" >> all_blast_hits_summary.tsv
    fi
    """
}