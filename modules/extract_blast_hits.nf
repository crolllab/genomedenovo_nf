// Process to extract top BLAST hit regions from assemblies
process ExtractBlastHits {
    tag "Extract BLAST hits: ${sample_id}"
    errorStrategy 'ignore'
    memory '2GB'
    cpus 1
    publishDir "${params.outdir}/blast_hit_sequences", mode: 'copy'
    
    input:
    tuple val(sample_id), path(blast_results)
    tuple val(sample_id_scaffolds), path(scaffolds_fasta)
    
    output:
    tuple val(sample_id), path("${sample_id}_top_hits.fasta"), emit: sequences
    path "${sample_id}_top_hits_info.tsv", emit: hit_info
    
    when:
    sample_id == sample_id_scaffolds
    
    script:
    """
    #!/bin/bash
    
    # Check if BLAST results file exists and is not empty
    if [ ! -s "${blast_results}" ]; then
        echo "No BLAST hits found for ${sample_id}" > ${sample_id}_top_hits_info.tsv
        touch ${sample_id}_top_hits.fasta
        exit 0
    fi
    
    # Create header for hit info file
    echo -e "Sample_ID\\tQuery_ID\\tSubject_ID\\tPercent_Identity\\tAlignment_Length\\tStart\\tEnd\\tE_value\\tBit_Score" > ${sample_id}_top_hits_info.tsv
    
    # Extract top hit per query (highest bit score)
    awk -F'\\t' '
    {
        # BLAST format: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
        query = \$1
        subject = \$2
        pident = \$3
        length = \$4
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
    
    # Process each top hit
    > ${sample_id}_top_hits.fasta
    
    while IFS=\$'\\t' read -r qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore; do
        # Determine start and end coordinates (ensure start < end)
        if [ "\$sstart" -le "\$send" ]; then
            start=\$sstart
            end=\$send
        else
            start=\$send
            end=\$sstart
        fi
        
        # Extract sequence using seqkit
        seqkit subseq -r \${start}:\${end} ${scaffolds_fasta} -n "\$sseqid" >> temp_hit.fasta
        
        # Rename the extracted sequence with informative header
        if [ -s temp_hit.fasta ]; then
            sed "s/>.*/>${sample_id}_\${qseqid}_hit_\${sseqid}:\${start}-\${end}/" temp_hit.fasta >> ${sample_id}_top_hits.fasta
            
            # Add hit information to TSV
            echo -e "${sample_id}\\t\${qseqid}\\t\${sseqid}\\t\${pident}\\t\${length}\\t\${start}\\t\${end}\\t\${evalue}\\t\${bitscore}" >> ${sample_id}_top_hits_info.tsv
        fi
        
        rm -f temp_hit.fasta
        
    done < top_hits.tsv
    
    # Clean up
    rm -f top_hits.tsv
    
    # If no sequences were extracted, create empty files
    if [ ! -s ${sample_id}_top_hits.fasta ]; then
        echo "No sequences could be extracted for ${sample_id}" >> ${sample_id}_top_hits_info.tsv
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
