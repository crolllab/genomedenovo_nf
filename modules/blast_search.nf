// Process to split multi-sequence query file into individual sequences
process SplitQuerySequences {
    tag "Splitting query sequences"
    cpus 4
    memory '2GB' 

    publishDir "${params.outdir}/6_blast_hits/query_sequences", mode: 'copy'
    
    input:
    path blast_query
    
    output:
    path "query_*.fasta", emit: individual_queries
    path "query_info.tsv", emit: query_info
    path "query_type.txt", emit: query_type
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    # Create info file with header
    echo -e "Query_ID\tFilename\tSequence_Length\tSequence_Type" > query_info.tsv
    
    # Function to detect sequence type (nucleotide or protein)
    detect_sequence_type() {
        local fasta_file=\$1
        # Extract sequence (remove header lines)
        seq=\$(grep -v "^>" "\$fasta_file" | tr -d '\n' | tr '[:lower:]' '[:upper:]')
        
        # Count nucleotide characters (A, T, G, C, U, N)
        nuc_count=\$(echo "\$seq" | grep -o '[ATGCUN]' | wc -l)
        total_count=\$(echo "\$seq" | wc -c)
        
        # If >95% are valid nucleotides, consider it nucleotide sequence
        # Otherwise, it's protein
        if [ \$total_count -gt 0 ]; then
            nuc_percent=\$((nuc_count * 100 / total_count))
            if [ \$nuc_percent -ge 95 ]; then
                echo "nucleotide"
            else
                echo "protein"
            fi
        else
            echo "unknown"
        fi
    }
    
    # Split FASTA file into individual sequences
    awk '
    BEGIN { 
        file_count = 0 
        seq_id = ""
        seq = ""
    }
    /^>/ { 
        # If we have a previous sequence, write it out
        if (seq_id != "" && seq != "") {
            # Clean up sequence ID for filename (remove special characters)
            clean_id = seq_id
            gsub(/[^a-zA-Z0-9_-]/, "_", clean_id)
            filename = "query_" file_count "_" clean_id ".fasta"
            print ">" seq_id > filename
            print seq > filename
            close(filename)
        }
        # Start new sequence
        seq_id = substr(\$0, 2)  # Remove ">" character
        seq = ""
        file_count++
    }
    !/^>/ { 
        seq = seq \$0 
    }
    END {
        # Write the last sequence
        if (seq_id != "" && seq != "") {
            clean_id = seq_id
            gsub(/[^a-zA-Z0-9_-]/, "_", clean_id)
            filename = "query_" file_count "_" clean_id ".fasta"
            print ">" seq_id > filename
            print seq > filename
            close(filename)
        }
    }' ${blast_query}
    
    # Check if any sequences were found
    query_count=\$(ls query_*.fasta 2>/dev/null | wc -l || echo "0")
    if [ "\$query_count" -eq 0 ]; then
        echo "ERROR: No valid sequences found in query file"
        exit 1
    fi
    
    echo "Split \$query_count sequences from ${blast_query}"
    
    # Detect sequence type for each query file and update query_info.tsv
    first_seq_type=""
    all_same_type=true
    
    for query_file in query_*.fasta; do
        if [ -f "\$query_file" ]; then
            seq_id=\$(grep "^>" "\$query_file" | sed 's|^>||' | head -n 1)
            seq_len=\$(grep -v "^>" "\$query_file" | tr -d '\n' | wc -c)
            seq_type=\$(detect_sequence_type "\$query_file")
            
            # Track if all sequences are the same type
            if [ -z "\$first_seq_type" ]; then
                first_seq_type="\$seq_type"
            elif [ "\$seq_type" != "\$first_seq_type" ]; then
                all_same_type=false
            fi
            
            echo -e "\${seq_id}\t\${query_file}\t\${seq_len}\t\${seq_type}" >> query_info.tsv
        fi
    done
    
    # Write dominant sequence type to query_type.txt for downstream use
    if [ "\$all_same_type" = true ]; then
        echo "\$first_seq_type" > query_type.txt
        echo "All queries are \$first_seq_type sequences"
    else
        # If mixed, default to nucleotide (safer choice)
        echo "mixed" > query_type.txt
        echo "WARNING: Mixed sequence types detected. Using nucleotide BLAST by default."
    fi
    """
}

// BLAST process definition - supports both blastn (nucleotide) and tblastn (protein)
process BLASTsearch {
    tag "BLAST search"
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    memory '4GB'
    cpus 4
    publishDir { "${params.outdir}/6_blast_hits/${query_id}/blast_results" }, mode: 'copy'
    
    input:
    tuple val(sample_id), path(scaffolds_fasta), val(query_id), path(query_fasta), val(seq_type)
    
    output:
    tuple val(query_id), val(sample_id), path("${sample_id}_blast_results.tsv"), emit: results
    path "${sample_id}_blast_results.tsv"
    
    script:
    // Choose BLAST tool based on sequence type
    def blast_tool = (seq_type == "protein") ? "tblastn" : "blastn"
    // Sanitise sample_id for safe use as a shell filename (replace chars unsafe in unquoted contexts)
    def safe_id = sample_id.replaceAll(/[^A-Za-z0-9._-]/, '_')
    """
    # Create BLAST database from scaffolds
    makeblastdb -in ${scaffolds_fasta} -dbtype nucl -out ${safe_id}_db

    # Run appropriate BLAST search based on query type
    echo "Running ${blast_tool} for ${seq_type} query: ${query_id}"

    ${blast_tool} -query ${query_fasta} \\
           -db ${safe_id}_db \\
           -out ${safe_id}_blast_results.tsv \\
           -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \\
           -num_threads ${task.cpus} \\
           -evalue 1e-5

    # Rename to the real sample_id in the output filename, using a safe temp name
    mv ${safe_id}_blast_results.tsv "${sample_id}_blast_results.tsv"

    # Clean up database files
    rm -f ${safe_id}_db.*
    """
}