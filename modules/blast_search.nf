// Process to split multi-sequence query file into individual sequences
process SplitQuerySequences {
    tag "Splitting query sequences"
    cpus 4
    memory '2GB' 

    publishDir "${params.outdir}/query_sequences", mode: 'copy'
    
    input:
    path blast_query
    
    output:
    path "query_*.fasta", emit: individual_queries
    path "query_info.tsv", emit: query_info
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    # Create info file with header
    echo -e "Query_ID\tFilename\tSequence_Length" > query_info.tsv
    
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
            print seq_id "\t" filename "\t" length(seq) >> "query_info.tsv"
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
            print seq_id "\t" filename "\t" length(seq) >> "query_info.tsv"
        }
    }' ${blast_query}
    
    # Check if any sequences were found
    query_count=\$(ls query_*.fasta 2>/dev/null | wc -l || echo "0")
    if [ "\$query_count" -eq 0 ]; then
        echo "ERROR: No valid sequences found in query file"
        exit 1
    fi
    
    echo "Split \$query_count sequences from ${blast_query}"
    """
}

// BLAST process definition (create this in modules/blast_search.nf)
process BLASTn {
    tag "BLAST search: ${sample_id} vs ${query_id}"
    errorStrategy 'ignore'
    memory '4GB'
    cpus 4
    publishDir "${params.outdir}/${query_id}/blast_results", mode: 'copy'
    
    input:
    tuple val(sample_id), path(scaffolds_fasta), val(query_id), path(query_fasta)
    
    output:
    tuple val(query_id), val(sample_id), path("${sample_id}_blast_results.tsv"), emit: results
    path "${sample_id}_blast_results.tsv"
    
    script:
    """
    # Create BLAST database from scaffolds
    makeblastdb -in ${scaffolds_fasta} -dbtype nucl -out ${sample_id}_db
    
    # Run BLAST search
    blastn -query ${query_fasta} \\
           -db ${sample_id}_db \\
           -out ${sample_id}_blast_results.tsv \\
           -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \\
           -num_threads ${task.cpus} \\
           -evalue 1e-5
    
    # Clean up database files
    rm -f ${sample_id}_db.*
    """
}

// Updated workflow section (add to your main workflow)
workflow denovo_query {
    // ... your existing code ...
    
    // Input validation - add this to your existing input checks
    if (params.blast_query && !file(params.blast_query).exists()) {
        exit 1, "ERROR: BLAST query file does not exist: ${params.blast_query}"
    }
    
    // ... your existing SRA and SPAdes code ...
    
    // Run SPAdes
    SPAdes(combined_pe_ch)
    
    // Run BLAST if query file is provided
    if (params.blast_query) {
        blast_query_file = file(params.blast_query)
        BLASTn(SPAdes.out.scaffolds, blast_query_file)
    }
}