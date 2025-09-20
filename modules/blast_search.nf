// BLAST process definition (create this in modules/blast_search.nf)
process BLASTn {
    tag "BLAST search: ${sample_id}"
    errorStrategy 'ignore'
    memory '4GB'
    cpus 4
    publishDir "${params.outdir}/blast_results", mode: 'copy'
    
    input:
    tuple val(sample_id), path(scaffolds_fasta)
    path blast_query
    
    output:
    tuple val(sample_id), path("${sample_id}_blast_results.tsv"), emit: results
    path "${sample_id}_blast_results.tsv"
    
    script:
    """
    # Create BLAST database from scaffolds
    makeblastdb -in ${scaffolds_fasta} -dbtype nucl -out ${sample_id}_db
    
    # Run BLAST search
    blastn -query ${blast_query} \\
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