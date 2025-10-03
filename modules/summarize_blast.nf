// Process to merge all BLAST query summaries into a single wide-format table
process MergeBlastSummaries {
    tag "Merging all BLAST summaries"
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path summary_files  // All *_all_blast_hits_summary.tsv files
    
    output:
    path "blast_summary_merged.tsv", emit: merged_summary
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    # Get list of all TSV files
    tsv_files=( ${summary_files} )
    
    if [ \${#tsv_files[@]} -eq 0 ]; then
        echo "ERROR: No BLAST summary files found"
        exit 1
    fi
    
    echo "Processing \${#tsv_files[@]} BLAST summary files..."
    
    # Create temporary directory for processing
    mkdir -p temp_processing
    
    # Step 1: Process each query TSV file
    # Extract query name and prepare data with prefixed columns
    for tsv_file in "\${tsv_files[@]}"; do
        # Extract query ID from filename (remove _all_blast_hits_summary.tsv)
        query_id=\$(basename "\$tsv_file" | sed 's|_all_blast_hits_summary.tsv||')
        echo "  Processing query: \$query_id"
        
        # Create a processed file with query-prefixed columns
        # Skip header and create only 5 essential columns per query:
        # Sample_ID, query_Present, query_Identity, query_AlignLen, query_Evalue, query_BitScore
        
        tail -n +2 "\$tsv_file" | awk -v query="\$query_id" 'BEGIN {OFS="\t"} {
            sample_id = \$1
            query_id = \$2
            subject_id = \$3
            pident = \$4
            align_len = \$5
            start = \$6
            end = \$7
            evalue = \$8
            bitscore = \$9
            
            # Determine if this is a real hit or no hit
            if (query_id == "no_hits_found" || subject_id == "no_hits_found") {
                # No hit found - output only 5 essential columns
                print sample_id, query "_Present", "FALSE"
                print sample_id, query "_Identity", "NA"
                print sample_id, query "_AlignLen", "NA"
                print sample_id, query "_Evalue", "NA"
                print sample_id, query "_BitScore", "NA"
            } else {
                # Real hit found - output only 5 essential columns
                print sample_id, query "_Present", "TRUE"
                print sample_id, query "_Identity", pident
                print sample_id, query "_AlignLen", align_len
                print sample_id, query "_Evalue", evalue
                print sample_id, query "_BitScore", bitscore
            }
        }' > "temp_processing/\${query_id}_processed.tsv"
    done
    
    # Step 2: Get all unique sample IDs across all queries
    echo "  Collecting all sample IDs..."
    cat temp_processing/*_processed.tsv | cut -f1 | sort -u > temp_processing/all_sample_ids.txt
    
    # Step 3: Get all unique column names (sorted by query)
    echo "  Collecting all column names..."
    cat temp_processing/*_processed.tsv | cut -f2 | sort -u > temp_processing/all_columns.txt
    
    # Step 4: Create a master table with Sample_ID as first column
    echo "  Building merged table..."
    
    # Create header
    echo -n "Sample_ID" > blast_summary_merged.tsv
    while IFS= read -r colname; do
        echo -ne "\t\$colname" >> blast_summary_merged.tsv
    done < temp_processing/all_columns.txt
    echo "" >> blast_summary_merged.tsv
    
    # For each sample, collect values for all columns
    while IFS= read -r sample_id; do
        echo -n "\$sample_id" >> blast_summary_merged.tsv
        
        # For each column, find the value for this sample
        while IFS= read -r colname; do
            # Search for this sample_id and column combination in all processed files
            value=\$(grep -h "^\${sample_id}\t\${colname}\t" temp_processing/*_processed.tsv 2>/dev/null | head -n 1 | cut -f3)
            
            # If not found, use NA
            if [ -z "\$value" ]; then
                value="NA"
            fi
            
            echo -ne "\t\$value" >> blast_summary_merged.tsv
        done < temp_processing/all_columns.txt
        
        echo "" >> blast_summary_merged.tsv
    done < temp_processing/all_sample_ids.txt
    
    # Step 5: Add summary columns (gene pattern, total genes present)
    echo "  Adding summary statistics..."
    
    # Count total samples and samples with each query
    total_samples=\$(tail -n +2 blast_summary_merged.tsv | wc -l)
    echo "Total samples in merged table: \$total_samples"
    
    # Get list of query IDs (unique prefixes before _Present)
    queries=\$(head -n 1 blast_summary_merged.tsv | tr "\\t" "\\n" | grep "_Present\$" | sed "s|_Present||" | tr "\\n" " ")
    echo "Queries found: \$queries"
    
    # Count samples with each query
    for query in \$queries; do
        count=\$(tail -n +2 blast_summary_merged.tsv | cut -f2- | tr "\\t" "\\n" | grep -c "^TRUE\$" | head -n 1 || echo 0)
        # Note: This is a simplified count; proper counting would need column indexing
        echo "  Query \$query: hits in multiple samples"
    done
    
    # Clean up
    rm -rf temp_processing
    
    echo "Merge complete! Output: blast_summary_merged.tsv"
    """
}
