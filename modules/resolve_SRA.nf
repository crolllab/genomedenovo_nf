process SRAresolve {
    tag "Querying $accessions_file for SRR ids"
    cpus 1
    memory '4GB'
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path accessions_file
    
    output:
    path "NCBI_SRR_PE_accessions.txt", emit: pe_file
    path "NCBI_SRR_SE_accessions.txt", emit: se_file

    script:
    """
    #!/bin/bash
    
    # Clear output files
    > NCBI_SRR_PE_accessions.txt
    > NCBI_SRR_SE_accessions.txt

    # Read all accessions into array
    accessions=( \$(grep -v '^\\s*\$' ${accessions_file}) )
    
    for accession in "\${accessions[@]}"; do
        echo "Processing \$accession"
        
        # Fetch runinfo and process each line
        esearch -db sra -query "\$accession" \
        | efetch -format runinfo 2>/dev/null \
        | tail -n +2 \
        | while IFS= read -r line; do
            # Extract run ID (first column) and layout (16th column)
            run_id=\$(echo "\$line" | cut -d',' -f1)
            layout=\$(echo "\$line" | cut -d',' -f16)
            
            echo "Debug: Processing \$run_id with layout '\$layout'" >&2
            
            # Check if layout field contains "PAIRED" or "SINGLE"
            if [[ "\$layout" == "PAIRED" ]]; then
                echo "\$run_id" >> NCBI_SRR_PE_accessions.txt
                echo "Added \$run_id to PE file" >&2
            elif [[ "\$layout" == "SINGLE" ]]; then
                echo "\$run_id" >> NCBI_SRR_SE_accessions.txt
                echo "Added \$run_id to SE file" >&2
            else
                echo "Warning: Unknown layout '\$layout' for \$run_id" >&2
            fi
        done
    done
    
    # Report results
    pe_count=\$(wc -l < NCBI_SRR_PE_accessions.txt)
    se_count=\$(wc -l < NCBI_SRR_SE_accessions.txt)

    echo "Done! Found \$pe_count paired-end and \$se_count single-end SRR accessions"
    echo "PE accessions saved in NCBI_SRR_PE_accessions.txt"
    echo "SE accessions saved in NCBI_SRR_SE_accessions.txt"
    """
}