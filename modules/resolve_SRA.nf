process SRAresolve {
    tag "Querying $accessions_file for SRR ids and ENA URLs"
    cpus 1
    memory '4GB'
    publishDir "${params.outdir}", mode: 'copy'
    errorStrategy 'ignore'

    input:
    path accessions_file
    
    output:
    path "NCBI_SRR_PE_accessions.txt", emit: pe_file
    path "NCBI_SRR_SE_accessions.txt", emit: se_file
    path "NCBI_download_urls.tsv", emit: url_file

    script:
    """
    #!/bin/bash
    
    # Clear output files
    > NCBI_SRR_PE_accessions.txt
    > NCBI_SRR_SE_accessions.txt
    > NCBI_download_urls.tsv
    
    # Header for URL file
    echo -e "SRR_ID\tLayout\tURL_1\tURL_2\tSource" > NCBI_download_urls.tsv

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
            
            echo "Processing \$run_id with layout '\$layout'" >&2
            
            # Try to get ENA FTP URLs (provides direct FASTQ links)
            ena_response=\$(wget -qO- "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=\$run_id&result=read_run&fields=fastq_ftp" | tail -n +2)
            
            url1=""
            url2=""
            source="NCBI"
            
            if [[ -n "\$ena_response" ]]; then
                # ENA provides semicolon-separated URLs for PE reads
                ena_urls=\$(echo "\$ena_response" | cut -f2)
                
                if [[ -n "\$ena_urls" && "\$ena_urls" != "fastq_ftp" ]]; then
                    url1=\$(echo "\$ena_urls" | cut -d';' -f1)
                    url2=\$(echo "\$ena_urls" | cut -d';' -f2)
                    
                    # Convert FTP to HTTP for better reliability
                    if [[ -n "\$url1" ]]; then
                        url1="http://\${url1#ftp://}"
                        source="ENA"
                        echo "Found ENA URL for \$run_id" >&2
                    fi
                    
                    if [[ -n "\$url2" ]]; then
                        url2="http://\${url2#ftp://}"
                    fi
                fi
            fi
            
            # Classify by layout and save
            if [[ "\$layout" == "PAIRED" ]]; then
                echo "\$run_id" >> NCBI_SRR_PE_accessions.txt
                echo -e "\$run_id\tPE\t\$url1\t\$url2\t\$source" >> NCBI_download_urls.tsv
                echo "Added \$run_id to PE file (source: \$source)" >&2
            elif [[ "\$layout" == "SINGLE" ]]; then
                echo "\$run_id" >> NCBI_SRR_SE_accessions.txt
                echo -e "\$run_id\tSE\t\$url1\t\t\$source" >> NCBI_download_urls.tsv
                echo "Added \$run_id to SE file (source: \$source)" >&2
            else
                echo "Warning: Unknown layout '\$layout' for \$run_id" >&2
            fi
        done
    done
    
    # Report results
    pe_count=\$(wc -l < NCBI_SRR_PE_accessions.txt)
    se_count=\$(wc -l < NCBI_SRR_SE_accessions.txt)
    ena_count=\$(grep -c "ENA" NCBI_download_urls.tsv || echo "0")

    echo "Done! Found \$pe_count paired-end and \$se_count single-end SRR accessions"
    echo "ENA URLs available for \$ena_count accessions"
    echo "PE accessions saved in NCBI_SRR_PE_accessions.txt"
    echo "SE accessions saved in NCBI_SRR_SE_accessions.txt"
    echo "Download URLs saved in NCBI_download_urls.tsv"
    """
}