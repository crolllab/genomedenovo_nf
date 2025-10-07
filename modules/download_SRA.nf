// Hybrid download for paired-end reads: tries ENA URL first, falls back to NCBI
process SRAdownloadPE {
    cpus 1
    memory '4GB'
    tag "Downloading PE: $srr"
    errorStrategy 'ignore'
    
    input:
    tuple val(srr), val(url1), val(url2), val(source)
    
    output:
    tuple val(srr), path("${srr}_1.fastq*"), path("${srr}_2.fastq*")
    
    script:
    """
    #!/bin/bash
    set -e
    
    echo "Downloading $srr from $source"
    
    # Try ENA direct download if URLs are available
    if [[ -n "$url1" && -n "$url2" && "$source" == "ENA" ]]; then
        echo "Attempting ENA direct download..."
        
        if wget -q -O ${srr}_1.fastq.gz "$url1" && wget -q -O ${srr}_2.fastq.gz "$url2"; then
            # Verify files are not empty
            if [ -s "${srr}_1.fastq.gz" ] && [ -s "${srr}_2.fastq.gz" ]; then
                echo "✓ Successfully downloaded from ENA"
                exit 0
            else
                echo "⚠ ENA files empty, falling back to NCBI"
                rm -f ${srr}_*.fastq.gz
            fi
        else
            echo "⚠ ENA download failed, falling back to NCBI"
            rm -f ${srr}_*.fastq.gz
        fi
    fi
    
    # Fallback to NCBI prefetch/fasterq-dump
    echo "Downloading from NCBI..."
    
    prefetch $srr
    fasterq-dump $srr --split-files -O .
    
    # Verify we got paired-end files
    if [ ! -f "${srr}_1.fastq" ] || [ ! -f "${srr}_2.fastq" ]; then
        echo "ERROR: Expected PE files not found for $srr"
        exit 1
    fi
    
    # Check files are not empty
    if [ ! -s "${srr}_1.fastq" ] || [ ! -s "${srr}_2.fastq" ]; then
        echo "ERROR: Empty fastq files for $srr"
        exit 1
    fi
    
    # Clean up SRA file
    rm -rf ${srr}/ ${srr}.sra
    
    echo "✓ Successfully downloaded from NCBI"
    """
}

// Hybrid download for single-end reads: tries ENA URL first, falls back to NCBI
process SRAdownloadSE {
    cpus 1
    memory '4GB'
    tag "Downloading SE: $srr"
    errorStrategy 'ignore'
    
    input:
    tuple val(srr), val(url1), val(source)
    
    output:
    tuple val(srr), path("${srr}.fastq*")
    
    script:
    """
    #!/bin/bash
    set -e
    
    echo "Downloading $srr from $source"
    
    # Try ENA direct download if URL is available
    if [[ -n "$url1" && "$source" == "ENA" ]]; then
        echo "Attempting ENA direct download..."
        
        if wget -q -O ${srr}.fastq.gz "$url1"; then
            # Verify file is not empty
            if [ -s "${srr}.fastq.gz" ]; then
                echo "✓ Successfully downloaded from ENA"
                exit 0
            else
                echo "⚠ ENA file empty, falling back to NCBI"
                rm -f ${srr}.fastq.gz
            fi
        else
            echo "⚠ ENA download failed, falling back to NCBI"
            rm -f ${srr}.fastq.gz
        fi
    fi
    
    # Fallback to NCBI prefetch/fasterq-dump
    echo "Downloading from NCBI..."
    
    prefetch $srr
    fasterq-dump $srr --split-files -O .
    
    # Handle SE naming variations (might be .fastq or _1.fastq)
    if [ -f "${srr}.fastq" ]; then
        echo "Found ${srr}.fastq"
    elif [ -f "${srr}_1.fastq" ]; then
        echo "Found ${srr}_1.fastq, renaming to ${srr}.fastq"
        mv ${srr}_1.fastq ${srr}.fastq
        rm -f ${srr}_2.fastq 2>/dev/null || true
    else
        echo "ERROR: No fastq file found for $srr"
        exit 1
    fi
    
    # Check file is not empty
    if [ ! -s "${srr}.fastq" ]; then
        echo "ERROR: Empty fastq file for $srr"
        exit 1
    fi
    
    # Clean up SRA file
    rm -rf ${srr}/ ${srr}.sra
    
    echo "✓ Successfully downloaded from NCBI"
    """
}
