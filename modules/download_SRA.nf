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
    
    # Configure SRA-tools to use current directory for temp files
    export NCBI_SETTINGS="\$PWD/.ncbi"
    mkdir -p "\$NCBI_SETTINGS"
    
    # Create config to use current directory for cache and temp
    cat > "\$NCBI_SETTINGS/user-settings.mkfg" << EOF
/repository/user/main/public/root = "\$PWD"
/repository/user/main/public/cache-enabled = "false"
EOF
    
    echo "Downloading $srr from $source"
    
    # Retry loop: attempt download up to 3 times if CRC validation fails
    max_attempts=3
    attempt=1
    success=false
    
    while [ \$attempt -le \$max_attempts ] && [ "\$success" = "false" ]; do
        if [ \$attempt -gt 1 ]; then
            echo "Retry attempt \$attempt of \$max_attempts..."
            # Clean up any partial files
            rm -f ${srr}_*.fastq* 2>/dev/null || true
            rm -rf ncbi_download 2>/dev/null || true
        fi
        
        # Try ENA direct download if URLs are available
        if [[ -n "$url1" && -n "$url2" && "$source" == "ENA" ]]; then
            echo "Attempting ENA direct download..."
            
            if wget -q -O ${srr}_1.fastq.gz "$url1" && wget -q -O ${srr}_2.fastq.gz "$url2"; then
                # Verify files are not empty
                if [ -s "${srr}_1.fastq.gz" ] && [ -s "${srr}_2.fastq.gz" ]; then
                    # Validate gzip integrity
                    if gzip -t ${srr}_1.fastq.gz 2>/dev/null && gzip -t ${srr}_2.fastq.gz 2>/dev/null; then
                        echo "✓ Successfully downloaded from ENA (CRC validated)"
                        success=true
                        break
                    else
                        echo "⚠ ENA files failed CRC validation (attempt \$attempt)"
                        rm -f ${srr}_*.fastq.gz
                    fi
                else
                    echo "⚠ ENA files empty, falling back to NCBI"
                    rm -f ${srr}_*.fastq.gz
                fi
            else
                echo "⚠ ENA download failed, falling back to NCBI"
                rm -f ${srr}_*.fastq.gz
            fi
        fi
        
        # Fallback to NCBI prefetch/fasterq-dump if ENA failed or was not attempted
        if [ "\$success" = "false" ]; then
            echo "Downloading from NCBI..."
            
            # Create explicit download directory to avoid /tmp/ quota issues
            mkdir -p ncbi_download
            mkdir -p fasterq_tmp
            
            # Use explicit output file path (not --output-directory which creates subdirs)
            # This ensures download happens in work directory, not /tmp/
            # Set max-size to 100GB to handle large sequencing runs
            if prefetch $srr --output-file ncbi_download/${srr}.sra --max-size 100G; then
                # Override TMPDIR to force fasterq-dump to use work directory
                # Even with --temp, fasterq-dump may use system TMPDIR for buffers
                export TMPDIR="\$PWD/fasterq_tmp"
                
                # Decompress with explicit temp directory in work location
                fasterq-dump ncbi_download/${srr}.sra --split-files -O . --temp fasterq_tmp
                
                # Clean up SRA file and temp directory immediately to save space
                rm -rf ncbi_download fasterq_tmp
                
                # Verify we got paired-end files
                if [ -f "${srr}_1.fastq" ] && [ -f "${srr}_2.fastq" ]; then
                    # Check files are not empty
                    if [ -s "${srr}_1.fastq" ] && [ -s "${srr}_2.fastq" ]; then
                        # Validate fastq format (check first read has 4 lines starting with @)
                        if head -n 1 ${srr}_1.fastq 2>/dev/null | grep -q '^@' && \
                           head -n 1 ${srr}_2.fastq 2>/dev/null | grep -q '^@'; then
                            echo "✓ Successfully downloaded from NCBI (format validated)"
                            success=true
                            break
                        else
                            echo "⚠ NCBI files failed format validation (attempt \$attempt)"
                        fi
                    else
                        echo "⚠ NCBI files empty (attempt \$attempt)"
                    fi
                else
                    echo "⚠ NCBI PE files not found (attempt \$attempt)"
                fi
            else
                echo "⚠ NCBI prefetch failed (attempt \$attempt)"
                rm -rf ncbi_download 2>/dev/null || true
            fi
        fi
        
        attempt=\$((attempt + 1))
    done
    
    # Final validation
    if [ "\$success" = "false" ]; then
        echo "ERROR: All download attempts failed for $srr after \$max_attempts tries"
        exit 1
    fi
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
    
    # Configure SRA-tools to use current directory for temp files
    export NCBI_SETTINGS="\$PWD/.ncbi"
    mkdir -p "\$NCBI_SETTINGS"
    
    # Create config to use current directory for cache and temp
    cat > "\$NCBI_SETTINGS/user-settings.mkfg" << EOF
/repository/user/main/public/root = "\$PWD"
/repository/user/main/public/cache-enabled = "false"
EOF
    
    echo "Downloading $srr from $source"
    
    # Retry loop: attempt download up to 3 times if CRC validation fails
    max_attempts=3
    attempt=1
    success=false
    
    while [ \$attempt -le \$max_attempts ] && [ "\$success" = "false" ]; do
        if [ \$attempt -gt 1 ]; then
            echo "Retry attempt \$attempt of \$max_attempts..."
            # Clean up any partial files
            rm -f ${srr}*.fastq* 2>/dev/null || true
            rm -rf ncbi_download 2>/dev/null || true
        fi
        
        # Try ENA direct download if URL is available
        if [[ -n "$url1" && "$source" == "ENA" ]]; then
            echo "Attempting ENA direct download..."
            
            if wget -q -O ${srr}.fastq.gz "$url1"; then
                # Verify file is not empty
                if [ -s "${srr}.fastq.gz" ]; then
                    # Validate gzip integrity
                    if gzip -t ${srr}.fastq.gz 2>/dev/null; then
                        echo "✓ Successfully downloaded from ENA (CRC validated)"
                        success=true
                        break
                    else
                        echo "⚠ ENA file failed CRC validation (attempt \$attempt)"
                        rm -f ${srr}.fastq.gz
                    fi
                else
                    echo "⚠ ENA file empty, falling back to NCBI"
                    rm -f ${srr}.fastq.gz
                fi
            else
                echo "⚠ ENA download failed, falling back to NCBI"
                rm -f ${srr}.fastq.gz
            fi
        fi
        
        # Fallback to NCBI prefetch/fasterq-dump if ENA failed or was not attempted
        if [ "\$success" = "false" ]; then
            echo "Downloading from NCBI..."
            
            # Create explicit download directory to avoid /tmp/ quota issues
            mkdir -p ncbi_download
            mkdir -p fasterq_tmp
            
            # Use explicit output file path to ensure download in work directory
            # Set max-size to 100GB to handle large sequencing runs
            if prefetch $srr --output-file ncbi_download/${srr}.sra --max-size 100G; then
                # Override TMPDIR to force fasterq-dump to use work directory
                export TMPDIR="\$PWD/fasterq_tmp"
                
                # Decompress with explicit temp directory in work location
                fasterq-dump ncbi_download/${srr}.sra --split-files -O . --temp fasterq_tmp
                
                # Clean up SRA file and temp directory immediately
                rm -rf ncbi_download fasterq_tmp
                
                # Handle SE naming variations (might be .fastq or _1.fastq)
                if [ -f "${srr}.fastq" ]; then
                    echo "Found ${srr}.fastq"
                    if [ -s "${srr}.fastq" ]; then
                        # Validate fastq format
                        if head -n 1 ${srr}.fastq 2>/dev/null | grep -q '^@'; then
                            echo "✓ Successfully downloaded from NCBI (format validated)"
                            success=true
                            break
                        else
                            echo "⚠ NCBI file failed format validation (attempt \$attempt)"
                        fi
                    else
                        echo "⚠ NCBI file empty (attempt \$attempt)"
                    fi
                elif [ -f "${srr}_1.fastq" ]; then
                    echo "Found ${srr}_1.fastq, renaming to ${srr}.fastq"
                    mv ${srr}_1.fastq ${srr}.fastq
                    rm -f ${srr}_2.fastq 2>/dev/null || true
                    if [ -s "${srr}.fastq" ]; then
                        # Validate fastq format
                        if head -n 1 ${srr}.fastq 2>/dev/null | grep -q '^@'; then
                            echo "✓ Successfully downloaded from NCBI (format validated)"
                            success=true
                            break
                        else
                            echo "⚠ NCBI file failed format validation (attempt \$attempt)"
                        fi
                    else
                        echo "⚠ NCBI file empty (attempt \$attempt)"
                    fi
                else
                    echo "⚠ NCBI file not found (attempt \$attempt)"
                fi
            else
                echo "⚠ NCBI prefetch failed (attempt \$attempt)"
                rm -rf ncbi_download 2>/dev/null || true
            fi
        fi
        
        attempt=\$((attempt + 1))
    done
    
    # Final validation
    if [ "\$success" = "false" ]; then
        echo "ERROR: All download attempts failed for $srr after \$max_attempts tries"
        exit 1
    fi
    """
}
