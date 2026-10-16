// Hybrid download for paired-end reads: tries ENA URL first, falls back to NCBI
process SRAdownloadPE {
    tag "Downloading PE: $srr"
    // A 404 from prefetch means the accession does not exist; retrying it just
    // burns hours of backoff, so exit 66 is treated as permanent.
    errorStrategy { task.exitStatus == 66 ? 'ignore' : (task.attempt <= 4 ? 'retry' : 'ignore') }
    maxRetries 4

    input:
    tuple val(srr), val(url1), val(url2), val(source)

    output:
    tuple val(srr), path("${srr}_1.fastq*"), path("${srr}_2.fastq*")

    script:
    """
    #!/bin/bash

    # Configure SRA-tools to use current directory for temp files
    export NCBI_SETTINGS="\$PWD/.ncbi"
    mkdir -p "\$NCBI_SETTINGS"

    # Create config to use current directory for cache and temp
    cat > "\$NCBI_SETTINGS/user-settings.mkfg" << EOF
/repository/user/main/public/root = "\$PWD"
/repository/user/main/public/cache-enabled = "false"
EOF

    echo "Downloading $srr from $source (attempt ${task.attempt} of 4)"

    # Retry loop: attempt download up to 4 times
    # Includes in-process exponential backoff to avoid hammering remote servers
    max_attempts=4
    attempt=1
    success=false

    # Timeouts (seconds) for the two long sra-tools steps. A ~8 GB .sra needs
    # well over ten minutes both to fetch and to convert, more again when the
    # work directory sits on network storage, so keep these generous.
    prefetch_timeout=3600
    fasterq_timeout=3600

    while [ \$attempt -le \$max_attempts ] && [ "\$success" = "false" ]; do
        if [ \$attempt -gt 1 ]; then
            # Long backoff with jitter to tolerate server-side outages/throttling:
            # attempt 2: 10 min, attempt 3: 30 min, attempt 4: 60 min (+0-59s jitter)
            case "\$attempt" in
                2) base_delay=600 ;;
                3) base_delay=1800 ;;
                4) base_delay=3600 ;;
                *) base_delay=600 ;;
            esac
            jitter=\$((RANDOM % 60))
            sleep_seconds=\$((base_delay + jitter))
            echo "Retry attempt \$attempt of \$max_attempts..."
            echo "Waiting \$sleep_seconds seconds before retry (long server cooldown)..."
            sleep \$sleep_seconds
            # Clean up partial FASTQ output, but keep ncbi_download: a complete
            # .sra is reused below when only the extraction step failed.
            rm -f ${srr}_*.fastq* 2>/dev/null || true
            rm -rf fasterq_tmp 2>/dev/null || true
        fi

        # Try ENA direct download if URLs are available
        if [[ -n "$url1" && -n "$url2" && "$source" == "ENA" ]]; then
            echo "Attempting ENA direct download..."

            # Ensure HTTPS is used (http:// hangs on compute nodes; https:// is reliable)
            dl_url1=\$(echo "$url1" | sed 's|^http://|https://|')
            dl_url2=\$(echo "$url2" | sed 's|^http://|https://|')

            # BusyBox wget: -T for network read timeout (seconds), -t for tries.
            # -T is accepted by both BusyBox and GNU wget; --timeout= is not portable.
            if wget -T 300 -t 3 -q -O ${srr}_1.fastq.gz "\$dl_url1" && \
               wget -T 300 -t 3 -q -O ${srr}_2.fastq.gz "\$dl_url2"; then
                # Verify files are not empty
                if [ -s "${srr}_1.fastq.gz" ] && [ -s "${srr}_2.fastq.gz" ]; then
                    # Validate gzip integrity
                    if gzip -t ${srr}_1.fastq.gz 2>/dev/null && gzip -t ${srr}_2.fastq.gz 2>/dev/null; then
                        # Validate paired-end read count consistency
                        echo "Validating read counts..."
                        count1=\$(zcat ${srr}_1.fastq.gz | wc -l | awk '{print int(\$1/4)}')
                        count2=\$(zcat ${srr}_2.fastq.gz | wc -l | awk '{print int(\$1/4)}')

                        if [ "\$count1" -eq "\$count2" ]; then
                            echo "✓ Successfully downloaded from ENA (CRC validated, \$count1 paired reads)"
                            success=true
                            break
                        else
                            echo "⚠ ENA read count mismatch: R1=\$count1, R2=\$count2 (attempt \$attempt)"
                            rm -f ${srr}_*.fastq.gz
                        fi
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

            # Reuse a .sra left behind by an earlier attempt. Extraction can fail
            # on its own (timeout, disk), and re-fetching several GB every round
            # costs far more than it saves. prefetch only writes the final
            # .sra name once the transfer completes, so a partial download is
            # never picked up here.
            sra_file=\$(find ncbi_download -type f -name "${srr}.sra" 2>/dev/null | head -n 1)
            if [ -n "\$sra_file" ]; then
                echo "Reusing .sra fetched on an earlier attempt: \$sra_file"
            else
                # Download into a controlled directory; newer sra-tools deprecates --output-file.
                # Set max-size to 100GB to handle large sequencing runs.
                # Use timeout to prevent hangs (prefetch occasionally hangs on network timeouts).
                # Some sra-tools builds can return non-zero despite a complete .sra being written,
                # so gate on file presence rather than prefetch exit status alone.
                timeout \$prefetch_timeout prefetch $srr --output-directory ncbi_download --max-size 100G > prefetch_${srr}.log 2>&1 || true
                cat prefetch_${srr}.log
                sra_file=\$(find ncbi_download -type f -name "${srr}.sra" 2>/dev/null | head -n 1)

                # A 404 / unresolvable accession is permanent. Fail fast with a
                # dedicated exit code so errorStrategy skips the retry ladder.
                if [ -z "\$sra_file" ] && grep -qE "failed to resolve accession|no data.*404" prefetch_${srr}.log 2>/dev/null; then
                    echo "ERROR: accession $srr does not exist in SRA/ENA (404) -- not retrying"
                    exit 66
                fi
            fi
            if [ -n "\$sra_file" ]; then
                # Override TMPDIR to force fasterq-dump to use work directory
                # Even with --temp, fasterq-dump may use system TMPDIR for buffers
                export TMPDIR="\$PWD/fasterq_tmp"

                # Decompress with explicit temp directory in work location.
                # Capture the status rather than testing it inline: a timeout kill
                # exits 124 and writes nothing to stderr, which is indistinguishable
                # from a silent failure unless the code itself is reported.
                timeout \$fasterq_timeout fasterq-dump "\$sra_file" -e ${task.cpus} --split-files -O . --temp fasterq_tmp 2>fasterq_${srr}.err
                fq_status=\$?
                if [ \$fq_status -eq 0 ]; then
                    # Clean up SRA file and temp directory immediately to save space
                    rm -rf ncbi_download fasterq_tmp

                    # Verify we got paired-end files
                    if [ -f "${srr}_1.fastq" ] && [ -f "${srr}_2.fastq" ]; then
                        # Check files are not empty
                        if [ -s "${srr}_1.fastq" ] && [ -s "${srr}_2.fastq" ]; then
                            # Validate fastq format (check first read has 4 lines starting with @)
                            if head -n 1 ${srr}_1.fastq 2>/dev/null | grep -q '^@' && \
                               head -n 1 ${srr}_2.fastq 2>/dev/null | grep -q '^@'; then

                                # Validate paired-end read count consistency
                                echo "Validating read counts..."
                                count1=\$(cat ${srr}_1.fastq | wc -l | awk '{print int(\$1/4)}')
                                count2=\$(cat ${srr}_2.fastq | wc -l | awk '{print int(\$1/4)}')

                                if [ "\$count1" -eq "\$count2" ]; then
                                    echo "✓ Successfully downloaded from NCBI (format validated, \$count1 paired reads)"
                                    success=true
                                    break
                                else
                                    echo "⚠ NCBI read count mismatch: R1=\$count1, R2=\$count2 (attempt \$attempt)"
                                    echo "   This indicates data loss during fasterq-dump extraction"
                                    rm -f ${srr}_*.fastq
                                fi
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
                    if [ \$fq_status -eq 124 ]; then
                        echo "⚠ NCBI fasterq-dump timed out after \${fasterq_timeout}s (attempt \$attempt)"
                        echo "   The .sra is kept for the next attempt; raise fasterq_timeout if this repeats"
                    else
                        echo "⚠ NCBI fasterq-dump failed with exit status \$fq_status (attempt \$attempt)"
                    fi
                    if [ -s fasterq_${srr}.err ]; then
                        echo "   Error: \$(tail -2 fasterq_${srr}.err 2>/dev/null)"
                    fi
                fi
            else
                echo "⚠ NCBI prefetch failed to produce .sra file (attempt \$attempt)"
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
    tag "Downloading SE: $srr"
    errorStrategy { task.exitStatus == 66 ? 'ignore' : (task.attempt <= 4 ? 'retry' : 'ignore') }
    maxRetries 4

    input:
    tuple val(srr), val(url1), val(source)

    output:
    tuple val(srr), path("${srr}.fastq*")

    script:
    """
    #!/bin/bash

    # Configure SRA-tools to use current directory for temp files
    export NCBI_SETTINGS="\$PWD/.ncbi"
    mkdir -p "\$NCBI_SETTINGS"

    # Create config to use current directory for cache and temp
    cat > "\$NCBI_SETTINGS/user-settings.mkfg" << EOF
/repository/user/main/public/root = "\$PWD"
/repository/user/main/public/cache-enabled = "false"
EOF

    echo "Downloading $srr from $source (attempt ${task.attempt} of 4)"

    max_attempts=4
    attempt=1
    success=false

    # Timeouts (seconds) for the two long sra-tools steps. A ~8 GB .sra needs
    # well over ten minutes both to fetch and to convert, more again when the
    # work directory sits on network storage, so keep these generous.
    prefetch_timeout=3600
    fasterq_timeout=3600

    while [ \$attempt -le \$max_attempts ] && [ "\$success" = "false" ]; do
        if [ \$attempt -gt 1 ]; then
            case "\$attempt" in
                2) base_delay=600 ;;
                3) base_delay=1800 ;;
                4) base_delay=3600 ;;
                *) base_delay=600 ;;
            esac
            jitter=\$((RANDOM % 60))
            sleep_seconds=\$((base_delay + jitter))
            echo "Retry attempt \$attempt of \$max_attempts..."
            echo "Waiting \$sleep_seconds seconds before retry (long server cooldown)..."
            sleep \$sleep_seconds
            # Clean up partial FASTQ output, but keep ncbi_download: a complete
            # .sra is reused below when only the extraction step failed.
            rm -f ${srr}*.fastq* 2>/dev/null || true
            rm -rf fasterq_tmp 2>/dev/null || true
        fi

        # Try ENA direct download if URL is available
        if [[ -n "$url1" && "$source" == "ENA" ]]; then
            echo "Attempting ENA direct download..."

            # Ensure HTTPS is used (http:// hangs on compute nodes; https:// is reliable)
            dl_url1=\$(echo "$url1" | sed 's|^http://|https://|')

            # BusyBox wget: -T for network read timeout (seconds), -t for tries.
            # -T is accepted by both BusyBox and GNU wget; --timeout= is not portable.
            if wget -T 300 -t 3 -q -O ${srr}.fastq.gz "\$dl_url1"; then
                if [ -s "${srr}.fastq.gz" ]; then
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

            mkdir -p ncbi_download
            mkdir -p fasterq_tmp

            # Reuse a .sra left behind by an earlier attempt (see PE process).
            sra_file=\$(find ncbi_download -type f -name "${srr}.sra" 2>/dev/null | head -n 1)
            if [ -n "\$sra_file" ]; then
                echo "Reusing .sra fetched on an earlier attempt: \$sra_file"
            else
                timeout \$prefetch_timeout prefetch $srr --output-directory ncbi_download --max-size 100G > prefetch_${srr}.log 2>&1 || true
                cat prefetch_${srr}.log
                sra_file=\$(find ncbi_download -type f -name "${srr}.sra" 2>/dev/null | head -n 1)

                # A 404 / unresolvable accession is permanent. Fail fast with a
                # dedicated exit code so errorStrategy skips the retry ladder.
                if [ -z "\$sra_file" ] && grep -qE "failed to resolve accession|no data.*404" prefetch_${srr}.log 2>/dev/null; then
                    echo "ERROR: accession $srr does not exist in SRA/ENA (404) -- not retrying"
                    exit 66
                fi
            fi
            if [ -n "\$sra_file" ]; then
                export TMPDIR="\$PWD/fasterq_tmp"

                timeout \$fasterq_timeout fasterq-dump "\$sra_file" -e ${task.cpus} --split-files -O . --temp fasterq_tmp 2>fasterq_${srr}.err
                fq_status=\$?
                if [ \$fq_status -eq 0 ]; then
                    rm -rf ncbi_download fasterq_tmp

                    # Handle SE naming variations (might be .fastq or _1.fastq).
                    # These are separate `if`s, not if/elif: fasterq-dump can emit
                    # ${srr}_1.fastq even when ${srr}.fastq is absent or empty.
                    if [ -f "${srr}.fastq" ]; then
                        echo "Found ${srr}.fastq"
                        if [ -s "${srr}.fastq" ]; then
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
                    fi
                    if [ -f "${srr}_1.fastq" ]; then
                        echo "Found ${srr}_1.fastq, renaming to ${srr}.fastq"
                        mv ${srr}_1.fastq ${srr}.fastq
                        rm -f ${srr}_2.fastq 2>/dev/null || true
                        if [ -s "${srr}.fastq" ]; then
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
                    if [ \$fq_status -eq 124 ]; then
                        echo "⚠ NCBI fasterq-dump timed out after \${fasterq_timeout}s (attempt \$attempt)"
                        echo "   The .sra is kept for the next attempt; raise fasterq_timeout if this repeats"
                    else
                        echo "⚠ NCBI fasterq-dump failed with exit status \$fq_status (attempt \$attempt)"
                    fi
                    if [ -s fasterq_${srr}.err ]; then
                        echo "   Error: \$(tail -2 fasterq_${srr}.err 2>/dev/null)"
                    fi
                fi
            else
                echo "⚠ NCBI prefetch failed to produce .sra file (attempt \$attempt)"
                rm -rf ncbi_download 2>/dev/null || true
            fi
        fi

        attempt=\$((attempt + 1))
    done

    if [ "\$success" = "false" ]; then
        echo "ERROR: All download attempts failed for $srr after \$max_attempts tries"
        exit 1
    fi
    """
}
