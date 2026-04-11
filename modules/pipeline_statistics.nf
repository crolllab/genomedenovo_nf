process PipelineStatistics {
    tag "Generating pipeline stats"
    errorStrategy 'retry'
    maxRetries 3
    
    publishDir params.outdir, mode: 'copy'
    
    input:
    val ready
    path trace_file
    
    output:
    path "pipeline_execution_stats.txt"
    path "pipeline_execution_stats.tsv"
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    OUTFILE="pipeline_execution_stats.txt"
    TSV_FILE="pipeline_execution_stats.tsv"
    TRACE_FILE="${trace_file}"
    
    echo "=================================================================" > \$OUTFILE
    echo "  Nextflow Pipeline - Process Statistics" >> \$OUTFILE
    echo "  Generated: \$(date)" >> \$OUTFILE
    echo "=================================================================" >> \$OUTFILE
    echo "" >> \$OUTFILE
    
    printf "%-30s %6s %14s %14s %14s %10s %10s\\n" \\
        "Process" "Count" "Avg_time" "Min_time" "Max_time" "AvgPeakGB" "MaxPeakGB" >> \$OUTFILE
    echo "---------------------------------------------------------------------------------------------------" >> \$OUTFILE

    printf "Process\\tCount\\tAvg_time\\tMin_time\\tMax_time\\tAvgPeakGB\\tMaxPeakGB\\n" > \$TSV_FILE

    fmt_time() {
        local s=\$(printf '%.0f' "\$1")
        local h=\$((s/3600)); local m=\$((s%3600/60)); local sec=\$((s%60))
        if   [ \$h -gt 0 ]; then printf '%dh %dm %ds' \$h \$m \$sec
        elif [ \$m -gt 0 ]; then printf '%dm %ds' \$m \$sec
        else printf '%ds' \$sec; fi
    }

    fmt_gb() { awk -v mb="\$1" 'BEGIN{printf "%.3f", mb/1024}'; }
    
    if [ ! -s "\$TRACE_FILE" ]; then
        echo "WARNING: Trace file not found or empty: \$TRACE_FILE" >> \$OUTFILE
        echo "No process statistics were generated." >> \$OUTFILE
        echo "" >> \$OUTFILE
        echo "=================================================================" >> \$OUTFILE
        echo "Notes:" >> \$OUTFILE
        echo "  - Statistics are derived from pipeline_trace.txt" >> \$OUTFILE
        echo "  - Ensure trace is enabled and run has completed enough tasks" >> \$OUTFILE
        echo "=================================================================" >> \$OUTFILE
        echo "Statistics generated with warnings"
        exit 0
    fi

    parsed_file=\$(mktemp)

    # Parse completed successful tasks; convert duration/realtime to seconds and memory to MB.
    awk -F '\\t' '
    function trim(s) {
        gsub(/^[ \\t]+|[ \\t]+\$/, "", s)
        return s
    }

    function to_seconds(raw,    s, n, i, tok, val, unit, total) {
        s = trim(raw)
        if (s == "") return -1

        total = 0
        n = split(s, tok, / +/)
        for (i = 1; i <= n; i++) {
            if (tok[i] == "") continue

            unit = tok[i]
            gsub(/^[0-9]*\\.?[0-9]+/, "", unit)
            val = tok[i]
            gsub(/[a-zA-Z]+\$/, "", val)
            val = val + 0

            if (unit == "ms") total += (val / 1000)
            else if (unit == "s") total += val
            else if (unit == "m") total += (val * 60)
            else if (unit == "h") total += (val * 3600)
            else if (unit == "d") total += (val * 86400)
            else if (unit == "") total += val
            else return -1
        }

        return total
    }

    function to_megabytes(raw,    s, val, unit) {
        s = trim(raw)
        if (s == "") return -1

        val  = s; gsub(/[a-zA-Z ].*/, "", val); val = val + 0
        unit = s; gsub(/^[0-9. ]*/, "", unit);  gsub(/ .*/, "", unit)

        if      (unit == "B")  return val / (1024*1024)
        else if (unit == "KB") return val / 1024
        else if (unit == "MB") return val
        else if (unit == "GB") return val * 1024
        else if (unit == "TB") return val * 1024 * 1024
        else if (unit == "")   return 0
        else                   return -1
    }

    NR == 1 { next }

    {
        process_name = \$4
        status       = \$7
        exit_code    = \$8
        realtime     = \$13
        duration     = \$12
        peak_mem     = \$18

        sub(/^.*:/, "", process_name)

        sec     = to_seconds(realtime)
        if (sec < 0) sec = to_seconds(duration)
        peak_mb = to_megabytes(peak_mem)

        if (status == "COMPLETED" && exit_code == "0" && sec >= 0) {
            printf "%s\\t%.6f\\t%.2f\\n", process_name, sec, peak_mb
        }
    }
    ' "\$TRACE_FILE" > "\$parsed_file"

    known_processes=(
        SRAresolve SRAdownloadPE SRAdownloadSE
        SPAdes SPAdesSE
        QUAST SummarizeAssemblyMetrics
        SplitQuerySequences BLASTsearch
        ExtractBlastHits CombineBlastHits
        MergeBlastSummaries
        AlignMAFFT
    )

    # Append any processes found in trace that are not in the known list
    mapfile -t trace_processes < <(awk -F '\\t' '{print \$1}' "\$parsed_file" | sort | uniq)
    for tp in "\${trace_processes[@]}"; do
        found=0
        for kp in "\${known_processes[@]}"; do
            [[ "\$tp" == "\$kp" ]] && found=1 && break
        done
        [[ \$found -eq 0 ]] && known_processes+=("\$tp")
    done

    for proc in "\${known_processes[@]}"; do
        stats=\$(awk -F '\\t' -v p="\$proc" '
        \$1 == p {
            count++
            sum += \$2
            if (count == 1 || \$2 < min) min = \$2
            if (count == 1 || \$2 > max) max = \$2
            memsum += \$3
            if (count == 1 || \$3 > peakmax) peakmax = \$3
        }
        END {
            if (count > 0) {
                printf "%d\\t%.0f\\t%.0f\\t%.0f\\t%.2f\\t%.2f", \\
                    count, (sum / count), min, max, (memsum / count), peakmax
            }
        }
        ' "\$parsed_file")

        if [ -n "\$stats" ]; then
            IFS=\$'\\t' read -r count avg min max avgmem maxmem <<< "\$stats"
            t_avg=\$(fmt_time "\$avg"); t_min=\$(fmt_time "\$min"); t_max=\$(fmt_time "\$max")
            g_avg=\$(fmt_gb "\$avgmem");  g_max=\$(fmt_gb "\$maxmem")
            printf "%-30s %6d %14s %14s %14s %10s %10s\\n" \\
                "\$proc" "\$count" "\$t_avg" "\$t_min" "\$t_max" "\$g_avg" "\$g_max" >> "\$OUTFILE"
            printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
                "\$proc" "\$count" "\$t_avg" "\$t_min" "\$t_max" "\$g_avg" "\$g_max" >> "\$TSV_FILE"
        else
            printf "%-30s %6d %14s %14s %14s %10s %10s\\n" \\
                "\$proc" "0" "." "." "." "." "." >> "\$OUTFILE"
            printf "%s\\t0\\t.\\t.\\t.\\t.\\t.\\n" "\$proc" >> "\$TSV_FILE"
        fi
    done

    echo "" >> \$OUTFILE
    echo "=================================================================" >> \$OUTFILE
    echo "Notes:" >> \$OUTFILE
    echo "  - Statistics from completed tasks only" >> \$OUTFILE
    echo "  - Times and memory from Nextflow trace realtime/duration/memory fields" >> \$OUTFILE
    echo "  - Trace file: \$TRACE_FILE" >> \$OUTFILE
    echo "=================================================================" >> \$OUTFILE
    
    rm -f "\$parsed_file"

    echo "Statistics generated successfully"
    """
}
