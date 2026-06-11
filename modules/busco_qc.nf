// Download BUSCO lineage database once (cached via storeDir to avoid re-downloading on -resume)
process DownloadBUSCODatabase {
    tag "Downloading BUSCO database: ${params.busco_lineage}"
    storeDir "${params.outdir}/busco_db_cache"

    output:
    path "busco_downloads", emit: db_dir

    script:
    """
    set -euo pipefail
    busco --download ${params.busco_lineage} --download_path busco_downloads

    # Verify the lineage was actually downloaded
    lineage_dir="busco_downloads/lineages/${params.busco_lineage}"
    if [ ! -d "\$lineage_dir" ]; then
        echo "" >&2
        echo "ERROR: BUSCO download failed — lineage '${params.busco_lineage}' was not found." >&2
        echo "" >&2
        echo "Available lineages in the BUSCO database:" >&2
        ls busco_downloads/lineages/ 2>/dev/null | sort >&2 || echo "  (no lineages found)" >&2
        echo "" >&2
        echo "Common fungal lineages: fungi_odb10, ascomycota_odb10, basidiomycota_odb10," >&2
        echo "  dothideomycetes_odb10, sordariomycetes_odb10" >&2
        echo "Browse all lineages at: https://busco-data.ezlab.org/v5/data/lineages/" >&2
        exit 1
    fi
    """
}

// BUSCO process for genome completeness assessment
process BUSCO {
    tag "BUSCO assessment"
    // Exit code 2 from our pre-check signals a configuration error (e.g. missing lineage) — skip retries.
    // NOTE: bash syntax errors also produce exit code 2, so we use 'ignore' rather than 'finish'
    // to avoid aborting the entire pipeline when a single sample has a problematic name.
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    publishDir "${params.outdir}/4_busco", mode: 'copy', pattern: "*_busco"

    input:
    tuple val(sample_id), path(scaffolds)
    path db_dir

    output:
    tuple val(sample_id), path("${sample_id}_busco"), emit: result_dir
    path "short_summary.*.txt",                        emit: summary
    tuple val(sample_id), path("${sample_id}_scogs"), emit: single_copy_dir

    script:
    """
    set -euo pipefail

    # Pre-check: verify the lineage exists before attempting to run BUSCO.
    # Exit with code 2 (config error) to skip retries and surface a clear message.
    lineage_dir="${db_dir}/lineages/${params.busco_lineage}"
    if [ ! -d "\$lineage_dir" ]; then
        echo "" >&2
        echo "ERROR: BUSCO lineage '${params.busco_lineage}' not found in the database directory." >&2
        echo "Expected path: \$lineage_dir" >&2
        echo "" >&2
        echo "Available lineages:" >&2
        ls "${db_dir}/lineages/" 2>/dev/null | sort >&2 || echo "  (no lineages directory found in ${db_dir})" >&2
        echo "" >&2
        echo "Common fungal lineages: fungi_odb10, ascomycota_odb10, basidiomycota_odb10," >&2
        echo "  dothideomycetes_odb10, sordariomycetes_odb10" >&2
        echo "Browse all lineages at: https://busco-data.ezlab.org/v5/data/lineages/" >&2
        exit 2
    fi

    busco \\
        -m genome \\
        -i "${scaffolds}" \\
        -l ${params.busco_lineage} \\
        -o "${sample_id}_busco" \\
        -c ${task.cpus} \\
        --offline --download_path ${db_dir}

    cp "${sample_id}_busco"/short_summary.*.txt .

    # Extract single-copy BUSCO protein sequences for gene phylogenomics (--phylo_tree).
    # BUSCO v6 uses miniprot (default), which writes .faa protein sequences in
    # single_copy_busco_sequences/ — no nucleotide .fna files are produced.
    mkdir -p "${sample_id}_scogs"
    for faa in "${sample_id}_busco"/run_*/busco_sequences/single_copy_busco_sequences/*.faa; do
        [ -f "\$faa" ] || continue
        gene=\$(basename "\$faa" .faa)
        # Replace FASTA header with sample_id so sequences are identifiable after grouping
        awk -v sid="${sample_id}" '/^>/{print ">" sid; next} {print}' "\$faa" \\
            > "${sample_id}_scogs/\${gene}.faa"
    done
    """
}

// Collect and summarize all BUSCO short summaries into a single TSV
process SummarizeBuscoMetrics {
    tag "Summarizing BUSCO metrics"
    publishDir "${params.outdir}/4_busco", mode: 'copy'

    input:
    path summaries

    output:
    path "busco_summary_metrics.tsv"

    script:
    """
    #!/usr/bin/env python3
    import glob
    import re
    import os

    summary_files = sorted(glob.glob("short_summary.*.txt"))
    rows = []

    for f in summary_files:
        basename = os.path.basename(f)
        m = re.match(r'short_summary\\.\\w+\\.\\S+\\.(.+)_busco\\.txt', basename)
        sample_id = m.group(1) if m else basename

        data = {k: '' for k in ['Sample_ID', 'Lineage', 'Complete_pct', 'Single_copy_pct',
                                 'Duplicated_pct', 'Fragmented_pct', 'Missing_pct',
                                 'Complete', 'Single_copy', 'Duplicated', 'Fragmented',
                                 'Missing', 'Total']}
        data['Sample_ID'] = sample_id

        with open(f) as fh:
            for line in fh:
                line = line.strip()

                lm = re.search(r'Lineage dataset:\\s+(\\S+)', line)
                if lm:
                    data['Lineage'] = lm.group(1)

                sm = re.search(
                    r'C:(\\d+\\.\\d+)%\\[S:(\\d+\\.\\d+)%,D:(\\d+\\.\\d+)%\\],'
                    r'F:(\\d+\\.\\d+)%,M:(\\d+\\.\\d+)%,n:(\\d+)', line)
                if sm:
                    g = sm.groups()
                    data['Complete_pct']    = g[0]
                    data['Single_copy_pct'] = g[1]
                    data['Duplicated_pct']  = g[2]
                    data['Fragmented_pct']  = g[3]
                    data['Missing_pct']     = g[4]
                    data['Total']           = g[5]

                for pattern, key in [
                    (r'^(\\d+)\\s+Complete BUSCOs \\(C\\)',                 'Complete'),
                    (r'^(\\d+)\\s+Complete and single-copy BUSCOs \\(S\\)', 'Single_copy'),
                    (r'^(\\d+)\\s+Complete and duplicated BUSCOs \\(D\\)',  'Duplicated'),
                    (r'^(\\d+)\\s+Fragmented BUSCOs \\(F\\)',               'Fragmented'),
                    (r'^(\\d+)\\s+Missing BUSCOs \\(M\\)',                  'Missing'),
                ]:
                    pm = re.match(pattern, line)
                    if pm:
                        data[key] = pm.group(1)

        rows.append(data)

    headers = ['Sample_ID', 'Lineage', 'Complete_pct', 'Single_copy_pct', 'Duplicated_pct',
               'Fragmented_pct', 'Missing_pct', 'Complete', 'Single_copy', 'Duplicated',
               'Fragmented', 'Missing', 'Total']
    with open("busco_summary_metrics.tsv", "w") as out:
        out.write("\\t".join(headers) + "\\n")
        for row in rows:
            out.write("\\t".join(str(row.get(h, '')) for h in headers) + "\\n")
    """
}
