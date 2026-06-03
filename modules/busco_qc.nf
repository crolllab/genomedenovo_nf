// BUSCO process for genome completeness assessment
process BUSCO {
    tag "BUSCO assessment"
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    publishDir "${params.outdir}/4_busco", mode: 'copy', pattern: "*_busco"

    input:
    tuple val(sample_id), path(scaffolds)

    output:
    tuple val(sample_id), path("${sample_id}_busco"), emit: result_dir
    path "short_summary.*.txt",                        emit: summary
    tuple val(sample_id), path("${sample_id}_scogs"), emit: single_copy_dir

    script:
    def db_arg = params.busco_db_path
        ? "--offline --download_path ${params.busco_db_path}"
        : "--download_path ./busco_downloads"
    """
    set -euo pipefail

    busco \\
        -m genome \\
        -i ${scaffolds} \\
        -l ${params.busco_lineage} \\
        -o ${sample_id}_busco \\
        -c ${task.cpus} \\
        ${db_arg}

    cp ${sample_id}_busco/short_summary.*.txt .

    # Extract single-copy BUSCO nucleotide sequences for gene phylogenomics (--phylo_tree)
    mkdir -p ${sample_id}_scogs
    for fna in ${sample_id}_busco/run_*/busco_sequences/single_copy_busco_sequences/*.fna; do
        [ -f "\$fna" ] || continue
        gene=\$(basename "\$fna" .fna)
        # Replace FASTA header with sample_id so sequences are identifiable after grouping
        awk -v sid="${sample_id}" '/^>/{print ">" sid; next} {print}' "\$fna" \\
            > ${sample_id}_scogs/\${gene}.fna
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
