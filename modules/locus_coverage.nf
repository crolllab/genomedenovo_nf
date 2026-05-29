// ---------------------------------------------------------------------------
// Locus coverage module
//
// Workflow:
//   1. BlastLocusOnReference  – BLAST each query against the reference genome
//      to identify the best-hit locus (one run per query, not per sample).
//   2. ComputeLocusCoverage   – For each sample BAM × locus:
//        a. Extract per-base read depth at the padded locus (MQ >= 20).
//        b. Compute mean window depth in 1 000 bp intervals across the whole
//           genome (MQ >= 20) in a single streaming pass.
//        c. Report locus mean/median depth and genome-wide window mean/SD.
//   3. MergeLocusCoverage     – Concatenate per-sample TSVs per query.
// ---------------------------------------------------------------------------

process BlastLocusOnReference {
    tag "BLAST locus on reference: ${query_id}"
    errorStrategy 'ignore'
    publishDir { "${params.outdir}/locus_coverage/${query_id}" }, mode: 'copy'

    input:
    tuple val(query_id), path(query_fasta), val(seq_type), path(reference)

    output:
    tuple val(query_id), path("*_ref_locus.tsv"), emit: locus

    script:
    def qid_safe   = query_id.replaceAll(/[^A-Za-z0-9._-]/, '_')
    def blast_prog  = (seq_type == 'protein') ? 'tblastn' : 'blastn'
    def blast_extra = (seq_type == 'protein') ? '' : '-task blastn'
    """
    set -euo pipefail

    makeblastdb -in ${reference} -dbtype nucl -out ref_db -logfile /dev/null

    ${blast_prog} -query ${query_fasta} -db ref_db ${blast_extra} \\
        -outfmt '6 sseqid sstart send sstrand evalue bitscore' \\
        -max_target_seqs 1 -max_hsps 1 \\
        -out blast_locus.txt

    # Write header; append best hit if found
    printf 'query_id\\tchrom\\thit_start\\thit_end\\n' > ${qid_safe}_ref_locus.tsv

    if [ -s blast_locus.txt ]; then
        awk -v qid="${query_id}" 'NR==1 {
            s = \$2; e = \$3
            if (s > e) { tmp = s; s = e; e = tmp }
            print qid "\\t" \$1 "\\t" s "\\t" e
        }' blast_locus.txt >> ${qid_safe}_ref_locus.tsv
    else
        echo "WARNING: No BLAST hit found for ${query_id} on reference genome" >&2
    fi
    """
}


process ComputeLocusCoverage {
    tag "Locus coverage: ${sample_id} / ${query_id}"
    errorStrategy 'ignore'
    publishDir { "${params.outdir}/locus_coverage/${query_id}" }, mode: 'copy'

    input:
    tuple val(query_id), val(sample_id), path(bam), path(bai), path(locus_tsv)

    output:
    tuple val(query_id), val(sample_id),
          path("${sample_id}_*_locus_coverage.tsv"), emit: coverage

    script:
    def qid_safe = query_id.replaceAll(/[^A-Za-z0-9._-]/, '_')
    """
    set -euo pipefail

    OUT="${sample_id}_${qid_safe}_locus_coverage.tsv"

    CHROM=\$(awk -F'\\t' 'NR==2{print \$2}' ${locus_tsv})
    HIT_START=\$(awk -F'\\t' 'NR==2{print \$3}' ${locus_tsv})
    HIT_END=\$(awk   -F'\\t' 'NR==2{print \$4}' ${locus_tsv})

    # Write header unconditionally
    printf 'sample_id\\tquery_id\\tscaffold\\thit_start\\thit_end\\t' > "\$OUT"
    printf 'locus_start\\tlocus_end\\tlocus_mean_depth\\tlocus_median_depth\\t' >> "\$OUT"
    printf 'genome_mean_depth\tgenome_sd_depth\n'      >> "\$OUT"

    # No BLAST hit for this query on the reference — write NA row and exit
    if [ -z "\$CHROM" ]; then
        printf '%s\\t%s\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\n' \\
            "${sample_id}" "${query_id}" >> "\$OUT"
        exit 0
    fi

    # ── Padded locus coordinates (clamped to chromosome boundaries) ──────────
    CHROM_LEN=\$(samtools view -H ${bam} | awk -v c="\$CHROM" '
        /^@SQ/ {
            sn = ""; ln = 0
            for (i = 1; i <= NF; i++) {
                if (\$i ~ /^SN:/) sn = substr(\$i, 4)
                if (\$i ~ /^LN:/) ln = substr(\$i, 4)
            }
            if (sn == c) { print ln; exit }
        }')
    [ -z "\$CHROM_LEN" ] && CHROM_LEN=999999999

    LOCUS_START=\$(( HIT_START - ${params.padding} > 1          ? HIT_START - ${params.padding} : 1           ))
    LOCUS_END=\$(( HIT_END   + ${params.padding} < CHROM_LEN  ? HIT_END   + ${params.padding} : CHROM_LEN   ))

    # ── Per-base depth at padded locus (MQ >= 20) ────────────────────────────
    samtools depth -q 20 -a -r "\${CHROM}:\${LOCUS_START}-\${LOCUS_END}" ${bam} \\
        > locus_depth.txt

    LOCUS_MEAN=\$(awk '
        { s += \$3; n++ }
        END { if (n > 0) printf "%.2f", s/n; else print "NA" }
    ' locus_depth.txt)

    LOCUS_MEDIAN=\$(awk '{print \$3}' locus_depth.txt | sort -n | awk '
        { a[NR] = \$1 }
        END {
            n = NR
            if (n == 0) { print "NA"; exit }
            if (n % 2 == 1) print a[(n+1)/2]
            else printf "%.1f", (a[n/2] + a[n/2+1]) / 2
        }')

    # ── Genome-wide per-bp depth (MQ >= 20): mean and SD via single streaming pass ──
    samtools depth -q 20 -a ${bam} | awk '
    BEGIN { s=0; s2=0; n=0 }
    { d = \$3; s += d; s2 += d*d; n++ }
    END {
        if (n > 0) { mn = s/n } else { mn = 0 }
        if (n > 1) { sd = sqrt((s2 - n*mn*mn) / (n-1)) } else { sd = 0 }
        printf "%.2f\\t%.2f\\n", mn, sd
    }' > genome_window_stats.txt

    GENOME_MEAN=\$(awk '{print \$1}' genome_window_stats.txt)
    GENOME_SD=\$(  awk '{print \$2}' genome_window_stats.txt)

    # ── Write data row ────────────────────────────────────────────────────────
    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \\
        "${sample_id}" "${query_id}" \\
        "\$CHROM" "\$HIT_START" "\$HIT_END" \\
        "\$LOCUS_START" "\$LOCUS_END" \\
        "\$LOCUS_MEAN" "\$LOCUS_MEDIAN" \\
        "\$GENOME_MEAN" "\$GENOME_SD" \\
        >> "\$OUT"
    """
}


process MergeLocusCoverage {
    tag "Merge locus coverage: ${query_id}"
    errorStrategy 'ignore'
    publishDir { "${params.outdir}/locus_coverage/${query_id}" }, mode: 'copy'

    input:
    tuple val(query_id), path(tsvs)

    output:
    tuple val(query_id), path("*_all_locus_coverage.tsv"), emit: merged

    script:
    def qid_safe = query_id.replaceAll(/[^A-Za-z0-9._-]/, '_')
    """
    set -euo pipefail

    files=( ${tsvs.join(' ')} )

    # Header from first file; data rows (skip header) from all files, sorted
    head -1 "\${files[0]}" > ${qid_safe}_all_locus_coverage.tsv
    for f in "\${files[@]}"; do
        tail -n +2 "\$f"
    done | sort -k1,1 >> ${qid_safe}_all_locus_coverage.tsv
    """
}
