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
    tag "BLAST locus on reference"
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6
    publishDir { "${params.outdir}/11_locus_coverage/${query_id}" }, mode: 'copy'

    input:
    // paralog_fasta: the query's paralogs (--query_paralog), or [] without
    tuple val(query_id), path(query_fasta), val(seq_type), path(paralog_fasta), path(reference)

    output:
    tuple val(query_id), path("*_ref_locus.tsv"), emit: locus

    script:
    def qid_safe   = query_id.replaceAll(/[^A-Za-z0-9._-]/, '_')
    def blast_prog  = (seq_type == 'protein') ? 'tblastn' : 'blastn'
    def blast_extra = (seq_type == 'protein') ? '' : '-task blastn'
    // --query_paralog: the locus must belong to the query. If a paralog hit
    // overlapping it scores higher, the reference has no copy of the query
    // there and no locus is reported (coverage NA in every sample).
    def paralog_check = paralog_fasta ? """
    if [ -s blast_locus.txt ]; then
        ${blast_prog} -query ${paralog_fasta} -db ref_db ${blast_extra} \\
            -outfmt '6 sseqid sstart send sstrand evalue bitscore' \\
            -num_threads ${task.cpus} \\
            -out paralog_locus.txt
        if awk -F'\\t' '
            NR == FNR { if (FNR == 1) { c = \$1; s = (\$2 < \$3) ? \$2 : \$3; e = (\$2 < \$3) ? \$3 : \$2; b = \$6 }; next }
            \$1 == c && ((\$2 < \$3) ? \$2 : \$3) <= e && ((\$2 < \$3) ? \$3 : \$2) >= s && \$6 > b { f = 1 }
            END { exit !f }' blast_locus.txt paralog_locus.txt; then
            echo "WARNING: The best hit of ${query_id} on the reference belongs to a paralog; no locus" >&2
            : > blast_locus.txt
        fi
    fi
""" : ''
    """
    set -euo pipefail

    makeblastdb -in ${reference} -dbtype nucl -out ref_db -logfile /dev/null

    ${blast_prog} -query ${query_fasta} -db ref_db ${blast_extra} \\
        -outfmt '6 sseqid sstart send sstrand evalue bitscore' \\
        -max_target_seqs 1 -max_hsps 1 \\
        -num_threads ${task.cpus} \\
        -out blast_locus.txt
${paralog_check}
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
    tag "Locus coverage"
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6
    publishDir { "${params.outdir}/11_locus_coverage/${query_id}" }, mode: 'copy'

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
    printf 'genome_mean_depth\\tgenome_sd_depth\\thit_mean_depth\\thit_breadth_pct\\n' >> "\$OUT"

    # No BLAST hit for this query on the reference — write NA row and exit
    if [ -z "\$CHROM" ]; then
        printf '%s\\t%s\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\n' \\
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

    # ── Depth and breadth (% of bases with >= 1 read) over the unpadded hit ──
    read -r HIT_MEAN HIT_BREADTH < <(awk -v s="\$HIT_START" -v e="\$HIT_END" '
        \$2 >= s && \$2 <= e { t += \$3; n++; if (\$3 > 0) c++ }
        END { if (n > 0) printf "%.2f %.1f\\n", t / n, 100 * c / n; else print "NA NA" }
    ' locus_depth.txt)

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
    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \\
        "${sample_id}" "${query_id}" \\
        "\$CHROM" "\$HIT_START" "\$HIT_END" \\
        "\$LOCUS_START" "\$LOCUS_END" \\
        "\$LOCUS_MEAN" "\$LOCUS_MEDIAN" \\
        "\$GENOME_MEAN" "\$GENOME_SD" \\
        "\$HIT_MEAN" "\$HIT_BREADTH" \\
        >> "\$OUT"
    """
}


process MergeLocusCoverage {
    tag "Merge locus coverage"
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6
    publishDir { "${params.outdir}/11_locus_coverage/${query_id}" }, mode: 'copy'

    input:
    // blast_summary: the query's BLAST hit summary (CombineBlastHits), for the
    // assembly status of every sample
    tuple val(query_id), path(tsvs), path(blast_summary, stageAs: 'blast/*')

    output:
    tuple val(query_id), path("*_all_locus_coverage.tsv"), emit: merged

    script:
    def qid_safe = query_id.replaceAll(/[^A-Za-z0-9._-]/, '_')
    """
    set -euo pipefail

    files=( ${tsvs.join(' ')} )

    # Coverage conflict: the assembly has no copy of the query (no BLAST hit,
    # or with --query_paralog only paralog hits) while the reads cover the
    # query's locus on the reference, e.g. an allele that failed to assemble
    head -1 "\${files[0]}" | awk '{ print \$0 "\\thit_to_genome_depth_ratio\\tassembly_hit\\tcoverage_conflict" }' \
        > ${qid_safe}_all_locus_coverage.tsv
    for f in "\${files[@]}"; do
        tail -n +2 "\$f"
    done | sort -k1,1 | awk -F'\\t' -v OFS='\\t' \
        -v min_ratio=${params.locus_conflict_min_depth_ratio} -v min_breadth=${params.locus_conflict_min_breadth} '
        NR == FNR {
            if (FNR > 1 && st[\$1] != "hit")
                st[\$1] = (\$3 == "paralog_only") ? "paralog_only" : (\$3 == "no_hits_found") ? "no_hit" : "hit"
            next
        }
        {
            # 10 genome_mean_depth, 12 hit_mean_depth, 13 hit_breadth_pct
            ratio = (\$10 + 0 > 0 && \$12 != "NA") ? sprintf("%.3f", \$12 / \$10) : "NA"
            hit = (\$1 in st) ? st[\$1] : "NA"
            conflict = "NA"
            if (ratio != "NA" && hit != "NA") {
                conflict = "none"
                if (hit != "hit" && ratio + 0 >= min_ratio && \$13 + 0 >= min_breadth)
                    conflict = hit "_locus_covered"
            }
            print \$0, ratio, hit, conflict
        }' ${blast_summary} - >> ${qid_safe}_all_locus_coverage.tsv
    """
}
