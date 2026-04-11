// GenomeThreader spliced alignment of a protein or cDNA query against a
// de novo assembled genome. Handles multi-exon eukaryotic gene structures.
//
// GTHsearch     -- per-sample: run gth, extract CDS + protein via extract_gth_cds.py
// CombineGTHHits -- per-query: merge all sample outputs into combined FASTA + TSV

// ---------------------------------------------------------------------------
// GTHsearch
// ---------------------------------------------------------------------------
process GTHsearch {
    tag "GTH: ${sample_id} vs ${query_id}"
    errorStrategy 'ignore'
    memory '8GB'
    cpus 1
    publishDir "${params.outdir}/${query_id}/gth_results", mode: 'copy'

    input:
    tuple val(sample_id), path(scaffolds_fasta), val(query_id), path(query_fasta), val(seq_type)

    output:
    tuple val(query_id), val(sample_id), path("${sample_id}_gth_cds.fasta"),     emit: cds_sequences
    tuple val(query_id), val(sample_id), path("${sample_id}_gth_protein.fasta"), emit: protein_sequences
    tuple val(query_id), path("${sample_id}_gth_summary.tsv"),                   emit: summary

    script:
    // Use -protein for protein queries (tblastn equivalent),
    // -cdna for nucleotide/cDNA queries (blastn equivalent)
    def gth_mode = (seq_type == 'protein') ? '-protein' : '-cdna'
    """
    #!/bin/bash
    set -euo pipefail

    # Initialise empty output files (written if GTH finds no gene models)
    touch ${sample_id}_gth_cds.fasta
    touch ${sample_id}_gth_protein.fasta
    printf 'Sample_ID\tQuery_ID\tScaffold\tStrand\tExon_Count\tCDS_Length\tProtein_Length\tExon_Coordinates\n' \
        > ${sample_id}_gth_summary.tsv

    # Run GenomeThreader: protein-to-genome or cDNA-to-genome spliced alignment
    # -gff3out      : output GFF3 (exon/CDS coordinates)
    # -skipalignmentout : suppress verbose alignment text
    # -force        : overwrite output if it already exists
    gth -genomic ${scaffolds_fasta} ${gth_mode} ${query_fasta} \
        -gff3out -skipalignmentout \
        > ${sample_id}_gth.gff3 \
        2> ${sample_id}_gth.log || true

    # Extract CDS exons, concatenate respecting strand, and translate to protein.
    # extract_gth_cds.py is staged from the pipeline bin/ directory.
    if grep -qP '\tCDS\t' ${sample_id}_gth.gff3 2>/dev/null; then
        QUERY_ID='${query_id}'
        python3 extract_gth_cds.py \
            ${scaffolds_fasta} \
            ${sample_id}_gth.gff3 \
            ${sample_id} \
            "\${QUERY_ID}"
    else
        echo "No CDS features found in GTH output for ${sample_id} vs ${query_id}" >> ${sample_id}_gth.log
    fi
    """
}

// ---------------------------------------------------------------------------
// CombineGTHHits
// ---------------------------------------------------------------------------
process CombineGTHHits {
    tag "Combining GTH hits for ${query_id}"
    publishDir "${params.outdir}/${query_id}", mode: 'copy'

    input:
    // query_fasta    : original query file (to prepend to combined output)
    // cds_fastas     : list of per-sample *_gth_cds.fasta files
    // protein_fastas : list of per-sample *_gth_protein.fasta files
    // summary_tsvs   : list of per-sample *_gth_summary.tsv files
    tuple val(query_id), path(query_fasta), path(cds_fastas), path(protein_fastas), path(summary_tsvs)

    output:
    // combined_cds     : for AlignMAFFT (nucleotide); query prepended if nucleotide query
    // combined_protein : for AlignProteinMAFFT; query prepended if protein query
    // summary          : merged per-sample hit table for this query
    tuple val(query_id), path("*_all_gth_cds.fasta"),     emit: combined_cds
    tuple val(query_id), path("*_all_gth_protein.fasta"), emit: combined_protein
    tuple val(query_id), path("*_all_gth_summary.tsv"),   emit: summary

    script:
    """
    #!/bin/bash
    set -euo pipefail

    clean_query_id="\$(echo '${query_id}' | sed 's|[^a-zA-Z0-9_.-]|_|g')"

    # -------------------------------------------------------------------------
    # Detect whether the original query is nucleotide or protein
    # (so we prepend it to the appropriate combined FASTA for alignment)
    # -------------------------------------------------------------------------
    query_type="unknown"
    if [ -f "${query_fasta}" ]; then
        seq_chars=\$(grep -v '^>' "${query_fasta}" | tr -d '\n[:space:]' | tr '[:lower:]' '[:upper:]')
        nuc_count=\$(echo "\${seq_chars}" | grep -o '[ATGCUN]' | wc -l)
        total_count=\$(echo "\${seq_chars}" | wc -c)
        if [ "\${total_count}" -gt 0 ]; then
            pct=\$((nuc_count * 100 / total_count))
            [ "\${pct}" -ge 95 ] && query_type="nucleotide" || query_type="protein"
        fi
    fi

    # -------------------------------------------------------------------------
    # Combined CDS FASTA (nucleotide)
    # Prepend original query only when it is a nucleotide sequence so that
    # MAFFT can align all sequences (query + hit CDS) together.
    # -------------------------------------------------------------------------
    if [ "\${query_type}" = "nucleotide" ] && [ -f "${query_fasta}" ]; then
        cat "${query_fasta}" > "\${clean_query_id}_all_gth_cds.fasta"
        echo "" >> "\${clean_query_id}_all_gth_cds.fasta"
    else
        : > "\${clean_query_id}_all_gth_cds.fasta"
    fi

    for f in *_gth_cds.fasta; do
        [[ "\${f}" == "\${clean_query_id}_all_gth_cds.fasta" ]] && continue
        [ -s "\${f}" ] && cat "\${f}" >> "\${clean_query_id}_all_gth_cds.fasta" || true
    done

    # -------------------------------------------------------------------------
    # Combined protein FASTA
    # Prepend original query only when it is a protein sequence.
    # -------------------------------------------------------------------------
    if [ "\${query_type}" = "protein" ] && [ -f "${query_fasta}" ]; then
        cat "${query_fasta}" > "\${clean_query_id}_all_gth_protein.fasta"
        echo "" >> "\${clean_query_id}_all_gth_protein.fasta"
    else
        : > "\${clean_query_id}_all_gth_protein.fasta"
    fi

    for f in *_gth_protein.fasta; do
        [[ "\${f}" == "\${clean_query_id}_all_gth_protein.fasta" ]] && continue
        [ -s "\${f}" ] && cat "\${f}" >> "\${clean_query_id}_all_gth_protein.fasta" || true
    done

    # -------------------------------------------------------------------------
    # Combined summary TSV — merge all per-sample summaries, add no-hit rows
    # -------------------------------------------------------------------------
    first_tsv=\$(ls *_gth_summary.tsv 2>/dev/null | grep -v '_all_gth_summary' | head -n 1 || echo "")

    if [ -n "\${first_tsv}" ]; then
        head -n 1 "\${first_tsv}" > "\${clean_query_id}_all_gth_summary.tsv"
        for tsv in *_gth_summary.tsv; do
            [[ "\${tsv}" == "\${clean_query_id}_all_gth_summary.tsv" ]] && continue
            sample_id=\$(basename "\${tsv}" "_gth_summary.tsv")
            line_count=\$(wc -l < "\${tsv}")
            if [ "\${line_count}" -gt 1 ]; then
                tail -n +2 "\${tsv}" >> "\${clean_query_id}_all_gth_summary.tsv"
            else
                echo -e "\${sample_id}\t${query_id}\tno_hits_found\t.\t0\t0\t0\t." \
                    >> "\${clean_query_id}_all_gth_summary.tsv"
            fi
        done
    else
        printf 'Sample_ID\tQuery_ID\tScaffold\tStrand\tExon_Count\tCDS_Length\tProtein_Length\tExon_Coordinates\n' \
            > "\${clean_query_id}_all_gth_summary.tsv"
    fi
    """
}
