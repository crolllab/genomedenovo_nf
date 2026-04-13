// GenomeThreader spliced alignment of a protein or cDNA query against a
// de novo assembled genome. Handles multi-exon eukaryotic gene structures.
//
// GTHsearch      -- per-sample: run gth, output raw GFF3 (genomethreader container, no python3)
// ExtractGTHCds  -- per-sample: parse GFF3, extract CDS + translate protein (spades container, has python3)
// CombineGTHHits -- per-query: merge all sample outputs into combined FASTA + TSV

// ---------------------------------------------------------------------------
// GTHsearch  (genomethreader container — no python3, BusyBox shell)
// ---------------------------------------------------------------------------
process GTHsearch {
    tag "GTH search in assemblies"
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${query_id}/gth_results", mode: 'copy', pattern: "*_gth.gff3"

    input:
    tuple val(sample_id), path(scaffolds_fasta), val(query_id), path(query_fasta), val(seq_type)

    output:
    // Pass scaffolds + GFF3 forward to ExtractGTHCds for Python-based CDS extraction
    tuple val(query_id), val(sample_id), path(scaffolds_fasta), path("${sample_id}_gth.gff3"), emit: gff3

    script:
    def gth_mode = (seq_type == 'protein') ? '-protein' : '-cdna'
    """
    #!/bin/bash
    set -euo pipefail

    # Run GenomeThreader: protein-to-genome or cDNA-to-genome spliced alignment
    # -gff3out      : output GFF3 (exon/CDS coordinates)
    # -skipalignmentout : suppress verbose alignment text
    gth -genomic ${scaffolds_fasta} ${gth_mode} ${query_fasta} \
        -gff3out -skipalignmentout \
        > "${sample_id}_gth.gff3" \
        2> "${sample_id}_gth.log" || true
    """
}

// ---------------------------------------------------------------------------
// ExtractGTHCds  (spades container — has python3)
// Parses the GTH GFF3 and extracts CDS + protein sequences per sample.
// ---------------------------------------------------------------------------
process ExtractGTHCds {
    tag "Extract GTH CDS"
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${query_id}/gth_results", mode: 'copy'

    input:
    tuple val(query_id), val(sample_id), path(scaffolds_fasta), path(gff3_file)

    output:
    tuple val(query_id), val(sample_id), path("${sample_id}_gth_cds.fasta"),     emit: cds_sequences
    tuple val(query_id), val(sample_id), path("${sample_id}_gth_protein.fasta"), emit: protein_sequences
    tuple val(query_id), path("${sample_id}_gth_summary.tsv"),                   emit: summary

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # Initialise empty output files
    touch "${sample_id}_gth_cds.fasta"
    touch "${sample_id}_gth_protein.fasta"
    printf 'Sample_ID\tQuery_ID\tScaffold\tStrand\tExon_Count\tCDS_Length\tProtein_Length\tExon_Coordinates\n' \
        > "${sample_id}_gth_summary.tsv"

    # Use grep -F with a literal tab (works in BusyBox and GNU grep)
    TAB=\$(printf '\t')
    if grep -qF "\${TAB}CDS\${TAB}" ${gff3_file} 2>/dev/null; then
        QUERY_ID='${query_id}'
        extract_gth_cds.py \
            ${scaffolds_fasta} \
            ${gff3_file} \
            "${sample_id}" \
            "\${QUERY_ID}"
    else
        echo "No CDS features found in GTH output for ${sample_id}" >> "${sample_id}_gth_extract.log"
    fi
    """
}

// ---------------------------------------------------------------------------
// CombineGTHHits
// ---------------------------------------------------------------------------
process CombineGTHHits {
    tag "Combining GTH hits"
    publishDir "${params.outdir}/${query_id}", mode: 'copy'

    input:
    tuple val(query_id), path(query_fasta), path(cds_fastas), path(protein_fastas), path(summary_tsvs)

    output:
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
    # If query is nucleotide, translate it and prepend to the protein FASTA so
    # the protein alignment includes the reference sequence in protein form.
    # -------------------------------------------------------------------------
    if [ "\${query_type}" = "nucleotide" ] && [ -f "${query_fasta}" ]; then
        python3 - "${query_fasta}" "\${clean_query_id}_all_gth_protein.fasta" << 'PYEOF'
import sys, os

CODON = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf).upper()
            name = line[1:].split()[0]
            buf = []
        else:
            buf.append(line)
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs

def translate(seq):
    prot = []
    for i in range(0, len(seq) - 2, 3):
        aa = CODON.get(seq[i:i+3], "X")
        if aa == "*":
            break
        prot.append(aa)
    return "".join(prot)

fasta_path, out_path = sys.argv[1], sys.argv[2]
seqs = read_fasta(fasta_path)
tmp = out_path + ".tmp"
with open(tmp, "w") as out:
    for name, seq in seqs.items():
        prot = translate(seq)
        out.write(">" + name + "_query_prot\\n")
        for i in range(0, len(prot), 60):
            out.write(prot[i:i+60] + "\\n")
    if os.path.exists(out_path):
        out.write(open(out_path).read())
os.replace(tmp, out_path)
PYEOF
    fi

    # -------------------------------------------------------------------------
    # Combined summary TSV
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
