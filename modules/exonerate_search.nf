// Exonerate spliced alignment of a protein or CDS query against a
// de novo assembled genome. Handles multi-exon eukaryotic gene structures.
//
// ExonerateSearch   -- per-sample: run exonerate, output target GFF
// ExtractExonerateCds  -- per-sample: parse GFF, extract CDS + translate protein (spades container, has python3)
// CombineExonerateHits -- per-query: merge all sample outputs into combined FASTA + TSV

// ---------------------------------------------------------------------------
// ExonerateSearch  (exonerate container)
// ---------------------------------------------------------------------------
process ExonerateSearch {
    tag "Exonerate search in assemblies"
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6
    publishDir { "${params.outdir}/7_exonerate_hits/${query_id}/gff3_results" }, mode: 'copy', pattern: "*_exonerate.gff3"

    input:
    tuple val(sample_id), path(scaffolds_fasta), val(query_id), path(query_fasta), val(seq_type)

    output:
    // Pass scaffolds + GFF forward to ExtractExonerateCds for Python-based CDS extraction
    tuple val(query_id), val(sample_id), path(scaffolds_fasta), path("${sample_id}_exonerate.gff3"), emit: gff3

    script:
    // Nucleotide queries are mature CDS, so use coding2genome (codon-aware, spliced).
    // est2genome labels the aligned region as UTR/exon and emits no cds features,
    // which leaves nothing for ExtractExonerateCds to extract.
    def exonerate_model = (seq_type == 'protein') ? 'protein2genome' : 'coding2genome'
    // --refine region re-aligns the heuristic hit region with full dynamic
    // programming; without it divergent alleles often lose whole exons and
    // look truncated. (--exhaustive hangs in this exonerate build.)
    def refine = (params.exonerate_refine in ['region', 'full'])
        ? "--refine ${params.exonerate_refine} --refineboundary 1000" : ''
    """
    #!/bin/bash
    set -euo pipefail

    # Run Exonerate: protein-to-genome or CDS-to-genome spliced alignment.
    # Use target-side GFF output for downstream CDS extraction.
    exonerate \
        --model ${exonerate_model} ${refine} \
        --query ${query_fasta} \
        --target ${scaffolds_fasta} \
        --showalignment no \
        --showvulgar yes \
        --showsugar no \
        --showquerygff no \
        --showtargetgff yes \
        > "${sample_id}_exonerate.gff3" \
        2> "${sample_id}_exonerate.log" || true
    """
}

// ---------------------------------------------------------------------------
// ExtractExonerateCds  (spades container — has python3)
// Parses the Exonerate GFF and extracts CDS + protein sequences per sample.
// ---------------------------------------------------------------------------
process ExtractExonerateCds {
    tag "Extract Exonerate CDS"
    errorStrategy { task.attempt <= 6 ? 'retry' : 'ignore' }
    maxRetries 6
    publishDir { "${params.outdir}/7_exonerate_hits/${query_id}/hit_sequences" }, mode: 'copy'

    input:
    tuple val(query_id), val(sample_id), path(scaffolds_fasta), path(gff3_file)

    output:
    tuple val(query_id), val(sample_id), path("${sample_id}_exonerate_cds.fasta"),     emit: cds_sequences
    tuple val(query_id), val(sample_id), path("${sample_id}_exonerate_protein.fasta"), emit: protein_sequences
    tuple val(query_id), path("${sample_id}_exonerate_summary.tsv"),                   emit: summary

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # Initialise empty output files
    touch "${sample_id}_exonerate_cds.fasta"
    touch "${sample_id}_exonerate_protein.fasta"
    printf 'sample_id\tquery_id\tsubject_id\tstrand\tquery_start\tquery_end\tsubject_start\tsubject_end\texon_count\tcds_length\tprotein_length\texon_coordinates\n' \
        > "${sample_id}_exonerate_summary.tsv"

    # Exonerate can emit CDS/cds; detect case-insensitively.
    if awk -F'\t' 'NF>=3 && tolower(\$3)=="cds" {found=1; exit} END{exit !found}' ${gff3_file}; then
        QUERY_ID='${query_id}'
        extract_exonerate_cds.py \
            ${scaffolds_fasta} \
            ${gff3_file} \
            "${sample_id}" \
            "\${QUERY_ID}"
    else
        echo "No CDS features found in Exonerate output for ${sample_id}" >> "${sample_id}_exonerate_extract.log"
    fi
    """
}

// ---------------------------------------------------------------------------
// CombineExonerateHits
// ---------------------------------------------------------------------------
process CombineExonerateHits {
    tag "Combining Exonerate hits"
    publishDir { "${params.outdir}/7_exonerate_hits/${query_id}" }, mode: 'copy'

    input:
    tuple val(query_id), path(query_fasta), path(cds_fastas), path(protein_fastas), path(summary_tsvs)

    output:
    tuple val(query_id), path("*_all_exonerate_cds.fasta"),     emit: combined_cds
    tuple val(query_id), path("*_all_exonerate_protein.fasta"), emit: combined_protein
    tuple val(query_id), path("*_all_exonerate_summary.tsv"),   emit: summary

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
        cat "${query_fasta}" > "\${clean_query_id}_all_exonerate_cds.fasta"
        echo "" >> "\${clean_query_id}_all_exonerate_cds.fasta"
    else
        : > "\${clean_query_id}_all_exonerate_cds.fasta"
    fi

    for f in *_exonerate_cds.fasta; do
        [[ "\${f}" == "\${clean_query_id}_all_exonerate_cds.fasta" ]] && continue
        [ -s "\${f}" ] && cat "\${f}" >> "\${clean_query_id}_all_exonerate_cds.fasta" || true
    done

    # -------------------------------------------------------------------------
    # Combined protein FASTA
    # -------------------------------------------------------------------------
    if [ "\${query_type}" = "protein" ] && [ -f "${query_fasta}" ]; then
        cat "${query_fasta}" > "\${clean_query_id}_all_exonerate_protein.fasta"
        echo "" >> "\${clean_query_id}_all_exonerate_protein.fasta"
    else
        : > "\${clean_query_id}_all_exonerate_protein.fasta"
    fi

    for f in *_exonerate_protein.fasta; do
        [[ "\${f}" == "\${clean_query_id}_all_exonerate_protein.fasta" ]] && continue
        [ -s "\${f}" ] && cat "\${f}" >> "\${clean_query_id}_all_exonerate_protein.fasta" || true
    done

    # -------------------------------------------------------------------------
    # If query is nucleotide, translate it and prepend to the protein FASTA so
    # the protein alignment includes the reference sequence in protein form.
    # -------------------------------------------------------------------------
    if [ "\${query_type}" = "nucleotide" ] && [ -f "${query_fasta}" ]; then
        python3 - "${query_fasta}" "\${clean_query_id}_all_exonerate_protein.fasta" << 'PYEOF'
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
    first_tsv=\$(ls *_exonerate_summary.tsv 2>/dev/null | grep -v '_all_exonerate_summary' | head -n 1 || echo "")

    if [ -n "\${first_tsv}" ]; then
        head -n 1 "\${first_tsv}" > "\${clean_query_id}_all_exonerate_summary.tsv"
        for tsv in *_exonerate_summary.tsv; do
            [[ "\${tsv}" == "\${clean_query_id}_all_exonerate_summary.tsv" ]] && continue
            sample_id=\$(basename "\${tsv}" "_exonerate_summary.tsv")
            line_count=\$(wc -l < "\${tsv}")
            if [ "\${line_count}" -gt 1 ]; then
                tail -n +2 "\${tsv}" >> "\${clean_query_id}_all_exonerate_summary.tsv"
            else
                echo -e "\${sample_id}\t${query_id}\tno_hits_found\t.\tNA\tNA\tNA\tNA\t0\t0\t0\t." \
                    >> "\${clean_query_id}_all_exonerate_summary.tsv"
            fi
        done
    else
        printf 'sample_id\tquery_id\tsubject_id\tstrand\tquery_start\tquery_end\tsubject_start\tsubject_end\texon_count\tcds_length\tprotein_length\texon_coordinates\n' \
            > "\${clean_query_id}_all_exonerate_summary.tsv"
    fi
    """
}
