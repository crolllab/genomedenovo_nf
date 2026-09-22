#!/usr/bin/env python3
"""
Assign consecutive isoform IDs to unique CDS and protein haplotypes
derived from MAFFT-aligned FASTA files.

Usage:
    summarize_isoforms.py <query_id> <cds_aligned.fasta> <protein_aligned.fasta>

Output:
    <safe_query_id>_isoforms.tsv              (columns: sample_id, cds_isoform, protein_isoform)
    <safe_query_id>_cds_isoforms.fasta        one ungapped sequence per unique CDS isoform
    <safe_query_id>_protein_isoforms.fasta    one ungapped sequence per unique protein isoform

Haplotype identity is determined from the full aligned sequence (gaps included),
uppercased. MAFFT reverse-complement headers prefixed with '_R_' are stripped.
The reference entry in the protein alignment (any header not ending in '_prot')
is skipped automatically.
"""

import sys
import re


def read_fasta(path):
    """Return list of (header, seq_uppercase) preserving record order.

    Strips the leading '_R_' MAFFT reverse-complement flag from headers.
    """
    records = []
    header = None
    parts = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n\r')
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(parts).upper()))
                raw = line[1:].strip()
                header = raw[3:] if raw.startswith('_R_') else raw
                parts = []
            else:
                parts.append(line)
    if header is not None:
        records.append((header, ''.join(parts).upper()))
    return records


def assign_isoforms(records):
    """Map each record to a consecutive isoform integer (first-seen ordering).

    Returns {header: isoform_int}.
    """
    haplotype_to_num = {}
    result = {}
    for header, seq in records:
        if seq not in haplotype_to_num:
            haplotype_to_num[seq] = len(haplotype_to_num) + 1
        result[header] = haplotype_to_num[seq]
    return result


def write_isoform_fasta(path, prefix, records, isoforms):
    """Write one ungapped sequence per isoform, numbered as in the TSV.

    The header records how many samples carry the isoform, e.g.
    '>cds_isoform_1 n_samples=3'.
    """
    seqs, counts = {}, {}
    for header, seq in records:
        num = isoforms[header]
        seqs.setdefault(num, seq.replace('-', ''))
        counts[num] = counts.get(num, 0) + 1
    with open(path, 'w') as out:
        for num in sorted(seqs):
            seq = seqs[num]
            out.write('>{}_{} n_samples={}\n'.format(prefix, num, counts[num]))
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + '\n')


def main():
    if len(sys.argv) != 4:
        sys.exit("Usage: summarize_isoforms.py <query_id> <cds_fasta> <protein_fasta>")

    query_id = sys.argv[1]
    cds_path = sys.argv[2]
    prot_path = sys.argv[3]

    # --- CDS alignment ---------------------------------------------------------
    cds_records = read_fasta(cds_path)
    cds_isoforms = assign_isoforms(cds_records)

    # --- Protein alignment -----------------------------------------------------
    # Keep only records whose header ends with '_prot' (sample records);
    # the reference entry is included by CombineExonerateHits but should be skipped.
    prot_records = [(h, s) for h, s in read_fasta(prot_path) if h.endswith('_prot')]
    prot_isoforms_raw = assign_isoforms(prot_records)
    # Re-key by stripping '_prot' suffix so keys match sample_id
    prot_isoforms = {h[:-5]: v for h, v in prot_isoforms_raw.items()}

    # --- Write output ----------------------------------------------------------
    safe_qid = re.sub(r'[^\w.\-]', '_', query_id)
    out_path = '{}_isoforms.tsv'.format(safe_qid)
    with open(out_path, 'w') as out:
        out.write('sample_id\tcds_isoform\tprotein_isoform\n')
        for sample_id in cds_isoforms:
            cds_iso = cds_isoforms[sample_id]
            prot_iso = prot_isoforms.get(sample_id, 'NA')
            out.write('{}\t{}\t{}\n'.format(sample_id, cds_iso, prot_iso))

    write_isoform_fasta('{}_cds_isoforms.fasta'.format(safe_qid),
                        'cds_isoform', cds_records, cds_isoforms)
    write_isoform_fasta('{}_protein_isoforms.fasta'.format(safe_qid),
                        'protein_isoform', prot_records, prot_isoforms_raw)


if __name__ == '__main__':
    main()
