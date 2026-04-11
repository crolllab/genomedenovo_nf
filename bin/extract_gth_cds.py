#!/usr/bin/env python3
"""
Extract CDS exon sequences from GenomeThreader GFF3 output, concatenate
them per gene model (respecting strand), and translate to protein.

Usage:
    python3 extract_gth_cds.py <genome.fasta> <gth.gff3> <sample_id> <query_id>

Outputs (in current directory):
    <sample_id>_gth_cds.fasta      -- nucleotide CDS (exons concatenated)
    <sample_id>_gth_protein.fasta  -- translated protein
    <sample_id>_gth_summary.tsv    -- per-gene-model statistics

Strand handling:
    + strand: sort exons ascending by start, extract forward sequence, concatenate.
    - strand: sort exons descending by start (5'->3' of transcript), reverse-complement
              each exon individually, then concatenate. This correctly assembles the
              spliced mRNA regardless of the number of exons.
"""

import sys
from collections import defaultdict


# ---------------------------------------------------------------------------
# FASTA reader
# ---------------------------------------------------------------------------

def read_fasta(path):
    """Return {seq_id: uppercase_sequence} for all sequences in *path*."""
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if name:
                    seqs[name] = ''.join(buf).upper()
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if name:
        seqs[name] = ''.join(buf).upper()
    return seqs


# ---------------------------------------------------------------------------
# Reverse complement
# ---------------------------------------------------------------------------

_RC_TABLE = str.maketrans('ATGCatgcNnRrYySsWwKkMmBbDdHhVv',
                           'TACGtacgNnYyRrSsWwMmKkVvHhDdBb')

def rc(seq):
    return seq.translate(_RC_TABLE)[::-1]


# ---------------------------------------------------------------------------
# Standard genetic code (NCBI table 1)
# ---------------------------------------------------------------------------

CODON = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def translate(seq):
    """Translate nucleotide *seq* to protein; stop at first stop codon."""
    prot = []
    for i in range(0, len(seq) - 2, 3):
        aa = CODON.get(seq[i:i + 3].upper(), 'X')
        if aa == '*':
            break
        prot.append(aa)
    return ''.join(prot)


def wrap(seq, w=60):
    """Wrap *seq* to lines of width *w*."""
    return '\n'.join(seq[i:i + w] for i in range(0, len(seq), w))


# ---------------------------------------------------------------------------
# GFF3 parser
# ---------------------------------------------------------------------------

def parse_gff3_cds(gff_path):
    """
    Return (cds_by_parent, parent_strand, parent_scaffold).

    cds_by_parent : {parent_id: [(start_1based, end_1based, scaffold), ...]}
    parent_strand  : {parent_id: '+' or '-'}
    parent_scaffold: {parent_id: scaffold_name}
    """
    cds_by_parent   = defaultdict(list)
    parent_strand   = {}
    parent_scaffold = {}

    with open(gff_path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            seqid, _, feat, start, end, _, strand, _, attrs = parts
            if feat != 'CDS':
                continue

            parent = None
            for item in attrs.split(';'):
                item = item.strip()
                if item.upper().startswith('PARENT='):
                    parent = item.split('=', 1)[1].split(',')[0].strip()
                    break
            if parent:
                cds_by_parent[parent].append((int(start), int(end), seqid))
                parent_strand[parent]   = strand
                parent_scaffold[parent] = seqid

    return cds_by_parent, parent_strand, parent_scaffold


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) != 5:
        sys.exit(f'Usage: {sys.argv[0]} <genome.fasta> <gth.gff3> <sample_id> <query_id>')

    genome_file, gff_file, sample_id, query_id_val = sys.argv[1:]

    genome = read_fasta(genome_file)
    cds_by_parent, parent_strand, parent_scaffold = parse_gff3_cds(gff_file)

    if not cds_by_parent:
        sys.stderr.write(f'No CDS features found in {gff_file}\n')
        return

    with (open(f'{sample_id}_gth_cds.fasta', 'w')     as cds_out,
          open(f'{sample_id}_gth_protein.fasta', 'w')  as prot_out,
          open(f'{sample_id}_gth_summary.tsv', 'w')    as sum_out):

        sum_out.write(
            'Sample_ID\tQuery_ID\tScaffold\tStrand\t'
            'Exon_Count\tCDS_Length\tProtein_Length\tExon_Coordinates\n'
        )

        for parent_id, exons in sorted(cds_by_parent.items()):
            strand   = parent_strand[parent_id]
            scaffold = parent_scaffold[parent_id]

            if scaffold not in genome:
                sys.stderr.write(f'WARNING: scaffold {scaffold!r} not in genome FASTA — skipping {parent_id}\n')
                continue

            # ----------------------------------------------------------------
            # Sort exons so they are in transcript 5'->3' order:
            #   + strand: ascending genomic coordinate
            #   - strand: descending genomic coordinate (highest coord = 5' end)
            # ----------------------------------------------------------------
            exons_sorted = sorted(exons, key=lambda x: x[0], reverse=(strand == '-'))

            full_cds = ''
            coords   = []
            for start, end, scaf in exons_sorted:
                # GFF3: 1-based inclusive -> Python slice [start-1 : end]
                exon_seq = genome[scaf][start - 1 : end]
                if strand == '-':
                    # RC each exon individually so concatenation gives
                    # the correct 5'->3' transcript sequence
                    exon_seq = rc(exon_seq)
                full_cds += exon_seq
                coords.append(f'{scaf}:{start}-{end}')

            if not full_cds:
                continue

            prot   = translate(full_cds)
            header = f'{sample_id}_{parent_id}'

            cds_out.write(f'>{header}\n{wrap(full_cds)}\n')
            prot_out.write(f'>{header}_prot\n{wrap(prot)}\n')
            sum_out.write(
                f'{sample_id}\t{query_id_val}\t{scaffold}\t{strand}\t'
                f'{len(exons_sorted)}\t{len(full_cds)}\t{len(prot)}\t'
                f'{";".join(coords)}\n'
            )


if __name__ == '__main__':
    main()
