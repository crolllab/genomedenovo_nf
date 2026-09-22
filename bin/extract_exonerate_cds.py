#!/usr/bin/env python3
"""
Extract CDS exon sequences from GenomeThreader/Exonerate GFF output,
concatenate them per gene model (respecting strand), and translate to protein.

Usage:
    python3 extract_exonerate_cds.py <genome.fasta> <spliced_aligner.gff> <sample_id> <query_id>

Outputs (in current directory):
    <sample_id>_exonerate_cds.fasta      -- nucleotide CDS (exons concatenated)
    <sample_id>_exonerate_protein.fasta  -- translated protein
    <sample_id>_exonerate_summary.tsv    -- per-gene-model statistics

Strand handling:
    + strand: sort exons ascending by start, extract forward sequence, concatenate.
    - strand: sort exons descending by start (5'->3' of transcript), reverse-complement
              each exon individually, then concatenate. This correctly assembles the
              spliced mRNA regardless of the number of exons.
"""

import re
import sys
from collections import defaultdict


# ---------------------------------------------------------------------------
# FASTA reader
# ---------------------------------------------------------------------------

def read_fasta(path):
    """Return {seq_id: uppercase_sequence} for all sequences in *path*.

    Some spliced-alignment outputs use seqids without trailing coordinate suffixes,
    so
    strip any ':start-end' coordinate suffix (added by samtools faidx or
    the extract_blast_hits step) to ensure GFF3 lookups succeed.
    """
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if name:
                    seqs[name] = ''.join(buf).upper()
                raw_name = line[1:].split()[0]
                # Strip samtools-style coordinate suffix, e.g. "SCAFFOLD:1-428"
                name = re.sub(r':\d+-\d+$', '', raw_name)
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
# Vulgar line parser — query coordinate extraction
# ---------------------------------------------------------------------------

def parse_vulgar_query_coords(path):
    """Parse vulgar lines from exonerate output (--showvulgar yes).

    Returns {target_id: (query_start_1based, query_end_1based)}.
    When multiple alignments hit the same target, keeps the highest-scoring one.
    Vulgar format: vulgar: query_id q_start q_end q_strand target_id t_start
                           t_end t_strand score operations...
    Coordinates are 0-based half-open; converted to 1-based inclusive here.
    """
    coords = {}
    scores = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith('vulgar:'):
                continue
            parts = line.split()
            if len(parts) < 10:
                continue
            q_start   = int(parts[2])
            q_end     = int(parts[3])
            target_id = parts[5]
            score     = int(parts[9])
            # 0-based half-open → 1-based inclusive
            q1 = q_start + 1
            q2 = q_end
            if target_id not in scores or score > scores[target_id]:
                coords[target_id] = (q1, q2)
                scores[target_id] = score
    return coords


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
# Exonerate GFF2 → GFF3 converter
# ---------------------------------------------------------------------------

def _parse_gff2_attrs(attr_field):
    """Parse space/semicolon-separated GFF2 attributes (key value ; key value).

    Returns a plain dict.  Values are stripped of surrounding quotes.
    """
    attrs = {}
    for token in attr_field.split(';'):
        token = token.strip()
        if not token:
            continue
        m = re.match(r'^(\S+)\s+"?(.*?)"?\s*$', token)
        if m:
            attrs[m.group(1)] = m.group(2)
    return attrs


def gff2_to_gff3_lines(gff_path):
    """Convert Exonerate GFF2 output to a list of GFF3-compliant strings.

    Exonerate's ``--showtargetgff`` produces GFF version 2:
      - ``gene`` lines carry ``gene_id N ; sequence QUERY_ID ; ...``
      - ``cds`` lines carry *no attributes at all* (8 columns)
      - ``exon`` and ``similarity`` lines are also present but not needed

    This function assigns each ``gene`` a stable GFF3 ID, then emits a
    properly attributed ``CDS`` line for every attribute-less ``cds`` line
    that follows it.  All other feature types are dropped.
    """
    gff3_lines = []
    current_id = None
    current_strand = None
    n_genes = 0

    with open(gff_path) as fh:
        for raw in fh:
            if raw.startswith('#') or not raw.strip():
                continue
            cols = raw.rstrip('\n').split('\t')
            if len(cols) < 8:
                continue
            seqid, source, feat, start, end, score, strand, frame = cols[:8]
            attr_col = cols[8] if len(cols) >= 9 else ''

            if feat.lower() == 'gene':
                gff2_attrs = _parse_gff2_attrs(attr_col)
                gid = gff2_attrs.get('gene_id', '')
                seq = gff2_attrs.get('sequence', '')
                # Exonerate restarts gene_id at 1 for every alignment, so add a
                # running index to keep separate alignments on one target apart.
                n_genes += 1
                current_id = f'{seqid}_gene{gid}_{seq}' if seq else f'{seqid}_gene{gid}'
                current_id = f'{current_id}_{n_genes}'
                current_strand = strand
                gff3_lines.append(
                    f'{seqid}\t{source}\tgene\t{start}\t{end}\t{score}\t{strand}\t.\t'
                    f'ID={current_id}\n'
                )

            elif feat.lower() == 'cds' and current_id is not None:
                # cds lines may have 8 or 9 columns; strand may be '.'
                eff_strand = strand if strand in ('+', '-') else current_strand
                gff3_lines.append(
                    f'{seqid}\t{source}\tCDS\t{start}\t{end}\t.\t{eff_strand}\t.\t'
                    f'Parent={current_id}\n'
                )

    return gff3_lines


# ---------------------------------------------------------------------------
# GFF3 parser
# ---------------------------------------------------------------------------

def parse_gff3_cds(gff_path):
    """Return (cds_by_parent, parent_strand, parent_scaffold).

    cds_by_parent : {parent_id: [(start_1based, end_1based, scaffold), ...]}
    parent_strand  : {parent_id: '+' or '-'}
    parent_scaffold: {parent_id: scaffold_name}

    The input is first normalised to GFF3 via :func:`gff2_to_gff3_lines` so
    that Exonerate GFF2 and standard GFF3 are handled uniformly.
    """
    cds_by_parent   = defaultdict(list)
    parent_strand   = {}
    parent_scaffold = {}
    gene_span       = {}

    for line in gff2_to_gff3_lines(gff_path):
        if line.startswith('#') or not line.strip():
            continue
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 9:
            continue
        seqid, _, feat, start, end, score, strand, _, attrs = parts
        if feat == 'gene':
            gid = attrs.split('=', 1)[1]
            gene_span[gid] = (seqid, int(start), int(end), float(score) if score != '.' else 0.0)
            continue
        if feat.upper() != 'CDS':
            continue

        parent = None
        for kv in attrs.split(';'):
            kv = kv.strip()
            if kv.lower().startswith('parent='):
                parent = kv.split('=', 1)[1].strip()
                break

        if parent:
            cds_by_parent[parent].append((int(start), int(end), seqid))
            parent_strand[parent]   = strand
            parent_scaffold[parent] = seqid

    drop_overlapping_models(cds_by_parent, gene_span)
    return cds_by_parent, parent_strand, parent_scaffold


def drop_overlapping_models(cds_by_parent, gene_span):
    """Keep only the highest-scoring gene model among models that overlap on
    the same target sequence.

    coding2genome (CDS queries) also reports a weaker alignment of the
    reverse-complemented query over the same locus; without this filter it
    would be emitted as a second, spurious CDS/protein for the sample.
    """
    kept = []
    for gid in sorted(cds_by_parent, key=lambda g: -gene_span.get(g, ('', 0, 0, 0.0))[3]):
        if gid not in gene_span:
            kept.append(gid)
            continue
        seqid, start, end, _ = gene_span[gid]
        if any(gene_span[k][0] == seqid and start <= gene_span[k][2] and gene_span[k][1] <= end
               for k in kept if k in gene_span):
            del cds_by_parent[gid]
        else:
            kept.append(gid)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) != 5:
        sys.exit(f'Usage: {sys.argv[0]} <genome.fasta> <spliced_aligner.gff> <sample_id> <query_id>')

    genome_file, gff_file, sample_id, query_id_val = sys.argv[1:]

    genome = read_fasta(genome_file)
    cds_by_parent, parent_strand, parent_scaffold = parse_gff3_cds(gff_file)
    query_coords = parse_vulgar_query_coords(gff_file)

    if not cds_by_parent:
        sys.stderr.write(f'No CDS features found in {gff_file}\n')
        return

    with (open(f'{sample_id}_exonerate_cds.fasta', 'w')     as cds_out,
          open(f'{sample_id}_exonerate_protein.fasta', 'w')  as prot_out,
          open(f'{sample_id}_exonerate_summary.tsv', 'w')    as sum_out):

        sum_out.write(
            'sample_id\tquery_id\tsubject_id\tstrand\t'
            'query_start\tquery_end\t'
            'subject_start\tsubject_end\t'
            'exon_count\tcds_length\tprotein_length\texon_coordinates\n'
        )

        n_models = len(cds_by_parent)

        for model_idx, (parent_id, exons) in enumerate(sorted(cds_by_parent.items()), start=1):
            strand      = parent_strand[parent_id]
            scaffold_raw = parent_scaffold[parent_id]

            # Parse extraction offset encoded in the FASTA header by ExtractBlastHits.
            # Header format: "<contig>:<padded_start>-<padded_end>" (1-based, samtools style).
            # The offset lets us convert GFF window-relative coords to absolute assembly coords.
            m_off = re.match(r'^(.+):([0-9]+)-([0-9]+)$', scaffold_raw)
            if m_off:
                scaffold = m_off.group(1)          # real contig name
                offset   = int(m_off.group(2)) - 1 # 0-based offset into assembly
            else:
                scaffold = scaffold_raw
                offset   = 0

            # Query start/end from vulgar lines (keyed by the full scaffold_raw
            # which matches the exonerate target ID including :start-end suffix)
            q_start_val, q_end_val = query_coords.get(scaffold_raw, ('NA', 'NA'))

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
            abs_positions = []
            for start, end, scaf in exons_sorted:
                # GFF3: 1-based inclusive -> Python slice [start-1 : end]
                # Use the stripped scaffold name for genome lookup (all exons share the same target)
                exon_seq = genome[scaffold][start - 1 : end]
                if strand == '-':
                    # RC each exon individually so concatenation gives
                    # the correct 5'->3' transcript sequence
                    exon_seq = rc(exon_seq)
                full_cds += exon_seq
                # Convert window-relative (1-based) coords to absolute assembly coords
                abs_start = offset + start
                abs_end   = offset + end
                coords.append(f'{scaffold}:{abs_start}-{abs_end}')
                abs_positions.extend([abs_start, abs_end])

            if not full_cds:
                continue

            prot = translate(full_cds)
            # Name by sample_id for the common single-model case; add a
            # numeric suffix only when multiple gene models are reported.
            header = sample_id if n_models == 1 else f'{sample_id}_{model_idx}'

            # Report in original orientation: start > end for minus-strand hits (BLAST convention)
            if abs_positions:
                contig_start = max(abs_positions) if strand == '-' else min(abs_positions)
                contig_end   = min(abs_positions) if strand == '-' else max(abs_positions)
            else:
                contig_start = 'NA'
                contig_end   = 'NA'

            cds_out.write(f'>{header}\n{wrap(full_cds)}\n')
            prot_out.write(f'>{header}_prot\n{wrap(prot)}\n')
            sum_out.write(
                f'{sample_id}\t{query_id_val}\t{scaffold}\t{strand}\t'
                f'{q_start_val}\t{q_end_val}\t'
                f'{contig_start}\t{contig_end}\t'
                f'{len(exons_sorted)}\t{len(full_cds)}\t{len(prot)}\t'
                f'{";".join(coords)}\n'
            )


if __name__ == '__main__':
    main()
