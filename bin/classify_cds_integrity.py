#!/usr/bin/env python3
"""
Classify every Exonerate gene model of a sample as full-length, truncated or
partial, and describe the mutations responsible.

Four causes of an incomplete gene model are distinguished. The first two are
limits of the assembly (the gene may well be intact in the genome), the last
two are evidence of a truncated gene:

  contig_edge         -- the contig ends before the missing query segment (or
        the stop codon) could fit: the gene runs off the contig.

  assembly_gap        -- an assembly gap (N run) next to the missing segment.

  nucleotide_deletion -- part of the query CDS is absent from the contig: a
        terminal query segment is missing although the contig continues with
        unrelated sequence, or a large in-frame internal deletion.

  ORF_mutation        -- the contig carries the query CDS (or the missing part
        is still found nearby) but the reading frame is disrupted: frameshift
        indels, premature stop codons, lost start or stop codons, non-canonical
        splice sites, or a query segment present on the contig but outside the
        ORF (unbridged frameshift / splice-site mutation / insertion).

A model with only contig_edge / assembly_gap causes is 'partial', one with a
nucleotide_deletion or ORF_mutation 'truncated' (or 'extended' when the stop
codon is lost and nothing else). How the encoded ORF ends is reported apart
(orf_end): at the expected stop, cut short by a premature stop, read through to
a later stop, or running off the contig end without a stop codon.

The Exonerate vulgar string is the primary evidence. It is walked codon by
codon, mapping every query codon (nucleotide query) or residue (protein query)
onto the contig so that the frame-corrected codon can be translated and
compared with the reference. Query segments that Exonerate did not align are
checked in two steps: (1) is the ORF simply divergent at that end but still
intact (in-frame start/stop codon at the expected distance)? (2) if not, is
the missing segment present in the flanking contig sequence (rescue alignment
with exonerate est2genome / protein2genome)? A segment that is present means
the ORF is disrupted; a segment that is absent means a nucleotide-level loss.

Usage:
    classify_cds_integrity.py --sample S --query-id Q --query query.fasta \\
        --seq-type nucleotide|protein --exonerate S_exonerate.gff3 \\
        --target S_top_hits.fasta [--assembly scaffolds.fasta] [options]

Outputs (in current directory):
    <sample>_cds_integrity.tsv         one row per gene model (no_hit row if none)
    <sample>_cds_integrity_events.tsv  one row per detected mutation event
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile

from extract_exonerate_cds import CODON, parse_gff3_cds, rc, read_hit_spans


SUMMARY_COLUMNS = [
    'sample_id', 'query_id', 'sequence_id', 'status', 'truncation_cause',
    'contig_edge', 'orf_end', 'deletion_types', 'orf_mutation_types', 'contig',
    'contig_length', 'strand', 'gene_start', 'gene_end',
    'contig_flank_5prime_nt', 'contig_flank_3prime_nt', 'exon_count', 'query_length_aa',
    'aligned_query_start_aa', 'aligned_query_end_aa', 'query_coverage_pct',
    'aa_identity_pct', 'observed_protein_length', 'in_frame_query_aa', 'in_frame_query_pct',
    'start_codon', 'stop_codon', 'n_frameshifts', 'n_premature_stops',
    'n_noncanonical_splice_sites', 'n_inframe_indels', 'first_disruption_aa',
    'mutation_details', 'notes',
]

EVENT_COLUMNS = [
    'sample_id', 'query_id', 'sequence_id', 'event_class', 'event_type',
    'query_aa_start', 'query_aa_end', 'contig', 'contig_start', 'contig_end',
    'ref', 'alt', 'size_nt', 'description',
]

STOPS = {'TAA', 'TAG', 'TGA'}
# Event classes that make a gene model incomplete; the first two are limits of
# the assembly, the last two evidence of a truncated gene
ASSEMBLY_CLASSES = ('contig_edge', 'assembly_gap')
TRUNCATING_CLASSES = ASSEMBLY_CLASSES + ('nucleotide_deletion', 'ORF_mutation')
NEAR_EDGE_NT = 300        # note when the contig ends this close beyond a missing segment
LOCAL_WINDOW = 10         # codons either side of a frameshift for local identity
LOW_LOCAL_IDENTITY = 50.0
CANONICAL_INTRONS = {('GT', 'AG'), ('GC', 'AG'), ('AT', 'AC')}
TELOMERE_MOTIFS = ('TTAGGG', 'CCCTAA', 'TTTAGGG', 'CCCTAAA')


# ---------------------------------------------------------------------------
# Sequence helpers
# ---------------------------------------------------------------------------

def read_fasta_full(path, wanted=None):
    """Return {full_header_id: sequence}; keep only ids in *wanted* if given."""
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(buf).upper()
                name = line[1:].split()[0]
                if wanted is not None and name not in wanted:
                    name = None
                buf = []
            elif name is not None:
                buf.append(line)
    if name is not None:
        seqs[name] = ''.join(buf).upper()
    return seqs


def aa(codon):
    return CODON.get(codon.upper(), 'X')


def fmt_pct(x):
    return f'{x:.1f}' if x is not None else 'NA'


# ---------------------------------------------------------------------------
# Exonerate output parsing
# ---------------------------------------------------------------------------

def parse_vulgar(line):
    parts = line.split()
    fields = parts[1:]
    blk = {
        'q_id': fields[0], 'q_start': int(fields[1]), 'q_end': int(fields[2]),
        'q_strand': fields[3], 't_id': fields[4], 't_start': int(fields[5]),
        't_end': int(fields[6]), 't_strand': fields[7], 'score': int(fields[8]),
    }
    rest = fields[9:]
    blk['ops'] = [(rest[i], int(rest[i + 1]), int(rest[i + 2]))
                  for i in range(0, len(rest) - 2, 3)]
    return blk


def parse_exonerate_blocks(path):
    """Return vulgar blocks in output order (index i <-> i-th GFF gene)."""
    blocks = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('vulgar:'):
                blocks.append(parse_vulgar(line))
    return blocks


def selected_models(gff_path, blocks, sample_id, target_path):
    """Mirror extract_exonerate_cds.py: same models, same FASTA headers."""
    cds_by_parent, _, _ = parse_gff3_cds(gff_path, read_hit_spans(target_path))
    parents = sorted(cds_by_parent)
    models = []
    for idx, pid in enumerate(parents, start=1):
        header = sample_id if len(parents) == 1 else f'{sample_id}_{idx}'
        block_no = int(pid.rsplit('_', 1)[1])
        if 1 <= block_no <= len(blocks):
            models.append((header, blocks[block_no - 1]))
    return models


# ---------------------------------------------------------------------------
# Coordinate context: oriented contig sequence around one gene model
# ---------------------------------------------------------------------------

class Context:
    """Gene-oriented view of the contig carrying a gene model.

    ``seq`` is the contig (or, without assembly, the Exonerate target window)
    oriented 5'->3' along the gene; ``to_contig`` maps a 0-based position in
    ``seq`` to a 1-based contig coordinate.
    """

    def __init__(self, blk, window_seqs, assembly):
        t_id = blk['t_id']
        window = window_seqs[t_id]
        m = re.match(r'^(.+):(\d+)-(\d+)$', t_id)
        contig, ws, we = (m.group(1), int(m.group(2)), int(m.group(3))) if m else (t_id, 1, len(window))

        self.exact = False
        cseq = assembly.get(contig) if assembly else None
        win_rc = False
        if cseq is not None and 0 < ws <= we <= len(cseq):
            fwd = cseq[ws - 1:we]
            if fwd == window:
                self.exact = True
            elif rc(fwd) == window:
                self.exact, win_rc = True, True
        if not self.exact:
            # Without the assembly the window orientation cannot be verified,
            # so coordinates stay relative to the Exonerate target window.
            cseq, ws, we, win_rc = window, 1, len(window), False
            contig = t_id

        self.contig = contig
        self.length = len(cseq)
        # Gene orientation on the contig
        self.minus = (blk['t_strand'] == '-') != win_rc
        self.seq = rc(cseq) if self.minus else cseq

        # First aligned window base (0-based) -> contig -> oriented position
        w0 = blk['t_start'] if blk['t_strand'] == '+' else blk['t_start'] - 1
        c0 = (we - 1 - w0) if win_rc else (ws - 1 + w0)
        self.start = (self.length - 1 - c0) if self.minus else c0

    @property
    def strand(self):
        return '-' if self.minus else '+'

    def to_contig(self, p):
        return self.length - p if self.minus else p + 1

    def span(self, p0, p1):
        """Contig coordinates (start, end) of oriented half-open [p0, p1)."""
        n = len(self.seq)
        p0 = min(max(p0, 0), n - 1)
        p1 = min(max(p1, p0 + 1), n)
        a, b = self.to_contig(p0), self.to_contig(p1 - 1)
        return (min(a, b), max(a, b))


# ---------------------------------------------------------------------------
# Vulgar walk
# ---------------------------------------------------------------------------

class Walk:
    """Walk vulgar operations along the gene-oriented contig.

    codons  : {query_aa_index: [ctx_pos, ctx_pos, ctx_pos]} frame-corrected
              (codon bases placed by their aligned query position)
    raw     : ctx positions of every coding base in transcript order (what the
              gene model would actually translate, frameshift bases included)
    indels  : merged G/F runs -> dicts (q_aa, net_nt, q_nt, t_nt, ctx_pos)
    introns : (ctx_start, ctx_end, donor, acceptor, q_aa)
    """

    def __init__(self, blk, ctx, nucleotide_query):
        self.codons = {}
        self.raw = []
        self.indels = []
        self.introns = []
        self.exons = []
        self.unaligned = []
        seq = ctx.seq
        tpos = ctx.start
        nt = nucleotide_query
        q = blk['q_start']              # nt (nucleotide query) or aa (protein)
        split_buf, split_aa = [], None
        pending = None                  # current merged indel run
        intron_start = donor = None
        exon_start = tpos

        def flush_indel():
            nonlocal pending
            if pending is not None:
                self.indels.append(pending)
                pending = None

        for op, ql, tl in blk['ops']:
            q_nt = ql if nt else 3 * ql
            if op in ('G', 'F'):
                if pending is None:
                    q_aa = (q // 3) if nt else q
                    pending = {'q_aa': q_aa, 'q_nt': 0, 't_nt': 0, 'ctx_pos': tpos,
                               'ops': ''}
                pending['q_nt'] += q_nt
                pending['t_nt'] += tl
                pending['ops'] += op
                self.raw.extend(range(tpos, tpos + tl))
                q += ql
                tpos += tl
                continue
            flush_indel()

            if op in ('M', 'C'):
                if nt:
                    for i in range(min(ql, tl)):
                        self.codons.setdefault((q + i) // 3, [None] * 3)[(q + i) % 3] = tpos + i
                else:
                    for i in range(min(ql, tl // 3)):
                        self.codons[q + i] = [tpos + 3 * i, tpos + 3 * i + 1, tpos + 3 * i + 2]
                self.raw.extend(range(tpos, tpos + tl))
                q += ql
                tpos += tl
            elif op == 'S':
                if nt:
                    for i in range(min(ql, tl)):
                        self.codons.setdefault((q + i) // 3, [None] * 3)[(q + i) % 3] = tpos + i
                else:
                    if not split_buf:
                        split_aa = q
                    split_buf.extend(range(tpos, tpos + tl))
                    if len(split_buf) >= 3:
                        self.codons[split_aa] = split_buf[:3]
                        split_buf = []
                self.raw.extend(range(tpos, tpos + tl))
                q += ql
                tpos += tl
            elif op == '5':
                self.exons.append((exon_start, tpos))
                intron_start = tpos
                donor = seq[tpos:tpos + 2]
                tpos += tl
            elif op == 'I':
                tpos += tl
            elif op == '3':
                acceptor = seq[tpos:tpos + 2]
                tpos += tl
                q_aa = (q // 3) if nt else q
                self.introns.append((intron_start, tpos, donor, acceptor, q_aa))
                exon_start = tpos
            elif op == 'N':
                self.unaligned.append(((q // 3) if nt else q, ql, tl, tpos))
                self.raw.extend(range(tpos, tpos + tl))
                q += ql
                tpos += tl
            else:
                # Unknown op: advance conservatively
                q += ql
                tpos += tl
        flush_indel()
        self.exons.append((exon_start, tpos))
        self.end = tpos                  # oriented position after last base
        # Drop incomplete codons (frameshift / gap boundaries)
        self.codons = {k: v for k, v in self.codons.items() if None not in v}


# ---------------------------------------------------------------------------
# Rescue alignment of an unaligned query segment against flanking sequence
# ---------------------------------------------------------------------------

MIN_SEGMENT_AA = 4       # shorter unaligned termini are not re-aligned
SHORT_SEGMENT_AA = 30    # up to this length, re-align locally in Python first


def local_align(a, b, match, mismatch, gap):
    """Smith-Waterman (linear gaps). Return (score, a0, a1, b0, b1, n_ident, n_cols)."""
    n, m = len(a), len(b)
    H = [[0] * (m + 1) for _ in range(n + 1)]
    best = (0, 0, 0)
    for i in range(1, n + 1):
        ai, Hi, Hp = a[i - 1], H[i], H[i - 1]
        for j in range(1, m + 1):
            v = Hp[j - 1] + (match if ai == b[j - 1] else mismatch)
            u = Hp[j] + gap
            if u > v:
                v = u
            u = Hi[j - 1] + gap
            if u > v:
                v = u
            if v < 0:
                v = 0
            Hi[j] = v
            if v > best[0]:
                best = (v, i, j)
    score, i, j = best
    a1, b1, ident, cols = i, j, 0, 0
    while i > 0 and j > 0 and H[i][j] > 0:
        cols += 1
        if H[i][j] == H[i - 1][j - 1] + (match if a[i - 1] == b[j - 1] else mismatch):
            ident += a[i - 1] == b[j - 1]
            i, j = i - 1, j - 1
        elif H[i][j] == H[i - 1][j] + gap:
            i -= 1
        else:
            j -= 1
    return score, i, a1, j, b1, ident, cols


def local_rescue(segment, near, nucleotide_query, min_frac=0.5, min_ident=0.7):
    """Re-align a short missing segment to sequence adjacent to the alignment.

    Nucleotide segments are aligned directly; protein segments against the
    three forward-frame translations. Returns a hit dict like rescue_segment.
    """
    if not segment or not near:
        return None
    cands = []
    if nucleotide_query:
        r = local_align(segment, near, 2, -3, -5)
        cands.append((r, 0, 1))
    else:
        for frame in range(3):
            tr = ''.join(aa(near[i:i + 3]) for i in range(frame, len(near) - 2, 3))
            cands.append((local_align(segment, tr, 3, -1, -4), frame, 3))
    best = None
    for (score, a0, a1, b0, b1, ident, cols), frame, unit in cands:
        if cols == 0:
            continue
        frac = (a1 - a0) / len(segment)
        min_len = 8 if nucleotide_query else 3
        if frac < min_frac or (a1 - a0) < min_len or ident / cols < min_ident:
            continue
        if best is None or score > best['score']:
            best = {'q0': a0, 'q1': a1, 't0': frame + unit * b0, 't1': frame + unit * b1,
                    'score': score, 'frac': frac}
    return best


def rescue_segment(segment, flank, nucleotide_query, min_frac):
    """Align *segment* to *flank* with Exonerate; return best forward hit or None.

    Hit dict: q0, q1 (segment units), t0, t1 (flank nt), score, frac.
    """
    exe = shutil.which('exonerate')
    if exe is None or not segment or len(flank) < 30:
        return None
    model = 'est2genome' if nucleotide_query else 'protein2genome'
    with tempfile.TemporaryDirectory() as tmp:
        qf, tf = os.path.join(tmp, 'q.fa'), os.path.join(tmp, 't.fa')
        with open(qf, 'w') as fh:
            fh.write(f'>segment\n{segment}\n')
        with open(tf, 'w') as fh:
            fh.write(f'>flank\n{flank}\n')
        cmd = [exe, '--model', model, '--query', qf, '--target', tf,
               '--refine', 'full', '--score', '20', '--percent', '25',
               '--showalignment', 'no', '--showvulgar', 'yes', '--bestn', '3']
        try:
            out = subprocess.run(cmd, capture_output=True, text=True, timeout=300).stdout
        except (subprocess.SubprocessError, OSError):
            return None
    best = None
    for line in out.splitlines():
        if not line.startswith('vulgar:'):
            continue
        b = parse_vulgar(line)
        if b['t_strand'] != '+' or b['q_strand'] == '-':
            continue
        frac = (b['q_end'] - b['q_start']) / len(segment)
        if frac < min_frac:
            continue
        if best is None or b['score'] > best['score']:
            best = {'q0': b['q_start'], 'q1': b['q_end'], 't0': b['t_start'],
                    't1': b['t_end'], 'score': b['score'], 'frac': frac}
    return best


# ---------------------------------------------------------------------------
# Classification of one gene model
# ---------------------------------------------------------------------------

class Model:

    def __init__(self, args, header, blk, ctx, query_seq, nucleotide_query):
        self.args = args
        self.header = header
        self.blk = blk
        self.ctx = ctx
        self.nt = nucleotide_query
        self.events = []
        self.notes = []
        self.edge_ends = set()          # '5prime' / '3prime' ends running off the contig
        self.query = query_seq

        if nucleotide_query:
            n_codons = len(query_seq) // 3
            self.ref_aa = ''.join(aa(query_seq[3 * i:3 * i + 3]) for i in range(n_codons))
            self.ref_has_stop = self.ref_aa.endswith('*')
            self.expected = n_codons - 1 if self.ref_has_stop else n_codons
            if len(query_seq) % 3:
                self.notes.append('query_length_not_multiple_of_3')
        else:
            self.ref_aa = query_seq.rstrip('*')
            self.ref_has_stop = False
            self.expected = len(self.ref_aa)
        self.walk = Walk(blk, ctx, nucleotide_query)

    # -- helpers -----------------------------------------------------------

    def codon_seq(self, k):
        return ''.join(self.ctx.seq[p] for p in self.walk.codons[k])

    def add(self, cls, etype, q0, q1, p0, p1, ref, alt, size, desc):
        if p0 is not None:
            c0, c1 = self.ctx.span(p0, p1)
        else:
            c0 = c1 = 'NA'
        self.events.append({
            'event_class': cls, 'event_type': etype,
            'query_aa_start': q0 + 1 if q0 is not None else 'NA',
            'query_aa_end': q1 if q1 is not None else 'NA',
            'contig': self.ctx.contig, 'contig_start': c0, 'contig_end': c1,
            'ref': ref or '.', 'alt': alt or '.', 'size_nt': size if size is not None else 'NA',
            'description': desc,
        })

    def local_identity(self, k):
        """aa identity (%) of aligned codons within LOCAL_WINDOW of residue k."""
        ks = [j for j in self.walk.codons if abs(j - k) <= LOCAL_WINDOW and j < self.expected]
        if not ks:
            return None
        return 100.0 * sum(aa(self.codon_seq(j)) == self.ref_aa[j] for j in ks) / len(ks)

    def ref_label(self, k):
        return f'{self.ref_aa[k]}{k + 1}' if 0 <= k < len(self.ref_aa) else f'{k + 1}'

    def aa_at_ctx(self, p):
        """Query residue index whose codon contains (or precedes) ctx position p."""
        best = None
        for k, pos in self.walk.codons.items():
            if pos[0] <= p and (best is None or k > best):
                best = k
        return best

    def read_frame(self, p):
        """Translate ctx forward from p in frame; return (n_sense, stop_pos or None)."""
        seq = self.ctx.seq
        n = 0
        while p + 3 <= len(seq):
            if seq[p:p + 3] in STOPS:
                return n, p
            n += 1
            p += 3
        return n, None

    # -- analysis ------------------------------------------------------------

    def analyse(self):
        a = self.args
        w = self.walk
        seq = self.ctx.seq
        codons = w.codons
        aligned = sorted(k for k in codons if k < self.expected)
        self.first_aa = aligned[0] if aligned else None
        self.last_aa = aligned[-1] if aligned else None

        # Aligned query range in residues (from vulgar)
        if self.nt:
            qa0 = -(-self.blk['q_start'] // 3)
            qa1 = min(self.blk['q_end'] // 3, self.expected)
        else:
            qa0, qa1 = self.blk['q_start'], min(self.blk['q_end'], self.expected)
        self.qa0, self.qa1 = qa0, qa1

        # --- codon-by-codon comparison ------------------------------------
        same = 0
        n_amb = 0
        for k in aligned:
            cod = self.codon_seq(k)
            alt = aa(cod)
            ref = self.ref_aa[k]
            if 'N' in cod:
                n_amb += 1
            if alt == ref:
                same += 1
            if alt == '*' and ref != '*':
                ref_cod = self.query[3 * k:3 * k + 3] if self.nt else None
                p = codons[k][0]
                self.add('ORF_mutation', 'premature_stop', k, k + 1, p, codons[k][2] + 1,
                         f'{ref_cod}({ref})' if ref_cod else ref, f'{cod}(*)', 3,
                         f'p.{self.ref_label(k)}* nonsense codon {cod}')
        self.identity = 100.0 * same / len(aligned) if aligned else None
        if n_amb:
            self.notes.append(f'ambiguous_bases_in_{n_amb}_codons')

        # --- indels (G/F runs) ---------------------------------------------
        for ind in w.indels:
            net = ind['t_nt'] - ind['q_nt']
            k = ind['q_aa']
            p0 = ind['ctx_pos']
            p1 = p0 + max(ind['t_nt'], 1)
            if net % 3:
                local = self.local_identity(k)
                desc = (f'p.{self.ref_label(k)}fs frameshift {net:+d} nt '
                        f'({"insertion" if net > 0 else "deletion"} in contig)')
                if local is not None and local < LOW_LOCAL_IDENTITY:
                    desc += (f' [low confidence: {local:.0f}% aa identity within '
                             f'{LOCAL_WINDOW} codons; alignment may be forced through divergent sequence]')
                    self.notes.append(f'low_confidence_frameshift_at_aa_{k + 1}')
                self.add('ORF_mutation', 'frameshift', k, k + 1, p0, p1,
                         '.', '.', net, desc)
            elif net < 0:
                n_aa = -net // 3
                big = n_aa >= a.min_deletion_aa
                self.add('nucleotide_deletion' if big else 'in_frame_indel',
                         'internal_deletion' if big else 'inframe_deletion',
                         k, k + n_aa, p0, p0 + 1, self.ref_aa[k:k + n_aa], '-', net,
                         f'in-frame deletion of {n_aa} aa (query aa {k + 1}-{k + n_aa})')
            elif net > 0:
                self.add('in_frame_indel', 'inframe_insertion', k, k, p0, p0 + net,
                         '-', seq[p0:p0 + net], net,
                         f'in-frame insertion of {net // 3} codon(s) after query aa {k}')

        for (qk, ql, tl, p) in w.unaligned:
            self.add('in_frame_indel', 'unaligned_block', qk, qk + (ql if not self.nt else ql // 3),
                     p, p + max(tl, 1), '.', '.', tl - (ql if self.nt else 3 * ql),
                     'non-equivalenced alignment block')

        # --- splice sites ------------------------------------------------------
        for (i0, i1, donor, acceptor, qk) in w.introns:
            if (donor, acceptor) not in CANONICAL_INTRONS:
                which = []
                if donor not in ('GT', 'GC', 'AT'):
                    which.append(f"donor GT>{donor}")
                if acceptor not in ('AG', 'AC'):
                    which.append(f"acceptor AG>{acceptor}")
                if not which:
                    which.append(f'non-canonical pair {donor}-{acceptor}')
                if len(which) == 2:
                    which.append('both ends non-canonical: possibly an insertion '
                                 '(e.g. transposable element) modelled as an intron')
                self.add('ORF_mutation', 'noncanonical_splice_site', qk, qk + 1, i0, i1,
                         'GT..AG', f'{donor}..{acceptor}', None,
                         f'intron between query aa {qk} and {qk + 1} ({i1 - i0} nt): ' + ', '.join(which))

        # --- start codon -------------------------------------------------------
        self.start_status = 'NA' if qa0 == 0 else 'not_aligned'
        if qa0 == 0 and 0 in codons:
            cod = self.codon_seq(0)
            if self.ref_aa[:1] == 'M':
                if cod == 'ATG':
                    self.start_status = 'ATG'
                else:
                    self.start_status = f'lost({cod})'
                    nxt = next((k for k in sorted(codons) if 0 < k < self.expected
                                and self.codon_seq(k) == 'ATG'), None)
                    alt_start = (f'; next in-frame ATG at query aa {nxt + 1}' if nxt is not None
                                 else '; no downstream in-frame ATG')
                    self.add('ORF_mutation', 'start_lost', 0, 1, codons[0][0], codons[0][2] + 1,
                             'ATG(M)', f'{cod}({aa(cod)})', 3, f'start codon ATG>{cod}{alt_start}')
            else:
                self.start_status = f'{cod}(query_not_Met)'
        # --- stop codon / 3' end ------------------------------------------------
        self.stop_status = 'not_aligned'
        self.ext_downstream = None

        # --- terminal truncations ---------------------------------------------
        self.check_five_prime()
        self.check_three_prime()

        # --- observed protein (actual spliced model incl. frameshift bases) ----
        self.observed_protein()

    # ..................................................................

    def short_terminus(self, m):
        """Unaligned termini up to this size may be accepted as divergent-but-intact."""
        return m <= max(10, self.args.max_unaligned_terminus * self.expected)

    def tol(self, m):
        return max(3, int(round(self.args.terminal_tolerance * m)))

    def flank_notes(self, flank):
        notes = []
        if re.search(r'N{10,}', flank):
            notes.append('assembly_gap')
        if any(flank.count(motif) >= 4 for motif in TELOMERE_MOTIFS):
            notes.append('telomeric_repeats')
        return notes

    def check_five_prime(self):
        # m = residues missing before the first complete aligned codon
        codons = self.walk.codons
        m = self.first_aa
        if not m:
            return
        seq = self.ctx.seq
        f = codons[m][0]
        # In-frame upstream scan for ATG before the first upstream stop
        atgs, k, p = [], 1, f - 3
        while p >= 0:
            cod = seq[p:p + 3]
            if cod in STOPS:
                break
            if cod == 'ATG':
                atgs.append(k)
            k += 1
            p -= 3
        near = [x for x in atgs if abs(x - m) <= self.tol(m)]
        if near and self.short_terminus(m):
            best = min(near, key=lambda x: abs(x - m))
            self.start_status = 'ATG(unaligned_5prime_end)'
            self.notes.append(f'5prime_end_unaligned_but_ORF_intact(ATG {best} codons upstream, expected {m})')
            return
        self.classify_missing('5prime', 0, m, f)

    def check_three_prime(self):
        codons = self.walk.codons
        if self.last_aa is None:
            return
        last = self.last_aa
        l = codons[last][2] + 1          # first base after the last aligned codon
        m = self.expected - 1 - last     # sense codons still expected
        if m <= 0:
            n, stop = self.read_frame(l)
            if l + 3 > len(self.ctx.seq):
                # The coding sequence reaches the contig end: the stop codon
                # position itself is not in the assembly
                self.stop_status = 'beyond_contig_end'
                self.add_edge('3prime', self.expected, self.expected + 1, l, len(self.ctx.seq),
                              '(stop)', f'C-terminal stop codon missing: {self.edge_reason(len(self.ctx.seq) - l)}')
            elif stop == l:
                self.stop_status = self.ctx.seq[l:l + 3]
            else:
                self.stop_status = 'lost'
                ext = f'+{n} aa until {self.ctx.seq[stop:stop + 3]}' if stop is not None \
                    else f'no stop codon within the remaining {len(self.ctx.seq) - l} nt'
                ref_cod = self.query[-3:] if self.ref_has_stop else '(stop)'
                self.add('ORF_mutation', 'stop_lost', self.expected, self.expected + 1,
                         l, l + 3, ref_cod, self.ctx.seq[l:l + 3], 3,
                         f'stop codon lost, read-through {ext}')
            return
        n, stop = self.read_frame(l)
        if stop is not None and abs(n - m) <= self.tol(m) and self.short_terminus(m):
            self.stop_status = f'{self.ctx.seq[stop:stop + 3]}(unaligned_3prime_end)'
            self.notes.append(f'3prime_end_unaligned_but_ORF_intact(stop {n} codons downstream, expected {m})')
            return
        self.ext_downstream = (n, stop)
        self.classify_missing('3prime', last + 1, self.expected, l)

    def classify_missing(self, end, k0, k1, anchor):
        """Query residues [k0, k1) are not aligned at the *end* terminus.

        anchor is the oriented contig position adjacent to the alignment
        (first aligned base for 5prime, first base after for 3prime).
        """
        a = self.args
        seq = self.ctx.seq
        m = k1 - k0
        need = 3 * m
        span = need + a.rescue_flank
        if end == '5prime':
            flank_len = anchor
            fl0 = max(0, anchor - span)
            flank = seq[fl0:anchor]
        else:
            flank_len = len(seq) - anchor
            fl0 = anchor
            flank = seq[anchor:anchor + span]
        notes = self.flank_notes(flank[:need + 500] if end == '3prime' else flank[-(need + 500):])
        label = 'N-terminal' if end == '5prime' else 'C-terminal'
        where = f'query aa {k0 + 1}-{k1} ({m} aa, {need} nt)'
        # The 3' segment also needs room for the stop codon
        edge_need = need + 3 if end == '3prime' else need
        p0, p1 = (max(0, anchor - 1), anchor) if end == '5prime' else (anchor, anchor + 1)

        if m < MIN_SEGMENT_AA:
            # Too short to re-align. If the contig ends before the missing
            # codons (and stop codon) could fit, it is a contig edge, not a
            # lost start / stop codon.
            if flank_len < edge_need:
                ref = self.ref_aa[k0:k1] if end == '5prime' else self.ref_aa[k0:k1] + '(stop)'
                self.add_edge(end, k0, k1, p0, p1, ref, f'{label} {where} missing: {self.edge_reason(flank_len)}')
                return
            # Too short to re-align: only the terminal codon is unaccounted for
            if end == '5prime':
                self.add('ORF_mutation', 'start_lost', k0, k1, max(0, anchor - 3), anchor,
                         'ATG(M)', '.', None,
                         f'{label} {where} unaligned and no in-frame ATG near the expected start')
            else:
                n, sp = self.ext_downstream
                ext = f'ORF continues {n} codon(s) until {seq[sp:sp + 3]}' if sp is not None \
                    else 'no in-frame stop codon before the contig end'
                self.add('ORF_mutation', 'stop_lost', k0, k1, anchor, anchor + 3,
                         '(stop)', '.', None,
                         f'{label} {where} unaligned and no stop codon at the expected position: {ext}')
            return

        hit = None
        if not a.no_rescue:
            if self.nt:
                # 3' segment keeps the reference stop codon as extra anchor
                q1 = len(self.query) if (end == '3prime' and self.ref_has_stop) else 3 * k1
                segment = self.query[3 * k0:q1]
            else:
                segment = self.ref_aa[k0:k1]
            if m <= SHORT_SEGMENT_AA:
                # Small indels only: search just beyond the expected position
                near = seq[max(0, anchor - need - 30):anchor] if end == '5prime' \
                    else seq[anchor:anchor + need + 30]
                off = max(0, anchor - need - 30) if end == '5prime' else anchor
                hit = local_rescue(segment, near, self.nt)
                if hit is not None:
                    hit['t0'] += off - fl0
                    hit['t1'] += off - fl0
            if hit is None and m >= a.min_rescue_aa:
                hit = rescue_segment(segment, flank, self.nt, a.rescue_min_fraction)

        if hit is not None:
            unit = 1 if self.nt else 3
            rq0, rq1 = 3 * k0 + hit['q0'] * unit, 3 * k0 + hit['q1'] * unit
            rt0, rt1 = fl0 + hit['t0'], fl0 + hit['t1']
            if end == '3prime':
                d = (rt0 - anchor) - (rq0 - 3 * k0)
            else:
                d = (anchor - rt1) - (3 * k1 - rq1)
            if abs(d) <= a.max_indel_nt and d % 3:
                mech = f'frameshift (net {d:+d} nt) not bridged by the alignment'
            elif abs(d) <= a.max_indel_nt:
                mech = 'in-frame disruption (premature stop / small indel)'
            elif d > 0:
                mech = (f'{d} nt of extra sequence between aligned parts '
                        '(retained intron / splice-site mutation or insertion)')
            else:
                mech = 'segment re-located upstream/overlapping (rearrangement or duplication)'
            stop_txt = ''
            if end == '3prime' and self.ext_downstream and self.ext_downstream[1] is not None:
                n, sp = self.ext_downstream
                stop_txt = f'; ORF ends {n} codon(s) after the last aligned residue ({seq[sp:sp + 3]})'
            self.add('ORF_mutation', 'segment_outside_ORF', k0, k1, rt0, rt1,
                     '.', '.', d,
                     f'{label} {where} present on contig but outside the ORF '
                     f'({100 * hit["frac"]:.0f}% of segment re-aligned): {mech}{stop_txt}')
            return

        ref = self.ref_aa[k0:k1] if m <= 30 else f'{m}aa'
        if flank_len < need:
            self.add_edge(end, k0, k1, p0, p1, ref, f'{label} {where} missing: {self.edge_reason(flank_len)}')
            return
        if 'assembly_gap' in notes:
            ecls, etype, why = 'assembly_gap', 'assembly_gap', 'assembly gap (N run) in the flanking sequence'
        else:
            ecls, etype, why = ('nucleotide_deletion', f'{end}_deletion',
                                'segment absent; contig continues with unrelated sequence')
            if flank_len - need < NEAR_EDGE_NT:
                why += f' (contig ends {flank_len} nt beyond the alignment)'
                self.notes.append(f'{end}_contig_end_{flank_len}nt_beyond_alignment')
        if 'telomeric_repeats' in notes:
            why += ' (telomeric repeats in flank)'
        if end == '3prime' and self.ext_downstream and self.ext_downstream[1] is not None:
            n, sp = self.ext_downstream
            why += f'; ORF continues {n} codon(s) into flank then {seq[sp:sp + 3]}'
        self.add(ecls, etype, k0, k1, p0, p1, ref, '-', -need,
                 f'{label} {where} missing: {why}')

    def edge_reason(self, flank_len):
        if self.ctx.exact:
            return f'contig ends {flank_len} nt beyond the alignment'
        return f'extracted window ends {flank_len} nt beyond the alignment (assembly not supplied)'

    def add_edge(self, end, k0, k1, p0, p1, ref, desc):
        """The gene runs off the contig (or the extracted window) at *end*."""
        self.edge_ends.add(end)
        self.add('contig_edge', 'contig_end' if self.ctx.exact else 'window_end',
                 k0, k1, p0, p1, ref, '-', None, desc)

    def observed_protein(self):
        """Translate the actual model (raw bases, then flank) from its first codon."""
        seq = self.ctx.seq
        raw = self.walk.raw
        codons = self.walk.codons
        self.obs_len = None
        self.obs_stop_aa = None
        self.orf_end = 'NA'
        self.in_frame = 0
        if self.first_aa is None:
            return
        start = codons[self.first_aa][0]
        try:
            i = raw.index(start)
        except ValueError:
            return
        bases = raw[i:]
        n = 0
        j = 0
        translated = set()      # codons (as position triples) actually read in frame
        while j + 3 <= len(bases):
            trip = tuple(bases[j:j + 3])
            cod = ''.join(seq[p] for p in trip)
            if cod in STOPS:
                self.obs_len = n
                self.obs_stop_aa = self.aa_at_ctx(bases[j])
                stop_q = self.obs_stop_aa + 1 if self.obs_stop_aa is not None else None
                break
            translated.add(trip)
            n += 1
            j += 3
        else:
            # No stop within the model: continue in frame into the flank
            p = bases[j] if j < len(bases) else bases[-1] + 1
            ext, stop = self.read_frame(p)
            self.obs_len = n + ext
            last = self.aa_at_ctx(p - 1)
            stop_q = last + 1 + ext if (stop is not None and last is not None) else None
            if stop is None:
                self.orf_end = 'contig_end'
                self.edge_ends.add('3prime')
        # Where the ORF stops, relative to the expected stop codon (query
        # residue index); the tolerance follows the divergent-terminus rule
        if stop_q is not None:
            tol = self.tol(max(0, self.expected - 1 - self.last_aa))
            if stop_q < self.expected - tol:
                self.orf_end = 'premature_stop'
            elif stop_q > self.expected + tol:
                self.orf_end = 'late_stop'
            else:
                self.orf_end = 'expected_stop'
        # Reference residues encoded in the correct frame before the stop
        self.in_frame = sum(1 for k, pos in codons.items()
                            if k < self.expected and tuple(pos) in translated)

    # -- summary -------------------------------------------------------------

    def summary(self, sample_id, query_id):
        ev = self.events
        del_types = sorted({e['event_type'] for e in ev if e['event_class'] == 'nucleotide_deletion'})
        orf_types = sorted({e['event_type'] for e in ev if e['event_class'] == 'ORF_mutation'})

        # Frameshift / splice defects truncate the product even when the
        # observed length is long; stop_lost alone means an extended protein.
        # A model incomplete only because of the assembly (contig end, gap) is
        # partial: there is no evidence that the gene itself is truncated.
        causes = [c for c in TRUNCATING_CLASSES if any(e['event_class'] == c for e in ev)]
        # A non-canonical splice site says where the model put an exon boundary, not that the
        # product is broken: a boundary placed a base or two off reads as a non-canonical acceptor
        # while the protein is untouched. On its own, with the reading frame still reaching the
        # expected stop, it is reported but does not truncate. (Without this, one shared boundary
        # call made every isolate of a 4,068-genome panel 'truncated' for the same effector.)
        hard_orf = sorted(t for t in orf_types if t != 'noncanonical_splice_site')
        splice_only = bool(orf_types) and not hard_orf
        if not causes:
            status = 'full_length'
        elif not del_types and not orf_types:
            status = 'partial'
        elif not del_types and splice_only and self.orf_end == 'expected_stop':
            status = 'full_length'
            causes = [c for c in causes if c != 'ORF_mutation']
            self.notes.append('noncanonical_splice_site_only_ORF_intact')
        elif not del_types and set(hard_orf) == {'stop_lost'}:
            status = 'extended'
        else:
            status = 'truncated'

        # Observed protein: only meaningful from the true start codon
        obs = self.obs_len
        full_start = self.qa0 == 0 or self.start_status.startswith('ATG(unaligned')
        obs_txt = str(obs) if obs is not None else 'NA'
        if obs is not None and (not full_start or self.orf_end == 'contig_end'):
            obs_txt = f'>={obs}'

        if status == 'truncated' and self.obs_stop_aa is not None and obs is not None \
                and full_start and obs < self.expected:
            self.notes.append(f'protein_ends_at_query_aa_{self.obs_stop_aa + 1}')

        disrupt = [e['query_aa_start'] for e in ev
                   if e['event_class'] in TRUNCATING_CLASSES and e['query_aa_start'] != 'NA']
        first_disrupt = min(disrupt) if disrupt else 'NA'

        detail = '; '.join(e['description'] for e in ev if e['event_class'] in TRUNCATING_CLASSES)
        minor = [e for e in ev if e['event_class'] == 'in_frame_indel']
        if minor:
            detail = '; '.join(filter(None, [detail, f'{len(minor)} minor in-frame indel(s): ' +
                                             ', '.join(e['description'] for e in minor)]))

        raw = self.walk.raw
        spans = [self.ctx.to_contig(p) for p in (raw[0], raw[-1])] if raw else ['NA', 'NA']
        # Contig sequence beyond the coding part of the model, in gene orientation
        flanks = [raw[0], len(self.ctx.seq) - 1 - raw[-1]] if raw and self.ctx.exact else ['NA', 'NA']
        edge = ','.join(e for e in ('5prime', '3prime') if e in self.edge_ends) or 'none'
        coverage = 100.0 * max(0, self.qa1 - self.qa0) / self.expected if self.expected else None

        return {
            'sample_id': sample_id, 'query_id': query_id, 'sequence_id': self.header,
            'status': status, 'truncation_cause': '+'.join(causes) or 'none',
            'contig_edge': edge, 'orf_end': self.orf_end,
            'deletion_types': ','.join(del_types) or 'none',
            'orf_mutation_types': ','.join(orf_types) or 'none',
            'contig': self.ctx.contig,
            'contig_length': self.ctx.length if self.ctx.exact else 'NA',
            'strand': self.ctx.strand if self.ctx.exact else self.blk['t_strand'],
            'gene_start': spans[0], 'gene_end': spans[1],
            'contig_flank_5prime_nt': flanks[0], 'contig_flank_3prime_nt': flanks[1],
            'exon_count': len(self.walk.introns) + 1,
            'query_length_aa': self.expected,
            'aligned_query_start_aa': self.qa0 + 1, 'aligned_query_end_aa': self.qa1,
            'query_coverage_pct': fmt_pct(coverage),
            'aa_identity_pct': fmt_pct(self.identity),
            'observed_protein_length': obs_txt,
            'in_frame_query_aa': self.in_frame,
            'in_frame_query_pct': fmt_pct(100.0 * self.in_frame / self.expected if self.expected else None),
            'start_codon': self.start_status,
            'stop_codon': self.stop_status,
            'n_frameshifts': sum(e['event_type'] == 'frameshift' for e in ev),
            'n_premature_stops': sum(e['event_type'] == 'premature_stop' for e in ev),
            'n_noncanonical_splice_sites': sum(e['event_type'] == 'noncanonical_splice_site' for e in ev),
            'n_inframe_indels': len(minor) + sum(e['event_type'] == 'internal_deletion' for e in ev),
            'first_disruption_aa': first_disrupt,
            'mutation_details': detail or 'none',
            'notes': ','.join(self.notes) or '.',
        }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--sample', required=True)
    ap.add_argument('--query-id', required=True)
    ap.add_argument('--query', required=True, help='query FASTA (single sequence)')
    ap.add_argument('--seq-type', required=True, choices=['nucleotide', 'protein'])
    ap.add_argument('--exonerate', required=True, help='raw Exonerate output (vulgar + target GFF)')
    ap.add_argument('--target', required=True, help='Exonerate target FASTA (BLAST hit windows)')
    ap.add_argument('--assembly', help='full assembly FASTA (exact contig ends and coordinates)')
    ap.add_argument('--terminal-tolerance', type=float, default=0.25,
                    help='fractional tolerance on the distance to an in-frame start/stop codon '
                         'when an unaligned terminus is accepted as divergent-but-intact (min 3 codons)')
    ap.add_argument('--max-unaligned-terminus', type=float, default=0.15,
                    help='largest unaligned terminus (fraction of query length, min 10 residues) '
                         'that can be accepted as divergent-but-intact when the ORF start/stop '
                         'codon lies at the expected distance')
    ap.add_argument('--min-deletion-aa', type=int, default=10,
                    help='in-frame internal deletions of at least this many residues count as '
                         'nucleotide-level deletions (smaller ones are reported as minor indels)')
    ap.add_argument('--min-rescue-aa', type=int, default=5,
                    help='minimum unaligned terminal segment (residues) to attempt a rescue alignment')
    ap.add_argument('--rescue-flank', type=int, default=3000,
                    help='extra flanking nt searched beyond the missing segment length')
    ap.add_argument('--rescue-min-fraction', type=float, default=0.3,
                    help='fraction of the missing segment a rescue alignment must cover')
    ap.add_argument('--max-indel-nt', type=int, default=50,
                    help='offsets up to this size between aligned and rescued parts are read as small indels')
    ap.add_argument('--no-rescue', action='store_true')
    ap.add_argument('--out-prefix', help='output file prefix (default: sample id)')
    args = ap.parse_args()

    qseqs = read_fasta_full(args.query)
    query_seq = next(iter(qseqs.values()), '').replace('U', 'T')
    nucleotide = args.seq_type == 'nucleotide'

    prefix = args.out_prefix or args.sample
    sum_path = f'{prefix}_cds_integrity.tsv'
    ev_path = f'{prefix}_cds_integrity_events.tsv'
    rows, events = [], []

    blocks = parse_exonerate_blocks(args.exonerate) if os.path.exists(args.exonerate) else []
    models = selected_models(args.exonerate, blocks, args.sample, args.target) if blocks else []

    if models:
        window_seqs = read_fasta_full(args.target)
        assembly = {}
        if args.assembly and os.path.exists(args.assembly):
            wanted = {re.sub(r':\d+-\d+$', '', b['t_id']) for _, b in models}
            assembly = read_fasta_full(args.assembly, wanted)
        for header, blk in models:
            if blk['t_id'] not in window_seqs:
                sys.stderr.write(f'WARNING: target {blk["t_id"]} not in {args.target}\n')
                continue
            ctx = Context(blk, window_seqs, assembly)
            model = Model(args, header, blk, ctx, query_seq, nucleotide)
            model.analyse()
            rows.append(model.summary(args.sample, args.query_id))
            for e in model.events:
                events.append(dict(e, sample_id=args.sample, query_id=args.query_id,
                                   sequence_id=header))

    if not rows:
        rows.append(dict({c: 'NA' for c in SUMMARY_COLUMNS},
                         sample_id=args.sample, query_id=args.query_id, sequence_id='.',
                         status='no_hit', truncation_cause='NA', notes='no Exonerate alignment'))

    with open(sum_path, 'w') as out:
        out.write('\t'.join(SUMMARY_COLUMNS) + '\n')
        for r in rows:
            out.write('\t'.join(str(r.get(c, 'NA')) for c in SUMMARY_COLUMNS) + '\n')
    with open(ev_path, 'w') as out:
        out.write('\t'.join(EVENT_COLUMNS) + '\n')
        for e in events:
            out.write('\t'.join(str(e.get(c, 'NA')) for c in EVENT_COLUMNS) + '\n')


if __name__ == '__main__':
    main()
