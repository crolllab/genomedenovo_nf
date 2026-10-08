#!/usr/bin/env python3
"""
List the transposable elements (TEs) on the scaffold that carries the primary
hit of a query in each sample, with their position and distance relative to
the hit in the hit's own 5'->3' orientation.

Usage:
    summarize_hit_tes.py <query_id> <hits.tsv> --rm-out <*_repeatmasker.out> ...

Input:
    hits.tsv   *_all_blast_hits_summary.tsv (CombineBlastHits) or
               *_all_exonerate_summary.tsv (CombineExonerateHits). Columns used:
               sample_id, subject_id, strand, subject_start, subject_end
               (minus-strand hits have subject_start > subject_end; samples
               without a hit have subject_id 'no_hits_found').
    --rm-out   RepeatMasker .out tables named <safe_sample_id>_repeatmasker.out,
               where safe_sample_id replaces characters outside [A-Za-z0-9._-]
               with '_' (as in the RepeatMasker process).

Output:
    <safe_query_id>_scaffold_TEs.tsv          one row per TE on the hit scaffold
    <safe_query_id>_scaffold_TEs_summary.tsv  one row per hit: number of TEs on
                                              the scaffold, overlapping TEs and the
                                              nearest upstream / downstream TE

Hit IDs: with Exonerate a sample can carry several gene models; they are named
sample_id_N (N in file order) as in the Exonerate FASTA headers, and plain
sample_id when there is a single model (the same rule applies to BLAST, which
reports one hit per sample). BLAST hits are the merged span of all
HSPs on the best scaffold, so their ends are not necessarily the start/stop
codons.

Distances (edge to edge, in bp) follow the orientation of the hit:
    negative  TE upstream (5') of the hit
    positive  TE downstream (3') of the hit
    0         TE overlaps the hit
For a + strand hit at [start, end], a TE ending before start has distance
te_end - start and a TE starting after end has distance te_start - end. For a
- strand hit the hit's 5' end is its highest coordinate, so a TE at higher
coordinates is upstream (-(te_start - end)) and a TE at lower coordinates is
downstream (+(start - te_end)). te_orientation gives the TE strand relative to
the hit (same / opposite).

Simple_repeat and Low_complexity annotations are not TEs and are skipped.
"""

import argparse
import os
import re


# RepeatMasker classes that are not transposable elements
NON_TE_CLASSES = ('Simple_repeat', 'Low_complexity')

TE_COLUMNS = [
    'sample_id', 'hit_id', 'query_id', 'scaffold', 'hit_start', 'hit_end', 'hit_strand',
    'te_name', 'te_class', 'te_start', 'te_end', 'te_strand', 'te_orientation', 'te_length',
    'position', 'distance_bp', 'perc_div', 'sw_score', 'te_id',
]

SUMMARY_COLUMNS = [
    'sample_id', 'hit_id', 'query_id', 'scaffold', 'hit_start', 'hit_end', 'hit_strand',
    'n_TEs_on_scaffold', 'n_overlapping_TEs', 'overlapping_TEs',
    'nearest_upstream_TE', 'nearest_upstream_distance_bp',
    'nearest_downstream_TE', 'nearest_downstream_distance_bp',
]


def safe_name(name):
    """Filesystem-safe sample ID, as used for the RepeatMasker file names."""
    return re.sub(r'[^A-Za-z0-9._-]', '_', name)


def parse_rm_out(path):
    """Return {scaffold: [TE dict]} from a RepeatMasker .out table.

    Header lines, the 'There were no repetitive sequences detected' message
    and any other line that does not parse are skipped.
    """
    tes = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 15:
                continue
            try:
                sw_score = int(f[0])
                float(f[1])
                start, end = int(f[5]), int(f[6])
            except ValueError:
                continue
            te_class = f[10]
            if te_class.split('/')[0] in NON_TE_CLASSES:
                continue
            tes.setdefault(f[4], []).append({
                'te_name': f[9],
                'te_class': te_class,
                'te_start': min(start, end),
                'te_end': max(start, end),
                'te_strand': '-' if f[8] == 'C' else '+',
                'perc_div': f[1],
                'sw_score': sw_score,
                'te_id': f[14],
            })
    return tes


def read_hits(path):
    """Return [(sample_id, hit dict or None)] in file order.

    The hit dict has hit_id, scaffold, start <= end and strand. Samples without
    a hit are returned with None.
    """
    rows = []
    with open(path) as fh:
        head = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:                 # e.g. 'No BLAST hit files found'
                continue
            rows.append(dict(zip(head, parts)))

    def is_hit(r):
        try:
            int(r.get('subject_start', '')), int(r.get('subject_end', ''))
        except ValueError:
            return False
        return r.get('subject_id', 'no_hits_found') not in ('no_hits_found', 'NA', '')

    n_per_sample = {}
    for r in rows:
        if is_hit(r):
            n_per_sample[r['sample_id']] = n_per_sample.get(r['sample_id'], 0) + 1

    hits, idx = [], {}
    for r in rows:
        sid = r['sample_id']
        if not is_hit(r):
            if sid not in n_per_sample:
                hits.append((sid, None))
                n_per_sample[sid] = 0          # report each no-hit sample once
            continue
        idx[sid] = idx.get(sid, 0) + 1
        # Exonerate: one FASTA record per gene model (sample_id_N when several)
        hit_id = sid if n_per_sample[sid] == 1 else '{}_{}'.format(sid, idx[sid])
        s, e = int(r['subject_start']), int(r['subject_end'])
        strand = r.get('strand', '+')
        if strand not in ('+', '-'):
            strand = '-' if s > e else '+'
        hits.append((sid, {'hit_id': hit_id, 'scaffold': r['subject_id'],
                           'start': min(s, e), 'end': max(s, e), 'strand': strand}))
    return hits


def relate(hit, te):
    """Return (position, distance_bp) of a TE relative to a hit (see docstring)."""
    gs, ge = hit['start'], hit['end']
    ts, te_end = te['te_start'], te['te_end']
    if te_end >= gs and ts <= ge:
        return 'overlap', 0
    if hit['strand'] == '-':
        if ts > ge:
            return 'upstream', -(ts - ge)
        return 'downstream', gs - te_end
    if te_end < gs:
        return 'upstream', te_end - gs
    return 'downstream', ts - ge


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('query_id')
    ap.add_argument('hits_tsv')
    ap.add_argument('--rm-out', nargs='*', default=[],
                    help='RepeatMasker .out tables (<safe_sample_id>_repeatmasker.out)')
    args = ap.parse_args()

    rm_files = {os.path.basename(p)[:-len('_repeatmasker.out')]: p
                for p in args.rm_out if p.endswith('_repeatmasker.out')}
    rm_cache = {}

    te_rows, summary_rows = [], []
    for sid, hit in read_hits(args.hits_tsv):
        base = {'sample_id': sid, 'query_id': args.query_id}
        if hit is None:
            row = dict(base, hit_id='NA', position='no_hit')
            te_rows.append(row)
            summary_rows.append(dict(base, hit_id='NA'))
            continue
        base.update(hit_id=hit['hit_id'], scaffold=hit['scaffold'], hit_start=hit['start'],
                    hit_end=hit['end'], hit_strand=hit['strand'])

        rm_path = rm_files.get(safe_name(sid))
        if rm_path is None:
            # RepeatMasker did not produce a table for this sample
            te_rows.append(dict(base, position='no_repeatmasker_output'))
            summary_rows.append(dict(base))
            continue
        if rm_path not in rm_cache:
            rm_cache[rm_path] = parse_rm_out(rm_path)
        tes = rm_cache[rm_path].get(hit['scaffold'], [])

        if not tes:
            te_rows.append(dict(base, position='none'))
            summary_rows.append(dict(base, n_TEs_on_scaffold=0, n_overlapping_TEs=0))
            continue

        hit_rows = []
        for te in tes:
            position, dist = relate(hit, te)
            hit_rows.append(dict(base, **te,
                                 te_orientation='same' if te['te_strand'] == hit['strand'] else 'opposite',
                                 te_length=te['te_end'] - te['te_start'] + 1,
                                 position=position, distance_bp=dist))
        hit_rows.sort(key=lambda r: (r['distance_bp'], r['te_start']))
        te_rows.extend(hit_rows)

        overlap = [r for r in hit_rows if r['position'] == 'overlap']
        up = [r for r in hit_rows if r['position'] == 'upstream']
        down = [r for r in hit_rows if r['position'] == 'downstream']
        nearest_up = max(up, key=lambda r: r['distance_bp']) if up else {}
        nearest_down = min(down, key=lambda r: r['distance_bp']) if down else {}
        summary_rows.append(dict(
            base, n_TEs_on_scaffold=len(hit_rows), n_overlapping_TEs=len(overlap),
            overlapping_TEs=','.join(r['te_name'] for r in overlap) or 'NA',
            nearest_upstream_TE=nearest_up.get('te_name', 'NA'),
            nearest_upstream_distance_bp=nearest_up.get('distance_bp', 'NA'),
            nearest_downstream_TE=nearest_down.get('te_name', 'NA'),
            nearest_downstream_distance_bp=nearest_down.get('distance_bp', 'NA')))

    # Samples in name order; within a hit, TEs from 5' to 3' (stable sort)
    te_rows.sort(key=lambda r: (r['sample_id'], str(r.get('hit_id'))))
    summary_rows.sort(key=lambda r: (r['sample_id'], str(r.get('hit_id'))))

    safe_qid = re.sub(r'[^\w.\-]', '_', args.query_id)
    for path, cols, rows in (('{}_scaffold_TEs.tsv'.format(safe_qid), TE_COLUMNS, te_rows),
                             ('{}_scaffold_TEs_summary.tsv'.format(safe_qid), SUMMARY_COLUMNS,
                              summary_rows)):
        with open(path, 'w') as out:
            out.write('\t'.join(cols) + '\n')
            for r in rows:
                out.write('\t'.join(str(r.get(c, 'NA')) for c in cols) + '\n')


if __name__ == '__main__':
    main()
