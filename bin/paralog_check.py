#!/usr/bin/env python3
"""
Paralog check: is the hit reported for a query in each sample closer to the
query or to one of its known paralogs?

The hit locus of a sample is the span reported in the query's BLAST hit
summary (best scaffold, all query hits on it merged; the region extracted for
the alignments and Exonerate). The query and every paralog were BLASTed
against the same assembly with the same program and settings; their hits
overlapping the locus are scored per reference:

    bitscore      sum of the bit scores of a set of non-overlapping hits
                  (greedy by bit score; hits overlapping an accepted one by
                  more than 20 % in query or scaffold coordinates are skipped)
    identity_pct  alignment-length-weighted identity of these hits
    coverage_pct  share of the reference covered by them

The closest reference is the one with the highest bit score. The call is
'query' or 'paralog', or 'ambiguous' when the second best reaches
--ambiguous-ratio of the best score. Each paralog's own best locus in the
sample (same rule as for the query) is listed as well, which shows whether the
sample carries a separate copy for it.

Usage:
    paralog_check.py <query_id> <query.fasta> <paralogs.fasta> <hits_summary.tsv>
        --query-blast <sample>_blast_results.tsv ...
        --paralog-blast <sample>_paralog_blast_results.tsv ...
        [--ambiguous-ratio 0.95]

Output:
    <safe_query_id>_paralog_check.tsv   one row per sample
"""

import argparse
import os
import re

MAX_OVERLAP = 0.2


def fasta_names(path):
    """Return [(first word of header, sequence length)] in file order."""
    records, name, n = [], None, 0
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                if name is not None:
                    records.append((name, n))
                name, n = line[1:].split()[0], 0
            elif name is not None:
                n += len(line.rstrip('*'))
    if name is not None:
        records.append((name, n))
    return records


def read_hsps(paths, suffix):
    """{sample_id: [hsp dict]} from BLAST outfmt-6 tables named <sample><suffix>."""
    hsps = {}
    for path in paths:
        base = os.path.basename(path)
        if not base.endswith(suffix):
            continue
        sample = base[:-len(suffix)]
        rows = hsps.setdefault(sample, [])
        with open(path) as fh:
            for line in fh:
                f = line.rstrip('\n').split('\t')
                if len(f) < 12:
                    continue
                s0, s1 = sorted((int(f[8]), int(f[9])))
                q0, q1 = sorted((int(f[6]), int(f[7])))
                rows.append({'ref': f[0], 'scaffold': f[1], 'pident': float(f[2]),
                             'length': int(f[3]), 'q0': q0, 'q1': q1, 's0': s0, 's1': s1,
                             'bitscore': float(f[11])})
    return hsps


def overlap(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0) + 1)


def score_locus(hsps, scaffold, lo, hi, ref_len):
    """Score one reference's hits overlapping scaffold:lo-hi (see module doc)."""
    cand = sorted((h for h in hsps if h['scaffold'] == scaffold and overlap(h['s0'], h['s1'], lo, hi)),
                  key=lambda h: -h['bitscore'])
    kept = []
    for h in cand:
        qlen, slen = h['q1'] - h['q0'] + 1, h['s1'] - h['s0'] + 1
        if any(overlap(h['q0'], h['q1'], k['q0'], k['q1']) > MAX_OVERLAP * qlen or
               overlap(h['s0'], h['s1'], k['s0'], k['s1']) > MAX_OVERLAP * slen for k in kept):
            continue
        kept.append(h)
    if not kept:
        return {'bitscore': 0.0, 'identity': None, 'coverage': 0.0}
    aln = sum(h['length'] for h in kept)
    covered = set()
    for h in kept:
        covered.update(range(h['q0'], h['q1'] + 1))
    return {'bitscore': sum(h['bitscore'] for h in kept),
            'identity': sum(h['pident'] * h['length'] for h in kept) / aln if aln else None,
            'coverage': 100.0 * min(len(covered), ref_len) / ref_len if ref_len else None}


def best_locus(hsps, ref):
    """Best scaffold (top bit score) and merged span of all hits of *ref* on it."""
    rows = [h for h in hsps if h['ref'] == ref]
    if not rows:
        return None
    top = max(rows, key=lambda h: h['bitscore'])['scaffold']
    on = [h for h in rows if h['scaffold'] == top]
    return top, min(h['s0'] for h in on), max(h['s1'] for h in on)


def fmt(x, digits=1):
    return 'NA' if x is None else f'{x:.{digits}f}'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('query_id')
    ap.add_argument('query_fasta')
    ap.add_argument('paralog_fasta')
    ap.add_argument('hits_summary')
    ap.add_argument('--query-blast', nargs='*', default=[])
    ap.add_argument('--paralog-blast', nargs='*', default=[])
    ap.add_argument('--ambiguous-ratio', type=float, default=0.95,
                    help="call 'ambiguous' when second best score >= ratio x best score")
    args = ap.parse_args()

    query_name, query_len = fasta_names(args.query_fasta)[0]
    paralogs = fasta_names(args.paralog_fasta)
    paralogs.sort(key=lambda r: int(re.search(r'_paralog(\d+)$', r[0]).group(1))
                  if re.search(r'_paralog(\d+)$', r[0]) else 0)
    refs = [(query_name, query_len, 'query')] + \
        [(name, n, re.sub(r'^.*_(paralog\d+)$', r'\1', name)) for name, n in paralogs]

    query_hsps = read_hsps(args.query_blast, '_blast_results.tsv')
    paralog_hsps = read_hsps(args.paralog_blast, '_paralog_blast_results.tsv')

    head = ['sample_id', 'query_id', 'hit_scaffold', 'hit_start', 'hit_end',
            'closest_reference', 'call', 'best_bitscore', 'second_reference',
            'second_bitscore', 'score_ratio']
    for _, _, label in refs:
        head += [f'{label}_bitscore', f'{label}_identity_pct', f'{label}_coverage_pct']
    for _, _, label in refs[1:]:
        head += [f'{label}_best_locus', f'{label}_best_locus_is_hit']

    rows = []
    with open(args.hits_summary) as fh:
        cols = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            r = dict(zip(cols, line.rstrip('\n').split('\t')))
            if 'sample_id' not in r or 'subject_start' not in r:
                continue
            rows.append(r)

    safe_qid = re.sub(r'[^\w.\-]', '_', args.query_id)
    with open(f'{safe_qid}_paralog_check.tsv', 'w') as out:
        out.write('\t'.join(head) + '\n')
        for r in rows:
            sid = r['sample_id']
            hsps = query_hsps.get(sid, []) + paralog_hsps.get(sid, [])
            row = {c: 'NA' for c in head}
            row.update(sample_id=sid, query_id=args.query_id)
            try:
                lo, hi = sorted((int(r['subject_start']), int(r['subject_end'])))
            except ValueError:
                row['call'] = 'no_hit'
                out.write('\t'.join(row[c] for c in head) + '\n')
                continue
            scaffold = r['subject_id']
            row.update(hit_scaffold=scaffold, hit_start=str(lo), hit_end=str(hi))

            scores = []
            for name, n, label in refs:
                s = score_locus([h for h in hsps if h['ref'] == name], scaffold, lo, hi, n)
                scores.append((s['bitscore'], name))
                row[f'{label}_bitscore'] = fmt(s['bitscore'])
                row[f'{label}_identity_pct'] = fmt(s['identity'])
                row[f'{label}_coverage_pct'] = fmt(s['coverage'])
            for name, _, label in refs[1:]:
                loc = best_locus(paralog_hsps.get(sid, []), name)
                if loc is None:
                    row[f'{label}_best_locus'] = 'no_hit'
                    continue
                row[f'{label}_best_locus'] = f'{loc[0]}:{loc[1]}-{loc[2]}'
                row[f'{label}_best_locus_is_hit'] = \
                    'yes' if loc[0] == scaffold and overlap(loc[1], loc[2], lo, hi) else 'no'

            # Highest score first; ties keep the query ahead of the paralogs
            ranked = sorted(scores, key=lambda x: -x[0])
            (b1, n1), (b2, n2) = ranked[0], ranked[1] if len(ranked) > 1 else (0.0, 'NA')
            ratio = b2 / b1 if b1 else None
            row.update(closest_reference=n1, best_bitscore=fmt(b1), second_reference=n2,
                       second_bitscore=fmt(b2), score_ratio=fmt(ratio, 3))
            if b1 <= 0:
                row['call'] = 'no_hit'
            elif ratio is not None and ratio >= args.ambiguous_ratio:
                row['call'] = 'ambiguous'
            else:
                row['call'] = 'query' if n1 == query_name else 'paralog'
            out.write('\t'.join(row[c] for c in head) + '\n')


if __name__ == '__main__':
    main()
