#!/usr/bin/env python3
"""
Paralog check of one sample: which BLAST hits of the query belong to the query
itself, and which to one of its known paralogs?

The query and its paralogs were BLASTed against the sample's assembly with the
same program and settings. Every hit (HSP) is assigned to the reference that
scores best there:

    ownership   an HSP belongs to its reference unless an HSP of another
                reference overlaps it on the scaffold with a higher bit score
                (ties go to the query, then to the lower paralog number)
    copy        the owned HSPs of one reference on one scaffold, split where an
                HSP owned by another reference lies between them (tandem
                arrays); the extracted hit span is one query copy

The query copy that the pipeline extracts is chosen as before (the scaffold of
the best-scoring HSP), but among the query-owned HSPs only, and it must also
win at the locus level: summed over all HSPs overlapping its span, no paralog
may score higher than the query. Copies that fail are skipped (listed in
skipped_paralog_loci) and the next query copy is tried.

    call = query      a query copy was found and is extracted
           ambiguous  as 'query', but the second-best reference reaches
                      --ambiguous-ratio of the query's score at the locus
           paralog    the query hits only loci that belong to a paralog;
                      nothing is extracted for this sample
           no_hit     the query has no BLAST hit

Scores at a locus: bit scores of non-overlapping HSPs summed (greedy by bit
score; an HSP overlapping an accepted one by more than 20 % in query or
scaffold coordinates is skipped), their length-weighted identity, and the
share of the reference they cover.

Usage:
    paralog_check.py <sample_id> <query_id> <query.fasta> <paralogs.fasta>
        <query_blast.tsv> <paralog_blast.tsv> [--ambiguous-ratio 0.95]

Outputs:
    <sample>_query_blast_results.tsv   the query HSPs of the extracted copy, in
                                       BLAST tabular format (input of the hit
                                       extraction), then one line
                                       '# paralog_copy <label> <scaffold>
                                       <start> <end>' per paralog copy next to
                                       it (see 'paralog copies next to the
                                       extracted query copy'); for a 'paralog'
                                       call it holds only a '# paralog_only'
                                       line
    <sample>_paralog_check.tsv         one row: the call and every reference's
                                       score at the locus and its own best copy
    <sample>_paralog_loci.tsv          one row per copy of every reference
"""

import argparse
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


def read_hsps(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if line.startswith('#') or len(f) < 12:
                continue
            s0, s1 = sorted((int(f[8]), int(f[9])))
            q0, q1 = sorted((int(f[6]), int(f[7])))
            rows.append({'ref': f[0], 'scaffold': f[1], 'pident': float(f[2]),
                         'length': int(f[3]), 'q0': q0, 'q1': q1, 's0': s0, 's1': s1,
                         'bitscore': float(f[11]), 'line': line.rstrip('\n')})
    return rows


def overlap(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0) + 1)


def score(hsps, ref_len):
    """Score a set of HSPs of one reference (see module doc)."""
    kept = []
    for h in sorted(hsps, key=lambda h: -h['bitscore']):
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


def fmt(x, digits=1):
    return 'NA' if x is None else f'{x:.{digits}f}'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('sample_id')
    ap.add_argument('query_id')
    ap.add_argument('query_fasta')
    ap.add_argument('paralog_fasta')
    ap.add_argument('query_blast')
    ap.add_argument('paralog_blast')
    ap.add_argument('--ambiguous-ratio', type=float, default=0.95,
                    help="call 'ambiguous' when second best score >= ratio x query score")
    args = ap.parse_args()
    sid = args.sample_id

    query_name, query_len = fasta_names(args.query_fasta)[0]
    paralogs = [(name, n) for name, n in fasta_names(args.paralog_fasta)]
    paralogs.sort(key=lambda r: int(re.search(r'_paralog(\d+)$', r[0]).group(1)))
    # (name, length, column label); list order is the tie-break order
    refs = [(query_name, query_len, 'query')] + \
        [(name, n, re.sub(r'^.*_(paralog\d+)$', r'\1', name)) for name, n in paralogs]
    rank = {name: i for i, (name, _, _) in enumerate(refs)}
    length = {name: n for name, n, _ in refs}
    label = {name: lab for name, _, lab in refs}

    hsps = [h for h in read_hsps(args.query_blast) if h['ref'] == query_name] + \
        [h for h in read_hsps(args.paralog_blast) if h['ref'] in rank and h['ref'] != query_name]

    # --- ownership of every HSP ---------------------------------------------
    by_scaffold = {}
    for h in hsps:
        by_scaffold.setdefault(h['scaffold'], []).append(h)
    for h in hsps:
        h['owner'] = h['ref']
        for c in by_scaffold[h['scaffold']]:
            if c['ref'] == h['ref'] or not overlap(h['s0'], h['s1'], c['s0'], c['s1']):
                continue
            if (c['bitscore'], -rank[c['ref']]) > (h['bitscore'], -rank[h['ref']]):
                h['owner'] = None
                break

    # --- copies: owned HSPs of one reference on one scaffold, split where an
    # HSP owned by another reference lies in between (tandem arrays) ---------
    copies = {}                      # ref -> [copy dict]
    for name, _, _ in refs:
        groups = {}
        for h in hsps:
            if h['owner'] == name:
                groups.setdefault(h['scaffold'], []).append(h)
        cl = []
        for scaffold, hs in groups.items():
            others = [h for h in by_scaffold[scaffold] if h['owner'] not in (None, name)]
            parts = []
            for h in sorted(hs, key=lambda h: h['s0']):
                if parts:
                    end = max(x['s1'] for x in parts[-1])
                    if not any(o['s0'] > end and o['s1'] < h['s0'] for o in others):
                        parts[-1].append(h)
                        continue
                parts.append([h])
            for part in parts:
                cl.append(dict(score(part, length[name]), ref=name, scaffold=scaffold, hsps=part,
                               start=min(h['s0'] for h in part), end=max(h['s1'] for h in part),
                               top=max(h['bitscore'] for h in part)))
        # Ordered as the hit extraction picks: best single HSP first
        cl.sort(key=lambda c: (-c['top'], -c['bitscore']))
        copies[name] = cl

    def locus_scores(scaffold, lo, hi):
        """Score of every reference over all its HSPs overlapping scaffold:lo-hi."""
        out = {}
        for name, n, _ in refs:
            out[name] = score([h for h in hsps if h['ref'] == name and h['scaffold'] == scaffold
                               and overlap(h['s0'], h['s1'], lo, hi)], n)
        return out

    def ranked(scores):
        return sorted(scores, key=lambda name: (-scores[name]['bitscore'], rank[name]))

    # --- the query copy to extract ------------------------------------------
    chosen, at_hit, skipped = None, None, []
    for c in copies[query_name]:
        sc = locus_scores(c['scaffold'], c['start'], c['end'])
        best = ranked(sc)[0]
        if best == query_name:
            chosen, at_hit = c, sc
            break
        skipped.append(f"{c['scaffold']}:{c['start']}-{c['end']}({label[best]})")

    query_hsps = [h for h in hsps if h['ref'] == query_name]

    # --- paralog copies next to the extracted query copy --------------------
    # The hit extraction keeps its padding off them and the CDS integrity check
    # does not search them for missing query segments, so that a neighbouring
    # paralog is neither modelled (alone or fused to a truncated query copy)
    # nor taken for part of the query gene.
    flanking = []
    if chosen is not None:
        for name, _, lab in refs[1:]:
            flanking += [(lab, c) for c in copies[name] if c['scaffold'] == chosen['scaffold']
                         and not overlap(c['start'], c['end'], chosen['start'], chosen['end'])]

    if chosen is not None:
        locus = (chosen['scaffold'], chosen['start'], chosen['end'])
    elif query_hsps:
        # Paralog-only sample: report the locus of the query's best HSP
        top = max(query_hsps, key=lambda h: h['bitscore'])
        on = [h for h in query_hsps if h['scaffold'] == top['scaffold']]
        locus = (top['scaffold'], min(h['s0'] for h in on), max(h['s1'] for h in on))
        at_hit = locus_scores(*locus)
    else:
        locus = None

    row = {'sample_id': sid, 'query_id': args.query_id}
    if locus is None:
        row['call'] = 'no_hit'
    else:
        order = ranked(at_hit)
        b1, b2 = at_hit[order[0]]['bitscore'], at_hit[order[1]]['bitscore']
        ratio = b2 / b1 if b1 else None
        if chosen is None:
            row['call'] = 'paralog'
        elif ratio is not None and ratio >= args.ambiguous_ratio:
            row['call'] = 'ambiguous'
        else:
            row['call'] = 'query'
        row.update(hit_scaffold=locus[0], hit_start=locus[1], hit_end=locus[2],
                   closest_reference=order[0], best_bitscore=fmt(b1),
                   second_reference=order[1], second_bitscore=fmt(b2), score_ratio=fmt(ratio, 3))
        for name, _, lab in refs:
            s = at_hit[name]
            row.update({f'{lab}_bitscore_at_hit': fmt(s['bitscore']),
                        f'{lab}_identity_pct_at_hit': fmt(s['identity']),
                        f'{lab}_coverage_pct_at_hit': fmt(s['coverage'])})
    row['skipped_paralog_loci'] = ';'.join(skipped) or 'none'
    for name, _, lab in refs:
        cl = copies[name]
        row[f'{lab}_n_copies'] = len(cl)
        if cl:
            c = max(cl, key=lambda c: c['bitscore'])
            row.update({f'{lab}_best_copy': f"{c['scaffold']}:{c['start']}-{c['end']}",
                        f'{lab}_best_copy_bitscore': fmt(c['bitscore']),
                        f'{lab}_best_copy_identity_pct': fmt(c['identity']),
                        f'{lab}_best_copy_coverage_pct': fmt(c['coverage'])})
        else:
            row[f'{lab}_best_copy'] = 'no_hit'

    head = ['sample_id', 'query_id', 'call', 'hit_scaffold', 'hit_start', 'hit_end',
            'closest_reference', 'best_bitscore', 'second_reference', 'second_bitscore',
            'score_ratio']
    for _, _, lab in refs:
        head += [f'{lab}_bitscore_at_hit', f'{lab}_identity_pct_at_hit', f'{lab}_coverage_pct_at_hit']
    head += ['skipped_paralog_loci']
    for _, _, lab in refs:
        head += [f'{lab}_n_copies', f'{lab}_best_copy', f'{lab}_best_copy_bitscore',
                 f'{lab}_best_copy_identity_pct', f'{lab}_best_copy_coverage_pct']
    with open(f'{sid}_paralog_check.tsv', 'w') as out:
        out.write('\t'.join(head) + '\n')
        out.write('\t'.join(str(row.get(c, 'NA')) for c in head) + '\n')

    with open(f'{sid}_paralog_loci.tsv', 'w') as out:
        out.write('\t'.join(['sample_id', 'query_id', 'reference', 'reference_label', 'scaffold',
                             'start', 'end', 'n_hsps', 'bitscore', 'identity_pct', 'coverage_pct',
                             'extracted']) + '\n')
        for name, _, lab in refs:
            for c in copies[name]:
                out.write('\t'.join(map(str, [
                    sid, args.query_id, name, lab, c['scaffold'], c['start'], c['end'],
                    len(c['hsps']), fmt(c['bitscore']), fmt(c['identity']), fmt(c['coverage']),
                    'yes' if c is chosen else 'no'])) + '\n')

    with open(f'{sid}_query_blast_results.tsv', 'w') as out:
        if chosen is not None:
            for h in chosen['hsps']:
                out.write(h['line'] + '\n')
            for lab, c in sorted(flanking, key=lambda x: x[1]['start']):
                out.write(f"# paralog_copy\t{lab}\t{c['scaffold']}\t{c['start']}\t{c['end']}\n")
        elif row['call'] == 'paralog':
            out.write(f"# paralog_only {row['closest_reference']}\n")


if __name__ == '__main__':
    main()
