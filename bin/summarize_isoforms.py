#!/usr/bin/env python3
"""
Assign isoform IDs to unique CDS and protein haplotypes derived from
MAFFT-aligned FASTA files, and report which isoforms are full-length or
truncated (from the CDS integrity tables).

Usage:
    summarize_isoforms.py <query_id> <cds_aligned.fasta> <protein_aligned.fasta>
        [--integrity <*_cds_integrity.tsv> ...] [--events <*_cds_integrity_events.tsv> ...]

Output:
    <safe_query_id>_isoforms.tsv                   intact isoforms only
                                                   (columns: sample_id, cds_isoform, protein_isoform,
                                                   cds_sequence, protein_sequence)
    <safe_query_id>_isoforms_incl_truncated.tsv    all sequences, with integrity status and
                                                   mutation description; samples without an
                                                   Exonerate hit are listed with status no_hit
    <safe_query_id>_isoform_mutations.tsv          one row per mutation event, with isoform IDs
    <safe_query_id>_cds_isoforms.fasta             one ungapped sequence per unique CDS isoform
    <safe_query_id>_protein_isoforms.fasta         one ungapped sequence per unique protein isoform
    <safe_query_id>_cds_isoforms_aligned.fasta     one aligned row per CDS isoform (all-gap
                                                   columns removed), input for the isoform tree
    <safe_query_id>_protein_isoforms_aligned.fasta same for the protein isoforms

Numbering is shared by all outputs and follows isoform frequency:
  1. isoforms carried by at least one intact sample sequence, by the number of
     sample sequences carrying them (descending; the reference is not counted);
  2. the query reference, when its sequence is not found in any intact sample;
  3. isoforms found only in truncated (or stop-lost, 'extended') sequences,
     again by number of sample sequences.
Ties keep the order of first appearance in the alignment. The intact-only
table (reference row included) therefore has consecutive numbers, and
truncated isoforms continue the series in the combined table.

Haplotype identity is determined from the full aligned sequence (gaps included),
uppercased. MAFFT reverse-complement headers prefixed with '_R_' are stripped.
Protein and CDS records of a sequence share the same header (sample_id, or
sample_id_N for several gene models). The query reference is the record whose
header equals the query ID or its first word (a translated nucleotide query is
named by the first word of its header).
"""

import argparse
import re


# Statuses excluded from the intact-only table
NOT_INTACT = {'truncated', 'extended'}


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


def read_tsv(paths):
    """Concatenate TSV files that share a header; return (columns, rows as dicts)."""
    columns, rows = None, []
    for path in paths:
        with open(path) as fh:
            head = fh.readline().rstrip('\n').split('\t')
            columns = columns or head
            for line in fh:
                if line.strip():
                    rows.append(dict(zip(head, line.rstrip('\n').split('\t'))))
    return columns or [], rows


def assign_isoforms(records, status_of):
    """Map each record to an isoform number (see module docstring for the order).

    records   : [(record_id, aligned_seq)]
    status_of : {record_id: status}
    Returns ({record_id: isoform_int}, {isoform_int: [statuses]}).
    """
    groups = {}          # seq -> [record_ids], in first-seen order
    for rid, seq in records:
        groups.setdefault(seq, []).append(rid)

    def rank(item):
        first_seen, ids = item
        statuses = [status_of.get(r, 'NA') for r in ids]
        samples = [s for s in statuses if s != 'reference']
        if any(s not in NOT_INTACT for s in samples):
            category = 0                 # carried by an intact sample sequence
        elif 'reference' in statuses:
            category = 1                 # reference not found in any intact sample
        else:
            category = 2                 # truncated / extended sequences only
        return (category, -len(samples), first_seen)

    ordered = sorted(enumerate(groups.values()), key=rank)
    numbers, statuses = {}, {}
    for num, (_, ids) in enumerate(ordered, start=1):
        statuses[num] = [status_of.get(r, 'NA') for r in ids]
        for r in ids:
            numbers[r] = num
    return numbers, statuses


def write_isoform_fasta(path, prefix, records, isoforms, statuses):
    """Write one ungapped sequence per isoform, numbered as in the TSVs.

    The header records how many sequences carry the isoform and their status,
    e.g. '>cds_isoform_1 n_samples=3 status=full_length'.
    """
    seqs = {}
    for rid, seq in records:
        seqs.setdefault(isoforms[rid], seq.replace('-', ''))
    with open(path, 'w') as out:
        for num in sorted(seqs):
            status = ','.join(sorted(set(statuses[num])))
            n = sum(s != 'reference' for s in statuses[num])
            out.write('>{}_{} n_samples={} status={}\n'.format(prefix, num, n, status))
            seq = seqs[num]
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + '\n')


def write_aligned_isoform_fasta(path, prefix, records, isoforms):
    """Write one aligned row per isoform (input for the isoform tree).

    Rows come from the existing MAFFT alignment; columns that are gaps in
    every isoform are dropped. Names are '<prefix>_N' without annotations,
    as IQ-TREE uses the whole header as the taxon name.
    """
    rows = {}
    for rid, seq in records:
        rows.setdefault(isoforms[rid], seq)
    nums = sorted(rows)
    keep = [i for i in range(max((len(rows[n]) for n in nums), default=0))
            if any(i < len(rows[n]) and rows[n][i] != '-' for n in nums)]
    with open(path, 'w') as out:
        for num in nums:
            seq = ''.join(rows[num][i] if i < len(rows[num]) else '-' for i in keep)
            out.write('>{}_{}\n'.format(prefix, num))
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + '\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('query_id')
    ap.add_argument('cds_fasta')
    ap.add_argument('protein_fasta')
    ap.add_argument('--integrity', nargs='*', default=[],
                    help='per-sample *_cds_integrity.tsv files for this query')
    ap.add_argument('--events', nargs='*', default=[],
                    help='per-sample *_cds_integrity_events.tsv files for this query')
    args = ap.parse_args()

    cds_records = read_fasta(args.cds_fasta)
    prot_all = read_fasta(args.protein_fasta)

    # --- sample vs reference records ------------------------------------------
    # Protein records carry the same header as their CDS (sample_id or
    # sample_id_N). The query reference keeps the full query header (protein
    # query, or nucleotide query in the CDS file) or is named by its first word
    # (translated nucleotide query in the protein file).
    query_words = args.query_id.split()
    ref_names = {args.query_id.strip()} | ({query_words[0]} if query_words else set())

    def is_ref(h):
        return h in ref_names

    ref_headers = [h for h, _ in cds_records + prot_all if is_ref(h)]
    ref_id = ref_headers[0] if ref_headers else None

    samples = list(dict.fromkeys(h for h, _ in cds_records + prot_all if not is_ref(h)))

    # Records keyed by sample id; both reference records map to one reference id
    cds_records = [(ref_id if is_ref(h) else h, seq) for h, seq in cds_records]
    prot_records = [(ref_id if is_ref(h) else h, seq) for h, seq in prot_all]

    # --- integrity status per sequence ----------------------------------------
    int_cols, int_rows = read_tsv(args.integrity)
    integrity = {r['sequence_id']: r for r in int_rows if r.get('status') != 'no_hit'}
    no_hit = [r for r in int_rows if r.get('status') == 'no_hit']

    status_of = {sid: integrity.get(sid, {}).get('status', 'NA') for sid in samples}
    if ref_id is not None:
        status_of[ref_id] = 'reference'

    cds_iso, cds_status = assign_isoforms(cds_records, status_of)
    prot_iso, prot_status = assign_isoforms(prot_records, status_of)

    # Row order: reference, then samples as they appear in the CDS alignment
    order = [ref_id] if ref_id is not None else []
    for rid, _ in cds_records + prot_records:
        if rid not in order:
            order.append(rid)

    safe_qid = re.sub(r'[^\w.\-]', '_', args.query_id)

    # Ungapped sequence per record, for the sequence columns of both tables
    cds_seq = {rid: seq.replace('-', '') for rid, seq in cds_records}
    prot_seq = {rid: seq.replace('-', '') for rid, seq in prot_records}

    def seq_cols(sid):
        return [cds_seq.get(sid) or 'NA', prot_seq.get(sid) or 'NA']

    # --- intact-only table ------------------------------------------------------
    with open('{}_isoforms.tsv'.format(safe_qid), 'w') as out:
        out.write('sample_id\tcds_isoform\tprotein_isoform\tcds_sequence\tprotein_sequence\n')
        for sid in order:
            if status_of.get(sid) in NOT_INTACT:
                continue
            out.write('\t'.join([sid, str(cds_iso.get(sid, 'NA')), str(prot_iso.get(sid, 'NA'))]
                                + seq_cols(sid)) + '\n')

    # --- combined table with integrity information ----------------------------
    extra = [c for c in int_cols if c not in ('sample_id', 'query_id', 'sequence_id', 'status')]
    with open('{}_isoforms_incl_truncated.tsv'.format(safe_qid), 'w') as out:
        out.write('\t'.join(['sample_id', 'cds_isoform', 'protein_isoform', 'status'] + extra
                            + ['cds_sequence', 'protein_sequence']) + '\n')
        for sid in order:
            row = integrity.get(sid, {})
            out.write('\t'.join([sid, str(cds_iso.get(sid, 'NA')), str(prot_iso.get(sid, 'NA')),
                                 status_of.get(sid, 'NA')]
                                + [row.get(c, 'NA') for c in extra] + seq_cols(sid)) + '\n')
        for row in no_hit:
            out.write('\t'.join([row['sample_id'], 'NA', 'NA', 'no_hit']
                                + ['NA'] * (len(extra) + 2)) + '\n')

    # --- mutation events with isoform IDs --------------------------------------
    ev_cols, ev_rows = read_tsv(args.events)
    head = ['sample_id', 'cds_isoform', 'protein_isoform'] + \
        [c for c in ev_cols if c not in ('sample_id', 'query_id', 'sequence_id')]
    rank = {sid: i for i, sid in enumerate(order)}
    ev_rows.sort(key=lambda r: rank.get(r['sequence_id'], len(rank)))
    with open('{}_isoform_mutations.tsv'.format(safe_qid), 'w') as out:
        out.write('\t'.join(head) + '\n')
        for r in ev_rows:
            sid = r['sequence_id']
            out.write('\t'.join([sid, str(cds_iso.get(sid, 'NA')), str(prot_iso.get(sid, 'NA'))]
                                + [r.get(c, 'NA') for c in head[3:]]) + '\n')

    write_isoform_fasta('{}_cds_isoforms.fasta'.format(safe_qid),
                        'cds_isoform', cds_records, cds_iso, cds_status)
    write_isoform_fasta('{}_protein_isoforms.fasta'.format(safe_qid),
                        'protein_isoform', prot_records, prot_iso, prot_status)
    write_aligned_isoform_fasta('{}_cds_isoforms_aligned.fasta'.format(safe_qid),
                                'cds_isoform', cds_records, cds_iso)
    write_aligned_isoform_fasta('{}_protein_isoforms_aligned.fasta'.format(safe_qid),
                                'protein_isoform', prot_records, prot_iso)


if __name__ == '__main__':
    main()
