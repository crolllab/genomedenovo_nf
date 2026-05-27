#!/usr/bin/env python3
"""
concat_busco_alignments.py

Concatenate per-gene MAFFT alignments into a supermatrix FASTA and
a NEXUS partition file for IQ-tree. Each gene alignment file must be
named <gene_id>_aligned.fasta and contain one sequence per sample
(FASTA header = sample_id, as produced by the BUSCO module).

Outputs:
  supermatrix.fasta  — concatenated alignment, one sequence per sample
  partitions.nex     — NEXUS charset file for per-gene model selection

Samples missing from a gene are filled with gap characters (-).
Only genes present in ≥ 2 samples are included.
"""

import glob
import os
import sys


def parse_fasta(filepath):
    """Return dict of {header: sequence} from a FASTA file."""
    seqs = {}
    current_id = None
    current_seq = []
    with open(filepath) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if current_id is not None:
                    seqs[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id is not None:
        seqs[current_id] = ''.join(current_seq)
    return seqs


def main():
    aln_files = sorted(glob.glob('*_aligned.fasta'))
    if not aln_files:
        sys.exit("ERROR: No *_aligned.fasta files found in working directory.")

    # Parse all gene alignments
    gene_seqs = {}   # {gene_id: {sample_id: sequence}}
    for fpath in aln_files:
        gene_id = os.path.basename(fpath).replace('_aligned.fasta', '')
        seqs = parse_fasta(fpath)
        if len(seqs) < 2:
            print(f"  Skipping {gene_id}: only {len(seqs)} sample(s) — need ≥ 2", flush=True)
            continue
        gene_seqs[gene_id] = seqs

    if not gene_seqs:
        sys.exit("ERROR: No genes with ≥ 2 samples found. Cannot build supermatrix.")

    # Collect all samples across all valid genes
    all_samples = sorted({s for seqs in gene_seqs.values() for s in seqs})
    valid_genes = sorted(gene_seqs.keys())

    # Build supermatrix: fill missing genes with gap characters
    supermatrix = {s: [] for s in all_samples}
    partitions = []
    position = 1

    for gene_id in valid_genes:
        seqs = gene_seqs[gene_id]
        aln_len = len(next(iter(seqs.values())))
        gap_seq = '-' * aln_len

        for sample in all_samples:
            supermatrix[sample].append(seqs.get(sample, gap_seq))

        partitions.append((gene_id, position, position + aln_len - 1))
        position += aln_len

    # Write supermatrix FASTA
    with open('supermatrix.fasta', 'w') as out:
        for sample in all_samples:
            seq = ''.join(supermatrix[sample])
            out.write(f'>{sample}\n{seq}\n')

    # Write NEXUS partition file
    with open('partitions.nex', 'w') as out:
        out.write('#nexus\nbegin sets;\n')
        for gene_id, start, end in partitions:
            out.write(f'    charset {gene_id} = {start}-{end};\n')
        out.write('end;\n')

    total_sites = position - 1
    print(
        f"Supermatrix written: {len(all_samples)} taxa, "
        f"{total_sites} sites, {len(valid_genes)} genes",
        flush=True,
    )


if __name__ == '__main__':
    main()
