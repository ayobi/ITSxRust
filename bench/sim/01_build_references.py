#!/usr/bin/env python3
"""
01_build_references.py

Parse ITSx positions file, select reads with all 5 regions found,
extract their sequences from the original FASTQ, and write:
  - references.fasta  (the "ground truth" reference sequences)
  - truth.tsv         (known boundary coordinates per reference)

Usage:
    python 01_build_references.py \
        --positions data/itsx_out_positions.txt \
        --fastq data/SRR21494940.fastq \
        --outdir sim_refs \
        --n-refs 2000
"""

import argparse
import os
import random
import re
import sys


def parse_itsx_positions(path):
    """
    Parse ITSx positions file.
    Returns list of dicts with read_id, length, and region boundaries.
    Only returns reads where ALL 5 regions are found.
    """
    records = []
    region_pat = re.compile(r'(\S+):\s+(\d+)-(\d+)')

    with open(path) as f:
        for line in f:
            line = line.rstrip('\n\r')
            if not line:
                continue
            fields = line.split('\t')
            if len(fields) < 7:
                continue

            read_id = fields[0]
            # Skip non-sequence lines (e.g., --END--)
            if read_id.startswith('-') or read_id.startswith('#'):
                continue
            # fields[1] is like "863 bp."
            try:
                length = int(fields[1].strip().split()[0])
            except (ValueError, IndexError):
                continue

            # Parse each region field (fields 2-6): "SSU: 1-30" or "SSU: Not found"
            regions = {}
            skip = False
            for i, name in enumerate(['SSU', 'ITS1', '5.8S', 'ITS2', 'LSU'], start=2):
                if i >= len(fields):
                    skip = True
                    break
                m = region_pat.search(fields[i])
                if m is None:
                    skip = True
                    break
                regions[name] = (int(m.group(2)), int(m.group(3)))

            if skip:
                continue

            records.append({
                'read_id': read_id,
                'length': length,
                'SSU': regions['SSU'],
                'ITS1': regions['ITS1'],
                '5.8S': regions['5.8S'],
                'ITS2': regions['ITS2'],
                'LSU': regions['LSU'],
            })

    return records


def read_fastq_subset(fastq_path, read_ids):
    """
    Read specific reads from a FASTQ file.
    Returns dict: read_id -> (seq, qual)
    """
    wanted = set(read_ids)
    found = {}

    # Handle gzipped
    if fastq_path.endswith('.gz'):
        import gzip
        fh = gzip.open(fastq_path, 'rt')
    else:
        fh = open(fastq_path)

    try:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().rstrip('\n')
            plus = fh.readline()
            qual = fh.readline().rstrip('\n')

            rid = header[1:].split()[0]
            if rid in wanted:
                found[rid] = (seq, qual)
                if len(found) == len(wanted):
                    break
    finally:
        fh.close()

    return found


def main():
    parser = argparse.ArgumentParser(description='Build reference sequences with known ITS boundaries')
    parser.add_argument('--positions', required=True, help='ITSx positions file')
    parser.add_argument('--fastq', required=True, help='Original FASTQ file')
    parser.add_argument('--outdir', default='sim_refs', help='Output directory')
    parser.add_argument('--n-refs', type=int, default=2000,
                        help='Number of reference sequences to sample (default: 2000)')
    parser.add_argument('--seed', type=int, default=42, help='Random seed')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Parse positions
    print(f"Parsing ITSx positions from {args.positions}...")
    records = parse_itsx_positions(args.positions)
    print(f"  Reads with all 5 regions: {len(records)}")

    if len(records) < args.n_refs:
        print(f"  Warning: only {len(records)} available, using all")
        sampled = records
    else:
        random.seed(args.seed)
        sampled = random.sample(records, args.n_refs)
    print(f"  Sampled: {len(sampled)}")

    # Extract sequences from FASTQ
    print(f"Extracting sequences from {args.fastq}...")
    read_ids = [r['read_id'] for r in sampled]
    seqs = read_fastq_subset(args.fastq, read_ids)
    print(f"  Found: {len(seqs)}/{len(read_ids)}")

    missing = [rid for rid in read_ids if rid not in seqs]
    if missing:
        print(f"  Warning: {len(missing)} reads not found in FASTQ, skipping them")
        sampled = [r for r in sampled if r['read_id'] in seqs]

    # Write reference FASTA
    ref_path = os.path.join(args.outdir, 'references.fasta')
    with open(ref_path, 'w') as f:
        for rec in sampled:
            rid = rec['read_id']
            seq = seqs[rid][0]
            f.write(f">{rid}\n{seq}\n")
    print(f"  Wrote {len(sampled)} references to {ref_path}")

    # Write truth table
    truth_path = os.path.join(args.outdir, 'truth.tsv')
    with open(truth_path, 'w') as f:
        header = [
            'ref_id', 'ref_len',
            'ssu_start', 'ssu_end',
            'its1_start', 'its1_end',
            's58_start', 's58_end',
            'its2_start', 'its2_end',
            'lsu_start', 'lsu_end',
        ]
        f.write('\t'.join(header) + '\n')
        for rec in sampled:
            rid = rec['read_id']
            seq_len = len(seqs[rid][0])
            cols = [
                rid, str(seq_len),
                str(rec['SSU'][0]), str(rec['SSU'][1]),
                str(rec['ITS1'][0]), str(rec['ITS1'][1]),
                str(rec['5.8S'][0]), str(rec['5.8S'][1]),
                str(rec['ITS2'][0]), str(rec['ITS2'][1]),
                str(rec['LSU'][0]), str(rec['LSU'][1]),
            ]
            f.write('\t'.join(cols) + '\n')
    print(f"  Wrote truth table to {truth_path}")

    # Print summary stats
    its1_lens = [r['ITS1'][1] - r['ITS1'][0] + 1 for r in sampled]
    its2_lens = [r['ITS2'][1] - r['ITS2'][0] + 1 for r in sampled]
    full_lens = [r['ITS2'][1] - r['ITS1'][0] + 1 for r in sampled]  # ITS1 through ITS2 inclusive
    print(f"\n  Reference stats:")
    print(f"    ITS1 length: median={sorted(its1_lens)[len(its1_lens)//2]}, "
          f"range=[{min(its1_lens)},{max(its1_lens)}]")
    print(f"    ITS2 length: median={sorted(its2_lens)[len(its2_lens)//2]}, "
          f"range=[{min(its2_lens)},{max(its2_lens)}]")
    print(f"    Full ITS:    median={sorted(full_lens)[len(full_lens)//2]}, "
          f"range=[{min(full_lens)},{max(full_lens)}]")


if __name__ == '__main__':
    main()