#!/usr/bin/env python3
"""
04_eval_boundaries.py

Evaluate boundary accuracy by aligning each tool's extracted sequences
back to the original references with minimap2, then comparing alignment
boundaries on the reference to the known true ITS boundaries.

For each extracted sequence:
  1. Align to reference FASTA with minimap2
  2. The alignment start/end on the reference = the tool's called boundary
  3. Compare to truth.tsv for that reference
  4. Boundary error = |called_boundary - true_boundary|

Usage:
    python 04_eval_boundaries.py \
        --truth sim_refs/truth.tsv \
        --refs sim_refs/references.fasta \
        --itsxrust-dir tool_outputs/itsxrust \
        --itsx-dir tool_outputs/itsx \
        --itsxpress-dir tool_outputs/itsxpress \
        --outdir results
"""

import argparse
import csv
import os
import subprocess
import sys
import tempfile
from collections import defaultdict

import numpy as np


def load_truth(path):
    """Load ground truth TSV. Returns dict: ref_id -> {its1_start, its1_end, ...}"""
    truth = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            rid = row['ref_id']
            truth[rid] = {
                'its1_start': int(row['its1_start']),
                'its1_end': int(row['its1_end']),
                'its2_start': int(row['its2_start']),
                'its2_end': int(row['its2_end']),
                's58_start': int(row['s58_start']),
                's58_end': int(row['s58_end']),
                'ssu_end': int(row['ssu_end']),
                'lsu_start': int(row['lsu_start']),
            }
    return truth


def align_to_refs(query_fasta, ref_fasta, threads=4):
    """
    Align query sequences to references with minimap2.
    Returns list of (query_id, ref_id, ref_start_1based, ref_end_1based, strand).
    Uses primary alignments only.
    """
    alignments = []
    cmd = [
        'minimap2', '-a',
        '--secondary=no',  # primary only
        '-t', str(threads),
        ref_fasta, query_fasta
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"minimap2 warning: {result.stderr[:500]}", file=sys.stderr)

    for line in result.stdout.split('\n'):
        if line.startswith('@') or not line.strip():
            continue
        fields = line.split('\t')
        if len(fields) < 11:
            continue

        qname = fields[0]
        flag = int(fields[1])
        rname = fields[2]
        pos = int(fields[3])  # 1-based leftmost position
        cigar = fields[5]

        if rname == '*' or flag & 4:  # unmapped
            continue
        if flag & 256 or flag & 2048:  # secondary or supplementary
            continue

        strand = '-' if (flag & 16) else '+'

        # Calculate alignment end from CIGAR
        ref_consumed = 0
        num = ''
        for ch in cigar:
            if ch.isdigit():
                num += ch
            else:
                n = int(num) if num else 0
                num = ''
                if ch in ('M', 'D', 'N', '=', 'X'):
                    ref_consumed += n

        ref_start = pos           # 1-based inclusive
        ref_end = pos + ref_consumed - 1  # 1-based inclusive

        alignments.append((qname, rname, ref_start, ref_end, strand))

    return alignments


def fasta_to_tmp(sequences_dict):
    """Write dict of {id: seq} to a temp FASTA file. Returns path."""
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
    for sid, seq in sequences_dict.items():
        tmp.write(f">{sid}\n{seq}\n")
    tmp.close()
    return tmp.name


def read_fasta(path):
    """Read FASTA file. Returns dict: id -> seq."""
    seqs = {}
    cur_id = None
    cur_seq = []
    if not os.path.exists(path):
        return seqs
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if cur_id is not None:
                    seqs[cur_id] = ''.join(cur_seq)
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line)
    if cur_id is not None:
        seqs[cur_id] = ''.join(cur_seq)
    return seqs


def read_fastq(path):
    """Read FASTQ file. Returns dict: id -> seq."""
    seqs = {}
    if not os.path.exists(path):
        return seqs
    with open(path) as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().rstrip()
            f.readline()  # +
            f.readline()  # qual
            rid = header[1:].split()[0]
            seqs[rid] = seq
    return seqs


def trace_sim_read_to_ref(sim_read_id):
    """
    Extract the source reference ID from a simulated read's extracted header.
    ITSxRust headers: "sim_read_id|region:start-end"
    ITSx headers vary. We strip the |region:... suffix if present.
    """
    # Handle ITSxRust-style: "SRR21494940.123 0-862 +|full:31-802"
    base = sim_read_id.split('|')[0].strip()
    # Badread headers: "ref_name start-end strand" → first token is ref-derived
    # The sim read ID from Badread often has the ref name embedded
    # We'll rely on minimap2 alignment instead of parsing names
    return base


def collect_itsxrust_extractions(itsxrust_dir):
    """
    Collect ITSxRust extracted sequences by region.
    Returns dict: {region: {extracted_id: seq}}
    """
    regions = {}
    for region in ['full', 'its1', 'its2']:
        path = os.path.join(itsxrust_dir, f'out.{region}.fasta')
        if os.path.exists(path):
            regions[region] = read_fasta(path)
        else:
            regions[region] = {}
    return regions


def collect_itsx_extractions(itsx_dir):
    """
    Collect ITSx extracted sequences by region.
    ITSx outputs: itsx_out_ITS1.fasta, itsx_out_ITS2.fasta, itsx_out.full.fasta
    """
    regions = {}
    for region, pattern in [('its1', 'ITS1'), ('its2', 'ITS2'), ('full', 'full')]:
        path = os.path.join(itsx_dir, f'itsx_out_{pattern}.fasta')
        if not os.path.exists(path):
            path = os.path.join(itsx_dir, f'itsx_out.{pattern}.fasta')
        if os.path.exists(path):
            regions[region] = read_fasta(path)
        else:
            regions[region] = {}
    return regions


def collect_itsxpress_extractions(itsxpress_dir):
    """
    Collect ITSxpress extracted sequences.
    Files: its1.fastq, its2.fastq, all.fastq
    """
    regions = {}
    for region, fname in [('its1', 'its1.fastq'), ('its2', 'its2.fastq'), ('full', 'all.fastq')]:
        path = os.path.join(itsxpress_dir, fname)
        if os.path.exists(path):
            regions[region] = read_fastq(path)
        else:
            regions[region] = {}
    return regions


def eval_tool(tool_name, extractions, ref_fasta, truth, threads=4):
    """
    Evaluate boundary accuracy for one tool.
    For each region (its1, its2, full), align extracted sequences to references,
    then compare alignment boundaries to truth.

    Returns list of dicts with per-read boundary errors.
    """
    results = []

    for region in ['its1', 'its2', 'full']:
        seqs = extractions.get(region, {})
        if not seqs:
            print(f"  {tool_name} {region}: no extractions found, skipping")
            continue

        # Write to temp FASTA for alignment
        tmp_fa = fasta_to_tmp(seqs)

        try:
            alns = align_to_refs(tmp_fa, ref_fasta, threads=threads)
        finally:
            os.unlink(tmp_fa)

        n_aligned = 0
        for qname, rname, ref_start, ref_end, strand in alns:
            if rname not in truth:
                continue
            t = truth[rname]
            n_aligned += 1

            # Determine true boundaries for this region
            if region == 'its1':
                true_start = t['its1_start']
                true_end = t['its1_end']
            elif region == 'its2':
                true_start = t['its2_start']
                true_end = t['its2_end']
            elif region == 'full':
                true_start = t['its1_start']  # full = SSU_end+1 to LSU_start-1
                true_end = t['its2_end']       # ... which is its1_start to its2_end
                # Actually full ITS = from after SSU to before LSU
                true_start = t['ssu_end'] + 1  # = its1_start
                true_end = t['lsu_start'] - 1   # = its2_end typically

            start_err = ref_start - true_start
            end_err = ref_end - true_end

            results.append({
                'tool': tool_name,
                'region': region,
                'query_id': qname,
                'ref_id': rname,
                'called_start': ref_start,
                'called_end': ref_end,
                'true_start': true_start,
                'true_end': true_end,
                'start_error': start_err,
                'end_error': end_err,
                'abs_start_error': abs(start_err),
                'abs_end_error': abs(end_err),
            })

        print(f"  {tool_name} {region}: {len(seqs)} extracted, {n_aligned} aligned to refs")

    return results


def compute_stats(results, tool, region):
    """Compute summary statistics for one tool + region."""
    subset = [r for r in results if r['tool'] == tool and r['region'] == region]
    if not subset:
        return None

    abs_errs = [max(r['abs_start_error'], r['abs_end_error']) for r in subset]
    abs_errs = np.array(abs_errs)

    start_errs = np.array([r['abs_start_error'] for r in subset])
    end_errs = np.array([r['abs_end_error'] for r in subset])

    return {
        'tool': tool,
        'region': region,
        'n': len(subset),
        'median_abs_err': float(np.median(abs_errs)),
        'mean_abs_err': float(np.mean(abs_errs)),
        'pct95_abs_err': float(np.percentile(abs_errs, 95)),
        'pct_within_5bp': float(np.mean(abs_errs <= 5) * 100),
        'pct_within_10bp': float(np.mean(abs_errs <= 10) * 100),
        'median_start_err': float(np.median(start_errs)),
        'median_end_err': float(np.median(end_errs)),
    }


def main():
    parser = argparse.ArgumentParser(description='Evaluate boundary accuracy on simulated reads')
    parser.add_argument('--truth', required=True, help='Ground truth TSV')
    parser.add_argument('--refs', required=True, help='Reference FASTA')
    parser.add_argument('--itsxrust-dir', required=True, help='ITSxRust output directory')
    parser.add_argument('--itsx-dir', required=True, help='ITSx output directory')
    parser.add_argument('--itsxpress-dir', required=True, help='ITSxpress output directory')
    parser.add_argument('--outdir', default='results', help='Output directory')
    parser.add_argument('--threads', type=int, default=4, help='minimap2 threads')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Load truth
    print("Loading ground truth...")
    truth = load_truth(args.truth)
    print(f"  {len(truth)} reference sequences with known boundaries")

    # Collect extractions from each tool
    print("\nCollecting extracted sequences...")
    tools = {
        'ITSxRust': collect_itsxrust_extractions(args.itsxrust_dir),
        'ITSx': collect_itsx_extractions(args.itsx_dir),
        'ITSxpress': collect_itsxpress_extractions(args.itsxpress_dir),
    }

    # Evaluate each tool
    all_results = []
    for tool_name, extractions in tools.items():
        print(f"\nEvaluating {tool_name}...")
        results = eval_tool(tool_name, extractions, args.refs, truth, threads=args.threads)
        all_results.extend(results)

    # Write per-read errors
    errors_path = os.path.join(args.outdir, 'boundary_errors.tsv')
    with open(errors_path, 'w') as f:
        if all_results:
            writer = csv.DictWriter(f, fieldnames=all_results[0].keys(), delimiter='\t')
            writer.writeheader()
            writer.writerows(all_results)
    print(f"\nPer-read errors: {errors_path} ({len(all_results)} rows)")

    # Compute summary statistics
    print("\n" + "=" * 70)
    print("BOUNDARY ACCURACY SUMMARY (Table 4)")
    print("=" * 70)

    header = f"{'Tool':<15} {'Region':<8} {'N':>6} {'Med err':>8} {'Mean':>8} {'95th':>8} {'≤5bp':>8} {'≤10bp':>8}"
    print(header)
    print("-" * len(header))

    summary_rows = []
    for tool in ['ITSxRust', 'ITSx', 'ITSxpress']:
        for region in ['its1', 'its2', 'full']:
            stats = compute_stats(all_results, tool, region)
            if stats is None:
                print(f"{tool:<15} {region:<8} {'---':>6}")
                continue
            print(f"{tool:<15} {region:<8} {stats['n']:>6} "
                  f"{stats['median_abs_err']:>8.1f} {stats['mean_abs_err']:>8.1f} "
                  f"{stats['pct95_abs_err']:>8.1f} {stats['pct_within_5bp']:>7.1f}% "
                  f"{stats['pct_within_10bp']:>7.1f}%")
            summary_rows.append(stats)

    # Write summary TSV
    summary_path = os.path.join(args.outdir, 'table4_boundary_accuracy.tsv')
    with open(summary_path, 'w') as f:
        if summary_rows:
            writer = csv.DictWriter(f, fieldnames=summary_rows[0].keys(), delimiter='\t')
            writer.writeheader()
            writer.writerows(summary_rows)
    print(f"\nSummary table: {summary_path}")


if __name__ == '__main__':
    main()
