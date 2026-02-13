#!/usr/bin/env python3
"""
Figure 1: Throughput scaling with thread count.
Reads thread_scaling.tsv produced by fig1_bench_threads.sh.

Usage:
    python fig1_thread_scaling.py --input thread_scaling.tsv
"""
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='thread_scaling.tsv')
    parser.add_argument('--reads', type=int, default=54659,
                        help='Total reads in dataset (for reads/sec calculation)')
    parser.add_argument('--outfile', default='fig1_thread_scaling')
    args = parser.parse_args()

    # ── Parse TSV ──
    # tool  threads  rep  wall_sec  peak_rss_kb
    raw = defaultdict(lambda: defaultdict(list))  # raw[tool][threads] = [wall_sec, ...]
    rss = defaultdict(lambda: defaultdict(list))

    with open(args.input) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            tool = parts[0]
            threads = int(parts[1])
            wall = float(parts[3])
            peak_rss = int(parts[4])
            raw[tool][threads].append(wall)
            rss[tool][threads].append(peak_rss)

    tools = sorted(raw.keys())
    tool_colors = {'itsxrust': '#4C78A8', 'itsx': '#F58518', 'itsxpress': '#54A24B'}
    tool_labels = {'itsxrust': 'ITSxRust', 'itsx': 'ITSx', 'itsxpress': 'ITSxpress v2'}
    tool_markers = {'itsxrust': 'o', 'itsx': 's', 'itsxpress': '^'}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.2))
    fig.subplots_adjust(wspace=0.32, left=0.08, right=0.96, top=0.88, bottom=0.14)

    # ── Panel A: Throughput (reads/sec) vs threads ──
    for tool in tools:
        threads_sorted = sorted(raw[tool].keys())
        means = []
        stds = []
        for t in threads_sorted:
            times = raw[tool][t]
            rates = [args.reads / w for w in times]
            means.append(np.mean(rates))
            stds.append(np.std(rates))

        c = tool_colors.get(tool, '#999')
        m = tool_markers.get(tool, 'o')
        lbl = tool_labels.get(tool, tool)

        ax1.errorbar(threads_sorted, means, yerr=stds, marker=m, color=c,
                     label=lbl, capsize=3, linewidth=1.5, markersize=6)

    ax1.set_xlabel('Threads', fontsize=10)
    ax1.set_ylabel('Throughput (reads/sec)', fontsize=10)
    ax1.set_title('A) Throughput scaling', fontsize=11, fontweight='bold', loc='left')
    ax1.legend(fontsize=8.5)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    if raw:
        all_threads = sorted(set(t for tool in raw for t in raw[tool]))
        ax1.set_xticks(all_threads)

    # ── Panel B: Wall time vs threads ──
    for tool in tools:
        threads_sorted = sorted(raw[tool].keys())
        means = [np.mean(raw[tool][t]) / 60 for t in threads_sorted]  # minutes
        stds = [np.std(raw[tool][t]) / 60 for t in threads_sorted]

        c = tool_colors.get(tool, '#999')
        m = tool_markers.get(tool, 'o')
        lbl = tool_labels.get(tool, tool)

        ax2.errorbar(threads_sorted, means, yerr=stds, marker=m, color=c,
                     label=lbl, capsize=3, linewidth=1.5, markersize=6)

    ax2.set_xlabel('Threads', fontsize=10)
    ax2.set_ylabel('Wall time (minutes)', fontsize=10)
    ax2.set_title('B) Wall time', fontsize=11, fontweight='bold', loc='left')
    ax2.legend(fontsize=8.5)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    if raw:
        ax2.set_xticks(all_threads)

    fig.suptitle(f'Performance scaling ({args.reads:,} ONT reads)',
                 fontsize=12, fontweight='bold', y=0.98)

    plt.savefig(f'{args.outfile}.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{args.outfile}.png', dpi=300, bbox_inches='tight')
    print(f"Saved {args.outfile}.pdf / .png")


if __name__ == '__main__':
    main()
