#!/usr/bin/env python3
"""
Figure 3: Boundary error distributions on simulated ONT reads.
Reads boundary_errors.tsv from the simulation pipeline.

Usage:
    python fig3_boundary_errors.py --input /path/to/results/boundary_errors.tsv
"""
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np

def load_data(path):
    """Parse boundary_errors.tsv -> list of dicts."""
    rows = []
    with open(path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < len(header):
                continue
            d = dict(zip(header, fields))
            try:
                d['max_err'] = int(d['max_err'])
            except (KeyError, ValueError):
                # fallback: try abs_start_error + abs_end_error -> max
                try:
                    se = int(d.get('abs_start_error', d.get('start_error', '0')))
                    ee = int(d.get('abs_end_error', d.get('end_error', '0')))
                    d['max_err'] = max(abs(se), abs(ee))
                except ValueError:
                    continue
            rows.append(d)
    return rows

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='boundary_errors.tsv path')
    parser.add_argument('--outfile', default='fig3_boundary_errors')
    args = parser.parse_args()

    rows = load_data(args.input)
    print(f"Loaded {len(rows)} rows from {args.input}")

    tools = ['ITSxRust', 'ITSx', 'ITSxpress']
    tool_colors = {'ITSxRust': '#4C78A8', 'ITSx': '#F58518', 'ITSxpress': '#54A24B'}
    regions = ['its1', 'its2', 'full']
    region_labels = {'its1': 'ITS1', 'its2': 'ITS2', 'full': 'Full ITS'}

    # ── Collect data per (tool, region) ──
    data = {}
    for t in tools:
        for reg in regions:
            key = (t, reg)
            # tool name in TSV may be lowercase
            errs = [r['max_err'] for r in rows
                    if r.get('tool', '').lower().startswith(t.lower())
                    and r.get('region', '').lower() == reg]
            data[key] = np.array(errs) if errs else np.array([])

    # ── Figure: 2 rows × 3 cols (top: violin, bottom: ECDF) ──
    fig, axes = plt.subplots(2, 3, figsize=(12, 6.5),
                              gridspec_kw={'height_ratios': [1, 1]})
    fig.subplots_adjust(hspace=0.35, wspace=0.28, left=0.07, right=0.97,
                        top=0.88, bottom=0.08)

    # ── Row 1: Box plots (capped at 100bp for visibility) ──
    CAP = 100
    for j, reg in enumerate(regions):
        ax = axes[0][j]
        box_data = []
        labels = []
        colors = []
        for t in tools:
            arr = data[(t, reg)]
            if len(arr) == 0:
                continue
            capped = np.clip(arr, 0, CAP)
            box_data.append(capped)
            labels.append(f'{t}\n(n={len(arr):,})')
            colors.append(tool_colors[t])

        if box_data:
            bp = ax.boxplot(box_data, tick_labels=labels, patch_artist=True,
                           showfliers=False, widths=0.55,
                           medianprops={'color': 'black', 'linewidth': 1.5})
            for patch, color in zip(bp['boxes'], colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.65)

        ax.set_ylabel('Max boundary error (bp)' if j == 0 else '', fontsize=9)
        ax.set_ylim(-2, CAP + 5)
        ax.set_title(region_labels[reg], fontsize=10, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis='x', labelsize=7.5)
        # reference lines
        ax.axhline(5, color='#999', linewidth=0.5, linestyle=':', zorder=0)
        ax.axhline(10, color='#bbb', linewidth=0.5, linestyle=':', zorder=0)
        if j == 2:
            ax.text(3.6, 6, '5 bp', fontsize=6.5, color='#777')
            ax.text(3.6, 11, '10 bp', fontsize=6.5, color='#999')

    # ── Row 2: ECDF up to 50 bp ──
    ECDF_MAX = min(500, max(arr))
    for j, reg in enumerate(regions):
        ax = axes[1][j]
        for t in tools:
            arr = data[(t, reg)]
            if len(arr) == 0:
                continue
            sorted_e = np.sort(arr)
            ecdf_y = np.arange(1, len(sorted_e) + 1) / len(sorted_e) * 100
            # plot only up to ECDF_MAX
            mask = sorted_e <= ECDF_MAX
            if mask.any():
                ax.plot(sorted_e[mask], ecdf_y[mask], label=f'{t} (n={len(arr):,})',
                        color=tool_colors[t], linewidth=1.5)
                # extend to edge
                last_y = ecdf_y[mask][-1]
                ax.plot([sorted_e[mask][-1], ECDF_MAX], [last_y, last_y],
                        color=tool_colors[t], linewidth=1.5, alpha=0.3)

        ax.set_xlim(0, ECDF_MAX)
        ax.set_ylim(0, 105)
        ax.set_xlabel('Max boundary error (bp)', fontsize=9)
        ax.set_ylabel('Cumulative % of reads' if j == 0 else '', fontsize=9)
        ax.yaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        # reference lines
        ax.axvline(5, color='#ccc', linewidth=0.5, linestyle=':', zorder=0)
        ax.axvline(10, color='#ddd', linewidth=0.5, linestyle=':', zorder=0)
        if j == 0:
            ax.legend(fontsize=7, loc='lower right', framealpha=0.9)

    fig.suptitle('Boundary accuracy on simulated ONT reads', fontsize=14,
                 fontweight='bold', y=0.98)

    plt.savefig(f'{args.outfile}.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{args.outfile}.png', dpi=300, bbox_inches='tight')
    print(f"Saved {args.outfile}.pdf / .png")


if __name__ == '__main__':
    main()
