#!/usr/bin/env python3
"""
06_tabulate.py

Read results TSVs and generate LaTeX table fragments for the paper.

Usage:
    python 06_tabulate.py <results_dir>
"""

import csv
import os
import sys


def load_tsv(path):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))


def gen_table4(rows, outpath):
    """Generate Table 4 (boundary accuracy) LaTeX."""
    # Group by tool, aggregate across regions (or show per-region)
    # Paper format: Simulation | Tool | Med abs err | Mean | 95th | % ≤5bp

    # We only have ONT-like simulation for now
    with open(outpath, 'w') as f:
        f.write(r"""% Table 4: Boundary accuracy on simulated ONT reads
\begin{table}[ht]
\centering
\caption{Boundary accuracy on simulated ONT reads with known true ITS boundaries.
Maximum of start and end absolute errors is reported per read.}
\label{tab:boundary}
\small
\begin{tabular}{@{}llrrrr@{}}
\toprule
Tool & Region & Med.\ abs.\ err.\ (bp) & Mean & 95th \%ile & \% $\leq$5\,bp \\
\midrule
""")

        tool_order = ['ITSxRust', 'ITSx', 'ITSxpress']
        region_order = ['its1', 'its2', 'full']
        region_labels = {'its1': 'ITS1', 'its2': 'ITS2', 'full': 'Full ITS'}

        prev_tool = None
        for tool in tool_order:
            for region in region_order:
                matching = [r for r in rows
                            if r['tool'] == tool and r['region'] == region]
                if not matching:
                    continue
                r = matching[0]

                # Add spacing between tools
                if prev_tool is not None and tool != prev_tool:
                    f.write(r'\addlinespace' + '\n')
                prev_tool = tool

                med = float(r['median_abs_err'])
                mean = float(r['mean_abs_err'])
                p95 = float(r['pct95_abs_err'])
                within5 = float(r['pct_within_5bp'])

                f.write(f"{tool} & {region_labels[region]} & "
                        f"{med:.1f} & {mean:.1f} & {p95:.1f} & {within5:.1f}\\% \\\\\n")

        f.write(r"""\bottomrule
\end{tabular}
\end{table}
""")

    print(f"  Table 4 LaTeX: {outpath}")


def gen_table5(rows, outpath):
    """Generate Table 5 (classification accuracy) LaTeX."""
    with open(outpath, 'w') as f:
        f.write(r"""% Table 5: Classification accuracy using UNITE + VSEARCH
\begin{table}[ht]
\centering
\caption{Classification accuracy on simulated ONT reads using UNITE and VSEARCH.
Accuracy is measured as agreement with the classification of the error-free reference sequence.}
\label{tab:classification}
\small
\begin{tabular}{@{}lrrrl@{}}
\toprule
Input to classifier & $N$ & Genus acc.\ & Species acc.\ & Notes \\
\midrule
""")

        for r in rows:
            label = r['label']
            n = int(r['n_classified'])
            g_acc = float(r['genus_accuracy'])
            s_acc = float(r['species_accuracy'])

            notes = ""
            if 'Raw' in label:
                notes = "Conserved flanks included"

            f.write(f"{label} & {n:,} & {g_acc:.1f}\\% & {s_acc:.1f}\\% & {notes} \\\\\n")

        f.write(r"""\bottomrule
\end{tabular}
\end{table}
""")

    print(f"  Table 5 LaTeX: {outpath}")


def main():
    results_dir = sys.argv[1] if len(sys.argv) > 1 else 'results'

    print("Generating LaTeX tables...")

    # Table 4
    t4_rows = load_tsv(os.path.join(results_dir, 'table4_boundary_accuracy.tsv'))
    if t4_rows:
        gen_table4(t4_rows, os.path.join(results_dir, 'table4.tex'))
    else:
        print("  Table 4: no data found")

    # Table 5
    t5_rows = load_tsv(os.path.join(results_dir, 'table5_classification.tsv'))
    if t5_rows:
        gen_table5(t5_rows, os.path.join(results_dir, 'table5.tex'))
    else:
        print("  Table 5: no data found")

    print("Done.")


if __name__ == '__main__':
    main()
