#!/usr/bin/env python3
"""bench/collect_results.py

Parse benchmark outputs and produce summary TSVs for Tables 2 and 3.

Usage:
    python3 bench/collect_results.py <results_dir> <total_reads>
"""

import json
import os
import re
import sys
from pathlib import Path


def parse_wall_seconds(timefile: Path) -> float:
    """Parse wall time from a .time file (written by run_benchmark.sh)."""
    text = timefile.read_text()
    # Format from our script: "wall=<value> rss_kb=<value>"
    m = re.search(r"wall=([\d:.]+)", text)
    if not m:
        return float("nan")
    val = m.group(1)
    # GNU time: h:mm:ss or m:ss.ss  |  macOS: seconds as float
    if ":" in val:
        parts = val.split(":")
        if len(parts) == 3:
            return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
        elif len(parts) == 2:
            return int(parts[0]) * 60 + float(parts[1])
    return float(val)


def parse_rss_kb(timefile: Path) -> float:
    """Parse peak RSS (KB) from a .time file."""
    text = timefile.read_text()
    m = re.search(r"rss_kb=(\d+)", text)
    return int(m.group(1)) if m else float("nan")


def count_fasta_seqs(path: Path) -> int:
    """Count '>' headers in a FASTA file."""
    if not path.exists():
        return 0
    count = 0
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count


def count_fastq_seqs(path: Path) -> int:
    """Count records in a FASTQ file (every 4th line starting with @)."""
    if not path.exists():
        return 0
    count = 0
    with open(path) as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                count += 1
    return count


def load_qc_json(path: Path) -> dict:
    """Load an ITSxRust QC JSON."""
    if not path.exists():
        return {}
    with open(path) as f:
        return json.load(f)


def count_itsx_outputs(outdir: Path, prefix: str) -> dict:
    """Count sequences in ITSx output files."""
    counts = {}
    for region, suffix in [("ITS1", ".ITS1.fasta"), ("ITS2", ".ITS2.fasta"),
                           ("full", ".full.fasta"), ("5_8S", ".5_8S.fasta")]:
        p = outdir / f"{prefix}{suffix}"
        counts[region] = count_fasta_seqs(p)
    return counts


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <results_dir> <total_reads>")
        sys.exit(1)

    outdir = Path(sys.argv[1])
    total_reads = int(sys.argv[2])

    # -----------------------------------------------------------------------
    # Table 2: Throughput and resource use
    # -----------------------------------------------------------------------
    table2_rows = []

    tools = [
        ("ITSxRust (full)", "itsxrust_full"),
        ("ITSxRust (all)", "itsxrust_all"),
        ("ITSxRust (derep)", "itsxrust_derep"),
        ("ITSxRust (ONT preset)", "itsxrust_ont"),
        ("ITSx", "itsx"),
        ("ITSxpress (ITS1)", "itsxpress_its1"),
        ("ITSxpress (ITS2)", "itsxpress_its2"),
        ("ITSxpress (ALL)", "itsxpress_all"),
    ]

    for label, prefix in tools:
        timefile = outdir / f"{prefix}.time"
        if not timefile.exists():
            continue

        wall_sec = parse_wall_seconds(timefile)
        rss_kb = parse_rss_kb(timefile)
        reads_per_sec = total_reads / wall_sec if wall_sec > 0 else 0
        rss_gb = rss_kb / (1024 * 1024)

        table2_rows.append({
            "tool": label,
            "total_reads": total_reads,
            "wall_sec": f"{wall_sec:.1f}",
            "wall_min": f"{wall_sec / 60:.1f}",
            "reads_per_sec": f"{reads_per_sec:.0f}",
            "peak_rss_gb": f"{rss_gb:.2f}",
        })

    # Write Table 2
    t2_path = outdir / "table2_throughput.tsv"
    with open(t2_path, "w") as f:
        header = ["tool", "total_reads", "wall_sec", "wall_min", "reads_per_sec", "peak_rss_gb"]
        f.write("\t".join(header) + "\n")
        for row in table2_rows:
            f.write("\t".join(str(row[h]) for h in header) + "\n")
    print(f"Wrote: {t2_path}")

    # -----------------------------------------------------------------------
    # Table 3: Extraction success rates
    # -----------------------------------------------------------------------
    table3_rows = []

    # ITSxRust (from QC JSONs)
    for label, prefix in [("ITSxRust (default)", "itsxrust_all"),
                          ("ITSxRust (ONT)", "itsxrust_ont")]:
        qc = load_qc_json(outdir / f"{prefix}_qc.json")
        if not qc:
            continue

        kept = qc.get("kept", {})
        regions = qc.get("regions", {})

        n_full = regions.get("full", kept.get("total", 0))
        n_its1 = regions.get("its1", 0)
        n_its2 = regions.get("its2", 0)
        n_ambig = kept.get("ambiguous", 0)

        table3_rows.append({
            "tool": label,
            "pct_its1": f"{100 * n_its1 / total_reads:.1f}" if total_reads else "0",
            "pct_its2": f"{100 * n_its2 / total_reads:.1f}" if total_reads else "0",
            "pct_full": f"{100 * n_full / total_reads:.1f}" if total_reads else "0",
            "pct_ambig": f"{100 * n_ambig / total_reads:.1f}" if total_reads else "—",
            "n_its1": n_its1,
            "n_its2": n_its2,
            "n_full": n_full,
            "n_ambig": n_ambig,
        })

    # ITSx
    itsx_counts = count_itsx_outputs(outdir, "itsx_out")
    if any(v > 0 for v in itsx_counts.values()):
        table3_rows.append({
            "tool": "ITSx",
            "pct_its1": f"{100 * itsx_counts['ITS1'] / total_reads:.1f}",
            "pct_its2": f"{100 * itsx_counts['ITS2'] / total_reads:.1f}",
            "pct_full": f"{100 * itsx_counts['full'] / total_reads:.1f}",
            "pct_ambig": "—",
            "n_its1": itsx_counts["ITS1"],
            "n_its2": itsx_counts["ITS2"],
            "n_full": itsx_counts["full"],
            "n_ambig": "—",
        })

    # ITSxpress
    for label, prefix, region in [
        ("ITSxpress (ITS1)", "itsxpress_its1", "ITS1"),
        ("ITSxpress (ITS2)", "itsxpress_its2", "ITS2"),
        ("ITSxpress (ALL)", "itsxpress_all", "full"),
    ]:
        fq_path = outdir / f"{prefix}.fastq"
        n = count_fastq_seqs(fq_path)
        if n > 0:
            row = {
                "tool": label,
                "pct_its1": "", "pct_its2": "", "pct_full": "",
                "pct_ambig": "—",
                "n_its1": "", "n_its2": "", "n_full": "",
                "n_ambig": "—",
            }
            pct = f"{100 * n / total_reads:.1f}"
            if region == "ITS1":
                row["pct_its1"] = pct
                row["n_its1"] = n
            elif region == "ITS2":
                row["pct_its2"] = pct
                row["n_its2"] = n
            else:
                row["pct_full"] = pct
                row["n_full"] = n
            table3_rows.append(row)

    # Write Table 3
    t3_path = outdir / "table3_extraction.tsv"
    with open(t3_path, "w") as f:
        header = ["tool", "pct_its1", "n_its1", "pct_its2", "n_its2",
                  "pct_full", "n_full", "pct_ambig", "n_ambig"]
        f.write("\t".join(header) + "\n")
        for row in table3_rows:
            f.write("\t".join(str(row.get(h, "")) for h in header) + "\n")
    print(f"Wrote: {t3_path}")

    # -----------------------------------------------------------------------
    # Print summary
    # -----------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("  TABLE 2: Throughput and Resource Use")
    print("=" * 60)
    print(f"{'Tool':<30} {'Wall(s)':>8} {'Reads/s':>9} {'RSS(GB)':>8}")
    print("-" * 60)
    for r in table2_rows:
        print(f"{r['tool']:<30} {r['wall_sec']:>8} {r['reads_per_sec']:>9} {r['peak_rss_gb']:>8}")

    print("\n" + "=" * 60)
    print("  TABLE 3: Extraction Success Rates")
    print("=" * 60)
    print(f"{'Tool':<30} {'ITS1%':>7} {'ITS2%':>7} {'Full%':>7} {'Ambig%':>7}")
    print("-" * 60)
    for r in table3_rows:
        print(f"{r['tool']:<30} {r.get('pct_its1',''):>7} {r.get('pct_its2',''):>7} "
              f"{r.get('pct_full',''):>7} {r.get('pct_ambig','—'):>7}")


if __name__ == "__main__":
    main()