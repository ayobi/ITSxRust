#!/usr/bin/env python3
"""Plot boundary error distributions (boxplots + cumulative curves) in the same
visual style as fig3_boundary_errors.py.

Input: boundary_errors.tsv from the simulation evaluation.
It tolerates different column names by computing max absolute error if needed.
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick


def _to_float(x):
    try:
        if x is None:
            return None
        x = str(x).strip()
        if x == "" or x.lower() in {"na", "nan", "none"}:
            return None
        return float(x)
    except Exception:
        return None


def _norm_tool(t):
    t = (t or "").strip()
    tl = t.lower()
    if tl.startswith("itsxrust"):
        return "ITSxRust"
    if tl.startswith("itsxpress"):
        return "ITSxpress"
    if tl.startswith("itsx"):
        return "ITSx"
    # fall back to original
    return t or "unknown"


def _norm_region(r):
    r = (r or "").strip().lower()
    if r in {"its1", "its_1", "its-1", "its1_region", "region_its1", "its1"}:
        return "its1"
    if r in {"its2", "its_2", "its-2", "its2_region", "region_its2", "its2"}:
        return "its2"
    if r in {"full", "fullits", "full_its", "full its", "all", "its", "full_its_region"}:
        return "full"
    return r


def _extract_max_err(d):
    # Direct max columns (common possibilities)
    for k in [
        "max_err", "max_error", "max_abs_err", "max_abs_error",
        "max_abs_err_bp", "max_abs_error_bp",
    ]:
        if k in d:
            v = _to_float(d.get(k))
            if v is not None:
                return abs(v)

    # Start/end columns: take max(|start|, |end|)
    start_keys = [
        "abs_start_error", "start_abs_error", "start_abs_err", "start_abs_err_bp",
        "start_err", "start_error", "start_error_bp", "start_delta", "start_diff",
    ]
    end_keys = [
        "abs_end_error", "end_abs_error", "end_abs_err", "end_abs_err_bp",
        "end_err", "end_error", "end_error_bp", "end_delta", "end_diff",
    ]
    s = None
    for k in start_keys:
        if k in d:
            s = _to_float(d.get(k))
            if s is not None:
                break
    e = None
    for k in end_keys:
        if k in d:
            e = _to_float(d.get(k))
            if e is not None:
                break

    if s is None and e is None:
        return None
    if s is None:
        return abs(e)
    if e is None:
        return abs(s)
    return max(abs(s), abs(e))


def load_data(path: Path):
    data = defaultdict(list)

    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Input looks empty or has no header row")
        # normalize keys to lowercase for robust lookup
        for row in reader:
            d = {k.strip().lower(): (v.strip() if isinstance(v, str) else v)
                 for k, v in row.items() if k is not None}
            tool = _norm_tool(d.get("tool"))
            region = _norm_region(d.get("region"))
            if tool == "unknown" or region not in {"its1", "its2", "full"}:
                continue

            max_err = _extract_max_err(d)
            if max_err is None:
                continue
            data[(tool, region)].append(max_err)

    # convert to numpy arrays
    out = {}
    for k, vals in data.items():
        out[k] = np.asarray(vals, dtype=float)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="boundary_errors.tsv")
    ap.add_argument("--outfile", required=True, help="output prefix (no extension)")
    ap.add_argument("--title", default="Boundary error distributions", help="figure title")
    ap.add_argument("--cap", type=float, default=100, help="cap for boxplots (bp)")
    ap.add_argument("--ecdf-max", type=float, default=None, help="x-axis max for cumulative curves (bp)")
    args = ap.parse_args()

    in_path = Path(args.input)
    out_prefix = args.outfile

    data = load_data(in_path)

    tools = ["ITSxRust", "ITSx", "ITSxpress"]
    regions = ["its1", "its2", "full"]
    region_labels = {"its1": "ITS1", "its2": "ITS2", "full": "Full ITS"}

    # Color palette to match the ONT figure style
    tool_colors = {
        "ITSxRust": "#4C78A8",  # blue
        "ITSx": "#F58518",      # orange
        "ITSxpress": "#54A24B", # green
    }

    # Determine ECDF max if not provided
    if args.ecdf_max is None:
        all_vals = []
        for k in data:
            all_vals.extend(data[k].tolist())
        ecdf_max = 500
        if all_vals:
            ecdf_max = min(500, float(np.max(all_vals)))
    else:
        ecdf_max = float(args.ecdf_max)

    CAP = float(args.cap)

    # ── Figure: 2 rows × 3 cols (top: boxplot, bottom: cumulative %) ──
    fig, axes = plt.subplots(
        2, 3, figsize=(12, 6.5),
        gridspec_kw={"height_ratios": [1, 1]}
    )
    fig.subplots_adjust(
        left=0.06, right=0.97, top=0.88, bottom=0.09, wspace=0.20, hspace=0.28
    )

    # ── Row 1: Box plots (capped for visibility) ──
    for j, reg in enumerate(regions):
        ax = axes[0][j]
        box_data = []
        labels = []
        colors = []
        for t in tools:
            arr = data.get((t, reg), np.array([]))
            if len(arr) == 0:
                continue
            capped = np.clip(arr, 0, CAP)
            box_data.append(capped)
            labels.append(f"{t}\n(n={len(arr):,})")
            colors.append(tool_colors[t])

        if box_data:
            bp = ax.boxplot(
                box_data, tick_labels=labels, patch_artist=True,
                showfliers=False, widths=0.55,
                medianprops={"color": "black", "linewidth": 1.5}
            )
            for patch, color in zip(bp["boxes"], colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.65)

        ax.set_ylabel("Max boundary error (bp)" if j == 0 else "", fontsize=9)
        ax.set_ylim(-2, CAP + 5)
        ax.set_title(region_labels[reg], fontsize=10, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="x", labelsize=7.5)

        # keep both lines
        ax.axhline(5,  ls=":", lw=1, color="0.7")
        ax.axhline(10, ls=":", lw=1, color="0.7")

        # label only 10 bp
        ax.text(1.02, 10, "10 bp", transform=ax.get_yaxis_transform(),
                ha="left", va="center", fontsize=9, color="0.5", clip_on=False)

           
    # ── Row 2: Cumulative % curves up to ECDF_MAX ──
    for j, reg in enumerate(regions):
        ax = axes[1][j]

        for t in tools:
            arr = data.get((t, reg), np.array([]))
            if len(arr) == 0:
                continue
            sorted_e = np.sort(arr)
            y = np.arange(1, len(sorted_e) + 1) / len(sorted_e) * 100.0
            mask = sorted_e <= ecdf_max
            if mask.any():
                ax.plot(
                    sorted_e[mask], y[mask],
                    label=f"{t} (n={len(arr):,})",
                    color=tool_colors[t], linewidth=1.5
                )
                # extend horizontally to the right edge (like ONT figure)
                last_y = y[mask][-1]
                ax.plot(
                    [sorted_e[mask][-1], ecdf_max], [last_y, last_y],
                    color=tool_colors[t], linewidth=1.5, alpha=0.3
                )

        ax.set_xlim(0, ecdf_max)
        ax.set_ylim(0, 105)
        ax.set_xlabel("Max boundary error (bp)", fontsize=9)
        ax.set_ylabel("Cumulative % of reads" if j == 0 else "", fontsize=9)
        ax.yaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # reference lines
        ax.axvline(5, color="#ccc", linewidth=0.5, linestyle=":", zorder=0)
        ax.axvline(10, color="#ddd", linewidth=0.5, linestyle=":", zorder=0)

        if j == 0:
            ax.legend(fontsize=7, loc="lower right", framealpha=0.9)

    fig.suptitle(args.title, fontsize=14, fontweight="bold", y=0.985)

    plt.savefig(f"{out_prefix}.pdf", dpi=300, bbox_inches="tight", pad_inches=0.12)
    plt.savefig(f"{out_prefix}.png", dpi=300, bbox_inches="tight", pad_inches=0.12)
    print(f"Saved {out_prefix}.pdf / .png")


if __name__ == "__main__":
    main()
