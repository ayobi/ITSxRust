#!/usr/bin/env python3
"""
Figure 4: ITSxRust pipeline overview schematic (pastel theme).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import matplotlib.patches as mpatches


# ── Pastel palette (soft blue / peach / mint + neutrals) ────────────────
C_IN   = "#E5E7EB"   # light gray (input)
C_HMM  = "#BFD7EA"   # pastel blue
C_PROC = "#FAD4B8"   # pastel peach
C_DEC  = "#FAD4B8"   # decision (same family)
C_OUT  = "#CBE7D3"   # pastel mint
C_QC   = "#D6D3F0"   # pastel lavender (distinct)
C_SKIP = "#F7B7B2"   # pastel salmon

TXT_DARK = "#1F2937"  # dark gray text
EDGE     = "#6B7280"  # softer border than #555
ARROW    = "#6B7280"  # arrows
LABEL    = "#6B7280"  # small italic labels


def rbox(
    ax, cx, cy, w, h, text, color,
    text_color=TXT_DARK, fontsize=10,
    subtext=None, subtextsize=8,
    alpha=1.0, edgecolor=EDGE, lw=1.0
):
    """Rounded rectangle with centered text."""
    box = FancyBboxPatch(
        (cx - w/2, cy - h/2), w, h,
        boxstyle="round,pad=0.015",
        linewidth=lw,
        edgecolor=edgecolor,
        facecolor=color,
        alpha=alpha,
        zorder=3
    )
    ax.add_patch(box)

    if subtext:
        ax.text(
            cx, cy + h * 0.14, text,
            ha="center", va="bottom",
            fontsize=fontsize, fontweight="bold",
            color=text_color, zorder=4
        )
        ax.text(
            cx, cy - h * 0.14, subtext,
            ha="center", va="top",
            fontsize=subtextsize,
            color=text_color, alpha=0.95, zorder=4
        )
    else:
        ax.text(
            cx, cy, text,
            ha="center", va="center",
            fontsize=fontsize, fontweight="bold",
            color=text_color, zorder=4
        )


def arr(ax, x1, y1, x2, y2, color=ARROW, lw=1.2):
    ax.annotate(
        "", xy=(x2, y2), xytext=(x1, y1),
        arrowprops=dict(
            arrowstyle="->",
            color=color, lw=lw,
            shrinkA=1, shrinkB=1
        ),
        zorder=2
    )


def main():
    fig, ax = plt.subplots(figsize=(10, 8.5))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # ── Layout ──────────────────────────────────────────────────────────
    XL = 0.33           # main column (left-of-center)
    XR = 0.77           # right column
    BW = 0.32           # main box width
    BWR = 0.22          # right box width
    BH = 0.065          # box height
    VGAP = 0.105        # vertical gap between centers
    Y0 = 0.93           # top

    # ════════════════════════════════════════════
    # 1. Input reads
    # ════════════════════════════════════════════
    Y1 = Y0
    rbox(ax, XL, Y1, BW, BH, "Input reads", C_IN,
         subtext="FASTQ / FASTA (.gz)")

    # ════════════════════════════════════════════
    # 2. nhmmer search
    # ════════════════════════════════════════════
    Y2 = Y1 - VGAP
    rbox(ax, XL, Y2, BW, BH, "nhmmer search", C_HMM,
         subtext="vs SSU / 5.8S / LSU profiles")
    arr(ax, XL, Y1 - BH/2, XL, Y2 + BH/2)

    # HMM profiles (right)
    rbox(ax, XR, Y2, BWR, BH * 0.8, "HMM profiles", "#DCECF7",
         text_color=TXT_DARK, fontsize=9)
    arr(ax, XR - BWR/2, Y2, XL + BW/2, Y2, color=ARROW)

    # ════════════════════════════════════════════
    # 3. Stream tblout
    # ════════════════════════════════════════════
    Y3 = Y2 - VGAP
    rbox(ax, XL, Y3, BW, BH, "Stream tblout", C_PROC,
         subtext="Top-K per anchor × strand")
    arr(ax, XL, Y2 - BH/2, XL, Y3 + BH/2)

    # Arrow label
    mid_23 = (Y2 - BH/2 + Y3 + BH/2) / 2
    ax.text(XL + BW/2 + 0.02, mid_23, "--tblout",
            fontsize=7.5, color=LABEL, fontstyle="italic", va="center")

    # ════════════════════════════════════════════
    # 4. Chain selection
    # ════════════════════════════════════════════
    Y4 = Y3 - VGAP
    rbox(ax, XL, Y4, BW, BH, "Chain selection", C_PROC,
         subtext="4-anchor or 2-anchor fallback")
    arr(ax, XL, Y3 - BH/2, XL, Y4 + BH/2)

    # Confidence annotation (right side)
    ax.text(
        XR, Y4, "Confidence labels:\nconfident / ambiguous",
        fontsize=8, color=TXT_DARK, va="center", ha="center",
        bbox=dict(
            boxstyle="round,pad=0.4",
            facecolor="#EEEAFB", edgecolor="#C7BFEA",
            linewidth=0.9
        )
    )

    # ════════════════════════════════════════════
    # 5. Bounds valid? (diamond)
    # ════════════════════════════════════════════
    Y5 = Y4 - VGAP
    DS = 0.038
    diamond = plt.Polygon(
        [(XL, Y5 + DS), (XL + DS * 1.5, Y5), (XL, Y5 - DS), (XL - DS * 1.5, Y5)],
        closed=True, facecolor=C_DEC, edgecolor=EDGE, lw=1.0, zorder=3
    )
    ax.add_patch(diamond)
    ax.text(
        XL, Y5, "Bounds\nvalid?",
        ha="center", va="center",
        fontsize=8, fontweight="bold",
        color=TXT_DARK, zorder=4, linespacing=1.1
    )
    arr(ax, XL, Y4 - BH/2, XL, Y5 + DS)

    # ── "no" → Skipped ──
    rbox(ax, XR, Y5, BWR * 0.8, BH * 0.75, "Skipped", C_SKIP,
         fontsize=9.5, text_color=TXT_DARK)
    arr(ax, XL + DS * 1.5, Y5, XR - (BWR * 0.8) / 2, Y5, color=ARROW)
    ax.text(
        (XL + DS * 1.5 + XR - (BWR * 0.8) / 2) / 2, Y5 - 0.022,
        "no", fontsize=8, color=LABEL, ha="center", fontstyle="italic"
    )

    # ════════════════════════════════════════════
    # 6. Trim & orient
    # ════════════════════════════════════════════
    Y6 = Y5 - VGAP - 0.075   # extra space for the yes label
    rbox(ax, XL, Y6, BW, BH, "Trim & orient", C_PROC,
         subtext="Reverse-complement if minus strand")
    arr(ax, XL, Y5 - DS, XL, Y6 + BH/2, color=ARROW)

    # "yes" label beside the arrow, vertically centered in the gap
    mid_56 = (Y5 - DS + Y6 + BH/2) / 2
    ax.text(
        XL + 0.025, mid_56, "yes",
        fontsize=8, color=LABEL, va="center", ha="left",
        fontstyle="italic"
    )

    # ── QC outputs (right of trim) ──
    rbox(ax, XR, Y6, BWR, BH, "QC outputs", C_QC,
         subtext="TSV / JSONL / summary", text_color=TXT_DARK)
    # Arrow from chain selection → QC
    arr(ax, XL + BW/2, Y4 - BH * 0.3, XR - BWR/2, Y6 + BH/2, color=ARROW)
    # Arrow from skipped → QC
    arr(ax, XR, Y5 - (BH * 0.75)/2, XR, Y6 + BH/2, color=ARROW)

    # ════════════════════════════════════════════
    # 7. Output files
    # ════════════════════════════════════════════
    Y7 = Y6 - VGAP
    OW = 0.15
    OH = BH * 0.78
    SPREAD = 0.19
    for label, xoff in [("ITS1", -SPREAD), ("Full ITS", 0.0), ("ITS2", SPREAD)]:
        rbox(ax, XL + xoff, Y7, OW, OH, label, C_OUT,
             fontsize=9.5, text_color=TXT_DARK)

    arr(ax, XL - 0.10, Y6 - BH/2, XL - SPREAD, Y7 + OH/2, color=ARROW)
    arr(ax, XL,        Y6 - BH/2, XL,         Y7 + OH/2, color=ARROW)
    arr(ax, XL + 0.10, Y6 - BH/2, XL + SPREAD, Y7 + OH/2, color=ARROW)

    # ── Legend ──
    handles = [
        mpatches.Patch(color=C_IN,   label="Input"),
        mpatches.Patch(color=C_HMM,  label="HMMER search"),
        mpatches.Patch(color=C_PROC, label="Processing"),
        mpatches.Patch(color=C_DEC,  label="Decision"),
        mpatches.Patch(color=C_OUT,  label="Output sequences"),
        mpatches.Patch(color=C_QC,   label="QC / diagnostics"),
        mpatches.Patch(color=C_SKIP, label="Skipped reads"),
    ]
    leg = ax.legend(
        handles=handles, loc="lower right",
        fontsize=8, framealpha=0.95,
        ncol=2, bbox_to_anchor=(1.0, 0.0),
        handlelength=1.2, handleheight=1.0,
        columnspacing=1.0
    )
    leg.get_frame().set_edgecolor(EDGE)
    leg.get_frame().set_linewidth(0.8)

    plt.savefig("fig4_pipeline.pdf", dpi=300, bbox_inches="tight", pad_inches=0.18)
    plt.savefig("fig4_pipeline.png", dpi=300, bbox_inches="tight", pad_inches=0.18)
    print("Saved fig4_pipeline.pdf / .png")


if __name__ == "__main__":
    main()
