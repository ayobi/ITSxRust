#!/usr/bin/env bash
# bench/smoke_test.sh
#
# Quick test: runs ITSxRust benchmarks only (no ITSx/ITSxpress)
# using data/subset.fasta and existing tblout to verify the pipeline.
#
# Usage:
#   bash bench/smoke_test.sh
#
set -euo pipefail

OUTDIR="bench/results_smoke"
mkdir -p "$OUTDIR"

INPUT="data/subset.fasta"
TBLOUT="data/subset.tblout"
NREADS=1000

echo "=== Smoke test: ITSxRust only, using existing tblout ==="
echo "Reads: $NREADS"
echo "$NREADS" > "$OUTDIR/nreads.txt"

# Detect time command
if command -v gtime &>/dev/null; then
    TIMECMD="gtime -v"
elif [ -x /usr/bin/time ]; then
    TIMECMD="/usr/bin/time -l"
else
    echo "WARNING: No GNU time found, timing will be approximate"
    TIMECMD="time"
fi

# --- ITSxRust: full region (default params) ---
echo ""
echo "--- itsxrust_full (default) ---"
$TIMECMD cargo run --release -- extract \
    --input "$INPUT" \
    --tblout-existing "$TBLOUT" \
    --output "$OUTDIR/itsxrust_full" \
    --region full \
    --qc-json "$OUTDIR/itsxrust_full_qc.json" \
    2> "$OUTDIR/itsxrust_full.time.raw"

# --- ITSxRust: all regions ---
echo ""
echo "--- itsxrust_all ---"
$TIMECMD cargo run --release -- extract \
    --input "$INPUT" \
    --tblout-existing "$TBLOUT" \
    --output "$OUTDIR/itsxrust_all" \
    --region all \
    --qc-json "$OUTDIR/itsxrust_all_qc.json" \
    2> "$OUTDIR/itsxrust_all.time.raw"

# --- ITSxRust: ONT preset ---
echo ""
echo "--- itsxrust_ont ---"
$TIMECMD cargo run --release -- extract \
    --input "$INPUT" \
    --tblout-existing "$TBLOUT" \
    --output "$OUTDIR/itsxrust_ont" \
    --region full \
    --preset ont \
    --qc-json "$OUTDIR/itsxrust_ont_qc.json" \
    2> "$OUTDIR/itsxrust_ont.time.raw"

# --- ITSxRust: with ambiguous filtering ---
echo ""
echo "--- itsxrust with ambiguous filtering (min_score=50) ---"
$TIMECMD cargo run --release -- extract \
    --input "$INPUT" \
    --tblout-existing "$TBLOUT" \
    --output "$OUTDIR/itsxrust_strict" \
    --region full \
    --min-anchor-score 50 \
    --write-ambiguous "$OUTDIR/itsxrust_ambiguous.fasta" \
    --qc-json "$OUTDIR/itsxrust_strict_qc.json" \
    2> "$OUTDIR/itsxrust_strict.time.raw"

# --- Summary ---
echo ""
echo "========================================="
echo "  QC JSON summaries"
echo "========================================="
for f in "$OUTDIR"/*_qc.json; do
    echo ""
    echo "--- $(basename $f) ---"
    python3 -c "
import json, sys
with open('$f') as fh:
    d = json.load(fh)
kept = d['kept']
skip = d['skipped']
print(f\"  total={d['total_reads']} kept={kept['total']} confident={kept['confident']} ambig={kept['ambiguous']} diverted={kept['ambiguous_diverted']} skipped={skip['total']}\")
print(f\"  skip reasons: {skip['by_reason']}\")
"
done

echo ""
echo "Done. Results in: $OUTDIR/"