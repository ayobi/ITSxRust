#!/usr/bin/env bash
# 03_run_tools.sh
#
# Run ITSxRust, ITSx, and ITSxpress on simulated reads.
#
# Usage:
#   bash 03_run_tools.sh <simulated.fastq> <outdir> [hmm_path] [threads]

set -euo pipefail

SIM_FQ="${1:?Usage: $0 <simulated.fastq> <outdir> [hmm_path] [threads]}"
OUTDIR="${2:?Usage: $0 <simulated.fastq> <outdir> [hmm_path] [threads]}"
HMM="${3:-F.hmm}"
THREADS="${4:-8}"

SIM_FA="${SIM_FQ%.fastq}.fasta"

mkdir -p "$OUTDIR/itsxrust" "$OUTDIR/itsx" "$OUTDIR/itsxpress"

echo "=== Running tools on simulated reads ==="
echo "  Input: $SIM_FQ"
echo "  HMM: $HMM"
echo "  Threads: $THREADS"

# -----------------------------------------------------------------------
# ITSxRust (ONT preset, all regions, with anchors TSV)
# -----------------------------------------------------------------------
echo ""
echo "--- ITSxRust ---"

ITSXRUST_BIN="${ITSXRUST_BIN:-itsxrust}"

"$ITSXRUST_BIN" extract \
    --input "$SIM_FQ" \
    --hmm "$HMM" \
    --output "$OUTDIR/itsxrust/out" \
    --region all \
    --input-format fastq \
    --output-format fasta \
    --hmmer-cpu "$THREADS" \
    --inc-e 1e-3 \
    --min-its1 30 --max-its1 1800 \
    --min-its2 30 --max-its2 2500 \
    --min-full 100 --max-full 5000 \
    --anchors-tsv "$OUTDIR/itsxrust/anchors.tsv" \
    --anchors-jsonl "$OUTDIR/itsxrust/anchors.jsonl" \
    2>&1 | tee "$OUTDIR/itsxrust/stderr.log"

echo "  ITSxRust outputs:"
ls -lh "$OUTDIR"/itsxrust/out.*.fasta 2>/dev/null || echo "    (no output files)"

# -----------------------------------------------------------------------
# ITSx
# -----------------------------------------------------------------------
echo ""
echo "--- ITSx ---"

# ITSx needs FASTA input
ITSx \
    -i "$SIM_FA" \
    -o "$OUTDIR/itsx/itsx_out" \
    -t F \
    --save_regions all \
    --detailed_results T \
    --cpu "$THREADS" \
    2>&1 | tee "$OUTDIR/itsx/stderr.log"

echo "  ITSx outputs:"
ls -lh "$OUTDIR"/itsx/itsx_out*.fasta "$OUTDIR"/itsx/itsx_out_positions.txt 2>/dev/null || true

# -----------------------------------------------------------------------
# ITSxpress (ITS1, ITS2, and ALL)
# -----------------------------------------------------------------------
echo ""
echo "--- ITSxpress ---"

for REGION in ITS1 ITS2 ALL; do
    REGION_LOWER=$(echo "$REGION" | tr '[:upper:]' '[:lower:]')
    echo "  ITSxpress --region $REGION ..."
    itsxpress \
        --fastq "$SIM_FQ" \
        --single_end \
        --region "$REGION" \
        --taxa Fungi \
        --threads "$THREADS" \
        --outfile "$OUTDIR/itsxpress/${REGION_LOWER}.fastq" \
        2>&1 | tee "$OUTDIR/itsxpress/${REGION_LOWER}.log" || true
done

echo "  ITSxpress outputs:"
ls -lh "$OUTDIR"/itsxpress/*.fastq 2>/dev/null || echo "    (no output files)"

echo ""
echo "=== All tools complete ==="