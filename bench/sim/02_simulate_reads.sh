#!/usr/bin/env bash
# 02_simulate_reads.sh
#
# Simulate ONT-like reads from reference sequences using Badread.
#
# Usage:
#   bash 02_simulate_reads.sh <references.fasta> <outdir> [n_reads]
#
# Badread encodes the source reference and coordinates in read names:
#   @<ref_id>,<ref_start>,<ref_end>,<strand>  (0-based coords)
# This lets us trace each simulated read back to its source.

set -euo pipefail

REFS="${1:?Usage: $0 <references.fasta> <outdir> [n_reads]}"
OUTDIR="${2:?Usage: $0 <references.fasta> <outdir> [n_reads]}"
N_READS="${3:-10000}"

mkdir -p "$OUTDIR"

# Badread's --quantity is in TARGET BASES (not read count).
# Calculate: desired_reads × mean_read_length = target_bases
MEAN_LEN=800
TARGET_BASES=$(( N_READS * MEAN_LEN ))

echo "=== Simulating ~$N_READS ONT reads from $(grep -c '^>' "$REFS") references ==="
echo "  Target bases: ${TARGET_BASES} (${N_READS} reads × ${MEAN_LEN} bp mean length)"

# Badread simulate parameters:
#   --quantity:  target total bases of output
#   --identity:  min,max,shape for read accuracy (ONT R10: ~92-98%)
#   --length:    mean,stdev for read lengths
#   --start_adapter_seq / --end_adapter_seq: "" to skip adapters

badread simulate \
    --reference "$REFS" \
    --quantity "$TARGET_BASES" \
    --length 800,200 \
    --identity 92,98,4 \
    --start_adapter_seq "" \
    --end_adapter_seq "" \
    --junk_reads 0 \
    --random_reads 0 \
    --chimeras 0 \
    --seed 42 \
    > "$OUTDIR/simulated.fastq" \
    2> "$OUTDIR/badread.log"

N_SIM=$(grep -c '^@' "$OUTDIR/simulated.fastq" || true)
echo "  Generated $N_SIM reads -> $OUTDIR/simulated.fastq"

# Also create FASTA version (needed by ITSx)
awk 'NR%4==1{sub(/^@/,">",$0); print} NR%4==2{print}' \
    "$OUTDIR/simulated.fastq" > "$OUTDIR/simulated.fasta"

echo "  FASTA copy -> $OUTDIR/simulated.fasta"

# Parse Badread read names to create a read-to-reference mapping
# Also create an alignment-based mapping (more robust across Badread versions)
echo "  Building read-to-reference mapping via minimap2..."

minimap2 -a --secondary=no -t 4 \
    "$REFS" "$OUTDIR/simulated.fasta" 2>/dev/null \
    | awk -F'\t' '!/^@/ && $2 != 4 { print $1 "\t" $3 }' \
    > "$OUTDIR/read_map.tsv"

# Add header (macOS-compatible: use temp file instead of sed -i)
{ echo -e "sim_read_id\tref_id"; cat "$OUTDIR/read_map.tsv"; } > "$OUTDIR/read_map.tmp"
mv "$OUTDIR/read_map.tmp" "$OUTDIR/read_map.tsv"

N_MAPPED=$(tail -n+2 "$OUTDIR/read_map.tsv" | wc -l)
echo "  Mapped $N_MAPPED / $N_SIM reads to references -> $OUTDIR/read_map.tsv"

echo "  Read map -> $OUTDIR/read_map.tsv"
echo "=== Simulation complete ==="