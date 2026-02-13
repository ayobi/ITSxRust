#!/usr/bin/env bash
# bench/run_benchmark.sh
#
# Benchmark ITSxRust vs ITSx vs ITSxpress v2.
# Collects data for paper Tables 2 (throughput/resource) and 3 (extraction success).
#
# Prerequisites:
#   brew install gnu-time          # macOS: provides gtime
#   conda install -c bioconda itsx itsxpress
#   cargo build --release          # ITSxRust binary
#
# Usage:
#   bash bench/run_benchmark.sh <input.fastq> <hmm_profile> <threads> <out_dir>
#
# Example:
#   bash bench/run_benchmark.sh data/ont_5M.fastq.gz $CONDA_PREFIX/bin/ITSx_db/HMMs/F.hmm 16 bench/results
#
set -euo pipefail

# ---------------------------------------------------------------------------
# Args
# ---------------------------------------------------------------------------
INPUT="${1:?Usage: $0 <input.fasta|fastq[.gz]> <hmm> <threads> <outdir>}"
HMM="${2:?Provide HMM profile path (e.g. F.hmm)}"
THREADS="${3:-16}"
OUTDIR="${4:-bench/results}"

mkdir -p "$OUTDIR"

# ---------------------------------------------------------------------------
# Detect input format
# ---------------------------------------------------------------------------
BASENAME=$(basename "$INPUT")
EXT="${BASENAME##*.}"

IS_GZ=false
if [ "$EXT" = "gz" ]; then
    IS_GZ=true
    INNER="${BASENAME%.gz}"
    EXT="${INNER##*.}"
fi

IS_FASTQ=false
if [[ "$EXT" == "fq" || "$EXT" == "fastq" ]]; then
    IS_FASTQ=true
fi

# ---------------------------------------------------------------------------
# Count reads
# ---------------------------------------------------------------------------
echo "=== Counting reads in $INPUT ==="
if $IS_GZ; then
    if $IS_FASTQ; then
        NREADS=$(zcat < "$INPUT" | awk 'NR%4==1' | wc -l | tr -d ' ')
    else
        NREADS=$(zcat < "$INPUT" | grep -c '^>' || true)
    fi
else
    if $IS_FASTQ; then
        NREADS=$(awk 'NR%4==1' "$INPUT" | wc -l | tr -d ' ')
    else
        NREADS=$(grep -c '^>' "$INPUT" || true)
    fi
fi
echo "Total reads: $NREADS"
echo "$NREADS" > "$OUTDIR/nreads.txt"

# ---------------------------------------------------------------------------
# Detect time command
# ---------------------------------------------------------------------------
if command -v gtime &>/dev/null; then
    TIMECMD="gtime"
elif command -v /usr/bin/time &>/dev/null; then
    TIMECMD="/usr/bin/time"
else
    echo "ERROR: Neither gtime nor /usr/bin/time found. Install gnu-time."
    exit 1
fi

# Helper: run a command, capture wall time + peak RSS
# Writes: <tool>.time containing wall_sec and peak_rss_kb
run_timed() {
    local label="$1"
    shift
    local timefile="$OUTDIR/${label}.time"

    # Skip completed steps unless FORCE=1
    if [[ -s "$timefile" && "${FORCE:-0}" != "1" ]]; then
        echo "--- SKIP: $label (found $timefile)"
        echo ""
        return 0
    fi

    echo "--- Running: $label ---"
    echo "CMD: $*"

    # GNU time -v gives "Maximum resident set size (kbytes)" and "Elapsed (wall clock)"
    # macOS /usr/bin/time -l gives "maximum resident set size" in bytes
    if [ "$TIMECMD" = "gtime" ]; then
        $TIMECMD -v "$@" 2> "$timefile.raw"
        WALL=$(grep "Elapsed (wall clock)" "$timefile.raw" | sed 's/.*: //')
        RSS_KB=$(grep "Maximum resident set size" "$timefile.raw" | awk '{print $NF}')
    else
        # macOS /usr/bin/time -l
        $TIMECMD -l "$@" 2> "$timefile.raw"
        WALL=$(grep "real" "$timefile.raw" | awk '{print $1}')
        RSS_BYTES=$(grep "maximum resident set size" "$timefile.raw" | awk '{print $1}')
        RSS_KB=$((RSS_BYTES / 1024))
    fi

    echo "wall=$WALL rss_kb=$RSS_KB" | tee "$timefile"
    echo ""
}

# ---------------------------------------------------------------------------
# Prepare FASTA for ITSx (it prefers FASTA)
# ---------------------------------------------------------------------------
if $IS_FASTQ; then
    echo "=== Converting FASTQ -> FASTA for ITSx ==="
    FASTA_INPUT="$OUTDIR/input_for_itsx.fasta"
    if $IS_GZ; then
        zcat < "$INPUT" | awk 'NR%4==1{sub(/^@/,">")} NR%4==1||NR%4==2' > "$FASTA_INPUT"
    else
        awk 'NR%4==1{sub(/^@/,">")} NR%4==1||NR%4==2' "$INPUT" > "$FASTA_INPUT"
    fi
else
    if $IS_GZ; then
        FASTA_INPUT="$OUTDIR/input_for_itsx.fasta"
        zcat < "$INPUT" > "$FASTA_INPUT"
    else
        FASTA_INPUT="$INPUT"
    fi
fi

# Also need uncompressed FASTQ for ITSxpress
if $IS_FASTQ; then
    if $IS_GZ; then
        FASTQ_INPUT="$OUTDIR/input_uncompressed.fastq"
        zcat < "$INPUT" > "$FASTQ_INPUT"
    else
        FASTQ_INPUT="$INPUT"
    fi
fi

# ---------------------------------------------------------------------------
# 1. ITSxRust
# ---------------------------------------------------------------------------
echo ""
echo "========================================="
echo "  ITSxRust"
echo "========================================="

# Full ITS
run_timed "itsxrust_full" \
    cargo run --release -- extract \
    --input "$INPUT" \
    --hmm "$HMM" \
    --output "$OUTDIR/itsxrust_full" \
    --region full \
    --hmmer-cpu "$THREADS" \
    --qc-json "$OUTDIR/itsxrust_full_qc.json"

# All regions (ITS1 + ITS2 + full)
run_timed "itsxrust_all" \
    cargo run --release -- extract \
    --input "$INPUT" \
    --hmm "$HMM" \
    --output "$OUTDIR/itsxrust_all" \
    --region all \
    --hmmer-cpu "$THREADS" \
    --qc-json "$OUTDIR/itsxrust_all_qc.json"

# With dereplication
run_timed "itsxrust_derep" \
    cargo run --release -- extract \
    --input "$INPUT" \
    --hmm "$HMM" \
    --output "$OUTDIR/itsxrust_derep" \
    --region full \
    --hmmer-cpu "$THREADS" \
    --derep \
    --qc-json "$OUTDIR/itsxrust_derep_qc.json"

# ONT preset
run_timed "itsxrust_ont" \
    cargo run --release -- extract \
    --input "$INPUT" \
    --hmm "$HMM" \
    --output "$OUTDIR/itsxrust_ont" \
    --region full \
    --hmmer-cpu "$THREADS" \
    --preset ont \
    --qc-json "$OUTDIR/itsxrust_ont_qc.json"

# ---------------------------------------------------------------------------
# 2. ITSx
# ---------------------------------------------------------------------------
echo ""
echo "========================================="
echo "  ITSx"
echo "========================================="

run_timed "itsx" \
    ITSx -i "$FASTA_INPUT" \
    -o "$OUTDIR/itsx_out" \
    --cpu "$THREADS" \
    -t F \
    --save_regions all \
    --detailed_results T

# ---------------------------------------------------------------------------
# 3. ITSxpress v2
# ---------------------------------------------------------------------------
echo ""
echo "========================================="
echo "  ITSxpress v2"
echo "========================================="

if $IS_FASTQ; then
    # ITSxpress ITS1
    run_timed "itsxpress_its1" \
        itsxpress --fastq "$FASTQ_INPUT" \
        --single_end \
        --outfile "$OUTDIR/itsxpress_its1.fastq" \
        --region ITS1 \
        --taxa Fungi \
        --threads "$THREADS"

    # ITSxpress ITS2
    run_timed "itsxpress_its2" \
        itsxpress --fastq "$FASTQ_INPUT" \
        --single_end \
        --outfile "$OUTDIR/itsxpress_its2.fastq" \
        --region ITS2 \
        --taxa Fungi \
        --threads "$THREADS"

    # ITSxpress ALL (full ITS)
    run_timed "itsxpress_all" \
        itsxpress --fastq "$FASTQ_INPUT" \
        --single_end \
        --outfile "$OUTDIR/itsxpress_all.fastq" \
        --region ALL \
        --taxa Fungi \
        --threads "$THREADS"
else
    echo "SKIP: ITSxpress requires FASTQ input. Run with FASTQ to benchmark ITSxpress."
fi

# ---------------------------------------------------------------------------
# Collect results
# ---------------------------------------------------------------------------
echo ""
echo "========================================="
echo "  Collecting results"
echo "========================================="

python3 bench/collect_results.py "$OUTDIR" "$NREADS"

echo ""
echo "Done. Results in: $OUTDIR/"