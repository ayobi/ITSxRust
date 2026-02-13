#!/usr/bin/env bash
set -euo pipefail
#
# fig1_bench_threads.sh  —  Benchmark ITSxRust and ITSx at various thread counts.
#
# Usage:
#   bash fig1_bench_threads.sh \
#     --input  data/SRR21494940.fastq \
#     --hmm    bench/sim/hmmer/F.hmm \
#     [--threads "1 2 4 8 16"] \
#     [--reps 3] \
#     [--outfile results/thread_scaling.tsv]
#
# Requires: itsxrust (in PATH or target/release/itsxrust), ITSx, /usr/bin/time or gtime

# ── Defaults ──
THREADS_LIST="1 2 4 8"
REPS=3
OUTFILE="thread_scaling.tsv"
INPUT=""
HMM=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)    INPUT="$2";        shift 2 ;;
        --hmm)      HMM="$2";         shift 2 ;;
        --threads)  THREADS_LIST="$2"; shift 2 ;;
        --reps)     REPS="$2";         shift 2 ;;
        --outfile)  OUTFILE="$2";      shift 2 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

if [[ -z "$INPUT" || -z "$HMM" ]]; then
    echo "Usage: $0 --input <fastq> --hmm <F.hmm> [--threads '1 2 4 8'] [--reps 3]"
    exit 1
fi

# Find GNU time
if command -v gtime &>/dev/null; then
    GTIME=gtime
elif /usr/bin/time --version 2>&1 | grep -q GNU; then
    GTIME=/usr/bin/time
else
    echo "ERROR: GNU time required (brew install gnu-time on macOS)"
    exit 1
fi

# Find itsxrust
ITSXRUST=""
if command -v itsxrust &>/dev/null; then
    ITSXRUST=itsxrust
elif [[ -x target/release/itsxrust ]]; then
    ITSXRUST=target/release/itsxrust
else
    echo "ERROR: itsxrust not found in PATH or target/release/"
    exit 1
fi

WORKDIR=$(mktemp -d)
trap "rm -rf $WORKDIR" EXIT

echo "Thread scaling benchmark"
echo "  Input:   $INPUT"
echo "  HMM:     $HMM"
echo "  Threads: $THREADS_LIST"
echo "  Reps:    $REPS"
echo "  Output:  $OUTFILE"
echo ""

# ── Write header ──
{
    echo -e "tool\tthreads\trep\twall_sec\tpeak_rss_kb"
} > "$OUTFILE"

for T in $THREADS_LIST; do
    echo "=== Threads: $T ==="

    for REP in $(seq 1 "$REPS"); do
        echo "  ITSxRust  rep=$REP  threads=$T"
        OUT="$WORKDIR/itsxrust_t${T}_r${REP}"
        TIME_LOG="$WORKDIR/time_itsxrust_t${T}_r${REP}.txt"

        "$GTIME" -v -o "$TIME_LOG" \
            "$ITSXRUST" extract \
                --input "$INPUT" \
                --hmm "$HMM" \
                --output "$OUT" \
                --region all \
                --hmmer-cpu "$T" \
                --inc-e 1e-3 \
                --min-its1 30 --max-its1 1800 \
                --min-its2 30 --max-its2 2500 \
                --min-full 100 --max-full 5000 \
            2>/dev/null

        # Parse GNU time output
        WALL=$(grep "Elapsed (wall clock)" "$TIME_LOG" | sed 's/.*: //')
        # Convert M:SS or H:MM:SS to seconds
        WALL_SEC=$(echo "$WALL" | awk -F: '{
            if (NF==3) print $1*3600 + $2*60 + $3;
            else if (NF==2) print $1*60 + $2;
            else print $1
        }')
        RSS=$(grep "Maximum resident" "$TIME_LOG" | sed 's/[^0-9]//g')

        echo -e "itsxrust\t$T\t$REP\t$WALL_SEC\t$RSS" >> "$OUTFILE"
        echo "    wall=${WALL_SEC}s  rss=${RSS}KB"
    done

    for REP in $(seq 1 "$REPS"); do
        echo "  ITSx      rep=$REP  threads=$T"
        ITSX_OUT="$WORKDIR/itsx_t${T}_r${REP}"
        TIME_LOG="$WORKDIR/time_itsx_t${T}_r${REP}.txt"

        # Convert FASTQ -> FASTA if needed
        FASTA_INPUT="$WORKDIR/input.fasta"
        if [[ ! -f "$FASTA_INPUT" ]]; then
            awk 'NR%4==1{print ">"substr($0,2)} NR%4==2{print}' "$INPUT" > "$FASTA_INPUT"
        fi

        "$GTIME" -v -o "$TIME_LOG" \
            ITSx -i "$FASTA_INPUT" \
                 -o "${ITSX_OUT}" \
                 -t F \
                 --cpu "$T" \
                 --save_regions all \
                 --detailed_results T \
            2>/dev/null

        WALL=$(grep "Elapsed (wall clock)" "$TIME_LOG" | sed 's/.*: //')
        WALL_SEC=$(echo "$WALL" | awk -F: '{
            if (NF==3) print $1*3600 + $2*60 + $3;
            else if (NF==2) print $1*60 + $2;
            else print $1
        }')
        RSS=$(grep "Maximum resident" "$TIME_LOG" | sed 's/[^0-9]//g')

        echo -e "itsx\t$T\t$REP\t$WALL_SEC\t$RSS" >> "$OUTFILE"
        echo "    wall=${WALL_SEC}s  rss=${RSS}KB"
    done

    echo ""
done

echo "Done. Results in: $OUTFILE"
echo ""
echo "To plot:  python fig1_thread_scaling.py --input $OUTFILE"
