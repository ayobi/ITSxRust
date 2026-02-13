#!/usr/bin/env bash
# run_all.sh — Master orchestrator for simulation-based evaluation
#
# Run from the bench/sim/ directory:
#   cd bench/sim && bash run_all.sh
#
# Prerequisites:
#   - badread, minimap2, vsearch, seqkit in PATH
#   - Python 3 with numpy
#   - itsxrust binary (set ITSXRUST_BIN env var or have in PATH)
#   - ITSx and itsxpress in PATH
#   - F.hmm profile (set HMM env var, default: F.hmm)
#   - Original data: data/SRR21494940.fastq, data/itsx_out_positions.txt
#     (or set DATA_DIR, POSITIONS_FILE, and FASTQ_FILE)
#   - UNITE database FASTA (set UNITE_DB env var)

set -euo pipefail

# -----------------------------------------------------------------------
# Configuration (override via environment)
# -----------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../.." && pwd)}"

DATA_DIR="${DATA_DIR:-$PROJECT_ROOT/data}"
POSITIONS_FILE="${POSITIONS_FILE:-$DATA_DIR/itsx_out.positions.txt}"
FASTQ_FILE="${FASTQ_FILE:-$DATA_DIR/SRR21494940.fastq}"
HMM="${HMM:-$SCRIPT_DIR/hmmer/F.hmm}"
UNITE_DB="${UNITE_DB:-$SCRIPT_DIR/unite_db/sh_general_release_dynamic_19.02.2025.fasta}"

N_REFS="${N_REFS:-2000}"
N_READS="${N_READS:-10000}"
THREADS="${THREADS:-8}"

WORKDIR="$SCRIPT_DIR/workdir"
RESULTS="$SCRIPT_DIR/results"

echo "============================================================"
echo "  Simulation-based evaluation pipeline"
echo "============================================================"
echo ""
echo "Configuration:"
echo "  Positions:  $POSITIONS_FILE"
echo "  FASTQ:      $FASTQ_FILE"
echo "  HMM:        $HMM"
echo "  UNITE DB:   $UNITE_DB"
echo "  N refs:     $N_REFS"
echo "  N sim reads:$N_READS"
echo "  Threads:    $THREADS"
echo "  Work dir:   $WORKDIR"
echo "  Results:    $RESULTS"
echo ""

# -----------------------------------------------------------------------
# Preflight checks
# -----------------------------------------------------------------------
echo "--- Preflight checks ---"
for CMD in python3 badread minimap2; do
    if ! command -v "$CMD" &>/dev/null; then
        echo "ERROR: $CMD not found. Install it first."
        exit 1
    fi
done

# Check optional tools (warn but don't fail)
for CMD in ITSx itsxpress vsearch; do
    if ! command -v "$CMD" &>/dev/null; then
        echo "WARNING: $CMD not found — its evaluation will be skipped"
    fi
done

ITSXRUST_BIN="${ITSXRUST_BIN:-itsxrust}"
if ! command -v "$ITSXRUST_BIN" &>/dev/null; then
    # Try common locations
    for CANDIDATE in \
        "$PROJECT_ROOT/target/release/itsxrust" \
        "$PROJECT_ROOT/target/debug/itsxrust"; do
        if [ -x "$CANDIDATE" ]; then
            ITSXRUST_BIN="$CANDIDATE"
            break
        fi
    done
fi
export ITSXRUST_BIN
echo "  ITSxRust binary: $ITSXRUST_BIN"

for F in "$POSITIONS_FILE" "$FASTQ_FILE" "$HMM"; do
    if [ ! -f "$F" ]; then
        echo "ERROR: Required file not found: $F"
        exit 1
    fi
done

python3 -c "import numpy" 2>/dev/null || {
    echo "ERROR: Python numpy not installed. Run: pip install numpy"
    exit 1
}

echo "  All checks passed."
echo ""

mkdir -p "$WORKDIR" "$RESULTS"

# -----------------------------------------------------------------------
# Step 1: Build references
# -----------------------------------------------------------------------
echo "============================================================"
echo "  Step 1: Build reference sequences + ground truth"
echo "============================================================"

python3 "$SCRIPT_DIR/01_build_references.py" \
    --positions "$POSITIONS_FILE" \
    --fastq "$FASTQ_FILE" \
    --outdir "$WORKDIR/sim_refs" \
    --n-refs "$N_REFS"

echo ""

# -----------------------------------------------------------------------
# Step 2: Simulate ONT reads
# -----------------------------------------------------------------------
echo "============================================================"
echo "  Step 2: Simulate ONT reads with Badread"
echo "============================================================"

bash "$SCRIPT_DIR/02_simulate_reads.sh" \
    "$WORKDIR/sim_refs/references.fasta" \
    "$WORKDIR/sim_reads" \
    "$N_READS"

echo ""

# -----------------------------------------------------------------------
# Step 3: Run all tools
# -----------------------------------------------------------------------
echo "============================================================"
echo "  Step 3: Run ITSxRust, ITSx, ITSxpress"
echo "============================================================"

bash "$SCRIPT_DIR/03_run_tools.sh" \
    "$WORKDIR/sim_reads/simulated.fastq" \
    "$WORKDIR/tool_outputs" \
    "$HMM" \
    "$THREADS"

echo ""

# -----------------------------------------------------------------------
# Step 4: Evaluate boundary accuracy (Table 4)
# -----------------------------------------------------------------------
echo "============================================================"
echo "  Step 4: Boundary accuracy evaluation (Table 4)"
echo "============================================================"

python3 "$SCRIPT_DIR/04_eval_boundaries.py" \
    --truth "$WORKDIR/sim_refs/truth.tsv" \
    --refs "$WORKDIR/sim_refs/references.fasta" \
    --itsxrust-dir "$WORKDIR/tool_outputs/itsxrust" \
    --itsx-dir "$WORKDIR/tool_outputs/itsx" \
    --itsxpress-dir "$WORKDIR/tool_outputs/itsxpress" \
    --outdir "$RESULTS" \
    --threads "$THREADS"

echo ""

# -----------------------------------------------------------------------
# Step 5: Classification evaluation (Table 5)
# -----------------------------------------------------------------------
if [ -f "$UNITE_DB" ]; then
    echo "============================================================"
    echo "  Step 5: Classification accuracy (Table 5)"
    echo "============================================================"

    bash "$SCRIPT_DIR/05_eval_classification.sh" \
        "$WORKDIR/sim_refs/truth.tsv" \
        "$WORKDIR/tool_outputs" \
        "$UNITE_DB" \
        "$RESULTS" \
        "$THREADS"
else
    echo "============================================================"
    echo "  Step 5: SKIPPED (UNITE database not found at $UNITE_DB)"
    echo "  To run: download UNITE and set UNITE_DB=/path/to/unite.fasta"
    echo "============================================================"
fi

echo ""

# -----------------------------------------------------------------------
# Step 6: Generate LaTeX tables
# -----------------------------------------------------------------------
echo "============================================================"
echo "  Step 6: Generate LaTeX tables"
echo "============================================================"

python3 "$SCRIPT_DIR/06_tabulate.py" "$RESULTS"

echo ""

# -----------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------
echo "============================================================"
echo "  PIPELINE COMPLETE"
echo "============================================================"
echo ""
echo "Results in: $RESULTS/"
ls -lh "$RESULTS"/*.tsv "$RESULTS"/*.tex 2>/dev/null || true
echo ""
echo "To use in paper:"
echo "  - Copy table4.tex and table5.tex content into your .tex file"
echo "  - Or use the .tsv files to make your own tables"
echo ""
echo "For plots (optional):"
echo "  - boundary_errors.tsv has per-read data for violin/ECDF plots"
