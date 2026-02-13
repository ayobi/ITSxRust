#!/usr/bin/env bash
# run_hifi_eval.sh
#
# Run the full HiFi simulation evaluation pipeline:
#   1. Subsample 165K reads to ~15K (to match ONT sim size)
#   2. Convert FASTQ → FASTA (for ITSx/nhmmer)
#   3. Run ITSxRust, ITSx, ITSxpress
#   4. Evaluate boundary accuracy (Table 4)
#   5. Evaluate classification accuracy (Table 5)
#
# Prerequisites:
#   - HiFi sim already generated: bench/sim/workdir/hifi_sim.fastq
#   - ONT sim references exist:   bench/sim/workdir/sim_refs/{references.fasta, truth.tsv}
#   - UNITE DB:                   bench/sim/unite_db/unite.fasta
#   - Tools in PATH: itsxrust, ITSx, itsxpress, minimap2, vsearch
#
# Usage (from project root):
#   bash bench/sim/run_hifi_eval.sh [threads]

set -euo pipefail

THREADS="${1:-8}"

# ---------- paths (adjust if your layout differs) ----------
PROJ_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
SIM_DIR="$PROJ_ROOT/bench/sim"
WORKDIR="$SIM_DIR/workdir"

HIFI_FULL="$WORKDIR/hifi_sim.fastq"
HIFI_SUB="$WORKDIR/hifi_sim_15k.fastq"
HIFI_FA="$WORKDIR/hifi_sim_15k.fasta"

REFS="$WORKDIR/sim_refs/references.fasta"
TRUTH="$WORKDIR/sim_refs/truth.tsv"

TOOL_OUT="$WORKDIR/hifi_tool_outputs"
RESULTS="$WORKDIR/hifi_results"

HMM="${HMM:-F.hmm}"
ITSXRUST_BIN="${ITSXRUST_BIN:-itsxrust}"
UNITE_DB="$SIM_DIR/unite_db/unite.fasta"

# ---------- sanity checks ----------
for F in "$HIFI_FULL" "$REFS" "$TRUTH"; do
    if [ ! -f "$F" ]; then
        echo "ERROR: missing required file: $F" >&2
        exit 1
    fi
done

echo "========================================================"
echo " HiFi simulation evaluation pipeline"
echo "========================================================"
echo "  Threads:    $THREADS"
echo "  HiFi FASTQ: $HIFI_FULL ($(grep -c '^@' "$HIFI_FULL" || true) reads)"
echo "  References: $REFS"
echo "  Truth:      $TRUTH"
echo "  HMM:        $HMM"
echo "  UNITE DB:   $UNITE_DB"
echo ""

# =======================================================================
# Step 1: Subsample to ~15,000 reads
# =======================================================================
echo "--- Step 1: Subsample to ~15K reads ---"

if [ -f "$HIFI_SUB" ]; then
    N_SUB=$(awk 'NR%4==1' "$HIFI_SUB" | wc -l)
    echo "  Already exists: $HIFI_SUB ($N_SUB reads), skipping"
else
    # 15000 reads × 4 lines = 60000 lines
    head -60000 "$HIFI_FULL" > "$HIFI_SUB"
    N_SUB=$(awk 'NR%4==1' "$HIFI_SUB" | wc -l)
    echo "  Created: $HIFI_SUB ($N_SUB reads)"
fi

# =======================================================================
# Step 2: Convert FASTQ → FASTA (for ITSx)
# =======================================================================
echo ""
echo "--- Step 2: Convert FASTQ → FASTA ---"

if [ -f "$HIFI_FA" ]; then
    echo "  Already exists: $HIFI_FA, skipping"
else
    awk 'NR%4==1{sub(/^@/,">",$0); print} NR%4==2{print}' "$HIFI_SUB" > "$HIFI_FA"
    echo "  Created: $HIFI_FA"
fi

# =======================================================================
# Step 3: Run all three tools
# =======================================================================
mkdir -p "$TOOL_OUT/itsxrust" "$TOOL_OUT/itsx" "$TOOL_OUT/itsxpress"

# --- ITSxRust ---
echo ""
echo "--- Step 3a: ITSxRust ---"

if [ -f "$TOOL_OUT/itsxrust/out.full.fasta" ]; then
    echo "  Output exists, skipping. Delete $TOOL_OUT/itsxrust/ to re-run."
else
    "$ITSXRUST_BIN" extract \
        --input "$HIFI_SUB" \
        --hmm "$HMM" \
        --output "$TOOL_OUT/itsxrust/out" \
        --region all \
        --input-format fastq \
        --output-format fasta \
        --hmmer-cpu "$THREADS" \
        --inc-e 1e-3 \
        --min-its1 30 --max-its1 1800 \
        --min-its2 30 --max-its2 2500 \
        --min-full 100 --max-full 5000 \
        --anchors-tsv "$TOOL_OUT/itsxrust/anchors.tsv" \
        --anchors-jsonl "$TOOL_OUT/itsxrust/anchors.jsonl" \
        2>&1 | tee "$TOOL_OUT/itsxrust/stderr.log"
fi

echo "  ITSxRust outputs:"
ls -lh "$TOOL_OUT"/itsxrust/out.*.fasta 2>/dev/null || echo "    (no output files)"

# --- ITSx ---
echo ""
echo "--- Step 3b: ITSx ---"

if [ -f "$TOOL_OUT/itsx/itsx_out_positions.txt" ]; then
    echo "  Output exists, skipping. Delete $TOOL_OUT/itsx/ to re-run."
else
    ITSx \
        -i "$HIFI_FA" \
        -o "$TOOL_OUT/itsx/itsx_out" \
        -t F \
        --save_regions all \
        --detailed_results T \
        --cpu "$THREADS" \
        2>&1 | tee "$TOOL_OUT/itsx/stderr.log"
fi

echo "  ITSx outputs:"
ls -lh "$TOOL_OUT"/itsx/itsx_out*.fasta "$TOOL_OUT"/itsx/itsx_out_positions.txt 2>/dev/null || true

# --- ITSxpress ---
echo ""
echo "--- Step 3c: ITSxpress ---"

for REGION in ITS1 ITS2 ALL; do
    REGION_LC="$(printf '%s' "$REGION" | tr '[:upper:]' '[:lower:]')"
    OUTFILE="$TOOL_OUT/itsxpress/${REGION_LC}.fastq"

    if [ -f "$OUTFILE" ]; then
        echo "  ITSxpress $REGION: output exists, skipping"
        continue
    fi

    echo "  ITSxpress --region $REGION ..."
    itsxpress \
        --fastq "$HIFI_SUB" \
        --single_end \
        --region "$REGION" \
        --taxa Fungi \
        --threads "$THREADS" \
        --outfile "$OUTFILE" \
        2>&1 | tee "$TOOL_OUT/itsxpress/${REGION_LC}.log" || true
done


echo "  ITSxpress outputs:"
ls -lh "$TOOL_OUT"/itsxpress/*.fastq 2>/dev/null || echo "    (no output files)"

echo ""
echo "=== All tools complete ==="

# =======================================================================
# Step 4: Evaluate boundary accuracy
# =======================================================================
echo ""
echo "--- Step 4: Boundary accuracy evaluation ---"

mkdir -p "$RESULTS"

python3 "$SIM_DIR/04_eval_boundaries.py" \
    --truth "$TRUTH" \
    --refs "$REFS" \
    --itsxrust-dir "$TOOL_OUT/itsxrust" \
    --itsx-dir "$TOOL_OUT/itsx" \
    --itsxpress-dir "$TOOL_OUT/itsxpress" \
    --outdir "$RESULTS" \
    --threads "$THREADS"

# =======================================================================
# Step 5: Evaluate classification accuracy
# =======================================================================
echo ""
echo "--- Step 5: Classification accuracy evaluation ---"

if [ ! -f "$UNITE_DB" ]; then
    echo "  WARNING: UNITE DB not found at $UNITE_DB"
    echo "  Skipping classification evaluation."
    echo "  To enable: place UNITE FASTA at $UNITE_DB"
else
    # The 05 script expects a specific directory layout for SIM_DIR.
    # We'll set up the raw simulated FASTA it needs.
    RAW_SIM_FA="$WORKDIR/hifi_sim_reads/simulated.fasta"
    mkdir -p "$WORKDIR/hifi_sim_reads"
    if [ ! -f "$RAW_SIM_FA" ]; then
        cp "$HIFI_FA" "$RAW_SIM_FA"
    fi

    # The 05 script derives SIM_DIR from tool_outputs location.
    # We need to adjust paths or call inline. Let's call inline for clarity.
    # Re-use the ref_classify.b6 from the ONT run if it exists:
    ONT_REF_B6="$WORKDIR/results/ref_classify.b6"
    if [ -f "$ONT_REF_B6" ] && [ ! -f "$RESULTS/ref_classify.b6" ]; then
        echo "  Reusing reference classifications from ONT run"
        cp "$ONT_REF_B6" "$RESULTS/ref_classify.b6"
        [ -f "$WORKDIR/results/ref_its.fasta" ] && cp "$WORKDIR/results/ref_its.fasta" "$RESULTS/ref_its.fasta"
    fi

    # Classify raw + tool outputs
    classify_file() {
        local LABEL="$1"
        local INPUT="$2"
        local OUTPUT="$RESULTS/${LABEL}_classify.b6"

        if [ ! -f "$INPUT" ] || [ ! -s "$INPUT" ]; then
            echo "  $LABEL: no input file or empty, skipping"
            return
        fi

        local FA_INPUT="$INPUT"
        if [[ "$INPUT" == *.fastq ]] || [[ "$INPUT" == *.fq ]]; then
            FA_INPUT="$RESULTS/${LABEL}_input.fasta"
            awk 'NR%4==1{sub(/^@/,">",$0); print} NR%4==2{print}' "$INPUT" > "$FA_INPUT"
        fi

        vsearch \
            --usearch_global "$FA_INPUT" \
            --db "$UNITE_DB" \
            --id 0.8 \
            --maxaccepts 1 \
            --blast6out "$OUTPUT" \
            --threads "$THREADS" \
            2> "$RESULTS/${LABEL}_classify.log"

        local N_CLASS=$(wc -l < "$OUTPUT")
        echo "  $LABEL: $N_CLASS classified"
    }

    # Classify reference ITS if not already done
    if [ ! -f "$RESULTS/ref_classify.b6" ]; then
        echo "  Classifying reference sequences..."
        python3 -c "
import sys
truth_path, refs_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
refs = {}; cur_id = None; cur_seq = []
with open(refs_path) as f:
    for line in f:
        line = line.rstrip()
        if line.startswith('>'):
            if cur_id: refs[cur_id] = ''.join(cur_seq)
            cur_id = line[1:].split()[0]; cur_seq = []
        else: cur_seq.append(line)
    if cur_id: refs[cur_id] = ''.join(cur_seq)
with open(truth_path) as f, open(out_path, 'w') as out:
    header = f.readline()
    for line in f:
        fields = line.strip().split('\t')
        ref_id = fields[0]
        its1_start = int(fields[4]); its2_end = int(fields[9])
        if ref_id in refs:
            seq = refs[ref_id]
            out.write(f'>{ref_id}\n{seq[its1_start-1:its2_end]}\n')
print(f'Extracted ITS from {len(refs)} references')
" "$TRUTH" "$REFS" "$RESULTS/ref_its.fasta"

        vsearch \
            --usearch_global "$RESULTS/ref_its.fasta" \
            --db "$UNITE_DB" \
            --id 0.8 \
            --maxaccepts 1 \
            --blast6out "$RESULTS/ref_classify.b6" \
            --threads "$THREADS" \
            2> "$RESULTS/ref_classify.log"
    fi

    echo ""
    echo "  Classifying tool outputs..."

    # Raw untrimmed
    classify_file "raw" "$RAW_SIM_FA"

    # ITSxRust
    for REGION in full its1 its2; do
        F="$TOOL_OUT/itsxrust/out.${REGION}.fasta"
        [ -f "$F" ] && classify_file "itsxrust_${REGION}" "$F"
    done

    # ITSx
    for REGION_TAG in ITS1:its1 ITS2:its2 full:full; do
        ITSX_NAME="${REGION_TAG%%:*}"
        TAG="${REGION_TAG##*:}"
        F="$TOOL_OUT/itsx/itsx_out_${ITSX_NAME}.fasta"
        [ ! -f "$F" ] && F="$TOOL_OUT/itsx/itsx_out.${ITSX_NAME}.fasta"
        [ -f "$F" ] && classify_file "itsx_${TAG}" "$F"
    done

    # ITSxpress
    for REGION in its1 its2 all; do
        TAG="$REGION"
        [ "$REGION" = "all" ] && TAG="full"
        F="$TOOL_OUT/itsxpress/${REGION}.fastq"
        [ -f "$F" ] && classify_file "itsxpress_${TAG}" "$F"
    done

    # Compute accuracy
    echo ""
    echo "  Computing classification accuracy..."

    python3 - "$RESULTS" "$WORKDIR" << 'PYEOF'
import csv, os, sys

results_dir = sys.argv[1]
workdir = sys.argv[2]

def parse_unite_taxonomy(subject_id):
    parts = subject_id.split('|')
    tax_str = parts[-1] if len(parts) > 1 else subject_id
    genus = species = "unclassified"
    for field in tax_str.split(';'):
        field = field.strip()
        if field.startswith('g__') and len(field) > 3:
            genus = field[3:]
        elif field.startswith('s__') and len(field) > 3:
            species = field[3:]
    return genus, species

def load_b6(path):
    results = {}
    if not os.path.exists(path):
        return results
    with open(path) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            qid = fields[0].split('|')[0]
            sid = fields[1]
            pct_id = float(fields[2])
            genus, species = parse_unite_taxonomy(sid)
            if qid not in results:
                results[qid] = (genus, species, pct_id)
    return results

def trace_to_ref(query_id):
    base = query_id.split('|')[0]
    base = base.split()[0]
    return base

# Load reference classifications
ref_class = load_b6(os.path.join(results_dir, 'ref_classify.b6'))
print(f"Reference sequences classified: {len(ref_class)}")

# Load sim_read -> ref_id mapping if available
# Try both ONT and HiFi read map locations
sim_to_ref = {}
for subdir in ['sim_reads', 'hifi_sim_reads']:
    read_map_path = os.path.join(workdir, subdir, 'read_map.tsv')
    if os.path.exists(read_map_path):
        with open(read_map_path) as f:
            next(f)
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 2:
                    sim_to_ref[fields[0]] = fields[1]

labels = [
    ('raw', 'Raw (untrimmed)'),
    ('itsxrust_full', 'ITSxRust (full)'),
    ('itsxrust_its1', 'ITSxRust (ITS1)'),
    ('itsxrust_its2', 'ITSxRust (ITS2)'),
    ('itsx_full', 'ITSx (full)'),
    ('itsx_its1', 'ITSx (ITS1)'),
    ('itsx_its2', 'ITSx (ITS2)'),
    ('itsxpress_full', 'ITSxpress (full)'),
    ('itsxpress_its1', 'ITSxpress (ITS1)'),
    ('itsxpress_its2', 'ITSxpress (ITS2)'),
]

summary = []
for tag, label in labels:
    b6_path = os.path.join(results_dir, f'{tag}_classify.b6')
    tool_class = load_b6(b6_path)
    if not tool_class:
        continue

    n_total = 0
    genus_correct = 0
    species_correct = 0

    for qid, (g, s, pct) in tool_class.items():
        base_id = trace_to_ref(qid)
        ref_id = sim_to_ref.get(base_id, base_id)
        if ref_id not in ref_class:
            continue
        ref_g, ref_s, _ = ref_class[ref_id]
        n_total += 1
        if g == ref_g and g != "unclassified":
            genus_correct += 1
        if s == ref_s and s != "unclassified":
            species_correct += 1

    if n_total > 0:
        g_acc = genus_correct / n_total * 100
        s_acc = species_correct / n_total * 100
    else:
        g_acc = s_acc = 0.0

    print(f"  {label:<25} N={n_total:>6}  Genus={g_acc:.1f}%  Species={s_acc:.1f}%")
    summary.append({
        'label': label, 'tag': tag,
        'n_classified': n_total,
        'genus_correct': genus_correct,
        'species_correct': species_correct,
        'genus_accuracy': g_acc,
        'species_accuracy': s_acc,
    })

out_path = os.path.join(results_dir, 'table5_classification.tsv')
with open(out_path, 'w') as f:
    if summary:
        writer = csv.DictWriter(f, fieldnames=summary[0].keys(), delimiter='\t')
        writer.writeheader()
        writer.writerows(summary)
print(f"\nClassification summary: {out_path}")
PYEOF
fi

# =======================================================================
# Summary
# =======================================================================
echo ""
echo "========================================================"
echo " HiFi evaluation complete!"
echo "========================================================"
echo ""
echo "Results directory: $RESULTS"
echo ""
echo "Key output files:"
echo "  $RESULTS/boundary_errors.tsv           (per-read boundary errors)"
echo "  $RESULTS/table4_boundary_accuracy.tsv   (summary stats)"
echo "  $RESULTS/table5_classification.tsv      (classification accuracy)"
echo ""
echo "Tool outputs:  $TOOL_OUT/{itsxrust,itsx,itsxpress}/"
echo ""
echo "To compare with ONT results:"
echo "  ONT: $WORKDIR/results/"
echo "  HiFi: $RESULTS/"
