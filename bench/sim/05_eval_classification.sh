#!/usr/bin/env bash
# 05_eval_classification.sh
#
# Evaluate downstream classification accuracy (Table 5).
#
# Strategy:
#   1. Classify the raw (untrimmed) simulated reads against UNITE with VSEARCH
#   2. Classify each tool's extracted sequences against UNITE
#   3. Compare genus/species assignments to a "true" taxonomy derived by
#      classifying the original (error-free) reference sequences
#
# Usage:
#   bash 05_eval_classification.sh \
#       <truth.tsv> <tool_outputs_dir> <unite.fasta> <results_dir>

set -euo pipefail

TRUTH="${1:?Usage: $0 <truth.tsv> <tool_outputs_dir> <unite.fasta> <results_dir>}"
TOOL_DIR="${2:?}"
UNITE_DB="${3:?}"
RESULTS="${4:?}"
THREADS="${5:-4}"

# Paths derived from earlier pipeline steps
SIM_DIR="$(dirname "$TOOL_DIR")/sim_reads"
REF_DIR="$(dirname "$TRUTH")"
REFS="$REF_DIR/references.fasta"

mkdir -p "$RESULTS"

echo "=== Classification evaluation (Table 5) ==="

# -----------------------------------------------------------------------
# Step 1: Establish "true" taxonomy by classifying error-free references
# -----------------------------------------------------------------------
echo ""
echo "--- Classifying reference sequences (ground truth taxonomy) ---"

# Extract full ITS region from references using truth coordinates
python3 - "$TRUTH" "$REFS" "$RESULTS/ref_its.fasta" << 'PYEOF'
import sys

truth_path = sys.argv[1]
refs_path = sys.argv[2]
out_path = sys.argv[3]

# Load references
refs = {}
cur_id = None
cur_seq = []
with open(refs_path) as f:
    for line in f:
        line = line.rstrip()
        if line.startswith('>'):
            if cur_id:
                refs[cur_id] = ''.join(cur_seq)
            cur_id = line[1:].split()[0]
            cur_seq = []
        else:
            cur_seq.append(line)
    if cur_id:
        refs[cur_id] = ''.join(cur_seq)

# Extract ITS regions using truth coords
with open(truth_path) as f, open(out_path, 'w') as out:
    header = f.readline()
    for line in f:
        fields = line.strip().split('\t')
        ref_id = fields[0]
        its1_start = int(fields[4])  # its1_start
        its2_end = int(fields[9])    # its2_end
        if ref_id in refs:
            seq = refs[ref_id]
            its_seq = seq[its1_start-1:its2_end]  # 1-based to 0-based
            out.write(f">{ref_id}\n{its_seq}\n")

print(f"Extracted ITS from {len(refs)} references")
PYEOF

# Classify references against UNITE
vsearch \
    --usearch_global "$RESULTS/ref_its.fasta" \
    --db "$UNITE_DB" \
    --id 0.8 \
    --maxaccepts 1 \
    --blast6out "$RESULTS/ref_classify.b6" \
    --threads "$THREADS" \
    2> "$RESULTS/ref_classify.log"

echo "  Reference classifications: $RESULTS/ref_classify.b6"

# -----------------------------------------------------------------------
# Step 2: Classify raw simulated reads and each tool's extractions
# -----------------------------------------------------------------------

classify_file() {
    local LABEL="$1"
    local INPUT="$2"
    local OUTPUT="$RESULTS/${LABEL}_classify.b6"

    if [ ! -f "$INPUT" ] || [ ! -s "$INPUT" ]; then
        echo "  $LABEL: no input file or empty, skipping"
        return
    fi

    # Convert FASTQ to FASTA if needed
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

echo ""
echo "--- Classifying simulated reads and tool outputs ---"

# Raw untrimmed
SIM_FA="${SIM_DIR}/simulated.fasta"
if [ -f "$SIM_FA" ]; then
    classify_file "raw" "$SIM_FA"
fi

# ITSxRust
for REGION in full its1 its2; do
    F="$TOOL_DIR/itsxrust/out.${REGION}.fasta"
    [ -f "$F" ] && classify_file "itsxrust_${REGION}" "$F"
done

# ITSx
for REGION_TAG in ITS1:its1 ITS2:its2 full:full; do
    ITSX_NAME="${REGION_TAG%%:*}"
    TAG="${REGION_TAG##*:}"
    F="$TOOL_DIR/itsx/itsx_out_${ITSX_NAME}.fasta"
    [ ! -f "$F" ] && F="$TOOL_DIR/itsx/itsx_out.${ITSX_NAME}.fasta"
    [ -f "$F" ] && classify_file "itsx_${TAG}" "$F"
done

# ITSxpress
for REGION in its1 its2 all; do
    TAG="$REGION"
    [ "$REGION" = "all" ] && TAG="full"
    F="$TOOL_DIR/itsxpress/${REGION}.fastq"
    [ -f "$F" ] && classify_file "itsxpress_${TAG}" "$F"
done

# -----------------------------------------------------------------------
# Step 3: Compare classifications to ground truth
# -----------------------------------------------------------------------
echo ""
echo "--- Computing classification accuracy ---"

python3 - "$RESULTS" "$SIM_DIR/read_map.tsv" << 'PYEOF'
import csv
import os
import sys
from collections import defaultdict

results_dir = sys.argv[1]
read_map_path = sys.argv[2]


def parse_unite_taxonomy(subject_id):
    """
    UNITE subject IDs contain taxonomy like:
    SH1234567.09FU|reps|k__Fungi;p__Ascomycota;c__...;g__Genus;s__Genus_species
    Extract genus and species.
    """
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
    """Load BLAST6 output. Returns dict: query_id -> (genus, species, pct_identity)."""
    results = {}
    if not os.path.exists(path):
        return results
    with open(path) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            qid = fields[0].split('|')[0]  # strip ITSxRust region tags
            sid = fields[1]
            pct_id = float(fields[2])
            genus, species = parse_unite_taxonomy(sid)
            # Keep first (best) hit per query
            if qid not in results:
                results[qid] = (genus, species, pct_id)
    return results


def trace_to_ref(query_id):
    """
    Extract the source reference ID from an extracted sequence header.
    Handles formats like "sim_read_id|its1:30-400" → "sim_read_id"
    and Badread-style "REF_NAME start-end strand" → "REF_NAME"
    """
    base = query_id.split('|')[0]
    # Badread often uses space-separated parts; first token is the read name
    base = base.split()[0]
    return base


# Load reference (ground truth) classifications
ref_class = load_b6(os.path.join(results_dir, 'ref_classify.b6'))
print(f"Reference sequences classified: {len(ref_class)}")

# Build truth lookup: sim_read_id -> ref_id is encoded in Badread read names
# We need to map through the alignment or read names
# For now, we use a simpler approach: classify extracted sequences directly
# and compare to reference classifications by aligning sim reads to refs

# For each tool output, compare genus/species to reference classification
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

# We need sim_read -> ref_id mapping. Load from minimap2 alignment.
sim_to_ref = {}
if os.path.exists(read_map_path):
    with open(read_map_path) as f:
        next(f)  # skip header
        for line in f:
            fields = line.strip().split('\t')
            sim_to_ref[fields[0]] = fields[1]

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
        # Trace back to reference
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
        'label': label,
        'tag': tag,
        'n_classified': n_total,
        'genus_correct': genus_correct,
        'species_correct': species_correct,
        'genus_accuracy': g_acc,
        'species_accuracy': s_acc,
    })

# Write summary
out_path = os.path.join(results_dir, 'table5_classification.tsv')
with open(out_path, 'w') as f:
    if summary:
        writer = csv.DictWriter(f, fieldnames=summary[0].keys(), delimiter='\t')
        writer.writeheader()
        writer.writerows(summary)
print(f"\nClassification summary: {out_path}")
PYEOF

echo ""
echo "=== Classification evaluation complete ==="
