# Simulation-based evaluation pipeline (Tables 4 & 5)

## Overview

This pipeline generates the data for:
- **Table 4**: Boundary accuracy on simulated ONT reads
- **Table 5**: Mock-community classification accuracy

**Strategy**: Use ITSx boundary calls on real reads (38,541 reads with all 5 regions found)
as ground truth. Simulate new ONT reads from these "reference" sequences with Badread,
run all three tools, then compare called boundaries and taxonomy to truth.

## Prerequisites

```bash
# Install tools (conda recommended)
conda install -c bioconda badread minimap2 vsearch samtools seqkit

# Python packages
pip install pandas numpy pysam

# Download UNITE database for classification evaluation
# (UNITE general FASTA release, dynamic)
wget https://doi.plutof.ut.ee/doi/10.15156/BIO/2959336 -O sh_general_release_dynamic.tar.gz
tar xzf sh_general_release_dynamic.tar.gz
# The extracted .fasta file goes in bench/sim/unite_db/
mkdir -p unite_db
mv sh_general_release_dynamic_*.fasta unite_db/unite.fasta
```

You also need:
- `data/SRR21494940.fastq` (the original ONT dataset)
- `data/itsx_out_positions.txt` (ITSx positions output from your benchmark)
- ITSxRust binary (`itsxrust`) in PATH or project target/release/
- ITSx installed (with F.hmm profiles)
- ITSxpress v2 installed

## Quick run

```bash
cd bench/sim
bash run_all.sh
```

This produces:
- `results/table4_boundary_accuracy.tsv` — boundary error stats
- `results/table5_classification.tsv` — genus/species accuracy
- `results/table4.tex` and `results/table5.tex` — ready-to-paste LaTeX

## Step-by-step

```bash
# 1. Build references + ground truth from ITSx calls
python 01_build_references.py \
    --positions data/itsx_out_positions.txt \
    --fastq data/SRR21494940.fastq \
    --outdir sim_refs \
    --n-refs 2000

# 2. Simulate ONT reads
bash 02_simulate_reads.sh sim_refs/references.fasta sim_reads 10000

# 3. Run all three tools on simulated reads
bash 03_run_tools.sh sim_reads/simulated.fastq tool_outputs

# 4. Evaluate boundary accuracy (Table 4)
python 04_eval_boundaries.py \
    --truth sim_refs/truth.tsv \
    --refs sim_refs/references.fasta \
    --itsxrust-tsv tool_outputs/itsxrust/anchors.tsv \
    --itsx-positions tool_outputs/itsx/itsx_out_positions.txt \
    --itsxpress-its1 tool_outputs/itsxpress/its1.fastq \
    --itsxpress-its2 tool_outputs/itsxpress/its2.fastq \
    --simulated sim_reads/simulated.fastq \
    --outdir results

# 5. Evaluate classification accuracy (Table 5)
bash 05_eval_classification.sh \
    sim_refs/truth.tsv \
    tool_outputs \
    unite_db/unite.fasta \
    results

# 6. Generate LaTeX tables
python 06_tabulate.py results
```

## Output structure

```
sim_refs/
  references.fasta       # 2000 reference sequences
  truth.tsv              # Ground-truth boundaries per reference
sim_reads/
  simulated.fastq        # Badread-simulated ONT reads
tool_outputs/
  itsxrust/              # ITSxRust extraction outputs
  itsx/                  # ITSx extraction outputs
  itsxpress/             # ITSxpress extraction outputs
results/
  table4_boundary_accuracy.tsv
  table5_classification.tsv
  table4.tex
  table5.tex
  boundary_errors.tsv    # Per-read boundary errors (for plots)
```
