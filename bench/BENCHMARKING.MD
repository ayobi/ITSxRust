# ITSxRust Benchmarking Guide

## Overview

The paper has 4 results tables to fill plus figure placeholders. This guide covers
how to collect data for each using the `bench/` scripts.

---

## Prerequisites

```bash
# macOS: GNU time for memory measurement
brew install gnu-time

# Conda environments
conda create -n bench-itsx -c bioconda itsx hmmer
conda create -n bench-itsxpress -c bioconda itsxpress hmmer

# ITSxRust is already built
cargo build --release
```

Verify tools:
```bash
conda run -n bench-itsx ITSx --help
conda run -n bench-itsxpress itsxpress --help
which gtime
```

---

## Datasets

### For Tables 2 and 3 (throughput + extraction success)

You need real amplicon FASTQ data at scale. Good options:

| Dataset | Platform | Reads | Accession | Notes |
|---------|----------|-------|-----------|-------|
| Your subset | ONT | 1,000 | data/subset.fasta | Too small for timing, good for smoke test |
| SRR21494940 (full) | ONT ITS | ~300K+ | SRA | The full run your subset comes from |
| Mock community ONT | ONT | varies | Check BioProject | Known composition = Table 5 too |
| HiFi ITS mock | PacBio HiFi | varies | Check BioProject | |

**Recommended minimum**: 100K+ reads for meaningful timing. Ideally 500K–5M.

Download example:
```bash
# Install SRA toolkit if needed
conda install -c bioconda sra-tools

# Download the full run
fastq-dump --split-spot --skip-technical SRR21494940 -O data/
# Or use fasterq-dump for speed
fasterq-dump SRR21494940 -O data/
```

### For Table 4 (boundary accuracy)

Requires **simulated reads with known true ITS boundaries**. Options:

1. **Badread** (ONT-like): Simulate from reference ITS+flanks sequences
2. **PBSIM2** (HiFi-like): Similar approach

Workflow:
```bash
# 1. Get reference sequences with known ITS boundaries
#    e.g., from UNITE (full-length including SSU/LSU flanks)
#    Or manually construct: SSU_3end + ITS1 + 5.8S + ITS2 + LSU_5end

# 2. Simulate ONT reads
conda install -c bioconda badread
badread simulate \
    --reference reference_its_with_flanks.fasta \
    --quantity 50000 \
    --error_model nanopore2020 \
    --length 1500,300 \
    > simulated_ont.fastq

# 3. Simulate HiFi reads
conda install -c bioconda pbsim2
pbsim2 --strategy wgs \
    --method errhmm \
    --errhmm ERRHMM-SEQUEL.model \
    --depth 100 \
    --length-min 800 --length-max 2000 \
    reference_its_with_flanks.fasta

# 4. Run all tools on simulated data
# 5. Compare extracted boundaries to known truth
```

The key is that you know the exact bp positions of ITS1/ITS2 in your reference,
so you can compute coordinate error per read.

### For Table 5 (classification accuracy)

Requires a **mock community** with known taxonomic composition:

```bash
# After extraction with each tool, classify using VSEARCH + UNITE
# Download UNITE reference
wget https://doi.plutof.ut.ee/doi/10.15156/BIO/2959336  # Check for latest

# Classify
vsearch --usearch_global extracted_its.fasta \
    --db unite_reference.fasta \
    --id 0.97 \
    --blast6out hits.tsv \
    --maxaccepts 1

# Compare assigned taxa to known mock composition
```

---

## Running Tables 2 and 3

### Quick test (your existing subset)

```bash
# Smoke test with 1000 reads — verifies the pipeline works
bash bench/run_benchmark.sh \
    data/subset.fasta \
    $CONDA_PREFIX/bin/ITSx_db/HMMs/F.hmm \
    8 \
    bench/results_subset
```

Note: ITSxpress needs FASTQ input, so it will be skipped for FASTA.
For a complete comparison, use FASTQ data.

### Full benchmark (FASTQ data)

```bash
# ONT dataset
conda activate itsxrust-hmmer
bash bench/run_benchmark.sh \
    data/ont_reads.fastq \
    $CONDA_PREFIX/bin/ITSx_db/HMMs/F.hmm \
    16 \
    bench/results_ont

# HiFi dataset (if available)
bash bench/run_benchmark.sh \
    data/hifi_reads.fastq \
    $CONDA_PREFIX/bin/ITSx_db/HMMs/F.hmm \
    16 \
    bench/results_hifi
```

### Important notes on fair comparison

- **ITSx** runs its own hmmsearch internally; the `--cpu` flag controls threads.
- **ITSxpress** also runs HMMER internally + does its own clustering.
- **ITSxRust** runs nhmmer externally, so the timed run includes nhmmer wall time.
- For fair comparison, all tools should use the **same thread count**.
- Run each tool **3 times** and report the median to account for variance.
- Make sure nothing else is running on the machine during benchmarks.

### Thread scaling (for Figure 1)

```bash
for T in 1 2 4 8 16; do
    bash bench/run_benchmark.sh data/ont_reads.fastq $HMM $T bench/results_t${T}
done
```

---

## Output files

After running, `bench/results_*/` will contain:

```
table2_throughput.tsv    # Tab-separated: tool, wall_sec, reads/sec, RSS
table3_extraction.tsv    # Tab-separated: tool, %ITS1, %ITS2, %full, %ambig
itsxrust_*_qc.json       # Full QC summaries from ITSxRust
*.time                   # Raw timing data per tool
*.time.raw               # Full GNU time output
```

---

## Filling in the paper tables

Once you have TSV results:

1. **Table 2**: Copy from `table2_throughput.tsv` — one row per tool per dataset
2. **Table 3**: Copy from `table3_extraction.tsv` — note ITSxRust reports ambiguous %
3. **Table 4**: Requires simulation pipeline (see above); compute median/mean/95th abs error
4. **Table 5**: Requires mock + UNITE classification; compute genus/species accuracy

Tables 4 and 5 are the most labor-intensive. If you want to ship the paper
sooner, you could initially report only Tables 2 and 3 and add 4 and 5 in revision.

---

## Figures

- **Figure 1** (throughput vs threads): Plot reads/sec from thread scaling runs
- **Figure 2** (failure code breakdown): Use the `by_reason` field from QC JSON
- **Figure 3** (boundary error ECDF): From Table 4 simulation data
- **Figure 4** (pipeline schematic): Diagram, not data-driven

For Figures 1 and 2, the QC JSON has everything you need. A simple matplotlib
or R script can generate them from the TSV/JSON outputs.