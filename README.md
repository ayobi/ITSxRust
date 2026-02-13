# ITSxRust

Fast ITS region extraction in Rust (HMMER-based), designed for long-read amplicon data (ONT / PacBio HiFi) and general FASTA/FASTQ inputs.

## Features
- Extract ITS1, ITS2, and/or full ITS region(s)
- Works with FASTA and FASTQ inputs (optionally gzipped if supported in your build)
- Produces extracted sequences plus optional boundary/anchor reporting
- Designed to be fast and reproducible

## Install

### From source (developer install)
Requires Rust (stable) and Cargo.

```bash
cargo build --release
./target/release/itsxrust --help
```

Or install into your Cargo bin dir:

```bash
cargo install --path .
itsxrust --help
```

### Planned distribution
The manuscript version will provide:
- Bioconda recipe
- Prebuilt binaries (GitHub Releases)
- Container images (GHCR)

## Usage

Basic help:

```bash
itsxrust --help
itsxrust extract --help
```

Example extraction (adjust flags to match your CLI):

```bash
itsxrust extract \
  --input data/example.fastq \
  --hmm bench/sim/hmmer/F.hmm \
  --region its2 \
  --output out_dir/
```

## Inputs / Outputs
**Inputs**
- FASTA / FASTQ
- HMM model file (HMMER)

**Outputs**
- FASTA of extracted regions (e.g. ITS1 / ITS2 / full)
- Optional tables/JSONL with anchors/boundaries (if enabled)

## Development

Run checks:

```bash
cargo fmt
cargo clippy --all-targets --all-features -- -D warnings
cargo test
```

Benchmarks & scripts live in `bench/`.

## Project layout
- `src/` Rust source
- `tests/` integration tests
- `bench/` benchmarking + simulation harness
- `testdata/` small fixtures used for tests
- `manuscript/figures/` final figure outputs

Large datasets and generated outputs should stay untracked.

## License
MIT (see `LICENSE`).

## Citation
A `CITATION.cff` will be added for the public release.
