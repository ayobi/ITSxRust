# ITSxRust

ITS subregion extraction for fungal metabarcoding at long-read scale.

As long-read amplicon sequencing (Oxford Nanopore and PacBio HiFi) becomes routine, extracting ITS subregions (ITS1, 5.8S, ITS2, full ITS) reliably at scale can become a throughput and robustness bottleneck. ITSxRust is a Rust-based ITS extractor that follows the standard approach of locating conserved ribosomal flanks using profile-HMMs (via HMMER), while adding long-read–oriented features for reproducible, high-throughput processing.

## Features
- HMMER/profile-HMM–based detection of conserved ribosomal flanks to extract ITS subregions
- Supports long-read workloads (ONT / HiFi) with built-in parameter presets
- Optional dereplication to reduce redundant HMMER searches
- Partial-chain fallback: recover subregions using two-anchor pairs when a full four-anchor chain is unavailable
- Structured failure diagnostics and QC summaries to help understand why reads were skipped or partially recovered
- Works with FASTA and FASTQ inputs

## Install

### Prebuilt binaries (recommended)
Download the appropriate binary for your OS from GitHub Releases:

- GitHub → Releases → `v0.1.0`

Then:

```bash
chmod +x itsxrust
./itsxrust --help
```

### From source
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

### Dependency: HMMER
ITSxRust coordinates HMMER searches (e.g., `hmmscan`) to locate ribosomal flanks. Ensure HMMER is available in your environment for typical extraction workflows.

## Usage

Help:

```bash
itsxrust --help
itsxrust extract --help
```

Example extraction:

```bash
itsxrust extract   --input reads.fastq.gz   --hmm path/to/F.hmm   --region its2   --output out_dir/   --hmmer-cpu 8
```

Presets (ONT / HiFi) are available via the CLI (see `itsxrust extract --help`).

## Inputs / Outputs

**Inputs**
- FASTA / FASTQ (optionally gzipped)
- HMM model file (profile-HMMs for ribosomal flanks)

**Outputs**
- FASTA of extracted regions (ITS1 / ITS2 / full)
- Optional anchor/boundary outputs (TSV/JSONL) and QC summaries (if enabled)

## Development

Run checks:

```bash
cargo fmt
cargo clippy --all-targets --all-features -- -D warnings
cargo test
```

Benchmarks and simulation scripts live in `bench/`.

## Project layout
- `src/` Rust source
- `tests/` integration tests
- `bench/` benchmarking + simulation harness
- `testdata/` small fixtures used for tests
- `manuscript/figures/` final figure outputs

Large datasets and generated outputs should stay untracked.

## Roadmap
- Container images (GHCR)
- Bioconda recipe

## License
MIT (see `LICENSE`).

## Citation
If you use ITSxRust, please cite the repository metadata via GitHub’s “Cite this repository” button (powered by `CITATION.cff`).
