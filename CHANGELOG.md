# Changelog

All notable changes to this project will be documented in this file.

The format is based on Keep a Changelog, and this project follows Semantic Versioning.

## [Unreleased]

## [0.3.0] - 2026-07-28

### Added
- `--region ssu` and `--region lsu` for flanking ribosomal gene export.
- `--plain-ids` to emit output with the input read ID unchanged, omitting the
  default `|<region>:<start>-<end>` suffix, for joining output back to input.
- Unrecognised profile names are counted, sampled and reported on stderr and in
  the QC JSON (`tblout.unclassified_hits`, `tblout.unclassified_models`), so a
  mismatched HMM set fails loudly rather than under-detecting anchors silently.

### Fixed
- Peak memory reduced ~6x on large inputs. `nhmmer`'s human-readable stdout
  report was being buffered in full and discarded; it scaled with reported hits
  at ~464 B/hit and accounted for 5.98 GB of a 7.07 GB peak RSS on a 230 Mbp
  input. Extraction results are unchanged.
- Empty input now reports "no sequence records" instead of an opaque `nhmmer`
  format-autodetection failure.

## [0.1.0] - 2026-02-13

### Added
- Initial release of ITSxRust.
- `extract` command for ITS region extraction using HMMER tblout anchors.
- Platform presets: `--preset ont` and `--preset hifi` with tuned E-values, length constraints, and confidence thresholds.
- Partial-chain fallback: recover subregions from 2-anchor pairs when the full 4-anchor chain is unavailable.
- Confidence classification: reads labelled `confident`, `ambiguous`, or `partial` based on per-anchor score/E-value thresholds.
- `--write-ambiguous` and `--write-skipped` flags to divert reads to separate files.
- Exact dereplication (`--derep`): hash identical sequences and search only unique representatives, projecting results back to duplicates.
- Multi-region output: `--region all` extracts ITS1, ITS2, and full ITS simultaneously, treating `--output` as a prefix.
- FASTQ support with quality score preservation; gzipped input support.
- `--qc-json` for per-sample QC summary (read counts, skip-reason breakdown, parameters).
- `--anchors-tsv` and `--anchors-jsonl` for per-read anchor coordinates and confidence labels.
- `--explain N` diagnostic flag for inspecting skip reasons.
- `--tblout-existing` to reuse a previous nhmmer tblout file.
- Docker image with bundled HMMER.
- Bioconda recipe.
- Deterministic test fixtures and CI checks (fmt, clippy, tests).