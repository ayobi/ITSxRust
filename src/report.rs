// src/report.rs
//
// Structured failure codes and per-sample QC summary.
//
// SkipReason is a machine-readable enum that replaces the old string-based
// diagnose() output. QcSummary aggregates all counters into a single JSON
// object suitable for MultiQC custom_data ingestion.

use serde::Serialize;
use std::collections::HashMap;
use std::fmt;

// ---------------------------------------------------------------------------
// Skip reason codes
// ---------------------------------------------------------------------------

/// Why a read was not extracted.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum SkipReason {
    /// No rows for this read in the tblout at all.
    NoHmmHits,
    /// Tblout rows exist but none match an anchor model (SSU_end/58S/LSU).
    NoAnchorHits,
    /// Some anchor types are missing entirely.
    MissingAnchors { missing: Vec<String> },
    /// All four anchors present but no valid chain under length constraints.
    NoValidChain,
    /// Chain found but the requested region yields invalid bounds (start > end).
    InvalidBounds,
    /// Bounds are valid logically but exceed the actual sequence length.
    TrimFailed,
}

impl SkipReason {
    /// Short machine-readable code for JSON/TSV output.
    pub fn code(&self) -> &'static str {
        match self {
            SkipReason::NoHmmHits => "no_hmm_hits",
            SkipReason::NoAnchorHits => "no_anchor_hits",
            SkipReason::MissingAnchors { .. } => "missing_anchors",
            SkipReason::NoValidChain => "no_valid_chain",
            SkipReason::InvalidBounds => "invalid_bounds",
            SkipReason::TrimFailed => "trim_failed",
        }
    }
}

impl fmt::Display for SkipReason {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SkipReason::NoHmmHits => write!(f, "no HMM hits"),
            SkipReason::NoAnchorHits => {
                write!(
                    f,
                    "no classified anchor hits (none of SSU_end/58S/LSU matched)"
                )
            }
            SkipReason::MissingAnchors { missing } => {
                write!(f, "missing anchors: {}", missing.join(","))
            }
            SkipReason::NoValidChain => {
                write!(f, "anchors present but no valid chain under constraints")
            }
            SkipReason::InvalidBounds => {
                write!(
                    f,
                    "chain found but requested region bounds invalid (start>end)"
                )
            }
            SkipReason::TrimFailed => {
                write!(
                    f,
                    "bounds valid but trimming failed (bounds exceed sequence length)"
                )
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Skip reason counter
// ---------------------------------------------------------------------------

/// Accumulates skip reason counts during trimming.
#[derive(Debug, Default)]
pub struct SkipCounter {
    counts: HashMap<String, usize>,
}

impl SkipCounter {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn record(&mut self, reason: &SkipReason) {
        *self.counts.entry(reason.code().to_string()).or_insert(0) += 1;
    }

    pub fn to_map(&self) -> HashMap<String, usize> {
        self.counts.clone()
    }

    pub fn total(&self) -> usize {
        self.counts.values().sum()
    }
}

// ---------------------------------------------------------------------------
// QC summary
// ---------------------------------------------------------------------------

#[derive(Serialize)]
pub struct QcSummary {
    pub total_reads: usize,

    pub kept: KeptSummary,
    pub skipped: SkippedSummary,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub regions: Option<RegionSummary>,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub derep: Option<DerepSummary>,

    pub tblout: TbloutSummary,
    pub params: ParamsSummary,
}

#[derive(Serialize)]
pub struct KeptSummary {
    pub total: usize,
    pub confident: usize,
    pub ambiguous: usize,
    pub ambiguous_diverted: usize,
    pub partial: usize,
}

#[derive(Serialize)]
pub struct SkippedSummary {
    pub total: usize,
    /// Breakdown by reason code.
    pub by_reason: HashMap<String, usize>,
}

#[derive(Serialize)]
pub struct RegionSummary {
    pub full: usize,
    pub its1: usize,
    pub its2: usize,
}

#[derive(Serialize)]
pub struct DerepSummary {
    pub total_seqs: usize,
    pub unique_seqs: usize,
    pub reduction_pct: f64,
}

#[derive(Serialize)]
pub struct TbloutSummary {
    pub total_hits: u64,
    pub anchor_hits: u64,
    pub stored_topk: u64,
    pub reads_with_hits: usize,
}

#[derive(Serialize)]
pub struct ParamsSummary {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub preset: Option<String>,
    pub region: String,
    pub inc_e: f64,
    pub max_per_anchor: usize,
    pub min_anchor_score: f64,
    pub max_anchor_evalue: f64,
    pub min_its1: i64,
    pub max_its1: i64,
    pub min_its2: i64,
    pub max_its2: i64,
    pub min_full: i64,
    pub max_full: i64,
}

impl QcSummary {
    /// Write the summary as pretty-printed JSON to the given path.
    pub fn write_json(&self, path: &std::path::Path) -> anyhow::Result<()> {
        let f = std::fs::File::create(path)?;
        let w = std::io::BufWriter::new(f);
        serde_json::to_writer_pretty(w, self)?;
        Ok(())
    }
}
