// src/preset.rs
//
// Platform presets and confidence classification for extracted reads.
//
// Presets bundle sensible defaults for ONT and HiFi long-read platforms.
// Users can override any individual parameter after applying a preset.
//
// Confidence classification assigns each chain one of:
//   - "confident"  : all anchors meet score/E-value thresholds
//   - "ambiguous"  : chain found but at least one anchor is weak
// Reads with no valid chain are simply "no_chain" (handled elsewhere).

use crate::select::{Constraints, HitIvl};

// ---------------------------------------------------------------------------
// Preset definitions
// ---------------------------------------------------------------------------

/// All tunable parameters that a preset can set.
/// `None` means "use the CLI default / don't override".
#[derive(Debug, Clone)]
pub struct PresetParams {
    pub inc_e: Option<f64>,
    pub max_per_anchor: Option<usize>,
    pub constraints: Option<Constraints>,
    /// Per-anchor minimum bitscore to be considered "confident"
    pub min_anchor_score: Option<f64>,
    /// Per-anchor maximum E-value to be considered "confident"
    pub max_anchor_evalue: Option<f64>,
}

/// ONT preset: tolerant of higher error rates and partial reads.
pub fn ont_preset() -> PresetParams {
    PresetParams {
        inc_e: Some(1e-3),
        max_per_anchor: Some(10),
        constraints: Some(Constraints {
            min_its1: 30,
            max_its1: 1800,
            min_its2: 30,
            max_its2: 2500,
            min_full: 100,
            max_full: 5000,
        }),
        min_anchor_score: Some(15.0),
        max_anchor_evalue: Some(1e-3),
    }
}

/// HiFi preset: stricter thresholds for high-accuracy reads.
pub fn hifi_preset() -> PresetParams {
    PresetParams {
        inc_e: Some(1e-10),
        max_per_anchor: Some(6),
        constraints: Some(Constraints {
            min_its1: 50,
            max_its1: 1500,
            min_its2: 50,
            max_its2: 2000,
            min_full: 150,
            max_full: 4000,
        }),
        min_anchor_score: Some(30.0),
        max_anchor_evalue: Some(1e-8),
    }
}

/// Default confidence thresholds when no preset is used.
pub fn default_confidence() -> (f64, f64) {
    // (min_anchor_score, max_anchor_evalue)
    (20.0, 1e-4)
}

// ---------------------------------------------------------------------------
// Confidence classification
// ---------------------------------------------------------------------------

/// Confidence level for a chain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Confidence {
    Confident,
    Ambiguous,
    /// From a 2-anchor pair fallback (no full 4-anchor chain available).
    Partial,
}

impl Confidence {
    pub fn as_str(&self) -> &'static str {
        match self {
            Confidence::Confident => "confident",
            Confidence::Ambiguous => "ambiguous",
            Confidence::Partial => "partial",
        }
    }
}

impl std::fmt::Display for Confidence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Classify a 4-anchor chain as confident or ambiguous.
///
/// A chain is "confident" when every anchor meets both thresholds:
///   - bitscore >= min_anchor_score
///   - E-value  <= max_anchor_evalue
///
/// Otherwise the chain is "ambiguous".
pub fn classify_chain(
    chain: &[HitIvl; 4],
    min_anchor_score: f64,
    max_anchor_evalue: f64,
) -> Confidence {
    for h in chain {
        if h.score < min_anchor_score || h.evalue > max_anchor_evalue {
            return Confidence::Ambiguous;
        }
    }
    Confidence::Confident
}

/// Return a short human-readable reason string for an ambiguous chain.
pub fn ambiguous_reason(
    chain: &[HitIvl; 4],
    min_anchor_score: f64,
    max_anchor_evalue: f64,
) -> String {
    let labels = ["SSU_end", "58S_start", "58S_end", "LSU_start"];
    let mut weak: Vec<String> = Vec::new();

    for (h, label) in chain.iter().zip(labels.iter()) {
        let mut issues = Vec::new();
        if h.score < min_anchor_score {
            issues.push(format!("score={:.1}<{:.1}", h.score, min_anchor_score));
        }
        if h.evalue > max_anchor_evalue {
            issues.push(format!("evalue={:.1e}>{:.1e}", h.evalue, max_anchor_evalue));
        }
        if !issues.is_empty() {
            weak.push(format!("{}({})", label, issues.join(",")));
        }
    }

    if weak.is_empty() {
        "none".to_string()
    } else {
        weak.join("; ")
    }
}
