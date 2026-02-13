use crate::report::SkipReason;
use crate::tblout::TblHit;

#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub enum Anchor {
    SsuEnd,
    S58Start,
    S58End,
    LsuStart,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Region {
    Full,
    Its1,
    Its2,
    All,
}

#[derive(Copy, Clone, Debug)]
pub struct Interval {
    pub start: i64, // 1-based inclusive
    pub end: i64,   // 1-based inclusive
}

#[derive(Clone, Debug)]
pub struct HitIvl {
    pub anchor: Anchor,
    pub ivl: Interval, // normalized to forward coordinates
    pub score: f64,
    pub evalue: f64,
    pub strand: char, // original '+' or '-'
    pub model: String,
    pub seq_len: i64,
}

#[derive(Copy, Clone, Debug)]
pub struct Constraints {
    pub min_its1: i64,
    pub max_its1: i64,
    pub min_its2: i64,
    pub max_its2: i64,
    pub min_full: i64,
    pub max_full: i64,
}

impl Default for Constraints {
    fn default() -> Self {
        Self {
            min_its1: 50,
            max_its1: 1500,
            min_its2: 50,
            max_its2: 2000,
            min_full: 150,
            max_full: 4000,
        }
    }
}

fn in_range(x: i64, lo: i64, hi: i64) -> bool {
    x >= lo && x <= hi
}

fn anchor_name(a: Anchor) -> &'static str {
    match a {
        Anchor::SsuEnd => "SSU_end",
        Anchor::S58Start => "58S_start",
        Anchor::S58End => "58S_end",
        Anchor::LsuStart => "LSU_start",
    }
}

/// Normalize a hit interval so it is always in "forward" coordinates.
/// For '-' strand hits, convert to reverse-complement coordinate system:
/// If original covered [a,b], RC covered [L-b+1, L-a+1].
pub fn normalized_env_interval(hit: &TblHit) -> Interval {
    let a = hit.env_from.min(hit.env_to);
    let b = hit.env_from.max(hit.env_to);

    if hit.strand == '+' {
        Interval { start: a, end: b }
    } else {
        let l = hit.seq_len;
        Interval {
            start: l - b + 1,
            end: l - a + 1,
        }
    }
}

/// Classify model names from ITSx fungal HMMs (F.hmm).
pub fn classify(model: &str) -> Option<Anchor> {
    if model.contains("SSU_end") || model.starts_with("1_SSU") {
        Some(Anchor::SsuEnd)
    } else if model.contains("58S_start") || model.starts_with("2_5.8") {
        Some(Anchor::S58Start)
    } else if model.contains("58S_end") || model.starts_with("3_End_") {
        Some(Anchor::S58End)
    } else if model.contains("LSU_start") || model.starts_with("4_LSU_") {
        Some(Anchor::LsuStart)
    } else {
        None
    }
}

fn to_anchor_hits(hits: &[TblHit]) -> Vec<HitIvl> {
    let mut out = Vec::new();
    for h in hits {
        if let Some(a) = classify(&h.model) {
            out.push(HitIvl {
                anchor: a,
                ivl: normalized_env_interval(h),
                score: h.score,
                evalue: h.evalue,
                strand: h.strand,
                model: h.model.clone(),
                seq_len: h.seq_len,
            });
        }
    }
    out
}

/// Sort by score desc, then evalue asc, keep top K
fn top_k<'a>(mut v: Vec<&'a HitIvl>, k: usize) -> Vec<&'a HitIvl> {
    v.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                a.evalue
                    .partial_cmp(&b.evalue)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });
    v.truncate(k);
    v
}

fn best_chain_one_strand(
    hits: &[HitIvl],
    strand: char,
    c: Constraints,
    max_per_anchor: usize,
) -> Option<[HitIvl; 4]> {
    let ssu = top_k(
        hits.iter()
            .filter(|h| h.strand == strand && h.anchor == Anchor::SsuEnd)
            .collect(),
        max_per_anchor,
    );
    let s58s = top_k(
        hits.iter()
            .filter(|h| h.strand == strand && h.anchor == Anchor::S58Start)
            .collect(),
        max_per_anchor,
    );
    let s58e = top_k(
        hits.iter()
            .filter(|h| h.strand == strand && h.anchor == Anchor::S58End)
            .collect(),
        max_per_anchor,
    );
    let lsu = top_k(
        hits.iter()
            .filter(|h| h.strand == strand && h.anchor == Anchor::LsuStart)
            .collect(),
        max_per_anchor,
    );

    if ssu.is_empty() || s58s.is_empty() || s58e.is_empty() || lsu.is_empty() {
        return None;
    }

    let mut best_chain: Option<[HitIvl; 4]> = None;
    let mut best_score: f64 = f64::NEG_INFINITY;

    for a in &ssu {
        for b in &s58s {
            let its1_len = b.ivl.start - a.ivl.end - 1;
            if its1_len <= 0 || !in_range(its1_len, c.min_its1, c.max_its1) {
                continue;
            }
            for c58 in &s58e {
                if !(b.ivl.start < c58.ivl.start) {
                    continue;
                }
                for d in &lsu {
                    let its2_len = d.ivl.start - c58.ivl.end - 1;
                    if its2_len <= 0 || !in_range(its2_len, c.min_its2, c.max_its2) {
                        continue;
                    }

                    let full_len = d.ivl.start - a.ivl.end - 1;
                    if full_len <= 0 || !in_range(full_len, c.min_full, c.max_full) {
                        continue;
                    }

                    let chain_score = a.score + b.score + c58.score + d.score;
                    if chain_score > best_score {
                        best_score = chain_score;
                        best_chain =
                            Some([(*a).clone(), (*b).clone(), (*c58).clone(), (*d).clone()]);
                    }
                }
            }
        }
    }

    best_chain
}

fn best_chain(
    hits: &[TblHit],
    constraints: Constraints,
    max_per_anchor: usize,
) -> Option<[HitIvl; 4]> {
    let ah = to_anchor_hits(hits);

    let plus = best_chain_one_strand(&ah, '+', constraints, max_per_anchor);
    let minus = best_chain_one_strand(&ah, '-', constraints, max_per_anchor);

    match (plus, minus) {
        (Some(p), None) => Some(p),
        (None, Some(m)) => Some(m),
        (Some(p), Some(m)) => {
            let sp: f64 = p.iter().map(|h| h.score).sum();
            let sm: f64 = m.iter().map(|h| h.score).sum();
            if sp >= sm { Some(p) } else { Some(m) }
        }
        (None, None) => None,
    }
}

/// Public: return the chosen 4-anchor chain (SSU_end, 58S_start, 58S_end, LSU_start)
pub fn compute_chain(
    hits: &[TblHit],
    max_per_anchor: usize,
    constraints: Constraints,
) -> Option<[HitIvl; 4]> {
    best_chain(hits, constraints, max_per_anchor)
}

// ---------------------------------------------------------------------------
// Partial-chain (2-anchor) fallback
// ---------------------------------------------------------------------------

/// Per-region bounds from a 2-anchor pair, with the anchor hits that produced them.
#[derive(Clone, Debug)]
pub struct PairBound {
    pub start: i64,
    pub end: i64,
    pub left: HitIvl,
    pub right: HitIvl,
    pub strand: char,
    pub seq_len: i64,
}

/// Partial bounds: per-region results from 2-anchor pairs.
/// Used as fallback when the full 4-anchor chain fails.
#[derive(Clone, Debug, Default)]
pub struct PartialBounds {
    pub its1: Option<PairBound>,
    pub its2: Option<PairBound>,
    pub full: Option<PairBound>,
}

impl PartialBounds {
    pub fn has_any(&self) -> bool {
        self.its1.is_some() || self.its2.is_some() || self.full.is_some()
    }
}

/// Find the best 2-anchor pair for a given (left_anchor, right_anchor) on one strand.
fn best_pair_one_strand(
    hits: &[HitIvl],
    strand: char,
    left_anchor: Anchor,
    right_anchor: Anchor,
    min_len: i64,
    max_len: i64,
    max_per_anchor: usize,
) -> Option<PairBound> {
    let lefts = top_k(
        hits.iter()
            .filter(|h| h.strand == strand && h.anchor == left_anchor)
            .collect(),
        max_per_anchor,
    );
    let rights = top_k(
        hits.iter()
            .filter(|h| h.strand == strand && h.anchor == right_anchor)
            .collect(),
        max_per_anchor,
    );

    let mut best: Option<PairBound> = None;
    let mut best_score: f64 = f64::NEG_INFINITY;

    for l in &lefts {
        for r in &rights {
            let region_len = r.ivl.start - l.ivl.end - 1;
            if region_len <= 0 || !in_range(region_len, min_len, max_len) {
                continue;
            }
            let score = l.score + r.score;
            if score > best_score {
                best_score = score;
                best = Some(PairBound {
                    start: l.ivl.end + 1,
                    end: r.ivl.start - 1,
                    left: (*l).clone(),
                    right: (*r).clone(),
                    strand,
                    seq_len: l.seq_len,
                });
            }
        }
    }

    best
}

/// Pick the better PairBound between + and - strand.
fn best_pair_both_strands(
    hits: &[HitIvl],
    left_anchor: Anchor,
    right_anchor: Anchor,
    min_len: i64,
    max_len: i64,
    max_per_anchor: usize,
) -> Option<PairBound> {
    let plus = best_pair_one_strand(
        hits,
        '+',
        left_anchor,
        right_anchor,
        min_len,
        max_len,
        max_per_anchor,
    );
    let minus = best_pair_one_strand(
        hits,
        '-',
        left_anchor,
        right_anchor,
        min_len,
        max_len,
        max_per_anchor,
    );

    match (plus, minus) {
        (Some(p), None) => Some(p),
        (None, Some(m)) => Some(m),
        (Some(p), Some(m)) => {
            let sp = p.left.score + p.right.score;
            let sm = m.left.score + m.right.score;
            if sp >= sm { Some(p) } else { Some(m) }
        }
        (None, None) => None,
    }
}

/// Public: compute per-region bounds from 2-anchor pairs.
/// This is called as a fallback when the full 4-anchor chain fails.
pub fn compute_partial_bounds(
    hits: &[TblHit],
    max_per_anchor: usize,
    constraints: Constraints,
) -> PartialBounds {
    let ah = to_anchor_hits(hits);

    let its1 = best_pair_both_strands(
        &ah,
        Anchor::SsuEnd,
        Anchor::S58Start,
        constraints.min_its1,
        constraints.max_its1,
        max_per_anchor,
    );

    let its2 = best_pair_both_strands(
        &ah,
        Anchor::S58End,
        Anchor::LsuStart,
        constraints.min_its2,
        constraints.max_its2,
        max_per_anchor,
    );

    let full = best_pair_both_strands(
        &ah,
        Anchor::SsuEnd,
        Anchor::LsuStart,
        constraints.min_full,
        constraints.max_full,
        max_per_anchor,
    );

    PartialBounds { its1, its2, full }
}

/// Public: compute bounds from a chain for a requested region.
/// Returns (start,end) 1-based inclusive in forward coords.
pub fn bounds_from_chain(chain: &[HitIvl; 4], region: Region) -> Option<(i64, i64)> {
    let ssu = &chain[0];
    let s58s = &chain[1];
    let s58e = &chain[2];
    let lsu = &chain[3];

    match region {
        Region::Its1 => {
            let start = ssu.ivl.end + 1;
            let end = s58s.ivl.start - 1;
            if start <= end {
                Some((start, end))
            } else {
                None
            }
        }
        Region::Its2 => {
            let start = s58e.ivl.end + 1;
            let end = lsu.ivl.start - 1;
            if start <= end {
                Some((start, end))
            } else {
                None
            }
        }
        Region::Full => {
            let start = ssu.ivl.end + 1;
            let end = lsu.ivl.start - 1;
            if start <= end {
                Some((start, end))
            } else {
                None
            }
        }
        Region::All => None,
    }
}

#[derive(Default)]
struct AnchorCounts {
    ssu_p: usize,
    ssu_m: usize,
    s58s_p: usize,
    s58s_m: usize,
    s58e_p: usize,
    s58e_m: usize,
    lsu_p: usize,
    lsu_m: usize,
}

fn count_anchors(ah: &[HitIvl]) -> AnchorCounts {
    let mut c = AnchorCounts::default();
    for h in ah {
        match (h.anchor, h.strand) {
            (Anchor::SsuEnd, '+') => c.ssu_p += 1,
            (Anchor::SsuEnd, '-') => c.ssu_m += 1,
            (Anchor::S58Start, '+') => c.s58s_p += 1,
            (Anchor::S58Start, '-') => c.s58s_m += 1,
            (Anchor::S58End, '+') => c.s58e_p += 1,
            (Anchor::S58End, '-') => c.s58e_m += 1,
            (Anchor::LsuStart, '+') => c.lsu_p += 1,
            (Anchor::LsuStart, '-') => c.lsu_m += 1,
            _ => {}
        }
    }
    c
}

/// Public: explain why a read likely failed to yield bounds.
pub fn diagnose(
    hits: &[TblHit],
    region: Region,
    max_per_anchor: usize,
    constraints: Constraints,
) -> String {
    let ah = to_anchor_hits(hits);
    if ah.is_empty() {
        return "no classified anchor hits (none of SSU_end/58S/LSU matched)".to_string();
    }

    let counts = count_anchors(&ah);

    let missing = {
        let mut v = Vec::new();
        let has_ssu = counts.ssu_p + counts.ssu_m > 0;
        let has_s58s = counts.s58s_p + counts.s58s_m > 0;
        let has_s58e = counts.s58e_p + counts.s58e_m > 0;
        let has_lsu = counts.lsu_p + counts.lsu_m > 0;

        if !has_ssu {
            v.push(anchor_name(Anchor::SsuEnd));
        }
        if !has_s58s {
            v.push(anchor_name(Anchor::S58Start));
        }
        if !has_s58e {
            v.push(anchor_name(Anchor::S58End));
        }
        if !has_lsu {
            v.push(anchor_name(Anchor::LsuStart));
        }
        v
    };

    let counts_str = format!(
        "counts(+/-): SSU_end {}/{}  58S_start {}/{}  58S_end {}/{}  LSU_start {}/{}",
        counts.ssu_p,
        counts.ssu_m,
        counts.s58s_p,
        counts.s58s_m,
        counts.s58e_p,
        counts.s58e_m,
        counts.lsu_p,
        counts.lsu_m
    );

    if !missing.is_empty() {
        return format!("missing anchors: {} | {}", missing.join(","), counts_str);
    }

    // anchors exist; do we have a valid chain under constraints?
    let chain = compute_chain(hits, max_per_anchor, constraints);
    if chain.is_none() {
        return format!(
            "anchors present but no valid SSU_endâ†’58S_startâ†’58S_endâ†’LSU_start chain under constraints | {}",
            counts_str
        );
    }

    let chain = chain.unwrap();
    if bounds_from_chain(&chain, region).is_none() {
        return format!(
            "chain found, but requested region {:?} bounds invalid (start>end) | {}",
            region, counts_str
        );
    }

    // If we reach here, bounds *should* exist; this typically means downstream trimming failed (seq shorter than bounds).
    format!(
        "bounds exist logically but trimming failed (likely bounds exceed sequence length) | {}",
        counts_str
    )
}

/// Structured version of `diagnose` returning a machine-readable `SkipReason`.
pub fn diagnose_structured(
    hits: &[TblHit],
    region: Region,
    max_per_anchor: usize,
    constraints: Constraints,
) -> SkipReason {
    let ah = to_anchor_hits(hits);
    if ah.is_empty() {
        return SkipReason::NoAnchorHits;
    }

    let counts = count_anchors(&ah);
    let has_ssu = counts.ssu_p + counts.ssu_m > 0;
    let has_s58s = counts.s58s_p + counts.s58s_m > 0;
    let has_s58e = counts.s58e_p + counts.s58e_m > 0;
    let has_lsu = counts.lsu_p + counts.lsu_m > 0;

    let mut missing = Vec::new();
    if !has_ssu {
        missing.push(anchor_name(Anchor::SsuEnd).to_string());
    }
    if !has_s58s {
        missing.push(anchor_name(Anchor::S58Start).to_string());
    }
    if !has_s58e {
        missing.push(anchor_name(Anchor::S58End).to_string());
    }
    if !has_lsu {
        missing.push(anchor_name(Anchor::LsuStart).to_string());
    }

    if !missing.is_empty() {
        return SkipReason::MissingAnchors { missing };
    }

    let chain = compute_chain(hits, max_per_anchor, constraints);
    if chain.is_none() {
        return SkipReason::NoValidChain;
    }

    let chain = chain.unwrap();
    if bounds_from_chain(&chain, region).is_none() {
        return SkipReason::InvalidBounds;
    }

    SkipReason::TrimFailed
}
