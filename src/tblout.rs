// src/tblout.rs

use anyhow::{Context, Result};
use std::cmp::Ordering;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::select::{self, Anchor};

#[derive(Debug, Clone)]
pub struct TblHit {
    pub read_id: String,
    pub model: String,

    pub env_from: i64,
    pub env_to: i64,
    pub seq_len: i64,
    pub strand: char,

    pub evalue: f64,
    pub score: f64,
}

fn parse_tblout_line(line: &str) -> Option<TblHit> {
    if line.is_empty() || line.starts_with('#') {
        return None;
    }
    let f: Vec<&str> = line.split_whitespace().collect();
    if f.len() < 15 {
        return None;
    }

    // nhmmer --tblout columns (typical):
    // 0 target name
    // 1 tacc
    // 2 query name
    // 3 qacc
    // 8 envfrom
    // 9 envto
    // 10 sqlen
    // 11 strand
    // 12 E-value
    // 13 score
    let read_id = f[0].to_string();
    let model = f[2].to_string();

    let env_from = f[8].parse().ok()?;
    let env_to = f[9].parse().ok()?;
    let seq_len = f[10].parse().ok()?;

    let strand_s = f[11];
    let strand = strand_s.chars().next()?;
    if strand != '+' && strand != '-' {
        return None;
    }

    let evalue = f[12].parse().ok()?;
    let score = f[13].parse().ok()?;

    Some(TblHit {
        read_id,
        model,
        env_from,
        env_to,
        seq_len,
        strand,
        evalue,
        score,
    })
}

/// 8 bins: 4 anchors × 2 strands
/// idx = anchor_index * 2 + strand_index, where strand_index: '+' -> 0, '-' -> 1
fn bin_index(anchor: Anchor, strand: char) -> usize {
    let a = match anchor {
        Anchor::SsuEnd => 0,
        Anchor::S58Start => 1,
        Anchor::S58End => 2,
        Anchor::LsuStart => 3,
    };
    let s = if strand == '+' { 0 } else { 1 };
    a * 2 + s
}

/// Sort: score desc, then evalue asc
fn cmp_hits(a: &TblHit, b: &TblHit) -> Ordering {
    b.score
        .partial_cmp(&a.score)
        .unwrap_or(Ordering::Equal)
        .then_with(|| a.evalue.partial_cmp(&b.evalue).unwrap_or(Ordering::Equal))
}

fn push_topk(vec: &mut Vec<TblHit>, hit: TblHit, k: usize) {
    vec.push(hit);
    vec.sort_by(cmp_hits);
    if vec.len() > k {
        vec.truncate(k);
    }
}

#[derive(Debug, Default, Clone)]
pub struct TopKHits {
    pub bins: [Vec<TblHit>; 8],
}

impl TopKHits {
    /// Compatibility helper: gather all stored hits into one small Vec
    /// (max size = 8 * K).
    pub fn flatten(&self) -> Vec<TblHit> {
        let mut out = Vec::new();
        for b in &self.bins {
            out.extend(b.iter().cloned());
        }
        out
    }

    pub fn stored_len(&self) -> usize {
        self.bins.iter().map(|b| b.len()).sum()
    }
}

#[derive(Debug, Clone)]
pub struct StreamStats {
    pub total_tblout_hits: u64, // non-# lines parsed into hits
    pub anchor_hits: u64,       // hits that matched an anchor + strand
    pub stored_hits: u64,       // final stored hits across all reads/bins
    pub reads: usize,           // reads with >=1 anchor hit
}

/// Stream tblout line-by-line and keep only Top-K per (anchor × strand) per read.
///
/// - `k`: typically your `--max-per-anchor`
/// - `evalue_cutoff`: optional extra guard (use Some(inc_e) if you want)
pub fn stream_topk_tblout(
    tblout_path: &Path,
    k: usize,
    evalue_cutoff: Option<f64>,
) -> Result<(HashMap<String, TopKHits>, StreamStats)> {
    let f = File::open(tblout_path)
        .with_context(|| format!("Failed to open tblout: {}", tblout_path.display()))?;
    let mut r = BufReader::with_capacity(1 << 20, f);

    let mut map: HashMap<String, TopKHits> = HashMap::new();
    let mut total_tblout_hits: u64 = 0;
    let mut anchor_hits: u64 = 0;

    let mut line = String::new();
    loop {
        line.clear();
        let n = r.read_line(&mut line)?;
        if n == 0 {
            break;
        }

        let Some(hit) = parse_tblout_line(line.trim_end()) else {
            continue;
        };
        total_tblout_hits += 1;

        if let Some(cut) = evalue_cutoff
            && hit.evalue > cut
        {
            continue;
        }

        let Some(anchor) = select::classify(&hit.model) else {
            continue;
        };

        anchor_hits += 1;
        let rid = hit.read_id.clone();
        let idx = bin_index(anchor, hit.strand);

        let entry = map.entry(rid).or_default();

        push_topk(&mut entry.bins[idx], hit, k);
    }

    let stored_hits: u64 = map.values().map(|t| t.stored_len() as u64).sum();
    let stats = StreamStats {
        total_tblout_hits,
        anchor_hits,
        stored_hits,
        reads: map.len(),
    };

    Ok((map, stats))
}
