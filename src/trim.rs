// src/trim.rs
//
// Unified record representation and trimming/writing logic.

use anyhow::{Result, anyhow};
use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::fasta;
use crate::fastq;
use crate::report::{SkipCounter, SkipReason};
use crate::select;
use crate::seq;
use crate::tblout::TopKHits;
use crate::{ExtractPlan, InputFormat, OutputFormat, PlanSet};

// ---------------------------------------------------------------------------
// Unified record: FASTA or FASTQ
// ---------------------------------------------------------------------------

pub struct Record {
    pub id: String,
    pub seq: String,
    pub qual: Option<String>,
}

// ---------------------------------------------------------------------------
// Unified reader
// ---------------------------------------------------------------------------

pub enum SeqReader {
    Fasta(fasta::FastaReader),
    Fastq(fastq::FastqReader),
}

impl SeqReader {
    pub fn open(path: &Path, fmt: InputFormat) -> Result<Self> {
        match fmt {
            InputFormat::Fasta => Ok(SeqReader::Fasta(fasta::FastaReader::new(path)?)),
            InputFormat::Fastq => Ok(SeqReader::Fastq(fastq::FastqReader::new(path)?)),
            InputFormat::Auto => Err(anyhow!("SeqReader::open called with Auto format")),
        }
    }

    pub fn next_record(&mut self) -> Result<Option<Record>> {
        match self {
            SeqReader::Fasta(r) => {
                if let Some(rec) = r.next_record()? {
                    Ok(Some(Record {
                        id: rec.id,
                        seq: rec.seq,
                        qual: None,
                    }))
                } else {
                    Ok(None)
                }
            }
            SeqReader::Fastq(r) => {
                if let Some(rec) = r.next_record()? {
                    Ok(Some(Record {
                        id: rec.id,
                        seq: rec.seq,
                        qual: Some(rec.qual),
                    }))
                } else {
                    Ok(None)
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Unified writer
// ---------------------------------------------------------------------------

pub enum SeqWriter {
    Fasta(BufWriter<std::fs::File>),
    Fastq(BufWriter<std::fs::File>),
}

impl SeqWriter {
    pub fn open(path: &Path, fmt: OutputFormat) -> Result<Self> {
        match fmt {
            OutputFormat::Fasta => Ok(SeqWriter::Fasta(fasta::open_writer(path)?)),
            OutputFormat::Fastq => Ok(SeqWriter::Fastq(fastq::open_writer(path)?)),
            OutputFormat::Auto => Err(anyhow!("SeqWriter::open called with Auto format")),
        }
    }

    pub fn write_record(&mut self, id: &str, seq_data: &str, qual: Option<&str>) -> Result<()> {
        match self {
            SeqWriter::Fasta(w) => fasta::write_fasta_record(w, id, seq_data),
            SeqWriter::Fastq(w) => {
                let q = qual.ok_or_else(|| {
                    anyhow!("FASTQ output requires quality scores (cannot convert FASTA to FASTQ)")
                })?;
                fastq::write_fastq_record(w, id, seq_data, q)
            }
        }
    }

    pub fn flush(&mut self) -> Result<()> {
        match self {
            SeqWriter::Fasta(w) => Ok(w.flush()?),
            SeqWriter::Fastq(w) => Ok(w.flush()?),
        }
    }
}

// ---------------------------------------------------------------------------
// Core extraction
// ---------------------------------------------------------------------------

fn extract_trimmed(rec: &Record, plan: &ExtractPlan) -> Option<(String, Option<String>)> {
    let s0 = (plan.orig_start - 1) as usize;
    let e0 = plan.orig_end as usize;

    if e0 > rec.seq.len() || s0 >= e0 {
        return None;
    }

    let seq_slice = &rec.seq[s0..e0];

    if plan.strand == '-' {
        let rc_seq = seq::revcomp_dna(seq_slice);
        let rc_qual = rec.qual.as_ref().map(|q| seq::reverse_qual(&q[s0..e0]));
        Some((rc_seq, rc_qual))
    } else {
        let out_seq = seq_slice.to_string();
        let out_qual = rec.qual.as_ref().map(|q| q[s0..e0].to_string());
        Some((out_seq, out_qual))
    }
}

fn trimmed_id(read_id: &str, region_tag: &str, plan: &ExtractPlan) -> String {
    format!(
        "{}|{}:{}-{}",
        read_id, region_tag, plan.norm_start, plan.norm_end
    )
}

// ---------------------------------------------------------------------------
// Counters
// ---------------------------------------------------------------------------

#[derive(Default)]
pub struct TrimStats {
    pub total_reads: usize,
    pub kept_any: usize,
    pub kept_full: usize,
    pub kept_its1: usize,
    pub kept_its2: usize,
    pub ambiguous_out: usize,
    pub skipped: usize,
    pub skip_reasons: SkipCounter,
}

// ---------------------------------------------------------------------------
// Region output helper
// ---------------------------------------------------------------------------

struct RegionOutput {
    tag: &'static str,
    writer: SeqWriter,
}

// ---------------------------------------------------------------------------
// Skip handling
// ---------------------------------------------------------------------------

fn handle_skip(
    rec: &Record,
    topk_map: &HashMap<String, TopKHits>,
    stats: &mut TrimStats,
    skipped_writer: &mut Option<SeqWriter>,
    diag_region: select::Region,
    max_per_anchor: usize,
    constraints: select::Constraints,
    explained: &mut usize,
    explain: usize,
) -> Result<()> {
    stats.skipped += 1;

    // Structured reason for QC summary
    let reason = if let Some(tk) = topk_map.get(&rec.id) {
        let hs = tk.flatten();
        select::diagnose_structured(&hs, diag_region, max_per_anchor, constraints)
    } else {
        SkipReason::NoHmmHits
    };
    stats.skip_reasons.record(&reason);

    if let Some(sw) = skipped_writer.as_mut() {
        sw.write_record(&rec.id, &rec.seq, rec.qual.as_deref())?;
    }

    if *explained < explain {
        // Use the Display impl for human-readable stderr output
        eprintln!("SKIP {}: {}", rec.id, reason);
        *explained += 1;
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Public: unified trimming loop — all regions
// ---------------------------------------------------------------------------

#[allow(clippy::too_many_arguments)]
pub fn trim_all(
    reader: &mut SeqReader,
    out_fmt: OutputFormat,
    out_full_path: &Path,
    out_its1_path: &Path,
    out_its2_path: &Path,
    bounds_map: &HashMap<String, PlanSet>,
    topk_map: &HashMap<String, TopKHits>,
    skipped_writer: &mut Option<SeqWriter>,
    ambiguous_writer: &mut Option<SeqWriter>,
    ambiguous_ids: Option<&HashMap<String, bool>>,
    diag_region: select::Region,
    max_per_anchor: usize,
    constraints: select::Constraints,
    explain: usize,
) -> Result<TrimStats> {
    let mut outputs = [
        RegionOutput {
            tag: "full",
            writer: SeqWriter::open(out_full_path, out_fmt)?,
        },
        RegionOutput {
            tag: "its1",
            writer: SeqWriter::open(out_its1_path, out_fmt)?,
        },
        RegionOutput {
            tag: "its2",
            writer: SeqWriter::open(out_its2_path, out_fmt)?,
        },
    ];

    let mut stats = TrimStats::default();
    let mut explained = 0usize;

    while let Some(rec) = reader.next_record()? {
        stats.total_reads += 1;

        if let Some(plans) = bounds_map.get(&rec.id) {
            let is_ambiguous = ambiguous_ids
                .map(|ids| ids.contains_key(&rec.id))
                .unwrap_or(false);

            if is_ambiguous {
                if let Some(aw) = ambiguous_writer.as_mut() {
                    aw.write_record(&rec.id, &rec.seq, rec.qual.as_deref())?;
                }
                stats.ambiguous_out += 1;
                continue;
            }

            let plan_array: [Option<ExtractPlan>; 3] = [plans.full, plans.its1, plans.its2];
            let mut wrote = false;

            for (i, maybe_plan) in plan_array.iter().enumerate() {
                if let Some(plan) = maybe_plan {
                    if let Some((out_seq, out_qual)) = extract_trimmed(&rec, plan) {
                        let out_id = trimmed_id(&rec.id, outputs[i].tag, plan);
                        outputs[i]
                            .writer
                            .write_record(&out_id, &out_seq, out_qual.as_deref())?;
                        wrote = true;
                        match i {
                            0 => stats.kept_full += 1,
                            1 => stats.kept_its1 += 1,
                            2 => stats.kept_its2 += 1,
                            _ => unreachable!(),
                        }
                    }
                }
            }

            if wrote {
                stats.kept_any += 1;
                continue;
            }
        }

        handle_skip(
            &rec,
            topk_map,
            &mut stats,
            skipped_writer,
            diag_region,
            max_per_anchor,
            constraints,
            &mut explained,
            explain,
        )?;
    }

    for o in &mut outputs {
        o.writer.flush()?;
    }

    Ok(stats)
}

// ---------------------------------------------------------------------------
// Public: unified trimming loop — single region
// ---------------------------------------------------------------------------

#[allow(clippy::too_many_arguments)]
pub fn trim_single(
    reader: &mut SeqReader,
    out_fmt: OutputFormat,
    out_path: &Path,
    region_tag: &str,
    bounds_map: &HashMap<String, ExtractPlan>,
    topk_map: &HashMap<String, TopKHits>,
    skipped_writer: &mut Option<SeqWriter>,
    ambiguous_writer: &mut Option<SeqWriter>,
    ambiguous_ids: Option<&HashMap<String, bool>>,
    diag_region: select::Region,
    max_per_anchor: usize,
    constraints: select::Constraints,
    explain: usize,
) -> Result<TrimStats> {
    let mut w = SeqWriter::open(out_path, out_fmt)?;
    let mut stats = TrimStats::default();
    let mut explained = 0usize;

    while let Some(rec) = reader.next_record()? {
        stats.total_reads += 1;

        if let Some(plan) = bounds_map.get(&rec.id) {
            let is_ambiguous = ambiguous_ids
                .map(|ids| ids.contains_key(&rec.id))
                .unwrap_or(false);

            if is_ambiguous {
                if let Some(aw) = ambiguous_writer.as_mut() {
                    aw.write_record(&rec.id, &rec.seq, rec.qual.as_deref())?;
                }
                stats.ambiguous_out += 1;
                continue;
            }

            if let Some((out_seq, out_qual)) = extract_trimmed(&rec, plan) {
                let out_id = trimmed_id(&rec.id, region_tag, plan);
                w.write_record(&out_id, &out_seq, out_qual.as_deref())?;
                stats.kept_any += 1;
                continue;
            }
        }

        handle_skip(
            &rec,
            topk_map,
            &mut stats,
            skipped_writer,
            diag_region,
            max_per_anchor,
            constraints,
            &mut explained,
            explain,
        )?;
    }

    w.flush()?;
    Ok(stats)
}
