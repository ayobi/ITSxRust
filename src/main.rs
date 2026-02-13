// src/main.rs

mod derep;
mod fasta;
mod fastq;
mod hmmer;
mod preset;
mod report;
mod select;
mod seq;
mod tblout;
mod trim;

use anyhow::{Result, anyhow};
use clap::{Parser, Subcommand, ValueEnum};
use rayon::prelude::*;
use serde::Serialize;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use tempfile::tempdir;

#[derive(Copy, Clone, Debug, ValueEnum)]
pub enum InputFormat {
    Auto,
    Fasta,
    Fastq,
}

#[derive(Copy, Clone, Debug, ValueEnum)]
pub enum OutputFormat {
    Auto,
    Fasta,
    Fastq,
}

#[derive(Copy, Clone, Debug, ValueEnum)]
pub enum Preset {
    Ont,
    Hifi,
}

#[derive(Parser, Debug)]
#[command(name = "itsxrust", version, about = "ONT ITS region extractor (v1)")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Run nhmmer -> tblout -> compute bounds -> write trimmed output (FASTA/FASTQ)
    Extract {
        /// Input sequence file (FASTA/FASTQ; optionally .gz)
        #[arg(short, long)]
        input: PathBuf,

        /// HMM profile file (SSU/5.8S/LSU)
        #[arg(long)]
        hmm: Option<PathBuf>,

        /// Reuse an existing nhmmer --tblout output (skips running nhmmer)
        #[arg(long, conflicts_with = "tblout")]
        tblout_existing: Option<PathBuf>,

        /// Output file path (single region) OR output prefix when --region all
        #[arg(short, long)]
        output: PathBuf,

        /// Input format (auto detects by extension)
        #[arg(long, default_value = "auto")]
        input_format: InputFormat,

        /// Output format (default auto: FASTQ if input is FASTQ else FASTA)
        #[arg(long, default_value = "auto")]
        output_format: OutputFormat,

        /// Output tblout path (optional; if omitted, a temp file is used)
        #[arg(long)]
        tblout: Option<PathBuf>,

        /// Platform preset: ont (tolerant) or hifi (strict). Explicit flags override preset values.
        #[arg(long)]
        preset: Option<Preset>,

        /// HMMER threads
        #[arg(long, default_value_t = 8)]
        hmmer_cpu: usize,

        /// Inclusion E-value threshold (preset default: ont=1e-3, hifi=1e-10, none=1e-5)
        #[arg(long)]
        inc_e: Option<f64>,

        /// Region to extract: full (default), its1, its2, all
        #[arg(long, default_value = "full")]
        region: String,

        /// Limit top hits per anchor per strand used to form chains
        #[arg(long)]
        max_per_anchor: Option<usize>,

        /// Constraints (bp) — preset provides defaults; explicit flags override
        #[arg(long)]
        min_its1: Option<i64>,
        #[arg(long)]
        max_its1: Option<i64>,
        #[arg(long)]
        min_its2: Option<i64>,
        #[arg(long)]
        max_its2: Option<i64>,
        #[arg(long)]
        min_full: Option<i64>,
        #[arg(long)]
        max_full: Option<i64>,

        /// Minimum per-anchor bitscore for "confident" classification (default: 20)
        #[arg(long)]
        min_anchor_score: Option<f64>,

        /// Maximum per-anchor E-value for "confident" classification (default: 1e-4)
        #[arg(long)]
        max_anchor_evalue: Option<f64>,

        /// Write chosen anchors per read as TSV (optional)
        #[arg(long)]
        anchors_tsv: Option<PathBuf>,

        /// Write chosen anchors per read as JSONL (optional)
        #[arg(long)]
        anchors_jsonl: Option<PathBuf>,

        /// Print skip reasons for first N skipped reads
        #[arg(long, default_value_t = 0)]
        explain: usize,

        /// Optionally write skipped reads to this file (same format as --output-format)
        #[arg(long)]
        write_skipped: Option<PathBuf>,

        /// Write ambiguous reads to a separate file (excludes them from main output)
        #[arg(long)]
        write_ambiguous: Option<PathBuf>,

        /// Write per-sample QC summary as JSON (for MultiQC or downstream aggregation)
        #[arg(long)]
        qc_json: Option<PathBuf>,

        /// Exact dereplication: search only unique sequences, project results to duplicates
        #[arg(long, default_value_t = false)]
        derep: bool,
    },
}

fn detect_format(path: &Path) -> Option<InputFormat> {
    let name = path.file_name()?.to_string_lossy().to_lowercase();
    let base = if name.ends_with(".gz") {
        name.trim_end_matches(".gz").to_string()
    } else {
        name
    };
    if base.ends_with(".fa") || base.ends_with(".fasta") || base.ends_with(".fna") {
        Some(InputFormat::Fasta)
    } else if base.ends_with(".fq") || base.ends_with(".fastq") {
        Some(InputFormat::Fastq)
    } else {
        None
    }
}

fn out_ext(out_fmt: OutputFormat) -> &'static str {
    match out_fmt {
        OutputFormat::Fasta => "fasta",
        OutputFormat::Fastq => "fastq",
        OutputFormat::Auto => unreachable!(),
    }
}

fn out_path(prefix: &Path, tag: &str, ext: &str) -> PathBuf {
    PathBuf::from(format!("{}.{}.{}", prefix.to_string_lossy(), tag, ext))
}

fn open_opt_writer(path: &Option<PathBuf>) -> Result<Option<BufWriter<File>>> {
    if let Some(p) = path {
        let f = File::create(p)?;
        Ok(Some(BufWriter::new(f)))
    } else {
        Ok(None)
    }
}

fn opt_pair_start(x: Option<(i64, i64)>) -> String {
    x.map(|p| p.0.to_string()).unwrap_or_default()
}
fn opt_pair_end(x: Option<(i64, i64)>) -> String {
    x.map(|p| p.1.to_string()).unwrap_or_default()
}

// ---------------------------------------------------------------------------
// Types shared with trim.rs
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug)]
pub struct ExtractPlan {
    pub norm_start: i64,
    pub norm_end: i64,
    pub orig_start: i64,
    pub orig_end: i64,
    pub strand: char,
    pub confidence: preset::Confidence,
}

#[derive(Clone, Copy, Debug)]
pub struct PlanSet {
    pub full: Option<ExtractPlan>,
    pub its1: Option<ExtractPlan>,
    pub its2: Option<ExtractPlan>,
    pub confidence: preset::Confidence,
}

fn mk_plan(
    bounds: Option<(i64, i64)>,
    strand: char,
    seq_len: i64,
    confidence: preset::Confidence,
) -> Option<ExtractPlan> {
    bounds.map(|(norm_start, norm_end)| {
        let (orig_start, orig_end) = if strand == '+' {
            (norm_start, norm_end)
        } else {
            (seq_len - norm_end + 1, seq_len - norm_start + 1)
        };
        ExtractPlan {
            norm_start,
            norm_end,
            orig_start,
            orig_end,
            strand,
            confidence,
        }
    })
}

// ---------------------------------------------------------------------------
// Anchor output types (for TSV/JSONL)
// ---------------------------------------------------------------------------

#[derive(Serialize)]
struct AnchorHitOut {
    model: String,
    strand: char,
    start: i64,
    end: i64,
    score: f64,
    evalue: f64,
}

#[derive(Serialize)]
struct BoundsOut {
    full: Option<(i64, i64)>,
    its1: Option<(i64, i64)>,
    its2: Option<(i64, i64)>,
    selected_region: String,
    selected: Option<(i64, i64)>,
}

#[derive(Serialize)]
struct AnchorRowJson {
    read_id: String,
    ssu_end: AnchorHitOut,
    s58_start: AnchorHitOut,
    s58_end: AnchorHitOut,
    lsu_start: AnchorHitOut,
    bounds: BoundsOut,
    confidence: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    ambiguous_reason: Option<String>,
}

impl AnchorHitOut {
    fn from_hit(h: &select::HitIvl) -> Self {
        AnchorHitOut {
            model: h.model.clone(),
            strand: h.strand,
            start: h.ivl.start,
            end: h.ivl.end,
            score: h.score,
            evalue: h.evalue,
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn write_anchor_tsv_row<W: Write>(
    w: &mut W,
    rid: &str,
    chain: &[select::HitIvl; 4],
    full: Option<(i64, i64)>,
    its1: Option<(i64, i64)>,
    its2: Option<(i64, i64)>,
    region: &str,
    selected: Option<(i64, i64)>,
    confidence: &preset::Confidence,
    ambiguous_reason: &str,
) -> Result<()> {
    let mut cols: Vec<String> = Vec::with_capacity(36);
    cols.push(rid.to_string());
    for h in chain {
        cols.push(h.model.clone());
        cols.push(h.strand.to_string());
        cols.push(h.ivl.start.to_string());
        cols.push(h.ivl.end.to_string());
        cols.push(h.score.to_string());
        cols.push(h.evalue.to_string());
    }
    cols.push(opt_pair_start(full));
    cols.push(opt_pair_end(full));
    cols.push(opt_pair_start(its1));
    cols.push(opt_pair_end(its1));
    cols.push(opt_pair_start(its2));
    cols.push(opt_pair_end(its2));
    cols.push(region.to_string());
    cols.push(opt_pair_start(selected));
    cols.push(opt_pair_end(selected));
    cols.push(confidence.to_string());
    cols.push(ambiguous_reason.to_string());
    debug_assert_eq!(cols.len(), 36);
    writeln!(w, "{}", cols.join("\t"))?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn write_anchor_jsonl_row<W: Write>(
    w: &mut W,
    rid: &str,
    chain: &[select::HitIvl; 4],
    full: Option<(i64, i64)>,
    its1: Option<(i64, i64)>,
    its2: Option<(i64, i64)>,
    region: &str,
    selected: Option<(i64, i64)>,
    confidence: &preset::Confidence,
    ambiguous_reason: Option<String>,
) -> Result<()> {
    let row = AnchorRowJson {
        read_id: rid.to_string(),
        ssu_end: AnchorHitOut::from_hit(&chain[0]),
        s58_start: AnchorHitOut::from_hit(&chain[1]),
        s58_end: AnchorHitOut::from_hit(&chain[2]),
        lsu_start: AnchorHitOut::from_hit(&chain[3]),
        bounds: BoundsOut {
            full,
            its1,
            its2,
            selected_region: region.to_string(),
            selected,
        },
        confidence: confidence.to_string(),
        ambiguous_reason,
    };
    let s = serde_json::to_string(&row)?;
    writeln!(w, "{}", s)?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Chain computation result
// ---------------------------------------------------------------------------

struct ChainResult {
    rid: String,
    chain: [select::HitIvl; 4],
    full: Option<(i64, i64)>,
    its1: Option<(i64, i64)>,
    its2: Option<(i64, i64)>,
    selected: Option<(i64, i64)>,
    strand: char,
    seq_len: i64,
    confidence: preset::Confidence,
    ambiguous_reason: Option<String>,
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Extract {
            input,
            hmm,
            output,
            input_format,
            output_format,
            tblout,
            tblout_existing,
            preset,
            hmmer_cpu,
            inc_e,
            region,
            max_per_anchor,
            min_its1,
            max_its1,
            min_its2,
            max_its2,
            min_full,
            max_full,
            min_anchor_score,
            max_anchor_evalue,
            anchors_tsv,
            anchors_jsonl,
            explain,
            write_skipped,
            write_ambiguous,
            qc_json,
            derep,
        } => {
            // ---------------------------------------------------------------
            // Resolve preset → effective parameters
            // ---------------------------------------------------------------
            let preset_params = match preset {
                Some(Preset::Ont) => {
                    eprintln!("Using preset: ont");
                    Some(preset::ont_preset())
                }
                Some(Preset::Hifi) => {
                    eprintln!("Using preset: hifi");
                    Some(preset::hifi_preset())
                }
                None => None,
            };

            let hardcoded_defaults = select::Constraints::default();
            let (default_conf_score, default_conf_evalue) = preset::default_confidence();

            // Helper: resolve a param with priority: explicit CLI > preset > hardcoded default
            let preset_constraints = preset_params.as_ref().and_then(|p| p.constraints.as_ref());

            let eff_inc_e = inc_e
                .or(preset_params.as_ref().and_then(|p| p.inc_e))
                .unwrap_or(1e-5);

            let eff_max_per_anchor = max_per_anchor
                .or(preset_params.as_ref().and_then(|p| p.max_per_anchor))
                .unwrap_or(8);

            let eff_constraints = select::Constraints {
                min_its1: min_its1
                    .or(preset_constraints.map(|c| c.min_its1))
                    .unwrap_or(hardcoded_defaults.min_its1),
                max_its1: max_its1
                    .or(preset_constraints.map(|c| c.max_its1))
                    .unwrap_or(hardcoded_defaults.max_its1),
                min_its2: min_its2
                    .or(preset_constraints.map(|c| c.min_its2))
                    .unwrap_or(hardcoded_defaults.min_its2),
                max_its2: max_its2
                    .or(preset_constraints.map(|c| c.max_its2))
                    .unwrap_or(hardcoded_defaults.max_its2),
                min_full: min_full
                    .or(preset_constraints.map(|c| c.min_full))
                    .unwrap_or(hardcoded_defaults.min_full),
                max_full: max_full
                    .or(preset_constraints.map(|c| c.max_full))
                    .unwrap_or(hardcoded_defaults.max_full),
            };

            let eff_min_anchor_score = min_anchor_score
                .or(preset_params.as_ref().and_then(|p| p.min_anchor_score))
                .unwrap_or(default_conf_score);

            let eff_max_anchor_evalue = max_anchor_evalue
                .or(preset_params.as_ref().and_then(|p| p.max_anchor_evalue))
                .unwrap_or(default_conf_evalue);

            // Are we filtering ambiguous reads out of the main output?
            let filter_ambiguous = write_ambiguous.is_some();

            // ---------------------------------------------------------------
            // Resolve formats
            // ---------------------------------------------------------------
            let in_fmt = match input_format {
                InputFormat::Auto => detect_format(&input).ok_or_else(|| {
                    anyhow!("Could not auto-detect input format. Use --input-format fasta|fastq")
                })?,
                other => other,
            };

            let out_fmt = match output_format {
                OutputFormat::Auto => match in_fmt {
                    InputFormat::Fastq => OutputFormat::Fastq,
                    _ => OutputFormat::Fasta,
                },
                other => other,
            };

            if matches!(in_fmt, InputFormat::Fasta) && matches!(out_fmt, OutputFormat::Fastq) {
                return Err(anyhow!(
                    "Cannot output FASTQ from FASTA input (no qualities). Use --output-format fasta."
                ));
            }

            // Convert region string -> enum
            let region_enum = match region.as_str() {
                "full" => select::Region::Full,
                "its1" => select::Region::Its1,
                "its2" => select::Region::Its2,
                "all" => select::Region::All,
                other => {
                    return Err(anyhow!(
                        "Invalid --region '{}'. Use full|its1|its2|all",
                        other
                    ));
                }
            };
            let is_all = matches!(region_enum, select::Region::All);

            eprintln!(
                "Params: in={:?} out={:?} region={} max_per_anchor={} derep={} inc_e={:.1e} \
                 its1=[{}..{}] its2=[{}..{}] full=[{}..{}] \
                 confidence: min_score={:.1} max_evalue={:.1e}",
                in_fmt,
                out_fmt,
                region,
                eff_max_per_anchor,
                derep,
                eff_inc_e,
                eff_constraints.min_its1,
                eff_constraints.max_its1,
                eff_constraints.min_its2,
                eff_constraints.max_its2,
                eff_constraints.min_full,
                eff_constraints.max_full,
                eff_min_anchor_score,
                eff_max_anchor_evalue,
            );

            // Anchor output writers
            let mut tsv_w = open_opt_writer(&anchors_tsv)?;
            let mut jsonl_w = open_opt_writer(&anchors_jsonl)?;

            // TSV header (36 columns: original 34 + confidence + ambiguous_reason)
            if let Some(w) = tsv_w.as_mut() {
                let header = [
                    "read_id",
                    "ssu_model",
                    "ssu_strand",
                    "ssu_start",
                    "ssu_end",
                    "ssu_score",
                    "ssu_evalue",
                    "s58s_model",
                    "s58s_strand",
                    "s58s_start",
                    "s58s_end",
                    "s58s_score",
                    "s58s_evalue",
                    "s58e_model",
                    "s58e_strand",
                    "s58e_start",
                    "s58e_end",
                    "s58e_score",
                    "s58e_evalue",
                    "lsu_model",
                    "lsu_strand",
                    "lsu_start",
                    "lsu_end",
                    "lsu_score",
                    "lsu_evalue",
                    "full_start",
                    "full_end",
                    "its1_start",
                    "its1_end",
                    "its2_start",
                    "its2_end",
                    "selected_region",
                    "selected_start",
                    "selected_end",
                    "confidence",
                    "ambiguous_reason",
                ]
                .join("\t");
                writeln!(w, "{}", header)?;
            }

            // ---------------------------------------------------------------
            // Resolve tblout source
            // ---------------------------------------------------------------
            let using_existing_tblout = tblout_existing.is_some();

            let (tbl_path, _tbl_tmp_guard) = if let Some(p) = tblout_existing.clone() {
                eprintln!("Using existing tblout: {}", p.display());
                (p, None)
            } else {
                if hmm.is_none() {
                    return Err(anyhow!(
                        "--hmm is required unless --tblout-existing is provided"
                    ));
                }
                if let Some(p) = tblout {
                    (p, None)
                } else {
                    let td = tempdir()?;
                    let p = td.path().join("hits.tblout");
                    (p, Some(td))
                }
            };

            // ---------------------------------------------------------------
            // Prepare FASTA target + optional dereplication
            // ---------------------------------------------------------------
            let mut derep_result: Option<derep::DerepResult> = None;

            if !using_existing_tblout {
                let hmm_path = hmm.as_ref().unwrap();

                let (fasta_target, _fasta_tmp_guard) = if derep {
                    let td = tempdir()?;
                    let tmp_fa = td.path().join("derep.fasta");

                    let dr = derep::derep_and_write_fasta(&input, in_fmt, &tmp_fa)?;
                    eprintln!(
                        "Dereplication: {} total → {} unique ({:.1}% reduction)",
                        dr.total_seqs,
                        dr.unique_seqs,
                        if dr.total_seqs > 0 {
                            (1.0 - dr.unique_seqs as f64 / dr.total_seqs as f64) * 100.0
                        } else {
                            0.0
                        }
                    );
                    derep_result = Some(dr);
                    (tmp_fa, Some(td))
                } else {
                    match in_fmt {
                        InputFormat::Fasta => (input.clone(), None),
                        InputFormat::Fastq => {
                            let td = tempdir()?;
                            let tmp_fa = td.path().join("target.fasta");
                            eprintln!("Converting FASTQ -> FASTA for nhmmer: {}", tmp_fa.display());

                            let mut r = fastq::FastqReader::new(&input)?;
                            let mut w = fasta::open_writer(&tmp_fa)?;
                            while let Some(rec) = r.next_record()? {
                                fasta::write_fasta_record(&mut w, &rec.id, &rec.seq)?;
                            }
                            (tmp_fa, Some(td))
                        }
                        InputFormat::Auto => unreachable!(),
                    }
                };

                hmmer::run_nhmmer(&hmmer::HmmerArgs {
                    hmm: hmm_path,
                    fasta: &fasta_target,
                    tblout: &tbl_path,
                    cpu: hmmer_cpu,
                    report_e: eff_inc_e,
                    inc_e: eff_inc_e,
                })?;
                eprintln!("Wrote tblout: {}", tbl_path.display());
            }

            // ---------------------------------------------------------------
            // Parse tblout
            // ---------------------------------------------------------------
            let (mut topk_map, stats) =
                tblout::stream_topk_tblout(&tbl_path, eff_max_per_anchor, Some(eff_inc_e))?;
            eprintln!(
                "tblout hits parsed: {} | anchor hits: {} | stored(topK): {} | reads w/anchor hits: {}",
                stats.total_tblout_hits, stats.anchor_hits, stats.stored_hits, stats.reads
            );

            // Expand dereplicated hits
            if let Some(ref dr) = derep_result {
                let before = topk_map.len();
                derep::expand_topk_map(&mut topk_map, dr);
                eprintln!(
                    "Derep expansion: {} representative entries → {} total entries",
                    before,
                    topk_map.len()
                );
            }

            // ---------------------------------------------------------------
            // Compute chains (parallel) + classify confidence
            // ---------------------------------------------------------------
            let chain_results: Vec<ChainResult> = topk_map
                .par_iter()
                .filter_map(|(rid, tk)| {
                    let hs = tk.flatten();
                    select::compute_chain(&hs, eff_max_per_anchor, eff_constraints).map(|chain| {
                        let strand = chain[0].strand;
                        let seq_len = chain[0].seq_len;

                        let full = select::bounds_from_chain(&chain, select::Region::Full);
                        let its1 = select::bounds_from_chain(&chain, select::Region::Its1);
                        let its2 = select::bounds_from_chain(&chain, select::Region::Its2);

                        let selected = if is_all {
                            None
                        } else {
                            select::bounds_from_chain(&chain, region_enum)
                        };

                        let confidence = preset::classify_chain(
                            &chain,
                            eff_min_anchor_score,
                            eff_max_anchor_evalue,
                        );

                        let ambiguous_reason = if confidence == preset::Confidence::Ambiguous {
                            Some(preset::ambiguous_reason(
                                &chain,
                                eff_min_anchor_score,
                                eff_max_anchor_evalue,
                            ))
                        } else {
                            None
                        };

                        ChainResult {
                            rid: rid.clone(),
                            chain,
                            full,
                            its1,
                            its2,
                            selected,
                            strand,
                            seq_len,
                            confidence,
                            ambiguous_reason,
                        }
                    })
                })
                .collect();

            // Sequential: build bounds maps + write anchor outputs
            let mut bounds_map_one: HashMap<String, ExtractPlan> = HashMap::new();
            let mut bounds_map_all: HashMap<String, PlanSet> = HashMap::new();
            let mut n_confident = 0usize;
            let mut n_ambiguous = 0usize;

            // Track ambiguous read IDs for filtering during trimming
            let mut ambiguous_ids: HashMap<String, bool> = HashMap::new();

            for cr in &chain_results {
                match cr.confidence {
                    preset::Confidence::Confident => n_confident += 1,
                    preset::Confidence::Ambiguous => n_ambiguous += 1,
                    preset::Confidence::Partial => {} // counted separately below
                }

                // If filtering ambiguous, record the ID but still build the plan
                // (we need it for --write-ambiguous output)
                if filter_ambiguous && cr.confidence == preset::Confidence::Ambiguous {
                    ambiguous_ids.insert(cr.rid.clone(), true);
                }

                if is_all {
                    let plans = PlanSet {
                        full: mk_plan(cr.full, cr.strand, cr.seq_len, cr.confidence),
                        its1: mk_plan(cr.its1, cr.strand, cr.seq_len, cr.confidence),
                        its2: mk_plan(cr.its2, cr.strand, cr.seq_len, cr.confidence),
                        confidence: cr.confidence,
                    };
                    if plans.full.is_some() || plans.its1.is_some() || plans.its2.is_some() {
                        bounds_map_all.insert(cr.rid.clone(), plans);
                    }
                } else if let Some(plan) =
                    mk_plan(cr.selected, cr.strand, cr.seq_len, cr.confidence)
                {
                    bounds_map_one.insert(cr.rid.clone(), plan);
                }

                // TSV
                if let Some(w) = tsv_w.as_mut() {
                    write_anchor_tsv_row(
                        w,
                        &cr.rid,
                        &cr.chain,
                        cr.full,
                        cr.its1,
                        cr.its2,
                        &region,
                        cr.selected,
                        &cr.confidence,
                        cr.ambiguous_reason.as_deref().unwrap_or(""),
                    )?;
                }
                // JSONL
                if let Some(w) = jsonl_w.as_mut() {
                    write_anchor_jsonl_row(
                        w,
                        &cr.rid,
                        &cr.chain,
                        cr.full,
                        cr.its1,
                        cr.its2,
                        &region,
                        cr.selected,
                        &cr.confidence,
                        cr.ambiguous_reason.clone(),
                    )?;
                }
            }

            let n_bounds = if is_all {
                bounds_map_all.len()
            } else {
                bounds_map_one.len()
            };
            let n_full_chain = n_bounds;

            // ---------------------------------------------------------------
            // Partial-chain fallback (2-anchor pairs)
            // For reads that failed the full 4-anchor chain, try extracting
            // individual regions from 2-anchor pairs:
            //   ITS1: SSU_end + 5.8S_start
            //   ITS2: 5.8S_end + LSU_start
            //   Full: SSU_end + LSU_start
            // ---------------------------------------------------------------
            let mut n_partial = 0usize;

            let full_chain_rids: std::collections::HashSet<&String> = if is_all {
                bounds_map_all.keys().collect()
            } else {
                bounds_map_one.keys().collect()
            };

            let partial_results: Vec<(String, select::PartialBounds)> = topk_map
                .par_iter()
                .filter(|(rid, _)| !full_chain_rids.contains(rid))
                .filter_map(|(rid, tk)| {
                    let hs = tk.flatten();
                    let pb =
                        select::compute_partial_bounds(&hs, eff_max_per_anchor, eff_constraints);
                    if pb.has_any() {
                        Some((rid.clone(), pb))
                    } else {
                        None
                    }
                })
                .collect();

            for (rid, pb) in &partial_results {
                n_partial += 1;
                let confidence = preset::Confidence::Partial;

                if is_all {
                    let plans = PlanSet {
                        full: pb.full.as_ref().map(|p| ExtractPlan {
                            norm_start: p.start,
                            norm_end: p.end,
                            orig_start: if p.strand == '+' {
                                p.start
                            } else {
                                p.seq_len - p.end + 1
                            },
                            orig_end: if p.strand == '+' {
                                p.end
                            } else {
                                p.seq_len - p.start + 1
                            },
                            strand: p.strand,
                            confidence,
                        }),
                        its1: pb.its1.as_ref().map(|p| ExtractPlan {
                            norm_start: p.start,
                            norm_end: p.end,
                            orig_start: if p.strand == '+' {
                                p.start
                            } else {
                                p.seq_len - p.end + 1
                            },
                            orig_end: if p.strand == '+' {
                                p.end
                            } else {
                                p.seq_len - p.start + 1
                            },
                            strand: p.strand,
                            confidence,
                        }),
                        its2: pb.its2.as_ref().map(|p| ExtractPlan {
                            norm_start: p.start,
                            norm_end: p.end,
                            orig_start: if p.strand == '+' {
                                p.start
                            } else {
                                p.seq_len - p.end + 1
                            },
                            orig_end: if p.strand == '+' {
                                p.end
                            } else {
                                p.seq_len - p.start + 1
                            },
                            strand: p.strand,
                            confidence,
                        }),
                        confidence,
                    };
                    if plans.full.is_some() || plans.its1.is_some() || plans.its2.is_some() {
                        bounds_map_all.insert(rid.clone(), plans);
                    }
                } else {
                    // For single-region mode, pick the matching partial bound
                    let pair = match region_enum {
                        select::Region::Its1 => pb.its1.as_ref(),
                        select::Region::Its2 => pb.its2.as_ref(),
                        select::Region::Full => pb.full.as_ref(),
                        select::Region::All => None, // handled by is_all branch
                    };
                    if let Some(p) = pair {
                        let plan = ExtractPlan {
                            norm_start: p.start,
                            norm_end: p.end,
                            orig_start: if p.strand == '+' {
                                p.start
                            } else {
                                p.seq_len - p.end + 1
                            },
                            orig_end: if p.strand == '+' {
                                p.end
                            } else {
                                p.seq_len - p.start + 1
                            },
                            strand: p.strand,
                            confidence,
                        };
                        bounds_map_one.insert(rid.clone(), plan);
                    }
                }
            }

            let n_bounds_total = if is_all {
                bounds_map_all.len()
            } else {
                bounds_map_one.len()
            };
            eprintln!(
                "Reads with computed bounds: {} (full-chain: {}, confident: {}, ambiguous: {}, partial: {})",
                n_bounds_total, n_full_chain, n_confident, n_ambiguous, n_partial
            );
            if filter_ambiguous {
                eprintln!(
                    "Ambiguous reads will be written separately ({} reads)",
                    ambiguous_ids.len()
                );
            }

            // ---------------------------------------------------------------
            // Writers
            // ---------------------------------------------------------------
            let mut skipped_writer: Option<trim::SeqWriter> =
                if let Some(p) = write_skipped.as_ref() {
                    eprintln!("Will write skipped reads to: {}", p.display());
                    Some(trim::SeqWriter::open(p, out_fmt)?)
                } else {
                    None
                };

            let mut ambiguous_writer: Option<trim::SeqWriter> =
                if let Some(p) = write_ambiguous.as_ref() {
                    eprintln!("Will write ambiguous reads to: {}", p.display());
                    Some(trim::SeqWriter::open(p, out_fmt)?)
                } else {
                    None
                };

            // ---------------------------------------------------------------
            // Trimming loop
            // ---------------------------------------------------------------
            let diag_region = if is_all {
                select::Region::Full
            } else {
                region_enum
            };

            let mut reader = trim::SeqReader::open(&input, in_fmt)?;

            let trim_stats = if is_all {
                let ext = out_ext(out_fmt);
                let p_full = out_path(&output, "full", ext);
                let p_its1 = out_path(&output, "its1", ext);
                let p_its2 = out_path(&output, "its2", ext);

                let s = trim::trim_all(
                    &mut reader,
                    out_fmt,
                    &p_full,
                    &p_its1,
                    &p_its2,
                    &bounds_map_all,
                    &topk_map,
                    &mut skipped_writer,
                    &mut ambiguous_writer,
                    if filter_ambiguous {
                        Some(&ambiguous_ids)
                    } else {
                        None
                    },
                    diag_region,
                    eff_max_per_anchor,
                    eff_constraints,
                    explain,
                )?;

                eprintln!("Wrote outputs:");
                eprintln!("  {}", p_full.display());
                eprintln!("  {}", p_its1.display());
                eprintln!("  {}", p_its2.display());
                eprintln!(
                    "Kept (reads with any region): {} | full: {} | its1: {} | its2: {} | partial: {} | ambiguous_out: {} | skipped: {}",
                    s.kept_any,
                    s.kept_full,
                    s.kept_its1,
                    s.kept_its2,
                    n_partial,
                    s.ambiguous_out,
                    s.skipped
                );
                s
            } else {
                let s = trim::trim_single(
                    &mut reader,
                    out_fmt,
                    &output,
                    &region,
                    &bounds_map_one,
                    &topk_map,
                    &mut skipped_writer,
                    &mut ambiguous_writer,
                    if filter_ambiguous {
                        Some(&ambiguous_ids)
                    } else {
                        None
                    },
                    diag_region,
                    eff_max_per_anchor,
                    eff_constraints,
                    explain,
                )?;

                eprintln!("Wrote output: {}", output.display());
                eprintln!(
                    "Kept: {} (partial: {})  Ambiguous (separate): {}  Skipped: {}",
                    s.kept_any, n_partial, s.ambiguous_out, s.skipped
                );
                s
            };

            // Flush
            if let Some(w) = tsv_w.as_mut() {
                w.flush()?;
            }
            if let Some(w) = jsonl_w.as_mut() {
                w.flush()?;
            }
            if let Some(w) = skipped_writer.as_mut() {
                w.flush()?;
            }
            if let Some(w) = ambiguous_writer.as_mut() {
                w.flush()?;
            }

            // ---------------------------------------------------------------
            // QC JSON summary
            // ---------------------------------------------------------------
            if let Some(ref qc_path) = qc_json {
                let summary = report::QcSummary {
                    total_reads: trim_stats.total_reads,
                    kept: report::KeptSummary {
                        total: trim_stats.kept_any,
                        confident: n_confident.min(trim_stats.kept_any),
                        ambiguous: n_ambiguous,
                        ambiguous_diverted: trim_stats.ambiguous_out,
                        partial: n_partial,
                    },
                    skipped: report::SkippedSummary {
                        total: trim_stats.skipped,
                        by_reason: trim_stats.skip_reasons.to_map(),
                    },
                    regions: if is_all {
                        Some(report::RegionSummary {
                            full: trim_stats.kept_full,
                            its1: trim_stats.kept_its1,
                            its2: trim_stats.kept_its2,
                        })
                    } else {
                        None
                    },
                    derep: derep_result.as_ref().map(|dr| report::DerepSummary {
                        total_seqs: dr.total_seqs,
                        unique_seqs: dr.unique_seqs,
                        reduction_pct: if dr.total_seqs > 0 {
                            (1.0 - dr.unique_seqs as f64 / dr.total_seqs as f64) * 100.0
                        } else {
                            0.0
                        },
                    }),
                    tblout: report::TbloutSummary {
                        total_hits: stats.total_tblout_hits,
                        anchor_hits: stats.anchor_hits,
                        stored_topk: stats.stored_hits,
                        reads_with_hits: stats.reads,
                    },
                    params: report::ParamsSummary {
                        preset: preset.map(|p| match p {
                            Preset::Ont => "ont".to_string(),
                            Preset::Hifi => "hifi".to_string(),
                        }),
                        region: region.clone(),
                        inc_e: eff_inc_e,
                        max_per_anchor: eff_max_per_anchor,
                        min_anchor_score: eff_min_anchor_score,
                        max_anchor_evalue: eff_max_anchor_evalue,
                        min_its1: eff_constraints.min_its1,
                        max_its1: eff_constraints.max_its1,
                        min_its2: eff_constraints.min_its2,
                        max_its2: eff_constraints.max_its2,
                        min_full: eff_constraints.min_full,
                        max_full: eff_constraints.max_full,
                    },
                };

                summary.write_json(qc_path)?;
                eprintln!("Wrote QC summary: {}", qc_path.display());
            }
        }
    }

    Ok(())
}
