use anyhow::{Context, Result, anyhow};
use std::path::Path;
use std::process::Command;

pub struct HmmerArgs<'a> {
    pub hmm: &'a Path,
    pub fasta: &'a Path,
    pub tblout: &'a Path,
    pub cpu: usize,

    /// Reporting threshold (nhmmer -E). Controls what gets printed to --tblout.
    pub report_e: f64,

    /// Inclusion threshold (nhmmer --incE). Used for “included” hits (mainly for pipelines),
    /// but not sufficient alone to clean up --tblout.
    pub inc_e: f64,
}

pub fn run_nhmmer(args: &HmmerArgs) -> Result<()> {
    let output = Command::new("nhmmer")
        .arg("--cpu")
        .arg(args.cpu.to_string())
        .arg("--noali")
        // Reporting threshold: this is the key for cleaning up --tblout.
        .arg("-E")
        .arg(format!("{}", args.report_e))
        // Inclusion threshold: keep same as report threshold by default.
        .arg("--incE")
        .arg(format!("{}", args.inc_e))
        .arg("--tblout")
        .arg(args.tblout)
        .arg(args.hmm)
        .arg(args.fasta)
        .output()
        .with_context(|| "Failed to spawn `nhmmer`. Is your conda env active?")?;

    if !output.status.success() {
        return Err(anyhow!(
            "nhmmer failed (exit={:?})\n--- stdout ---\n{}\n--- stderr ---\n{}",
            output.status.code(),
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr),
        ));
    }

    Ok(())
}
